#include <iostream>
#include <fstream>
#include <functional>
#include <random>

#include "absl/log/initialize.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/strings/str_join.h"
#include "absl/time/time.h"
#include "cxxopts.hpp"
#include "nlohmann/json.hpp"

#include "coal_sim.h"
#include "pop_model.h"
#include "tree.h"
#include "version.h"
#include "evo_hky.h"
#include "sequence_overlay.h"

using json = nlohmann::json;

namespace sapling {

auto fatal(std::string_view msg) -> void {
  std::cerr << "ERROR: " << msg << "\n";
  std::exit(EXIT_FAILURE);
}

auto parse_iso_date(std::string_view iso_date_str, const absl::Time& epoch) -> double {
  auto d = absl::Time{};
  auto err = std::string{};
  if (not absl::ParseTime("%Y-%m-%d", iso_date_str, &d, &err)) {
    fatal(absl::StrFormat("Badly formatted ISO date: %s (error: %s)", iso_date_str, err));
  }
  return absl::ToDoubleHours(d - epoch) / 24.0 / 365.0;
}

auto to_iso_date(double t, const absl::Time& epoch) -> std::string {
  auto abs_time = absl::Time{epoch + absl::Hours(t * 24.0 * 365.0 + 1e-5)};  // +1e-5 to avoid roundoff
  return absl::FormatTime("%Y-%m-%d", abs_time, absl::UTCTimeZone());
}


struct Options {
  std::vector<std::string> args;
  uint32_t seed;
  
  // Times are measured forward in years (= 365.0*24.0 hours) from a reference time t0, by default 2020-01-01T00:00:00Z
  absl::Time t0;

  // Sampling strategy
  std::unique_ptr<Pop_model> pop_model;
  int num_samples;
  double min_tip_t, max_tip_t;

  // Substitution model (fixed to HKY, but with the right parameters, can model JC and others)
  double mu;
  double hky_kappa;
  Seq_vector<double> hky_pi_a;

  // Genome
  int num_sites;
  
  std::optional<std::string> out_info_filename;
  std::optional<std::string> out_newick_filename;
  std::optional<std::string> out_nexus_filename;
  std::optional<std::string> out_fasta_filename;
  std::optional<std::string> out_maple_filename;
};

auto process_args(int argc, char** argv) -> Options {
  // Record CLI invocation
  auto args = std::vector<std::string>{argv, argv+argc};
  
  cxxopts::Options options("sapling", "Sapling - Simulation of pandemic-scale sequencing data");

  options.add_options("Generic")
      ("version", "print version string")
      ("h,help", "print usage")
      ("seed", "Seed to use for random number generation (default: random)",
       cxxopts::value<uint32_t>())
      ("t0", "Epoch for measuring time (default: 2020-01-01)",
       cxxopts::value<std::string>()->default_value("2020-01-01"))
      ;

  options.add_options("Constant population model [N_e(t) = n_0]")
      ("const-pop-n0", "Effective population size n_0 (years)",
       cxxopts::value<double>())
      ;
  options.add_options("Exponential population model [N_e(t) = n_0 * exp(g * (t - t0))]")
      ("exp-pop-n0", "Effective population size at t_0 (years)",
       cxxopts::value<double>())
      ("exp-pop-g", "Exponential growth rate (per year)",
       cxxopts::value<double>())
      ;

  options.add_options("Sampling strategy")
      ("n,num-samples", "Number of tips to sample",
       cxxopts::value<int>())
      ("min-tip-t", "Earliest possible tip date (e.g., 2019-06-15; default: 1M years before t0)",
       cxxopts::value<std::string>())
      ("max-tip-t", "Latest possible tip date (e.g., 2020-01-01; default: t0)",
       cxxopts::value<std::string>())
      ;

  options.add_options("Substitution model [currently fixed to HKY]")
      ("mu", "Mutation rate per site per year (default: 1e-3)",
       cxxopts::value<double>()->default_value("1e-3"))
      ("hky-kappa", "Transition vs transversion ratio (default: 1.0)",
       cxxopts::value<double>()->default_value("1.0"))
      ("hky-pi-A", "Stationary base frequency for A (default: 0.25)",
       cxxopts::value<double>()->default_value("0.25"))
      ("hky-pi-C", "Stationary base frequency for C (default: 0.25)",
       cxxopts::value<double>()->default_value("0.25"))
      ("hky-pi-G", "Stationary base frequency for G (default: 0.25)",
       cxxopts::value<double>()->default_value("0.25"))
      ("hky-pi-T", "Stationary base frequency for T (default: 0.25)",
       cxxopts::value<double>()->default_value("0.25"))
      ;

  options.add_options("Genome")
      ("L,num-sites", "Number of sites in the genome",
       cxxopts::value<int>())
      ;

  options.add_options("Output")
      ("out-prefix", "Default filename prefix for output files, e.g., 'foo' creates 'foo.info', 'foo.fasta', etc.",
       cxxopts::value<std::string>())
      ("out-info", "Filename of output info file (default: <prefix>_info.json)",
       cxxopts::value<std::string>())
      ("out-newick", "Filename of output Newick file (default: <prefix>.nwk)",
       cxxopts::value<std::string>())
      ("out-nexus", "Filename of output Nexus file (default: <prefix>.nexus)",
       cxxopts::value<std::string>())
      ("out-fasta", "Filename of output FASTA file (default: <prefix>.fasta)",
       cxxopts::value<std::string>())
      ("out-maple", "Filename of output Maple file (default: <prefix>.maple)",
       cxxopts::value<std::string>())
      ;

  try {
    auto opts = options.parse(argc, argv);

    if (opts.count("version")) {
      std::cout << absl::StreamFormat("Sapling Version %s (build %d, commit %s)",
                                      k_sapling_version_string,
                                      k_sapling_build_number,
                                      k_sapling_commit_string) << std::endl;
      std::exit(EXIT_SUCCESS);
    }
    if (opts.count("help")) {
      std::cout << options.help() << std::endl;
      std::exit(EXIT_SUCCESS);
    }

    // Validate options
    // ================

    // Generic

    auto seed = uint32_t{};
    if (opts.count("seed")) {
      seed = opts["seed"].as<uint32_t>();
    } else {
      seed = std::random_device{}();
    }

    auto t0 = absl::Time{};
    auto t0_err = std::string{};
    if (not absl::ParseTime("%Y-%m-%d", opts["t0"].as<std::string>(), &t0, &t0_err)) {
      fatal(absl::StrFormat("Badly formatted reference time (-t0): %s (error: %s)", opts["t0"].as<std::string>(), t0_err));
    }
    
    // Population model

    auto pop_model = std::unique_ptr<Pop_model>{};
    auto check_no_pop_model_yet = [&]() {
      if (pop_model) {
        fatal("Multiple population models specified; only one should be used.");
      }
    };

    if (opts.count("const-pop-n0")) {
      check_no_pop_model_yet();
      auto n0 = opts["const-pop-n0"].as<double>();
      if (n0 <= 0.0) {
        fatal("Effective population size n0 (--const-pop-n0) should be positive.");
      }
      pop_model = std::make_unique<Const_pop_model>(n0);
    }

    if (opts.count("exp-pop-n0")) {
      check_no_pop_model_yet();
      auto n0 = opts["exp-pop-n0"].as<double>();
      auto g = opts["exp-pop-g"].as<double>();
      if (n0 <= 0.0) {
        fatal("Effective population size n0 (--exp-pop-n0) should be positive.");
      }
      if (g < 0.0) {
        fatal("Exponential growth rate (--exp-pop-g) should be non-negative.");
      }
      pop_model = std::make_unique<Exp_pop_model>(0.0, n0, g);
    }

    if (not pop_model) {
      fatal("No population model specified.");
    }

    // Sampling strategy

    if (not opts.count("n")) {
      fatal("Number of samples (-n) not specified");
    }
    auto num_samples = opts["n"].as<int>();
    if (num_samples <= 0) {
      fatal("Number of samples (-n) should be positive");
    }
    auto min_tip_t = -1e6;
    if (opts.count("min-tip-t")) {
      min_tip_t = parse_iso_date(opts["min-tip-t"].as<std::string>(), t0);
    }
    auto max_tip_t = 0.0;
    if (opts.count("max-tip-t")) {
      max_tip_t = parse_iso_date(opts["max-tip-t"].as<std::string>(), t0);
    }

    // Substitution model

    auto mu = opts["mu"].as<double>();
    if (mu <= 0.0) {
      fatal("Mutation rate (--mu) must be positive");
    }
    auto hky_kappa = opts["hky-kappa"].as<double>();
    if (hky_kappa <= 0.0) {
      fatal("HKY kappa parameter (--hky-kappa) must be positive");
    }
    auto hky_pi_A = opts["hky-pi-A"].as<double>();
    if (hky_pi_A < 0.0) {
      fatal("HKY stationary base frequency for A (--hky-pi-A) must be nonnegative");
    }
    auto hky_pi_C = opts["hky-pi-C"].as<double>();
    if (hky_pi_C < 0.0) {
      fatal("HKY stationary base frequency for C (--hky-pi-C) must be nonnegative");
    }
    auto hky_pi_G = opts["hky-pi-G"].as<double>();
    if (hky_pi_G < 0.0) {
      fatal("HKY stationary base frequency for G (--hky-pi-G) must be nonnegative");
    }
    auto hky_pi_T = opts["hky-pi-T"].as<double>();
    if (hky_pi_T < 0.0) {
      fatal("HKY stationary base frequency for T (--hky-pi-T) must be nonnegative");
    }
    auto hky_pi_sum = hky_pi_A + hky_pi_C + hky_pi_G + hky_pi_T;
    if (std::abs(hky_pi_sum - 1.0) > 1e-6) {
      fatal("HKY stationary base frequencies (--hky-pi-{A,C,G,T}) should add to 1.0");
    }

    // Genome
    auto num_sites = opts["L"].as<int>();
    
    // Output files
    auto out_prefix = std::optional<std::string>{};
    if (opts.count("out-prefix")) {
      out_prefix = opts["out-prefix"].as<std::string>();
    }
    auto complete_filename =
        [&](const std::string& opt_name, const std::string& default_suffix) -> std::optional<std::string>
        {
          if (opts.count(opt_name)) {
            return opts[opt_name].as<std::string>();
          } else if (out_prefix.has_value()) {
            return out_prefix.value() + default_suffix;
          } else {
            return {};
          }
        };
    auto out_info = complete_filename("out-info", "_info.json");
    auto out_newick = complete_filename("out-newick", ".nwk");
    auto out_nexus = complete_filename("out-nexus", ".nexus");
    auto out_fasta = complete_filename("out-fasta", ".fasta");
    auto out_maple = complete_filename("out-maple", ".maple");
    
    // Done
    
    return {
      .args = std::move(args),
      .seed = seed,
      .t0 = t0,
      .pop_model = std::move(pop_model),
      .num_samples = num_samples,
      .min_tip_t = min_tip_t,
      .max_tip_t = max_tip_t,
      .mu = mu,
      .hky_kappa = hky_kappa,
      .hky_pi_a = Seq_vector{hky_pi_A, hky_pi_C, hky_pi_G, hky_pi_T},
      .num_sites = num_sites,
      .out_info_filename = std::move(out_info),
      .out_newick_filename = std::move(out_newick),
      .out_nexus_filename = std::move(out_nexus),
      .out_fasta_filename = std::move(out_fasta),
      .out_maple_filename = std::move(out_maple)
    };
    
  } catch (cxxopts::exceptions::exception& x) {
    std::cerr << "ERROR: " << x.what() << "\n" << options.help() << "\n";
    std::exit(EXIT_FAILURE);
  }
}

auto calc_total_branch_length(const Phylo_tree& tree) -> double {
  auto T = 0.0;
  for (const auto& node : index_order_traversal(tree)) {
    if (node != tree.root()) {
      T += tree.at(node).t - tree.at_parent_of(node).t;
    }
  }
  return T;
}

auto dump_info(const Options& opts, const Phylo_tree& tree) -> std::string {

  using enum Real_seq_letter;
  
  auto pop_model_j = json{};
  if (auto pop_model = dynamic_cast<Const_pop_model*>(opts.pop_model.get()); pop_model != nullptr) {
    pop_model_j = json{
      {"type", "const"},
      {"n0", pop_model->pop()}
    };
    
  } else if (auto pop_model = dynamic_cast<Exp_pop_model*>(opts.pop_model.get()); pop_model != nullptr) {
    pop_model_j = json{
      {"type", "exp"},
      {"n0", pop_model->pop_at_t0()},
      {"g", pop_model->growth_rate()}
    };
    
  } else {
    CHECK(false) << "Unknown population model type: " << *opts.pop_model;
  }

  auto sampling_j = json{
    {"num_samples", opts.num_samples},
    {"min_tip_t", to_iso_date(opts.min_tip_t, opts.t0)},
    {"max_tip_t", to_iso_date(opts.max_tip_t, opts.t0)}
  };

  auto subst_j = json{
    {"type", "hky"},
    {"mu", opts.mu},
    {"kappa", opts.hky_kappa},
    {"pi", {opts.hky_pi_a[A], opts.hky_pi_a[C], opts.hky_pi_a[G], opts.hky_pi_a[T]}}
  };

  auto result = json{
    {"sapling_version", k_sapling_version_string},
    {"sapling_build_number", k_sapling_build_number},
    {"sapling_commit", k_sapling_commit_string},
    {"args", absl::StrJoin(opts.args, " ")},
    {"seed", opts.seed},
    {"t0", absl::FormatTime("%Y-%m-%d", opts.t0, absl::UTCTimeZone())},
    {"pop_model", pop_model_j},
    {"sampling_strategy", sampling_j},
    {"subst_model", subst_j},
    {"num_sites", opts.num_sites},
    {"tree_stats", {
        {"num_mutations", std::ssize(tree.mutations)},
        {"total_branch_length", calc_total_branch_length(tree)}}}
  };
  return result.dump(2);
}

auto choose_tip_times(
    const Pop_model& pop_model,
    int num_samples,
    double min_tip_t,
    double max_tip_t,
    const absl::Time& t0,
    absl::BitGenRef rng)
    -> std::vector<double> {

  CHECK_GT(num_samples, 0);
  CHECK_LT(min_tip_t, max_tip_t);
  
  auto result = std::vector<double>{};
  result.reserve(num_samples);

  for (auto i = 0; i != num_samples; ++i) {
    auto raw_t = pop_model.sample(min_tip_t, max_tip_t, rng);
    // Peg to nearest ISO date
    auto t = parse_iso_date(to_iso_date(raw_t, t0), t0);
    result.push_back(t);
  }

  return result;
}

auto ladderize_tree(Phylo_tree& tree) -> void {
  auto num_descendants = count_all_descendants(tree);
  for (const auto& node : post_order_traversal(tree)) {
    if (tree.at(node).is_inner_node()) {
      auto left_child = tree.at(node).left_child();
      auto right_child = tree.at(node).right_child();
      if (num_descendants[left_child] > num_descendants[right_child]) {
        // Swap
        tree.at(node).set_children(right_child, left_child);
      }
    }
  }
}

auto name_nodes(Phylo_tree& tree, const absl::Time& t0) -> void {
  auto num_nodes = static_cast<Node_index>(std::ssize(tree));
  auto num_tips = (num_nodes + 1) / 2;
  
  auto next_tip_id = 1;
  auto next_inner_node_id = num_tips + 1;

  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      tree.at(node).name = absl::StrFormat("TIP_%d|%s", next_tip_id, to_iso_date(tree.at(node).t, t0));
      ++next_tip_id;
    } else {
      tree.at(node).name = absl::StrFormat("NODE_%d|%s", next_inner_node_id, to_iso_date(tree.at(node).t, t0));
      ++next_inner_node_id;
    }
  }
}

auto gen_random_sequence(
    int num_sites,
    const Seq_vector<double>& pi_a,
    absl::BitGenRef rng)
    -> Real_sequence {
  
  auto result = Real_sequence{};
  result.reserve(num_sites);
  for (auto l = 0; l != num_sites; ++l) {
    result.push_back(pick_state(pi_a, rng));
  }
  return result;
}

auto simulate_mutations(Phylo_tree& tree, const Global_evo_model& evo, absl::BitGenRef rng) -> void {
  tree.mutations.clear();

  struct Pending_node { Node_index node; double lambda; Sequence_overlay seq; };

  // Calculate lambda at the root
  auto lambda_root = 0.0;
  for (auto l = Site_index{0}; l != tree.num_sites(); ++l) {
    lambda_root += evo.Q_l_a(l, tree.ref_sequence[l]);
  }

  // Kick off pre-order traversal while updating lambda at every turn
  auto work_stack = std::stack<Pending_node, std::vector<Pending_node>>{};
  work_stack.push({tree.root(), lambda_root, Sequence_overlay{tree.ref_sequence}});
  while (not work_stack.empty()) {
    auto [node, lambda_at_node, seq_at_node] = std::move(work_stack.top());
    work_stack.pop();

    for (const auto& child : tree.at(node).children()) {
      auto lambda = lambda_at_node;
      auto t = tree.at(node).t;
      auto t_max = tree.at(child).t;
      auto seq = seq_at_node;

      // Evolve sequence from beginning to end of branch
      while (true) {
        // Pick a time
        auto next_t = t + absl::Exponential(rng, lambda);
        if (next_t >= t_max) { break; }

        // Pick a site (TODO: need to do this *way* more efficiently using a Fenwick Tree)
        // e.g.: https://medium.com/@sandeepsign/binary-index-tree-efficient-accumulative-sum-or-count-3cae9b83e41a
        // or here: https://en.wikipedia.org/wiki/Fenwick_tree
        auto target_cum_Q = absl::Uniform(absl::IntervalClosedOpen, rng, 0.0, lambda);
        auto l = Site_index{0};
        auto cum_Q = evo.Q_l_a(0, seq[0]);
        while (l != tree.num_sites() && cum_Q < target_cum_Q) {
          // Edge case when mutating last site but rounding errors make cum_Q < target_cum_Q
          if ((l+1) == tree.num_sites()) { break; }
          
          ++l;
          cum_Q += evo.Q_l_a(l, seq[l]);
        }
        CHECK_GE(l, 0);
        CHECK_LT(l, tree.num_sites());

        // Pick a target state
        auto a = seq[l];
        auto w = Seq_vector<double>{0.0};
        for (const auto& b : k_all_real_seq_letters) {
          w[b] = (a == b ? 0.0 : evo.Q_l_ab(l, a, b));
        }
        auto b = pick_state(w, rng);

        // Do it
        tree.mutations.insert(Mutation{{child, next_t}, {l, a, b}});
        seq[l] = b;
        t = next_t;
        lambda += evo.Q_l_a(l, b) - evo.Q_l_a(l, a);
      }

      if (not tree.at(child).is_tip()) {
        work_stack.push({child, lambda, seq});
      }
    }
  }
}

auto output_newick_tree(std::ostream& os, const Phylo_tree& tree, const absl::Time& t0, bool annotated) -> void {
  if (annotated) {
    // Mark as rooted tree and print out sequence at root
    os << "[&R] ";
    os << "[&root_seq=";
    for (const auto& a : tree.ref_sequence) {
     os << to_char(a);
    }
    os << "] ";
  }

  auto seq = Sequence_overlay{tree.ref_sequence};
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      if (tree.at(node).is_inner_node()) {
        os << "(";
      }

      // Apply mutations on branch leading up to here
      for (const auto& [loc, mm] : on_branch(node, tree.mutations)) {
        CHECK_EQ(seq[mm.site], mm.from);
        seq[mm.site] = mm.to;
      }
    }

    if (children_so_far > 0 && children_so_far < std::ssize(tree.at(node).children())) {
      // In between two children of an inner node
      os << ",";
    }

    if (children_so_far == std::ssize(tree.at(node).children())) {
      // Exiting `node`
      if (tree.at(node).is_inner_node()) {
        os << ")";
      }
      os << tree.at(node).name;

      if (annotated) {
        if (not tree.at(node).is_tip()) {
          os << "[&date=" << to_iso_date(tree.at(node).t, t0) << "]";
        }
      }
      
      os << ":";

      auto t_parent = (node == tree.root() ? tree.at_root().t : tree.at_parent_of(node).t);
      if (annotated) {
        // Add mutations, if any, as branch attribute
        auto mutations = std::vector<std::string>{};
        for (const auto& [loc, mm] : on_branch(node, tree.mutations)) {
          mutations.push_back(absl::StrFormat(
              "%c%d%c,%g", to_char(mm.from), mm.site+1, to_char(mm.to), loc.t - t_parent));
        }
        if (not mutations.empty()) {
          os << "[&mutations={" << absl::StrJoin(mutations, ",") << "}]";
        }
      }
      os << tree.at(node).t - t_parent;

      // Revert mutations on branch leading up to here
      for (const auto& [loc, mm] : on_branch(node, tree.mutations) | std::views::reverse) {
        CHECK_EQ(seq[mm.site], mm.to);
        seq[mm.site] = mm.from;
      }
    }
  }
  os << ";";
}

auto output_fasta(std::ostream& os, const Phylo_tree& tree) -> void {
  auto seq = Sequence_overlay{tree.ref_sequence};
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      // Apply mutations on branch leading up to here
      for (const auto& [loc, mm] : on_branch(node, tree.mutations)) {
        CHECK_EQ(seq[mm.site], mm.from);
        seq[mm.site] = mm.to;
      }
    }

    if (children_so_far == std::ssize(tree.at(node).children())) {
      // Exiting `node`
      if (tree.at(node).is_tip()) {
        os << ">" << tree.at(node).name << "\n";
        for (const auto& a : seq) {
          os << to_char(a);
        }
        os << "\n";
      }

      // Revert mutations on branch leading up to here
      for (const auto& [loc, mm] : on_branch(node, tree.mutations) | std::views::reverse) {
        CHECK_EQ(seq[mm.site], mm.to);
        seq[mm.site] = mm.from;
      }
    }
  }
}

auto output_maple(std::ostream& os, const Phylo_tree& tree) -> void {
  os << ">reference\n";
  for (const auto& a : tree.ref_sequence) {
    os << to_char(a);
  }
  os << "\n";
  
  auto seq = Sequence_overlay{tree.ref_sequence};
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      // Entering `node`
      // Apply mutations on branch leading up to here
      for (const auto& [loc, mm] : on_branch(node, tree.mutations)) {
        CHECK_EQ(seq[mm.site], mm.from);
        seq[mm.site] = mm.to;
      }
    }

    if (children_so_far == std::ssize(tree.at(node).children())) {
      // Exiting `node`
      if (tree.at(node).is_tip()) {
        os << ">" << tree.at(node).name << "\n";
        auto deltas = std::vector<std::pair<Site_index, Real_seq_letter>>{seq.deltas().begin(), seq.deltas().end()};
        std::ranges::sort(deltas);
        for (const auto& [l, b] : deltas) {
          os << to_char(b) << "\t" << (l+1) << "\n";
        }
      }

      // Revert mutations on branch leading up to here
      for (const auto& [loc, mm] : on_branch(node, tree.mutations) | std::views::reverse) {
        CHECK_EQ(seq[mm.site], mm.to);
        seq[mm.site] = mm.from;
      }
    }
  }
}

auto write_to(const std::optional<std::string> maybe_filename, std::string_view desc, auto writer) -> void {
  if (not maybe_filename.has_value()) {
    return;
  }
  
  auto os = std::ofstream{std::string{maybe_filename.value()}};
  if (not os) {
    fatal(absl::StrFormat("Could not write %s file '%s'", desc, maybe_filename.value()));
  }
  writer(os);
  os.close();
  std::cerr << "Wrote " << desc << " to " << maybe_filename.value() << '\n';
}

}  // namespace sapling

auto main(int argc, char** argv) -> int {
  using namespace sapling;
  
  absl::InitializeLog();

  auto opts = process_args(argc, argv);
  std::cerr << "Seed: " << opts.seed << "\n";

  auto rng = std::mt19937{};
  rng.seed(opts.seed);

  auto tip_times = choose_tip_times(*opts.pop_model, opts.num_samples, opts.min_tip_t, opts.max_tip_t, opts.t0, rng);

  auto tree = coal_sim(*opts.pop_model, tip_times, rng);
  ladderize_tree(tree);
  name_nodes(tree, opts.t0);
  tree.ref_sequence = gen_random_sequence(opts.num_sites, opts.hky_pi_a, rng);

  auto hky = Hky_model{};
  hky.mu = opts.mu;    // per site, per year
  hky.kappa = opts.hky_kappa;
  hky.pi_a = opts.hky_pi_a;
  
  auto evo = make_single_partition_global_evo_model(opts.num_sites);
  evo.site_evo_model = hky.derive_site_evo_model();

  // TODO: Fill in evo.nu_l[l] with something other than 1.0
  
  simulate_mutations(tree, evo, rng);
  
  write_to(opts.out_info_filename, "Info", [&](auto& os) {
    os << dump_info(opts, tree) << "\n";
  });

  write_to(opts.out_newick_filename, "Newick", [&](auto& os) {
    output_newick_tree(os, tree, opts.t0, false);
  });

  write_to(opts.out_nexus_filename, "Nexus", [&](auto& os) {
    os << "#NEXUS\n"
       << "\n"
       << "Begin trees;\n"
       << "tree TREE1 = ";
    output_newick_tree(os, tree, opts.t0, true);
    os << "\nEnd;\n";
  });

  write_to(opts.out_fasta_filename, "FASTA", [&](auto& os) {
    output_fasta(os, tree);
  });

  write_to(opts.out_maple_filename, "Maple", [&](auto& os) {
    output_maple(os, tree);
  });

  std::cerr << "Total number of mutations: " << tree.mutations.size() << "\n";
  std::cerr << "Total branch length: " << calc_total_branch_length(tree) << " years\n";
  
  return EXIT_SUCCESS;
}
