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

#include "alias_sampler.h"
#include "coal_sim.h"
#include "dates.h"
#include "pop_model.h"
#include "tip_file.h"
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


struct Options {
  std::vector<std::string> args;
  uint32_t seed;
  
  // Times are measured forward in years (= 365.0*24.0 hours) from a reference time t0, by default 2020-01-01T00:00:00Z
  absl::Time t0;

  // Sampling strategy
  std::unique_ptr<Pop_model> pop_model;
  std::optional<std::string> tip_file;
  int num_samples;
  double min_tip_t, max_tip_t;

  // Substitution model (fixed to HKY, but with the right parameters, can model JC and others)
  double mu;
  double hky_kappa;
  Seq_vector<double> hky_pi_a;

  // Site-rate heterogeneity
  double site_rate_heterogeneity_alpha;

  // Genome
  int num_sites;
  
  bool output_mutation_counts;

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
      ("tip-file", "File with tip names and dates, one per line (format: Name|YYYY-MM-DD)",
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
      ("site-rate-heterogeneity-alpha",
       "Shape and rate parameter alpha for Gamma-distributed site-rate heterogeneity. "
       "Each site l gets a rate modifier nu_l ~ Gamma(shape=alpha, rate=alpha) with mean 1 and variance 1/alpha. "
       "A value of 0.0 (default) means no site-rate heterogeneity.",
       cxxopts::value<double>()->default_value("0.0"))
      ;

  options.add_options("Genome")
      ("L,num-sites", "Number of sites in the genome",
       cxxopts::value<int>())
      ;

  options.add_options("Output")
      ("output-mutation-counts", "Include per-site mutation counts in the info JSON output")
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

    auto tip_file = std::optional<std::string>{};
    auto num_samples = 0;
    auto min_tip_t = -1e6;
    auto max_tip_t = 0.0;

    if (opts.count("tip-file")) {
      // --tip-file is mutually exclusive with -n, --min-tip-t, --max-tip-t
      if (opts.count("n")) {
        fatal("Cannot specify both --tip-file and --num-samples (-n)");
      }
      if (opts.count("min-tip-t")) {
        fatal("Cannot specify both --tip-file and --min-tip-t");
      }
      if (opts.count("max-tip-t")) {
        fatal("Cannot specify both --tip-file and --max-tip-t");
      }
      tip_file = opts["tip-file"].as<std::string>();
    } else {
      if (not opts.count("n")) {
        fatal("Number of samples (-n) not specified");
      }
      num_samples = opts["n"].as<int>();
      if (num_samples <= 0) {
        fatal("Number of samples (-n) should be positive");
      }
      if (opts.count("min-tip-t")) {
        min_tip_t = parse_iso_date(opts["min-tip-t"].as<std::string>(), t0);
      }
      if (opts.count("max-tip-t")) {
        max_tip_t = parse_iso_date(opts["max-tip-t"].as<std::string>(), t0);
      }
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

    auto site_rate_heterogeneity_alpha = opts["site-rate-heterogeneity-alpha"].as<double>();
    if (site_rate_heterogeneity_alpha < 0.0) {
      fatal("Site-rate heterogeneity alpha (--site-rate-heterogeneity-alpha) must be non-negative");
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
      .tip_file = std::move(tip_file),
      .num_samples = num_samples,
      .min_tip_t = min_tip_t,
      .max_tip_t = max_tip_t,
      .mu = mu,
      .hky_kappa = hky_kappa,
      .hky_pi_a = Seq_vector{hky_pi_A, hky_pi_C, hky_pi_G, hky_pi_T},
      .site_rate_heterogeneity_alpha = site_rate_heterogeneity_alpha,
      .num_sites = num_sites,
      .output_mutation_counts = opts.count("output-mutation-counts") > 0,
      .out_info_filename = std::move(out_info),
      .out_newick_filename = std::move(out_newick),
      .out_nexus_filename = std::move(out_nexus),
      .out_fasta_filename = std::move(out_fasta),
      .out_maple_filename = std::move(out_maple)
    };
    
  } catch (cxxopts::exceptions::exception& x) {
    std::cerr << "ERROR: " << x.what() << "\n" << options.help() << "\n";
    std::exit(EXIT_FAILURE);
  } catch (std::exception& x) {
    std::cerr << "ERROR: " << x.what() << "\n";
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

auto calc_tree_height(const Phylo_tree& tree) -> double {
  auto max_tip_t = -std::numeric_limits<double>::infinity();
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      max_tip_t = std::max(max_tip_t, tree.at(node).t);
    }
  }
  return max_tip_t - tree.at_root().t;
}

auto dump_info(const Options& opts, const Phylo_tree& tree, const Global_evo_model& evo) -> std::string {

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

  auto sampling_j = json{};
  if (opts.tip_file.has_value()) {
    sampling_j = json{
      {"tip_file", opts.tip_file.value()}
    };
  } else {
    sampling_j = json{
      {"num_samples", opts.num_samples},
      {"min_tip_t", to_iso_date(opts.min_tip_t, opts.t0)},
      {"max_tip_t", to_iso_date(opts.max_tip_t, opts.t0)}
    };
  }

  auto subst_j = json{
    {"type", "hky"},
    {"mu", opts.mu},
    {"kappa", opts.hky_kappa},
    {"pi", {opts.hky_pi_a[A], opts.hky_pi_a[C], opts.hky_pi_a[G], opts.hky_pi_a[T]}}
  };

  if (opts.site_rate_heterogeneity_alpha > 0.0) {
    subst_j["site_rate_heterogeneity_alpha"] = opts.site_rate_heterogeneity_alpha;
    subst_j["nu_l"] = evo.nu_l;
  }

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
        {"total_branch_length", calc_total_branch_length(tree)},
        {"t_mrca", tree.at_root().t},
        {"t_mrca_date", to_iso_date(tree.at_root().t, opts.t0)},
        {"tree_height", calc_tree_height(tree)}}}
  };

  if (opts.output_mutation_counts) {
    auto mutation_counts = Site_vector<int>(opts.num_sites, 0);
    for (const auto& [loc, mm] : tree.mutations) {
      ++mutation_counts[mm.site];
    }
    result["mutation_counts"] = mutation_counts;
  }

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

auto name_nodes(Phylo_tree& tree, const absl::Time& t0,
                const std::vector<std::string>& tip_names = {}) -> void {
  auto num_nodes = static_cast<Node_index>(std::ssize(tree));
  auto num_tips = (num_nodes + 1) / 2;

  auto next_tip_id = 1;
  auto next_inner_node_id = num_tips + 1;

  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      if (not tip_names.empty()) {
        // Tip nodes occupy indices 0..num_tips-1 in the tree, matching the order in tip_names
        tree.at(node).name = tip_names.at(node);
      } else {
        tree.at(node).name = absl::StrFormat("TIP_%d|%s", next_tip_id, to_iso_date(tree.at(node).t, t0));
      }
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

  // We use a thinning approach combined with the alias method for O(1) site selection.
  //
  // The exact per-site mutation rate is Q_l_a(l, seq[l]) = mu * nu_l * q_a(seq[l]),
  // which depends on the current base at site l.  A static alias table built from nu_l
  // alone cannot account for the base-dependent escape rate q_a(seq[l]).
  //
  // Instead, we use a constant upper-bound rate:
  //   lambda_upper = mu * q_a_max * sum_l nu_l
  // where q_a_max = max_a q_a(a) is the maximum escape rate over all 4 bases.
  // This rate is independent of the sequence state, so it never needs updating.
  //
  // At each event:
  //   1. Draw event time from Exponential(lambda_upper)
  //   2. Pick site l from the alias table (weight proportional to nu_l), O(1)
  //   3. Accept with probability q_a(seq[l]) / q_a_max (thinning)
  //   4. If accepted, pick target base from off-diagonal q_ab row and apply
  //
  // This is a valid thinning of the Poisson process: the effective per-site rate is
  //   lambda_upper * P(pick l) * P(accept at l)
  //   = (mu * q_a_max * sum_j nu_j) * (nu_l / sum_j nu_j) * (q_a(seq[l]) / q_a_max)
  //   = mu * nu_l * q_a(seq[l])
  // which is exactly Q_l_a(l, seq[l]).

  auto site_sampler = Alias_sampler{evo.nu_l};

  // Compute q_a_max and lambda_upper
  auto q_a_max = 0.0;
  for (const auto& a : k_all_real_seq_letters) {
    q_a_max = std::max(q_a_max, evo.site_evo_model.q_a(a));
  }
  auto sum_nu_l = std::accumulate(evo.nu_l.begin(), evo.nu_l.end(), 0.0);
  auto lambda_upper = evo.site_evo_model.mu * q_a_max * sum_nu_l;

  struct Pending_node { Node_index node; Sequence_overlay seq; };

  auto work_stack = std::stack<Pending_node, std::vector<Pending_node>>{};
  work_stack.push({tree.root(), Sequence_overlay{tree.ref_sequence}});
  while (not work_stack.empty()) {
    auto [node, seq_at_node] = std::move(work_stack.top());
    work_stack.pop();

    for (const auto& child : tree.at(node).children()) {
      auto t = tree.at(node).t;
      auto t_max = tree.at(child).t;
      auto seq = seq_at_node;

      // Evolve sequence from beginning to end of branch
      while (true) {
        auto next_t = t + absl::Exponential(rng, lambda_upper);
        if (next_t >= t_max) { break; }
        t = next_t;

        // Pick a site from the alias table (weight proportional to nu_l)
        auto l = static_cast<Site_index>(site_sampler.sample(rng));

        // Thinning: accept with probability q_a(seq[l]) / q_a_max
        auto a = seq[l];
        auto q_a = evo.site_evo_model.q_a(a);
        if (absl::Uniform(absl::IntervalClosedOpen, rng, 0.0, q_a_max) >= q_a) {
          continue;  // rejected
        }

        // Pick a target state from off-diagonal q_ab row
        auto w = Seq_vector<double>{0.0};
        for (const auto& b : k_all_real_seq_letters) {
          w[b] = (a == b ? 0.0 : evo.site_evo_model.q_ab[a][b]);
        }
        auto b = pick_state(w, rng);

        tree.mutations.insert(Mutation{{child, t}, {l, a, b}});
        seq[l] = b;
      }

      if (not tree.at(child).is_tip()) {
        work_stack.push({child, seq});
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

  try {
    auto opts = process_args(argc, argv);
    std::cerr << "Seed: " << opts.seed << "\n";

    auto rng = std::mt19937{};
    rng.seed(opts.seed);

    auto tip_names = std::vector<std::string>{};
    auto tip_times = std::vector<double>{};

    if (opts.tip_file.has_value()) {
      auto tips = parse_tip_file(opts.tip_file.value(), opts.t0);
      for (auto& [name, t] : tips) {
        tip_names.push_back(std::move(name));
        tip_times.push_back(t);
      }
    } else {
      tip_times = choose_tip_times(*opts.pop_model, opts.num_samples, opts.min_tip_t, opts.max_tip_t, opts.t0, rng);
    }

    auto tree = coal_sim(*opts.pop_model, tip_times, rng);
    ladderize_tree(tree);
    name_nodes(tree, opts.t0, tip_names);
    tree.ref_sequence = gen_random_sequence(opts.num_sites, opts.hky_pi_a, rng);

    auto hky = Hky_model{};
    hky.mu = opts.mu;    // per site, per year
    hky.kappa = opts.hky_kappa;
    hky.pi_a = opts.hky_pi_a;

    auto evo = make_single_partition_global_evo_model(opts.num_sites);
    evo.site_evo_model = hky.derive_site_evo_model();

    if (opts.site_rate_heterogeneity_alpha > 0.0) {
      auto alpha = opts.site_rate_heterogeneity_alpha;
      // nu_l ~ Gamma(shape=alpha, rate=alpha), i.e., Gamma(shape=alpha, scale=1/alpha)
      // Mean = 1, variance = 1/alpha
      auto gamma_dist = std::gamma_distribution<double>{alpha, 1.0 / alpha};
      for (auto l = Site_index{0}; l != opts.num_sites; ++l) {
        evo.nu_l[l] = gamma_dist(rng);
      }
    } else {
      for (auto l = Site_index{0}; l != opts.num_sites; ++l) {
        evo.nu_l[l] = 1.0;
      }
    }

    simulate_mutations(tree, evo, rng);

    write_to(opts.out_info_filename, "Info", [&](auto& os) {
      os << dump_info(opts, tree, evo) << "\n";
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
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
