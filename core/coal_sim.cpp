#include "coal_sim.h"

#include <ranges>

#include "absl/random/distributions.h"

namespace sapling {

auto coal_sim(const Pop_model& pop_model, std::span<const double> tip_times, absl::BitGenRef rng) -> Phylo_tree {
  auto num_tips = static_cast<Node_index>(std::ssize(tip_times));
  auto num_nodes = 2*num_tips - 1;
  auto tree = Phylo_tree{num_nodes};

  for (auto i = Node_index{0}; i != num_tips; ++i) {
    tree.at(i).t = tip_times[i];
  }

  // We advance from the maximum time backwards.  At any time t:
  // * active_branches contains branch i (i.e., ending at node i) if it is active
  // * pending_tips is a stack of tips having times earlier than t, where `top` returns the latest one
  // * c is the next free inner node index
  auto t = std::ranges::max(tip_times);
  auto active_branches = std::vector<Node_index>{};
  auto pending_tips = std::vector<Node_index>{};
  for (auto i = Node_index{0}; i != num_tips; ++i) {
    pending_tips.push_back(i);
  }
  std::ranges::sort(pending_tips, {}, [&tree](Node_index i) { return tree.at(i).t; });
  auto c = Node_index{num_tips};

  while (not (pending_tips.empty() && std::ssize(active_branches) < 2)) {
    
    // Should the latest node before t be a coalescence or a tip?
    
    auto next_t_tip = double{};
    if (pending_tips.empty()) {
      next_t_tip = -std::numeric_limits<double>::max();
    } else {
      next_t_tip = tree.at(pending_tips.back()).t;
    }
    
    auto next_t_coal = double{};
    if (std::ssize(active_branches) < 2) {
      next_t_coal = -std::numeric_limits<double>::max();
    } else {
      auto k = std::ssize(active_branches);
      auto k_choose_2 = (k * (k-1)) / 2;
      auto delta_intensity = absl::Exponential<double>(rng, k_choose_2);
      next_t_coal = pop_model.inverse_intensity(pop_model.intensity_at_time(t) - delta_intensity);
    }

    if (next_t_tip > next_t_coal) {
      // Tip
      auto i = pending_tips.back();
      CHECK_EQ(next_t_tip, tree.at(i).t);
      active_branches.push_back(i);
      pending_tips.pop_back();
      
      t = next_t_tip;
    } else {
      // Coalesce between two random active_branches
      CHECK_GE(std::ssize(active_branches), 2);
      auto ii = absl::Uniform<Node_index>(absl::IntervalClosedOpen, rng, 0, std::ssize(active_branches));
      auto jj = ii;
      while (jj == ii) {
        jj = absl::Uniform<Node_index>(absl::IntervalClosedOpen, rng, 0, std::ssize(active_branches));
      }
      if (ii > jj) { std::swap(ii, jj); }

      auto i = active_branches[ii];
      auto j = active_branches[jj];

      // Remove elements ii & jj from active_branches (order of the remaining items is irrelevant)
      // This is slightly tricky!
      CHECK_NE(ii, std::ssize(active_branches) - 1);
      active_branches[jj] = active_branches.back();
      active_branches.pop_back();
      active_branches[ii] = active_branches.back();
      active_branches.pop_back();

      tree.at(c).t = next_t_coal;
      tree.at(c).set_children(i, j);
      tree.at(i).set_parent(c);
      tree.at(j).set_parent(c);
      active_branches.push_back(c);
      ++c;

      t = next_t_coal;
    }
  }

  // Should be done and the remaining active branch should be the root of the tree
  CHECK(pending_tips.empty());
  CHECK_EQ(c, num_nodes);
  CHECK_EQ(std::ssize(active_branches), 1);
  
  auto root = active_branches[0];
  tree.at(root).set_parent(k_no_node);
  tree.set_root(root);

  return tree;
}

}  // namespace sapling
