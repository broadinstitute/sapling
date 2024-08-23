#ifndef SAPLING_PHYLO_TREE_H_
#define SAPLING_PHYLO_TREE_H_

#include "tree.h"
#include "sequence.h"
#include "mutations.h"

namespace sapling {

struct Phylo_node : public Binary_node {
  std::string name;
  double t;

  auto operator<=>(const Phylo_node& that) const = default;
};
inline auto operator<<(std::ostream& os, const Phylo_node& node) -> std::ostream& {
  return os << absl::StreamFormat(
      "Phylo_node{name=\"%s\", t=%f; parent=%d, children=[%s]}",
      node.name,
      node.t,
      node.parent(),
      absl::StrJoin(node.children(), ", "));
}

struct Phylo_tree : public Tree<Phylo_node> {
  using Tree::Tree;

  Real_sequence ref_sequence;
  Mutations mutations;

  auto num_sites() const -> Site_index { return std::ssize(ref_sequence); }
};
auto operator<<(std::ostream& os, const Phylo_tree& tree) -> std::ostream&;

// Assertions

auto assert_phylo_tree_integrity(const Phylo_tree& tree, bool force = false) -> void;

}  // namespace sapling

#endif // SAPLING_PHYLO_TREE_H_
