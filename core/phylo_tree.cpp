#include "phylo_tree.h"

namespace sapling {

auto operator<<(std::ostream& os, const Phylo_tree& tree) -> std::ostream& {
  os << "Ref sequence: '" << absl::StrJoin(tree.ref_sequence, "", absl::StreamFormatter()) << "'\n";
  for (const auto& node : index_order_traversal(tree)) {
    os << (node == tree.root() ? '*' : ' ') << absl::StreamFormat("[%3d] ", node) << tree.at(node) << "\n";
  }
  return os;
}

auto assert_phylo_tree_integrity(const Phylo_tree& tree, bool force) -> void {
  if (estd::is_debug_enabled || force) {
    assert_tree_integrity(tree, force);

    for (const auto& node : pre_order_traversal(tree)) {
      if (tree.at(node).is_inner_node()) {
        for (const auto& child : tree.at(node).children()) {
          CHECK_LE(tree.at(node).t, tree.at(child).t);
        }
      }
    }
  }
}

}  // namespace sapling
