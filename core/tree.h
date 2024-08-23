#ifndef SAPLING_TREE_H_
#define SAPLING_TREE_H_

#include <vector>
#include <ranges>
#include <stack>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/log/check.h"
#include "cppcoro/generator.hpp"

#include "estd.h"

namespace sapling {

// Trees and Nodes
// ===============

// Trees consist of a contiguous list of nodes.  Internally and externally,
// nodes are referred to by their index into such a list, using `Node_index`.  This representation
// is compact and easy to transfer across a process boundary or the network.  It also makes
// it easy to externally decorate all or some nodes on a tree (see `Node_vector<T>` and `Node_map<T>` below).
//
// By convention, variables holding indices have names like `node`, not `node_index`.
// By convention, the node objects themselves are almost never given a name; we prefer instead
// indirect access via a tree, e.g., to get the parent node of node i, write `tree.at(i).parent()`.
using Node_index = int;

// We refer to branches by their endpoint nodes.  Hence, branch `i` joins node `tree.at(i).parent()` to `i`.
using Branch_index = Node_index;

// Sentinel value for a blank reference to a node (like `nullptr`)
inline constexpr Node_index k_no_node = -1;

// Stores a `T` for every node in a tree, indexed by a node's index
template<typename T, typename Alloc = std::allocator<T>>
using Node_vector = std::vector<T, Alloc>;

// Stores a `T` for a small number of nodes in a tree, indexed by a node's index
template<typename T>
using Node_map = absl::flat_hash_map<Node_index, T>;

// Stores a set of nodes in a tree, represented by their index
using Node_set = absl::flat_hash_set<Node_index>;

// Nodes are the basic building block for trees.  Node<C> encodes topological info for navigating trees,
// using the container C for the child indices.  Concrete node classes derive from Node<C>.
// The `Node_like` concept is used to restrict template parameters to descendants of Node<C> for some C.
//
// Most commonly, concrete nodes subclass N one of two specializations:
//
// - Binary_node (compact representation for when inner nodes have at most two children)
// - Nary_node (no restrictions on number of children per node)
//
// You can then make trees of such nodes by deriving from Tree<N> (or using it directly if there is
// no tree-level information other than the root index).
//
// Example:
//
//  struct My_node : public Binary_node {
//    std::string name;
//    double coeff;
//  };
//
//  auto tree = Tree<My_node>{};
//  //...
//  auto node = tree.add_node();
//  tree.at(node).coeff = 42.0;
//
template<typename Child_indices>
class Node {
 public:
  auto parent() const -> Node_index { return parent_; }
  auto set_parent(Node_index parent) -> void { parent_ = parent; }

  auto children() -> Child_indices& { return children_; }
  auto children() const -> const Child_indices& { return children_; }
  auto clear_children() -> void { children_.clear(); }
  auto add_child(Node_index child) -> void { children_.push_back(child); }

  auto is_inner_node() const -> bool { return not children().empty(); }
  auto is_tip() const -> bool { return children().empty(); }

  auto operator<=>(const Node& that) const = default;

 private:
  Node_index parent_{k_no_node};
  Child_indices children_{};
};
template<typename N>
concept Node_like = requires(N node){
  // Adapted from https://stackoverflow.com/a/71921982 : N should be descended from Node<C> for some C
  []<typename C>(const Node<C>&){}(node);
};

auto operator<<(std::ostream& os, const Node_like auto& node) -> std::ostream& {
  return os << absl::StreamFormat(
      "Node{parent_index=%d, child_indices=[%s]}",
      node.parent_index(),
      absl::StrJoin(node.child_indices(), ", "));
}

// A container with up to 2 child indices.  Like a std::vector<Node_index>, but stored inline
// and restricted to up to 2 elements (neither of which is k_no_node).  Code is modelled on std::array.
// Casing of names follows std for integration with generic code
class Binary_child_indices {
  using Backing = std::array<Node_index, 2>;

 public:
  using value_type = Backing::value_type;
  using pointer = Backing::pointer;
  using const_pointer = Backing::const_pointer;
  using reference = Backing::reference;
  using const_reference = Backing::const_reference;
  using iterator = Backing::iterator;
  using const_iterator = Backing::const_iterator;
  using size_type = Backing::size_type;
  using difference_type = Backing::difference_type;
  using reverse_iterator = Backing::reverse_iterator;
  using const_reverse_iterator = Backing::const_reverse_iterator;

  // Constructors and assignments
  Binary_child_indices() : indices_{k_no_node, k_no_node} {}
  Binary_child_indices(Node_index left_index, Node_index right_index) { assign(left_index, right_index); }
  Binary_child_indices(const Binary_child_indices& other) = default;
  auto operator=(const Binary_child_indices& other) -> Binary_child_indices& = default;
  Binary_child_indices(Binary_child_indices&& other) noexcept = default;
  auto operator=(Binary_child_indices&& other) noexcept -> Binary_child_indices& = default;

  auto swap(Binary_child_indices& other) noexcept -> void { this->indices_.swap(other.indices_); }

  auto assign(Node_index left_index, Node_index right_index) -> void {
    DCHECK(left_index != k_no_node);
    DCHECK(right_index != k_no_node);
    indices_ = {left_index, right_index};
  }

  // Iterators.  Note that empty lists correctly have begin == end
  auto begin() noexcept -> iterator { return indices_.begin(); }
  auto begin() const noexcept -> const_iterator { return indices_.begin(); }
  auto end() noexcept -> iterator { return std::ranges::find(indices_, k_no_node); }
  auto end() const noexcept -> const_iterator { return std::ranges::find(indices_, k_no_node); }

  // Capacity
  auto size() const noexcept -> size_type { return end() - begin(); }
  constexpr auto max_size() const noexcept -> size_type { return indices_.max_size(); }
  auto empty() const -> bool { return begin() == end(); }

  // Element access
  auto operator[](size_type pos) noexcept -> reference { return indices_[pos]; }
  auto operator[](size_type pos) const noexcept -> const_reference { return indices_[pos]; }
  auto at(size_type pos) -> reference {
    if (pos >= size()) { throw std::out_of_range("Invalid pos"); }
    return indices_[pos];
  }
  auto at(size_type pos) const -> const_reference {
    if (pos >= size()) { throw std::out_of_range("Invalid pos"); }
    return indices_[pos];
  }
  auto data() const noexcept -> const_pointer { return indices_.data(); }

  auto clear() -> void { indices_ = {k_no_node, k_no_node}; }
  auto push_back(Node_index index) -> void {
    DCHECK(index != k_no_node);
    DCHECK_LT(std::ssize(*this), 2);
    indices_[size()] = index;  // This automatically increases the size (one fewer element has value `k_no_node`)
  }

  // Comparisons
  auto operator<=>(const Binary_child_indices& that) const -> std::strong_ordering {
    for (auto i = 0; i != std::ssize(indices_); ++i) {
      if (auto cmp = (this->indices_[i] <=> that.indices_[i]); cmp != 0) {
        return cmp;
      }
    }
    return std::strong_ordering::equal;
  }

 private:
  // Only the indices up to the first value k_no_node are valid
  Backing indices_;
};
inline auto swap(Binary_child_indices& one, Binary_child_indices& two) -> void { one.swap(two); }

class Binary_node : public Node<Binary_child_indices> {
 public:
  auto set_children(Node_index left_child, Node_index right_child) -> void {
    children().assign(left_child, right_child);
  }
  
  auto left_child() const -> Node_index {
    DCHECK_EQ(std::ssize(children()), 2);
    return children()[0];
  }
  auto right_child() const -> Node_index {
    DCHECK_EQ(std::ssize(children()), 2);
    return children()[1];
  }
  auto sibling_of(Node_index X) const -> Node_index {
    CHECK(X == left_child() || X == right_child()) << X << "," << left_child() << "," << right_child();
    return X == left_child() ? right_child() : left_child();
  }
  auto operator<=>(const Binary_node& that) const = default;
};

using Nary_node = Node<std::vector<Node_index>>;


// A Tree groups together a set of nodes, one of which is a root, and operations
// to access data for those nodes (`tree.at(node)`).  Concrete trees derive
// from Tree<N,A>, whose template parameters are the node class and an allocator.
// The `Tree_like` concept is used to restrict template parameters to descendants of Tree<N,A> for some N and A.
template<Node_like Node, template<typename> typename Alloc = std::allocator>
class Tree {
 public:
  explicit(false) Tree(const Alloc<Node>& alloc = {}) : root_{k_no_node}, nodes_{alloc} {}
  explicit Tree(Node_vector<Node, Alloc<Node>> nodes) : root_{k_no_node}, nodes_{std::move(nodes)} {}
  explicit Tree(Node_index num_nodes, const Alloc<Node>& alloc = {})
      : root_{k_no_node}, nodes_(num_nodes, Node{}, alloc) {}

  auto root() const -> Node_index { return root_; }
  auto set_root(Node_index root) -> void {
    DCHECK(root == k_no_node || (root >= 0 && root < std::ssize(nodes_)));
    root_ = root;
  }
  auto size() const -> Node_index { return std::ssize(nodes_); }

  auto at(Node_index i) -> Node& { return nodes_.at(i); }
  auto at(Node_index i) const -> const Node& { return nodes_.at(i); }

  // Often enough, we need to access the parent node of i.
  // Use tree.at_parent_of(i) instead of tree.at(tree.at(i).parent()) for clarity
  auto at_parent_of(Node_index i) -> Node& { return at(at(i).parent()); }
  auto at_parent_of(Node_index i) const -> const Node& { return at(at(i).parent()); }

  // Likewise for the root
  auto at_root() -> Node& { return at(root()); }
  auto at_root() const -> const Node& { return at(root()); }

  auto add_node(const Node& node) -> Node_index {
    auto new_node_index = std::ssize(nodes_);
    nodes_.push_back(node);
    return new_node_index;
  }
  auto add_node(Node&& node = {}) -> Node_index {
    auto new_node_index = std::ssize(nodes_);
    nodes_.push_back(std::move(node));
    return new_node_index;
  }

 private:
  Node_index root_;
  Node_vector<Node, Alloc<Node>> nodes_;
};
template<typename T>
concept Tree_like = requires(T tree){
  // Adapted from https://stackoverflow.com/a/71921982 : T should be descended from Tree<N,A> for some N,A
  []<typename N, template<typename> typename A>(const Tree<N,A>&){}(tree);
};


// Tree Traversals
// ===============

// Perform a depth-first walk through a tree, yielding node every time it is visited even transiently.
// The standard pre-order, in-order and post-order traversals are specializations of this walk (see below).
// The general walk is something useful if you need to update information on entering and/or exiting every node.
namespace detail {
struct Node_visitation { Node_index node; int children_so_far; };
};
auto traversal(const Tree_like auto& tree) -> cppcoro::generator<detail::Node_visitation> {
  if (std::ssize(tree) == 0) { co_return; }
  
  auto work_stack = std::stack<detail::Node_visitation>{};
  work_stack.emplace(tree.root(), -1);
  
  while (not work_stack.empty()) {
    auto [node, children_so_far] = work_stack.top();
    work_stack.pop();

    if (children_so_far != -1) {
      co_yield detail::Node_visitation{node, children_so_far};
    } else {
      const auto& children = tree.at(node).children();
      auto num_children = std::ssize(children);
      work_stack.emplace(node, num_children);
      for (auto i = num_children - 1; i >= 0; --i) {
        work_stack.emplace(children[i], -1);
        work_stack.emplace(node, i);
      }
    }
  }
}

auto pre_order_traversal(const Tree_like auto& tree) -> cppcoro::generator<Node_index> {
  // auto work = std::stack<Node_index>{};
  // work.push(tree.root());
  // while (not work.empty()) {
  //   auto node = work.top();
  //   work.pop();

  //   co_yield node;

  //   for (const auto& child : tree.at(node).children() | std::views::reverse) {
  //     work.push(child);
  //   }
  // }
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == 0) {
      co_yield Node_index{node};
    }
  }
}

auto in_order_traversal(const Tree_like auto& tree) -> cppcoro::generator<Node_index> {
  // struct Pending_node { Node_index node; bool visit_now; };
  // auto work_stack = std::stack<Pending_node>{};
  // work_stack.emplace(tree.root(), false);
  
  // while (not work_stack.empty()) {
  //   auto [node, visit_now] = work_stack.top();
  //   work_stack.pop();

  //   visit_now |= tree.at(node).is_tip();
  //   if (visit_now) { co_yield node; }
  //   else {
  //     const auto& children = tree.at(node).children();
  //     work_stack.emplace(children[1], false);
  //     work_stack.emplace(node, true);
  //     work_stack.emplace(children[0], false);
  //   }
  // }
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == std::ssize(tree.at(node).children()) / 2) {
      co_yield Node_index{node};
    }
  }
}

auto post_order_traversal(const Tree_like auto& tree) -> cppcoro::generator<Node_index> {
  for (const auto& [node, children_so_far] : traversal(tree)) {
    if (children_so_far == std::ssize(tree.at(node).children())) {
      co_yield Node_index{node};
    }
  }
}

auto index_order_traversal(const Tree_like auto& tree) {
  return std::views::iota(0, std::ssize(tree));
}


// Assertions
// ==========

inline auto assert_node_integrity(const Tree_like auto& tree, Node_index node, bool force = false) -> void {
  if (estd::is_debug_enabled || force) {
    auto parent = tree.at(node).parent();
    if (parent != k_no_node) {
      CHECK_GE(parent, 0);
      CHECK_LT(parent, std::ssize(tree));

      CHECK_GT(std::ranges::count(tree.at(parent).children(), node), 0);
    }

    for (const auto& child : tree.at(node).children()) {
      CHECK_GE(child, 0);
      CHECK_LT(child, std::ssize(tree));

      CHECK_EQ(tree.at(child).parent(), node);
    }
  }
}

inline auto assert_tree_integrity(const Tree_like auto& tree, bool force = false) -> void {
  if (estd::is_debug_enabled || force) {
    if (std::ssize(tree) == 0) {
      CHECK_EQ(tree.root(), k_no_node);
      return;
    }

    auto root = tree.root();
    CHECK_NE(root, k_no_node);
    CHECK_GE(root, 0);
    CHECK_LT(root, std::ssize(tree));

    auto visited = Node_vector<bool>(std::ssize(tree), false);
    for (const auto& node : pre_order_traversal(tree)) {
      CHECK(!visited[node]) << node;
      visited[node] = true;
      assert_node_integrity(tree, node, force);
    }

    CHECK(std::ranges::all_of(visited, [](bool v) { return v == true; }))
        << std::ranges::count(visited, true) << " vs " << std::ssize(visited);
  }
}

// Generic operations
// ==================

auto copy_topology(const Tree_like auto& src_tree, Tree_like auto& dst_tree) -> void {
  CHECK_EQ(std::ssize(src_tree), std::ssize(dst_tree));

  for (const auto& src_node : index_order_traversal(src_tree)) {
    auto dst_node = src_node;
    dst_tree.at(dst_node).clear_children();
    for (const auto& src_child : src_tree.at(src_node).children()) {
      auto dst_child = src_child;
      dst_tree.at(dst_node).add_child(dst_child);
    }
    auto dst_parent = src_tree.at(src_node).parent();
    dst_tree.at(dst_node).set_parent(dst_parent);
  }
  dst_tree.set_root(src_tree.root());
}

// result[i] is the number of descendants of i
auto count_all_descendants(const Tree_like auto& tree) -> Node_vector<int> {
  auto descendants = Node_vector<int>(std::ssize(tree), 0);
  for (const auto& node : post_order_traversal(tree)) {  // postorder because we aggregate from tips to root
    auto node_descendants = 0;
    for (const auto& child : tree.at(node).children()) {
      node_descendants += 1 + descendants[child];
    }
    descendants[node] = node_descendants;
  }
  return descendants;
}

}  // namespace sapling

#endif // SAPLING_TREE_H_
