#ifndef SAPLING_MUTATIONS_H_
#define SAPLING_MUTATIONS_H_

#include <iostream>

#include "absl/container/btree_map.h"
#include "absl/strings/str_format.h"

#include "sequence.h"
#include "tree.h"

namespace sapling {


// Tree annotations: info objects hanging off specific locations on a tree.
// A Tree_loc specifies a point on a tree in terms of a branch and a time on that branch.
// A given tree location can have either single annotations (e.g., location changes)
// or multiple annotations (e.g., mutations at different sites).  For the former, use a Tree_annotation_map;
// for the latter, use a Tree_annotation_multimap.
struct Tree_loc {
  Branch_index branch;
  double t;
  auto operator<=>(const Tree_loc& that) const = default;
};
inline auto operator<<(std::ostream& os, const Tree_loc& loc) -> std::ostream& {
  return os << absl::StreamFormat("[[%d,%g]]", loc.branch, loc.t);
}

template<typename X>
using Tree_annotation = std::pair<const Tree_loc, X>;

inline auto operator<<(std::ostream& os, const Tree_annotation<auto>& annotation) -> std::ostream& {
  const auto& [loc, info] = annotation;
  return os << absl::StreamFormat("%s:%s", absl::FormatStreamed(loc), absl::FormatStreamed(info));
}

// Tree_loc_cmp enables heterogeneous lookup of tree locations based on branch and on full tree location
struct Tree_loc_cmp {
  using is_transparent = void;  // Marker to enable heterogeneous lookup

  auto operator()(const Tree_loc& lhs_loc, Branch_index rhs_branch) const {
    return lhs_loc.branch < rhs_branch;
  }
  auto operator()(Branch_index lhs_branch, const Tree_loc& rhs_loc) const {
    return lhs_branch < rhs_loc.branch;
  }

  auto operator()(const Tree_loc& lhs, const Tree_loc& rhs) const { return lhs < rhs; }
};

template<typename X>
using Tree_annotation_map = absl::btree_map<Tree_loc, X, Tree_loc_cmp>;
template<typename X>
using Tree_annotation_multimap = absl::btree_multimap<Tree_loc, X, Tree_loc_cmp>;

template<typename X>
auto operator<<(std::ostream& os, const Tree_annotation_map<X>& atts) -> std::ostream& {
  return os << "{" << absl::StrJoin(atts, ", ", absl::StreamFormatter()) << "}";
}
template<typename X>
auto operator<<(std::ostream& os, const Tree_annotation_multimap<X>& atts) -> std::ostream& {
  return os << "{" << absl::StrJoin(atts, ", ", absl::StreamFormatter()) << "}";
}

template<typename As>
inline auto on_branch(Branch_index i, As& atts) {
  auto [first, last] = atts.equal_range(i);
  return std::ranges::subrange(first, last);
}
template<typename As>
inline auto on_branch(Branch_index i, const As& atts) {
  auto [first, last] = atts.equal_range(i);
  return std::ranges::subrange(first, last);
}


// Mutations are annotations at specific points on a tree where the sequence changes.
// Multiple sites can in principle change at the same time (though it's probably very rare); hence, mulimap.

struct Mutation_info {
  Site_index site;
  Real_seq_letter from;
  Real_seq_letter to;
  auto operator<=>(const Mutation_info& that) const = default;
};
inline auto operator<<(std::ostream& os, const Mutation_info& mm) -> std::ostream& {
  return os << absl::StreamFormat("%c%d%c", to_char(mm.from), mm.site, to_char(mm.to));
}

using Mutation = Tree_annotation<Mutation_info>;
using Mutations = Tree_annotation_multimap<Mutation_info>;

// Special method for erasing a (single) mutation at a particular location and site
// (we expect multiple mutations at the same location to be very rare, so we just do a linear search)
inline auto erase_mutation(Mutations& ms, const Tree_loc& loc, Site_index site) -> void {
  for (auto it = ms.find(loc); it != ms.end(); ++it) {
    const auto& [it_loc, it_mm] = *it;
    if (it_loc != loc) { return; }
    if (it_mm.site == site) {
      ms.erase(it);  // invalidates `it`, but ok because we're about to return
      return;
    }
  }
}
inline auto erase_mutation(Mutations& ms, const Mutation& m) -> void {
  const auto& [loc, mm] = m;
  erase_mutation(ms, loc, mm.site);
}

}  // namespace sapling

#endif // SAPLING_MUTATIONS_H_
