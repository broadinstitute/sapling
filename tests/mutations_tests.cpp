#include "mutations.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace sapling {

inline constexpr auto rA = Real_seq_letter::A;
inline constexpr auto rC = Real_seq_letter::C;
inline constexpr auto rG = Real_seq_letter::G;
inline constexpr auto rT = Real_seq_letter::T;

TEST(Mutations_test, print) {
  const auto m = Mutation{{4, 3.25}, {123, rA, rG}};

  auto ss = std::stringstream{};
  ss << m;

  EXPECT_THAT(ss.str(), testing::StrEq("[[4,3.25]]:A123G"));
}

TEST(Mutations_test, mutations_lookup) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},  // Purposely inserted unordered to test reordering below
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  auto ss = std::stringstream{};
  ss << ms;

  EXPECT_THAT(ss.str(), testing::StrEq(
      "{"
      "[[1,0.25]]:T123A, "
      "[[1,0.5]]:C345T, "
      "[[1,0.75]]:A456C, "
      "[[4,3.25]]:A123G"
      "}"));
}

template<std::ranges::range R>
static auto to_vec(R&& range) {
  auto result = std::vector<std::ranges::range_value_t<R>>{};
  for (const auto& elem : range) {
    result.push_back(elem);
  }
  return result;
}

TEST(Mutations_test, on_branch) {
  // ms is `const`; checking that on_branch can take a const Mutations
  const auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},  // Purposely inserted unordered to test reordering below
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  EXPECT_THAT(to_vec(on_branch(1, ms)), testing::ElementsAre(
      Mutation{{1, 0.25}, {123, rT, rA}},
      Mutation{{1, 0.50}, {345, rC, rT}},
      Mutation{{1, 0.75}, {456, rA, rC}}));
}

TEST(Mutations_test, on_branch_mut) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},  // Purposely inserted unordered to test reordering below
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  for (auto& [loc, mm] : on_branch(1, ms)) {
    mm.site += 9000;  // Check we can change underlying mutations
  }

  EXPECT_THAT(to_vec(on_branch(1, ms)), testing::ElementsAre(
      Mutation{{1, 0.25}, {9123, rT, rA}},
      Mutation{{1, 0.50}, {9345, rC, rT}},
      Mutation{{1, 0.75}, {9456, rA, rC}}));
}

TEST(Mutations_test, heterogeneous_lookup_loc) {
  const auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},  // Purposely inserted unordered to test reordering below
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  auto first = ms.lower_bound(Tree_loc{1, 0.45});
  auto last = ms.upper_bound(Tree_loc{1, 0.80});
  EXPECT_THAT(to_vec(std::ranges::subrange(first, last)), testing::ElementsAre(
      Mutation{{1, 0.50}, {345, rC, rT}},
      Mutation{{1, 0.75}, {456, rA, rC}}));
}

TEST(Mutations_test, heterogeneous_lookup_branch_erase) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},  // Purposely inserted unordered to test reordering below
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  ms.erase(1);

  EXPECT_THAT(ms, testing::ElementsAre(
      Mutation{{4, 3.25}, {123, rA, rG}}));
}

TEST(Mutations_test, heterogeneous_lookup_loc_erase) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},  // Purposely inserted unordered to test reordering below
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  ms.erase(Tree_loc{1, 0.50});

  EXPECT_THAT(ms, testing::ElementsAre(
      Mutation{{1, 0.25}, {123, rT, rA}},
      Mutation{{1, 0.75}, {456, rA, rC}},
      Mutation{{4, 3.25}, {123, rA, rG}}));
}

TEST(Mutations_test, erase_mutation_missing) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},
    {{1, 0.75}, {456, rA, rC}},
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  erase_mutation(ms, Tree_loc{2, 6.50}, 123);  // no mutation at this location
  erase_mutation(ms, Tree_loc{1, 0.50}, 456);  // no mutation at this location with that site

  EXPECT_THAT(std::ssize(ms), testing::Eq(4));
}

TEST(Mutations_test, erase_mutation_single) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},  // Delete this one
    {{1, 0.75}, {456, rA, rC}},
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  erase_mutation(ms, Tree_loc{1, 0.5}, 345);

  EXPECT_THAT(ms, testing::ElementsAre(
      Mutation{{1, 0.25}, {123, rT, rA}},
      Mutation{{1, 0.75}, {456, rA, rC}},
      Mutation{{4, 3.25}, {123, rA, rG}}));
}

TEST(Mutations_test, erase_mutation_multiple) {
  auto ms = Mutations{
    {{1, 0.50}, {123, rC, rT}},
    {{1, 0.50}, {456, rT, rA}},  // Delete this one
    {{1, 0.50}, {789, rA, rG}},
    {{1, 0.75}, {456, rA, rC}},
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  erase_mutation(ms, Tree_loc{1, 0.5}, 456);

  EXPECT_THAT(ms, testing::ElementsAre(
      Mutation{{1, 0.25}, {123, rT, rA}},
      Mutation{{1, 0.50}, {123, rC, rT}},
      Mutation{{1, 0.50}, {789, rA, rG}},
      Mutation{{1, 0.75}, {456, rA, rC}},
      Mutation{{4, 3.25}, {123, rA, rG}}));
}

TEST(Mutations_test, erase_mutation_convenience) {
  auto ms = Mutations{
    {{1, 0.50}, {345, rC, rT}},  // Delete this one
    {{1, 0.75}, {456, rA, rC}},
    {{1, 0.25}, {123, rT, rA}},
    {{4, 3.25}, {123, rA, rG}}
  };

  auto m = Mutation{{1, 0.50}, {345, rC, rT}};
  
  erase_mutation(ms, m);

  EXPECT_THAT(ms, testing::ElementsAre(
      Mutation{{1, 0.25}, {123, rT, rA}},
      Mutation{{1, 0.75}, {456, rA, rC}},
      Mutation{{4, 3.25}, {123, rA, rG}}));
}

}  // namespace sapling
