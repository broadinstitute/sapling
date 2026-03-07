#ifndef SAPLING_ALIAS_SAMPLER_H_
#define SAPLING_ALIAS_SAMPLER_H_

// Alias_sampler: O(1) weighted random sampling using Vose's alias method.
//
// References:
// - Walker, A.J. (1977). "An Efficient Method for Generating Discrete Random Variables
//   with General Distributions". ACM TOMS 3(3):253-256.
//   https://doi.org/10.1145/355744.355749
// - Vose, M.D. (1991). "A Linear Algorithm for Generating Random Numbers with a Given
//   Distribution". IEEE Trans. Software Eng. 17(9):972-975.
//   https://doi.org/10.1109/32.92917
// - Knuth, D.E. "The Art of Computer Programming", Vol. 2: Seminumerical Algorithms,
//   3rd ed., Section 3.4.1 A, pp. 120-121.
// - Keith Schwarz's tutorial: https://www.keithschwarz.com/darts-dice-coins/

#include <numeric>
#include <vector>

#include "absl/log/check.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/distributions.h"

namespace sapling {

struct Alias_sampler {
  std::vector<double> prob;
  std::vector<int> alias;

  // Build alias table from weights (need not be normalized).
  explicit Alias_sampler(const std::vector<double>& weights) {
    auto n = static_cast<int>(std::ssize(weights));
    CHECK_GT(n, 0);

    auto total = std::accumulate(weights.begin(), weights.end(), 0.0);
    CHECK_GT(total, 0.0);

    // Initialize prob with normalized weights (each scaled so they sum to n)
    prob.resize(n);
    alias.resize(n);
    auto scale = static_cast<double>(n) / total;
    for (auto i = 0; i < n; ++i) {
      prob[i] = weights[i] * scale;
    }

    // Partition indices into small (prob < 1) and large (prob >= 1).
    // Note: it is possible to eliminate these stacks using two scanning pointers
    // and a temporary variable (see Vose 1991, Section VI).  The stackless
    // algorithm is fully spelled out in claude_plans/2026-03-06-SiteRateHeterogeneity.md.
    // We prefer the explicit stacks: they are clearer and the O(n) allocation is negligible.
    auto small = std::vector<int>{};
    auto large = std::vector<int>{};
    for (auto i = 0; i < n; ++i) {
      if (prob[i] < 1.0) {
        small.push_back(i);
      } else {
        large.push_back(i);
      }
    }

    // Pair up small and large entries
    while (not small.empty() and not large.empty()) {
      auto s = small.back();
      small.pop_back();
      auto l = large.back();
      large.pop_back();

      alias[s] = l;
      // Transfer part of l's probability to fill up s's bucket.
      // Numerically stable version of prob[l] -= 1.0 - prob[s].
      // After this, prob[s] has its final alias-method meaning, while prob[l]
      // is still working space.  Since prob[l] >= 1.0 and prob[s] < 1.0, the
      // result is >= 0.0, but l's bucket may go from "large" to "small".
      prob[l] = (prob[l] + prob[s]) - 1.0;

      if (prob[l] < 1.0) {
        small.push_back(l);
      } else {
        large.push_back(l);
      }
    }

    // Remaining entries get probability 1.0 (numerical cleanup)
    for (auto i : large) {
      prob[i] = 1.0;
    }
    for (auto i : small) {
      prob[i] = 1.0;
    }
  }

  // Sample an index in [0, n) with probability proportional to weights[i].
  auto sample(absl::BitGenRef rng) const -> int {
    auto n = static_cast<int>(std::ssize(prob));
    auto i = absl::Uniform(absl::IntervalClosedOpen, rng, 0, n);
    if (prob[i] == 1.0) { return i; }
    auto u = absl::Uniform(absl::IntervalClosedOpen, rng, 0.0, 1.0);
    return (u < prob[i]) ? i : alias[i];
  }

  auto size() const -> int { return static_cast<int>(std::ssize(prob)); }
};

}  // namespace sapling

#endif // SAPLING_ALIAS_SAMPLER_H_
