# Site-Rate Heterogeneity

## Overview

Add a `--site-rate-heterogeneity-alpha` CLI parameter to enable Gamma-distributed site-rate
heterogeneity.  When provided with a positive value, each site `l` gets a rate modifier
`nu_l ~ Gamma(shape=alpha, rate=alpha)` (mean 1, variance 1/alpha).  A default value of
`0.0` means "no site-rate heterogeneity" (all `nu_l = 1.0`).  Site selection during the
Gillespie mutation simulation switches from a linear scan to the alias algorithm for O(1)
weighted sampling, combined with thinning to handle base-dependent escape rates.

## Current State

- `Global_evo_model` already has `nu_l` (`Site_vector<double>`), initialized to all 1.0 by
  `make_single_partition_global_evo_model()` in `core/evo_model.cpp`.
- `simulate_mutations()` in `sapling.cpp` already uses `evo.Q_l_a(l, seq[l])` which
  multiplies by `nu_l[l]`, so per-site rates are wired through the rate calculations.
- Site selection currently uses a linear O(L) scan over cumulative `Q_l_a` values.
- There is a TODO at `sapling.cpp:743`: `// TODO: Fill in evo.nu_l[l] with something other than 1.0`.

## Design

### 1. CLI Option

Add `--site-rate-heterogeneity-alpha` to the "Substitution model" option group in
`process_args()`.

- Type: `double`, default `0.0`.
- Value `0.0`: no site-rate heterogeneity (all `nu_l = 1.0`, current behavior).
- Positive value: enables Gamma-distributed site rates.
- Store as `double site_rate_heterogeneity_alpha` in `Options`.

### 2. Sampling `nu_l`

In `main()`, after `make_single_partition_global_evo_model()`:

- If `alpha > 0.0`, for each site `l` in `[0, L)`:
  - `nu_l[l] ~ Gamma(shape=alpha, rate=alpha)`, i.e., `Gamma(shape=alpha, scale=1/alpha)`.
  - Note: we use the (shape, rate) parameterization convention; `std::gamma_distribution`
    takes (shape, scale), so we use `std::gamma_distribution<double>{alpha, 1.0/alpha}`.
    (`absl::Gamma` does not exist.)
  - This gives mean 1 and variance 1/alpha.
- If `alpha == 0.0`, explicitly fill `nu_l` with 1.0.

### 3. Alias Table for Site Selection

**New file: `core/alias_sampler.h`**

Implement a reusable `Alias_sampler` class using Vose's alias method.  Builds an O(n) table
from a vector of non-negative weights and provides O(1) weighted random sampling.

References:
- Walker, A.J. (1977). "An Efficient Method for Generating Discrete Random Variables with
  General Distributions". ACM TOMS 3(3):253-256.
  https://doi.org/10.1145/355744.355749
- Vose, M.D. (1991). "A Linear Algorithm for Generating Random Numbers with a Given
  Distribution". IEEE Trans. Software Eng. 17(9):972-975.
  https://doi.org/10.1109/32.92917
- Knuth, D.E. "The Art of Computer Programming", Vol. 2: Seminumerical Algorithms,
  3rd ed., Section 3.4.1 A, pp. 120-121.
- Keith Schwarz's tutorial: https://www.keithschwarz.com/darts-dice-coins/

**Why not use an existing implementation?**

- `std::discrete_distribution` does not guarantee O(1) sampling — the C++ standard leaves
  the algorithm unspecified.  In practice, libstdc++ and libc++ both use O(log n) binary
  search on cumulative weights.
- `boost::random::discrete_distribution` also uses O(log n) binary search internally,
  not the alias method.
- Abseil does not provide a weighted discrete sampler.
- Standalone alias-method libraries exist on GitHub (e.g., `AliasMethod` by various
  authors), but they are not well-maintained or widely adopted, and adding a dependency
  for ~40 lines of straightforward code is not worthwhile.

```
struct Alias_sampler {
  // Build alias table from weights (need not be normalized)
  explicit Alias_sampler(const std::vector<double>& weights);

  // Sample an index in [0, n) with probability proportional to weights[i]
  auto sample(absl::BitGenRef rng) const -> int;

  // Number of entries
  auto size() const -> int;

  // Internals
  std::vector<double> prob;   // probability thresholds
  std::vector<int> alias;     // alias indices
};
```

The `sample()` method skips the second uniform draw when `prob[i] == 1.0` (i.e., the
bucket has no alias).

The alias table is built once from the `nu_l` values and remains static for the entire
simulation (since `nu_l` doesn't change).

### 4. Changes to `simulate_mutations()`

The core change: replace the O(L) linear scan for site selection with the alias sampler,
combined with thinning to handle base-dependent escape rates.  The `Alias_sampler` is
constructed inside `simulate_mutations()` from `evo.nu_l`, so the function signature
does not change.

**Current approach** (exact, O(L) per mutation):
- Total rate: `lambda = sum_l mu * nu_l * q_a(seq[l])`
- Site selection: linear scan over cumulative `Q_l_a` values
- Updates `lambda` incrementally when mutations occur

**New approach** (thinning + alias, O(1) per mutation attempt):

The per-site rate `Q_l_a(l, seq[l]) = mu * nu_l * q_a(seq[l])` depends on the current base
at site `l`, not just `nu_l`.  A static alias table built from `nu_l` alone cannot account
for the base-dependent escape rate `q_a`.  We handle this with thinning:

1. Precompute `q_a_max = max_a q_a(a)` (the maximum escape rate over all 4 bases).
2. Upper-bound total rate: `lambda_upper = mu * q_a_max * sum_l nu_l`.
   This is constant throughout the simulation (independent of sequence state).
3. Draw next event time: `t += Exponential(lambda_upper)`.
4. Pick site `l` from the alias table (weight proportional to `nu_l`), O(1).
5. Accept mutation with probability `q_a(seq[l]) / q_a_max`.
   - If rejected, discard this event and draw the next one (thinning).
   - If accepted, pick target base `b` from off-diagonal `q_ab[a][b]` row and apply.
6. `lambda_upper` never needs updating (it's independent of sequence state).

This is a valid thinning of the inhomogeneous Poisson process.  The expected acceptance
rate is `E[q_a(seq[l])] / q_a_max`, which is close to 1 when `q_a` doesn't vary much
across bases (exactly 1 when kappa=1 and equal base frequencies, i.e., all bases have
the same escape rate).

This approach replaces exact rate tracking with thinning.  The statistical output is
identical (thinning produces the exact same distribution of mutations), but the sequence
of random draws differs, so results won't be reproducible relative to the old code at
the same seed even when alpha is 0.0.  This is acceptable.

An alternative would be exact rate tracking using a Fenwick tree for O(log L) dynamic
weighted sampling, but the O(1) alias + thinning approach is simpler and more efficient.

Remove the existing TODO comment about Fenwick trees in `simulate_mutations()`.  Add a
detailed explanatory comment in the code describing the thinning + alias algorithm and
the reasoning behind it.

### 5. Info JSON Output

In `dump_info()`, add to the `"subst_model"` object (only when alpha > 0.0):

- `"site_rate_heterogeneity_alpha"`: the alpha value
- `"nu_l"`: array of all L values of `nu_l[l]`

When alpha is 0.0, these fields are omitted.

The `nu_l` values live in `Global_evo_model`, so `dump_info()` needs an additional
parameter.  Change its signature to also accept `const Global_evo_model& evo`.

Additionally, add an `--output-mutation-counts` boolean CLI flag.  When set, the info
JSON output includes a `"mutation_counts"` array (parallel to `"nu_l"`) containing the
total number of mutations observed at each site across the whole tree.  This is useful
for validating that the site-rate heterogeneity sampling is correct.

### 6. Unit Tests

**New file: `tests/alias_sampler_tests.cpp`**

- Test edge case: single element (always returns 0).
- Test that a single non-zero weight among zeros always returns that index.
- Test that uniform weights produce roughly uniform samples (3-sigma binomial tolerance).
- Test that known weights (1:2:3:4) produce the correct distribution (3-sigma binomial
  tolerance).
- Test extreme weight ratios (1000:1:1:1) produce the correct distribution.

No changes to existing tests needed.

## Files Changed

| File | Change |
|------|--------|
| `core/alias_sampler.h` | **New.** Alias_sampler class (header-only) |
| `tests/alias_sampler_tests.cpp` | **New.** Unit tests for Alias_sampler |
| `tests/CMakeLists.txt` | Add `alias_sampler_tests.cpp` |
| `sapling.cpp` | Add `--site-rate-heterogeneity-alpha` CLI option, `--output-mutation-counts` flag, sample `nu_l` from Gamma, build alias table, rewrite site selection in `simulate_mutations()` to use thinning + alias, add alpha/nu_l/mutation_counts to info JSON output |

## Implementation Order

1. Implement `Alias_sampler` in `core/alias_sampler.h` + tests
2. Add `--site-rate-heterogeneity-alpha` CLI option and Gamma sampling of `nu_l`
3. Rewrite `simulate_mutations()` to use thinning + alias table
4. Add alpha and nu_l to info JSON output
5. Build and run tests
6. Manual integration test with a sample run

## Endnote: Stackless Vose Algorithm

Vose (1991, Section VI) hints that the `small` and `large` stacks can be eliminated using
two forward-scanning pointers `j` (for small entries) and `k` (for large entries), plus a
single temporary variable.  Here is the fully spelled-out algorithm:

```
j = 0       // scans forward for small entries (prob < 1.0)
k = 0       // scans forward for large entries (prob >= 1.0)
temp = -1   // deferred small entry (at most one)

loop:
    // 1. Get next small entry
    if temp >= 0:
        s = temp
        temp = -1
    else:
        advance j until prob[j] < 1.0 (or j >= n → done)
        s = j; j++

    // 2. Get next large entry
    advance k until prob[k] >= 1.0 (or k >= n → done)
    l = k

    // 3. Pair them (same as standard Vose)
    alias[s] = l
    prob[l] = (prob[l] + prob[s]) - 1.0

    // 4. If l became small, handle it
    if prob[l] < 1.0:
        if l < j:     // j already passed l → must defer it
            temp = l
        // else:      // j hasn't reached l → j will find it naturally
```

**Why this works:** `j` and `k` scan the same `prob` array, which is mutated in place.
When a large entry at position `l` gets drained below 1.0, there are two cases:

1. **`l >= j`**: `j` hasn't reached `l` yet, so it will naturally find it as a small entry
   during a future scan.  No special handling needed.

2. **`l < j`**: `j` has already passed `l` and will never see it again.  We save `l` in
   `temp` so it is used as the small entry in the next iteration, *before* resuming the `j`
   scan.

There is never a conflict with `temp` already being set, because `temp` is always consumed
at the *start* of the very next iteration (before any new pairing could set it again).

**Why `k` and `j` don't interfere:**

- `k` skips entries with `prob < 1.0`, so it never picks up a small entry (including entries
  recently used as `s`).
- `j` skips entries with `prob >= 1.0`, so it never picks up a large entry.
- After pairing `(s, l)`, if `prob[l]` drops but stays `>= 1.0`, `k` stays at position `l`
  and reuses it next iteration (it still has excess probability to donate).  `k` only advances
  past `l` when `prob[l]` finally drops below 1.0.
- `j` never re-finds `s` because `j` was already incremented past `s`, and `prob[s]` is
  unchanged (only `prob[l]` is modified).

**Worked example:** Weights `[5, 5, 1, 1]`, normalized to `prob = [1.667, 1.667, 0.333, 0.333]`:

| Iter | Source              | s | l (k)                | alias[s]=l | new prob[l] | temp              |
|------|---------------------|---|----------------------|------------|-------------|-------------------|
| 1    | j→2 (skip 0,1)      | 2 | 0                    | alias[2]=0 | 1.0         | —                 |
| 2    | j→3                 | 3 | 0                    | alias[3]=0 | 0.333       | 0 (since 0 < j=4) |
| 3    | temp=0              | 0 | 1 (k advances past 0)| alias[0]=1 | 1.0         | —                 |

Result: `prob = [0.333, 1.0, 0.333, 0.333]`, `alias = [1, —, 0, 0]`.  Correct.

**Why we don't use it:** This eliminates two `vector` allocations but adds a conditional
branch and the `temp` variable.  The logic is trickier to reason about.  For ~40 lines of
code running once at startup on a table of size L, the explicit stacks are clearer and the
O(L) allocation is negligible.
