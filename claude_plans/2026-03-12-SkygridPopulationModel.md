# Skygrid Population Model Plan

## Context

Sapling currently supports constant and exponential population models. We want to support arbitrary population curves via the Skygrid model (Gill et al 2012), where `log(N(t))` is defined at a set of knots and interpolated between them (staircase or log-linear). This enables simulating coalescent trees under realistic, time-varying population dynamics.

The `Pop_model` interface requires `intensity_at_time(t)` and `inverse_intensity(I)`. For Skygrid, these have exact closed-form solutions within each interval — the same formulas already used for constant and exponential population models. The implementation walks through the Skygrid intervals, using precomputed cumulative values at interval boundaries for efficient lookup.

`coal_sim` works unchanged.

## `Pop_model` interface changes (in `pop_model.h/.cpp`)

Replace the virtual `sample(min_t, max_t, prng)` method with a virtual `inverse_cum_pop(P)` method. Sampling becomes model-agnostic and moves to `sapling.cpp`:

```cpp
// In choose_tip_times (sapling.cpp):
auto cum_pop_lo = pop_model.cum_pop_at_time(min_tip_t);
auto cum_pop_hi = pop_model.cum_pop_at_time(max_tip_t);
for (...) {
  auto u = absl::Uniform(absl::IntervalClosedOpen, rng, cum_pop_lo, cum_pop_hi);
  auto raw_t = pop_model.inverse_cum_pop(u);
  ...
}
```

This works for all population models: `inverse_cum_pop` is the exact inverse of `cum_pop_at_time`, so sampling proportional to `N(t)` is just uniform sampling in cum_pop space.

Update existing models:
- `Const_pop_model`: `inverse_cum_pop(P) = P / pop` (trivial inverse of `cum_pop_at_time(t) = t * pop`)
- `Exp_pop_model`: `inverse_cum_pop(P) = t0 + log1p(P * g / n0) / g` (inverse of `cum_pop_at_time(t) = (exp(g*(t-t0)) - 1) * n0/g`); `g == 0`: `P / n0`

Remove `distributions.h` entirely (no longer needed).

## `Skygrid_pop_model` class (in `pop_model.h/.cpp`)

Structure mirrors Delphy's `Skygrid_pop_model` in `pop_model.h`. Copy across Delphy's full explanatory comments verbatim (the ASCII diagram, interval numbering, `Interval(t)` definition, `log N(t)` formulas for both staircase and log-linear modes, and the relationship to Gill et al 2012 / BEAST notation).

Add `#include <vector>` to `pop_model.h`.

```cpp
class Skygrid_pop_model : public Pop_model {
 public:
  enum class Type { k_staircase, k_log_linear };

  Skygrid_pop_model(
      std::vector<double> x,      // x[k] = time of knot k
      std::vector<double> gamma,  // gamma[k] = log(N(x[k]))
      Type type);

  // Accessors (no setters, same pattern as other Pop_models)
  auto x() const -> const std::vector<double>&;
  auto gamma() const -> const std::vector<double>&;
  auto type() const -> Type;
  auto M() const -> int;  // = ssize(x) - 1

  // Derived convenience accessors matching Delphy's notation
  auto x(int k) const -> double;
  auto gamma(int k) const -> double;

  // Pop_model interface
  auto pop_at_time(double t) const -> double override;
  auto cum_pop_at_time(double t) const -> double override;
  auto intensity_at_time(double t) const -> double override;
  auto inverse_intensity(double I) const -> double override;
  auto inverse_cum_pop(double P) const -> double override;

  // Interval(t) in the top-level comment for Skygrid_pop_model
  auto interval_containing_t(double t) const -> int;

 private:
  std::vector<double> x_;
  std::vector<double> gamma_;
  Type type_;

  // Precomputed at knot boundaries (size = M+1, indices 0..M)
  // cum_intensity_[0] = 0, cum_intensity_[k] = integral of 1/N(t) from x_0 to x_k
  std::vector<double> cum_intensity_;
  // cum_pop_[0] = 0, cum_pop_[k] = integral of N(t) from x_0 to x_k
  std::vector<double> cum_pop_;

  auto print_to(std::ostream& os) const -> void override;
};
```

### Key method implementations

Within each Skygrid interval, `N(t)` is either constant (staircase) or exponential (log-linear).

**Interval `k` (`1 <= k <= M`)** runs from `x_{k-1}` to `x_k`, has width `w = x_k - x_{k-1}`, and has population size `N(x_k) = exp(gamma_k)` at its right endpoint.

```
                     |         |       |           |       |
                     |         |       |           |       |
        Int.0        |  Int.1  | Int.2 |    ...    | Int.M |   Int.{M+1}
                     |         |       |           |       |
                     |         |       |           |       |
   -INFINITY ...-----+---------+-------+----....---+-------+---... +INFINITY  --> t
                    x_0       x_1     x_2       x_{M-1}   x_M
```

**Staircase interval formulas:**
In staircase mode, `N(t) = exp(gamma_k)` (constant) throughout interval `k`. All integrals reduce to the constant-N case:
- Intensity from `x_{k-1}` to `t`: `(t - x_{k-1}) / exp(gamma_k)`
- Pop from `x_{k-1}` to `t`: `exp(gamma_k) * (t - x_{k-1})`
- Full interval intensity: `w / exp(gamma_k)`
- Full interval pop: `exp(gamma_k) * w`

Note: staircase uses `exp(gamma_k)` (the right endpoint's value), which differs from the log-linear `g == 0` case that uses `exp(gamma_{k-1})` when `gamma_{k-1} != gamma_k`.

**Log-linear interval parametrization:**
Following Delphy's `Exp_pop_model` convention, `N(t)` within interval `k` is parametrized using the interval's right endpoint:

```
N(t) = n0 * exp(g * (t - x_k))

where n0 = exp(gamma_k)                                 (population at interval end x_k)
      g  = (gamma_k - gamma_{k-1}) / (x_k - x_{k-1})    (local growth rate)
```

At `t = x_k`: `N(x_k) = exp(gamma_k)`. At `t = x_{k-1}`: `N(x_{k-1}) = exp(gamma_k - g*w) = exp(gamma_{k-1})`.

When `g == 0` (i.e., `gamma_{k-1} == gamma_k`), the formulas reduce to the constant-N case. This must be handled explicitly (same `if (g == 0.0)` pattern as `Exp_pop_model`).

**Key identity:** `n0 * exp(-g*w) = exp(gamma_k) * exp(-(gamma_k - gamma_{k-1})) = exp(gamma_{k-1})`. This simplification appears repeatedly in the inverse formulas below.

**Interval integrals (log-linear, following Delphy's `intensity_integral`/`pop_integral` style):**

For any sub-range `[a, b]` within interval `k`:

```
int_a^b 1/N(t) dt = (1/(g*n0)) * exp(-g*(b - x_k)) * expm1(g*(b - a))    [g != 0]
                  = (b - a) / n0                                         [g = 0]

int_a^b N(t) dt   = (n0/g) * exp(g*(a - x_k)) * expm1(g*(b - a))         [g != 0]
                  = n0 * (b - a)                                         [g = 0]
```

For the **full interval** (`a = x_{k-1}`, `b = x_k`, so `b - a = w`):
- intensity: `(1/(g*n0)) * expm1(g*w)` ; `g = 0`: `w / n0`
- pop: `-(n0/g) * expm1(-g*w)` ; `g = 0`: `n0 * w`
  (since `exp(g*(x_{k-1} - x_k)) * expm1(g*w) = exp(-g*w) * (exp(g*w) - 1) = 1 - exp(-g*w) = -expm1(-g*w)`)

For a **partial interval** from `x_{k-1}` to `t` (with `delta = t - x_{k-1}`):
- intensity: `(1/(g*n0)) * exp(-g*(t - x_k)) * expm1(g*delta)` ; `g = 0`: `delta / n0`
- pop: `(n0/g) * exp(-g*w) * expm1(g*delta) = (exp(gamma_{k-1})/g) * expm1(g*delta)` ; `g = 0`: `n0 * delta`
  (since `exp(g*(x_{k-1} - x_k)) = exp(-g*w)`, and by the key identity `n0 * exp(-g*w) = exp(gamma_{k-1})`)

**Inverse formulas** (solving for `t` given a target integral value from `x_{k-1}`):

For intensity inverse (given `remaining_I`, find `delta = t - x_{k-1}`):

Starting from the partial interval intensity formula:
```
remaining_I = (1/(g*n0)) * exp(-g*(t - x_k)) * expm1(g*delta)
```
Since `t - x_k = delta - w`: `exp(-g*(t - x_k)) = exp(-g*delta) * exp(g*w)`, so:
```
remaining_I = (exp(g*w) / (g*n0)) * exp(-g*delta) * (exp(g*delta) - 1)
            = (exp(g*w) / (g*n0)) * (1 - exp(-g*delta))
```
Using the key identity `exp(g*w) / n0 = exp(g*w) / exp(gamma_k) = 1 / exp(gamma_{k-1})`:
```
remaining_I = (1 - exp(-g*delta)) / (g * exp(gamma_{k-1}))
=> exp(-g*delta) = 1 - remaining_I * g * exp(gamma_{k-1})
=> delta = -log1p(-remaining_I * g * exp(gamma_{k-1})) / g    [g != 0]
   delta = remaining_I * exp(gamma_{k-1})                      [g = 0]
   t = x_{k-1} + delta
```

For cum_pop inverse (given `remaining_pop`, find `delta = t - x_{k-1}`):

Starting from the partial interval pop formula and applying the key identity `n0 * exp(-g*w) = exp(gamma_{k-1})`:
```
remaining_pop = (n0/g) * exp(-g*w) * expm1(g*delta)
              = (exp(gamma_{k-1})/g) * expm1(g*delta)
=> expm1(g*delta) = remaining_pop * g / exp(gamma_{k-1})
=> g*delta = log1p(remaining_pop * g / exp(gamma_{k-1}))
=> delta = log1p(remaining_pop * g / exp(gamma_{k-1})) / g    [g != 0]
   delta = remaining_pop / exp(gamma_{k-1})                    [g = 0]
   t = x_{k-1} + delta
```

---

**`pop_at_time(t)`**: Exact Skygrid formula. Find interval via `interval_containing_t`, apply staircase (constant) or log-linear (interpolation) formula. Same as Delphy's `log_N(t)` + `exp`.

**`intensity_at_time(t)`** (integral of `1/N` from reference point `x_0` to `t`):
- Find interval `k` containing `t`
- For intervals `1 <= k <= M`: look up `cum_intensity_[k-1]` for the contribution of all complete intervals before `k`, then add partial contribution from `x_{k-1}` to `t` using the intensity integral formula above
- For interval 0 (`t <= x_0`): `(t - x_0) / exp(gamma_0)` (negative or zero, constant `N = exp(gamma_0)`)
- For interval `M+1` (`t > x_M`): `cum_intensity_[M] + (t - x_M) / exp(gamma_M)` (constant `N = exp(gamma_M)`)

**`inverse_intensity(I)`**:
- `I <= 0`: in region at or below `x_0`, constant N. `t = x_0 + I * exp(gamma_0)`. Returns `x_0` when `I == 0`.
- `I >= cum_intensity_[M]`: in region at or above `x_M`, constant N. `t = x_M + (I - cum_intensity_[M]) * exp(gamma_M)`. Returns `x_M` when `I == cum_intensity_[M]`.
- Otherwise (`0 < I < cum_intensity_[M]`): binary search (`std::upper_bound`) in `cum_intensity_` to find the smallest `k` such that `cum_intensity_[k] > I`. Then the target is in interval `k`, between `x_{k-1}` and `x_k`:
  - `remaining_I = I - cum_intensity_[k-1]`
  - Staircase: `t = x_{k-1} + remaining_I * exp(gamma_k)`
  - Log-linear: `t = x_{k-1} - log1p(-remaining_I * g * exp(gamma_{k-1})) / g`; falls back to `t = x_{k-1} + remaining_I * exp(gamma_{k-1})` when `g == 0`

**`cum_pop_at_time(t)`**: Analogous to `intensity_at_time` but integrating `N(t)` instead of `1/N(t)`. Uses `cum_pop_` precomputed array and the pop integral formulas above.
- For interval 0: `exp(gamma_0) * (t - x_0)` (negative or zero)
- For interval `M+1`: `cum_pop_[M] + exp(gamma_M) * (t - x_M)`

**`inverse_cum_pop(P)`**: Find `t` such that `cum_pop_at_time(t) = P`. Exact inverse of `cum_pop_at_time`, paralleling `inverse_intensity`:
- `P <= 0`: in region at or below `x_0`. `t = x_0 + P / exp(gamma_0)`
- `P >= cum_pop_[M]`: in region at or above `x_M`. `t = x_M + (P - cum_pop_[M]) / exp(gamma_M)`
- Otherwise: binary search in `cum_pop_` to find interval `k`. Invert locally using cum_pop inverse formula:
  - `remaining = P - cum_pop_[k-1]`
  - Staircase: `t = x_{k-1} + remaining / exp(gamma_k)`
  - Log-linear: `t = x_{k-1} + log1p(remaining * g / exp(gamma_{k-1})) / g`; falls back to `t = x_{k-1} + remaining / exp(gamma_{k-1})` when `g == 0`

### Construction

Validate inputs (matching Delphy's checks):
- At least 2 knots
- `ssize(x) == ssize(gamma)`
- Strictly increasing knot times

Precompute `cum_intensity_` and `cum_pop_` arrays (both size `M+1`) by accumulating exact per-interval contributions across all `M` interior intervals (`k = 1..M`), using the full-interval integral formulas above.

## Sampling (in `sapling.cpp`)

`choose_tip_times` uses the model-agnostic `cum_pop_at_time` / `inverse_cum_pop` pair for inverse-CDF sampling:

```cpp
auto cum_pop_lo = pop_model.cum_pop_at_time(min_tip_t);
auto cum_pop_hi = pop_model.cum_pop_at_time(max_tip_t);
for (...) {
  auto u = absl::Uniform(absl::IntervalClosedOpen, rng, cum_pop_lo, cum_pop_hi);
  auto raw_t = pop_model.inverse_cum_pop(u);
  ...
}
```

This replaces the current per-model `sample` method. `distributions.h` is no longer needed and is deleted.

## CLI options (in `sapling.cpp`)

New option group:
```
Skygrid population model [log N(t) specified at evenly-spaced knots]:
  --skygrid-first-knot-date First knot date x_0 (YYYY-MM-DD)
  --skygrid-last-knot-date  Last knot date x_M (YYYY-MM-DD)
  --skygrid-gamma           Comma-separated log population sizes gamma_k = ln(N(x_k)) at knots,
                            where N(t) is measured in years.
                            Mutually exclusive with --skygrid-Ns.
  --skygrid-Ns              Comma-separated population sizes N(x_k) at knots, measured in years.
                            Converted to gamma_k = ln(N_k) internally.
                            Mutually exclusive with --skygrid-gamma.
  --skygrid-type            "staircase" (default) or "log-linear"
```

The number of parameters (`M+1`) is inferred from the count of comma-separated values in `--skygrid-gamma` or `--skygrid-Ns`. Knot times are evenly spaced: `x_k = x_0 + k * (x_M - x_0) / M` for `k = 0, ..., M`.

Exactly one of `--skygrid-gamma` or `--skygrid-Ns` must be provided. Mutually exclusive with `--const-pop-*` and `--exp-pop-*`.

Population model mutual exclusivity is enforced by detecting the presence of *any* option from each model group (constant: `--const-pop-n0`; exponential: `--exp-pop-n0`, `--exp-pop-g`; Skygrid: `--skygrid-first-knot-date`, `--skygrid-last-knot-date`, `--skygrid-gamma`, `--skygrid-Ns`, `--skygrid-type`) and requiring that exactly one group has options present.

## Info JSON output (in `sapling.cpp`)

Add Skygrid case to `dump_info`:
```json
{
  "type": "skygrid",
  "skygrid_type": "staircase",
  "x_k": [-1.0, -0.5, 0.0],
  "x_k_dates": ["2019-01-01", "2019-07-02", "2020-01-01"],
  "gamma_k": [3.0, 4.0, 3.5],
  "N_k": [20.09, 54.60, 33.12]
}
```

- `x_k`: knot times in years relative to t0
- `x_k_dates`: knot times as ISO dates
- `gamma_k`: log population sizes at knots (`ln(N)`, with N in years)
- `N_k`: population sizes at knots in years (`= exp(gamma_k)`)

## Tests (in `tests/pop_model_tests.cpp`)

1. **Construction validation**: bad inputs (mismatched sizes, non-increasing knots, equal knots, < 2 knots)
2. **Exact `pop_at_time`**: staircase and log-linear against hand-computed values, including at knot boundaries, between knots, and outside the knot range
3. **Intensity roundtrip**: `inverse_intensity(intensity_at_time(t)) == t` (within tolerance) for many `t` values: inside intervals, outside range, at knot boundaries — for both staircase and log-linear
4. **Cum_pop roundtrip**: `inverse_cum_pop(cum_pop_at_time(t)) == t` (within tolerance) for the same sweep of `t` values — for both staircase and log-linear
5. **Consistency with simpler models**: uniform-gamma Skygrid should match `Const_pop_model`; single-interval log-linear should match `Exp_pop_model`
6. **`cum_pop_at_time` numerical derivative**: verify `cum_pop_at_time(t + dt) - cum_pop_at_time(t) == dt * pop_at_time(t)` within `O(dt)` approximation error, for many `t` values across all interval types (staircase/log-linear, interior/exterior) — avoiding knot boundaries where the derivative is discontinuous (staircase)
7. **`intensity_at_time` numerical derivative**: verify `intensity_at_time(t + dt) - intensity_at_time(t) == dt / pop_at_time(t)` within `O(dt)` approximation error, for the same sweep of `t` values
8. **`inverse_cum_pop` for existing models**: verify `inverse_cum_pop(cum_pop_at_time(t)) == t` roundtrip for `Const_pop_model` and `Exp_pop_model` (including zero-growth `Exp_pop_model`)
9. **Integration with `coal_sim`**: Skygrid model produces valid trees (`assert_phylo_tree_integrity`) — for both staircase and log-linear
10. **Print**: verify debug output format contains expected substrings for both staircase and log-linear

## Files to modify

| File | Changes |
|------|---------|
| `core/pop_model.h` | Replace `sample` with `inverse_cum_pop` in `Pop_model` base class; add `#include <vector>`; add `Skygrid_pop_model` class with full Delphy-style comments |
| `core/pop_model.cpp` | Replace `sample` with `inverse_cum_pop` for `Const_pop_model` and `Exp_pop_model`; implement `Skygrid_pop_model` methods |
| `core/distributions.h` | Delete (no longer needed) |
| `core/CMakeLists.txt` | Remove `distributions.h` from source list |
| `sapling.cpp` | Sampling via `inverse_cum_pop` in `choose_tip_times`; Skygrid CLI options; robust mutual exclusivity check across all three model groups; `dump_info` Skygrid case |
| `tests/pop_model_tests.cpp` | Skygrid tests; `inverse_cum_pop` roundtrip tests for existing models |

## Verification

1. Build: `cmake --build --preset conan-debug`
2. Run tests: `./build/debug/tests/tests --gtest_filter="Pop_model*"`
3. End-to-end:
```bash
./build/debug/sapling \
  --skygrid-first-knot-date=2019-01-01 \
  --skygrid-last-knot-date=2020-01-01 \
  --skygrid-gamma="3.0,4.0,3.5" \
  --skygrid-type=staircase \
  -n 100 --min-tip-t 2019-01-01 --max-tip-t 2020-01-01 \
  -L 1000 --out-prefix test_skygrid
```
