# Tip-Date Uncertainty

## Context

Real sequencing metadata is often incomplete: some samples have only a month-level date (YYYY-MM),
others only a year-level date (YYYY).  To make Sapling's simulated data more realistic, we add the
ability to simulate tip-date uncertainty.  For each tip, the full date is independently masked to
month precision or year precision with configurable probabilities.

## Overview

For each tip, independently draw from a discrete distribution:
- With probability `p_month`: mask the day, outputting `YYYY-MM` instead of `YYYY-MM-DD`.
- With probability `p_year`: mask the day and month, outputting `YYYY` instead of `YYYY-MM-DD`.
- With probability `1 - p_month - p_year`: keep the full date `YYYY-MM-DD`.

When tip-date uncertainty is active (i.e., `p_month > 0` or `p_year > 0`):
- **Complete** FASTA/MAPLE files (with `-COMPLETE` suffix) preserve full dates in tip names.
- **Normal** FASTA/MAPLE files use uncertain (masked) dates in tip names.
- Newick/Nexus files use uncertain dates in tip names (no `-COMPLETE` versions).
- The info JSON includes a `tip_date_uncertainty` section.

When combined with missing data simulation, a single pair of `-COMPLETE` files is produced with
both full sequences and full dates.  The normal files have both masked sequences and uncertain dates.

When all uncertainty parameters are at their defaults (0.0), no uncertainty is simulated,
no `-COMPLETE` files are produced (unless missing data is active), and the info JSON has no
`tip_date_uncertainty` section.  The output is identical to the current behavior.

## CLI Options

Add a new option group `"Tip-date uncertainty"` with:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--p-tip-date-uncertain-upto-month` | `double` | `0.0` | Probability of masking tip date to month precision (YYYY-MM) |
| `--p-tip-date-uncertain-upto-year` | `double` | `0.0` | Probability of masking tip date to year precision (YYYY) |

Store in `Options` as two `double` fields: `p_tip_date_uncertain_upto_month` and
`p_tip_date_uncertain_upto_year`.

Validate:
- Both values are non-negative.
- Their sum is at most 1.0.

Tip-date uncertainty is considered **active** when `p_tip_date_uncertain_upto_month > 0.0` or
`p_tip_date_uncertain_upto_year > 0.0`.

## Data Structures

```cpp
// Per-tip date uncertainty level
enum class Date_uncertainty { none, month, year };

// All tip-date uncertainty, keyed by tip Node_index
using Tip_date_uncertainty = absl::flat_hash_map<Node_index, Date_uncertainty>;
```

These go directly in `sapling.cpp`.

## Algorithm: Simulating Tip-Date Uncertainty

```
Input: tree, p_month, p_year, rng
Output: Tip_date_uncertainty

1. For each tip node in the tree:
   a. Draw u ~ Uniform(0, 1).
   b. If u < p_month: assign Date_uncertainty::month.
      Else if u < p_month + p_year: assign Date_uncertainty::year.
      Else: assign Date_uncertainty::none.
   c. Store in the map.

Every tip gets an entry in the map when tip-date uncertainty is active,
even if Date_uncertainty::none was drawn.
```

## Output Strategy: Temporary Name Replacement

The tip dates are embedded in node names (e.g., `TIP_1|2020-06-15`).  All output functions
(FASTA, MAPLE, Newick, Nexus) read dates from these names.  To avoid threading uncertainty
through every output function, we temporarily swap node names when outputting files that need
uncertain dates.

```
apply_tip_date_uncertainty_to_names(tree, tip_date_uncertainty)
    -> flat_hash_map<Node_index, string>:
  Save each tip's current name, then for each tip with non-none uncertainty,
  truncate the date portion in tree.at(node).name.
  Return the saved names map.

restore_tip_names(tree, saved_names):
  For each entry in saved_names, restore tree.at(node).name.
```

The tree always holds full-precision names.  We swap to uncertain names only around specific
output calls, then restore immediately after.

**Helper function** to truncate the date in a tip name:

```cpp
auto truncate_date_in_name(const std::string& name, Date_uncertainty uncertainty) -> std::string;
```

Since `to_iso_date` always produces `YYYY-MM-DD` format, this simply truncates the date string
after the last `|`:
- `Date_uncertainty::month`: remove the last 3 characters (`-DD`), giving `YYYY-MM`.
- `Date_uncertainty::year`: remove the last 6 characters (`-MM-DD`), giving `YYYY`.
- `Date_uncertainty::none`: return unchanged.

Add explicit `CHECK`s that the input name conforms to assumptions:
- Name contains `|`.
- The date portion after the last `|` has at least 10 characters (`YYYY-MM-DD`).

Works for both auto-generated names (`TIP_1|2020-06-15`) and user-supplied names from tip files
(`SomeSample|2020-06-15`), since the tip file parser enforces the `Name|YYYY-MM-DD` format.

## Output Changes

### Complete Files

`-COMPLETE` files are now produced when **either** missing data or tip-date uncertainty is active.
Introduce a combined predicate:

```cpp
auto is_masking_active(const Options& opts) -> bool {
  return is_missing_data_active(opts) || is_tip_date_uncertainty_active(opts);
}
```

Replace uses of `missing_active` in `main()` with `masking_active` where it controls `-COMPLETE`
file production.  Missing data masking of sequences still only happens when missing data is
specifically active.

### Newick/Nexus Output

No changes to `output_newick_tree()`.  The Newick and Nexus outputs use uncertain dates because
the names are temporarily swapped before the output call.  No `-COMPLETE` Newick/Nexus files are
produced.

### Info JSON

When tip-date uncertainty is active, add a `"tip_date_uncertainty"` section to the info JSON:

```json
{
  "tip_date_uncertainty": {
    "p_uncertain_upto_month": 0.3,
    "p_uncertain_upto_year": 0.1,
    "p_certain": 0.6
  }
}
```

Just the simulation parameters (`p_certain` = `1 - p_uncertain_upto_month - p_uncertain_upto_year`,
included for completeness) — no per-tip breakdown (it can be reconstructed by comparing
normal and `-COMPLETE` output files if needed).

When masking is active (missing data, tip-date uncertainty, or both), add top-level fields for
the complete file paths:

```json
{
  "complete_fasta": "foo-COMPLETE.fasta",
  "complete_maple": "foo-COMPLETE.maple",
  "missing_data": { ... },
  "tip_date_uncertainty": { ... }
}
```

The `complete_fasta` / `complete_maple` fields are only present if those outputs were requested.
Move them out of the `missing_data` section to the top level, since they now serve both features.

**Breaking change:** The `complete_fasta` and `complete_maple` fields previously lived inside the
`missing_data` section.  They are now at the top level of the info JSON.  This is acceptable at
this early stage of development.

## Pipeline Integration

**Seed stability:** Tip-date uncertainty draws from the RNG *after* `simulate_missing_data()`
completes (which itself comes after `simulate_mutations()`).  This means:
- Adding or changing tip-date uncertainty parameters does not alter the coalescent tree, the
  mutations, the complete sequences, or the missing data masking.
- A run with `--p-tip-date-uncertain-upto-month 0.3` at seed 42 produces the same `-COMPLETE`
  files as a run without tip-date uncertainty at seed 42.

In `main()`, after missing data simulation:

```
1. Simulate tip-date uncertainty (if active).
2. Output info JSON (always uses full names, never swapped).
3. If masking is active (missing data or tip-date uncertainty):
   a. Output -COMPLETE FASTA/MAPLE files (full dates, full sequences).
   b. Swap to uncertain tip names.
   c. Output Newick file (uncertain dates in names).
   d. Output Nexus file (uncertain dates in names).
   e. If missing data is active:
      - Output normal FASTA with missing data (uncertain dates + masked sequences).
      - Output normal MAPLE with missing data (uncertain dates + masked sequences).
   f. Else:
      - Output normal FASTA (uncertain dates, full sequences).
      - Output normal MAPLE (uncertain dates, full sequences).
   g. Restore full tip names.
4. If masking is NOT active:
   - Output all files normally (current behavior, unchanged order).
```

The output order within the masking-active branch preserves the existing ordering as much as
possible: info first, then Newick, then Nexus, then FASTA, then MAPLE.  The only additions are
the `-COMPLETE` files (before Newick) and the name swap/restore around the outputs that need
uncertain dates.

## Files Changed

| File | Change |
|------|--------|
| `sapling.cpp` | Add CLI options, `Date_uncertainty`/`Tip_date_uncertainty` types, `is_tip_date_uncertainty_active()`, `is_masking_active()`, `truncate_date_in_name()`, `apply_tip_date_uncertainty_to_names()` (returns saved names map), `restore_tip_names()`, `simulate_tip_date_uncertainty()`, update `dump_info()` (add `tip_date_uncertainty` section, move `complete_*` to top level; no signature change needed), restructure `main()` output section |

No new core library files.  No test changes.  No changes to `output_newick_tree()`.

## Implementation Order

1. Add CLI options and `Options` fields
2. Add data structures (`Date_uncertainty`, `Tip_date_uncertainty`)
3. Implement `is_tip_date_uncertainty_active()` and `is_masking_active()`
4. Implement `truncate_date_in_name()`
5. Implement `apply_tip_date_uncertainty_to_names()` (returns saved names) and `restore_tip_names()`
6. Implement `simulate_tip_date_uncertainty()`
7. Update `dump_info()` to include `tip_date_uncertainty` section and move `complete_*` to top level
8. Restructure `main()` output section with temporary name swapping
9. Build and manual integration test

## Verification

1. Build: `cmake --build --preset conan-debug`
2. Run without uncertainty (should be identical to current behavior):
   ```
   ./build/debug/sapling --const-pop-n0 1.0 -n 10 -L 1000 --out-prefix /tmp/test
   ```
3. Run with missing data only (verify `complete_*` fields moved to top level in info JSON):
   ```
   ./build/debug/sapling --const-pop-n0 1.0 -n 10 -L 10000 \
     --missing-data-mean-num-gaps 2.0 --missing-data-mean-gap-length 500 \
     --out-prefix /tmp/test-missing
   ```
4. Run with tip-date uncertainty only:
   ```
   ./build/debug/sapling --const-pop-n0 1.0 -n 20 -L 1000 \
     --p-tip-date-uncertain-upto-month 0.3 --p-tip-date-uncertain-upto-year 0.1 \
     --out-prefix /tmp/test-dates
   ```
5. Run with both missing data and tip-date uncertainty:
   ```
   ./build/debug/sapling --const-pop-n0 1.0 -n 20 -L 10000 \
     --missing-data-mean-num-gaps 2.0 --missing-data-mean-gap-length 500 \
     --p-tip-date-uncertain-upto-month 0.3 --p-tip-date-uncertain-upto-year 0.1 \
     --out-prefix /tmp/test-both
   ```
6. Verify `-COMPLETE` files have full `YYYY-MM-DD` dates in all tip names
7. Verify normal FASTA/MAPLE files have a mix of `YYYY-MM-DD`, `YYYY-MM`, and `YYYY` dates
8. Verify Newick/Nexus files have uncertain dates in tip names
9. Verify `-COMPLETE` files from runs with and without uncertainty are identical (same seed)
10. Verify info JSON has correct `tip_date_uncertainty` section and top-level `complete_*` fields
11. Verify missing-data-only run has `complete_*` at top level (not inside `missing_data`)
