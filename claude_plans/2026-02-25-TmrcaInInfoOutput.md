# Plan: Output tMRCA of Root in Info JSON File

## Context

The `--out-info` JSON file contains metadata about a simulation run: version info,
model parameters, and tree statistics (number of mutations, total branch length).
It does not currently include the time of the most recent common ancestor (tMRCA),
i.e., the time of the root node.  This is a key summary statistic that users need
to characterize the simulated tree.

## Change

Add `"t_mrca"` and `"t_mrca_date"` fields to the `"tree_stats"` section of the
info JSON output.

### In `sapling.cpp`, function `dump_info()` (~line 385)

Currently, the `tree_stats` section is:

```cpp
{"tree_stats", {
    {"num_mutations", std::ssize(tree.mutations)},
    {"total_branch_length", calc_total_branch_length(tree)}}}
```

Change it to:

```cpp
{"tree_stats", {
    {"num_mutations", std::ssize(tree.mutations)},
    {"total_branch_length", calc_total_branch_length(tree)},
    {"t_mrca", tree.at_root().t},
    {"t_mrca_date", to_iso_date(tree.at_root().t, opts.t0)}}}
```

- `t_mrca` is the raw numeric time of the root (in years relative to `t0`,
  matching the internal time representation).
- `t_mrca_date` is the corresponding ISO date string (e.g., `"2020-03-15"`),
  using the same `to_iso_date()` function already used elsewhere in the file.

### Example output

Before:
```json
{
  ...
  "tree_stats": {
    "num_mutations": 42,
    "total_branch_length": 12.345
  }
}
```

After:
```json
{
  ...
  "tree_stats": {
    "num_mutations": 42,
    "total_branch_length": 12.345,
    "t_mrca": -3.456,
    "t_mrca_date": "2021-03-15"
  }
}
```

## Files modified

- `sapling.cpp` — two lines added in `dump_info()` (~line 387)

## No new tests needed

The info output is not currently unit-tested (it's a top-level function in
`sapling.cpp`, not in the core library).  This change adds two fields to an
existing JSON object using already-tested building blocks (`tree.at_root().t`
and `to_iso_date()`).

## Verification

1. Build: `cmake --build --preset conan-debug`
2. Run a simulation outputting both info and nexus:
   ```bash
   ./build/debug/sapling --const-pop 1 -n 10 --t0 2024-07-31 \
     --min-tip-t -0.5 --max-tip-t 0 --mu 0.001 -L 1000 \
     --out-info test.json --out-nexus test.nexus 2>/dev/null
   ```
3. Confirm `t_mrca` and `t_mrca_date` appear in the `tree_stats` section of `test.json`.
4. Confirm `t_mrca_date` matches the date annotation on the root node in `test.nexus`.
