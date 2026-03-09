# Missing Data Simulation

## Context

Real sequencing data contains missing data: regions of the genome where bases could not be
confidently called.  To make Sapling's simulated data more realistic, we add the ability to
simulate missing data for each tip sequence.  Missing data comes from two sources: (1) contiguous
"gaps" (e.g., from amplicon dropout) and (2) individual missing sites (e.g., from low coverage at
a single position).

## Overview

For each tip, independently simulate:
- A Poisson-distributed number of **gaps**, each with an exponentially-distributed length and
  a uniformly-distributed start position across the genome.
- A Poisson-distributed number of individual **missing sites**, each at a uniformly-distributed
  position across the genome.

Use a "level" approach to merge overlapping gaps and individual sites into a sorted list of
non-overlapping missing ranges per tip.  Then output:
- **Complete** FASTA/MAPLE files (no missing data) with a `-COMPLETE` suffix before the extension.
- **Normal** FASTA/MAPLE files with missing sites replaced by `N`.
- **Nexus** annotations with per-tip `missing_data_ranges`.
- An updated **info JSON** with a `missing_data` section.

When all missing data parameters are at their defaults (0.0), no missing data is simulated,
no `-COMPLETE` files are produced, and the info JSON has no `missing_data` section.  The output
is identical to the current behavior.

## Data Structures

```cpp
// A half-open range [start, end) of sites, 0-indexed
struct Site_range {
  Site_index start;
  Site_index end;
};

// Per-tip missing data: the raw simulation inputs (for info JSON) and merged ranges (for output)
struct Tip_missing_data {
  std::vector<Site_range> gaps;                // raw gaps as [start, end), for info
  std::vector<Site_index> missing_sites;       // raw individual sites, for info
  std::vector<Site_range> missing_ranges;      // merged non-overlapping [start, end)
};

// All missing data, keyed by tip Node_index
using Missing_data = absl::flat_hash_map<Node_index, Tip_missing_data>;
```

These go directly in `sapling.cpp` (no new core file needed).

## Algorithm: Simulating Missing Data for One Tip

```
Input: num_sites (L), mean_num_gaps, mean_gap_length, mean_num_missing_sites, rng
Output: Tip_missing_data

1. Draw num_gaps ~ Poisson(mean_num_gaps).
2. For each gap:
   a. Draw raw_length ~ Exponential(mean_gap_length), length = round(raw_length).
      Clamp to at most L.  If length == 0, skip this gap (do not insert it).
   b. Draw start ~ Uniform(0, L - length - 1), so the gap fits entirely within the genome.
   c. end = start + length.
   d. Record [start, end) in gaps list.
   e. Add level events: +1 at start, -1 at end.

3. Draw num_missing ~ Poisson(mean_num_missing_sites).
4. For each missing site:
   a. site ~ Uniform(0, L-1)
   b. Record site in missing_sites list.
   c. Add level events: +1 at site, -1 at site+1.

5. Sort level events by position (ties broken arbitrarily -- both orderings are correct).
6. Sweep through events, tracking cumulative level.  Merge contiguous regions where level > 0
   into missing_ranges.

Every tip gets an entry in the Missing_data map when missing data is active, even if 0 gaps
and 0 missing sites were drawn for that tip (the entry will have empty vectors).
```

## CLI Options

Add a new option group `"Missing data"` with:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--missing-data-mean-num-gaps` | `double` | `0.0` | Mean number of gaps per tip (Poisson) |
| `--missing-data-mean-gap-length` | `double` | `1000.0` | Mean gap length in sites (Exponential) |
| `--missing-data-mean-num-missing-sites` | `double` | `0.0` | Mean number of individually missing sites per tip (Poisson) |

Store in `Options` as three `double` fields.  Validate that all three are non-negative.
Additionally, validate that `mean_gap_length > 0` when `mean_num_gaps > 0` (an exponential
distribution with mean 0 is degenerate).

Missing data is considered **active** when `mean_num_gaps > 0.0` or `mean_num_missing_sites > 0.0`.

## Output Changes

### Complete Files

When missing data is active and FASTA/MAPLE output is requested, also write complete versions:
- `foo.fasta` -> also write `foo-COMPLETE.fasta`
- `foo.maple` -> also write `foo-COMPLETE.maple`

The complete files use the existing `output_fasta()` and `output_maple()` functions unchanged.
The normal files use new functions that incorporate missing data.

### FASTA with Missing Data

New function `output_fasta_with_missing_data(os, tree, missing_data)`:
- Same traversal as `output_fasta()`.
- When outputting a tip's sequence, check each site against the tip's `missing_ranges`.
  Since `missing_ranges` is sorted and non-overlapping, use a pointer/index that advances
  through ranges as we scan sites left-to-right.  Output `N` for missing sites, the real
  base otherwise.

### MAPLE with Missing Data

New function `output_maple_with_missing_data(os, tree, missing_data)`:
- Same traversal as `output_maple()`.
- When outputting a tip's deltas:
  1. Collect the sorted deltas (as currently done).
  2. Merge the deltas with the missing ranges.  For each missing range, output
     `N\t<1-indexed start>\t<length>`.  Skip any deltas that fall within a missing range.
     Output non-missing deltas normally as `<base>\t<1-indexed site>`.

### Nexus with Missing Data

When missing data is active, `output_newick_tree()` (in annotated mode) adds a
`missing_data_ranges` annotation to each tip node.  The ranges are encoded as a flat list
of numbers (pairs of start, end values) to avoid confusing existing Nexus parsers like
FigTree and `baltic`:

```
TIP_1|2020-06-15[&missing_data_ranges={100,1100,5000,5500}]:0.5
```

The values are 0-indexed half-open intervals, read as consecutive pairs: `[100,1100)` and
`[5000,5500)`.  If a tip has no missing ranges, the annotation is omitted.

**Why flat lists, not nested `{{start,end},{start,end}}`?**  Nested brace syntax in Nexus
annotations is poorly supported across the phylogenetics ecosystem:
- **baltic** (Gytis Dudas): can output nested braces but its input parser doesn't robustly
  handle them.
- **DendroPy**: has a documented bug where it truncates at the first closing brace in
  nested lists.
- **ETE2/ETE3**: cannot handle BEAST-style annotations at all without preprocessing.
- **FigTree**: primarily targets single-level BEAST-style lists; nested support is unclear.

Flat single-level lists `{a,b,c,d}` are the safe, widely-supported format.

This requires passing the `Missing_data` map to `output_newick_tree()` (only used when
`annotated=true`).

### Info JSON

When missing data is active, add a `"missing_data"` object to the info JSON:

```json
{
  "missing_data": {
    "mean_num_gaps": 2.0,
    "mean_gap_length": 1000.0,
    "mean_num_missing_sites": 5.0,
    "complete_fasta": "foo-COMPLETE.fasta",
    "complete_maple": "foo-COMPLETE.maple",
    "per_tip": {
      "TIP_1|2020-06-15": {
        "gaps": [[100, 1100], [5000, 5500]],
        "missing_sites": [200, 3000, 7500],
        "missing_ranges": [[100, 1100], [3000, 3001], [5000, 5500], [7500, 7501]]
      },
      ...
    }
  }
}
```

The `complete_fasta`/`complete_maple` fields are only present if those outputs were requested.
Gap and range coordinates are 0-indexed half-open intervals.  Individual missing sites are
0-indexed.

## Pipeline Integration

**Seed stability:** Missing data simulation draws from the RNG *after* `simulate_mutations()`
completes.  This means that adding or changing missing data parameters does not alter the
coalescent tree, the mutations, or the complete sequences — only the masking changes.  A run
with `--missing-data-mean-num-gaps 2.0` at seed 42 produces the same `-COMPLETE` files as a
run without missing data at seed 42.  This is a deliberate design choice that makes it easy
to study the effect of missing data in isolation.

In `main()`, after `simulate_mutations()` and before output:

```
1. If missing data is active:
   a. For each tip in the tree, call simulate_tip_missing_data(...).
   b. Collect results into Missing_data map.

2. Output:
   a. If missing data is active and FASTA output requested:
      - Write complete FASTA to foo-COMPLETE.fasta (using output_fasta)
      - Write masked FASTA to foo.fasta (using output_fasta_with_missing_data)
   b. If missing data is active and MAPLE output requested:
      - Write complete MAPLE to foo-COMPLETE.maple (using output_maple)
      - Write masked MAPLE to foo.maple (using output_maple_with_missing_data)
   c. If missing data is NOT active:
      - Write FASTA/MAPLE normally (current behavior, no -COMPLETE files)
   d. Info JSON: output as before, with missing_data section added if active.
   e. Newick: output unchanged (no missing data info).
   f. Nexus: pass missing_data to output_newick_tree for tip annotations.
```

## Files Changed

| File | Change |
|------|--------|
| `sapling.cpp` | Add CLI options, `Site_range`/`Tip_missing_data`/`Missing_data` types, `simulate_tip_missing_data()`, `output_fasta_with_missing_data()`, `output_maple_with_missing_data()`, `make_complete_filename()`, update `output_newick_tree()` (missing_data_ranges annotation on tips), update `dump_info()`, update `main()` |

No new core library files.  No test changes (missing data simulation is stochastic and best
validated by integration testing with a sample run).

## Implementation Order

1. Add CLI options and `Options` fields
2. Add data structures (`Site_range`, `Tip_missing_data`, `Missing_data`)
3. Implement `simulate_tip_missing_data()`
4. Implement `make_complete_filename()`
5. Implement `output_fasta_with_missing_data()`
6. Implement `output_maple_with_missing_data()`
7. Update `dump_info()` to include `missing_data` section
8. Update `main()` pipeline
9. Build and manual integration test

## Verification

1. Build: `cmake --build --preset conan-debug`
2. Run without missing data (should be identical to current behavior):
   ```
   ./build/debug/sapling --const-pop-n0 1.0 -n 10 -L 1000 --out-prefix /tmp/test
   ```
3. Run with missing data:
   ```
   ./build/debug/sapling --const-pop-n0 1.0 -n 10 -L 10000 \
     --missing-data-mean-num-gaps 2.0 --missing-data-mean-gap-length 500 \
     --missing-data-mean-num-missing-sites 5.0 \
     --out-prefix /tmp/test-missing
   ```
4. Verify `-COMPLETE` files match what would be output without missing data
5. Verify normal FASTA has `N` characters at expected positions
6. Verify normal MAPLE has `N` ranges at expected positions
7. Verify info JSON has correct `missing_data` section
