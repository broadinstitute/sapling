# Plan: Add --tip-file option for explicit tip names and dates

## Context

Currently, sapling generates tip times by randomly sampling from the population model's
distribution (`choose_tip_times()`), and auto-generates names like `TIP_1|2024-07-31`.
This feature adds a `--tip-file` option so users can supply exact tip names and dates,
which is useful when simulating data that matches a real-world sampling scheme.

## Input file format

One line per tip, each line is a tip name of the form `TipNameHere|YYYY-MM-DD`:

```
Alpha|2024-01-15
Beta|2024-03-22
Gamma|2024-07-01
```

## Changes

### 1. Add `core/tip_file.h` / `core/tip_file.cpp`

New core library files with a function to parse the tip file format:

```cpp
auto parse_tip_file(std::istream& is, const absl::Time& t0)
    -> std::vector<std::pair<std::string, double>>;
```

- Reads each non-empty line, splits on the last `|` character.
- Parses the date part into a double (years from epoch), using `absl::ParseTime`
  with format `%Y-%m-%d` (same logic as `parse_iso_date()` in `sapling.cpp`).
- Returns a vector of (name, time) pairs.
- Throws on malformed lines (no `|`, bad date format) or empty input.

Also add an overload that takes a filename:

```cpp
auto parse_tip_file(const std::string& filename, const absl::Time& t0)
    -> std::vector<std::pair<std::string, double>>;
```

### 2. Add `tests/tip_file_tests.cpp`

Unit tests for `parse_tip_file()` covering:

- **Valid input**: multiple tips, correct names and times extracted.
- **Single tip**: file with one line.
- **Empty lines**: blank lines are skipped.
- **Missing `|`**: line with no pipe character produces an error.
- **Bad date**: line with invalid date (e.g., `Foo|not-a-date`) produces an error.
- **Empty file**: no tips produces an error.
- **Name with special characters**: names containing spaces or multiple `|` characters
  (split on the *last* `|`).

### 3. Update `core/CMakeLists.txt`

Add `tip_file.h tip_file.cpp` to the `sapling_core` library source list.

### 4. Update `tests/CMakeLists.txt`

Add `tip_file_tests.cpp` to the test executable source list.

### 5. Add CLI option in `sapling.cpp`

Add `--tip-file` to the "Sampling strategy" option group (after `--max-tip-t`, ~line 103):

```cpp
("tip-file", "File with tip names and dates, one per line (format: Name|YYYY-MM-DD)",
 cxxopts::value<std::string>())
```

### 6. Add field to Options struct in `sapling.cpp`

Add to the Options struct (~line 53, alongside `num_samples`, `min_tip_t`, `max_tip_t`):

```cpp
std::optional<std::string> tip_file;
```

### 7. Mutual exclusion validation in process_args()

In the sampling strategy validation section (~line 210), add a check that `--tip-file`
is mutually exclusive with `-n`, `--min-tip-t`, and `--max-tip-t`. If `--tip-file` is
provided, error if any of the others are also specified. If `--tip-file` is NOT provided,
require `-n` as before.

### 8. Update main() flow in `sapling.cpp`

In `main()` (~line 678), branch on whether `opts.tip_file` is set:

- **If set:** call `parse_tip_file()` to get names and times. Extract the times
  into a `vector<double>` for `coal_sim()`. Store the names for later use.
- **If not set:** use the existing `choose_tip_times()` path (unchanged).

The rest of the pipeline (`coal_sim`, `ladderize_tree`, sequence generation,
mutation simulation, output) is unchanged -- `coal_sim()` just needs a
`span<const double>` of tip times regardless of source.

### 9. Update name_nodes() to accept optional supplied names

Change `name_nodes()` to accept an optional vector of tip names. When a tip
file is provided, the tip nodes (0..N-1) get their names from the file (in the order
they appear in the tip_times vector). `ladderize_tree()` does not reorder
node indices (it only swaps child pointers), so the mapping from node index to
supplied name is stable.

Inner nodes still get auto-generated names (`NODE_k|YYYY-MM-DD`).

Signature change:

```cpp
auto name_nodes(Phylo_tree& tree, const absl::Time& t0,
                const std::vector<std::string>& tip_names = {}) -> void;
```

If `tip_names` is non-empty, use `tip_names[node]` for tip nodes instead of
`TIP_k|YYYY-MM-DD`. If empty, use the existing auto-generated format.

### 10. Update dump_info()

Include the tip file path in the JSON info output when `--tip-file` is used
(in the `sampling_strategy` section).

## Files modified

- `core/tip_file.h` -- new file
- `core/tip_file.cpp` -- new file
- `core/CMakeLists.txt` -- add tip_file sources
- `tests/tip_file_tests.cpp` -- new file
- `tests/CMakeLists.txt` -- add tip_file_tests.cpp
- `sapling.cpp` -- CLI option, Options struct, main() flow, name_nodes(), dump_info()

## Verification

1. Build: `cd build/debug && cmake ../.. -DCMAKE_BUILD_TYPE=Debug && make -j 6`
2. Run unit tests: `cd build/debug && tests/tests`
3. Create a test tip file, e.g., `test_tips.txt`:
   ```
   Alpha|2024-01-15
   Beta|2024-03-22
   Gamma|2024-07-01
   ```
4. Run with tip file:
   ```
   ./sapling --tip-file test_tips.txt --exp-pop-n0 6.0 --exp-pop-g 10 \
     --t0 2024-07-31 --mu 0.1 -L 10 --out-prefix test
   ```
5. Verify tip names in output files (test.nwk, test.fasta, test.maple) match the
   names from the file (Alpha, Beta, Gamma), not auto-generated TIP_N names.
6. Verify error when combining `--tip-file` with `-n`:
   ```
   ./sapling --tip-file test_tips.txt -n 5 --exp-pop-n0 6.0 --exp-pop-g 10 ...
   ```
