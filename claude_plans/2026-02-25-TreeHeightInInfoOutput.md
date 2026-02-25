# Plan: Output Tree Height in Info JSON File

## Context

The `--out-info` JSON file's `tree_stats` section currently includes `num_mutations`,
`total_branch_length`, `t_mrca`, and `t_mrca_date`.  This feature adds the tree height,
defined as the time span from the root to the latest (most recent) tip.

## Change

Add a `calc_tree_height()` helper function and a `"tree_height"` field to the
`"tree_stats"` section of the info JSON output.

### 1. Add `calc_tree_height()` in `sapling.cpp` (after `calc_total_branch_length()`, ~line 330)

```cpp
auto calc_tree_height(const Phylo_tree& tree) -> double {
  auto max_tip_t = -std::numeric_limits<double>::infinity();
  for (const auto& node : index_order_traversal(tree)) {
    if (tree.at(node).is_tip()) {
      max_tip_t = std::max(max_tip_t, tree.at(node).t);
    }
  }
  return max_tip_t - tree.at_root().t;
}
```

### 2. Add `"tree_height"` to `dump_info()` (~line 389)

Add `{"tree_height", calc_tree_height(tree)}` to the `tree_stats` JSON object,
after the `t_mrca_date` entry.

### Example output

```json
{
  ...
  "tree_stats": {
    "num_mutations": 42,
    "total_branch_length": 12.345,
    "t_mrca": -3.456,
    "t_mrca_date": "2021-03-15",
    "tree_height": 3.456
  }
}
```

## Files modified

- `sapling.cpp` — new `calc_tree_height()` function (~8 lines) and one line added in `dump_info()`

## Verification

1. Build: `cmake --build --preset conan-debug`
2. Run all unit tests: `./build/debug/tests/tests`
3. Run a simulation:
   ```bash
   ./build/debug/sapling --const-pop-n0 1 -n 10 --t0 2024-07-31 \
     --min-tip-t 2024-01-01 --max-tip-t 2024-07-31 --mu 0.001 -L 1000 \
     --out-info test.json --out-nexus test.nexus 2>/dev/null
   ```
4. Confirm `tree_height` appears in `test.json` and equals the difference between
   the latest tip date and the root date in `test.nexus`.
