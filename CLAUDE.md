# Claude Code Instructions

## Overview

Sapling simulates pandemic-scale sequencing data with known parameters.  Given a population model, a sampling strategy, an HKY substitution model and a genome size, it produces a coalescent genealogy and evolves sequences along it using the Gillespie algorithm.  The output is an annotated Newick tree and sequences in FASTA or MAPLE format.  Because it uses an Explicit Mutation-Annotated Tree internally (never materializing the full alignment), it can generate very large datasets with modest resources.

---

## CRITICAL: Read CONTRIBUTING.md First

**Before writing any code, you MUST read [CONTRIBUTING.md](CONTRIBUTING.md) in full.**

It defines the coding style for this project. Failure to follow these conventions will result in rejected code.

---

## Key Coding Conventions

These are the most important conventions from CONTRIBUTING.md that you MUST follow:

| Convention | Correct | Incorrect |
|------------|---------|-----------|
| Variable/function names | `snake_case` | `camelCase` |
| Type names (classes, structs, aliases) | `Stroustrup_case` | `snake_case` |
| Boolean negation | `not x` | `!x` |
| Natural logarithm | `log(x)` | `ln(x)` |
| Const style | `const char* str` | `char const* str` |
| Pointer/reference | `int* ptr`, `int& ref` | `int *ptr`, `int &ref` |
| Container size | `std::ssize(c)` | `c.size()` |
| Hash maps/sets | `absl::flat_hash_map` | `std::unordered_map` |
| Data containers | structs with public fields | classes with getters/setters |
| Error handling | exceptions | error codes |
| Function style (AAA) | `auto func(int x) -> int;` | `int func(int x);` |

**Compiler flags:** Code must compile cleanly with `-Wall -Wextra -Werror`.

---

## Git Conventions

1. **Never amend pushed commits** - Always create fresh commits. Do not use `git commit --amend` or `git rebase` on commits that have already been pushed.
2. **Feature branches** - Develop each new feature in a separate branch and git worktree, so multiple features can be developed in parallel. Create a worktree with:
   ```bash
   git worktree add -b feature/FeatureName ../sapling-FeatureName main
   ```
   This creates a branch `feature/FeatureName` and checks it out in `../sapling-FeatureName`. Then initialize git submodules in the new worktree:
   ```bash
   cd ../sapling-FeatureName
   git submodule update --init
   ```
   When done, open a PR to merge into `main`.

---

## Planning

Feature plans go in `claude_plans/` with filenames of the form `YYYY-MM-DD-BriefFeatureDescription.md` (e.g., `2026-02-24-SiteRateHeterogeneity.md`).  Use plan mode to draft these before implementing a feature.

**Important:** Always write the plan to a file in `claude_plans/`, not just to the ephemeral plan mode file.  This ensures the plan is visible for direct editing, persists across sessions, and is committed alongside the code.  When in plan mode, write your plan to `claude_plans/` first, then reference it from the plan mode file.

---

## Directory layout

```
sapling/
  sapling.cpp              CLI entry point & simulation pipeline
  version.h / version.cpp.in   Version string (filled by CMake)
  CMakeLists.txt            Top-level build configuration
  CheckGit.cmake            Embeds git commit hash into version string
  conanfile.txt             Conan 2.x dependency list

  core/                     Core library (sapling_core)
    tree.h                    Generic tree data structure + traversals
    phylo_tree.h / .cpp       Phylogenetic tree (= Tree<Phylo_node> + ref_sequence + mutations)
    sequence.h                Seq_letter (IUPAC bitmask), Real_seq_letter (A/C/G/T enum),
                                Seq_vector<T>, Seq_matrix<T>, pick_state()
    sequence_overlay.h        Memory-efficient mutable view over a Real_sequence (base + sparse deltas)
    mutations.h               Tree_loc, Mutation_info, Mutations (btree_multimap), on_branch()
    pop_model.h / .cpp        Pop_model interface + Const_pop_model, Exp_pop_model
    coal_sim.h / .cpp         Coalescent simulator: pop model + tip times -> Phylo_tree
    evo_model.h / .cpp        Site_evo_model (mu, pi_a, q_ab), Global_evo_model (per-site rates)
    evo_hky.h / .cpp          HKY substitution model -> derives Site_evo_model using Eigen
    dates.h / .cpp            ISO date parsing and formatting (parse_iso_date, to_iso_date)
    tip_file.h / .cpp         Tip file parser (Name|YYYY-MM-DD format)
    distributions.h           Bounded_exponential_distribution (used by Exp_pop_model::sample)
    estd.h                    is_debug_enabled flag
    CMakeLists.txt            Builds sapling_core library

  tests/                    Unit tests (GoogleTest)
    tests_main.cpp            GoogleTest main
    tree_tests.cpp            Binary/Nary nodes, traversals, integrity checks
    sequence_tests.cpp        IUPAC conversions, Real_seq_letter round-trips
    sequence_overlay_tests.cpp  Overlay creation, delta tracking, materialize
    mutations_tests.cpp       Mutation insertion, on_branch lookup, erase_mutation
    pop_model_tests.cpp       Const/Exp population model functions, intensity, inverse_intensity
    dates_tests.cpp           ISO date parsing/formatting round-trips
    tip_file_tests.cpp        Tip file parsing (valid/invalid inputs, edge cases)
    CMakeLists.txt

  claude_plans/             Dated feature plans (YYYY-MM-DD-BriefFeatureDescription.md)
  demo-data/                Reference output files for a sample run
  third-party/              Git submodules: abseil-cpp, cxxopts, cppcoro
```

## Key types

| Type | File | Purpose |
|------|------|---------|
| `Node_index` (int) | `tree.h` | Index of a node in a tree's contiguous node array |
| `Branch_index` (int) | `tree.h` | Index of a branch (= index of its child endpoint) |
| `Node<C>` | `tree.h` | Generic tree node; parent index + children container `C` |
| `Binary_node` | `tree.h` | `Node<Binary_child_indices>` -- compact inline 2-child container |
| `Tree<N>` | `tree.h` | Vector of `N` nodes + root index + traversal support |
| `Phylo_node` | `phylo_tree.h` | `Binary_node` + `name` (string) + `t` (double, time) |
| `Phylo_tree` | `phylo_tree.h` | `Tree<Phylo_node>` + `ref_sequence` + `mutations` |
| `Seq_letter` (uint8_t) | `sequence.h` | IUPAC bitmask (A=0001, C=0010, G=0100, T=1000, N=1111) |
| `Real_seq_letter` (enum) | `sequence.h` | Unambiguous base: A=0, C=1, G=2, T=3 |
| `Seq_vector<T>` | `sequence.h` | 4-element array indexed by `Real_seq_letter` (e.g., base frequencies) |
| `Seq_matrix<T>` | `sequence.h` | 4x4 array indexed by pairs of `Real_seq_letter` (e.g., rate matrix) |
| `Real_sequence` | `sequence.h` | `vector<Real_seq_letter>` |
| `Sequence_overlay` | `sequence_overlay.h` | Non-owning view of a `Real_sequence` + sparse `flat_hash_map` of per-site deltas; avoids copying full genomes during tree traversal |
| `Tree_loc` | `mutations.h` | Point on a tree: `{branch, t}` |
| `Mutation_info` | `mutations.h` | `{site, from, to}` |
| `Mutations` | `mutations.h` | `btree_multimap<Tree_loc, Mutation_info>` keyed by branch then time; supports `on_branch(i, mutations)` to iterate over mutations on branch `i` |
| `Pop_model` | `pop_model.h` | Abstract: `pop_at_time(t)`, `intensity_at_time(t)`, `inverse_intensity(I)`, `sample(min_t, max_t, rng)` |
| `Const_pop_model` | `pop_model.h/.cpp` | Constant N_e; intensity I(t) = t/N |
| `Exp_pop_model` | `pop_model.h/.cpp` | N_e(t) = n0 * exp(g*(t-t0)); intensity via log formula |
| `Site_evo_model` | `evo_model.h` | `{mu, pi_a, q_ab}` -- single-site substitution rate matrix |
| `Global_evo_model` | `evo_model.h` | `{nu_l[], site_evo_model}` -- per-site relative rates (all 1.0 for now) |
| `Hky_model` | `evo_hky.h/.cpp` | `{mu, kappa, pi_a}` -> derives `Site_evo_model` via Eigen |
| `Bounded_exponential_distribution` | `distributions.h` | Samples from exponential pdf truncated to [a, b]; used for tip time sampling under exponential growth |

## Main simulation pipeline

All orchestrated in `main()` in `sapling.cpp`:

```
1. process_args(argc, argv)
   Parse CLI options with cxxopts.  Choose population model (const or exp).
   Parse ISO dates for tip time range.  Collect HKY parameters and output paths.
   -> Options struct

2. choose_tip_times(pop_model, num_samples, min_tip_t, max_tip_t, t0, rng)
   Draw N random sample times from the population model's sampling distribution
   (uniform for const, bounded-exponential for exp growth).
   Each time is rounded to the nearest ISO date.
   -> vector<double> tip_times
   Alternatively, if --tip-file is given, parse_tip_file() reads explicit
   tip names and dates from a file (format: Name|YYYY-MM-DD, one per line).
   -> vector<pair<string, double>> with names preserved for step 5

3. coal_sim(pop_model, tip_times, rng)
   Build a coalescent tree backwards in time:
     a. Allocate 2N-1 nodes; first N are tips with their sampling times.
     b. Sort tips by time.  Start at t = max(tip_times).
     c. Maintain active_branches (lineages alive at current t) and pending_tips.
     d. Loop: compare next tip activation time vs next coalescence time.
        - Coalescence time: draw exponential with rate k*(k-1)/2, convert
          through the population model's intensity function (inverse_intensity).
        - If a tip is next: activate it (add to active_branches).
        - If coalescence is next: pick two random active branches, create
          a new inner node joining them, add it as a new active branch.
     e. Last remaining active branch becomes the root.
   -> Phylo_tree (topology + times, no sequences yet)

4. ladderize_tree(tree)
   Post-order pass: swap children so the subtree with fewer descendants is on the left.

5. name_nodes(tree, t0, tip_names={})
   Index-order pass: tips get "TIP_k|YYYY-MM-DD" (or names from tip file if provided),
   inner nodes get "NODE_k|YYYY-MM-DD".

6. gen_random_sequence(num_sites, pi_a, rng)
   Draw each site i.i.d. from the stationary base frequencies pi_a.
   -> tree.ref_sequence (the MRCA sequence)

7. Hky_model::derive_site_evo_model()
   Build the HKY rate matrix q_ab using Eigen:
     r_ab = 1 for transversions, kappa for transitions
     q_ab = r_ab * pi_b / R,  where R = pi^T r pi  (normalization)
     diagonal: q_aa = -(sum of off-diagonal row)
   -> Site_evo_model {mu, pi_a, q_ab}

8. make_single_partition_global_evo_model(num_sites)
   Wrap the site model with per-site rate multipliers nu_l (all 1.0 for now).
   -> Global_evo_model

9. simulate_mutations(tree, evo, rng)
   Gillespie algorithm over the tree, pre-order:
     a. Calculate total event rate lambda = sum_l mu * nu_l * q_a(l) at root.
     b. Use a stack-based pre-order walk.  At each node, for each child branch:
        - t = parent time, t_max = child time
        - Draw next event time: t += Exponential(lambda).  If t >= t_max, done with branch.
        - Pick site: sample cumulative Q until target_cum_Q is reached (linear scan).
        - Pick target state: categorical from off-diagonal q_ab row.
        - Record mutation: tree.mutations.insert({{child, t}, {site, from, to}}).
        - Update seq[site] and adjust lambda incrementally.
     c. Only push inner children onto the work stack (tips are leaves).
   Mutations are stored globally in the tree's btree_multimap, keyed by (branch, time).

10. Output (each guarded by write_to, which skips if filename is absent):
     - Info (JSON):     dump_info() -- model params, tree stats, CLI invocation
     - Newick (.nwk):   output_newick_tree() -- plain topology + branch lengths
     - Nexus (.nexus):  output_newick_tree() with annotations -- root sequence,
                         per-branch mutation lists, inner node dates
     - FASTA (.fasta):  output_fasta() -- full tip sequences
     - MAPLE (.maple):  output_maple() -- reference sequence + per-tip deltas
```

## Output traversal pattern (Newick, FASTA, MAPLE)

All three output functions share the same traversal strategy, using `traversal(tree)` from `tree.h` which yields `(node, children_so_far)` pairs during a depth-first walk:

- **Entering a node** (`children_so_far == 0`): apply all mutations on the branch leading to this node to the current `Sequence_overlay`.
- **Exiting a node** (`children_so_far == num_children`): emit output for this node (for tips: sequence data; for Newick: close parens + name + branch length), then **revert** all mutations on the branch (iterating in reverse).

This apply-then-revert pattern means the traversal maintains a single `Sequence_overlay` that always reflects the sequence at the current node, without ever copying the full genome.

## Module dependency graph

```
sapling.cpp
  |
  +---> coal_sim.{h,cpp}       (coalescent tree generation)
  |       +---> phylo_tree.h    (Phylo_tree, Phylo_node)
  |       |       +---> tree.h  (Tree<N>, Binary_node, traversals)
  |       |       +---> sequence.h
  |       |       +---> mutations.h
  |       +---> pop_model.{h,cpp}
  |               +---> distributions.h
  |
  +---> dates.{h,cpp}          (ISO date parsing / formatting)
  |       +---> absl::time
  |
  +---> tip_file.{h,cpp}       (tip file parser: Name|YYYY-MM-DD)
  |       +---> dates.h
  |       +---> absl::time
  |
  +---> evo_hky.{h,cpp}        (HKY model -> Site_evo_model)
  |       +---> evo_model.h
  |       +---> Eigen3
  |
  +---> evo_model.{h,cpp}      (Global_evo_model)
  |       +---> sequence.h
  |
  +---> sequence_overlay.h     (memory-efficient sequence view)
  |       +---> sequence.h
  |       +---> absl::flat_hash_map
  |
  +---> tree.h                 (traversals, count_all_descendants)
  |       +---> cppcoro/generator.hpp
  |       +---> estd.h
  |       +---> absl::flat_hash_map, absl::check
  |
  +---> cxxopts                (CLI argument parsing)
  +---> nlohmann/json          (info file output)
  +---> absl::*                (logging, RNG, string formatting, time)
```

## Build system

- **CMake** (C++20, `-Wall -Wextra -Werror -ffast-math -msse3`)
- **Conan 2.25** for Eigen3 and Boost
- **FetchContent** for nlohmann/json and GoogleTest
- **Git submodules** for abseil-cpp, cxxopts, cppcoro
- Two build targets: `sapling_core` (static library) and `sapling` (executable linking `sapling_core` + nlohmann_json)
- Tests link `sapling_core` + GoogleTest

## Building and running tests

One-time setup (requires Python 3 for Conan):

```bash
# Create a virtualenv if needed (e.g., Ubuntu >= 24.04)
python3 -m venv delphy-venv && source delphy-venv/bin/activate

# Install Conan 2.25
pip3 install 'conan==2.25'
conan profile detect

# Check out git submodules
git submodule update --init
```

Build (from the repo root):

```bash
conan install . --output-folder=build/debug --build=missing --settings=build_type=Debug
cmake --preset conan-debug
cmake --build --preset conan-debug
```

This produces two executables in `build/debug`:
- `sapling` -- the main CLI tool
- `tests/tests` -- the unit test runner

Run all tests:

```bash
./build/debug/tests/tests
```

Or via CTest:

```bash
cd build/debug && ctest
```

For a release build, replace `debug` with `release` and `Debug` with `Release` in the commands above.

## Releasing a new version

In the top-level `CMakeLists.txt`, update three variables near the top of the file:

1. Commit all final changes and note the commit hash X.
2. Bump `VERSION` (in the `project(sapling VERSION ...)` line) and `BUILD_NUMBER`.
3. Set `BUILD_PREV_GIT_HASH` to X.
4. Commit the `CMakeLists.txt` change and push.
5. Tag the commit: `git tag ${VERSION} && git push origin ${VERSION}`.
6. Compile and distribute the release artifacts.

---

## Dependencies

### Via Conan
- boost/1.83.0 (header-only)
- eigen/3.4.0

### Via Git Submodules (third-party/)
- abseil-cpp - Hash maps, logging, random numbers, string formatting, time utilities
- cxxopts - Command-line option parsing

### Vendored in third-party/
- cppcoro - Coroutines/generators for tree traversals (header-only)

### Via CMake FetchContent
- nlohmann/json v3.11.3 - JSON serialization for info file output
- GoogleTest - Unit testing framework

---

## Testing Patterns

Tests use GoogleTest and GoogleMock, with a custom `tests_main.cpp` that initializes `absl::InitializeLog()` before running tests.  All tests are wrapped in `namespace sapling { ... }`.

```cpp
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "tree.h"

namespace sapling {

TEST(Tree_test, binary_simple) {
  auto tree = Tree<Binary_node>{2};
  auto i = 0;
  auto j = 1;
  auto k = tree.add_node();

  tree.at(i).set_children(j, k);
  tree.at(j).set_parent(i);
  tree.at(k).set_parent(i);

  tree.set_root(i);

  assert_tree_integrity(tree);

  EXPECT_EQ(std::ssize(tree.at(i).children()), 2);
  EXPECT_THAT(tree.at(i).left_child(), testing::Eq(j));
  EXPECT_THAT(tree.at(i).right_child(), testing::Eq(k));
}

}  // namespace sapling
```

Test files follow the pattern `{module}_tests.cpp` in the `tests/` directory.
