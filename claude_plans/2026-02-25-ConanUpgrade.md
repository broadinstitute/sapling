# Upgrade from Conan 1.65 to Conan 2.25

## Context

Sapling currently uses Conan 1.65.0, which is end-of-life. Delphy has already migrated to Conan 2.x. This upgrade brings Sapling in line with Delphy and removes Conan 1.x workarounds (e.g., the macOS Apple Clang version override).

All work in the `feature/ConanUpgrade` branch, worktree at `../sapling-ConanUpgrade`.

## Phase 1: Core build changes (do first, verify before proceeding)

### 1. `conanfile.txt` — update generators and options syntax

Replace Conan 1.x generators with Conan 2.x equivalents:

```
[requires]
boost/1.83.0
eigen/3.4.0

[options]
boost/*:header_only=True

[generators]
CMakeDeps
CMakeToolchain
```

Changes:
- `cmake_find_package` + `cmake_paths` → `CMakeDeps` + `CMakeToolchain`
- `boost:header_only=True` → `boost/*:header_only=True` (Conan 2.x syntax)

### 2. `CMakeLists.txt` — remove `conan_paths.cmake` include

Remove line 21:
```cmake
include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
```

Conan 2.x uses CMake presets and toolchain files instead. No manual include is needed.

### 3. `.gitignore` — ignore generated `CMakeUserPresets.json`

Conan 2.x generates a `CMakeUserPresets.json` in the source directory. Add it to `.gitignore`.

### Verification (Phase 1)

Use a venv to install Conan, as described in INSTALL.md:

```bash
cd ../sapling-ConanUpgrade
python3 -m venv conan-venv
source conan-venv/bin/activate
pip3 install 'conan==2.25'
conan profile detect --force
conan install . --output-folder=build/debug --build=missing --settings=build_type=Debug
cmake --preset conan-debug
cmake --build --preset conan-debug
./build/debug/tests/tests
```

Also run the demo command from INSTALL.md to verify the binary works. Fix any errors until the build is stable before proceeding to Phase 2.

## Phase 2: Update docs and CI (only after Phase 1 is verified)

### 4. `INSTALL.md` — update build instructions

- Conan version: `1.65.0` → `2.25`
- Profile setup: `conan profile new default --detect` + `conan profile update ...` → `conan profile detect`
- Install: `conan install ../..` → `conan install . --output-folder=build/debug --build=missing --settings=build_type=Debug`
- CMake configure: `cmake ../.. -DCMAKE_BUILD_TYPE=Debug` → `cmake --preset conan-debug`
- Build: `make -j 6` → `cmake --build --preset conan-debug`
- Remove the `libstdc++11` profile update step (Conan 2.x auto-detects this)

### 5. `CLAUDE.md` — update build instructions section

Mirror the same changes as INSTALL.md in the "Building and running tests" section.

### 6. `.github/workflows/ci.yml` — update CI pipeline

- `CONAN_VERSION: '1.65.0'` → `CONAN_VERSION: '2.25'`
- Cache path: `~/.conan` → `~/.conan2`
- Linux profile: replace `conan profile new` + `conan profile update` with `conan profile detect --force`
- macOS profile: replace with `conan profile detect --force` (remove Apple Clang version override)
- Build steps (both static and non-static):
  - `conan install . --output-folder=build/release --build=missing --settings=build_type=Release`
  - `cmake --preset conan-release` (with static flags appended for the static variant)
  - `cmake --build --preset conan-release`
