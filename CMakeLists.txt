# Sapling cmake build definition
# ==============================

cmake_minimum_required(VERSION 3.17.5)

# To make a release:
# 1. Commit final changes and note down commit hash X
# 2. Bump VERSION and BUILD_NUMBER below
# 3. Set BUILD_PREV_GIT_HASH to X
# 4. Commit these changes to CMakeLists.txt and push (git push)
# 5. Tag Git commit (git tag ${VERSION} && git push origin ${VERSION})
# 6. Compile & distribute the release artifacts
project(sapling VERSION 0.0.0)
set(BUILD_NUMBER 0)
set(BUILD_PREV_GIT_HASH 1234567)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Dependency management
include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
cmake_policy(SET CMP0135 NEW)  # Idiosyncracies about downloaded file timestamps
include(FetchContent)

# Git hash in version string (adapted from https://jonathanhamberg.com/post/cmake-embedding-git-hash/)
include(${CMAKE_SOURCE_DIR}/CheckGit.cmake)
CheckGitSetup()

# Set up common options
add_library(common INTERFACE)

# Set up compiler and linker flags
target_compile_options(common INTERFACE -Wall -Wextra -Werror)  # Treat all warnings as errors
add_compile_options(-fno-omit-frame-pointer -mno-omit-leaf-frame-pointer)  # Sane stack traces
add_compile_options(-ffast-math -msse3)
# add_compile_options(
#   -fsanitize=address
#   -fno-optimize-sibling-calls
#   -fsanitize-address-use-after-scope
#   -fsanitize=undefined)
# add_link_options(
#   -fsanitize=address
#   -fno-optimize-sibling-calls
#   -fsanitize-address-use-after-scope
#   -fsanitize=undefined)

# Find required and optional libraries
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
set(ABSL_PROPAGATE_CXX_STD ON)
add_subdirectory(third-party/abseil-cpp EXCLUDE_FROM_ALL)

# To update to latest nlohmann_json release:
# 1. Lookup latest release, e.g., at the top of https://github.com/nlohmann/json/releases
# 2. Adapt the url below to use the latest release version
# 3. Update the SHA256 of the artifact to match what's in the release description
FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
  URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
  )
FetchContent_MakeAvailable(json)

# Core library sources
add_subdirectory(core)

# Tests
# To update to latest GoogleTest release:
# 1. Lookup latest commit SHA, e.g., at the top of https://github.com/google/googletest/commits/main
# 2. Adapt the url below to use the commit SHA & bump timestamp in comment
# 3. Recalculate and update the SHA256 of the artifact by running:
#      curl -sL --output - https://github.com/google/googletest/archive/<COMMIT_SHA>.zip | sha256sum -
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/b3a9ba2b8e975550799838332803d468797ae2e1.zip # 2023-12-07
  URL_HASH SHA256=93b52882cff4e8fb344070106facee0df5e61be91b1948a69c8581b3bcbe0c61
  )
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(tests)

# cli executable
add_executable(sapling sapling.cpp)
target_link_libraries(sapling
  sapling_core
  nlohmann_json::nlohmann_json
  )

