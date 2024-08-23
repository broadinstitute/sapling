add_library(sapling_core
  # Source files in alphabetical order
  estd.h
  coal_sim.h coal_sim.cpp
  distributions.h
  evo_hky.h evo_hky.cpp
  evo_model.h evo_model.cpp
  mutations.h
  phylo_tree.h phylo_tree.cpp
  pop_model.h pop_model.cpp
  sequence.h
  sequence_overlay.h
  tree.h
  )

target_link_libraries(sapling_core
  common
  version
  Boost::headers
  Eigen3::Eigen
  absl::check
  absl::flat_hash_map
  absl::log
  absl::log_initialize
  absl::random_random
  absl::random_bit_gen_ref
  absl::str_format
  )

target_include_directories(sapling_core INTERFACE .)
target_include_directories(sapling_core PUBLIC ../third-party/cppcoro/include)
target_include_directories(sapling_core PUBLIC ../third-party/cxxopts/include)

# Uncomment to leave assembler intermediates around
#set_target_properties(sapling_core PROPERTIES COMPILE_FLAGS "-save-temps")
