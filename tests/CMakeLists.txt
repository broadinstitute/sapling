project(tests)

add_executable(tests
  # Source files in alphabetical order
  mutations_tests.cpp
  pop_model_tests.cpp
  sequence_overlay_tests.cpp
  sequence_tests.cpp
  tests_main.cpp
  tree_tests.cpp
  )

target_link_libraries(tests
  common
  sapling_core
  GTest::gtest
  GTest::gmock)

# Uncomment to leave assembler intermediates around
#set_target_properties(tests PROPERTIES COMPILE_FLAGS "-save-temps")
