cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

# Header-only: link the `lrsnash` target and include "lib.h". The arithmetic
# type is chosen in code (LRSNash::Solver64, Solver128, SolverMP, SolverGMP
# or Solver). GMP is optional; without it there is no SolverGMP, and Solver
# falls back to lrsmp.
add_library(lrsnash INTERFACE)
target_include_directories(lrsnash INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(lrsnash INTERFACE cxx_std_17)

find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
  target_include_directories(lrsnash INTERFACE ${GMP_INCLUDE_DIR})
  target_link_libraries(lrsnash INTERFACE ${GMP_LIBRARY})
endif()

# The benchmark and the self-contained tests in tools/ (they need GMP). The
# tests against frozen references and against the legacy library are run by
# hand; see CONTEXT.md, "How to re-verify quickly".
option(LRSNASH_BUILD_TOOLS "Build the benchmark and tests in tools/" OFF)
if(LRSNASH_BUILD_TOOLS AND GMP_LIBRARY)
  enable_testing()
  add_executable(lrsnash_bench tools/bench.cc)
  target_link_libraries(lrsnash_bench PRIVATE lrsnash)
  foreach(t test_kernel test_lh test_lp)
    add_executable(${t} tools/${t}.cc)
    target_link_libraries(${t} PRIVATE lrsnash)
  endforeach()
  add_test(NAME kernel COMMAND test_kernel 1000 1)
  add_test(NAME lemke_howson COMMAND test_lh 200 7 1)
  add_test(NAME lp COMMAND test_lp 300 9 1)
  set_tests_properties(kernel lemke_howson lp PROPERTIES
    FAIL_REGULAR_EXPRESSION "(mismatch[a-z_]*|not_eq|unsound_bound)=[1-9]|MISMATCH|NOT AN EQUILIBRIUM")
endif()
