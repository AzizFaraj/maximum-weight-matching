cmake_minimum_required(VERSION 3.20)

project(MaximumWeightMatchingAlgorithms VERSION 1.0.0 LANGUAGES CXX)

include(CTest)

option(MWM_NATIVE_OPTIMIZATION "Compile with -march=native" OFF)
option(MWM_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library(mwm_algorithms STATIC
  src/Graph.cc
  src/MatchingAlgorithms.cc)

target_compile_features(mwm_algorithms PUBLIC cxx_std_17)
target_compile_definitions(mwm_algorithms PRIVATE
  MBX_NO_HOTPROF
  MBX_NO_DIAGNOSTIC_OUTPUT
  SCALING_S5_FULL_KEY_DRIVER_DISABLE=1)
target_include_directories(mwm_algorithms
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/classic"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/scaling")
target_link_libraries(mwm_algorithms PUBLIC OpenMP::OpenMP_CXX)

if(MWM_NATIVE_OPTIMIZATION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(mwm_algorithms PRIVATE -march=native)
  else()
    message(FATAL_ERROR "MWM_NATIVE_OPTIMIZATION requires GCC or Clang")
  endif()
endif()

if(MWM_ENABLE_SANITIZERS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(mwm_algorithms PUBLIC
      -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(mwm_algorithms PUBLIC
      -fno-omit-frame-pointer -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "Sanitizers require GCC or Clang")
  endif()
endif()

add_executable(mwm_run app/mwm_run.cc)
target_link_libraries(mwm_run PRIVATE mwm_algorithms)
set_target_properties(mwm_run PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(BUILD_TESTING)
  add_executable(test_exact_algorithms tests/test_exact_algorithms.cc)
  target_link_libraries(test_exact_algorithms PRIVATE mwm_algorithms)
  add_test(NAME exact_algorithms COMMAND test_exact_algorithms)

  add_executable(test_randomized tests/test_randomized.cc)
  target_link_libraries(test_randomized PRIVATE mwm_algorithms)
  add_test(NAME randomized_bruteforce COMMAND test_randomized)

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Keep assertion-based correctness gates active in Release builds.
    target_compile_options(test_exact_algorithms PRIVATE -UNDEBUG)
    target_compile_options(test_randomized PRIVATE -UNDEBUG)
  endif()
endif()
