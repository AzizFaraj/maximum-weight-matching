cmake_minimum_required(VERSION 3.20)

project(MaximumWeightMatchingAlgorithms VERSION 1.0.0 LANGUAGES CXX)

include(CTest)

option(MWM_NATIVE_OPTIMIZATION "Compile with -march=native" OFF)
option(MWM_PAPER_CPU_TARGET "Compile with the fixed x86-64-v3 paper CPU target" OFF)
option(MWM_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library(mwm_algorithms STATIC
  src/Graph.cc
  src/MatchingAlgorithms.cc)

target_compile_features(mwm_algorithms PUBLIC cxx_std_17)
target_compile_definitions(mwm_algorithms PRIVATE
  MBX_NO_HOTPROF
  MBX_NO_DIAGNOSTIC_OUTPUT
  MBX_SHARED_D3_TIE_ORDER=1
  MBX_GMG_D3_OWNERSHIP_EXPERIMENT=1
  SCALING_S5_FULL_KEY_DRIVER_DISABLE=1
  SCALING_EVENT_INITIAL_DISSOLUTION_FRONTIER=1
  SCALING_ZERO_ROOT_CLEANUP_FRONTIER=1)
target_include_directories(mwm_algorithms
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/classic"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/scaling")
target_link_libraries(mwm_algorithms PUBLIC OpenMP::OpenMP_CXX)

if(MWM_NATIVE_OPTIMIZATION AND MWM_PAPER_CPU_TARGET)
  message(FATAL_ERROR "Choose only one of MWM_NATIVE_OPTIMIZATION and MWM_PAPER_CPU_TARGET")
endif()

if(MWM_NATIVE_OPTIMIZATION OR MWM_PAPER_CPU_TARGET)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(MWM_PAPER_CPU_TARGET)
      target_compile_options(mwm_algorithms PRIVATE -march=x86-64-v3 -mtune=generic)
    else()
      target_compile_options(mwm_algorithms PRIVATE -march=native)
    endif()
  else()
    message(FATAL_ERROR "CPU optimization options require GCC or Clang")
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
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME frozen_source_identity
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_source_identity.py")
  endif()

  add_executable(test_exact_algorithms tests/test_exact_algorithms.cc)
  target_link_libraries(test_exact_algorithms PRIVATE mwm_algorithms)
  add_test(NAME exact_algorithms COMMAND test_exact_algorithms)

  add_executable(test_randomized tests/test_randomized.cc)
  target_link_libraries(test_randomized PRIVATE mwm_algorithms)
  add_test(NAME randomized_bruteforce COMMAND test_randomized)

  add_executable(test_membership_regression tests/test_membership_regression.cc)
  target_link_libraries(test_membership_regression PRIVATE mwm_algorithms)
  add_test(NAME membership_regression
    COMMAND test_membership_regression
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/membership_regression_676.mtx")

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Keep assertion-based correctness gates active in Release builds.
    target_compile_options(test_exact_algorithms PRIVATE -UNDEBUG)
    target_compile_options(test_randomized PRIVATE -UNDEBUG)
    target_compile_options(test_membership_regression PRIVATE -UNDEBUG)
  endif()
endif()
