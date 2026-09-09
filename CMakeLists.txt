cmake_minimum_required(VERSION 3.20)

project(MWMExperimentalEvaluation VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

option(MWM_NATIVE_OPTIMIZATION "Compile with -march=native" OFF)
option(MWM_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)
option(MWM_BUILD_EXTENDED_TESTS "Build the longer randomized certificate gate" OFF)
option(MWM_BUILD_LEMON_ORACLE "Build the optional external LEMON oracle" OFF)
set(MWM_LEMON_ROOT "" CACHE PATH "LEMON installation prefix")
set(MWM_SOURCE_REVISION "" CACHE STRING "Source identity embedded in result records")

if(NOT MWM_SOURCE_REVISION)
  execute_process(
    COMMAND git rev-parse --verify HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE MWM_SOURCE_REVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
endif()
if(NOT MWM_SOURCE_REVISION)
  set(MWM_SOURCE_REVISION "unknown")
endif()

file(GLOB_RECURSE MWM_IDENTITY_FILES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/drivers/*.cc")
list(SORT MWM_IDENTITY_FILES)
set(MWM_SOURCE_HASH_INPUT "")
foreach(source_file IN LISTS MWM_IDENTITY_FILES)
  file(SHA256 "${source_file}" source_hash)
  file(RELATIVE_PATH source_rel "${CMAKE_CURRENT_SOURCE_DIR}" "${source_file}")
  string(APPEND MWM_SOURCE_HASH_INPUT "${source_rel}:${source_hash}\n")
endforeach()
string(SHA256 MWM_SOURCE_TREE_DIGEST "${MWM_SOURCE_HASH_INPUT}")
string(SUBSTRING "${MWM_SOURCE_TREE_DIGEST}" 0 12 MWM_SOURCE_TREE_SHORT)
set(MWM_SOURCE_STATE "${MWM_SOURCE_REVISION}+tree-${MWM_SOURCE_TREE_SHORT}")

set(MWM_BUILD_CONFIG "cmake-${CMAKE_BUILD_TYPE}")
if(MWM_NATIVE_OPTIMIZATION)
  string(APPEND MWM_BUILD_CONFIG "-native")
else()
  string(APPEND MWM_BUILD_CONFIG "-portable")
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library(mwm_build_options INTERFACE)
target_compile_features(mwm_build_options INTERFACE cxx_std_17)
target_compile_definitions(mwm_build_options INTERFACE
  MBX_NO_HOTPROF
  MBX_NO_DIAGNOSTIC_OUTPUT
  HPC_SOURCE_COMMIT="${MWM_SOURCE_REVISION}"
  HPC_SOURCE_STATE="${MWM_SOURCE_STATE}"
  HPC_BUILD_CONFIG="${MWM_BUILD_CONFIG}")
target_link_libraries(mwm_build_options INTERFACE OpenMP::OpenMP_CXX)

if(MWM_NATIVE_OPTIMIZATION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(mwm_build_options INTERFACE -march=native)
  else()
    message(FATAL_ERROR "MWM_NATIVE_OPTIMIZATION is supported only by GNU or Clang")
  endif()
endif()

message(STATUS "MWM source identity: ${MWM_SOURCE_STATE}")
message(STATUS "MWM build configuration: ${MWM_BUILD_CONFIG}")

if(MWM_ENABLE_SANITIZERS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(mwm_build_options INTERFACE
      -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(mwm_build_options INTERFACE
      -fno-omit-frame-pointer -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "Sanitizers are supported only by GNU or Clang")
  endif()
endif()

add_library(mwm_core STATIC
  src/Graph.cc
  src/MatchingEngine.cc
  src/Stats.cc)
target_include_directories(mwm_core PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(mwm_core PUBLIC mwm_build_options)

add_executable(mwm_benchmark benchmarks/drivers/mwm_benchmark.cc)
target_link_libraries(mwm_benchmark PRIVATE mwm_core)
# This is the exact practical Scaling configuration used for headline timings.
target_compile_definitions(mwm_benchmark PRIVATE
  SCALING_S5_FULL_KEY_DRIVER_DISABLE)

set_target_properties(mwm_benchmark PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(MWM_BUILD_LEMON_ORACLE)
  find_path(LEMON_INCLUDE_DIR lemon/matching.h
    HINTS "${MWM_LEMON_ROOT}/include")
  find_library(LEMON_LIBRARY lemon
    HINTS "${MWM_LEMON_ROOT}/lib")
  if(NOT LEMON_INCLUDE_DIR)
    message(FATAL_ERROR "LEMON headers were not found; set MWM_LEMON_ROOT")
  endif()
  add_executable(lemon_oracle benchmarks/drivers/lemon_oracle.cc)
  target_include_directories(lemon_oracle PRIVATE "${LEMON_INCLUDE_DIR}")
  target_link_libraries(lemon_oracle PRIVATE mwm_build_options)
  if(LEMON_LIBRARY)
    target_link_libraries(lemon_oracle PRIVATE "${LEMON_LIBRARY}")
  endif()
  set_target_properties(lemon_oracle PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

if(BUILD_TESTING)
  function(add_mwm_test target source)
    add_executable("${target}" "${source}")
    target_link_libraries("${target}" PRIVATE mwm_core)
    add_test(NAME "${target}" COMMAND "${target}")
  endfunction()

  add_mwm_test(test_gabow_backend tests/testGabowBackend.cc)
  add_mwm_test(test_hybrid_identity tests/testHybridIdentity.cc)
  add_mwm_test(test_scaling_input_contract tests/testScalingInputContract.cc)
  add_mwm_test(test_scaling_scale_endpoint tests/testScalingScaleEndpoint.cc)
  add_mwm_test(test_two_copy_reduction tests/testTwoCopyReduction.cc)

  foreach(engine edmonds gabow scaling)
    add_test(
      NAME "smoke_${engine}"
      COMMAND mwm_benchmark "${engine}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/perfect_k4.mtx"
        general file 1000 1 1 1)
    set_tests_properties("smoke_${engine}" PROPERTIES
      PASS_REGULAR_EXPRESSION "RESULT status=success .*exact-weight=19")
  endforeach()

  if(MWM_BUILD_EXTENDED_TESTS)
    add_executable(test_randomized_certificates tests/xvalidateCertificates.cc)
    target_link_libraries(test_randomized_certificates PRIVATE mwm_core)
    add_test(NAME randomized_certificates
      COMMAND test_randomized_certificates 10000 80
        "${CMAKE_BINARY_DIR}/certificate_failure.mtx" 4 1000)
    set_tests_properties(randomized_certificates PROPERTIES LABELS extended)
  endif()
endif()
