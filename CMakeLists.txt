cmake_minimum_required(VERSION 2.8.12)

project(dune-fempy CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR ${PROJECT_BINARY_DIR})
endif()

find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${dune-common_MODULE_PATH})

include(DuneMacros)

dune_project()
dune_enable_all_packages()

add_subdirectory("cmake/modules")
add_subdirectory("dune")
add_subdirectory("python")
add_subdirectory("doc")
add_subdirectory("data")
add_subdirectory("navierstokes")
add_subdirectory("phasefield")

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

if(UNIX)
  message("=============================================================================")
  message("  To use this build of dunefempy within python, simply type")
  message("    export PYTHONPATH=\"${CMAKE_BINARY_DIR}/python:\$PYTHONPATH\"")
  message("  at the command prompt.")
  message("=============================================================================")
endif()
