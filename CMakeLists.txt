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

add_subdirectory("bin")
add_subdirectory("cmake/modules")
add_subdirectory("dune")
add_subdirectory("python")
add_subdirectory("doc")
add_subdirectory("notebooks")
add_subdirectory("data")
add_subdirectory("forchheimer")
add_subdirectory("navierstokes")
add_subdirectory("phasefield")
add_subdirectory("models")

if(NOT (${DUNE_COMMON_VERSION} VERSION_LESS 2.6))
  dune_python_install_package(PATH "python")
endif()

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
