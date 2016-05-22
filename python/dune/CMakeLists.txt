add_subdirectory(generated)
add_subdirectory(generator)

add_library(common SHARED common.cc)
set_target_properties(common PROPERTIES PREFIX "")

add_library(mpihelper SHARED mpihelper.cc)
set_target_properties(mpihelper PROPERTIES PREFIX "")

add_python_target(grid)
add_python_target(function)

configure_file(__init__.py __init__.py COPYONLY)
