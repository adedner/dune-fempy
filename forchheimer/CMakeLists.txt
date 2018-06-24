dune_fem_add_elliptic_models(forchheimer.ufl)

add_executable( forchheimer_main main.cc forchheimer.hh )
target_compile_definitions( forchheimer_main PRIVATE YASPGRID;GRIDDIM=21;POLORDER=2 )
add_dune_mpi_flags( forchheimer_main )
target_link_libraries( forchheimer_main )

configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(unitcube-2d.dgf ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
