dune_fem_add_integrands_models(forchheimer.ufl)

add_executable( forchheimer EXCLUDE_FROM_ALL main.cc forchheimer.hh )
target_compile_definitions( forchheimer PRIVATE YASPGRID;GRIDDIM=21;POLORDER=2 )
add_dune_mpi_flags( forchheimer )
target_link_libraries( forchheimer )

configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(unitcube-2d.dgf ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
