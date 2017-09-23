set(uflfiles
  heat
)

dune_fem_add_elliptic_models(${uflfiles})

set(headers)
foreach(uflfile ${uflfiles})
  get_filename_component(base ${uflfile} NAME_WE)
  list(APPEND headers ${base}.hh)
endforeach()

add_custom_target(compare_models DEPENDS ${headers})

dune_symlink_to_source_files(FILES unitcube.dgf)
dune_symlink_to_source_files(FILES parameter)
dune_symlink_to_source_files(FILES heat.py)

add_executable( heatCPP main.cc)
target_compile_definitions( heatCPP PRIVATE YASPGRID;GRIDDIM=2;WANT_ISTL=1;POLORDER=1 )
add_dune_mpi_flags( heatCPP )
target_link_libraries( heatCPP ${DUNE_LIBS} )
