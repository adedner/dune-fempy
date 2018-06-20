dune_fem_add_elliptic_models(forchheimer.ufl)

set(headers)
get_filename_component(base forchheimer.ufl NAME_WE)
list(APPEND headers ${base}.hh)

add_custom_target(forchheimer DEPENDS ${headers})
