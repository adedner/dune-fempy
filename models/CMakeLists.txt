set(uflfiles
  boundary.ufl
  laplacemass.ufl
  meancurvature.ufl
  nonlinear.ufl
  system.ufl
  transport.ufl
  varyingcoeff.ufl
)

dune_fem_add_elliptic_models(${uflfiles})

set(headers)
foreach(uflfile ${uflfiles})
  get_filename_component(base ${uflfile} NAME_WE)
  list(APPEND headers ${base}.hh)
endforeach()

add_custom_target(models DEPENDS ${headers})
