set(DGFFILES
  circle.dgf
  corner.dgf
  sphere.dgf
  spiral-2d.dgf
  crystal-2d.dgf
  unitcube-1d.dgf
  unitcube-2d.dgf
  unitcube-3d.dgf
  hole2_larger.dgf
  hole2.dgf
)
# add_gmshgeo_target(karmanvortexstreet)

configure_file(parameter parameter-navier ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

foreach(DGFFILE ${DGFFILES})
  configure_file(${DGFFILE} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()
