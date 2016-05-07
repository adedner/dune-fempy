set(DGFFILES
  circle.dgf
  corner.dgf
  sphere.dgf
  unitcube-1d.dgf
  unitcube-2d.dgf
  unitcube-3d.dgf
)

foreach(DGFFILE ${DGFFILES})
  configure_file(${DGFFILE} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()
