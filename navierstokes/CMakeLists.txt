#add_library(stokes SHARED stokes.cc)
#add_library(burgers SHARED burgers.cc)
#target_link_libraries(stokes ${Boost_LIBRARIES})
#set_target_properties(stokes PROPERTIES PREFIX "")
#target_link_libraries(burgers ${Boost_LIBRARIES})
#set_target_properties(burgers PROPERTIES PREFIX "")
#install(FILES
#  conversion.hh
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/construction)
add_python_target(naviermain)
