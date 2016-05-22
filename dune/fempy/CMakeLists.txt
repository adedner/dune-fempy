add_subdirectory(function)
add_subdirectory(py)
add_subdirectory(pybind11)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy)
