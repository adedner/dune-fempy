
foreach(script dune-elliptic-model-compiler dune-integrands-model-compiler)
  configure_file(${script}.in ${script})
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()
