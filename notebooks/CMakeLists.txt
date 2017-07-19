include(DuneSphinxDoc)
set(SPHINXDIR ${CMAKE_BINARY_DIR}/doc/sphinx)
file(COPY unitcube-2d.dgf DESTINATION ${SPHINXDIR})
file(COPY sphere.dgf DESTINATION ${SPHINXDIR})
file(COPY parameter DESTINATION ${SPHINXDIR})
file(COPY mcf.gif DESTINATION ${SPHINXDIR})

set(JUPYTER_FILES
  create.ipynb
  crystal.ipynb
  extending.ipynb
  grid-construction.ipynb
  gridfunctions.ipynb
  laplace-adaptive.ipynb
  laplace-coefficients.ipynb
  laplace-dg.ipynb
  laplace-dirichlet.ipynb
  laplace-intro.ipynb
  laplace-la.ipynb
  mcf.ipynb
  parameters.ipynb
  spiral.ipynb
  uflbindings.ipynb
)

add_sphinx_files(notebooks ${JUPYTER_FILES})

install(FILES ${JUPYTER_FILES} DESTINATION ${CMAKE_INSTALL_FULL_DOCDIR}/jupyter)
