include(DuneSphinxDoc)
set(SPHINXDIR ${CMAKE_BINARY_DIR}/doc/sphinx)
file(COPY unitcube-2d.dgf DESTINATION ${SPHINXDIR})
file(COPY sphere.dgf DESTINATION ${SPHINXDIR})
file(COPY parameter DESTINATION ${SPHINXDIR})
file(COPY mcf.gif DESTINATION ${SPHINXDIR})

add_sphinx_files(notebooks
  create.ipynb
  crystal.ipynb
  extending.ipynb
  grid-construction.ipynb
  gridfunctions.ipynb
  laplace-intro.ipynb
  laplace-dirichlet.ipynb
  laplace-coefficients.ipynb
  laplace-la.ipynb
  laplace-adaptive.ipynb
  laplace-dg.ipynb
  mcf.ipynb
  parameters.ipynb
  spiral.ipynb
  uflbindings.ipynb
)
