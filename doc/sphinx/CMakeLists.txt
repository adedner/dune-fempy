include(DuneSphinxDoc)
configure_file(sphere.dgf sphere.dgf)
configure_file(unitcube-2d.dgf unitcube-2d.dgf)
configure_file(parameter parameter)
configure_file(mcf.gif mcf.gif)
add_sphinx_targets(sphinx
  create.ipynb
  crystal.ipynb
  demos.rst
  details.rst
  #duneufl.rst
  extending.ipynb
  grid-construction.ipynb
  gridfunctions.ipynb
  index.rst
  installation.rst
  laplace.rst
  laplace-intro.ipynb
  laplace-dirichlet.ipynb
  laplace-coefficients.ipynb
  laplace-la.ipynb
  laplace-adaptive.ipynb
  laplace-dg.ipynb
  mcf.ipynb
  parameter
  parameters.ipynb
  spiral.ipynb
  #uflbindings.ipynb
  unitcube-2d.dgf
)
