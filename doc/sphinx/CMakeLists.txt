include(DuneSphinxDoc)
configure_file(sphere.dgf sphere.dgf)
configure_file(mcf.gif mcf.gif)
add_sphinx_targets(sphinx
  advanced.rst
  installation.rst
  usage.rst
  index.rst
  laplace.rst
  tutorial.rst
  laplace-intro.ipynb
  laplace-dirichlet.ipynb
  laplace-coefficients.ipynb
  laplace-la.ipynb
  laplace-adaptive.ipynb
  laplace-dg.ipynb
  crystal-jupyter.ipynb
  mcf.ipynb
  spiral.ipynb
)
