include(DuneSphinxDoc)
configure_file(sphere.dgf sphere.dgf)
configure_file(mcf.gif mcf.gif)
add_sphinx_targets(sphinx
  advanced.rst
  installation.rst
  usage.rst
  index.rst
  tutorial.rst
  crystal-jupyter.ipynb
  laplace-jupyter.ipynb
  mcf.ipynb
)
