include(DuneSphinxDoc)
dune_sphinx_doc()

add_sphinx_targets(sphinx
  advanced.rst
  installation.rst
  usage.rst
  index.rst
  tutorial.rst
  crystal-jupyter.ipynb
  laplace-jupyter.ipynb
  mcf.ipynb
  sphere.dgf
  mcf.gif
)
