include(DuneSphinxDoc)
add_sphinx_targets(sphinx
                   additional.rst      index.rst             twophaseflow_descr.rst contributions.rst   installation.rst      twophaseflow.rst crystal.rst         laplace-adaptive.rst  uzawa-scipy.rst dune-corepy.rst     mcf-algorithm.rst     vemdemo_descr.rst dune-fempy.rst      mcf.rst               vemdemo.rst elasticity.rst      spiral.rst            wave.rst gettingstarted.rst  topics.rst
                   .special.rst
                  )

dune_symlink_to_source_files(FILES
crystal_files      laplace-adaptive_files  twophaseflow_files
dune-corepy_files  mcf-algorithm_files     uzawa-scipy_files
dune-fempy_files   mcf_files               vemdemo_files
elasticity_files   spiral_files            wave_files
sphere.dgf  unitcube-2d.dgf  wave_tank.msh
nbscript2rst.sh
3dexample.py          laplace-adaptive.py  pandoc-formatting.py  uzawa.py
crystal.py            laplace-dual.py      parcrystal.py         uzawa-scipy.py
dune-corepy.py        limit.py             parcrystalTest.py     vemdemo.py
dune-fempy.py         lineplot.py          py2ipynb.py           wave.py
elasticity.py         mcf-algorithm.py     spiral.py
gitlab-formatting.py  mcf_cmp_plot.py      svg2pdf.py
interpolation.py      mcf.py               twophaseflow.py
crystal_nb.ipynb      laplace-adaptive_nb.ipynb  spiral_nb.ipynb        wave_nb.ipynb
dune-corepy_nb.ipynb  lineplot_nb.ipynb          twophaseflow_nb.ipynb
dune-fempy_nb.ipynb   mcf-algorithm_nb.ipynb     uzawa-scipy_nb.ipynb
elasticity_nb.ipynb   mcf_nb.ipynb               vemdemo_nb.ipynb
figures _static
pointfunctional.hh  radius.hh  utility.hh  uzawa.hh
dune-fempy.bib
)
