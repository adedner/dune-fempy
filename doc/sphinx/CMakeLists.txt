include(DuneSphinxDoc)
add_sphinx_targets(sphinx)
dune_symlink_to_source_files(FILES
additional.rst index.rst twophaseflow_descr.rst contributions.rst
installation.rst vemdemo_descr.rst gettingstarted.rst  topics.rst
.special.rst
3dexample.py          limit.py              spiral.py
crystal.py            lineplot.py           svg2pdf.py
dune-corepy.py        mcf-algorithm.py      twophaseflow.py
dune-fempy.py         mcf_cmp_plot.py       uzawa.py
elasticity.py         mcf.py                uzawa-scipy.py
gitlab-formatting.py  pandoc-formatting.py  vemdemo.py
interpolation.py      parcrystal.py         wave.py
laplace-adaptive.py   parcrystalTest.py
laplace-dual.py       py2ipynb.py
figures tables generated _static
dune-fempy.bib
sphere.dgf  unitcube-2d.dgf  wave_tank.msh
nbscript2rst.sh
pointfunctional.hh  radius.hh  utility.hh  uzawa.hh
)
