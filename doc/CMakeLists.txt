add_subdirectory(geoview_nb_files)
include(DuneSphinxDoc)
add_sphinx_targets(sphinx)
dune_symlink_to_source_files(FILES
additional.rst index.rst twophaseflow_descr.rst contributions.rst
installation.rst vemdemo_descr.rst gettingstarted.rst  nextsteps.rst topics.rst developers.rst
3dexample.py          limit.py              spiral.py
crystal.py            lineplot.py           svg2pdf.py
dune-corepy.py        mcf-algorithm.py      twophaseflow.py
concepts.py           dune-fempy.py         mcf_cmp_plot.py
elasticity.py         mcf.py                uzawa-scipy.py        othergrids.py
gitlab-formatting.py  pandoc-formatting.py  vemdemo.py            discontinuousgalerkin.py
interpolation.py      wave.py               laplace-adaptive.py   laplace-dwr.py
py2ipynb.py
crystal_nb.ipynb           mcf_nb.ipynb
dune-corepy_nb.ipynb       spiral_nb.ipynb
concepts_nb.ipynb          dune-fempy_nb.ipynb        twophaseflow_nb.ipynb
elasticity_nb.ipynb        uzawa-scipy_nb.ipynb       othergrids_nb.ipynb
laplace-adaptive_nb.ipynb  vemdemo_nb.ipynb           discontinuousgalerkin_nb.ipynb
laplace-dwr_nb.ipynb       lineplot_nb.ipynb          wave_nb.ipynb
mcf-algorithm_nb.ipynb
nbscript2rst.sh laplace-dwr.hh  utility.hh  uzawa.hh
figures _static dune-fempy.bib .special.rst
soap.dgf sphere.dgf unitcube-2d.dgf wave_tank.msh
config.opts
levelgridview_nb.rst
geoview_nb.rst geoview_nb_files/geoview_nb_1_0.png
)


#add_custom_target(compilegrids2 ${PYTHON_EXECUTABLE} compilegrids.py 2)
#add_custom_target(compilegrids3 ${PYTHON_EXECUTABLE} compilegrids.py 3)
