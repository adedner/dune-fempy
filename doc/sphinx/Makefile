SHELL := /bin/bash
PATH := bin:$(PATH)

FIGURES = generated/figures/3dexample.png generated/figures/interpolation_discrete.png generated/figures/interpolation_exact.png generated/figures/interpolation_error.png
RST = generated/lineplot.rst generated/spiral.rst generated/vemdemo.rst generated/uzawa-scipy.rst generated/laplace-adaptive.rst generated/crystal.rst generated/elasticity.rst generated/mcf.rst generated/mcf-algorithm.rst generated/dune-fempy.rst generated/dune-corepy.rst generated/wave.rst generated/twophaseflow.rst
NBS = generated/spiral_nb.ipynb generated/uzawa-scipy_nb.ipynb generated/elasticity_nb.ipynb generated/wave_nb.ipynb \
      generated/laplace-adaptive_nb.ipynb generated/crystal_nb.ipynb generated/mcf_nb.ipynb \
			generated/mcf-algorithm_nb.ipynb generated/lineplot_nb.ipynb \
      generated/vemdemo_nb.ipynb generated/twophaseflow_nb.ipynb \
			generated/dune-fempy_nb.ipynb generated/dune-corepy_nb.ipynb

# .EXPORT_ALL_VARIABLES:

unexport DUNE_LOG_LEVEL

all: $(RST) $(NBS) $(FIGURES) index.rst additional.rst installation.rst gettingstarted.rst topics.rst \
	   contributions.rst vemdemo_descr.rst twophaseflow_descr.rst
	@rm -rf html
	@sphinx-build -b html . html

.PHONY: clean distclean linkcheck rst

linkcheck:
	@sphinx-build -b linkcheck . html

clean:
	@rm -f *.vtu *.pvtu *.p *.png $(TABLES) $(FIGURES) $(RST) $(NBS)
	@rm -rf spiral_files  vemdemo_files battery_files  elasticity_files mcf_files wave_files  uzawa-scipy_files \
         crystal_files  laplace-adaptive_files mcf-algorithm_files \
				 html
distclean: clean
	@rm -f *.bbl
	@rm -rf html

rst: spiral.rst vemdemo.rst uzawa-scipy.rst laplace-adaptive.rst crystal.rst elasticity.rst mcf.rst mcf-algorithm.rst dune-fempy.rst dune-corepy.rst wave.rst twophaseflow.rst


.SECONDARY:

generated/%_nb.ipynb: %.py
	@python3 py2ipynb.py $< generated/$*_nb.ipynb --image="png"
generated/%.rst: %_nb.ipynb
	@bash nbscript2rst.sh $*

generated/figures/3dexample.png: 3dexample.py dune-fempy_nb.ipynb
	@pvpython 3dexample.py
generated/figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png: interpolation.py dune-corepy_nb.ipynb
	@pvpython interpolation.py
