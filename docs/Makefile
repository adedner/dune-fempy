SHELL := /bin/bash
PATH := bin:$(PATH)

TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space
FIGURES = figures/3dexample.png figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png
RST = spiral.rst vemdemo.rst uzawa-scipy.rst laplace-adaptive.rst laplace-dual.rst crystal.rst elasticity.rst mcf.rst mcf-algorithm.rst dune-fempy.rst dune-corepy.rst wave.rst twophaseflow.rst
NBS = spiral_nb.ipynb uzawa-scipy_nb.ipynb elasticity_nb.ipynb wave_nb.ipynb \
      laplace-adaptive_nb.ipynb laplace_dual_nb.ipynb crystal_nb.ipynb mcf_nb.ipynb \
			mcf-algorithm_nb.ipynb lineplot_nb.ipynb \
      vemdemo_nb.ipynb twophaseflow_nb.ipynb \
			dune-fempy_nb.ipynb dune-corepy_nb.ipynb

# .EXPORT_ALL_VARIABLES:

unexport DUNE_LOG_LEVEL

all: $(NBS) $(TABLE) $(FIGURES) index.rst additional.rst installation.rst gettingstarted.rst topics.rst \
	   contributions.rst vemdemo_descr.rst twophaseflow_descr.rst
	@rm -rf html
	@sphinx-build -b html . html

.PHONY: clean distclean linkcheck

linkcheck:
	@sphinx-build -b linkcheck . html

clean:
	@rm -f *.vtu *.pvtu *.p *.png $(TABLES) $(FIGURES) $(RST) $(NBS)
	@rm -rf spiral_files  vemdemo_files battery_files  elasticity_files mcf_files wave_files  uzawa-scipy_files \
         crystal_files  laplace-adaptive_files laplace-dual_files  mcf-algorithm_files \
				 html
distclean: clean
	@rm -f *.bbl
	@rm -rf html

.SECONDARY:

%_nb.ipynb: %.py
	@python3 py2ipynb.py $< $*_nb.ipynb --image="png"
%.rst: %_nb.ipynb
	@bash nbscript2rst.sh $*

figures/3dexample.png: 3dexample.py dune-fempy_nb.ipynb
	@pvpython 3dexample.py
figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png: interpolation.py dune-corepy_nb.ipynb
	@pvpython interpolation.py
