SHELL := /bin/bash
PATH := bin:$(PATH)

TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space
FIGURES = figures/3dexample.png figures/mcf-comparison.png figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png
RST = spiral.rst vemdemo.rst uzawa-scipy.rst laplace-adaptive.rst crystal.rst elasticity.rst mcf.rst mcf-algorithm.rst dune-fempy.rst dune-corepy.rst wave.rst twophaseflow.rst

all: $(RST) $(TABLE) $(FIGURES) index.rst additional.rst installation.rst gettingstarted.rst topics.rst \
	   contributions.rst vemdemo_descr.rst twophaseflow_descr.rst
	@rm -rf html
	@sphinx-build -b html . html

.PHONY: clean distclean
clean:
	@rm -f *.ipynb *.vtu *.pvtu *.p *.png $(TABLES) $(FIGURES) $(RST)
	@rm -rf spiral_files  vemdemo_files battery_files  elasticity_files mcf_files wave_files  uzawa-scipy_files \
         crystal_files  laplace-adaptive_files  mcf-algorithm_files \
				 html
distclean: clean
	@rm -f *.bbl
	@rm -rf html

.SECONDARY:

%_nb.ipynb: %.py
	@python3 py2ipynb.py $< $*_nb.ipynb --image="png"
%.rst: %_nb.ipynb
	@bash nbscript2rst.sh $*


cpp_time.p: mcf-algorithm_nb.ipynb
python_time.p: mcf-algorithm_nb.ipynb

figures/mcf-comparison.png: cpp_time.p python_time.p
	@python3 mcf-comparison-plot.py
figures/3dexample.png: 3dexample.py dune-fempy_nb.ipynb
	@pvpython 3dexample.py
figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png: interpolation.py dune-corepy_nb.ipynb
	@pvpython interpolation.py
