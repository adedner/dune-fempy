SHELL := /bin/bash
PATH := bin:$(PATH)

FIGURES = figures/3dexample.png figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png
NBS = spiral_nb.ipynb uzawa-scipy_nb.ipynb elasticity_nb.ipynb wave_nb.ipynb \
      laplace-adaptive_nb.ipynb crystal_nb.ipynb \
			mcf_nb.ipynb \
			mcf-algorithm_nb.ipynb lineplot_nb.ipynb laplace-dwr_nb.ipynb \
			dune-fempy_nb.ipynb dune-corepy_nb.ipynb \
      twophaseflow_nb.ipynb # vemdemo_nb.ipynb
			# spiral_adapt_nb.ipynb

# .EXPORT_ALL_VARIABLES:

unexport DUNE_LOG_LEVEL

all: $(NBS) $(FIGURES)
	# @rm -rf html
	# @sphinx-build -b html . html

.PHONY: clean distclean linkcheck

linkcheck:
	@sphinx-build -b linkcheck . html

clean:
	@rm -f *.vtu *.pvtu *.p *.png $(FIGURES) $(NBS)
	@rm -rf *_files
distclean: clean
	@rm -f *.bbl
	@rm -rf html

.SECONDARY:

%_nb.ipynb: %.py
	@python3 py2ipynb.py $< $*_nb.ipynb --image="png"

figures/3dexample.png: 3dexample.py dune-fempy_nb.ipynb
	@pvpython 3dexample.py
figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png: interpolation.py dune-corepy_nb.ipynb
	@pvpython interpolation.py
