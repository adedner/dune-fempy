SHELL := /bin/bash
PATH := bin:$(PATH)

PDF = dune-fempy.pdf
PY = vemdemo.ipynb uzawa-scipy.ipynb laplace-adaptive.ipynb crystal.ipynb elasticity.ipynb mcf.ipynb mcf-algorithm.ipynb dune-fempy.ipynb wave.ipynb twophaseflow.ipynb dune-fempy.py
TEX = vemdemo.tex uzawa-scipy.tex laplace-adaptive.tex crystal.tex elasticity.tex mcf.tex mcf-algorithm.tex wave.tex dune-fempy.tex twophaseflow.tex
TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space
FIGURES = figures/3dexample.png figures/mcf-comparison.png figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png

.PHONY: all
all: $(TABLES) $(FIGURES) $(PDF) $(PY)

.PHONY: clean distclean
clean:
	@rm -f *.vtu *.pvtu *.p *.aux *.blg *.fdb_latexmk *.fls *.log *.out *.png *.tex *.md *.bbl *.toc $(TABLES) $(FIGURES) $(PY)
	@rm -rf vemdemo_files battery_files  elasticity_files mcf_files wave_files  uzawa-scipy_files \
         crystal_files  laplace-adaptive_files  mcf-algorithm_files
distclean: clean
	@rm -f *.bbl $(PDF) $(PY) $(TEX) *.vtu

.PRECIOUS: %.tex

%.tex: %.pmd
	@$(PYTHON_ENV) pweave -f texpweave $<

dune-fempy-doc.tex: dune-fempy.pmd
	@pweave -f texpweave -d -o dune-fempy-doc.tex $<

dune-fempy.pdf: $(TEX) dune-fempy.pmd $(TABLE) $(FIGURES)
	@pdflatex --interaction=nonstopmode dune-fempy
	@bibtex dune-fempy.aux
	@pdflatex --interaction=nonstopmode dune-fempy
	@pdflatex --interaction=nonstopmode dune-fempy

dune-fempy-doc.pdf: $(TEX) dune-fempy-doc.tex dune-fempy.pmd $(TABLE) $(FIGURES)
	@pdflatex --interaction=nonstopmode dune-fempy-doc
	@bibtex dune-fempy-doc.aux
	@pdflatex --interaction=nonstopmode dune-fempy-doc
	@pdflatex --interaction=nonstopmode dune-fempy-doc

dune-fempy.ipynb: dune-fempy.py
	@python3 py2ipynb.py dune-fempy.py dune-fempy.ipynb

dune-fempy.py: dune-fempy.pmd
	@ptangle $<

%.ipynb: %.py
	@python3 py2ipynb.py $< $@
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --to notebook --inplace $@
%.md: %.ipynb
	@jupyter nbconvert --to markdown $<
%.tex: %.md
	@pandoc --filter svg2pdf.py --listings -f markdown -t latex $< -o $@
	@python3 pandoc-formatting.py $@


cpp_time.p: mcf-algorithm.ipynb
python_time.p: mcf-algorithm.ipynb

figures/mcf-comparison.png: cpp_time.p python_time.p
	@python3 mcf-comparison-plot.py
figures/3dexample.png: 3dexample.py $(TEX)
	@pvpython 3dexample.py

figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png: interpolation.py $(TEX)
	@pvpython $<
