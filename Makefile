SHELL := /bin/bash
PATH := bin:$(PATH)

PDF = dune-fempy.pdf
PY = dune-fempy.py laplace-adaptive.ipynb crystal.ipynb mcf.ipynb mcf-algorithm.ipynb dune-fempy.ipynb
TEX = laplace-adaptive.tex crystal.tex mcf.tex mcf-algorithm.tex
TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space
FIGURES = figures/mcf-comparison.png

.PHONY: all
all: $(TABLES) $(FIGURES) $(PDF) $(PY)

.PHONY: clean distclean
clean:
	@rm -f *.vtu *.pvtu *.p *.aux *.blg *.fdb_latexmk *.fls *.log *.out *.png *.tex *.md *.bbl *.toc $(TABLES) $(FIGURES) $(PDF) $(PY)
distclean: clean
	@rm -f *.bbl $(PDF) $(PY) $(TEX) *.vtu

.PRECIOUS: %.tex

dune-fempy.tex: dune-fempy.pmd
	@pweave -f texpweave $<

dune-fempy-doc.tex: dune-fempy.pmd
	@pweave -f texpweave -d -o dune-fempy-doc.tex $<

dune-fempy.pdf: $(TEX) dune-fempy.tex dune-fempy.pmd $(TABLE) $(FIGURES)
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

laplace-adaptive.ipynb: laplace-adaptive.py
	@python3 py2ipynb.py laplace-adaptive.py laplace-adaptive.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to notebook --inplace laplace-adaptive.ipynb
laplace-adaptive.md: laplace-adaptive.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --allow-errors --to markdown laplace-adaptive.ipynb
laplace-adaptive.tex: laplace-adaptive.md
	@pandoc --listings -f markdown -t latex laplace-adaptive.md -o laplace-adaptive.tex
	@python3 pandoc-formatting.py laplace-adaptive.tex

crystal.ipynb: crystal.py
	@python3 py2ipynb.py crystal.py crystal.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to notebook --inplace crystal.ipynb
crystal.md: crystal.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --allow-errors --to markdown crystal.ipynb
crystal.tex: crystal.md
	@pandoc --listings -f markdown -t latex crystal.md -o crystal.tex --biblatex --bibliography=dune-fempy.bib
	@python3 pandoc-formatting.py crystal.tex

mcf.ipynb: mcf.py
	@python3 py2ipynb.py mcf.py mcf.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to notebook --inplace mcf.ipynb
mcf.md: mcf.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --allow-errors --to markdown mcf.ipynb
mcf.tex: mcf.md
	@pandoc --listings -f markdown -t latex mcf.md -o mcf.tex
	@python3 pandoc-formatting.py mcf.tex

mcf-algorithm.ipynb: mcf-algorithm.py
	@python3 py2ipynb.py mcf-algorithm.py mcf-algorithm.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to notebook --inplace mcf-algorithm.ipynb
mcf-algorithm.md: mcf-algorithm.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --allow-errors --to markdown mcf-algorithm.ipynb
mcf-algorithm.tex: mcf-algorithm.md
	@pandoc --listings -f markdown -t latex mcf-algorithm.md -o mcf-algorithm.tex
	@python3 pandoc-formatting.py mcf-algorithm.tex

battery.md: battery.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to markdown battery.ipynb
battery.tex: battery.md
	@pandoc --listings -f markdown -t latex battery.md -o battery.tex --biblatex --bibliography=dune-fempy.bib
	@python3 pandoc-formatting.py battery.tex

cpp_time.p: mcf-algorithm.md
python_time.p: mcf-algorithm.md
figures/mcf-comparison.png: cpp_time.p python_time.p
	@python3 mcf-comparison-plot.py
