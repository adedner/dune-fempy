SHELL := /bin/bash
PATH := bin:$(PATH)

PDF = dune-fempy.pdf
PY = dune-fempy.py
TEX = dune-fempy.tex laplace-adaptive.tex crystal.tex mcf.tex mcf-algorithm.tex
TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space
FIGURES = figures/mcf-comparison.png

.PHONY: all
all: $(TABLES) $(FIGURES) $(PDF) $(PY) dune-fempy.ipynb

.PHONY: clean distclean
clean:
	@rm -f *.aux *.blg *.fdb_latexmk *.fls *.log *.out *.png *.tex *.md *.bbl *.toc
distclean: clean
	@rm -f *.bbl $(PDF) $(PY) $(TEX) *.vtu dune-fempy.py

.PRECIOUS: %.tex

%.tex: %.pmd
	@pweave -f texpweave $<

$(PY): dune-fempy.pmd
	@ptangle $<

laplace-adaptive.md: laplace-adaptive.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to markdown laplace-adaptive.ipynb

laplace-adaptive.tex: laplace-adaptive.md
	@pandoc --listings -f markdown -t latex laplace-adaptive.md -o laplace-adaptive.tex
	@python3 pandoc-formatting.py laplace-adaptive.tex

crystal.md: crystal.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to markdown crystal.ipynb

crystal.tex: crystal.md
	@pandoc --listings -f markdown -t latex crystal.md -o crystal.tex --biblatex --bibliography=dune-fempy.bib
	@python3 pandoc-formatting.py crystal.tex

mcf.md: mcf.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to markdown mcf.ipynb

mcf.tex: mcf.md
	@pandoc --listings -f markdown -t latex mcf.md -o mcf.tex
	@python3 pandoc-formatting.py mcf.tex

mcf-algorithm.md: mcf-algorithm.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to markdown mcf-algorithm.ipynb

mcf-algorithm.tex: mcf-algorithm.md
	@pandoc --listings -f markdown -t latex mcf-algorithm.md -o mcf-algorithm.tex
	@python3 pandoc-formatting.py mcf-algorithm.tex

battery.md: battery.ipynb
	@jupyter nbconvert --ExecutePreprocessor.timeout=-1 --execute --allow-errors --to markdown battery.ipynb

battery.tex: battery.md
	@pandoc --listings -f markdown -t latex battery.md -o battery.tex --biblatex --bibliography=dune-fempy.bib
	@python3 pandoc-formatting.py battery.tex

$(PDF): $(TEX) dune-fempy.pmd $(TABLE) $(FIGURES)
	@pdflatex --interaction=nonstopmode $<
	@bibtex dune-fempy.aux
	@pdflatex --interaction=nonstopmode $<
	@pdflatex --interaction=nonstopmode $<

dune-fempy.ipynb: $(PY)
	@python3 py2ipynb.py $(PY) dune-fempy.ipynb

cpp_time.p: mcf-algorithm.md
python_time.p: mcf-algorithm.md
figures/mcf-comparison.png: cpp_time.p python_time.p
	@python3 mcf-comparison-plot.py
