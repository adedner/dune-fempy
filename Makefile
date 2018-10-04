SHELL := /bin/bash
PATH := bin:$(PATH)

PDF = dune-fempy.pdf
PY = dune-fempy.py
TEX = dune-fempy.tex battery.tex laplace-adaptive.tex crystal.tex mcf.tex
TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space

.PHONY: all
all: $(PDF) $(PY) dune-fempy.ipynb

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
	@jupyter nbconvert --to markdown laplace-adaptive.ipynb

laplace-adaptive.tex: laplace-adaptive.md
	@pandoc --listings -f markdown -t latex laplace-adaptive.md -o laplace-adaptive.tex
	@python3 pandoc-formatting.py laplace-adaptive.tex

crystal.md: crystal.ipynb
	@jupyter nbconvert --to markdown crystal.ipynb

crystal.tex: crystal.md
	@pandoc --listings -f markdown -t latex crystal.md -o crystal.tex --biblatex --bibliography=dune-fempy.bib
	@python3 pandoc-formatting.py crystal.tex

mcf.md: mcf.ipynb
	@jupyter nbconvert --to markdown mcf.ipynb

mcf.tex: mcf.md
	@pandoc --listings -f markdown -t latex mcf.md -o mcf.tex
	@python3 pandoc-formatting.py mcf.tex

battery.md: battery.ipynb
	@jupyter nbconvert --to markdown battery.ipynb

battery.tex: battery.md
	@pandoc --listings -f markdown -t latex battery.md -o battery.tex --biblatex --bibliography=dune-fempy.bib
	@python3 pandoc-formatting.py battery.tex

$(PDF): $(TEX) dune-fempy.pmd $(TABLE)
	@pdflatex --interaction=nonstopmode $<
	@bibtex dune-fempy.aux
	@pdflatex --interaction=nonstopmode $<
	@pdflatex --interaction=nonstopmode $<

dune-fempy.ipynb: $(PY)
	@python3 py2ipynb.py $(PY) dune-fempy.ipynb
