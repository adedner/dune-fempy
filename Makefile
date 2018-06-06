PDF = dune-fempy.pdf
PY = dune-fempy.py
TEX = dune-fempy.tex
AUX = dune-fempy.aux

.PHONY: all
all: $(PDF) $(PY)

.PHONY: clean distclean
clean:
	@rm -f *.aux *.blg *.fdb_latexmk *.fls *.log *.out *.png
distclean: clean
	@rm -f *.bbl $(PDF) $(PY) $(TEX) *.vtu

.PRECIOUS: %.tex

%.tex: %.pmd
	@$(PYTHON_ENV) pweave -f texpweave $<

%.py: %.pmd
	@ptangle $<

laplace-la.md: laplace-la.ipynb
	@jupyter nbconvert --to markdown laplace-la.ipynb

laplace-la.tex: laplace-la.md
	@pandoc --listings -f markdown -t latex laplace-la.md -o laplace-la.tex
	@python python-highlight.py laplace-la.tex

crystal.md: crystal.ipynb
	@jupyter nbconvert --to markdown crystal.ipynb

crystal.tex: crystal.md
	@pandoc --listings -f markdown -t latex crystal.md -o crystal.tex --biblatex --bibliography=dune-fempy.bib
	@python python-highlight.py crystal.tex

mcf.md: mcf.ipynb
	@jupyter nbconvert --to markdown mcf.ipynb

mcf.tex: mcf.md
	@pandoc --listings -f markdown -t latex mcf.md -o mcf.tex
	@python python-highlight.py mcf.tex

battery.md: battery.ipynb
	@jupyter nbconvert --to markdown battery.ipynb

battery.tex: battery.md
	@pandoc --listings -f markdown -t latex battery.md -o battery.tex --biblatex --bibliography=dune-fempy.bib
	@python python-highlight.py battery.tex

%.pdf: %.tex laplace-la.tex crystal.tex mcf.tex battery.tex
	@pdflatex --interaction=nonstopmode $<
	@bibtex $(AUX)
	@pdflatex --interaction=nonstopmode $<
	@pdflatex --interaction=nonstopmode $<
