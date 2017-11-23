PDF = dune-fempy.pdf
PY = dune-fempy.py
TEX = dune-fempy.tex
AUX = dune-fempy.aux

#FIGURES =

.PHONY: all
all: $(PDF) $(PY)

.PHONY: clean distclean
clean:
	@rm -f *.aux *.blg *.fdb_latexmk *.fls *.log *.out *.png
distclean: clean
	@rm -f *.bbl $(PDF) $(PY) $(TEX) *.vtu

.PRECIOUS: %.tex

%.tex: %.pmd
	@pweave -f tex $<

%.py: %.pmd
	@ptangle $<

%.pdf: %.tex
	@pdflatex --interaction=nonstopmode $<
	@bibtex $(AUX)
	@pdflatex --interaction=nonstopmode $<
	@pdflatex --interaction=nonstopmode $<
