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
	@pweave -f texminted $<

%.py: %.pmd
	@ptangle $<

%.pdf: %.tex
	@pdflatex -shell-escape --interaction=nonstopmode $<
	@bibtex $(AUX)
	@pdflatex -shell-escape --interaction=nonstopmode $<
	@pdflatex -shell-escape --interaction=nonstopmode $<
