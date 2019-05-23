SHELL := /bin/bash
PATH := bin:$(PATH)

PDF = dune-fempy.pdf
TEX = spiral.tex vemdemo.tex uzawa-scipy.tex laplace-adaptive.tex crystal.tex elasticity.tex mcf.tex mcf-algorithm.tex wave.tex dune-fempy.tex twophaseflow.tex
TABLE = tables/features_discretefunction tables/features_grid tables/features_operator tables/features_solver tables/features_view tables/features_function tables/features_model tables/features_scheme tables/features_space
FIGURES = figures/3dexample.png figures/mcf-comparison.png figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png
GL = spiral_gl.md vemdemo_gl.md uzawa-scipy_gl.md laplace-adaptive_gl.md crystal_gl.md elasticity_gl.md mcf_gl.md mcf-algorithm_gl.md dune-fempy_gl.md wave_gl.md twophaseflow_gl.md
RST = spiral.rst vemdemo.rst uzawa-scipy.rst laplace-adaptive.rst crystal.rst elasticity.rst mcf.rst mcf-algorithm.rst dune-fempy.rst dune-corepy.rst wave.rst twophaseflow.rst

.PHONY: all gitlab tex rst sphinx
all: $(TABLES) $(FIGURES) $(PDF)
tex: $(TEX)
gitlab: $(GL)
rst: $(RST)
sphinx-html: $(RST) dune-fempy.py dune-corepy.py $(TABLE) $(FIGURES) index.rst installation.rst gettingstarted.rst adaptivity.rst moving.rst furtherexamples.rst furtherprojects.rst
	@sphinx-build -b html . html

.PHONY: clean distclean
clean:
	@rm -f *.ipynb $(RST) *.vtu *.pvtu *.p *.aux *.blg *.fdb_latexmk *.fls *.log *.out *.png *.tex *.md *.bbl *.toc $(TABLES) $(FIGURES)
	@rm -rf spiral_files  vemdemo_files battery_files  elasticity_files mcf_files wave_files  uzawa-scipy_files \
         crystal_files  laplace-adaptive_files  mcf-algorithm_files \
				 html
distclean: clean
	@rm -f *.bbl $(PDF) $(TEX)

.PRECIOUS: %.tex

%.tex: %.pmd
	@$(PYTHON_ENV) pweave -f texpweave $<
%.md: %.pmd
	@$(PYTHON_ENV) pweave -f markdown $<
%.py: %.pmd
	@ptangle $<

dune-corepy-doc.tex: dune-corepy.pm
	@pweave -f texpweave -d -o dune-corepy-doc.tex $<
dune-algorithm-doc.tex: dune-algorithm.pmd
	@pweave -f texpweave -d -o dune-algorithm-doc.tex $<
dune-fempy-doc.tex: dune-fempy.pmd
	@pweave -f texpweave -d -o dune-fempy-doc.tex $<

dune-fempy.pdf: $(TEX) dune-fempy.tex dune-corepy.tex dune-algorithm.tex $(TABLE) $(FIGURES)
	@pdflatex --interaction=nonstopmode dune-fempy
	@bibtex dune-fempy.aux
	@pdflatex --interaction=nonstopmode dune-fempy
	@pdflatex --interaction=nonstopmode dune-fempy

dune-fempy-doc.pdf: $(TEX) dune-fempy-doc.tex dune-corepy-doc.tex dune-algorithm-doc.tex $(TABLE) $(FIGURES)
	@pdflatex --interaction=nonstopmode dune-fempy-doc
	@bibtex dune-fempy-doc.aux
	@pdflatex --interaction=nonstopmode dune-fempy-doc
	@pdflatex --interaction=nonstopmode dune-fempy-doc

%.ipynb: %.py
	@python3 py2ipynb.py $< $*_nb.ipynb --image="png"
%.md: %.ipynb
	@jupyter nbconvert --to markdown $<
%.tex: %.ipynb
	@jupyter nbconvert --to markdown $*_nb.ipynb --output $*
	@pandoc --filter svg2pdf.py --listings -f markdown -t latex $*.md -o $@
	@python3 pandoc-formatting.py $@
%.rst: %.ipynb
	@jupyter nbconvert --to rst $*_nb.ipynb --output $*
	@sed -i '/^----------/a :download:`(notebook) <$*_nb.ipynb>` :download:`(script) <$*.py>`' $*.rst
	@sed -i "s/raw:: latex/math::/g" $*.rst
	@sed -i "s/raw-latex/math/g" $*.rst


cpp_time.p: mcf-algorithm_nb.ipynb
python_time.p: mcf-algorithm_nb.ipynb

figures/mcf-comparison.png: cpp_time.p python_time.p
	@python3 mcf-comparison-plot.py
figures/3dexample.png: 3dexample.py $(TEX)
	@pvpython 3dexample.py

figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png: interpolation.py $(TEX)
	@pvpython $<




%_gl.md: %.py
	@python3 py2ipynb.py $< $*_gl.ipynb --image="png"
	@jupyter nbconvert --to markdown $*_gl.ipynb --output $*_gl
	@python3 gitlab-formatting.py $@
