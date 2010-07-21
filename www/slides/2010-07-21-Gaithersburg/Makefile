R_PROG = R
RNWSRC = sparse-intro-10 sparse-model-matrices
DRIVERS= MaechlerBates
TEXSRC =
SWEAVE = SweaveSlides.sty # SweaveSlides.cfg

XtrTARGETS = $(DRIVERS:=-4.pdf) $(DRIVERS:=-4a4.pdf) $(RNWSRC:_R.txt)

## 'make' should make the  slides file only
slides: $(DRIVERS:=D.pdf)

all: slides $(XtrTARGETS)

allR: $(RNWSRC:=.R)
	touch allR
allTeX: $(RNWSRC:=.tex)
	touch allTeX
%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex
	${R_PROG} CMD texi2pdf $<

dist: clean
## ^ currently not more
clean:
	$(RM) $(RNWSRC:=.tex) \
          Rplots.ps  Rplots.pdf $(DRIVERS:=*.out) allR figs/*.pdf figs/*.eps \
          *.aux *.log *.toc *.nav *.snm *.vrb *.dvi

realclean: clean
	$(RM) $(DRIVERS:=D.pdf) $(XtrTARGETS)

%D.tex: %.tex allTeX
	@echo "\\documentclass[dvipsnames,pdflatex,beamer]{beamer}" > $@
	@echo "\\input{"$*"}" >> $@
	@echo "\\end{document}" >> $@

%H.tex: %.tex allTeX
	@echo "\\documentclass[dvipsnames,pdflatex,handout]{beamer}" > $@
	@echo "\\input{"$*"}" >> $@
	@echo "\\end{document}" >> $@

%H.pdf: %H.tex
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<
#	texclean $<

%D.pdf: %D.tex
	TEXINPUTS=${TPUTS} texi2dvi --pdf $<
#	texclean $<

%-4.tex: %H.pdf
	@echo "\\documentclass[letterpaper,landscape]{article}" > $@
	@echo "\\usepackage{pdfpages}" >> $@
	@echo "\\begin{document}" >> $@
	@echo "\\includepdf[nup=2x2,pages=-]{"$*"H.pdf}" >> $@
	@echo "\\end{document}" >> $@

%-4.pdf: %-4.tex
	texi2dvi --pdf $<

%-2.tex: %H.pdf
	@echo "\\documentclass[letterpaper]{article}" > $@
	@echo "\\usepackage{pdfpages}" >> $@
	@echo "\\begin{document}" >> $@
	@echo "\\includepdf[nup=1x2,pages=-]{"$*"H.pdf}" >> $@
	@echo "\\end{document}" >> $@

%-2.pdf: %-2.tex
	texi2dvi --pdf $<
#	texclean $<


# A4 paper version: is   echo "-4" | sed s/letter/a4/
%-4a4.tex: %H.pdf
	@echo "\\documentclass[a4paper,landscape]{article}" > $@
	@echo "\\usepackage{pdfpages}" >> $@
	@echo "\\begin{document}" >> $@
	@echo "\\includepdf[nup=2x2,pages=-]{"$*"H.pdf}" >> $@
	@echo "\\end{document}" >> $@

%-4a4.pdf: %-4a4.tex
	texi2dvi --pdf $<
