R_PROG = R
RNWSRC = BatesMaechler
DRIVERS= $(RNWSRC)
TEXSRC =
SWEAVE = SweaveSlides.sty # SweaveSlides.cfg

## 'make' should make the  slides file only
slides: $(DRIVERS:=D.pdf)

all: slides $(DRIVERS:=-4.pdf) $(DRIVERS:=-4a4.pdf) $(RNWSRC:_R.txt)

allR: $(RNWSRC:=.R)
	touch allR

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex
	${R_PROG} CMD texi2pdf $<

clean:  realclean
realclean: distclean
	$(RM) $(DRIVERS:=D.pdf) $(DRIVERS:=-4.pdf) $(DRIVERS:=-4a4.pdf)

dist: distclean
## ^ currently not more
distclean:
	$(RM) $(RNWSRC:=.tex) \
          Rplots.ps  Rplots.pdf $(DRIVERS:=.out) allR figs/*.pdf figs/*.eps \
          *.aux *.log *.toc *.nav *.snm *.vrb *.dvi

%D.tex: %.tex
	@echo "\\documentclass[dvipsnames,pdflatex,beamer]{beamer}" > $@
	@echo "\\input{"$*"}" >> $@
	@echo "\\end{document}" >> $@

%H.tex: %.tex
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


# A4 paper version
%-4a4.tex: %H.pdf
	@echo "\\documentclass[a4paper,landscape]{article}" > $@
	@echo "\\usepackage{pdfpages}" >> $@
	@echo "\\\\begin{document}" >> $@
	@echo "\\includepdf[nup=2x2,pages=-]{"$*"H.pdf}" >> $@
	@echo "\\end{document}" >> $@

%-4a4.pdf: %-4.tex
	texi2dvi --pdf $<
