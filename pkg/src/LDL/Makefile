#-------------------------------------------------------------------------------
# Makefile for the LDL mexFunction
#-------------------------------------------------------------------------------

#default: all
default: lib

include ../UFconfig/UFconfig.mk

all: intro libldl.a ldlsimple ldlmain ldlamd

#-------------------------------------------------------------------------------
# the ldl library:
#-------------------------------------------------------------------------------

libldl.a: ldl.c ldl.h
	$(CC) -c -O ldl.c -o ldl.o
	$(AR) libldl.a ldl.o
	- $(RANLIB) libldl.a

lib: libldl.a
	- cp libldl.a ../LDL.a

#-------------------------------------------------------------------------------
# clean-up:
#-------------------------------------------------------------------------------

distclean: purge

purge: clean
	- 'rm' ldl.mex* ldl.dll libldl.a ldlmain ldlamd ldlsimple
	- 'rm' my_ldlmain.out my_ldlamd.out my_ldlsimple.out
	- 'rm' ldlmain.mex* ldlamd.mex* ldlmain.dll ldlamd.dll
	- 'rm' ldlsymbol.mex* ldlsymbol.dll
	- 'rm' *.dvi *.aux *.log *.bak *.bbl *.blg


#-------------------------------------------------------------------------------
# user guide:
#-------------------------------------------------------------------------------

doc: ldl_userguide.tex ldl.bib
	latex ldl_userguide
	- bibtex ldl_userguide
	latex ldl_userguide
	latex ldl_userguide
	dvips ldl_userguide -o ldl_userguide.ps
	pdflatex ldl_userguide
	pdflatex ldl_userguide

