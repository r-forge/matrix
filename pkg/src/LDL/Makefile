#-------------------------------------------------------------------------------
# Makefile for the LDL mexFunction
#-------------------------------------------------------------------------------

LIB = ../LDL.a

#default: all
default: lib

include ../UFconfig/UFconfig.mk

all: intro $(LIB) ldlsimple ldlmain ldlamd

#-------------------------------------------------------------------------------
# the ldl library:
#-------------------------------------------------------------------------------

$(LIB): ldl.c ldl.h
	$(CC) -c -O ldl.c -o ldl.o
	$(AR) $(LIB) ldl.o
	- $(RANLIB) $(LIB)
lib: $(LIB)
	@touch lib

#-------------------------------------------------------------------------------
# clean-up:
#-------------------------------------------------------------------------------

distclean: purge
	- 'rm' $(LIB)

purge: clean
	- 'rm' ldl.mex* ldl.dll ldlmain ldlamd ldlsimple
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

