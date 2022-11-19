#ifndef MATRIX_SUBSCRIPT_H
#define MATRIX_SUBSCRIPT_H

#include "Mutils.h"

SEXP R_subscript_1ary    (SEXP x, SEXP i);
SEXP R_subscript_1ary_mat(SEXP x, SEXP i);
SEXP R_subscript_2ary    (SEXP x, SEXP i, SEXP j, SEXP drop);

/* TODO: compare with macros in ./Mdefines.h */

#define VALID_NONVIRTUAL						\
/*  0 */ "dgCMatrix", "dgRMatrix", "dgTMatrix", "dgeMatrix",		\
/*  4 */ "dsCMatrix", "dsRMatrix", "dsTMatrix", "dspMatrix", "dsyMatrix", \
/*  9 */ "dtCMatrix", "dtRMatrix", "dtTMatrix", "dtpMatrix", "dtrMatrix", \
/* 14 */ "ddiMatrix",							\
/* 15 */ "lgCMatrix", "lgRMatrix", "lgTMatrix", "lgeMatrix",		\
/* 19 */ "lsCMatrix", "lsRMatrix", "lsTMatrix", "lspMatrix", "lsyMatrix", \
/* 24 */ "ltCMatrix", "ltRMatrix", "ltTMatrix", "ltpMatrix", "ltrMatrix", \
/* 29 */ "ldiMatrix",							\
/* 30 */ "ngCMatrix", "ngRMatrix", "ngTMatrix", "ngeMatrix",		\
/* 34 */ "nsCMatrix", "nsRMatrix", "nsTMatrix", "nspMatrix", "nsyMatrix", \
/* 39 */ "ntCMatrix", "ntRMatrix", "ntTMatrix", "ntpMatrix", "ntrMatrix", \
/* 44 */ /* "ndiMatrix", */						\
/* 44 */ "igCMatrix", "igRMatrix", "igTMatrix", "igeMatrix",		\
/* 48 */ "isCMatrix", "isRMatrix", "isTMatrix", "ispMatrix", "isyMatrix", \
/* 53 */ "itCMatrix", "itRMatrix", "itTMatrix", "itpMatrix", "itrMatrix", \
/* 58 */ "idiMatrix",							\
/* 59 */ "zgCMatrix", "zgRMatrix", "zgTMatrix", "zgeMatrix",		\
/* 63 */ "zsCMatrix", "zsRMatrix", "zsTMatrix", "zspMatrix", "zsyMatrix", \
/* 68 */ "ztCMatrix", "ztRMatrix", "ztTMatrix", "ztpMatrix", "ztrMatrix", \
/* 73 */ "zdiMatrix",							\
/* 74 */ "indMatrix", "pMatrix"

#endif
