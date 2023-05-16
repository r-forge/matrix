#ifndef MATRIX_PPMATRIX_H
#define MATRIX_PPMATRIX_H

#include "dspMatrix.h"

/* defined in factorizations.c : */
SEXP dppMatrix_trf_(SEXP, int);

SEXP dppMatrix_rcond(SEXP obj);

#endif
