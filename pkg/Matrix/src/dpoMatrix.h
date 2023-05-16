#ifndef MATRIX_POMATRIX_H
#define MATRIX_POMATRIX_H

#include "dsyMatrix.h"

/* defined in factorizations.c : */
SEXP dpoMatrix_trf_(SEXP, int);

SEXP dpoMatrix_rcond(SEXP obj);

#endif
