#ifndef MATRIX_PPMATRIX_H
#define MATRIX_PPMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"
#include "dspMatrix.h"

SEXP dppMatrix_chol(SEXP obj);

SEXP dppMatrix_rcond(SEXP obj);
SEXP dppMatrix_solve(SEXP a);
SEXP dppMatrix_matrix_solve(SEXP a, SEXP b);

#endif
