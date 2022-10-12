#ifndef MATRIX_POMATRIX_H
#define MATRIX_POMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"
#include "dsyMatrix.h"

SEXP dpoMatrix_chol(SEXP obj);

SEXP dpoMatrix_rcond(SEXP obj);
SEXP dpoMatrix_solve(SEXP a);
SEXP dpoMatrix_matrix_solve(SEXP a, SEXP b);

#endif
