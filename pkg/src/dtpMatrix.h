#ifndef MATRIX_TPMATRIX_H
#define MATRIX_TPMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP dtpMatrix_validate(SEXP obj);
SEXP dtpMatrix_norm(SEXP obj, SEXP type);
SEXP dtpMatrix_rcond(SEXP obj, SEXP type);
SEXP dtpMatrix_solve(SEXP a);
SEXP dtpMatrix_matrix_solve(SEXP a, SEXP b);
SEXP dtpMatrix_as_dtrMatrix(SEXP from);
SEXP dtrMatrix_as_dtpMatrix(SEXP from);

#endif /* MATRIX_TPMATRIX_H */
