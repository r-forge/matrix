#ifndef MATRIX_DETERMINANT_H
#define MATRIX_DETERMINANT_H

#include <Rinternals.h>

SEXP           denseLU_determinant(SEXP, SEXP);
SEXP denseBunchKaufman_determinant(SEXP, SEXP);
SEXP     denseCholesky_determinant(SEXP, SEXP);
SEXP          sparseLU_determinant(SEXP, SEXP);
SEXP          sparseQR_determinant(SEXP, SEXP);
SEXP    sparseCholesky_determinant(SEXP, SEXP, SEXP);

#endif /* MATRIX_DETERMINANT_H */
