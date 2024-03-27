#ifndef MATRIX_SOLVE_H
#define MATRIX_SOLVE_H

#include <Rinternals.h>

SEXP           denseLU_solve(SEXP, SEXP);
SEXP denseBunchKaufman_solve(SEXP, SEXP);
SEXP     denseCholesky_solve(SEXP, SEXP);
SEXP          trMatrix_solve(SEXP, SEXP);
SEXP          sparseLU_solve(SEXP, SEXP, SEXP);
SEXP    sparseCholesky_solve(SEXP, SEXP, SEXP, SEXP);
SEXP          tCMatrix_solve(SEXP, SEXP, SEXP);

SEXP sparseQR_matmult(SEXP, SEXP, SEXP, SEXP, SEXP);

#endif /* MATRIX_SOLVE_H */
