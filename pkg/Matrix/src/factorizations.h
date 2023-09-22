#ifndef MATRIX_FACTORIZATIONS_H
#define MATRIX_FACTORIZATIONS_H

#include <Rinternals.h>

SEXP dgeMatrix_trf (SEXP, SEXP);
SEXP dsyMatrix_trf (SEXP, SEXP);
SEXP dspMatrix_trf (SEXP, SEXP);
SEXP dpoMatrix_trf (SEXP, SEXP, SEXP, SEXP);
SEXP dppMatrix_trf (SEXP, SEXP);

SEXP dgCMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP dgCMatrix_orf(SEXP, SEXP, SEXP);
SEXP dpCMatrix_trf(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP BunchKaufman_expand(SEXP, SEXP);

SEXP      denseLU_determinant(SEXP, SEXP);
SEXP BunchKaufman_determinant(SEXP, SEXP, SEXP);
SEXP     Cholesky_determinant(SEXP, SEXP, SEXP);
SEXP     sparseLU_determinant(SEXP, SEXP);
SEXP     sparseQR_determinant(SEXP, SEXP);
SEXP    CHMfactor_determinant(SEXP, SEXP, SEXP);

SEXP      denseLU_solve(SEXP, SEXP);
SEXP BunchKaufman_solve(SEXP, SEXP);
SEXP     Cholesky_solve(SEXP, SEXP);
SEXP    dtrMatrix_solve(SEXP, SEXP);
SEXP     sparseLU_solve(SEXP, SEXP, SEXP);
SEXP    CHMfactor_solve(SEXP, SEXP, SEXP, SEXP);
SEXP    dtCMatrix_solve(SEXP, SEXP, SEXP);

SEXP sparseQR_matmult(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP CHMfactor_diag_get(SEXP, SEXP);
SEXP CHMfactor_update(SEXP, SEXP, SEXP);
SEXP CHMfactor_updown(SEXP, SEXP, SEXP);

#endif /* MATRIX_FACTORIZATIONS_H */
