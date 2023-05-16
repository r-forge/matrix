#ifndef MATRIX_FACTORS_H
#define MATRIX_FACTORS_H

#include "cs.h"
#include "chm_common.h"
#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dgeMatrix_trf_(SEXP obj,  int warn);
SEXP dsyMatrix_trf_(SEXP obj,  int warn);
SEXP dspMatrix_trf_(SEXP obj,  int warn);
SEXP dpoMatrix_trf_(SEXP obj,  int warn);
SEXP dppMatrix_trf_(SEXP obj,  int warn);

SEXP dgeMatrix_trf (SEXP obj, SEXP warn);
SEXP dsyMatrix_trf (SEXP obj, SEXP warn);
SEXP dspMatrix_trf (SEXP obj, SEXP warn);
SEXP dpoMatrix_trf (SEXP obj, SEXP warn);
SEXP dppMatrix_trf (SEXP obj, SEXP warn);

int dgCMatrix_trf_(cs *A, css **S, csn **N, int order, double tol);
SEXP dgCMatrix_trf(SEXP obj, SEXP order, SEXP tol, SEXP doError);

int dgCMatrix_orf_(cs *A, css **S, csn **N, int order);
SEXP dgCMatrix_orf(SEXP obj, SEXP order, SEXP doError);

int dpCMatrix_trf_(CHM_SP A, CHM_FR *L,
		   int perm, int ldl, int super, double mult);
SEXP dpCMatrix_trf(SEXP obj,
		   SEXP perm, SEXP ldl, SEXP super, SEXP mult);

SEXP denseLU_expand(SEXP obj);
SEXP BunchKaufman_expand(SEXP obj, SEXP packed);

SEXP      denseLU_determinant(SEXP obj, SEXP logarithm);
SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm, SEXP packed);
SEXP     Cholesky_determinant(SEXP obj, SEXP logarithm, SEXP packed);
SEXP     sparseLU_determinant(SEXP obj, SEXP logarithm);
SEXP     sparseQR_determinant(SEXP obj, SEXP logarithm);
SEXP    CHMfactor_determinant(SEXP obj, SEXP logarithm);

SEXP      denseLU_solve(SEXP a, SEXP b);
SEXP BunchKaufman_solve(SEXP a, SEXP b, SEXP packed);
SEXP     Cholesky_solve(SEXP a, SEXP b, SEXP packed);
SEXP     sparseLU_solve(SEXP a, SEXP b, SEXP sparse);
/* MJ: not needed since we have 'sparseQR_coef' : */
#if 0
SEXP     sparseQR_solve(SEXP a, SEXP b, SEXP sparse);
#endif /* MJ */
SEXP    CHMfactor_solve(SEXP a, SEXP b, SEXP sparse, SEXP system);

SEXP    dtrMatrix_solve(SEXP a, SEXP b, SEXP packed);
SEXP    dtCMatrix_solve(SEXP a, SEXP b, SEXP sparse);

#endif
