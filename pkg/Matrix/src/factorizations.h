#ifndef MATRIX_FACTORS_H
#define MATRIX_FACTORS_H

#include "cs.h"
#include "Lapack-etc.h"
#include "Mutils.h"
#include "chm_common.h"

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

int dgCMatrix_trf_(cs *A, css **S, csn **N, int doError,
		   int order, double tol);
SEXP dgCMatrix_trf(SEXP obj, SEXP doError, SEXP keepDimNames,
		   SEXP order, SEXP tol);

int dgCMatrix_orf_(cs *A, css **S, csn **N, int doError,
		   int order);
SEXP dgCMatrix_orf(SEXP obj, SEXP doError, SEXP keepDimNames,
		   SEXP order);

int dpCMatrix_trf_(CHM_SP A, CHM_FR *L,
		   int perm, int ldl, int super, double mult);
SEXP dpCMatrix_trf(SEXP obj,
		   SEXP perm, SEXP ldl, SEXP super, SEXP mult);

SEXP denseLU_expand(SEXP obj);
SEXP BunchKaufman_expand(SEXP obj);

SEXP denseLU_determinant(SEXP obj, SEXP logarithm);
SEXP sparseLU_determinant(SEXP obj, SEXP logarithm);
SEXP sparseQR_determinant(SEXP obj, SEXP logarithm);
SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm);

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP LU_validate(SEXP obj);
#endif /* MJ */

/* MJ: no longer needed ... prefer denseLU_expand() */
#if 0
SEXP LU_expand(SEXP x);
#endif /* MJ */

#endif
