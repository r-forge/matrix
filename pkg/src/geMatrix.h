#ifndef MATRIX_GEMATRIX_H
#define MATRIX_GEMATRIX_H

#include <R_ext/Lapack.h>
#include "Mutils.h"

SEXP dgeMatrix_validate(SEXP obj);
SEXP dgeMatrix_norm(SEXP obj, SEXP norm);
SEXP dgeMatrix_rcond(SEXP obj, SEXP type);
SEXP dgeMatrix_crossprod(SEXP x);
SEXP dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y);
SEXP dgeMatrix_matrix_crossprod(SEXP x, SEXP y);
SEXP dgeMatrix_getDiag(SEXP x);
SEXP dgeMatrix_LU(SEXP x);
SEXP dgeMatrix_determinant(SEXP x, SEXP logarithm);
SEXP dgeMatrix_solve(SEXP a);
SEXP dgeMatrix_dgeMatrix_mm(SEXP a, SEXP b);
SEXP dgeMatrix_svd(SEXP x, SEXP nu, SEXP nv);
SEXP dgeMatrix_exp(SEXP x);

/* DGESDD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or      */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.                                   */
void F77_NAME(dgesdd)(const char *jobz,
		      const int *m, const int *n,
		      double *a, const int *lda, double *s,
		      double *u, const int *ldu,
		      double *vt, const int *ldvt,
		      double *work, const int *lwork, int *iwork, int *info);


#endif
