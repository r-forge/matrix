#ifndef MATRIX_TRMATRIX_H
#define MATRIX_TRMATRIX_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Lapack-etc.h"
#include "Mutils.h"
#endif /* MJ */

/* MJ: no longer needed ... prefer Cholesky_solve() */
#if 0
SEXP dtrMatrix_chol2inv(SEXP a);
#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_diag_[gs]et() */
#if 0
SEXP dtrMatrix_getDiag(SEXP x);
SEXP ltrMatrix_getDiag(SEXP x);
SEXP dtrMatrix_setDiag(SEXP x, SEXP d);
SEXP ltrMatrix_setDiag(SEXP x, SEXP d);
SEXP dtrMatrix_addDiag(SEXP x, SEXP d);
#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0
SEXP dtrMatrix_as_dtpMatrix(SEXP from);
#endif

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0
SEXP dtrMatrix_as_matrix(SEXP from, SEXP keep_dimnames);
#endif

#endif
