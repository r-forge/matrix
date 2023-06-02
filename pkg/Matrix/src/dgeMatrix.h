#ifndef MATRIX_GEMATRIX_H
#define MATRIX_GEMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP dgeMatrix_Schur(SEXP x, SEXP vectors, SEXP isDGE);
SEXP dgeMatrix_svd(SEXP x, SEXP nu, SEXP nv);
SEXP dgeMatrix_exp(SEXP x);

/* MJ: no longer needed ... prefer more general unpackedMatrix_diag_[gs]et() */
#if 0
SEXP dgeMatrix_getDiag(SEXP x);
SEXP lgeMatrix_getDiag(SEXP x);
SEXP dgeMatrix_setDiag(SEXP x, SEXP d);
SEXP lgeMatrix_setDiag(SEXP x, SEXP d);
/* was unused, not replaced: */
SEXP dgeMatrix_addDiag(SEXP x, SEXP d);
#endif /* MJ */

/* MJ: no longer needed ... prefer more general R_dense_(col|row)Sums() */
#if 0
SEXP dgeMatrix_colsums(SEXP x, SEXP naRmP, SEXP cols, SEXP mean);
#endif /* MJ */

#endif
