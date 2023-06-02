#ifndef MATRIX_TPMATRIX_H
#define MATRIX_TPMATRIX_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Lapack-etc.h"
#include "Mutils.h"
#endif /* MJ */

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0
SEXP dtpMatrix_getDiag(SEXP x);
SEXP ltpMatrix_getDiag(SEXP x);
SEXP dtpMatrix_setDiag(SEXP x, SEXP d);
SEXP ltpMatrix_setDiag(SEXP x, SEXP d);
/* was unused, not replaced: */
SEXP dtpMatrix_addDiag(SEXP x, SEXP d);
#endif

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0
SEXP dtpMatrix_as_dtrMatrix(SEXP from);
#endif

#endif
