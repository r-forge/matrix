#ifndef MATRIX_SYMATRIX_H
#define MATRIX_SYMATRIX_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Lapack-etc.h"
#include "Mutils.h"
#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0
SEXP dsyMatrix_as_dspMatrix(SEXP from);
#endif

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0
SEXP dsyMatrix_as_matrix(SEXP from, SEXP keep_dimnames);
#endif

#endif
