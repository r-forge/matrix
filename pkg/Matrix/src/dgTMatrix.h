#ifndef MATRIX_TRIPLET_H
#define MATRIX_TRIPLET_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Mutils.h"
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP xTMatrix_validate(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_dense(), R_sparse_as_matrix() */
#if 0
SEXP dgTMatrix_to_dgeMatrix(SEXP x);
SEXP lgTMatrix_to_lgeMatrix(SEXP x);
SEXP dgTMatrix_to_matrix(SEXP x);
SEXP lgTMatrix_to_matrix(SEXP x);
#endif /* MJ */

#endif
