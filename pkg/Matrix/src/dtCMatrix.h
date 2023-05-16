#ifndef MATRIX_TSC_H
#define MATRIX_TSC_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Mutils.h"
#include "dgCMatrix.h"
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP tCMatrix_validate(SEXP x);
SEXP tRMatrix_validate(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./factorizations.c */
#if 0
SEXP dtCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed);
SEXP dtCMatrix_sparse_solve(SEXP a, SEXP b);
#endif /* MJ */

#endif
