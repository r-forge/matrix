#ifndef MATRIX_TRS_H
#define MATRIX_TRS_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Mutils.h"
#endif

/* MJ: no longer needed ... prefer R_sparse_as_dense() */
#if 0
SEXP dsTMatrix_as_dsyMatrix(SEXP x);
SEXP lsTMatrix_as_lsyMatrix(SEXP x);
SEXP dtTMatrix_as_dtrMatrix(SEXP x);
SEXP ltTMatrix_as_ltrMatrix(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_general */
#if 0
SEXP dsTMatrix_as_dgTMatrix(SEXP x);
SEXP lsTMatrix_as_lgTMatrix(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_dense() */
#if 0
SEXP nsTMatrix_as_nsyMatrix(SEXP x);
SEXP ntTMatrix_as_ntrMatrix(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_general */
#if 0
SEXP nsTMatrix_as_ngTMatrix(SEXP x);
#endif /* MJ */

#endif
