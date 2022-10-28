#ifndef MATRIX_TSPARSE_H
#define MATRIX_TSPARSE_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Mutils.h"
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP Tsparse_validate(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer Tsparse_as_CRsparse() */
#if 0
SEXP Tsparse_to_Csparse(SEXP x, SEXP tri);
#endif /* MJ */

/* MJ: unused */
#if 0
SEXP Tsparse_to_tCsparse(SEXP x, SEXP uplo, SEXP diag);
#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_diag_U2N() */
#if 0
SEXP Tsparse_diagU2N(SEXP x);
#endif /* MJ */

#endif
