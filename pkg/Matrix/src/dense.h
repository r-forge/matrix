#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP R_dense_as_kind(SEXP from, SEXP kind);
SEXP R_dense_as_matrix(SEXP from);
SEXP R_geMatrix_as_matrix(SEXP from);
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2);

SEXP lsq_dense_Chol(SEXP X, SEXP y);
SEXP lsq_dense_QR(SEXP X, SEXP y);
SEXP lapack_qr(SEXP Xin, SEXP tl);
SEXP dense_to_Csparse(SEXP x);
SEXP ddense_symmpart(SEXP x);
SEXP ddense_skewpart(SEXP x);

/* MJ: no longer needed ... prefer (un)?packedMatrix_force_symmetric() */
#if 0
SEXP dense_to_symmetric(SEXP x, SEXP uplo, SEXP symm_test);
#endif

/* MJ: no longer needed ... prefer R_dense_band() above */
#if 0
SEXP dense_band(SEXP x, SEXP k1, SEXP k2);
#endif /* MJ */

#endif
