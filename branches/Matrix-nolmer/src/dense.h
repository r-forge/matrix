#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Rdefines.h"
#include "R_ext/Lapack.h"

SEXP lsq_dense_Chol(SEXP X, SEXP y);
SEXP lsq_dense_QR(SEXP X, SEXP y);
SEXP lapack_qr(SEXP Xin, SEXP tl);
SEXP dense_to_Csparse(SEXP x);
SEXP ddense_band(SEXP x, SEXP k1, SEXP k2);

#endif
