#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP R_sparse_drop0(SEXP from);
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0);
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0);
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);
SEXP R_sparse_transpose(SEXP from);
SEXP tCRsparse_as_RCsparse(SEXP from);

#endif
