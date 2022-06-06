#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0);
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0);

SEXP Csparse_drop0(SEXP from);
SEXP Rsparse_drop0(SEXP from);
SEXP Tsparse_drop0(SEXP from);

SEXP tCsparse_as_Rsparse(SEXP from);
SEXP tRsparse_as_Csparse(SEXP from);

SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);

#endif
