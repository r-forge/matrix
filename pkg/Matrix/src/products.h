#ifndef MATRIX_PRODUCTS_H
#define MATRIX_PRODUCTS_H

#include "chm_common.h"
#include "Lapack-etc.h"
#include "Mutils.h"

SEXP R_dense_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans);

SEXP R_sparse_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans, SEXP ztrans,
                      SEXP boolean);

SEXP R_diagonal_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans,
                        SEXP boolean);

#endif
