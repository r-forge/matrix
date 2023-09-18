#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Mutils.h"

SEXP sparse_drop0(SEXP from, const char *class, double tol);
SEXP R_sparse_drop0(SEXP from, SEXP tol);

SEXP sparse_diag_U2N(SEXP from, const char *class);
SEXP R_sparse_diag_U2N(SEXP from);

SEXP sparse_diag_N2U(SEXP from, const char *class);
SEXP R_sparse_diag_N2U(SEXP from);

SEXP sparse_band(SEXP from, const char *class, int a, int b);
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);

SEXP sparse_diag_get(SEXP obj, const char *class, int names);
SEXP R_sparse_diag_get(SEXP obj, SEXP names);

SEXP sparse_diag_set(SEXP from, const char *class, SEXP value);
SEXP R_sparse_diag_set(SEXP from, SEXP value);

SEXP sparse_transpose(SEXP from, const char *class, int lazy);
SEXP R_sparse_transpose(SEXP from, SEXP lazy);

SEXP sparse_force_symmetric(SEXP from, const char *class, char ul);
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo);

SEXP sparse_symmpart(SEXP from, const char *class);
SEXP R_sparse_symmpart(SEXP from);

SEXP sparse_skewpart(SEXP from, const char *class);
SEXP R_sparse_skewpart(SEXP from);

int sparse_is_symmetric(SEXP obj, const char *class, int checkDN);
SEXP R_sparse_is_symmetric(SEXP obj, SEXP checkDN);

int sparse_is_triangular(SEXP obj, const char *class, int upper);
SEXP R_sparse_is_triangular(SEXP obj, SEXP upper);

int sparse_is_diagonal(SEXP obj, const char *class);
SEXP R_sparse_is_diagonal(SEXP obj);

SEXP sparse_marginsum(SEXP obj, const char *class, int margin,
                      int narm, int mean, int sparse);
SEXP R_sparse_marginsum(SEXP obj, SEXP margin,
                        SEXP narm, SEXP mean, SEXP sparse);

SEXP sparse_sum(SEXP obj, const char *class, int narm);
SEXP R_sparse_sum(SEXP obj, SEXP narm);

SEXP sparse_prod(SEXP obj, const char *class, int narm);
SEXP R_sparse_prod(SEXP obj, SEXP narm);

SEXP Tsparse_aggregate(SEXP from);

#endif /* MATRIX_SPARSE_H */
