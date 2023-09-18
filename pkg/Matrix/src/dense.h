#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include "Mutils.h"

SEXP dense_band(SEXP from, const char *class, int a, int b);
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2);

SEXP dense_diag_get(SEXP obj, const char *class, int names);
SEXP R_dense_diag_get(SEXP obj, SEXP names);

SEXP dense_diag_set(SEXP from, const char *class, SEXP value, int new);
SEXP R_dense_diag_set(SEXP from, SEXP value);

SEXP dense_transpose(SEXP from, const char *class);
SEXP R_dense_transpose(SEXP from);

SEXP dense_force_symmetric(SEXP from, const char *class, char ul);
SEXP R_dense_force_symmetric(SEXP from, SEXP uplo);

SEXP dense_symmpart(SEXP from, const char *class);
SEXP R_dense_symmpart(SEXP from);

SEXP dense_skewpart(SEXP from, const char *class);
SEXP R_dense_skewpart(SEXP from);

int dense_is_symmetric(SEXP obj, const char *class, int checkDN);
SEXP R_dense_is_symmetric(SEXP obj, SEXP checkDN);

int dense_is_triangular(SEXP obj, const char *class, int upper);
SEXP R_dense_is_triangular(SEXP obj, SEXP upper);

int dense_is_diagonal(SEXP obj, const char *class);
SEXP R_dense_is_diagonal(SEXP obj);

SEXP dense_marginsum(SEXP obj, const char *class, int margin,
                     int narm, int mean);
SEXP R_dense_marginsum(SEXP obj, SEXP margin,
                       SEXP narm, SEXP mean);

SEXP dense_sum(SEXP obj, const char *class, int narm);
SEXP R_dense_sum(SEXP obj, SEXP narm);

SEXP dense_prod(SEXP obj, const char *class, int narm);
SEXP R_dense_prod(SEXP obj, SEXP narm);

#endif /* MATRIX_DENSE_H */
