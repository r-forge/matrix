#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include <Rinternals.h>

SEXP sparse_aggregate(SEXP, const char *);
SEXP R_sparse_aggregate(SEXP);

SEXP sparse_drop0(SEXP, const char *, double);
SEXP R_sparse_drop0(SEXP, SEXP);

SEXP sparse_diag_U2N(SEXP, const char *);
SEXP R_sparse_diag_U2N(SEXP);

SEXP sparse_diag_N2U(SEXP, const char *);
SEXP R_sparse_diag_N2U(SEXP);

SEXP sparse_marginsum(SEXP, const char *, int, int, int, int);
SEXP R_sparse_marginsum(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP sparse_sum(SEXP, const char *, int);
SEXP R_sparse_sum(SEXP, SEXP);

SEXP sparse_prod(SEXP, const char *, int);
SEXP R_sparse_prod(SEXP, SEXP);

#endif /* MATRIX_SPARSE_H */
