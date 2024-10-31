#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include <Rinternals.h>

SEXP dense_marginsum(SEXP, const char *, int, int, int);
SEXP R_dense_marginsum(SEXP, SEXP, SEXP, SEXP);

SEXP dense_sum(SEXP, const char *, int);
SEXP R_dense_sum(SEXP, SEXP);

SEXP dense_prod(SEXP, const char *, int);
SEXP R_dense_prod(SEXP, SEXP);

#endif /* MATRIX_DENSE_H */
