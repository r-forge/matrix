#ifndef MATRIX_DENSE_H
#define MATRIX_DENSE_H

#include <Rinternals.h>

SEXP dense_force_symmetric(SEXP, const char *, char, char);
SEXP R_dense_force_symmetric(SEXP, SEXP, SEXP);

SEXP dense_symmpart(SEXP, const char *, char, char);
SEXP R_dense_symmpart(SEXP, SEXP, SEXP);

SEXP dense_skewpart(SEXP, const char *, char);
SEXP R_dense_skewpart(SEXP, SEXP);

SEXP dense_marginsum(SEXP, const char *, int, int, int);
SEXP R_dense_marginsum(SEXP, SEXP, SEXP, SEXP);

SEXP dense_sum(SEXP, const char *, int);
SEXP R_dense_sum(SEXP, SEXP);

SEXP dense_prod(SEXP, const char *, int);
SEXP R_dense_prod(SEXP, SEXP);

#endif /* MATRIX_DENSE_H */
