#ifndef MATRIX_COERCE_H
#define MATRIX_COERCE_H

#include <Rinternals.h>

SEXP vector_as_dense(SEXP, const char *, char, char, char, int, int, int, SEXP);
SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
SEXP sparse_as_dense(SEXP, const char *, int);
SEXP diagonal_as_dense(SEXP, const char *, char, char, int, char, char);
SEXP index_as_dense(SEXP, const char *, char);

SEXP Vector_as_sparse(SEXP, const char *, char, char, char, int, int, int, SEXP);
SEXP matrix_as_sparse(SEXP, const char *, char, char, char, int);
SEXP dense_as_sparse(SEXP, const char *, char);
SEXP diagonal_as_sparse(SEXP, const char *, char, char, char, char, char);
SEXP index_as_sparse(SEXP, const char *, char, char);

SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP sparse_as_kind(SEXP, const char *, char);
SEXP diagonal_as_kind(SEXP, const char *, char);
SEXP index_as_kind(SEXP, const char *, char);

SEXP dense_as_general(SEXP, const char *, int);
SEXP sparse_as_general(SEXP, const char *);

SEXP dense_as_unpacked(SEXP, const char *);
SEXP dense_as_packed(SEXP, const char *, char, char, char);
SEXP sparse_as_Csparse(SEXP, const char *);
SEXP sparse_as_Rsparse(SEXP, const char *);
SEXP sparse_as_Tsparse(SEXP, const char *);

SEXP vector_as_Vector(SEXP, char);
SEXP sparse_as_Vector(SEXP, const char *);
SEXP diagonal_as_Vector(SEXP, const char *);
SEXP index_as_Vector(SEXP, const char *);

#endif /* MATRIX_COERCE_H */
