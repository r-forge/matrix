#ifndef MATRIX_SUBSCRIPT_H
#define MATRIX_SUBSCRIPT_H

#include <Rinternals.h>

SEXP R_subscript_1ary     (SEXP, SEXP, SEXP);
SEXP R_subscript_1ary_2col(SEXP, SEXP, SEXP);
SEXP R_subscript_2ary     (SEXP, SEXP, SEXP);

#endif /* MATRIX_SUBSCRIPT_H */
