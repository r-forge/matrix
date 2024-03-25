#ifndef MATRIX_KAPPA_H
#define MATRIX_KAPPA_H

#include <Rinternals.h>

SEXP geMatrix_norm(SEXP, SEXP);
SEXP syMatrix_norm(SEXP, SEXP);
SEXP spMatrix_norm(SEXP, SEXP);
SEXP trMatrix_norm(SEXP, SEXP);
SEXP tpMatrix_norm(SEXP, SEXP);

SEXP geMatrix_rcond(SEXP, SEXP, SEXP);
SEXP syMatrix_rcond(SEXP, SEXP, SEXP);
SEXP spMatrix_rcond(SEXP, SEXP, SEXP);
SEXP poMatrix_rcond(SEXP, SEXP, SEXP);
SEXP ppMatrix_rcond(SEXP, SEXP, SEXP);
SEXP trMatrix_rcond(SEXP, SEXP);
SEXP tpMatrix_rcond(SEXP, SEXP);

#endif /* MATRIX_KAPPA_H */
