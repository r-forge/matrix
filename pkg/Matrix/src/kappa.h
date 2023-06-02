#ifndef MATRIX_KAPPA_H
#define MATRIX_KAPPA_H

#include "Lapack-etc.h"
#include "Mutils.h"

/* defined in factorizations.c : */
SEXP dgeMatrix_trf_(SEXP, int);
SEXP dsyMatrix_trf_(SEXP, int);
SEXP dspMatrix_trf_(SEXP, int);
SEXP dpoMatrix_trf_(SEXP, int, int, double);
SEXP dppMatrix_trf_(SEXP, int);

SEXP dgeMatrix_norm(SEXP obj, SEXP type);
SEXP dgeMatrix_rcond(SEXP obj, SEXP type);

SEXP dtrMatrix_norm(SEXP obj, SEXP type);
SEXP dtrMatrix_rcond(SEXP obj, SEXP type);

SEXP dtpMatrix_norm(SEXP obj, SEXP type);
SEXP dtpMatrix_rcond(SEXP obj, SEXP type);

SEXP dsyMatrix_norm(SEXP obj, SEXP type);
SEXP dsyMatrix_rcond(SEXP obj);

SEXP dspMatrix_norm(SEXP obj, SEXP type);
SEXP dspMatrix_rcond(SEXP obj);

SEXP dpoMatrix_rcond(SEXP obj);
SEXP dppMatrix_rcond(SEXP obj);

#endif
