#ifndef MATRIX_FACTOR_H
#define MATRIX_FACTOR_H

#include <Rinternals.h>

SEXP geMatrix_scf(SEXP, SEXP, SEXP);
SEXP syMatrix_scf(SEXP, SEXP, SEXP);
SEXP spMatrix_scf(SEXP, SEXP, SEXP);
SEXP geMatrix_trf(SEXP, SEXP);
SEXP syMatrix_trf(SEXP, SEXP);
SEXP spMatrix_trf(SEXP, SEXP);
SEXP poMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP ppMatrix_trf(SEXP, SEXP);

SEXP dgCMatrix_orf(SEXP, SEXP, SEXP);
SEXP dgCMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP dpCMatrix_trf(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP denseBunchKaufman_expand(SEXP);

SEXP CHMfactor_diag_get(SEXP, SEXP);
SEXP CHMfactor_update(SEXP, SEXP, SEXP);
SEXP CHMfactor_updown(SEXP, SEXP, SEXP);

#endif /* MATRIX_FACTOR_H */
