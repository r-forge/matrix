#ifndef MATRIX_CSPARSE_H
#define MATRIX_CSPARSE_H

#include <Rinternals.h>

SEXP Csparse_validate2(SEXP, SEXP);

SEXP dCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP lCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP iCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP nCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP zCsparse_subassign(SEXP, SEXP, SEXP, SEXP);

SEXP Csparse_MatrixMarket(SEXP, SEXP);
SEXP Csparse_dmperm(SEXP, SEXP, SEXP);

SEXP diag_tC(SEXP, SEXP);

#endif /* MATRIX_CSPARSE_H */
