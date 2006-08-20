#ifndef MATRIX_SPARSEQR_H
#define MATRIX_SPARSEQR_H

#include "Mutils.h"
#include "cs_utils.h"

SEXP sparseQR_validate(SEXP x);
SEXP sparseQR_qty(SEXP qr, SEXP y, SEXP classed, SEXP trans);
SEXP sparseQR_coef(SEXP qr, SEXP y, SEXP classed);
SEXP sparseQR_resid_fitted(SEXP qr, SEXP y, SEXP classed, SEXP resid);

#endif
