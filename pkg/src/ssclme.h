#ifndef MATRIX_SSCLME_H
#define MATRIX_SSCLME_H

#include "sscCrosstab.h"
#include "ldl.h"
#include "Metis_utils.h"
#include <R_ext/Lapack.h>
#include <R_ext/Constants.h>

SEXP ctab_permute(SEXP ctab);
SEXP ssclme_create(SEXP facs, SEXP ncv, SEXP threshhold);
SEXP ssclme_transfer_dimnames(SEXP x, SEXP facs, SEXP mmats);
SEXP ssclme_update_mm(SEXP x, SEXP facs, SEXP mmats);
SEXP ssclme_inflate_and_factor(SEXP x);
SEXP ssclme_factor(SEXP x);
SEXP ssclme_invert(SEXP x);
SEXP ssclme_initial(SEXP x);
SEXP ssclme_fixef(SEXP x);
SEXP ssclme_ranef(SEXP x);
SEXP ssclme_sigma(SEXP x, SEXP REML);
SEXP ssclme_coef(SEXP x);
SEXP ssclme_coefUnc(SEXP x);
SEXP ssclme_coefGetsUnc(SEXP x, SEXP coef);
SEXP ssclme_coefGets(SEXP x, SEXP coef);
SEXP ssclme_EMsteps(SEXP x, SEXP nsteps, SEXP REMLp, SEXP verb);
SEXP ssclme_fitted(SEXP x, SEXP facs, SEXP mmats);
SEXP ssclme_variances(SEXP x, SEXP REML);
SEXP ssclme_gradient(SEXP x, SEXP REMLp, SEXP Uncp);

#endif
