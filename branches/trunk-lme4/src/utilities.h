#ifndef NLME_UTILITIES_H
#define NLME_UTILITIES_H
#include <Rdefines.h>

SEXP nlme_replaceSlot(SEXP obj, SEXP names, SEXP value);
SEXP nlme_weight_matrix_list(SEXP MLin, SEXP wts, SEXP adjst, SEXP MLout);

#endif
