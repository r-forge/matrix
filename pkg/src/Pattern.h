#ifndef MATRIX_PATTERN_H
#define MATRIX_PATTERN_H

#include <Rinternals.h>
#include <R.h>
#include "chm_common.h"

cholmod_sparse *factor_to_pattern(SEXP fact);

SEXP factor_prod(SEXP f1, SEXP f2);

#endif
