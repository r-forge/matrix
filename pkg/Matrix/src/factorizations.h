#ifndef MATRIX_FACTORS_H
#define MATRIX_FACTORS_H

#include "Mutils.h"

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP LU_validate(SEXP obj);
#endif /* MJ */

SEXP LU_expand(SEXP x);

#endif
