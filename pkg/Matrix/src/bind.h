#ifndef MATRIX_BIND_H
#define MATRIX_BIND_H

#include "Mutils.h"

/* defined in ./sparse.c : */
SEXP Tsparse_aggregate(SEXP);

SEXP R_bind(SEXP args);

#endif
