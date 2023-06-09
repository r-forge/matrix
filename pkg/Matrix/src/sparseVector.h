#ifndef MATRIX_SPVECTOR_H
#define MATRIX_SPVECTOR_H

#include "Mutils.h"

/* defined in ./sparse.c : */
SEXP R_sparse_as_general(SEXP);

SEXP v2spV(SEXP from);
SEXP CR2spV(SEXP from);

#endif /* MATRIX_SPVECTOR_H */
