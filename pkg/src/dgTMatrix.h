#ifndef MATRIX_TRIPLET_H
#define MATRIX_TRIPLET_H

#include "Mutils.h"

SEXP dgTMatrix_validate(SEXP x);
SEXP dgTMatrix_to_dgeMatrix(SEXP x);
SEXP dgTMatrix_to_matrix(SEXP x);

#endif
