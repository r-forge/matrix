#ifndef MATRIX_BCROSSTAB_H
#define MATRIX_BCROSSTAB_H

#include "Mutils.h"
#include "cscMatrix.h"
#include "Metis_utils.h"

SEXP bCrosstab_project(SEXP ctab, SEXP j);
SEXP bCrosstab_permute(SEXP ctab, SEXP i, SEXP perm);

#endif
