#ifndef MATRIX_BCROSSTAB_H
#define MATRIX_BCROSSTAB_H

#include "Mutils.h"
#include "cscMatrix.h"
#include "Metis_utils.h"
#include "triplet_to_col.h"
#include "ldl.h"
#include "sscMatrix.h"
#include "tscMatrix.h"

SEXP bCrosstab_project(SEXP ctab, SEXP j);
SEXP bCrosstab_convert(SEXP bCtab);

#endif
