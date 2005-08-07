#ifndef HBMM_H
#define HBMM_H

#include <Rdefines.h>
#include "Mutils.h"
#include "triplet_to_col.h"

SEXP Matrix_readHarwellBoeing(SEXP filename);
SEXP Matrix_readMatrixMarket(SEXP filename);

#endif
