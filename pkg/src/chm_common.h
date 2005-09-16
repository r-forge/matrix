#ifndef CHM_COMMON_H
#define CHM_COMMON_H

#include "CHOLMOD/Include/cholmod.h"
#include <Rinternals.h>

cholmod_common c;

cholmod_sparse *as_cholmod_sparse(SEXP x);

#endif
