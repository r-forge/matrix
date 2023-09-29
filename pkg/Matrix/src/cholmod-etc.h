#ifndef MATRIX_CHOLMOD_ETC_H
#define MATRIX_CHOLMOD_ETC_H

#include <Rinternals.h>
#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"

extern cholmod_common c ;
extern cholmod_common cl;

cholmod_factor *M2CF(SEXP, int);
cholmod_sparse *M2CS(SEXP, int);
cholmod_dense  *M2CD(SEXP, int);

SEXP CF2M(cholmod_factor *, int);
SEXP CS2M(cholmod_sparse *, int, char);
SEXP CD2M(cholmod_dense  *, int, char);

#endif /* MATRIX_CHOLMOD_ETC_H */
