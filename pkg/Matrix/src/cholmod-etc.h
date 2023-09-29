#ifndef MATRIX_CHOLMOD_ETC_H
#define MATRIX_CHOLMOD_ETC_H

#include <Rinternals.h>
#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"

extern cholmod_common c ;
extern cholmod_common cl;

cholmod_sparse *M2CS(SEXP obj, int values);
cholmod_dense  *M2CD(SEXP obj, int  trans);
cholmod_factor *M2CF(SEXP obj);

SEXP CS2M(cholmod_sparse *A, int values, char shape);
SEXP CD2M(cholmod_dense  *A, int  trans, char shape);
SEXP CF2M(cholmod_factor *L);

#endif /* MATRIX_CHOLMOD_ETC_H */
