#ifndef MATRIX_CHOLMOD_ETC_H
#define MATRIX_CHOLMOD_ETC_H

#include <Rinternals.h>
#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"

extern cholmod_common c ;
extern cholmod_common cl;

cholmod_sparse *dgC2cholmod(SEXP obj, int values);
cholmod_dense  *dge2cholmod(SEXP obj, int  trans);
cholmod_factor * mf2cholmod(SEXP obj);

SEXP cholmod2dgC(cholmod_sparse *A, int values, char shape);
SEXP cholmod2dge(cholmod_dense  *A, int  trans, char shape);
SEXP cholmod2mf (cholmod_factor *L);

#endif /* MATRIX_CHOLMOD_ETC_H */
