#ifndef MATRIX_CHOLMOD_ETC_H
#define MATRIX_CHOLMOD_ETC_H

#include "SuiteSparse/CHOLMOD/Include/cholmod.h"
#include <Rinternals.h>

extern cholmod_common c ;
extern cholmod_common cl;

int Matrix_cholmod_start (cholmod_common *);
int Matrix_cholmod_finish(cholmod_common *);

cholmod_factor *M2CHF(SEXP,  int);
cholmod_sparse *M2CHS(SEXP,  int);
cholmod_dense  *M2CHD(SEXP, char);

SEXP CHF2M(cholmod_factor *,  int);
SEXP CHS2M(cholmod_sparse *,  int, char);
SEXP CHD2M(cholmod_dense  *, char, char);

#endif /* MATRIX_CHOLMOD_ETC_H */
