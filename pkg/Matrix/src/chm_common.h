#ifndef MATRIX_CHM_COMMON_H
#define MATRIX_CHM_COMMON_H

#include "cholmod-etc.h"

typedef cholmod_common  * CHM_CM;
typedef cholmod_factor  * CHM_FR;
typedef cholmod_sparse  * CHM_SP;
typedef cholmod_triplet * CHM_TR;
typedef cholmod_dense   * CHM_DN;

SEXP get_SuiteSparse_version(void);

int R_cholmod_start(CHM_CM);
SEXP R_cholmod_common_envini(SEXP);
void R_cholmod_common_envset(void);
void R_cholmod_common_envget(void);

CHM_FR sexp_as_cholmod_factor   (CHM_FR, SEXP);
CHM_SP sexp_as_cholmod_sparse   (CHM_SP, SEXP, Rboolean, Rboolean);
CHM_DN sexp_as_cholmod_dense    (CHM_DN, SEXP);

SEXP cholmod_factor_as_sexp (CHM_FR, int);
SEXP cholmod_sparse_as_sexp (CHM_SP, int, int, int, const char *, SEXP);

double chm_factor_ldetL2(CHM_FR);
CHM_FR chm_factor_update(CHM_FR, CHM_SP, double);
CHM_DN numeric_as_chm_dense(CHM_DN, double *, int, int);

#endif /* MATRIX_CHM_COMMON_H */
