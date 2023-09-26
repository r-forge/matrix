#ifndef MATRIX_CHM_COMMON_H
#define MATRIX_CHM_COMMON_H

#include "cholmod-etc.h"

typedef cholmod_common  * CHM_CM;
typedef cholmod_dense   * CHM_DN;
typedef cholmod_factor  * CHM_FR;
typedef cholmod_sparse  * CHM_SP;
typedef cholmod_triplet * CHM_TR;

SEXP get_SuiteSparse_version(void);

int R_cholmod_start(CHM_CM);
SEXP CHM_set_common_env(SEXP);
void CHM_store_common(void);
void CHM_restore_common(void);

CHM_SP as_cholmod_sparse   (CHM_SP, SEXP, Rboolean, Rboolean);
CHM_TR as_cholmod_triplet  (CHM_TR, SEXP, Rboolean);
CHM_DN as_cholmod_dense    (CHM_DN, SEXP);
CHM_FR as_cholmod_factor   (CHM_FR, SEXP);
CHM_FR as_cholmod_factor3  (CHM_FR, SEXP, Rboolean);
CHM_DN numeric_as_chm_dense(CHM_DN, double *, int, int);

SEXP chm_factor_to_SEXP (CHM_FR, int);
SEXP chm_sparse_to_SEXP (CHM_SP, int, int, int, const char *, SEXP);
SEXP chm_triplet_to_SEXP(CHM_TR, int, int, int, const char *, SEXP);

double chm_factor_ldetL2(CHM_FR);
CHM_FR chm_factor_update(CHM_FR, CHM_SP, double);

#endif /* MATRIX_CHM_COMMON_H */
