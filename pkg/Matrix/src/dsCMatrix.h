#ifndef MATRIX_SSC_H
#define MATRIX_SSC_H

/* MJ: no longer needed ... nothing below */
#if 0
#include "Mutils.h"
#include "chm_common.h"
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./factorizations.c */
#if 0
SEXP R_chkName_Cholesky(SEXP nm, SEXP perm, SEXP LDL, SEXP super);
SEXP R_chm_factor_name (         SEXP perm, SEXP LDL, SEXP super);
SEXP dsCMatrix_chol(SEXP x, SEXP pivot);
SEXP dsCMatrix_Cholesky(SEXP Ap, SEXP perm, SEXP LDL, SEXP super, SEXP Imult);
SEXP dsCMatrix_LDL_D(SEXP Ap, SEXP permP, SEXP resultKind);
SEXP dsCMatrix_Csparse_solve(SEXP a, SEXP b, SEXP LDL);
SEXP dsCMatrix_matrix_solve (SEXP a, SEXP b, SEXP LDL);
#endif /* MJ */

/* MJ: no longer used ... prefer R_sparse_as_general(), CRsparse_as_Tsparse() */
#if 0
SEXP dsCMatrix_to_dgTMatrix(SEXP x);
#endif /* MJ */

#endif
