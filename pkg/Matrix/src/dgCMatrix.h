#ifndef MATRIX_DGCMATRIX_H
#define MATRIX_DGCMATRIX_H

#include "Mutils.h"
#include "cs_utils.h"

/* defined in factorizations.c : */
SEXP dgCMatrix_trf(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP dgCMatrix_orf(SEXP, SEXP, SEXP, SEXP);

SEXP compressed_non_0_ij(SEXP x, SEXP colP);
SEXP dgCMatrix_qrsol(SEXP x, SEXP y, SEXP ord);
SEXP dgCMatrix_cholsol(SEXP x, SEXP y);
SEXP dgCMatrix_matrix_solve(SEXP Ap, SEXP bp, SEXP give_sparse);

/* MJ: unused */
#if 0
SEXP dgCMatrix_lusol(SEXP x, SEXP y);
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./factorizations.c */
#if 0
SEXP dgCMatrix_QR(SEXP Ap, SEXP order, SEXP keep_dimnames);

#ifdef Matrix_WithSPQR
SEXP dgCMatrix_SPQR(SEXP Ap, SEXP ordering, SEXP econ, SEXP tol);
#endif

SEXP dgCMatrix_LU(SEXP Ap, SEXP orderp, SEXP tolp,
		  SEXP error_on_sing, SEXP keep_dimnames);
#endif /* MJ */

/* MJ: no longer needed ... prefer CRsparse_(col|row)Sums() */
#if 0
SEXP dgCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
SEXP igCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
SEXP lgCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
SEXP ngCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means);
#endif /* MJ */

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0
SEXP xCMatrix_validate(SEXP x);
SEXP xRMatrix_validate(SEXP x);
#endif /* MJ */

/* MJ: no longer needed ... prefer CRsparse_as_Tsparse() */
#if 0
SEXP compressed_to_TMatrix(SEXP x, SEXP colP);
#endif /* MJ */

/* MJ: no longer needed ... 
   now done via R_sparse_transpose(), tCRsparse_as_RCsparse() */
#if 0
SEXP R_to_CMatrix(SEXP x);
#endif /* MJ */

#endif
