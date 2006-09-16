#include "Mutils.h"
#include "HBMM.h"
#include "chm_common.h"
#include "CHMfactor.h"
#include "Csparse.h"
#include "Tsparse.h"
#include "dense.h"
#include "dgCMatrix.h"
#include "dgTMatrix.h"
#include "dgeMatrix.h"
#include "dpoMatrix.h"
#include "dppMatrix.h"
#include "dsCMatrix.h"
#include "dsTMatrix.h"
#include "dspMatrix.h"
#include "dsyMatrix.h"
#include "dtCMatrix.h"
#include "dtTMatrix.h"
#include "dtrMatrix.h"
#include "dtpMatrix.h"
#include "factorizations.h"
#include "ldense.h"
#include "lgCMatrix.h"
#include "lgTMatrix.h"
#include "lsCMatrix.h"
#include "ltCMatrix.h"
#include "sparseQR.h"
#include <R_ext/Rdynload.h>

#include "Syms.h"

static R_CallMethodDef CallEntries[] = {
    {"BunchKaufman_validate", (DL_FUNC) &BunchKaufman_validate, 1},
    {"pBunchKaufman_validate", (DL_FUNC) &pBunchKaufman_validate, 1},
    {"CHMfactor_to_sparse", (DL_FUNC) &CHMfactor_to_sparse, 1},
    {"Cholesky_validate", (DL_FUNC) &Cholesky_validate, 1},
    {"Csparse_Csparse_prod", (DL_FUNC) &Csparse_Csparse_prod, 2},
    {"Csparse_band", (DL_FUNC) &Csparse_band, 3},
    {"Csparse_crossprod", (DL_FUNC) &Csparse_crossprod, 3},
    {"Csparse_dense_crossprod", (DL_FUNC) &Csparse_dense_crossprod, 2},
    {"Csparse_dense_prod", (DL_FUNC) &Csparse_dense_prod, 2},
    {"Csparse_diagU2N", (DL_FUNC) &Csparse_diagU2N, 1},
    {"Csparse_horzcat", (DL_FUNC) &Csparse_horzcat, 2},
    {"Csparse_to_Tsparse", (DL_FUNC) &Csparse_to_Tsparse, 2},
    {"Csparse_to_dense", (DL_FUNC) &Csparse_to_dense, 1},
    {"Csparse_to_nz_pattern", (DL_FUNC) &Csparse_to_nz_pattern, 2},
    {"Csparse_to_matrix", (DL_FUNC) &Csparse_to_matrix, 1},
    {"Csparse_submatrix", (DL_FUNC) &Csparse_submatrix, 3},
    {"Csparse_symmetric_to_general",
     (DL_FUNC) &Csparse_symmetric_to_general, 1},
    {"Csparse_transpose", (DL_FUNC) &Csparse_transpose, 2},
    {"Csparse_validate", (DL_FUNC) &Csparse_validate, 1},
    {"Csparse_vertcat", (DL_FUNC) &Csparse_vertcat, 2},
    {"pCholesky_validate", (DL_FUNC) &pCholesky_validate, 1},
#ifdef _valid_only_for_old_graph_package
    {"graphNEL_as_dgTMatrix", (DL_FUNC) &graphNEL_as_dgTMatrix, 2},
#endif
    {"LU_expand", (DL_FUNC) &LU_expand, 1},
    {"LU_validate", (DL_FUNC) &LU_validate, 1},
    {"Matrix_expand_pointers", (DL_FUNC) &Matrix_expand_pointers, 1},
    {"Matrix_writeHarwellBoeing", (DL_FUNC) &Matrix_writeHarwellBoeing, 3},
    {"Matrix_writeMatrixMarket", (DL_FUNC) &Matrix_writeMatrixMarket, 3},
    {"SVD_validate", (DL_FUNC) &SVD_validate, 1},
    {"Tsparse_validate", (DL_FUNC) &Tsparse_validate, 1},
    {"Tsparse_to_Csparse", (DL_FUNC) &Tsparse_to_Csparse, 2},
/*     {"csc_check_column_sorting", (DL_FUNC) &csc_check_column_sorting, 1}, */
    {"compressed_to_dgTMatrix", (DL_FUNC) &compressed_to_dgTMatrix, 2},
    {"compressed_non_0_ij", (DL_FUNC) &compressed_non_0_ij, 2},
    {"dense_to_Csparse", (DL_FUNC) &dense_to_Csparse, 1},
    {"dense_nonpacked_validate", (DL_FUNC) &dense_nonpacked_validate, 1},
    {"ddense_band", (DL_FUNC) &ddense_band, 3},
    {"dMatrix_validate", (DL_FUNC) &dMatrix_validate, 1},

    {"dgCMatrix_LU", (DL_FUNC) &dgCMatrix_LU, 3},
    {"dgCMatrix_QR", (DL_FUNC) &dgCMatrix_QR, 2},
    {"dgCMatrix_validate", (DL_FUNC) &dgCMatrix_validate, 1},
    {"dgCMatrix_lusol", (DL_FUNC) &dgCMatrix_lusol, 2},
    {"dgCMatrix_matrix_solve", (DL_FUNC) &dgCMatrix_matrix_solve, 2},
    {"dgCMatrix_qrsol", (DL_FUNC) &dgCMatrix_qrsol, 2},
    {"dgTMatrix_to_dgeMatrix", (DL_FUNC) &dgTMatrix_to_dgeMatrix, 1},
    {"dgTMatrix_to_matrix", (DL_FUNC) &dgTMatrix_to_matrix, 1},
    {"dgTMatrix_validate", (DL_FUNC) &dgTMatrix_validate, 1},
    {"dgeMatrix_LU", (DL_FUNC) &dgeMatrix_LU, 1},
    {"dgeMatrix_Schur", (DL_FUNC) &dgeMatrix_Schur, 2},
    {"dgeMatrix_colsums", (DL_FUNC) &dgeMatrix_colsums, 4},
    {"dgeMatrix_crossprod", (DL_FUNC) &dgeMatrix_crossprod, 2},
    {"dgeMatrix_determinant", (DL_FUNC) &dgeMatrix_determinant, 2},
    {"dgeMatrix_dgeMatrix_crossprod", (DL_FUNC) &dgeMatrix_dgeMatrix_crossprod, 3},
    {"dgeMatrix_matrix_mm", (DL_FUNC) &dgeMatrix_matrix_mm, 3},
    {"dgeMatrix_matrix_solve", (DL_FUNC) &dgeMatrix_matrix_solve, 2},
    {"dgeMatrix_dtpMatrix_mm", (DL_FUNC) &dgeMatrix_dtpMatrix_mm, 2},
    {"dgeMatrix_exp", (DL_FUNC) &dgeMatrix_exp, 1},
    {"dgeMatrix_getDiag", (DL_FUNC) &dgeMatrix_getDiag, 1},
    {"dgeMatrix_matrix_crossprod", (DL_FUNC) &dgeMatrix_matrix_crossprod, 3},
    {"dgeMatrix_norm", (DL_FUNC) &dgeMatrix_norm, 2},
    {"dgeMatrix_rcond", (DL_FUNC) &dgeMatrix_rcond, 2},
    {"dgeMatrix_solve", (DL_FUNC) &dgeMatrix_solve, 1},
    {"dgeMatrix_validate", (DL_FUNC) &dgeMatrix_validate, 1},
    {"dpoMatrix_chol", (DL_FUNC) &dpoMatrix_chol, 1},
    {"dpoMatrix_dgeMatrix_solve", (DL_FUNC) &dpoMatrix_dgeMatrix_solve, 2},
    {"dpoMatrix_matrix_solve", (DL_FUNC) &dpoMatrix_matrix_solve, 2},
    {"dpoMatrix_rcond", (DL_FUNC) &dpoMatrix_rcond, 2},
    {"dpoMatrix_solve", (DL_FUNC) &dpoMatrix_solve, 1},
    {"dpoMatrix_validate", (DL_FUNC) &dpoMatrix_validate, 1},
    {"dppMatrix_chol", (DL_FUNC) &dppMatrix_chol, 1},
    {"dppMatrix_matrix_solve", (DL_FUNC) &dppMatrix_matrix_solve, 2},
    {"dppMatrix_rcond", (DL_FUNC) &dppMatrix_rcond, 2},
    {"dppMatrix_solve", (DL_FUNC) &dppMatrix_solve, 1},
    {"dppMatrix_validate", (DL_FUNC) &dppMatrix_validate, 1},
    {"dsCMatrix_Cholesky", (DL_FUNC) &dsCMatrix_Cholesky, 4},
    {"dsCMatrix_chol", (DL_FUNC) &dsCMatrix_chol, 2},
    {"dsCMatrix_matrix_solve", (DL_FUNC) &dsCMatrix_matrix_solve, 2},
    {"dsCMatrix_to_dgTMatrix", (DL_FUNC) &dsCMatrix_to_dgTMatrix, 1},
    {"dsCMatrix_validate", (DL_FUNC) &dsCMatrix_validate, 1},
    {"dsTMatrix_as_dgTMatrix", (DL_FUNC) &dsTMatrix_as_dgTMatrix, 1},
    {"dsTMatrix_as_dsyMatrix", (DL_FUNC) &dsTMatrix_as_dsyMatrix, 1},
    {"dsTMatrix_validate", (DL_FUNC) &dsTMatrix_validate, 1},
    {"dsyMatrix_as_dspMatrix", (DL_FUNC) &dsyMatrix_as_dspMatrix, 1},
    {"dsyMatrix_as_matrix", (DL_FUNC) &dsyMatrix_as_matrix, 1},
    {"dsyMatrix_matrix_mm", (DL_FUNC) &dsyMatrix_matrix_mm, 3},
    {"dsyMatrix_matrix_solve", (DL_FUNC) &dsyMatrix_matrix_solve, 2},
    {"dsyMatrix_norm", (DL_FUNC) &dsyMatrix_norm, 2},
    {"dsyMatrix_rcond", (DL_FUNC) &dsyMatrix_rcond, 2},
    {"dsyMatrix_solve", (DL_FUNC) &dsyMatrix_solve, 1},
    {"dsyMatrix_validate", (DL_FUNC) &dsyMatrix_validate, 1},
    {"dspMatrix_as_dsyMatrix", (DL_FUNC) &dspMatrix_as_dsyMatrix, 1},
    {"dspMatrix_matrix_mm", (DL_FUNC) &dspMatrix_matrix_mm, 2},
    {"dspMatrix_matrix_solve", (DL_FUNC) &dspMatrix_matrix_solve, 2},
    {"dspMatrix_norm", (DL_FUNC) &dspMatrix_norm, 2},
    {"dspMatrix_rcond", (DL_FUNC) &dspMatrix_rcond, 2},
    {"dspMatrix_solve", (DL_FUNC) &dspMatrix_solve, 1},
    {"dspMatrix_trf", (DL_FUNC) &dspMatrix_trf, 1},
    {"dspMatrix_validate", (DL_FUNC) &dspMatrix_validate, 1},
    {"dtCMatrix_solve", (DL_FUNC) &dtCMatrix_solve, 1},
    {"dtCMatrix_matrix_solve", (DL_FUNC) &dtCMatrix_matrix_solve, 3},
    {"dtCMatrix_upper_solve", (DL_FUNC) &dtCMatrix_upper_solve, 1},
    {"dtCMatrix_validate", (DL_FUNC) &dtCMatrix_validate, 1},
    {"dtTMatrix_as_dtrMatrix", (DL_FUNC) &dtTMatrix_as_dtrMatrix, 1},
    {"dtTMatrix_as_dgCMatrix", (DL_FUNC) &dtTMatrix_as_dgCMatrix, 1},
    {"dtTMatrix_validate", (DL_FUNC) &dtTMatrix_validate, 1},
    {"dtpMatrix_as_dtrMatrix", (DL_FUNC) &dtpMatrix_as_dtrMatrix, 1},
    {"dtpMatrix_getDiag", (DL_FUNC) &dtpMatrix_getDiag, 1},
    {"dtpMatrix_matrix_mm", (DL_FUNC) &dtpMatrix_matrix_mm, 2},
    {"dtpMatrix_matrix_solve", (DL_FUNC) &dtpMatrix_matrix_solve, 2},
    {"dtpMatrix_norm", (DL_FUNC) &dtpMatrix_norm, 2},
    {"dtpMatrix_rcond", (DL_FUNC) &dtpMatrix_rcond, 2},
    {"dtpMatrix_solve", (DL_FUNC) &dtpMatrix_solve, 1},
    {"dtpMatrix_validate", (DL_FUNC) &dtpMatrix_validate, 1},
    {"dtrMatrix_as_dtpMatrix", (DL_FUNC) &dtrMatrix_as_dtpMatrix, 1},
    {"dtrMatrix_as_matrix", (DL_FUNC) &dtrMatrix_as_matrix, 1},
    {"dtrMatrix_matrix_mm", (DL_FUNC) &dtrMatrix_matrix_mm, 3},
    {"dtrMatrix_getDiag", (DL_FUNC) &dtrMatrix_getDiag, 1},
    {"dtrMatrix_matrix_solve", (DL_FUNC) &dtrMatrix_matrix_solve, 2},
    {"dtrMatrix_norm", (DL_FUNC) &dtrMatrix_norm, 2},
    {"dtrMatrix_rcond", (DL_FUNC) &dtrMatrix_rcond, 2},
    {"dtrMatrix_solve", (DL_FUNC) &dtrMatrix_solve, 1},
    {"dtrMatrix_validate", (DL_FUNC) &dtrMatrix_validate, 1},
    {"dup_mMatrix_as_dgeMatrix", (DL_FUNC) &dup_mMatrix_as_dgeMatrix, 1},

    {"lapack_qr", (DL_FUNC) &lapack_qr, 2},

    {"lcsc_to_matrix", (DL_FUNC) &lcsc_to_matrix, 1},
    {"ncsc_to_matrix", (DL_FUNC) &ncsc_to_matrix, 1},
    {"lgCMatrix_validate", (DL_FUNC) &lgCMatrix_validate, 1},
    {"lgTMatrix_validate", (DL_FUNC) &lgTMatrix_validate, 1},

    {"lspMatrix_as_lsyMatrix", (DL_FUNC) &lspMatrix_as_lsyMatrix, 2},
    {"lsyMatrix_as_lspMatrix", (DL_FUNC) &lsyMatrix_as_lspMatrix, 2},
    {"lsyMatrix_as_lgeMatrix", (DL_FUNC) &lsyMatrix_as_lgeMatrix, 2},
    {"ltpMatrix_as_ltrMatrix", (DL_FUNC) &ltpMatrix_as_ltrMatrix, 2},
    {"ltrMatrix_as_lgeMatrix", (DL_FUNC) &ltrMatrix_as_lgeMatrix, 2},
    {"ltrMatrix_as_ltpMatrix", (DL_FUNC) &ltrMatrix_as_ltpMatrix, 2},

    {"lsCMatrix_validate", (DL_FUNC) &lsCMatrix_validate, 1},
    {"ltCMatrix_validate", (DL_FUNC) &ltCMatrix_validate, 1},
    {"lsq_dense_Chol", (DL_FUNC) &lsq_dense_Chol, 2},
    {"lsq_dense_QR", (DL_FUNC) &lsq_dense_QR, 2},
    {"sparseQR_validate", (DL_FUNC) &sparseQR_validate, 1},
    {"sparseQR_qty", (DL_FUNC) &sparseQR_qty, 3},
    {"sparseQR_coef", (DL_FUNC) &sparseQR_coef, 2},
    {"sparseQR_resid_fitted", (DL_FUNC) &sparseQR_resid_fitted, 3},
/*     {"tsc_to_dgTMatrix", (DL_FUNC) &tsc_to_dgTMatrix, 1}, */
    {"triangularMatrix_validate", (DL_FUNC) &triangularMatrix_validate, 1},
    {"symmetricMatrix_validate", (DL_FUNC) &symmetricMatrix_validate, 1},
    {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_Matrix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    R_RegisterCCallable("Matrix", "as_cholmod_dense", (DL_FUNC)as_cholmod_dense);
    R_RegisterCCallable("Matrix", "as_cholmod_factor", (DL_FUNC)as_cholmod_factor);
    R_RegisterCCallable("Matrix", "as_cholmod_sparse", (DL_FUNC)as_cholmod_sparse);
    R_RegisterCCallable("Matrix", "chm_factor_to_SEXP", (DL_FUNC)chm_factor_to_SEXP);

    R_RegisterCCallable("Matrix", "cholmod_aat", (DL_FUNC)cholmod_aat);
    R_RegisterCCallable("Matrix", "cholmod_add", (DL_FUNC)cholmod_add);
    R_RegisterCCallable("Matrix", "cholmod_allocate_dense", (DL_FUNC)cholmod_allocate_dense);
    R_RegisterCCallable("Matrix", "cholmod_allocate_sparse", (DL_FUNC)cholmod_allocate_sparse);
    R_RegisterCCallable("Matrix", "cholmod_analyze", (DL_FUNC)cholmod_analyze);
    R_RegisterCCallable("Matrix", "cholmod_copy", (DL_FUNC)cholmod_copy);
    R_RegisterCCallable("Matrix", "cholmod_copy_dense", (DL_FUNC)cholmod_copy_dense);
    R_RegisterCCallable("Matrix", "cholmod_copy_factor", (DL_FUNC)cholmod_copy_factor);
    R_RegisterCCallable("Matrix", "cholmod_copy_sparse", (DL_FUNC)cholmod_copy_sparse);
    R_RegisterCCallable("Matrix", "cholmod_factor_to_sparse", (DL_FUNC)cholmod_factor_to_sparse);
    R_RegisterCCallable("Matrix", "cholmod_factorize", (DL_FUNC)cholmod_factorize);
    R_RegisterCCallable("Matrix", "cholmod_finish", (DL_FUNC)cholmod_finish);
    R_RegisterCCallable("Matrix", "cholmod_free_dense", (DL_FUNC)cholmod_free_dense);
    R_RegisterCCallable("Matrix", "cholmod_free_factor", (DL_FUNC)cholmod_free_sparse);
    R_RegisterCCallable("Matrix", "cholmod_free_sparse", (DL_FUNC)cholmod_free_sparse);
    R_RegisterCCallable("Matrix", "cholmod_nnz", (DL_FUNC)cholmod_nnz);
    R_RegisterCCallable("Matrix", "cholmod_sdmult", (DL_FUNC)cholmod_sdmult);
    R_RegisterCCallable("Matrix", "cholmod_solve", (DL_FUNC)cholmod_solve);
    R_RegisterCCallable("Matrix", "cholmod_speye", (DL_FUNC)cholmod_speye);
    R_RegisterCCallable("Matrix", "cholmod_spsolve", (DL_FUNC)cholmod_spsolve);
    R_RegisterCCallable("Matrix", "cholmod_start", (DL_FUNC)cholmod_start);
    R_RegisterCCallable("Matrix", "cholmod_transpose", (DL_FUNC)cholmod_transpose);

    R_RegisterCCallable("Matrix", "dpoMatrix_chol", (DL_FUNC)dpoMatrix_chol);
    R_RegisterCCallable("Matrix", "numeric_as_chm_dense", (DL_FUNC)numeric_as_chm_dense);

    cholmod_start(&c);

    Matrix_DimNamesSym = install("Dimnames");
    Matrix_DimSym = install("Dim");
    Matrix_diagSym = install("diag");
    Matrix_factorSym = install("factors");
    Matrix_iSym = install("i");
    Matrix_jSym = install("j");
    Matrix_pSym = install("p");
    Matrix_permSym = install("perm");
    Matrix_uploSym = install("uplo");
    Matrix_xSym = install("x");
}

void R_unload_Matrix(DllInfo *dll)
{
    cholmod_finish(&c);
}
