#include <Rinternals.h> /* SEXP, Rcomplex */
#include "abIndex.h"
#include "chm_common.h"
#include "CHMfactor.h"
#include "Csparse.h"
#include "dense.h"
#include "dgCMatrix.h"
#include "dgeMatrix.h"
#include "dpoMatrix.h"
#include "dppMatrix.h"
#include "dsCMatrix.h"
#include "dspMatrix.h"
#include "dsyMatrix.h"
#include "dtCMatrix.h"
#include "dtrMatrix.h"
#include "dtpMatrix.h"
#include "factorizations.h"
#include "sparseQR.h"
#include "packedMatrix.h"
#include "unpackedMatrix.h"
#include "sparse.h"
#include "validity.h"
#include <R_ext/Rdynload.h>

#include "Syms.h"
Rcomplex Matrix_zzero, Matrix_zone;

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define EXTDEF(name, n)   {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
    CALLDEF(CHMfactor_to_sparse, 1),
    CALLDEF(CHMfactor_solve, 3),
    CALLDEF(CHMfactor_spsolve, 3),
    CALLDEF(CHMfactor_ldetL2, 1),
    CALLDEF(CHMfactor_ldetL2up, 3),
    CALLDEF(CHMfactor_update, 3),
    CALLDEF(CHMfactor_updown,3),
    CALLDEF(destructive_CHM_update, 3),
    CALLDEF(Csparse_Csparse_prod, 3),
    CALLDEF(Csparse_Csparse_crossprod, 4),
    CALLDEF(Csparse_MatrixMarket, 2),
    
/* MJ: no longer needed ... prefer R_sparse_band() */
/* MJ: however, some reverse dependencies built with Matrix < 1.5-0 need it */
#ifdef Matrix_SupportingCachedMethods
    CALLDEF(Csparse_band, 3),
#endif

    CALLDEF(Csparse_crossprod, 4),
    CALLDEF(Csparse_dense_crossprod, 3),
    CALLDEF(Csparse_dense_prod, 3),
    CALLDEF(Csparse_drop, 2),
    CALLDEF(Csparse_horzcat, 2),
    CALLDEF(Csparse_sort, 1),
    
/* MJ: no longer needed ... prefer CRsparse_to_Tsparse() */
#if 0
    CALLDEF(compressed_to_TMatrix, 2),
    CALLDEF(Csparse_to_Tsparse, 2),
#endif
    
/* MJ: unused */
#if 0
    CALLDEF(Csparse_to_tCsparse, 3),
    CALLDEF(Csparse_to_tTsparse, 3),
#endif

/* MJ: no longer needed ... prefer R_sparse_as_dense() */
#if 0
    CALLDEF(Csparse_to_dense, 2),
#endif

/* MJ: no longer needed ... prefer R_sparse_as_kind() */
#if 0
    CALLDEF(Csparse_to_nz_pattern, 2),
    CALLDEF(nz_pattern_to_Csparse, 2),	
#endif

/* MJ: no longer needed ... prefer R_sparse_as_matrix() */
#if 0
    CALLDEF(Csparse_to_matrix, 3),
#endif

/* MJ: no longer needed ... prefer R_sparse_as_vector() */
#if 0
    CALLDEF(Csparse_to_vector, 1),
#endif
    
    CALLDEF(Csparse_submatrix, 3),
    CALLDEF(dCsparse_subassign, 4),
    CALLDEF(lCsparse_subassign, 4),
    CALLDEF(iCsparse_subassign, 4),
    CALLDEF(nCsparse_subassign, 4),
    CALLDEF(zCsparse_subassign, 4),

/* MJ: no longer needed ... prefer R_sparse_force_symmetric() */
#if 0    
    CALLDEF(Csparse_general_to_symmetric, 3),
#endif

/* MJ: no longer needed ... prefer R_sparse_as_general() */
#if 0
    CALLDEF(Csparse_symmetric_to_general, 1),
#endif    

/* MJ: no longer needed ... prefer R_sparse_transpose() */
/* MJ: however, some reverse dependencies built with Matrix < 1.5-0 need it */
#ifdef Matrix_SupportingCachedMethods
    CALLDEF(Csparse_transpose, 2),
#endif
    
    CALLDEF(Csparse_validate2, 2),
    CALLDEF(Csparse_vertcat, 2),
    CALLDEF(Csparse_dmperm, 3),
    CALLDEF(diag_tC, 2),

/* MJ: no longer needed ... prefer denseLU_expand() */
#if 0
    CALLDEF(LU_expand, 1),
#endif

/* MJ: no longer needed ... prefer R_dense_as_sparse() */
#if 0
    CALLDEF(matrix_to_Csparse, 2),
#endif
    
    CALLDEF(Matrix_expand_pointers, 1),
    CALLDEF(R_rbind2_vector, 2),
    CALLDEF(R_all0, 1),
    CALLDEF(R_any0, 1),

/* MJ: no longer needed ... 
   now done via R_sparse_transpose(), tCRsparse_as_RCsparse() */
#if 0
    CALLDEF(R_to_CMatrix, 1),
#endif

/* MJ: no longer needed ... prefer R_sparse_diag_(U2N|N2U)() */
#if 0
    CALLDEF(Csparse_diagU2N, 1),
    CALLDEF(Csparse_diagN2U, 1),
    CALLDEF(Tsparse_diagU2N, 1),
#endif
    
/* MJ: no longer needed ... prefer Tsparse_as_CRsparse() */
#if 0
    CALLDEF(Tsparse_to_Csparse, 2),
#endif

/* MJ: unused */
#if 0    
    CALLDEF(Tsparse_to_tCsparse, 3),
#endif
    
    CALLDEF(compressed_non_0_ij, 2),

/* MJ: no longer needed ... prefer R_dense_as_sparse() */
#if 0
    CALLDEF(dense_to_Csparse, 1),
#endif
    
/* MJ: no longer needed ... prefer R_dense_band() */
#if 0
    CALLDEF(dense_band, 3),
#endif

/* MJ: no longer needed ... prefer (un)?packedMatrix_force_symmetric() */
#if 0
    CALLDEF(dense_to_symmetric, 3),
#endif

/* MJ: no longer needed ... prefer (un)?packedMatrix_(symm|skew)part() */
#if 0
    CALLDEF(ddense_symmpart, 1),
    CALLDEF(ddense_skewpart, 1),
#endif
    
    CALLDEF(dgCMatrix_LU, 5),
    CALLDEF(dgCMatrix_QR, 3),
#ifdef Matrix_with_SPQR
    CALLDEF(dgCMatrix_SPQR, 4),
#endif

/* MJ: no longer needed ... prefer CRsparse_(col|row)Sums() */
#if 0
    CALLDEF(dgCMatrix_colSums, 5),
    CALLDEF(igCMatrix_colSums, 5),
    CALLDEF(lgCMatrix_colSums, 5),
    CALLDEF(ngCMatrix_colSums, 5),
#endif
    
    CALLDEF(dgCMatrix_cholsol, 2),
    /* CALLDEF(dgCMatrix_lusol, 2), */
    CALLDEF(dgCMatrix_matrix_solve, 3),
    CALLDEF(dgCMatrix_qrsol, 3),

/* MJ: no longer needed ... prefer R_sparse_as_dense() */
#if 0
    CALLDEF(dgTMatrix_to_dgeMatrix, 1),
    CALLDEF(lgTMatrix_to_lgeMatrix, 1),
#endif

/* MJ: no longer needed ... prefer R_sparse_as_matrix() */
#if 0
    CALLDEF(dgTMatrix_to_matrix, 1),
    CALLDEF(lgTMatrix_to_matrix, 1),
#endif

/* MJ: no longer needed ... prefer R_dense_(col|row)Sums() */
#if 0
    CALLDEF(dgeMatrix_colsums, 4),
#endif
    
    CALLDEF(dgeMatrix_trf, 2),
    CALLDEF(dgeMatrix_norm, 2),
    CALLDEF(dgeMatrix_rcond, 2),
    CALLDEF(dgeMatrix_determinant, 2),
    CALLDEF(dgeMatrix_solve, 1),
    CALLDEF(dgeMatrix_matrix_solve, 2),
    CALLDEF(dgeMatrix_crossprod, 2),
    CALLDEF (geMatrix_crossprod, 2),
    CALLDEF(dgeMatrix_dgeMatrix_crossprod, 3),
    CALLDEF (geMatrix_geMatrix_crossprod, 3),
    CALLDEF(dgeMatrix_matrix_crossprod, 3),
    CALLDEF (geMatrix_matrix_crossprod, 3),
    CALLDEF(dgeMatrix_matrix_mm, 3),
    CALLDEF (geMatrix_matrix_mm, 3),
    CALLDEF(dgeMatrix_Schur, 3),
    CALLDEF(dgeMatrix_exp, 1),
    
/* MJ: no longer needed ... prefer unpackedMatrix_diag_[gs]et() */
#if 0
    CALLDEF(dgeMatrix_getDiag, 1),
    CALLDEF(lgeMatrix_getDiag, 1),
    CALLDEF(dgeMatrix_setDiag, 2),
    CALLDEF(lgeMatrix_setDiag, 2),
    /* was unused, not replaced: */
    CALLDEF(dgeMatrix_addDiag, 2),
#endif
    
    CALLDEF(dpoMatrix_trf, 2),
    CALLDEF(dpoMatrix_rcond, 1),
    CALLDEF(dpoMatrix_solve, 1),
    CALLDEF(dpoMatrix_matrix_solve, 2),

    CALLDEF(dppMatrix_trf, 2),
    CALLDEF(dppMatrix_rcond, 1),
    CALLDEF(dppMatrix_solve, 1),
    CALLDEF(dppMatrix_matrix_solve, 2),

    CALLDEF(R_chkName_Cholesky, 4),
    CALLDEF(R_chm_factor_name, 3),
    CALLDEF(dsCMatrix_Cholesky, 5),
    CALLDEF(dsCMatrix_LDL_D, 3),
    CALLDEF(dsCMatrix_chol, 2),
    CALLDEF(dsCMatrix_Csparse_solve, 3),
    CALLDEF(dsCMatrix_matrix_solve,  3),
    
/* MJ: no longer needed ... 
   prefer R_sparse_as_general(), R_sparse_as_dense(), etc. */
#if 0
    CALLDEF(dsCMatrix_to_dgTMatrix, 1),
    CALLDEF(dsTMatrix_as_dgTMatrix, 1),
    CALLDEF(lsTMatrix_as_lgTMatrix, 1),
    CALLDEF(nsTMatrix_as_ngTMatrix, 1),
    CALLDEF(dsTMatrix_as_dsyMatrix, 1),
    CALLDEF(lsTMatrix_as_lsyMatrix, 1),
    CALLDEF(nsTMatrix_as_nsyMatrix, 1),
#endif

/* MJ: no longer needed ... prefer unpackedMatrix_pack() */
#if 0    
    CALLDEF(dsyMatrix_as_dspMatrix, 1),
#endif
    
/* MJ: no longer needed ... prefer R_dense_as_matrix() */
#if 0
    CALLDEF(dsyMatrix_as_matrix, 2),
#endif
    
    CALLDEF(dsyMatrix_trf, 2),
    CALLDEF(dsyMatrix_norm, 2),
    CALLDEF(dsyMatrix_rcond, 1),
    CALLDEF(dsyMatrix_determinant, 2),
    CALLDEF(dsyMatrix_solve, 1),
    CALLDEF(dsyMatrix_matrix_solve, 2),
    CALLDEF(dsyMatrix_matrix_mm, 3),
    
/* MJ: no longer needed ... prefer packedMatrix_unpack() */
#if 0
    CALLDEF(dspMatrix_as_dsyMatrix, 1),
    CALLDEF(dtpMatrix_as_dtrMatrix, 1),
#endif
    
    CALLDEF(dspMatrix_trf, 2),
    CALLDEF(dspMatrix_norm, 2),
    CALLDEF(dspMatrix_rcond, 1),
    CALLDEF(dspMatrix_determinant, 2),
    CALLDEF(dspMatrix_solve, 1),
    CALLDEF(dspMatrix_matrix_solve, 2),
    CALLDEF(dspMatrix_matrix_mm, 2),
    
/* MJ: no longer needed ... prefer packedMatrix_diag_[gs]et() */
#if 0
    CALLDEF(dspMatrix_getDiag, 1),
    CALLDEF(lspMatrix_getDiag, 1),
    CALLDEF(dspMatrix_setDiag, 2),
    CALLDEF(lspMatrix_setDiag, 2),
    /* was unused, not replaced: */
    CALLDEF(dtpMatrix_addDiag, 2),
#endif
    
    /* CALLDEF(dtCMatrix_solve, 1), */
    CALLDEF(dtCMatrix_matrix_solve, 3),
    CALLDEF(dtCMatrix_sparse_solve, 2),

/* MJ: no longer needed ... prefer R_sparse_as_dense() */
#if 0
    CALLDEF(dtTMatrix_as_dtrMatrix, 1),
    CALLDEF(ltTMatrix_as_ltrMatrix, 1),
    CALLDEF(ntTMatrix_as_ntrMatrix, 1),
#endif
    
/* MJ: no longer needed ... prefer packedMatrix_diag_[gs]et() */
#if 0
    CALLDEF(dtpMatrix_getDiag, 1),
    CALLDEF(ltpMatrix_getDiag, 1),
    CALLDEF(dtpMatrix_setDiag, 2),
    CALLDEF(ltpMatrix_setDiag, 2),
#endif
    
/* MJ: no longer needed ... prefer unpackedMatrix_diag_[gs]et() */
#if 0
    CALLDEF(dtrMatrix_getDiag, 1),
    CALLDEF(ltrMatrix_getDiag, 1),
    CALLDEF(dtrMatrix_setDiag, 2),
    CALLDEF(ltrMatrix_setDiag, 2),
#endif
    
    CALLDEF(dtpMatrix_norm, 2),
    CALLDEF(dtpMatrix_rcond, 2),
    CALLDEF(dtpMatrix_solve, 1),
    CALLDEF(dtpMatrix_matrix_solve, 2),
    CALLDEF(dtpMatrix_matrix_mm, 4),
    CALLDEF(dgeMatrix_dtpMatrix_mm, 2),

/* MJ: no longer needed ... prefer unpackedMatrix_pack() */
#if 0
    CALLDEF(dtrMatrix_as_dtpMatrix, 1),
#endif
    
/* MJ: no longer needed ... prefer R_dense_as_matrix() */
#if 0
    CALLDEF(dtrMatrix_as_matrix, 2),
#endif
    
    CALLDEF(dtrMatrix_norm, 2),
    CALLDEF(dtrMatrix_rcond, 2),
    CALLDEF(dtrMatrix_solve, 1),
    CALLDEF(dtrMatrix_matrix_solve, 2),
    CALLDEF(dtrMatrix_dtrMatrix_mm, 4),
    CALLDEF(dtrMatrix_matrix_mm, 4),
    CALLDEF(dtrMatrix_chol2inv, 1),
    CALLDEF(dtrMatrix_addDiag, 2),
    
/* MJ: no longer needed ... prefer R_sparse_as_matrix() */
#if 0
    CALLDEF(lgC_to_matrix, 1),
    CALLDEF(ngC_to_matrix, 1),
#endif
    
/* MJ: no longer needed ... prefer (un)?packedMatrix_(un)?pack() */
#if 0
    CALLDEF(lspMatrix_as_lsyMatrix, 2),
    CALLDEF(lsyMatrix_as_lspMatrix, 2),
    CALLDEF(ltpMatrix_as_ltrMatrix, 2),
    CALLDEF(ltrMatrix_as_ltpMatrix, 2),
#endif

/* MJ: no longer needed ... prefer R_dense_as_general() */
#if 0
    CALLDEF(ltrMatrix_as_lgeMatrix, 2),
    CALLDEF(lsyMatrix_as_lgeMatrix, 2),
#endif

    CALLDEF(lapack_qr, 2),
    CALLDEF(lsq_dense_Chol, 2),
    CALLDEF(lsq_dense_QR, 2),
    CALLDEF(sparseQR_qty, 4),
    CALLDEF(sparseQR_coef, 2),
    CALLDEF(sparseQR_resid_fitted, 3),

    CALLDEF(Matrix_validate, 1),
    CALLDEF(MatrixFactorization_validate, 1),
    CALLDEF(compMatrix_validate, 1),
    CALLDEF(dMatrix_validate, 1),
    CALLDEF(lMatrix_validate, 1),
    CALLDEF(ndenseMatrix_validate, 1),
    CALLDEF(iMatrix_validate, 1),
    CALLDEF(zMatrix_validate, 1),

    CALLDEF(symmetricMatrix_validate, 1),
    CALLDEF(triangularMatrix_validate, 1),
    CALLDEF(diagonalMatrix_validate, 1),
    CALLDEF(indMatrix_validate, 1),
    CALLDEF(pMatrix_validate, 1),
    CALLDEF(CsparseMatrix_validate, 1),
    CALLDEF(RsparseMatrix_validate, 1),
    CALLDEF(TsparseMatrix_validate, 1),

    CALLDEF(sCMatrix_validate, 1),
    CALLDEF(tCMatrix_validate, 1),
    CALLDEF(sRMatrix_validate, 1),
    CALLDEF(tRMatrix_validate, 1),
    CALLDEF(sTMatrix_validate, 1),
    CALLDEF(tTMatrix_validate, 1),
    
    CALLDEF(xgCMatrix_validate, 1),
    CALLDEF(xsCMatrix_validate, 1),
    CALLDEF(xtCMatrix_validate, 1),
    CALLDEF(xgRMatrix_validate, 1),
    CALLDEF(xsRMatrix_validate, 1),
    CALLDEF(xtRMatrix_validate, 1),
    CALLDEF(xgTMatrix_validate, 1),
    CALLDEF(xsTMatrix_validate, 1),
    CALLDEF(xtTMatrix_validate, 1),
    
    CALLDEF(unpackedMatrix_validate, 1),
    CALLDEF(packedMatrix_validate, 1),
    CALLDEF(dpoMatrix_validate, 1),
    CALLDEF(dppMatrix_validate, 1),
    CALLDEF(corMatrix_validate, 1),
    CALLDEF(Cholesky_validate, 1),
    CALLDEF(pCholesky_validate, 1),
    CALLDEF(BunchKaufman_validate, 1),
    CALLDEF(pBunchKaufman_validate, 1),
    CALLDEF(Schur_validate, 1),
    CALLDEF(denseLU_validate, 1),
    CALLDEF(sparseLU_validate, 1),
    CALLDEF(sparseQR_validate, 1),
    CALLDEF(CHMfactor_validate, 1),
    CALLDEF(CHMsimpl_validate, 1),
    CALLDEF(CHMsuper_validate, 1),

/* MJ: some reverse dependencies built with Matrix < 1.5-0 need these */
#ifdef Matrix_SupportingCachedMethods
    { "Csparse_validate", (DL_FUNC) &CsparseMatrix_validate, 1},
    { "Rsparse_validate", (DL_FUNC) &RsparseMatrix_validate, 1},
    { "Tsparse_validate", (DL_FUNC) &TsparseMatrix_validate, 1},
    {"xCMatrix_validate", (DL_FUNC)     &xgCMatrix_validate, 1},
    {"xRMatrix_validate", (DL_FUNC)     &xgRMatrix_validate, 1},
    {"xTMatrix_validate", (DL_FUNC)     &xgTMatrix_validate, 1},
    {      "LU_validate", (DL_FUNC)       &denseLU_validate, 1},
#endif

    CALLDEF(R_Dim_validate, 1),
    CALLDEF(R_DimNames_validate, 2),

/* MJ: some reverse dependencies built with Matrix < 1.5-0 need these */
#ifdef Matrix_SupportingCachedMethods
    {     "Dim_validate", (DL_FUNC)      &R_Dim_validate_old, 2},
    {"dimNames_validate", (DL_FUNC) &R_DimNames_validate_old, 1},
#endif
    
    CALLDEF(R_DimNames_fixup, 1),
    CALLDEF(R_DimNames_is_symmetric, 1),
    CALLDEF(R_symmDN, 1),
    CALLDEF(R_revDN, 1),
    CALLDEF(R_Matrix_kind, 2),
    CALLDEF(R_Matrix_shape, 1),
    CALLDEF(R_Matrix_repr, 1),
    CALLDEF(R_index_triangle, 4),
    CALLDEF(R_index_diagonal, 3),
    CALLDEF(R_nnz, 3),

    CALLDEF(R_sparse_as_dense, 2),
    CALLDEF(R_sparse_as_matrix, 1),
    CALLDEF(R_sparse_as_vector, 1),
    CALLDEF(R_sparse_as_kind, 3),
    CALLDEF(R_sparse_as_general, 1),

    CALLDEF(R_diagonal_as_sparse, 4),
    CALLDEF(R_diagonal_as_dense, 3),
    CALLDEF(R_diagonal_as_kind, 2),

    CALLDEF(R_sparse_drop0, 1),
    CALLDEF(R_sparse_band, 3),
    CALLDEF(R_sparse_diag_get, 2),
    CALLDEF(R_sparse_diag_set, 2),
    CALLDEF(R_sparse_diag_U2N, 1),
    CALLDEF(R_sparse_diag_N2U, 1),
    CALLDEF(R_sparse_transpose, 1),
    CALLDEF(R_sparse_force_symmetric, 2),
    CALLDEF(R_sparse_symmpart, 1),
    CALLDEF(R_sparse_skewpart, 1),

    CALLDEF(CRsparse_as_Tsparse, 1),
    CALLDEF(Tsparse_as_CRsparse, 2),
    CALLDEF(Tsparse_aggregate, 1),
    CALLDEF(tCRsparse_as_RCsparse, 1),
    CALLDEF(Csparse_is_diagonal, 1),
    CALLDEF(Rsparse_is_diagonal, 1),
    CALLDEF(Tsparse_is_diagonal, 1),
    CALLDEF(Csparse_is_triangular, 2),
    CALLDEF(Rsparse_is_triangular, 2),
    CALLDEF(Tsparse_is_triangular, 2),
    CALLDEF(Csparse_is_symmetric, 2),
    CALLDEF(Rsparse_is_symmetric, 2),
    CALLDEF(CRsparse_colSums, 4),
    CALLDEF(CRsparse_rowSums, 4),
    
    CALLDEF(R_matrix_as_dense, 4),
    CALLDEF(R_dense_as_general, 2),
    CALLDEF(R_dense_as_sparse, 4),
    CALLDEF(R_dense_as_kind, 2),
    CALLDEF(R_dense_as_matrix, 2),
    CALLDEF(R_geMatrix_as_matrix, 2),
    CALLDEF(R_dense_as_vector, 2),
    CALLDEF(R_geMatrix_as_vector, 2),
    CALLDEF(R_dense_band, 3),
    CALLDEF(R_dense_colSums, 3),
    CALLDEF(R_dense_rowSums, 3),
    
    CALLDEF(matrix_is_symmetric, 2),
    CALLDEF(matrix_is_triangular, 2),
    CALLDEF(matrix_is_diagonal, 1),
    CALLDEF(matrix_symmpart, 1),
    CALLDEF(matrix_skewpart, 1),
    CALLDEF(matrix_trf, 3),
    
    CALLDEF(unpackedMatrix_pack, 4),
    CALLDEF(unpackedMatrix_force_symmetric, 2),
    CALLDEF(unpackedMatrix_is_symmetric, 2),
    CALLDEF(unpackedMatrix_is_triangular, 2),
    CALLDEF(unpackedMatrix_is_diagonal, 1),
    CALLDEF(unpackedMatrix_transpose, 1),
    CALLDEF(unpackedMatrix_diag_get, 2),
    CALLDEF(unpackedMatrix_diag_set, 2),
    CALLDEF(unpackedMatrix_symmpart, 1),
    CALLDEF(unpackedMatrix_skewpart, 1),
    
    CALLDEF(packedMatrix_unpack, 2),
    CALLDEF(packedMatrix_force_symmetric, 2),
    CALLDEF(packedMatrix_is_symmetric, 2),
    CALLDEF(packedMatrix_is_triangular, 2),
    CALLDEF(packedMatrix_is_diagonal, 1),
    CALLDEF(packedMatrix_transpose, 1),
    CALLDEF(packedMatrix_diag_get, 2),
    CALLDEF(packedMatrix_diag_set, 2),
    CALLDEF(packedMatrix_symmpart, 1),
    CALLDEF(packedMatrix_skewpart, 1),
    CALLDEF(packedMatrix_sub1, 2),
    CALLDEF(packedMatrix_sub1_mat, 2),
    CALLDEF(packedMatrix_sub2, 4),

    CALLDEF(denseLU_expand, 1),
    CALLDEF(BunchKaufman_expand, 1),
    CALLDEF(denseLU_determinant, 2),
    CALLDEF(BunchKaufman_determinant, 2),
    
    CALLDEF(CHM_set_common_env, 1),

    CALLDEF(inv_permutation, 3),
    CALLDEF(m_encodeInd,  4),
    CALLDEF(m_encodeInd2, 5),

    CALLDEF(Matrix_rle_i, 2),
    CALLDEF(Matrix_rle_d, 2),

    CALLDEF(R_set_factor, 4),
    CALLDEF(R_empty_factors, 2),

    CALLDEF(get_SuiteSparse_version, 0),
    {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExtEntries[] = {
    EXTDEF(Mmatrix, 7),
    {NULL, NULL, 0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_Matrix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, ExtEntries);
    R_useDynamicSymbols(dll, FALSE);

/* These are callable from other packages' C code: */

#define RREGDEF(name)  R_RegisterCCallable("Matrix", #name, (DL_FUNC) name)

#if 0
    RREGDEF(Csparse_diagU2N);
    RREGDEF(dpoMatrix_chol);
#else
    R_RegisterCCallable(
	"Matrix", "Csparse_diagU2N", (DL_FUNC) R_sparse_diag_U2N);
    R_RegisterCCallable(
	"Matrix", "dpoMatrix_chol", (DL_FUNC) dpoMatrix_trf);
#endif
    
    RREGDEF(as_cholmod_dense);
    RREGDEF(as_cholmod_factor);
    RREGDEF(as_cholmod_factor3);
    RREGDEF(as_cholmod_sparse);
    RREGDEF(as_cholmod_triplet);
    RREGDEF(chm_factor_to_SEXP);
    RREGDEF(chm_factor_ldetL2);
    RREGDEF(chm_factor_update);
    RREGDEF(chm_sparse_to_SEXP);
    RREGDEF(chm_triplet_to_SEXP);

    RREGDEF(cholmod_aat);
    RREGDEF(cholmod_add);
    RREGDEF(cholmod_allocate_dense);
    RREGDEF(cholmod_allocate_sparse);
    RREGDEF(cholmod_allocate_triplet);
    RREGDEF(cholmod_analyze);
    RREGDEF(cholmod_analyze_p);
    RREGDEF(cholmod_band_inplace);
    RREGDEF(cholmod_change_factor);
    RREGDEF(cholmod_copy);
    RREGDEF(cholmod_copy_dense);
    RREGDEF(cholmod_copy_factor);
    RREGDEF(cholmod_copy_sparse);
    RREGDEF(cholmod_dense_to_sparse);
    RREGDEF(cholmod_factor_to_sparse);
    RREGDEF(cholmod_factorize);
    RREGDEF(cholmod_factorize_p);
    RREGDEF(cholmod_finish);
    RREGDEF(cholmod_free_dense);
    RREGDEF(cholmod_free_factor);
    RREGDEF(cholmod_free_sparse);
    RREGDEF(cholmod_free_triplet);
    RREGDEF(cholmod_nnz);
    RREGDEF(cholmod_scale);
    RREGDEF(cholmod_sdmult);
    RREGDEF(cholmod_solve);
    RREGDEF(cholmod_solve2);
    RREGDEF(cholmod_sort);
    RREGDEF(cholmod_sparse_to_dense);
    RREGDEF(cholmod_sparse_to_triplet);
    RREGDEF(cholmod_speye);
    RREGDEF(cholmod_spsolve);
    RREGDEF(cholmod_ssmult);
    RREGDEF(cholmod_start);
    RREGDEF(cholmod_submatrix);
    RREGDEF(cholmod_transpose);
    RREGDEF(cholmod_triplet_to_sparse);
    RREGDEF(cholmod_vertcat);
    RREGDEF(cholmod_updown);

    RREGDEF(numeric_as_chm_dense);

    R_cholmod_start(&c);
//    R_cholmod_start(&cl); << TODO; needs more work in ./chm_common.c etc

    Matrix_DimNamesSym = install("Dimnames");
    Matrix_DimSym      = install("Dim");
    Matrix_LSym        = install("L");
    Matrix_QSym        = install("Q");
    Matrix_RSym        = install("R");
    Matrix_TSym        = install("T");
    Matrix_USym        = install("U");
    Matrix_VSym        = install("V");
    Matrix_betaSym     = install("beta");
    Matrix_diagSym     = install("diag");
    Matrix_factorSym   = install("factors");
    Matrix_iSym        = install("i");
    Matrix_jSym        = install("j");
    Matrix_lengthSym   = install("length");
    Matrix_pSym        = install("p");
    Matrix_permSym     = install("perm");
    Matrix_qSym        = install("q");
    Matrix_sdSym       = install("sd");
    Matrix_uploSym     = install("uplo");
    Matrix_xSym        = install("x");
    
    Matrix_NS = R_FindNamespace(mkString("Matrix"));
    if(Matrix_NS == R_UnboundValue)
	error(_("missing 'Matrix' namespace: should never happen"));
#ifdef Matrix_Debug
    if(isEnvironment(Matrix_NS))
	Rprintf("Matrix_NS: %s\n",
		CHAR(asChar(eval(lang2(install("format"), Matrix_NS),
				 R_GlobalEnv))));
    else
#else
    if(!isEnvironment(Matrix_NS))
#endif
	error(_("Matrix namespace not determined correctly"));

    Matrix_zzero.r = 0.0; Matrix_zone.r = 1.0;
    Matrix_zzero.i = 0.0; Matrix_zone.i = 0.0;
}

void R_unload_Matrix(DllInfo *dll)
{
    cholmod_finish(&c);
}
