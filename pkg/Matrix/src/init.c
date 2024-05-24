#include "Mdefines.h"
#include "Csparse.h"
#include "attrib.h"
#include "bind.h"
#include "cholmod-api.h"
#include "coerce.h"
#include "dense.h"
#include "determinant.h"
#include "expm.h"
#include "factor.h"
#include "kappa.h"
#include "matmult.h"
#include "objects.h"
#include "perm.h"
#include "solve.h"
#include "sparse.h"
#include "subassign.h"
#include "subscript.h"
#include "utils-R.h"
#include "validity.h"
#include "vector.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "Msymbols.h"
Rcomplex Matrix_zzero, Matrix_zone, Matrix_zna;

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define  EXTDEF(name, n)  {#name, (DL_FUNC) &name, n}
#define RREGDEF(name   )  R_RegisterCCallable("Matrix", #name, (DL_FUNC) name)

static R_CallMethodDef CallEntries[] = {
	CALLDEF(R_all0, 1),
	CALLDEF(R_any0, 1),

	CALLDEF(compressed_non_0_ij, 2),
	CALLDEF(Matrix_expand_pointers, 1),
	CALLDEF(m_encodeInd,  4),
	CALLDEF(m_encodeInd2, 5),

	CALLDEF(Matrix_rle_i, 2),
	CALLDEF(Matrix_rle_d, 2),

	CALLDEF(CsparseMatrix_validate_maybe_sorting, 1),

	CALLDEF(dgCMatrix_lusol, 2),
	CALLDEF(dgCMatrix_qrsol, 3),
	CALLDEF(dgCMatrix_cholsol, 2),
	CALLDEF(dtCMatrix_diag, 2),

	CALLDEF(Csparse_dmperm, 3),
	CALLDEF(Csparse_writeMM, 2),

	CALLDEF(nCsparse_subassign, 4),
	CALLDEF(lCsparse_subassign, 4),
	CALLDEF(iCsparse_subassign, 4),
	CALLDEF(dCsparse_subassign, 4),
	CALLDEF(zCsparse_subassign, 4),

	CALLDEF(Matrix_validate, 1),

	CALLDEF(nMatrix_validate, 1),
	CALLDEF(lMatrix_validate, 1),
	CALLDEF(iMatrix_validate, 1),
	CALLDEF(dMatrix_validate, 1),
	CALLDEF(zMatrix_validate, 1),

	CALLDEF(generalMatrix_validate, 1),
	CALLDEF(symmetricMatrix_validate, 1),
	CALLDEF(triangularMatrix_validate, 1),
	CALLDEF(diagonalMatrix_validate, 1),
	CALLDEF(indMatrix_validate, 1),
	CALLDEF(pMatrix_validate, 1),

	CALLDEF(unpackedMatrix_validate, 1),
	CALLDEF(packedMatrix_validate, 1),
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

	CALLDEF(xpoMatrix_validate, 1),
	CALLDEF(xppMatrix_validate, 1),
	CALLDEF(xpCMatrix_validate, 1),
	CALLDEF(xpRMatrix_validate, 1),
	CALLDEF(xpTMatrix_validate, 1),

	CALLDEF(corMatrix_validate, 1),
	CALLDEF(copMatrix_validate, 1),

	CALLDEF(sparseVector_validate, 1),
	CALLDEF(lsparseVector_validate, 1),
	CALLDEF(isparseVector_validate, 1),
	CALLDEF(dsparseVector_validate, 1),
	CALLDEF(zsparseVector_validate, 1),

	CALLDEF(MatrixFactorization_validate, 1),

	CALLDEF(denseSchur_validate, 1),
	CALLDEF(denseQR_validate, 1),
	CALLDEF(denseLU_validate, 1),
	CALLDEF(denseBunchKaufman_validate, 1),
	CALLDEF(denseCholesky_validate, 1),

	CALLDEF(sparseQR_validate, 1),
	CALLDEF(sparseLU_validate, 1),
	CALLDEF(sparseCholesky_validate, 1),
	CALLDEF(simplicialCholesky_validate, 1),
	CALLDEF(supernodalCholesky_validate, 1),

	CALLDEF(R_Dim_validate, 1),
	CALLDEF(R_DimNames_validate, 2),

	CALLDEF(R_DimNames_fixup, 1),
	CALLDEF(R_DimNames_is_symmetric, 1),
	CALLDEF(R_symDN, 1),
	CALLDEF(R_revDN, 1),
	CALLDEF(R_Matrix_nonvirtual, 2),
	CALLDEF(R_Matrix_kind, 1),
	CALLDEF(R_Matrix_shape, 1),
	CALLDEF(R_Matrix_repr, 1),
	CALLDEF(R_index_triangle, 4),
	CALLDEF(R_index_diagonal, 3),
	CALLDEF(R_nnz, 3),
	CALLDEF(R_isPerm, 2),
	CALLDEF(R_signPerm, 2),
	CALLDEF(R_invertPerm, 3),
	CALLDEF(R_asPerm, 4),
	CALLDEF(R_set_factor, 4),

	CALLDEF(R_subscript_1ary, 2),
	CALLDEF(R_subscript_1ary_mat, 2),
	CALLDEF(R_subscript_2ary, 3),

	CALLDEF(R_dense_band, 3),
	CALLDEF(R_dense_diag_get, 2),
	CALLDEF(R_dense_diag_set, 2),
	CALLDEF(R_dense_transpose, 2),
	CALLDEF(R_dense_force_symmetric, 3),
	CALLDEF(R_dense_symmpart, 2),
	CALLDEF(R_dense_skewpart, 2),
	CALLDEF(R_dense_is_symmetric, 3),
	CALLDEF(R_dense_is_triangular, 2),
	CALLDEF(R_dense_is_diagonal, 1),
	CALLDEF(R_dense_marginsum, 4),
	CALLDEF(R_dense_sum, 2),
	CALLDEF(R_dense_prod, 2),

	CALLDEF(R_sparse_drop0, 2),
	CALLDEF(R_sparse_diag_U2N, 1),
	CALLDEF(R_sparse_diag_N2U, 1),
	CALLDEF(R_sparse_band, 3),
	CALLDEF(R_sparse_diag_get, 2),
	CALLDEF(R_sparse_diag_set, 2),
	CALLDEF(R_sparse_transpose, 2),
	CALLDEF(R_sparse_force_symmetric, 3),
	CALLDEF(R_sparse_symmpart, 1),
	CALLDEF(R_sparse_skewpart, 1),
	CALLDEF(R_sparse_is_symmetric, 2),
	CALLDEF(R_sparse_is_triangular, 2),
	CALLDEF(R_sparse_is_diagonal, 1),
	CALLDEF(R_sparse_marginsum, 5),
	CALLDEF(R_sparse_sum, 2),
	CALLDEF(R_sparse_prod, 2),
	CALLDEF(Tsparse_aggregate, 1),

	CALLDEF(R_dense_matmult, 4),
	CALLDEF(R_sparse_matmult, 6),
	CALLDEF(R_diagonal_matmult, 5),

	CALLDEF(geMatrix_expm, 1),

	CALLDEF(geMatrix_scf, 3),
	CALLDEF(syMatrix_scf, 3),
	CALLDEF(spMatrix_scf, 3),
	CALLDEF(geMatrix_trf, 2),
	CALLDEF(syMatrix_trf, 2),
	CALLDEF(spMatrix_trf, 2),
	CALLDEF(poMatrix_trf, 4),
	CALLDEF(ppMatrix_trf, 2),
	CALLDEF(gCMatrix_orf, 3),
	CALLDEF(gCMatrix_trf, 4),
	CALLDEF(pCMatrix_trf, 6),

	CALLDEF(sparseCholesky_update, 3),
	CALLDEF(sparseCholesky_updown, 3),
	CALLDEF(sparseCholesky_diag_get, 2),
	CALLDEF(denseBunchKaufman_expand, 1),

	CALLDEF(denseLU_solve, 2),
	CALLDEF(denseBunchKaufman_solve, 2),
	CALLDEF(denseCholesky_solve, 2),
	CALLDEF(trMatrix_solve, 2),
	CALLDEF(sparseLU_solve, 3),
	CALLDEF(sparseCholesky_solve, 4),
	CALLDEF(tCMatrix_solve, 3),
	CALLDEF(sparseQR_matmult, 5),

	CALLDEF(denseLU_determinant, 2),
	CALLDEF(denseBunchKaufman_determinant, 2),
	CALLDEF(denseCholesky_determinant, 2),
	CALLDEF(sparseQR_determinant, 2),
	CALLDEF(sparseLU_determinant, 2),
	CALLDEF(sparseCholesky_determinant, 3),

	CALLDEF(geMatrix_norm, 2),
	CALLDEF(syMatrix_norm, 2),
	CALLDEF(spMatrix_norm, 2),
	CALLDEF(trMatrix_norm, 2),
	CALLDEF(tpMatrix_norm, 2),

	CALLDEF(geMatrix_rcond, 3),
	CALLDEF(syMatrix_rcond, 3),
	CALLDEF(spMatrix_rcond, 3),
	CALLDEF(poMatrix_rcond, 3),
	CALLDEF(ppMatrix_rcond, 3),
	CALLDEF(trMatrix_rcond, 2),
	CALLDEF(tpMatrix_rcond, 2),

	CALLDEF(v2spV, 1),
	CALLDEF(CR2spV, 1),

	CALLDEF(R_vector_as_dense, 9),
	CALLDEF(R_matrix_as_dense, 6),
	CALLDEF(R_sparse_as_dense, 2),
	CALLDEF(R_diagonal_as_dense, 6),
	CALLDEF(R_index_as_dense, 2),
	CALLDEF(R_vector_as_sparse, 9),
	CALLDEF(R_matrix_as_sparse, 6),
	CALLDEF(R_dense_as_sparse, 2),
	CALLDEF(R_diagonal_as_sparse, 6),
	CALLDEF(R_index_as_sparse, 3),
	CALLDEF(R_dense_as_kind, 2),
	CALLDEF(R_sparse_as_kind, 2),
	CALLDEF(R_diagonal_as_kind, 2),
	CALLDEF(R_index_as_kind, 2),
	CALLDEF(R_dense_as_general, 1),
	CALLDEF(R_sparse_as_general, 1),
	CALLDEF(R_dense_as_unpacked, 1),
	CALLDEF(R_dense_as_packed, 3),
	CALLDEF(R_sparse_as_Csparse, 1),
	CALLDEF(R_sparse_as_Rsparse, 1),
	CALLDEF(R_sparse_as_Tsparse, 1),
	CALLDEF(R_Matrix_as_vector, 1),
	CALLDEF(R_Matrix_as_matrix, 1),
	CALLDEF(R_Matrix_as_unpacked, 1),
	CALLDEF(R_Matrix_as_packed, 1),
	CALLDEF(R_Matrix_as_Csparse, 1),
	CALLDEF(R_Matrix_as_Rsparse, 1),
	CALLDEF(R_Matrix_as_Tsparse, 1),
	CALLDEF(R_Matrix_as_kind, 3),
	CALLDEF(R_Matrix_as_general, 2),

	CALLDEF(R_Matrix_version, 0),

	{NULL, NULL, 0}
};

static const R_ExternalMethodDef ExtEntries[] = {
	EXTDEF(Mmatrix, 7),
	EXTDEF(R_bind, -1),
	{NULL, NULL, 0}
};

void attribute_visible R_init_Matrix(DllInfo *info)
{
	R_registerRoutines(info, NULL, CallEntries, NULL, ExtEntries);
	R_useDynamicSymbols(info, FALSE);

/* These are callable from other packages' C code: */

	/* CHOLMOD: */
	RREGDEF(cholmod_aat);
	RREGDEF(cholmod_add);
	RREGDEF(cholmod_allocate_dense);
	RREGDEF(cholmod_allocate_factor);
	RREGDEF(cholmod_allocate_sparse);
	RREGDEF(cholmod_allocate_triplet);
	RREGDEF(cholmod_analyze);
	RREGDEF(cholmod_analyze_p);
	RREGDEF(cholmod_band_inplace);
	RREGDEF(cholmod_change_factor);
	RREGDEF(cholmod_check_common);
	RREGDEF(cholmod_check_dense);
	RREGDEF(cholmod_check_factor);
	RREGDEF(cholmod_check_sparse);
	RREGDEF(cholmod_check_triplet);
	RREGDEF(cholmod_copy);
	RREGDEF(cholmod_copy_dense);
	RREGDEF(cholmod_copy_factor);
	RREGDEF(cholmod_copy_sparse);
	RREGDEF(cholmod_copy_triplet);
	RREGDEF(cholmod_defaults);
	RREGDEF(cholmod_dense_to_sparse);
	RREGDEF(cholmod_factor_to_sparse);
	RREGDEF(cholmod_factorize);
	RREGDEF(cholmod_factorize_p);
	RREGDEF(cholmod_finish);
	RREGDEF(cholmod_free_dense);
	RREGDEF(cholmod_free_factor);
	RREGDEF(cholmod_free_sparse);
	RREGDEF(cholmod_free_triplet);
	RREGDEF(cholmod_horzcat);
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
	RREGDEF(cholmod_updown);
	RREGDEF(cholmod_vertcat);
	/* Matrix: */
	RREGDEF(sexp_as_cholmod_factor);
	RREGDEF(sexp_as_cholmod_sparse);
	RREGDEF(sexp_as_cholmod_triplet);
	RREGDEF(sexp_as_cholmod_dense);
	RREGDEF(numeric_as_cholmod_dense);
	RREGDEF(cholmod_factor_as_sexp);
	RREGDEF(cholmod_sparse_as_sexp);
	RREGDEF(cholmod_triplet_as_sexp);
	RREGDEF(cholmod_dense_as_sexp);
	RREGDEF(cholmod_factor_ldetA);
	RREGDEF(cholmod_factor_update);

	Matrix_DimNamesSym = install("Dimnames");
	Matrix_DimSym      = install("Dim");
	Matrix_LSym        = install("L");
	Matrix_RSym        = install("R");
	Matrix_USym        = install("U");
	Matrix_VSym        = install("V");
	Matrix_betaSym     = install("beta");
	Matrix_colcountSym = install("colcount");
	Matrix_diagSym     = install("diag");
	Matrix_factorsSym  = install("factors");
	Matrix_iSym        = install("i");
	Matrix_isllSym     = install("is_ll");
	Matrix_ismtSym     = install("is_monotonic");
	Matrix_jSym        = install("j");
	Matrix_lengthSym   = install("length");
	Matrix_marginSym   = install("margin");
	Matrix_maxcsizeSym = install("maxcsize");
	Matrix_maxesizeSym = install("maxesize");
	Matrix_minorSym    = install("minor");
	Matrix_nextSym     = install("next");
	Matrix_nzSym       = install("nz");
	Matrix_orderingSym = install("ordering");
	Matrix_pSym        = install("p");
	Matrix_permSym     = install("perm");
	Matrix_piSym       = install("pi");
	Matrix_prevSym     = install("prev");
	Matrix_pxSym       = install("px");
	Matrix_qSym        = install("q");
	Matrix_sSym        = install("s");
	Matrix_sdSym       = install("sd");
	Matrix_superSym    = install("super");
	Matrix_transSym    = install("trans");
	Matrix_uploSym     = install("uplo");
	Matrix_valuesSym   = install("values");
	Matrix_vectorsSym  = install("vectors");
	Matrix_xSym        = install("x");

	Matrix_zzero.r = 0.0; Matrix_zone.r = 1.0; Matrix_zna.r = NA_REAL;
	Matrix_zzero.i = 0.0; Matrix_zone.i = 0.0; Matrix_zna.i = NA_REAL;

	Matrix_cholmod_start(&c);
	return;
}

void R_unload_Matrix(DllInfo *info)
{
	Matrix_cholmod_finish(&c);
	return;
}
