#include "cholmod-api.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

/* Csparse.c : */
SEXP CsparseMatrix_validate_maybe_sorting(SEXP);
SEXP dgCMatrix_lusol(SEXP, SEXP);
SEXP dgCMatrix_qrsol(SEXP, SEXP, SEXP);
SEXP dgCMatrix_cholsol(SEXP, SEXP);
SEXP dtCMatrix_diag(SEXP, SEXP);
SEXP Csparse_dmperm(SEXP, SEXP, SEXP);
SEXP Csparse_writeMM(SEXP, SEXP);

/* Summary.c : */
SEXP R_dense_sum(SEXP, SEXP);
SEXP R_sparse_sum(SEXP, SEXP);
SEXP R_dense_prod(SEXP, SEXP);
SEXP R_sparse_prod(SEXP, SEXP);

/* aggregate.c : */
SEXP R_sparse_aggregate(SEXP);

/* attrib.c : */
SEXP R_DimNames_is_symmetric(SEXP);
SEXP R_symDN(SEXP);
SEXP R_set_factor(SEXP, SEXP, SEXP, SEXP);

/* band.c : */
SEXP R_dense_band(SEXP, SEXP, SEXP);
SEXP R_sparse_band(SEXP, SEXP, SEXP);

/* bind.c : */
SEXP R_bind(SEXP);

/* coerce.c : */
SEXP R_vector_as_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_matrix_as_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_sparse_as_dense(SEXP, SEXP);
SEXP R_diagonal_as_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_index_as_dense(SEXP, SEXP);
SEXP R_Vector_as_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_matrix_as_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_dense_as_sparse(SEXP, SEXP);
SEXP R_diagonal_as_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_index_as_sparse(SEXP, SEXP, SEXP);
SEXP R_dense_as_kind(SEXP, SEXP);
SEXP R_sparse_as_kind(SEXP, SEXP);
SEXP R_diagonal_as_kind(SEXP, SEXP);
SEXP R_index_as_kind(SEXP, SEXP);
SEXP R_dense_as_general(SEXP);
SEXP R_sparse_as_general(SEXP);
SEXP R_dense_as_unpacked(SEXP);
SEXP R_dense_as_packed(SEXP, SEXP, SEXP, SEXP);
SEXP R_sparse_as_Csparse(SEXP);
SEXP R_sparse_as_Rsparse(SEXP);
SEXP R_sparse_as_Tsparse(SEXP);
SEXP R_vector_as_Vector(SEXP, SEXP);
SEXP R_sparse_as_Vector(SEXP);
SEXP R_diagonal_as_Vector(SEXP);
SEXP R_index_as_Vector(SEXP);
SEXP R_Matrix_as_vector(SEXP);
SEXP R_Matrix_as_matrix(SEXP);
SEXP R_Matrix_as_unpacked(SEXP);
SEXP R_Matrix_as_packed(SEXP);
SEXP R_Matrix_as_Csparse(SEXP);
SEXP R_Matrix_as_Rsparse(SEXP);
SEXP R_Matrix_as_Tsparse(SEXP);
SEXP R_Matrix_as_Vector(SEXP);
SEXP R_Matrix_as_kind(SEXP, SEXP, SEXP);
SEXP R_Matrix_as_general(SEXP, SEXP);

/* colSums.c : */
SEXP R_dense_marginsum(SEXP, SEXP, SEXP, SEXP);
SEXP R_sparse_marginsum(SEXP, SEXP, SEXP, SEXP, SEXP);

/* determinant.c : */
SEXP denseLU_determinant(SEXP, SEXP);
SEXP denseBunchKaufman_determinant(SEXP, SEXP);
SEXP denseCholesky_determinant(SEXP, SEXP);
SEXP sparseQR_determinant(SEXP, SEXP);
SEXP sparseLU_determinant(SEXP, SEXP);
SEXP sparseCholesky_determinant(SEXP, SEXP, SEXP);

/* diag.c : */
SEXP R_dense_diag_get(SEXP, SEXP);
SEXP R_sparse_diag_get(SEXP, SEXP);
SEXP R_dense_diag_set(SEXP, SEXP);
SEXP R_sparse_diag_set(SEXP, SEXP);
SEXP R_sparse_diag_U2N(SEXP);
SEXP R_sparse_diag_N2U(SEXP);
SEXP denseCholesky_diag_get(SEXP, SEXP);
SEXP sparseCholesky_diag_get(SEXP, SEXP);

/* dropzero.c : */
SEXP R_sparse_dropzero(SEXP, SEXP);

/* expand.c : */
SEXP denseBunchKaufman_expand(SEXP);

/* expm.c : */
SEXP geMatrix_expm(SEXP);

/* factor.c : */
SEXP geMatrix_scf(SEXP, SEXP, SEXP);
SEXP syMatrix_scf(SEXP, SEXP, SEXP);
SEXP spMatrix_scf(SEXP, SEXP, SEXP);
SEXP geMatrix_trf(SEXP, SEXP);
SEXP syMatrix_trf(SEXP, SEXP);
SEXP spMatrix_trf(SEXP, SEXP);
SEXP poMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP ppMatrix_trf(SEXP, SEXP);
SEXP gCMatrix_orf(SEXP, SEXP, SEXP);
SEXP gCMatrix_trf(SEXP, SEXP, SEXP, SEXP);
SEXP pCMatrix_trf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* forceCanonical.c : */
SEXP R_dense_force_canonical(SEXP, SEXP);
SEXP R_sparse_force_canonical(SEXP, SEXP);

/* forceSymmetric.c : */
SEXP R_dense_force_symmetric(SEXP, SEXP, SEXP);
SEXP R_sparse_force_symmetric(SEXP, SEXP, SEXP);

/* isCanonical.c : */
SEXP R_dense_is_canonical(SEXP);
SEXP R_sparse_is_canonical(SEXP);

/* isDiagonal.c : */
SEXP R_dense_is_diagonal(SEXP);
SEXP R_sparse_is_diagonal(SEXP);

/* isSymmetric.c : */
SEXP R_dense_is_symmetric(SEXP, SEXP, SEXP, SEXP);
SEXP R_sparse_is_symmetric(SEXP, SEXP, SEXP, SEXP);

/* isTriangular.c : */
SEXP R_dense_is_triangular(SEXP, SEXP);
SEXP R_sparse_is_triangular(SEXP, SEXP);

/* matmult.c : */
SEXP R_dense_matmult(SEXP, SEXP, SEXP, SEXP);
SEXP R_sparse_matmult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_diagonal_matmult(SEXP, SEXP, SEXP, SEXP, SEXP);

/* norm.c : */
SEXP geMatrix_norm(SEXP, SEXP);
SEXP syMatrix_norm(SEXP, SEXP);
SEXP spMatrix_norm(SEXP, SEXP);
SEXP trMatrix_norm(SEXP, SEXP);
SEXP tpMatrix_norm(SEXP, SEXP);

/* objects.c : */
SEXP R_Matrix_class(SEXP, SEXP);
SEXP R_Matrix_kind(SEXP);
SEXP R_Matrix_shape(SEXP);
SEXP R_Matrix_repr(SEXP);

/* perm.c : */
SEXP R_isPerm(SEXP, SEXP);
SEXP R_signPerm(SEXP, SEXP);
SEXP R_invertPerm(SEXP, SEXP, SEXP);
SEXP R_asPerm(SEXP, SEXP, SEXP, SEXP);

/* rcond.c : */
SEXP geMatrix_rcond(SEXP, SEXP, SEXP);
SEXP syMatrix_rcond(SEXP, SEXP, SEXP);
SEXP spMatrix_rcond(SEXP, SEXP, SEXP);
SEXP poMatrix_rcond(SEXP, SEXP, SEXP);
SEXP ppMatrix_rcond(SEXP, SEXP, SEXP);
SEXP trMatrix_rcond(SEXP, SEXP);
SEXP tpMatrix_rcond(SEXP, SEXP);

/* skewpart.c : */
SEXP R_dense_skewpart(SEXP, SEXP);
SEXP R_sparse_skewpart(SEXP, SEXP);

/* solve.c : */
SEXP denseLU_solve(SEXP, SEXP);
SEXP denseBunchKaufman_solve(SEXP, SEXP);
SEXP denseCholesky_solve(SEXP, SEXP);
SEXP trMatrix_solve(SEXP, SEXP);
SEXP sparseLU_solve(SEXP, SEXP, SEXP);
SEXP sparseCholesky_solve(SEXP, SEXP, SEXP, SEXP);
SEXP tCMatrix_solve(SEXP, SEXP, SEXP);
SEXP sparseQR_matmult(SEXP, SEXP, SEXP, SEXP, SEXP);

/* subassign.c : */
SEXP nCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP lCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP iCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP dCsparse_subassign(SEXP, SEXP, SEXP, SEXP);
SEXP zCsparse_subassign(SEXP, SEXP, SEXP, SEXP);

/* subscript.c : */
SEXP R_subscript_1ary(SEXP, SEXP, SEXP);
SEXP R_subscript_1ary_2col(SEXP, SEXP, SEXP);
SEXP R_subscript_2ary(SEXP, SEXP, SEXP);

/* symmpart.c : */
SEXP R_dense_symmpart(SEXP, SEXP, SEXP);
SEXP R_sparse_symmpart(SEXP, SEXP, SEXP);

/* t.c : */
SEXP R_dense_transpose(SEXP, SEXP);
SEXP R_sparse_transpose(SEXP, SEXP, SEXP);

/* updown.c : */
SEXP sparseCholesky_updown(SEXP, SEXP, SEXP);
SEXP sparseCholesky_update(SEXP, SEXP, SEXP);

/* utils-R.c : */
SEXP R_index_triangle(SEXP, SEXP, SEXP, SEXP);
SEXP R_index_diagonal(SEXP, SEXP, SEXP);
SEXP R_nnz(SEXP, SEXP, SEXP);
SEXP R_all0(SEXP);
SEXP R_any0(SEXP);
SEXP Mmatrix(SEXP);
SEXP compressed_non_0_ij(SEXP, SEXP);
SEXP Matrix_expand_pointers(SEXP);
SEXP m_encodeInd(SEXP, SEXP, SEXP, SEXP);
SEXP m_encodeInd2(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Matrix_rle_d(SEXP, SEXP);
SEXP Matrix_rle_i(SEXP, SEXP);

/* validity.c : */
SEXP R_Dim_validate(SEXP);
SEXP R_DimNames_validate(SEXP, SEXP);
SEXP R_DimNames_fixup(SEXP);
SEXP Matrix_validate(SEXP);
SEXP nMatrix_validate(SEXP);
SEXP lMatrix_validate(SEXP);
SEXP iMatrix_validate(SEXP);
SEXP dMatrix_validate(SEXP);
SEXP zMatrix_validate(SEXP);
SEXP generalMatrix_validate(SEXP);
SEXP symmetricMatrix_validate(SEXP);
SEXP triangularMatrix_validate(SEXP);
SEXP unpackedMatrix_validate(SEXP);
SEXP packedMatrix_validate(SEXP);
SEXP CsparseMatrix_validate(SEXP);
SEXP RsparseMatrix_validate(SEXP);
SEXP TsparseMatrix_validate(SEXP);
SEXP diagonalMatrix_validate(SEXP);
SEXP indMatrix_validate(SEXP);
SEXP pMatrix_validate(SEXP);
SEXP sCMatrix_validate(SEXP);
SEXP tCMatrix_validate(SEXP);
SEXP sRMatrix_validate(SEXP);
SEXP tRMatrix_validate(SEXP);
SEXP sTMatrix_validate(SEXP);
SEXP tTMatrix_validate(SEXP);
SEXP xgCMatrix_validate(SEXP);
SEXP xsCMatrix_validate(SEXP);
SEXP xtCMatrix_validate(SEXP);
SEXP xgRMatrix_validate(SEXP);
SEXP xsRMatrix_validate(SEXP);
SEXP xtRMatrix_validate(SEXP);
SEXP xgTMatrix_validate(SEXP);
SEXP xsTMatrix_validate(SEXP);
SEXP xtTMatrix_validate(SEXP);
SEXP xpoMatrix_validate(SEXP);
SEXP xppMatrix_validate(SEXP);
SEXP xpCMatrix_validate(SEXP);
SEXP xpRMatrix_validate(SEXP);
SEXP xpTMatrix_validate(SEXP);
SEXP corMatrix_validate(SEXP);
SEXP copMatrix_validate(SEXP);
SEXP sparseVector_validate(SEXP);
SEXP lsparseVector_validate(SEXP);
SEXP isparseVector_validate(SEXP);
SEXP dsparseVector_validate(SEXP);
SEXP zsparseVector_validate(SEXP);
SEXP MatrixFactorization_validate(SEXP);
SEXP denseSchur_validate(SEXP);
SEXP denseQR_validate(SEXP);
SEXP denseLU_validate(SEXP);
SEXP denseBunchKaufman_validate(SEXP);
SEXP denseCholesky_validate(SEXP);
SEXP sparseQR_validate(SEXP);
SEXP sparseLU_validate(SEXP);
SEXP sparseCholesky_validate(SEXP);
SEXP simplicialCholesky_validate(SEXP);
SEXP supernodalCholesky_validate(SEXP);

/* version.c : */
SEXP R_Matrix_version(void);

#define     CALL_METHOD(name, n) {#name, (DL_FUNC) &name, n}
#define EXTERNAL_METHOD(name, n) {#name, (DL_FUNC) &name, n}
#define        REGISTER(name   ) R_RegisterCCallable("Matrix", #name, (DL_FUNC) name)

static R_CallMethodDef CallMethodTable[] = {
	/* Csparse.c : */
	CALL_METHOD(CsparseMatrix_validate_maybe_sorting, 1),
	CALL_METHOD(dgCMatrix_lusol, 2),
	CALL_METHOD(dgCMatrix_qrsol, 3),
	CALL_METHOD(dgCMatrix_cholsol, 2),
	CALL_METHOD(dtCMatrix_diag, 2),
	CALL_METHOD(Csparse_dmperm, 3),
	CALL_METHOD(Csparse_writeMM, 2),

	/* Summary.c : */
	CALL_METHOD(R_dense_sum, 2),
	CALL_METHOD(R_sparse_sum, 2),
	CALL_METHOD(R_dense_prod, 2),
	CALL_METHOD(R_sparse_prod, 2),

	/* aggregate.c : */
	CALL_METHOD(R_sparse_aggregate, 1),

	/* attrib.c : */
	CALL_METHOD(R_DimNames_is_symmetric, 1),
	CALL_METHOD(R_symDN, 1),
	CALL_METHOD(R_set_factor, 4),

	/* band.c : */
	CALL_METHOD(R_dense_band, 3),
	CALL_METHOD(R_sparse_band, 3),

	/* coerce.c : */
	CALL_METHOD(R_vector_as_dense, 9),
	CALL_METHOD(R_matrix_as_dense, 6),
	CALL_METHOD(R_sparse_as_dense, 2),
	CALL_METHOD(R_diagonal_as_dense, 6),
	CALL_METHOD(R_index_as_dense, 2),
	CALL_METHOD(R_Vector_as_sparse, 9),
	CALL_METHOD(R_matrix_as_sparse, 6),
	CALL_METHOD(R_dense_as_sparse, 2),
	CALL_METHOD(R_diagonal_as_sparse, 6),
	CALL_METHOD(R_index_as_sparse, 3),
	CALL_METHOD(R_dense_as_kind, 2),
	CALL_METHOD(R_sparse_as_kind, 2),
	CALL_METHOD(R_diagonal_as_kind, 2),
	CALL_METHOD(R_index_as_kind, 2),
	CALL_METHOD(R_dense_as_general, 1),
	CALL_METHOD(R_sparse_as_general, 1),
	CALL_METHOD(R_dense_as_unpacked, 1),
	CALL_METHOD(R_dense_as_packed, 4),
	CALL_METHOD(R_sparse_as_Csparse, 1),
	CALL_METHOD(R_sparse_as_Rsparse, 1),
	CALL_METHOD(R_sparse_as_Tsparse, 1),
	CALL_METHOD(R_vector_as_Vector, 2),
	CALL_METHOD(R_sparse_as_Vector, 1),
	CALL_METHOD(R_diagonal_as_Vector, 1),
	CALL_METHOD(R_index_as_Vector, 1),
	CALL_METHOD(R_Matrix_as_vector, 1),
	CALL_METHOD(R_Matrix_as_matrix, 1),
	CALL_METHOD(R_Matrix_as_unpacked, 1),
	CALL_METHOD(R_Matrix_as_packed, 1),
	CALL_METHOD(R_Matrix_as_Csparse, 1),
	CALL_METHOD(R_Matrix_as_Rsparse, 1),
	CALL_METHOD(R_Matrix_as_Tsparse, 1),
	CALL_METHOD(R_Matrix_as_Vector, 1),
	CALL_METHOD(R_Matrix_as_kind, 3),
	CALL_METHOD(R_Matrix_as_general, 2),

	/* colSums.c : */
	CALL_METHOD(R_dense_marginsum, 4),
	CALL_METHOD(R_sparse_marginsum, 5),

	/* determinant.c : */
	CALL_METHOD(denseLU_determinant, 2),
	CALL_METHOD(denseBunchKaufman_determinant, 2),
	CALL_METHOD(denseCholesky_determinant, 2),
	CALL_METHOD(sparseQR_determinant, 2),
	CALL_METHOD(sparseLU_determinant, 2),
	CALL_METHOD(sparseCholesky_determinant, 3),

	/* diag.c : */
	CALL_METHOD(R_dense_diag_get, 2),
	CALL_METHOD(R_sparse_diag_get, 2),
	CALL_METHOD(R_dense_diag_set, 2),
	CALL_METHOD(R_sparse_diag_set, 2),
	CALL_METHOD(R_sparse_diag_U2N, 1),
	CALL_METHOD(R_sparse_diag_N2U, 1),
	CALL_METHOD(denseCholesky_diag_get, 2),
	CALL_METHOD(sparseCholesky_diag_get, 2),

	/* dropzero.c : */
	CALL_METHOD(R_sparse_dropzero, 2),

	/* expand.c : */
	CALL_METHOD(denseBunchKaufman_expand, 1),

	/* expm.c : */
	CALL_METHOD(geMatrix_expm, 1),

	/* factor.c : */
	CALL_METHOD(geMatrix_scf, 3),
	CALL_METHOD(syMatrix_scf, 3),
	CALL_METHOD(spMatrix_scf, 3),
	CALL_METHOD(geMatrix_trf, 2),
	CALL_METHOD(syMatrix_trf, 2),
	CALL_METHOD(spMatrix_trf, 2),
	CALL_METHOD(poMatrix_trf, 4),
	CALL_METHOD(ppMatrix_trf, 2),
	CALL_METHOD(gCMatrix_orf, 3),
	CALL_METHOD(gCMatrix_trf, 4),
	CALL_METHOD(pCMatrix_trf, 6),

	/* forceCanonical.c : */
	CALL_METHOD(R_dense_force_canonical, 2),
	CALL_METHOD(R_sparse_force_canonical, 2),

	/* forceSymmetric.c : */
	CALL_METHOD(R_dense_force_symmetric, 3),
	CALL_METHOD(R_sparse_force_symmetric, 3),

	/* isCanonical.c : */
	CALL_METHOD(R_dense_is_canonical, 1),
	CALL_METHOD(R_sparse_is_canonical, 1),

	/* isDiagonal.c : */
	CALL_METHOD(R_dense_is_diagonal, 1),
	CALL_METHOD(R_sparse_is_diagonal, 1),

	/* isSymmetric.c : */
	CALL_METHOD(R_dense_is_symmetric, 4),
	CALL_METHOD(R_sparse_is_symmetric, 4),

	/* isTriangular.c : */
	CALL_METHOD(R_dense_is_triangular, 2),
	CALL_METHOD(R_sparse_is_triangular, 2),

	/* matmult.c : */
	CALL_METHOD(R_dense_matmult, 4),
	CALL_METHOD(R_sparse_matmult, 6),
	CALL_METHOD(R_diagonal_matmult, 5),

	/* norm.c : */
	CALL_METHOD(geMatrix_norm, 2),
	CALL_METHOD(syMatrix_norm, 2),
	CALL_METHOD(spMatrix_norm, 2),
	CALL_METHOD(trMatrix_norm, 2),
	CALL_METHOD(tpMatrix_norm, 2),

	/* objects.c : */
	CALL_METHOD(R_Matrix_class, 2),
	CALL_METHOD(R_Matrix_kind, 1),
	CALL_METHOD(R_Matrix_shape, 1),
	CALL_METHOD(R_Matrix_repr, 1),

	/* perm.c : */
	CALL_METHOD(R_isPerm, 2),
	CALL_METHOD(R_signPerm, 2),
	CALL_METHOD(R_invertPerm, 3),
	CALL_METHOD(R_asPerm, 4),

	/* rcond.c : */
	CALL_METHOD(geMatrix_rcond, 3),
	CALL_METHOD(syMatrix_rcond, 3),
	CALL_METHOD(spMatrix_rcond, 3),
	CALL_METHOD(poMatrix_rcond, 3),
	CALL_METHOD(ppMatrix_rcond, 3),
	CALL_METHOD(trMatrix_rcond, 2),
	CALL_METHOD(tpMatrix_rcond, 2),

	/* skewpart.c : */
	CALL_METHOD(R_dense_skewpart, 2),
	CALL_METHOD(R_sparse_skewpart, 2),

	/* solve.c : */
	CALL_METHOD(denseLU_solve, 2),
	CALL_METHOD(denseBunchKaufman_solve, 2),
	CALL_METHOD(denseCholesky_solve, 2),
	CALL_METHOD(trMatrix_solve, 2),
	CALL_METHOD(sparseLU_solve, 3),
	CALL_METHOD(sparseCholesky_solve, 4),
	CALL_METHOD(tCMatrix_solve, 3),
	CALL_METHOD(sparseQR_matmult, 5),

	/* subassign.c : */
	CALL_METHOD(nCsparse_subassign, 4),
	CALL_METHOD(lCsparse_subassign, 4),
	CALL_METHOD(iCsparse_subassign, 4),
	CALL_METHOD(dCsparse_subassign, 4),
	CALL_METHOD(zCsparse_subassign, 4),

	/* subscript.c : */
	CALL_METHOD(R_subscript_1ary, 3),
	CALL_METHOD(R_subscript_1ary_2col, 3),
	CALL_METHOD(R_subscript_2ary, 3),

	/* symmpart.c : */
	CALL_METHOD(R_dense_symmpart, 3),
	CALL_METHOD(R_sparse_symmpart, 3),

	/* t.c : */
	CALL_METHOD(R_dense_transpose, 2),
	CALL_METHOD(R_sparse_transpose, 3),

	/* updown.c : */
	CALL_METHOD(sparseCholesky_updown, 3),
	CALL_METHOD(sparseCholesky_update, 3),

	/* utils-R.c : */
	CALL_METHOD(R_index_triangle, 4),
	CALL_METHOD(R_index_diagonal, 3),
	CALL_METHOD(R_nnz, 3),
	CALL_METHOD(R_all0, 1),
	CALL_METHOD(R_any0, 1),
	CALL_METHOD(compressed_non_0_ij, 2),
	CALL_METHOD(Matrix_expand_pointers, 1),
	CALL_METHOD(m_encodeInd,  4),
	CALL_METHOD(m_encodeInd2, 5),
	CALL_METHOD(Matrix_rle_i, 2),
	CALL_METHOD(Matrix_rle_d, 2),

	/* validity.c : */
	CALL_METHOD(R_Dim_validate, 1),
	CALL_METHOD(R_DimNames_validate, 2),
	CALL_METHOD(R_DimNames_fixup, 1),
	CALL_METHOD(Matrix_validate, 1),
	CALL_METHOD(nMatrix_validate, 1),
	CALL_METHOD(lMatrix_validate, 1),
	CALL_METHOD(iMatrix_validate, 1),
	CALL_METHOD(dMatrix_validate, 1),
	CALL_METHOD(zMatrix_validate, 1),
	CALL_METHOD(generalMatrix_validate, 1),
	CALL_METHOD(symmetricMatrix_validate, 1),
	CALL_METHOD(triangularMatrix_validate, 1),
	CALL_METHOD(unpackedMatrix_validate, 1),
	CALL_METHOD(packedMatrix_validate, 1),
	CALL_METHOD(CsparseMatrix_validate, 1),
	CALL_METHOD(RsparseMatrix_validate, 1),
	CALL_METHOD(TsparseMatrix_validate, 1),
	CALL_METHOD(diagonalMatrix_validate, 1),
	CALL_METHOD(indMatrix_validate, 1),
	CALL_METHOD(pMatrix_validate, 1),
	CALL_METHOD(sCMatrix_validate, 1),
	CALL_METHOD(tCMatrix_validate, 1),
	CALL_METHOD(sRMatrix_validate, 1),
	CALL_METHOD(tRMatrix_validate, 1),
	CALL_METHOD(sTMatrix_validate, 1),
	CALL_METHOD(tTMatrix_validate, 1),
	CALL_METHOD(xgCMatrix_validate, 1),
	CALL_METHOD(xsCMatrix_validate, 1),
	CALL_METHOD(xtCMatrix_validate, 1),
	CALL_METHOD(xgRMatrix_validate, 1),
	CALL_METHOD(xsRMatrix_validate, 1),
	CALL_METHOD(xtRMatrix_validate, 1),
	CALL_METHOD(xgTMatrix_validate, 1),
	CALL_METHOD(xsTMatrix_validate, 1),
	CALL_METHOD(xtTMatrix_validate, 1),
	CALL_METHOD(xpoMatrix_validate, 1),
	CALL_METHOD(xppMatrix_validate, 1),
	CALL_METHOD(xpCMatrix_validate, 1),
	CALL_METHOD(xpRMatrix_validate, 1),
	CALL_METHOD(xpTMatrix_validate, 1),
	CALL_METHOD(corMatrix_validate, 1),
	CALL_METHOD(copMatrix_validate, 1),
	CALL_METHOD(sparseVector_validate, 1),
	CALL_METHOD(lsparseVector_validate, 1),
	CALL_METHOD(isparseVector_validate, 1),
	CALL_METHOD(dsparseVector_validate, 1),
	CALL_METHOD(zsparseVector_validate, 1),
	CALL_METHOD(MatrixFactorization_validate, 1),
	CALL_METHOD(denseSchur_validate, 1),
	CALL_METHOD(denseQR_validate, 1),
	CALL_METHOD(denseLU_validate, 1),
	CALL_METHOD(denseBunchKaufman_validate, 1),
	CALL_METHOD(denseCholesky_validate, 1),
	CALL_METHOD(sparseQR_validate, 1),
	CALL_METHOD(sparseLU_validate, 1),
	CALL_METHOD(sparseCholesky_validate, 1),
	CALL_METHOD(simplicialCholesky_validate, 1),
	CALL_METHOD(supernodalCholesky_validate, 1),

	/* version.c */
	CALL_METHOD(R_Matrix_version, 0),

	{NULL, NULL, 0}
};

static const R_ExternalMethodDef ExternalMethodTable[] = {
	/* bind.c */
	EXTERNAL_METHOD(R_bind, -1),

	/* utils-R.c */
	EXTERNAL_METHOD(Mmatrix, 7),

	{NULL, NULL, 0}
};

#include "Msymbols.h"
Rcomplex Matrix_zzero, Matrix_zunit, Matrix_zna;

void attribute_visible R_init_Matrix(DllInfo *info)
{
	R_registerRoutines(info, NULL, CallMethodTable, NULL, ExternalMethodTable);
	R_useDynamicSymbols(info, FALSE);

	/* CHOLMOD : */
	REGISTER(cholmod_aat);
	REGISTER(cholmod_add);
	REGISTER(cholmod_allocate_dense);
	REGISTER(cholmod_allocate_factor);
	REGISTER(cholmod_allocate_sparse);
	REGISTER(cholmod_allocate_triplet);
	REGISTER(cholmod_analyze);
	REGISTER(cholmod_analyze_p);
	REGISTER(cholmod_band_inplace);
	REGISTER(cholmod_change_factor);
	REGISTER(cholmod_check_common);
	REGISTER(cholmod_check_dense);
	REGISTER(cholmod_check_factor);
	REGISTER(cholmod_check_sparse);
	REGISTER(cholmod_check_triplet);
	REGISTER(cholmod_copy);
	REGISTER(cholmod_copy_dense);
	REGISTER(cholmod_copy_factor);
	REGISTER(cholmod_copy_sparse);
	REGISTER(cholmod_copy_triplet);
	REGISTER(cholmod_defaults);
	REGISTER(cholmod_dense_to_sparse);
	REGISTER(cholmod_factor_to_sparse);
	REGISTER(cholmod_factorize);
	REGISTER(cholmod_factorize_p);
	REGISTER(cholmod_finish);
	REGISTER(cholmod_free_dense);
	REGISTER(cholmod_free_factor);
	REGISTER(cholmod_free_sparse);
	REGISTER(cholmod_free_triplet);
	REGISTER(cholmod_horzcat);
	REGISTER(cholmod_nnz);
	REGISTER(cholmod_scale);
	REGISTER(cholmod_sdmult);
	REGISTER(cholmod_solve);
	REGISTER(cholmod_solve2);
	REGISTER(cholmod_sort);
	REGISTER(cholmod_sparse_to_dense);
	REGISTER(cholmod_sparse_to_triplet);
	REGISTER(cholmod_speye);
	REGISTER(cholmod_spsolve);
	REGISTER(cholmod_ssmult);
	REGISTER(cholmod_start);
	REGISTER(cholmod_submatrix);
	REGISTER(cholmod_transpose);
	REGISTER(cholmod_triplet_to_sparse);
	REGISTER(cholmod_updown);
	REGISTER(cholmod_vertcat);

	/* Matrix : */
	REGISTER(sexp_as_cholmod_factor);
	REGISTER(sexp_as_cholmod_sparse);
	REGISTER(sexp_as_cholmod_triplet);
	REGISTER(sexp_as_cholmod_dense);
	REGISTER(numeric_as_cholmod_dense);
	REGISTER(cholmod_factor_as_sexp);
	REGISTER(cholmod_sparse_as_sexp);
	REGISTER(cholmod_triplet_as_sexp);
	REGISTER(cholmod_dense_as_sexp);
	REGISTER(cholmod_factor_ldetA);
	REGISTER(cholmod_factor_update);

	Matrix_DimNamesSym = Rf_install("Dimnames");
	Matrix_DimSym      = Rf_install("Dim");
	Matrix_LSym        = Rf_install("L");
	Matrix_RSym        = Rf_install("R");
	Matrix_USym        = Rf_install("U");
	Matrix_VSym        = Rf_install("V");
	Matrix_betaSym     = Rf_install("beta");
	Matrix_colcountSym = Rf_install("colcount");
	Matrix_diagSym     = Rf_install("diag");
	Matrix_factorsSym  = Rf_install("factors");
	Matrix_iSym        = Rf_install("i");
	Matrix_isllSym     = Rf_install("is_ll");
	Matrix_ismtSym     = Rf_install("is_monotonic");
	Matrix_jSym        = Rf_install("j");
	Matrix_lengthSym   = Rf_install("length");
	Matrix_marginSym   = Rf_install("margin");
	Matrix_maxcsizeSym = Rf_install("maxcsize");
	Matrix_maxesizeSym = Rf_install("maxesize");
	Matrix_minorSym    = Rf_install("minor");
	Matrix_nextSym     = Rf_install("next");
	Matrix_nzSym       = Rf_install("nz");
	Matrix_orderingSym = Rf_install("ordering");
	Matrix_pSym        = Rf_install("p");
	Matrix_permSym     = Rf_install("perm");
	Matrix_piSym       = Rf_install("pi");
	Matrix_prevSym     = Rf_install("prev");
	Matrix_pxSym       = Rf_install("px");
	Matrix_qSym        = Rf_install("q");
	Matrix_sSym        = Rf_install("s");
	Matrix_sdSym       = Rf_install("sd");
	Matrix_superSym    = Rf_install("super");
	Matrix_transSym    = Rf_install("trans");
	Matrix_uploSym     = Rf_install("uplo");
	Matrix_valuesSym   = Rf_install("values");
	Matrix_vectorsSym  = Rf_install("vectors");
	Matrix_xSym        = Rf_install("x");

	Matrix_LChar = Rf_mkChar("L");
	Matrix_TChar = Rf_mkChar("T");
	Matrix_UChar = Rf_mkChar("U");

	Matrix_zzero.r = 0.0; Matrix_zunit.r = 1.0; Matrix_zna.r = NA_REAL;
	Matrix_zzero.i = 0.0; Matrix_zunit.i = 0.0; Matrix_zna.i = NA_REAL;

	Matrix_cholmod_start(&c);
	return;
}

void R_unload_Matrix(DllInfo *info)
{
	Matrix_cholmod_finish(&c);
	return;
}
