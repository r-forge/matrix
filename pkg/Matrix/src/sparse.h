#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Mutils.h"

SEXP sparse_drop0(SEXP from, const char *class, double tol);
SEXP R_sparse_drop0(SEXP from, SEXP tol);

SEXP sparse_band(SEXP from, const char *class, int a, int b);
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);

SEXP R_sparse_diag_get(SEXP obj, SEXP nms);
SEXP R_sparse_diag_set(SEXP obj, SEXP val);

SEXP sparse_diag_U2N(SEXP from, const char *class);
SEXP R_sparse_diag_U2N(SEXP from);

SEXP sparse_diag_N2U(SEXP from, const char *class);
SEXP R_sparse_diag_N2U(SEXP from);

SEXP sparse_transpose(SEXP from, const char *class, int lazy);
SEXP R_sparse_transpose(SEXP from, SEXP lazy);

SEXP sparse_force_symmetric(SEXP from, const char *class, char ul);
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo);

SEXP R_sparse_symmpart(SEXP from);
SEXP R_sparse_skewpart(SEXP from);

SEXP Tsparse_aggregate(SEXP from);

SEXP Csparse_is_diagonal(SEXP obj);
SEXP Rsparse_is_diagonal(SEXP obj);
SEXP Tsparse_is_diagonal(SEXP obj);
SEXP Csparse_is_triangular(SEXP obj, SEXP upper);
SEXP Rsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Csparse_is_symmetric(SEXP obj, SEXP checkDN);
SEXP Rsparse_is_symmetric(SEXP obj, SEXP checkDN);
#if 0 /* unimplemented ... currently going via CsparseMatrix */
SEXP Tsparse_is_symmetric(SEXP obj, SEXP checkDN);
#endif

SEXP CRsparse_colSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse);
SEXP CRsparse_rowSums(SEXP obj, SEXP narm, SEXP mean, SEXP sparse);

#define SPARSE_CASES(_SEXPTYPE_, _DO_) \
do { \
	switch (_SEXPTYPE_) { \
	case LGLSXP: \
		_DO_(int, LOGICAL, 0, 1, ISNZ_LOGICAL); \
		break; \
	case INTSXP: \
		_DO_(int, INTEGER, 0, 1, ISNZ_INTEGER); \
		break; \
	case REALSXP: \
		_DO_(double, REAL, 0.0, 1.0, ISNZ_REAL); \
		break; \
	case CPLXSXP: \
		_DO_(Rcomplex, COMPLEX, Matrix_zzero, Matrix_zone, ISNZ_COMPLEX); \
		break; \
	default: \
		break; \
	} \
} while (0)

#endif
