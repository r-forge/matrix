#ifndef MATRIX_SPARSE_H
#define MATRIX_SPARSE_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP sparse_as_dense(SEXP from, int packed);
SEXP R_sparse_as_dense(SEXP from, SEXP packed);
SEXP R_sparse_as_matrix(SEXP from);
SEXP R_sparse_as_vector(SEXP from);
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0);
SEXP R_sparse_as_general(SEXP from);
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0);

SEXP R_sparse_drop0(SEXP from);
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2);
SEXP R_sparse_transpose(SEXP from);
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo_to);

SEXP CRsparse_as_Tsparse(SEXP from);
SEXP tCRsparse_as_RCsparse(SEXP from);

SEXP Csparse_is_diagonal(SEXP obj);
SEXP Rsparse_is_diagonal(SEXP obj);
SEXP Tsparse_is_diagonal(SEXP obj);
SEXP Csparse_is_triangular(SEXP obj, SEXP upper);
SEXP Rsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper);
SEXP Csparse_is_symmetric(SEXP obj, SEXP checkDN);
SEXP Rsparse_is_symmetric(SEXP obj, SEXP checkDN);

/* TODO: compare with macros in ./Mutils.h */

#define VALID_DSPARSE				\
"dgCMatrix", "dgRMatrix", "dgTMatrix",		\
"dtCMatrix", "dtRMatrix", "dtTMatrix",		\
"dsCMatrix", "dsRMatrix", "dsTMatrix"

#define VALID_LSPARSE				\
"lgCMatrix", "lgRMatrix", "lgTMatrix",		\
"ltCMatrix", "ltRMatrix", "ltTMatrix",		\
"lsCMatrix", "lsRMatrix", "lsTMatrix"

#define VALID_NSPARSE				\
"ngCMatrix", "ngRMatrix", "ngTMatrix",		\
"ntCMatrix", "ntRMatrix", "ntTMatrix",		\
"nsCMatrix", "nsRMatrix", "nsTMatrix"

#define VALID_CRSPARSE				\
"dgCMatrix", "dtCMatrix", "dsCMatrix",		\
"dgRMatrix", "dtRMatrix", "dsRMatrix",		\
"lgCMatrix", "ltCMatrix", "lsCMatrix",		\
"lgRMatrix", "ltRMatrix", "lsRMatrix",		\
"ngCMatrix", "ntCMatrix", "nsCMatrix",		\
"ngRMatrix", "ntRMatrix", "nsRMatrix"

#define VALID_CSPARSE				\
"dgCMatrix", "dtCMatrix", "dsCMatrix",		\
"lgCMatrix", "ltCMatrix", "lsCMatrix",		\
"ngCMatrix", "ntCMatrix", "nsCMatrix"

#define VALID_RSPARSE				\
"dgRMatrix", "dtRMatrix", "dsRMatrix",		\
"lgRMatrix", "ltRMatrix", "lsRMatrix",		\
"ngRMatrix", "ntRMatrix", "nsRMatrix"

#define VALID_TSPARSE				\
"dgTMatrix", "dtTMatrix", "dsTMatrix",		\
"lgTMatrix", "ltTMatrix", "lsTMatrix",		\
"ngTMatrix", "ntTMatrix", "nsTMatrix"
    
#define VALID_DIAGONAL				\
"ddiMatrix", "ldiMatrix"

#define SPARSE_CASES(_SEXPTYPE_, _DO_)			\
    do {						\
	switch (_SEXPTYPE_) {				\
	case REALSXP:					\
	    _DO_(double, REAL, 1.0);			\
	    break;					\
	case LGLSXP:					\
	    _DO_(int, LOGICAL, 1);			\
	    break;					\
	case INTSXP:					\
	    _DO_(int, INTEGER, 1);			\
	    break;					\
	case CPLXSXP:					\
	{						\
	    Rcomplex one; one.r = 1.0; one.i = 0.0;	\
	    _DO_(Rcomplex, COMPLEX, one);		\
	    break;					\
	}						\
	default:					\
	    break;					\
	}						\
    } while (0)

#define DROP0_CASES(_SEXPTYPE_, _DO_)			\
    do {						\
	switch (_SEXPTYPE_) {				\
	case REALSXP:					\
	    _DO_(double, REAL, NZ_REAL);		\
	    break;					\
	case LGLSXP:					\
	    _DO_(int, LOGICAL, NZ_INTEGER);		\
	    break;					\
	case INTSXP:					\
	    _DO_(int, INTEGER, NZ_INTEGER);		\
	    break;					\
	case CPLXSXP:					\
	    _DO_(Rcomplex, COMPLEX, NZ_COMPLEX);	\
	    break;					\
	default:					\
	    break;					\
	}						\
    } while (0)

#endif
