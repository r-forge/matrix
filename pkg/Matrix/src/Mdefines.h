#ifndef MATRIX_MDEFINES_H
#define MATRIX_MDEFINES_H

#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif

#include <float.h> /* DBL_EPSILON */
#include <limits.h> /* INT_MAX */
#include <math.h> /* fabs, hypot, sqrt, ... */
#include <stdarg.h> /* va_list, va_start, ... */
#include <stddef.h> /* size_t */
#include <stdint.h> /* int_fast64_t */
#include <stdio.h> /* vsnprintf */
#include <string.h> /* memcpy, memset */

#include <Rconfig.h> /* HAVE_ALLOCA_H, ENABLE_NLS */
#if !defined(__GNUC__) && defined(HAVE_ALLOCA_H)
# include <alloca.h> /* alloca */
#endif
#ifdef ENABLE_NLS
# include <libintl.h> /* dgettext, dngettext */
#endif

#include <R_ext/Arith.h> /* ISNAN, R_FINITE, ... */
#include <R_ext/Boolean.h> /* Rboolean */
#include <R_ext/Complex.h> /* Rcomplex */
#include <R_ext/Error.h> /* Rf_error, Rf_warning */
#include <R_ext/Memory.h> /* R_alloc */
#include <R_ext/Utils.h> /* R_CheckStack */
#include <R_ext/RS.h> /* R_Calloc, R_Free */
#include <Rinternals.h> /* SEXP, ... */
#include <Rversion.h> /* R_VERSION, ... */

#define Matrix_ErrorBufferSize   4096
#define Matrix_CallocThreshold   8192
#define Matrix_TranslationDomain "Matrix"

#ifdef __GNUC__
# define Matrix_alloca(x) __builtin_alloca((x))
#else
# define Matrix_alloca(x)           alloca((x))
#endif

#define Matrix_Calloc(p, n, t) \
do { \
	if ((n) >= Matrix_CallocThreshold) \
		(p) = R_Calloc((n), t); \
	else { \
		(p) = (t *) Matrix_alloca(sizeof(t) * (size_t) (n)); \
		R_CheckStack(); \
		memset((p), 0, sizeof(t) * (size_t) (n)); \
	} \
} while (0)

#define Matrix_Free(p, n) \
do { \
	if ((n) >= Matrix_CallocThreshold) \
		R_Free((p)); \
} while (0)

#ifndef ENABLE_NLS
# define dgettext(Domain, String) (String)
# define dngettext(Domain, String, StringP, N) (((N) == 1) ? String : StringP)
#endif
#define _(String) dgettext(Matrix_TranslationDomain, String)

#define errorChar(...)   Rf_mkChar  (Matrix_sprintf(__VA_ARGS__))
#define errorString(...) Rf_mkString(Matrix_sprintf(__VA_ARGS__))

#define HAS_SLOT(x, name)        R_has_slot      ((x), (name))
#define GET_SLOT(x, name)        R_do_slot       ((x), (name))
#define SET_SLOT(x, name, value) R_do_slot_assign((x), (name), (value))

#define COPY_SLOT(dest, src, name) \
	do { \
		SEXP value = GET_SLOT((src), (name)); \
		if (XLENGTH(value) > 0) { \
			PROTECT(value); \
			SET_SLOT((dest), (name), value); \
			UNPROTECT(1); \
		} \
	} while (0)

#define     DIM(x) \
	INTEGER(GET_SLOT((x), Matrix_DimSym))
#define SET_DIM(x, m, n) \
	do { \
		int __m__ = (m), __n__ = (n); \
		if (__m__ != __n__ || __n__ > 0) { \
			int *p = INTEGER(GET_SLOT((x), Matrix_DimSym)); \
			p[0] = (m); p[1] = (n); \
		} \
	} while (0)

#define     DIMNAMES(x, mode) \
	(DIMNAMES)((x), (mode))
#define SET_DIMNAMES(x, mode, value) \
	(SET_DIMNAMES)((x), (mode), (value))

#define     UPLO(x) \
	CHAR(STRING_ELT(GET_SLOT((x), Matrix_uploSym), 0))[0]
#define SET_UPLO(x) \
	SET_STRING_ELT(GET_SLOT((x), Matrix_uploSym), 0, Matrix_LChar)

#define     TRANS(x) \
	CHAR(STRING_ELT(GET_SLOT((x), Matrix_transSym), 0))[0]
#define SET_TRANS(x) \
	SET_STRING_ELT(GET_SLOT((x), Matrix_transSym), 0, Matrix_TChar)

#define     DIAG(x) \
	CHAR(STRING_ELT(GET_SLOT((x), Matrix_diagSym), 0))[0]
#define SET_DIAG(x) \
	SET_STRING_ELT(GET_SLOT((x), Matrix_diagSym), 0, Matrix_UChar)

#define     MARGIN(x) \
	(INTEGER(GET_SLOT((x), Matrix_marginSym))[0] - 1)
#define SET_MARGIN(x, j) \
	do { \
		INTEGER(GET_SLOT((x), Matrix_marginSym))[0] = (j) + 1; \
	} while (0)

#define TYPEOF(s) \
	((SEXPTYPE) (TYPEOF)((s)))

#define DENSE_INDEX_N(i, j, m) \
	((i) + (j) * (m))
#define DENSE_INDEX_U(i, j, m) \
	((i) + ((j) * (            (j) + 1U)) / 2U)
#define DENSE_INDEX_L(i, j, m) \
	((i) + ((j) * ((m) - (j) - 1U + (m))) / 2U)
#define PACKED_LENGTH(n) \
	((n) + ((n) * (            (n) - 1U)) / 2U)

#define ABS(i) \
	(((i) < 0) ? -(i) : (i))

#define SWAP(a, b, t, op) \
	do { t tmp = op(a); a = op(b); b = tmp; } while (0)

#define ERROR_OOM(_FUNC_) \
	Rf_error(_("out of memory in '%s'"), \
	         _FUNC_)

#define ERROR_INVALID_TYPE(_X_, _FUNC_) \
	Rf_error(_("invalid type \"%s\" in '%s'"), \
	         Rf_type2char(TYPEOF(_X_)), _FUNC_)

#define ERROR_INVALID_CLASS(_X_, _FUNC_) \
do { \
	if (!Rf_isObject(_X_)) \
		ERROR_INVALID_TYPE(_X_, _FUNC_); \
	else \
		Rf_error(_("invalid class \"%s\" in '%s'"), \
		         CHAR(STRING_ELT(Rf_getAttrib(_X_, R_ClassSymbol), 0)), _FUNC_); \
} while (0)

#define VALID_UPLO(s, c) \
do { \
	if (TYPEOF(s) != STRSXP || LENGTH(s) < 1 || \
		(s = STRING_ELT(s, 0)) == NA_STRING || \
		((c = CHAR(s)[0]) != 'U' && c != 'L')) \
		Rf_error(_("'%s' must be \"%c\" or \"%c\""),  "uplo", 'U', 'L'); \
} while (0)

#define VALID_TRANS(s, c) \
do { \
	if (TYPEOF(s) != STRSXP || LENGTH(s) < 1 || \
	    (s = STRING_ELT(s, 0)) == NA_STRING || \
	    ((c = CHAR(s)[0]) != 'C' && c != 'T')) \
		Rf_error(_("'%s' is not \"%c\" or \"%c\""), "trans", 'C', 'T'); \
} while (0)

#define VALID_DIAG(s, c) \
do { \
	if (TYPEOF(s) != STRSXP || LENGTH(s) < 1 || \
	    (s = STRING_ELT(s, 0)) == NA_STRING || \
	    ((c = CHAR(s)[0]) != 'N' && c != 'U')) \
		Rf_error(_("'%s' is not \"%c\" or \"%c\""),  "diag", 'N', 'U'); \
} while (0)

#define VALID_KIND(s, c) \
do { \
	if (TYPEOF(s) != STRSXP || LENGTH(s) < 1 || \
	    (s = STRING_ELT(s, 0)) == NA_STRING || \
	    ((c = CHAR(s)[0]) != 'n' && c != 'l' && c != 'i' && c != 'd' && c != 'z' && c != '.' && c != ',')) \
		Rf_error(_("'%s' is not \"%c\", \"%c\", \"%c\", \"%c\", or \"%c\""), \
		         "kind", 'n', 'l', 'i', 'd', 'z'); \
} while (0)

#define VALID_SHAPE(s, c) \
do { \
	if (TYPEOF(s) != STRSXP || LENGTH(s) < 1 || \
	    (s = STRING_ELT(s, 0)) == NA_STRING || \
	    ((c = CHAR(s)[0]) != 'g' && c != 's' && c != 'p' && c != 't')) \
		Rf_error(_("'%s' is not \"%c\", \"%c\", \"%c\", or \"%c\""), \
		         "shape", 'g', 's', 'p', 't'); \
} while (0)

#define VALID_REPR(s, c, dot) \
do { \
	if (TYPEOF(s) != STRSXP || LENGTH(s) < 1 || \
	    (s = STRING_ELT(s, 0)) == NA_STRING || \
	    ((c = CHAR(s)[0]) != 'C' && c != 'R' && c != 'T' && !(dot && c != '.'))) \
		Rf_error(_("'%s' is not \"%c\", \"%c\", or \"%c\""), \
		         "repr", 'C', 'R', 'T'); \
} while (0)

#define VALID_MARGIN(s, d) \
do { \
	if (TYPEOF(s) != INTSXP || LENGTH(s) < 1 || \
	    ((d = INTEGER(s)[0] - 1) != 0 && d != 1)) \
		Rf_error(_("'%s' is not %d or %d"), "margin", 1, 2); \
} while (0)

#define VALID_LOGIC2(s, d) \
do { \
	if (TYPEOF(s) != LGLSXP || LENGTH(s) < 1 || \
	    ((d = LOGICAL(s)[0]) == NA_LOGICAL)) \
		Rf_error(_("'%s' is not %s or %s"), #d, "TRUE", "FALSE"); \
} while (0)

#define VALID_LOGIC3(s, d) \
do { \
	if (TYPEOF(s) != LGLSXP || LENGTH(s) < 1) \
		Rf_error(_("'%s' is not %s, %s, or %s"), #d, "TRUE", "FALSE", "NA"); \
	else d = LOGICAL(s)[0]; \
} while (0)

#define VALID_DENSE \
"ngeMatrix", "lgeMatrix", "igeMatrix", "dgeMatrix", "zgeMatrix", \
"nsyMatrix", "lsyMatrix", "isyMatrix", "dsyMatrix", "zsyMatrix", \
                                       "dpoMatrix", "zpoMatrix", \
                                       "corMatrix",              \
"ntrMatrix", "ltrMatrix", "itrMatrix", "dtrMatrix", "ztrMatrix", \
"nspMatrix", "lspMatrix", "ispMatrix", "dspMatrix", "zspMatrix", \
                                       "dppMatrix", "zppMatrix", \
                                       "copMatrix",              \
"ntpMatrix", "ltpMatrix", "itpMatrix", "dtpMatrix", "ztpMatrix"

#define VALID_SPARSE_COMPRESSED \
"ngCMatrix", "lgCMatrix", "igCMatrix", "dgCMatrix", "zgCMatrix", \
"nsCMatrix", "lsCMatrix", "isCMatrix", "dsCMatrix", "zsCMatrix", \
                                       "dpCMatrix", "zpCMatrix", \
"ntCMatrix", "ltCMatrix", "itCMatrix", "dtCMatrix", "ztCMatrix", \
"ngRMatrix", "lgRMatrix", "igRMatrix", "dgRMatrix", "zgRMatrix", \
"nsRMatrix", "lsRMatrix", "isRMatrix", "dsRMatrix", "zsRMatrix", \
                                       "dpRMatrix", "zpRMatrix", \
"ntRMatrix", "ltRMatrix", "itRMatrix", "dtRMatrix", "ztRMatrix"

#define VALID_SPARSE_TRIPLET \
"ngTMatrix", "lgTMatrix", "igTMatrix", "dgTMatrix", "zgTMatrix", \
"nsTMatrix", "lsTMatrix", "isTMatrix", "dsTMatrix", "zsTMatrix", \
                                       "dpTMatrix", "zpTMatrix", \
"ntTMatrix", "ltTMatrix", "itTMatrix", "dtTMatrix", "ztTMatrix"

#define VALID_SPARSE \
VALID_SPARSE_COMPRESSED, VALID_SPARSE_TRIPLET

#define VALID_DIAGONAL \
"ndiMatrix", "ldiMatrix", "idiMatrix", "ddiMatrix", "zdiMatrix"

#define VALID_INDEX \
"indMatrix",   "pMatrix"

#define VALID_MATRIX \
VALID_DENSE, VALID_SPARSE, VALID_DIAGONAL, VALID_INDEX

#define VALID_VECTOR \
"nsparseVector", "lsparseVector", "isparseVector", "dsparseVector", "zsparseVector"

#define VALID_MATRIX_OR_VECTOR \
VALID_MATRIX, VALID_VECTOR

/* What we want declared "everywhere" : */

/* Often used symbols, declared in ./init.c */
extern
#include "Msymbols.h"

/* Often used numbers, declared in ./init.c */
extern
Rcomplex Matrix_zzero, Matrix_zunit, Matrix_zna; /* 0+0i, 1+0i, NA+NAi */

/* Often used arrays of nonvirtual class names, declared in ./objects.c */
extern
const char *valid_dense[], *valid_sparse[],
	*valid_sparse_compressed[], *valid_sparse_triplet[],
	*valid_diagonal[], *valid_index[],
	*valid_matrix[], *valid_vector[], *valid_matrix_or_vector[];

#if R_VERSION < R_Version(4, 5, 0)
int ANY_ATTRIB(SEXP);
void CLEAR_ATTRIB(SEXP);
#endif

char *Matrix_sprintf(const char *, ...);

int equalString(SEXP, SEXP, R_xlen_t);
SEXP duplicateVector(SEXP);
SEXP allocZero(SEXPTYPE, R_xlen_t);
SEXP allocUnit(SEXPTYPE, R_xlen_t);
SEXP allocSeqInt(int, R_xlen_t);
void naToUnit(SEXP);

SEXP newObject(const char *);
void validObject(SEXP, const char *);

char typeToKind(SEXPTYPE);
SEXPTYPE kindToType(char);
size_t kindToSize(char);

const char *Matrix_superclass(const char *, int);
const char *Matrix_class(SEXP, const char **, int, const char *);

int DimNames_is_trivial(SEXP);
int DimNames_is_symmetric(SEXP);

void symDN(SEXP, SEXP, int);
void cpyDN(SEXP, SEXP, int);

SEXP (DIMNAMES)(SEXP, int);
void (SET_DIMNAMES)(SEXP, int, SEXP);

#endif /* MATRIX_MDEFINES_H */
