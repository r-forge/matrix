#ifndef MATRIX_MDEFINES_H
#define MATRIX_MDEFINES_H

#include "version.h"

#define Matrix_Domain "Matrix"
#define Matrix_CallocThreshold 8192
#define Matrix_ErrorBufferSize 4096

/* NB: system headers should come before R headers */

#ifdef __GLIBC__
/* ensure that strdup() and others are declared when string.h is included : */
# define _POSIX_C_SOURCE 200809L
#endif

#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <float.h>

#ifndef STRICT_R_HEADERS
# define STRICT_R_HEADERS
#endif

#include <R.h>
#include <Rinternals.h>

/* Copy and paste from WRE : */
#ifdef ENABLE_NLS
# include <libintl.h>
# define _(String) dgettext(Matrix_Domain, String)
#else
# define _(String) (String)
# define dngettext(Domain, String, StringP, N) ((N == 1) ? String : StringP)
#endif

/* Copy and paste from Defn.h : */
/* 'alloca' is neither C99 nor POSIX */
#ifdef __GNUC__
/* This covers GNU, Clang and Intel compilers */
/* #undef needed in case some other header, e.g. malloc.h, already did this */
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
# ifdef HAVE_ALLOCA_H
/* This covers native compilers on Solaris and AIX */
#  include <alloca.h>
# endif
/* It might have been defined via some other standard header, e.g. stdlib.h */
# if !HAVE_DECL_ALLOCA
extern void *alloca(size_t);
# endif
#endif

#define Matrix_Calloc(_VAR_, _N_, _CTYPE_) \
do { \
	if (_N_ >= Matrix_CallocThreshold) \
		_VAR_ = R_Calloc(_N_, _CTYPE_); \
	else { \
		_VAR_ = (_CTYPE_ *) alloca(sizeof(_CTYPE_) * (size_t) (_N_)); \
		R_CheckStack(); \
		memset(_VAR_, 0, sizeof(_CTYPE_) * (size_t) (_N_)); \
	} \
} while (0)

#define Matrix_Free(_VAR_, _N_) \
do { \
	if (_N_ >= Matrix_CallocThreshold) \
		R_Free(_VAR_); \
} while (0)

#define errorChar(...)   mkChar  (Matrix_sprintf(__VA_ARGS__))
#define errorString(...) mkString(Matrix_sprintf(__VA_ARGS__))

#define HAS_SLOT(x, name)        R_has_slot      (x, name)
#define GET_SLOT(x, name)        R_do_slot       (x, name)
#define SET_SLOT(x, name, value) R_do_slot_assign(x, name, value)

/* TYPEOF returns int, not SEXPTYPE (unsigned int) => -Wsign-conversion */
#define TYPEOF(s) ((SEXPTYPE) (TYPEOF)(s))

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

#define LONGDOUBLE_AS_DOUBLE(x) \
	((x > DBL_MAX) ? R_PosInf : ((x < -DBL_MAX) ? R_NegInf : (double) x))

#define    MINOF(x, y) ((x < y) ? x : y)
#define    MAXOF(x, y) ((x < y) ? y : x)
#define  FIRSTOF(x, y) (x)
#define SECONDOF(x, y) (y)

#define ISNA_PATTERN(_X_) (0)
#define ISNA_LOGICAL(_X_) ((_X_) == NA_LOGICAL)
#define ISNA_INTEGER(_X_) ((_X_) == NA_INTEGER)
#define ISNA_REAL(_X_)    (ISNAN(_X_))
#define ISNA_COMPLEX(_X_) (ISNAN((_X_).r) || ISNAN((_X_).i))

#define NOTZERO_PATTERN(_X_) ((_X_) != 0)
#define NOTZERO_LOGICAL(_X_) ((_X_) != 0)
#define NOTZERO_INTEGER(_X_) ((_X_) != 0)
#define NOTZERO_REAL(_X_)    ((_X_) != 0.0)
#define NOTZERO_COMPLEX(_X_) ((_X_).r != 0.0 || (_X_).i != 0.0)

#define STRICTLY_NOTZERO_PATTERN(_X_) \
	(                      NOTZERO_PATTERN(_X_))
#define STRICTLY_NOTZERO_LOGICAL(_X_) \
	(!ISNA_LOGICAL(_X_) && NOTZERO_LOGICAL(_X_))
#define STRICTLY_NOTZERO_INTEGER(_X_) \
	(!ISNA_INTEGER(_X_) && NOTZERO_INTEGER(_X_))
#define STRICTLY_NOTZERO_REAL(_X_) \
	(!ISNA_REAL   (_X_) && NOTZERO_REAL   (_X_))
#define STRICTLY_NOTZERO_COMPLEX(_X_) \
	(!ISNA_COMPLEX(_X_) && NOTZERO_COMPLEX(_X_))

#define NOTREAL_PATTERN(_X_) 0
#define NOTREAL_LOGICAL(_X_) 0
#define NOTREAL_INTEGER(_X_) 0
#define NOTREAL_REAL(_X_)    0
#define NOTREAL_COMPLEX(_X_) (_X_.i != 0.0)

#define NOTEQUAL_PATTERN(_X_, _Y_) \
	((_X_ != 0) != (_Y_ != 0))
#define NOTEQUAL_LOGICAL(_X_, _Y_) \
	(_X_ != _Y_)
#define NOTEQUAL_INTEGER(_X_, _Y_) \
	(_X_ != _Y_)
#define NOTEQUAL_REAL(_X_, _Y_) \
	((ISNAN(_X_)) ? !ISNAN(_Y_) : ISNAN(_Y_) || _X_ != _Y_)
#define NOTEQUAL_COMPLEX(_X_, _Y_) \
	(((ISNAN(_X_.r)) ? !ISNAN(_Y_.r) : ISNAN(_Y_.r) || _X_.r != _Y_.r) || \
	 ((ISNAN(_X_.i)) ? !ISNAN(_Y_.i) : ISNAN(_Y_.i) || _X_.i != _Y_.i))

#define NOTCONJ_PATTERN(_X_, _Y_) \
	((_X_ != 0) != (_Y_ != 0))
#define NOTCONJ_LOGICAL(_X_, _Y_) \
	(_X_ != _Y_)
#define NOTCONJ_INTEGER(_X_, _Y_) \
	(_X_ != _Y_)
#define NOTCONJ_REAL(_X_, _Y_) \
	((ISNAN(_X_)) ? !ISNAN(_Y_) : ISNAN(_Y_) || _X_ != _Y_)
#define NOTCONJ_COMPLEX(_X_, _Y_) \
	(((ISNAN(_X_.r)) ? !ISNAN(_Y_.r) : ISNAN(_Y_.r) || _X_.r !=  _Y_.r) || \
	 ((ISNAN(_X_.i)) ? !ISNAN(_Y_.i) : ISNAN(_Y_.i) || _X_.i != -_Y_.i))

#define INCREMENT_PATTERN(_X_, _Y_) \
	do { \
		_X_ = 1; \
	} while (0)
#define INCREMENT_LOGICAL(_X_, _Y_) \
	do { \
		if (_Y_ == NA_LOGICAL) { \
			if (_X_ == 0) \
				_X_ = NA_LOGICAL; \
		} else if (_Y_ != 0) \
			_X_ = 1; \
	} while (0)
#define INCREMENT_INTEGER(_X_, _Y_) \
	do { \
		if (_X_ != NA_INTEGER) { \
			if (_Y_ == NA_INTEGER) \
				_X_ = NA_INTEGER; \
			else if ((_Y_ < 0) \
					 ? (_X_ <= INT_MIN - _Y_) \
					 : (_X_ >  INT_MAX - _Y_)) { \
				warning(_("NAs produced by integer overflow")); \
				_X_ = NA_INTEGER; \
			} else \
				_X_ += _Y_; \
		} \
	} while (0)
#define INCREMENT_REAL(_X_, _Y_) \
	do { \
		_X_ += _Y_; \
	} while (0)
#define INCREMENT_COMPLEX_ID(_X_, _Y_) \
	do { \
		_X_.r += _Y_.r; \
		_X_.i += _Y_.i; \
	} while (0)
#define INCREMENT_COMPLEX_CJ(_X_, _Y_) \
	do { \
		_X_.r += _Y_.r; \
		_X_.i -= _Y_.i; \
	} while (0)

#define ASSIGN1_REAL_RE(_X_, _Y_) _X_ = _Y_
#define ASSIGN1_REAL_IM(_X_, _Y_)

#define ASSIGN1_COMPLEX_RE(_X_, _Y_) _X_.r = _Y_
#define ASSIGN1_COMPLEX_IM(_X_, _Y_) _X_.i = _Y_

#define ASSIGN2_REAL_ID(_X_, _Y_) _X_ = _Y_
#define ASSIGN2_REAL_CJ(_X_, _Y_) _X_ = _Y_
#define ASSIGN2_REAL_RE(_X_, _Y_) _X_ = _Y_
#define ASSIGN2_REAL_IM(_X_, _Y_) _X_ = 0.0

#define ASSIGN2_COMPLEX_ID(_X_, _Y_) \
	do { _X_.r = _Y_.r; _X_.i =  _Y_.i; } while (0)
#define ASSIGN2_COMPLEX_CJ(_X_, _Y_) \
	do { _X_.r = _Y_.r; _X_.i = -_Y_.i; } while (0)
#define ASSIGN2_COMPLEX_RE(_X_, _Y_) \
	do { _X_.r = _Y_.r; _X_.i =    0.0; } while (0)
#define ASSIGN2_COMPLEX_IM(_X_, _Y_) \
	do { _X_.r =   0.0; _X_.i =  _Y_.i; } while (0)

#define SCALE1_REAL(_X_, _A_) \
	do { _X_   *= _A_;               } while (0)
#define SCALE1_COMPLEX(_X_, _A_) \
	do { _X_.r *= _A_; _X_.i *= _A_; } while (0)

#define SCALE2_REAL(_X_, _A_) \
	do { _X_   /= _A_;               } while (0)
#define SCALE2_COMPLEX(_X_, _A_) \
	do { _X_.r /= _A_; _X_.i /= _A_; } while (0)

#define PACKED_AR21_UP(i, j) \
	((i) + ((j) * (            (j) + 1U)) / 2U)
#define PACKED_AR21_LO(i, j, n) \
	((i) + ((j) * ((n) - (j) - 1U + (n))) / 2U)
#define PACKED_LENGTH(n) \
	((n) + ((n) * (            (n) - 1U)) / 2U)

#define SHOW(...) __VA_ARGS__
#define HIDE(...)

#define ERROR_OOM(_FUNC_) \
	error(_("out of memory in '%s'"), \
	      _FUNC_)

#define ERROR_INVALID_TYPE(_X_, _FUNC_) \
	error(_("invalid type \"%s\" in '%s'"), \
	      type2char(TYPEOF(_X_)), _FUNC_)

#define ERROR_INVALID_CLASS(_X_, _FUNC_) \
do { \
	if (!OBJECT(_X_)) \
		ERROR_INVALID_TYPE(_X_, _FUNC_); \
	else { \
		SEXP class = PROTECT(getAttrib(_X_, R_ClassSymbol)); \
		error(_("invalid class \"%s\" in '%s'"), \
		      CHAR(STRING_ELT(class, 0)), _FUNC_); \
		UNPROTECT(1); \
	} \
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

#include "utils.h"

SEXP newObject(const char *);
void validObject(SEXP, const char *);

char typeToKind(SEXPTYPE);
SEXPTYPE kindToType(char);
size_t kindToSize(char);

const char *Matrix_class(SEXP, const char **, int, const char *);

int DimNames_is_trivial(SEXP);
int DimNames_is_symmetric(SEXP);

void symDN(SEXP, SEXP, int);
void revDN(SEXP, SEXP);

SEXP get_symmetrized_DimNames(SEXP, int);
SEXP get_reversed_DimNames(SEXP);

void set_symmetrized_DimNames(SEXP, SEXP, int);
void set_reversed_DimNames(SEXP, SEXP);

#endif /* MATRIX_MDEFINES_H */
