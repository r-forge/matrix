#include "packedMatrix.h"

#define UNPACK(_PREFIX_, _CTYPE_, _ONE_)				\
void _PREFIX_ ## dense_unpack(_CTYPE_ *dest, const _CTYPE_ *src, int n, \
			      char uplo, char diag)			\
{									\
    int i, j;								\
    R_xlen_t dpos = 0, spos = 0;					\
    if (uplo == 'U') {							\
	for (j = 0; j < n; dpos += n-(++j))				\
	    for (i = 0; i <= j; ++i)					\
		dest[dpos++] = src[spos++];				\
    } else {								\
	for (j = 0; j < n; dpos += (++j))				\
	    for (i = j; i <  n; ++i)					\
		dest[dpos++] = src[spos++];				\
    }									\
    if (diag != 'N') {							\
	dpos = 0;							\
	R_xlen_t n1a = (R_xlen_t) n + 1;				\
	for (j = 0; j < n; ++j, dpos += n1a)				\
	    dest[dpos] = _ONE_;						\
    }									\
    return;								\
}

/**
 * @brief Unpack a `packedMatrix`.
 * 
 * Copies `src` to the upper or lower triangular part of `dest`, 
 * where it is stored _non_-contiguously ("unpacked"). Optionally 
 * resets the diagonal elements to 1.
 * 
 * @param dest,src Pointers to the first elements of length-`n*n` and
 *     length-`(n*(n+1))/2` (resp.) arrays, usually the "data" of the 
 *     `x` slot of an `n`-by-`n` `unpackedMatrix` and `packedMatrix` 
 *     (resp.).
 * @param n Size of matrix being unpacked.
 * @param uplo,diag `char` specifying whether to copy `src` to
 *     the upper (`'U'`) or lower (`'L'`) triangle of `dest` and
 *     whether to "force" a unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_unpack() */
UNPACK(d, double, 1.0)
/* idense_unpack() */
UNPACK(i, int, 1)
/* zdense_unpack() */
UNPACK(z, Rcomplex, Matrix_zone)

#undef UNPACK

#define MAKE_BANDED(_PREFIX_, _CTYPE_, _ZERO_, _ONE_)			\
void _PREFIX_ ## dense_packed_make_banded(_CTYPE_ *x, int n,		\
					  int a, int b,	char uplo, char diag) \
{									\
    if (n == 0)								\
	return;								\
    if (a > b || a >= n || b <= -n) {					\
	Memzero(x, PM_LENGTH(n));					\
	return;								\
    }									\
    if (uplo == 'U') {							\
	if (a <   0) a = 0;						\
	if (b >=  n) b = n-1;						\
    } else {								\
	if (b >   0) b = 0;						\
	if (a <= -n) a = 1-n;						\
    }									\
    									\
    int i, j, i0, i1,							\
	j0 = (a < 0) ? 0 : a,						\
	j1 = (b < 0) ? n+b : n;						\
    									\
    if (uplo == 'U') {							\
	if (j0 > 0) {							\
	    R_xlen_t dx;						\
	    Memzero(x, dx = PM_LENGTH(j0));				\
	    x += dx;							\
	}								\
	for (j = j0; j < j1; x += (++j)) {				\
	    i0 = j - b;							\
	    i1 = j - a + 1;						\
	    for (i = 0; i < i0; ++i)					\
		*(x + i) = _ZERO_;					\
	    for (i = i1; i <= j; ++i)					\
		*(x + i) = _ZERO_;					\
	}								\
	if (j1 < n)							\
	    Memzero(x, PM_LENGTH(n) - PM_LENGTH(j1));			\
	if (diag != 'N' && a == 0) {					\
	    x -= PM_LENGTH(j);						\
	    for (j = 0; j < n; x += (++j)+1)				\
		*x = _ONE_;						\
	}								\
    } else {								\
	if (j0 > 0) {							\
	    R_xlen_t dx;						\
	    Memzero(x, dx = PM_LENGTH(n) - PM_LENGTH(j0));		\
	    x += dx;							\
	}								\
	for (j = j0; j < j1; x += n-(j++)) {				\
	    i0 = j - b;							\
	    i1 = j - a + 1;						\
	    for (i = j; i < i0; ++i)					\
		*(x + i - j) = _ZERO_;					\
	    for (i = i1; i < n; ++i)					\
		*(x + i - j) = _ZERO_;					\
	}								\
	if (j1 < n)							\
	    Memzero(x, PM_LENGTH(n - j1));				\
	if (diag != 'N' && b == 0) {					\
	    x -= PM_LENGTH(n) - PM_LENGTH(j);				\
	    for (j = 0; j < n; x += n-(j++))				\
		*x = _ONE_;						\
	}								\
    }									\
    return;								\
}
MAKE_BANDED(d, double, 0.0, 1.0)
MAKE_BANDED(i, int, 0, 1)
MAKE_BANDED(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef MAKE_BANDED

#define COPY_DIAGONAL(_PREFIX_, _CTYPE_, _ONE_)				\
void _PREFIX_ ## dense_packed_copy_diagonal(_CTYPE_ *dest,		\
					    const _CTYPE_ *src,		\
					    int n, R_xlen_t len,	\
					    char uplo_dest,		\
					    char uplo_src,		\
					    char diag)			\
{									\
    int j;								\
    if (diag != 'N') {							\
	if (uplo_dest != 'L') {						\
	    for (j = 0; j < n; dest += (++j)+1)				\
		*dest = _ONE_;						\
	} else {							\
	    for (j = 0; j < n; dest += n-(j++))				\
		*dest = _ONE_;						\
	}								\
    } else if (len == n) {						\
	/* copying from diagonalMatrix */				\
	if (uplo_dest != 'L') {						\
	    for (j = 0; j < n; dest += (++j)+1, ++src)			\
		*dest = *src;						\
	} else {							\
	    for (j = 0; j < n; dest += n-(j++), ++src)			\
		*dest = *src;						\
	}								\
    } else if (len == PM_LENGTH(n)) {					\
	/* copying from packedMatrix */					\
	if (uplo_dest != 'L') {						\
       	    if (uplo_src != 'L') {					\
	    	for (j = 0; j < n; src += (++j)+1, dest += j+1)		\
		    *dest = *src;					\
	    } else {							\
		for (j = 0; j < n; src += n-j, dest += (++j)+1)		\
		    *dest = *src;					\
	    }								\
	} else {							\
	    if (uplo_src != 'L') {					\
	    	for (j = 0; j < n; dest += n-(j++), src += j+1)		\
		    *dest = *src;					\
	    } else {							\
		for (j = 0; j < n; dest += n-j, src += n-(j++))		\
		    *dest = *src;					\
	    }								\
	}								\
    } else if (len == (R_xlen_t) n * n) {				\
	/* copying from square unpackedMatrix */			\
	R_xlen_t n1a = (R_xlen_t) n + 1;				\
	if (uplo_dest != 'L') {						\
	    for (j = 0; j < n; dest += (++j)+1, src += n1a)		\
		*dest = *src;						\
	} else {							\
	    for (j = 0; j < n; dest += n-(j++), src += n1a)		\
		*dest = *src;						\
	}								\
    } else {								\
	error(_("incompatible 'n' and 'len' to '*_copy_diagonal()'"));	\
    }									\
    return;								\
}

/**
 * Copy a length-`n` diagonal to a length-`(n*(n+1))/2` array.
 * 
 * @param dest A pointer to the first element of a length-`(n*(n+1))/2` array,
 *     usually the "data" of the `x` slot of an `n`-by-`n` `packedMatrix`.
 * @param src A pointer to the first element of a length-`n`, 
 *     length-`(n*(n+1))/2`, or length-`n*n` array, usually the "data"
 *     of the `x` slot of an `n`-by-`n` `diagonalMatrix`, `packedMatrix`,
 *     or `unpackedMatrix`, respectively.
 * @param n Size of matrix being copied from and to.
 * @param len Length of `src` array.
 * @param uplo_dest,uplo_src,diag_src `char` constants specifying 
 *     whether `dest` stores an upper (`'U'`) or lower (`'L'`) triangle, 
 *     whether `src` stores an upper (`'U'`) or lower (`'L'`) triangle 
 *     when `len == (n*(n+1))/2`, and whether the matrix should have a
 *     unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_packed_copy_diagonal() */
COPY_DIAGONAL(d, double, 1.0)
/* idense_packed_copy_diagonal() */
COPY_DIAGONAL(i, int, 1)
/* zdense_packed_copy_diagonal() */
COPY_DIAGONAL(z, Rcomplex, Matrix_zone)

#undef COPY_DIAGONAL

#define IS_DIAGONAL(_PREFIX_, _CTYPE_,					\
		    _U_IS_NOT_ZERO_, _L_IS_NOT_ZERO_)			\
static Rboolean _PREFIX_ ## dense_packed_is_diagonal(const _CTYPE_ *x,	\
						     int n, char uplo)	\
{									\
    int i, j;								\
    if (uplo == 'U') {							\
    	for (j = 0; j < n; ++j, ++x)					\
	    for (i = 0; i < j; ++i)					\
		if (_U_IS_NOT_ZERO_)					\
		    return FALSE;					\
    } else {								\
	for (j = 0; j < n; ++j, ++x)					\
	    for (i = j+1; i < n; ++i)					\
		if (_L_IS_NOT_ZERO_)					\
		    return FALSE;					\
    }									\
    return TRUE;							\
}

/* ddense_packed_is_diagonal() */
IS_DIAGONAL(d, double,
		   ISNAN(*x) || *(x++) != 0.0,
		   ISNAN(*(++x)) || *x != 0.0)
/* idense_packed_is_diagonal() */
IS_DIAGONAL(i, int,
		   *(x++) != 0,
		   *(++x) != 0)
/* zdense_packed_is_diagonal() */
IS_DIAGONAL(z, Rcomplex,
		   ISNAN((*x).r) || (*x).r != 0.0 ||
		   ISNAN((*x).i) || (*(x++)).i != 0.0,
		   ISNAN((*(++x)).r) || (*x).r != 0.0 ||
		   ISNAN((*x).i) || (*x).i != 0.0)

#undef IS_DIAGONAL

#define TRANSPOSE(_PREFIX_, _CTYPE_)					\
static void _PREFIX_ ## dense_packed_transpose(_CTYPE_ *dest,		\
					       const _CTYPE_ *src,	\
					       int n, char uplo)	\
{									\
    int i, j;								\
    if (uplo == 'U') {							\
	for (j = 0; j < n; ++j)						\
	    for (i = j; i < n; ++i)					\
		*(dest++) = *(src + PM_AR21_UP(j, i));			\
    } else {								\
	R_xlen_t n2 = (R_xlen_t) n * 2;					\
	for (j = 0; j < n; ++j)						\
	    for (i = 0; i <= j; ++i)					\
		*(dest++) = *(src + PM_AR21_LO(j, i, n2));		\
    }									\
    return;								\
}

/* ddense_packed_transpose() */
TRANSPOSE(d, double)
/* idense_packed_transpose() */
TRANSPOSE(i, int)
/* zdense_packed_transpose() */
TRANSPOSE(z, Rcomplex)

#undef TRANSPOSE

SEXP packed_transpose(SEXP x, int n, char uplo)
{
    SEXPTYPE tx = TYPEOF(x);
    if (tx < LGLSXP || tx > CPLXSXP)
	ERROR_INVALID_TYPE("'x'", tx, "packed_transpose");
    R_xlen_t nx = XLENGTH(x);
    SEXP y = PROTECT(allocVector(tx, nx));

#define TRANSPOSE(_PREFIX_, _PTR_)					\
    _PREFIX_ ## dense_packed_transpose(_PTR_(y), _PTR_(x), n, uplo)

    switch (tx) {
    case LGLSXP:
	TRANSPOSE(i, LOGICAL);
	break;
    case INTSXP:
	TRANSPOSE(i, INTEGER);
	break;
    case REALSXP:
	TRANSPOSE(d, REAL);
	break;
    case CPLXSXP:
	TRANSPOSE(z, COMPLEX);
	break;
    default:
	break;
    }

#undef TRANSPOSE
    
    UNPROTECT(1);
    return y;
}

/* unpack(x), returning unpackedMatrix */
SEXP packedMatrix_unpack(SEXP from, SEXP strict)
{
    static const char *valid_from[] = {
	/* 0 */ "pCholesky", "pBunchKaufman", /* must match before dtpMatrix */
	/* 2 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 5 */ "dppMatrix", /* must match before dspMatrix */
	/* 6 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    static const char *valid_to[] = {
	/* 0 */ "Cholesky", "BunchKaufman",
	/* 2 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/* 5 */ "dpoMatrix",
	/* 6 */ "dsyMatrix", "lsyMatrix", "nsyMatrix", ""};
    int ivalid = R_check_class_etc(from, valid_from);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(from, "packedMatrix_unpack");
    if (asLogical(strict) == 0) {
	if (ivalid < 2)
	    ivalid = 2; /* pCholesky,pBunchKaufman->dtrMatrix */
	else if (ivalid == 5)
	    ivalid = 6; /* dppMatrix->dsyMatrix */
    }

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid_to[ivalid]));

    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
    if (n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */

    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */

    SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    if (ul != 'U')
	SET_SLOT(to, Matrix_uploSym, uplo);
    UNPROTECT(1); /* uplo */

    if (ivalid < 5) {
	/* .tpMatrix */
	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (di != 'N')
	    SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */
	if (ivalid == 1) {
	    /* pBunchKaufman */
	    SEXP perm = PROTECT(GET_SLOT(from, Matrix_permSym));
	    SET_SLOT(to, Matrix_permSym, perm);
	    UNPROTECT(1); /* perm */
	}
    } else {
	/* .spMatrix */
	SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
	if (LENGTH(factors) > 0)
	    SET_SLOT(to, Matrix_factorSym, factors);
	UNPROTECT(1); /* factors */
    }
    
    SEXPTYPE tx;
    R_xlen_t nx = (R_xlen_t) n * n;
    SEXP x_from = PROTECT(GET_SLOT(from, Matrix_xSym)),
	x_to = PROTECT(allocVector(tx = TYPEOF(x_from), nx));
    
#define UNPACK(_PREFIX_, _PTR_)						\
    do {								\
	Memzero(_PTR_(x_to), nx);					\
	_PREFIX_ ## dense_unpack(_PTR_(x_to), _PTR_(x_from), n, ul, 'N'); \
    } while (0)
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	UNPACK(d, REAL);
	break; 
    case LGLSXP: /* [ln]..Matrix */
	UNPACK(i, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	UNPACK(i, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
	UNPACK(z, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_unpack");
	break;
    }

#undef UNPACK
    
    SET_SLOT(to, Matrix_xSym, x_to);

    UNPROTECT(3); /* x_to, x_from, to */
    return to;
}

/* forceSymmetric(x, uplo), returning .spMatrix */
SEXP packedMatrix_force_symmetric(SEXP from, SEXP uplo_to)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(from, "packedMatrix_force_symmetric");
    const char *clf = valid[ivalid];

    SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
    char ulf = *CHAR(STRING_ELT(uplo_from, 0)), ult = ulf;
    UNPROTECT(1); /* uplo_from */

    if (!isNull(uplo_to) &&
	(TYPEOF(uplo_to) != STRSXP || LENGTH(uplo_to) < 1 ||
	 (uplo_to = STRING_ELT(uplo_to, 0)) == NA_STRING ||
	 ((ult = *CHAR(uplo_to)) != 'U' && ult != 'L')))
	error(_("invalid 'uplo' to 'packedMatrix_force_symmetric()'"));
    
    if (clf[1] == 's') {
	/* .spMatrix */
	if (ulf == ult)
	    return from;
	SEXP to = PROTECT(packedMatrix_transpose(from));
	if (clf[0] == 'z') {
	    /* Need _conjugate_ transpose */
	    SEXP x_to = PROTECT(GET_SLOT(to, Matrix_xSym));
	    conjugate(x_to);
	    UNPROTECT(1); /* x_to */
	}
	UNPROTECT(1); /* to */
	return to;
    }
    
    /* Now handling just .tpMatrix ... */
    
    char clt[] = ".spMatrix";
    clt[0] = clf[0];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	x_from = PROTECT(GET_SLOT(from, Matrix_xSym));
    
    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    if (n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */
    
    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    set_symmetrized_DimNames(to, dimnames, -1);
    UNPROTECT(1); /* dimnames */

    if (ult != 'U') {
	PROTECT(uplo_to = mkString("L"));
	SET_SLOT(to, Matrix_uploSym, uplo_to);
	UNPROTECT(1); /* uplo_to */
    }
    
    if (ulf == ult) {
	/* .tpMatrix with correct uplo */
	SET_SLOT(to, Matrix_xSym, x_from);
    } else {
	/* .tpMatrix with incorrect uplo */
	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	
	SEXPTYPE tx = TYPEOF(x_from);
	R_xlen_t nx = XLENGTH(x_from);
	SEXP x_to = PROTECT(allocVector(tx, nx));
	
#define COPY_DIAGONAL(_PREFIX_, _PTR_)					\
	do {								\
	    Memzero(_PTR_(x_to), nx);					\
	    _PREFIX_ ## dense_packed_copy_diagonal(			\
		_PTR_(x_to), _PTR_(x_from), n, nx, ult, ulf, di);	\
	} while (0)
	
	switch (tx) {
	case REALSXP: /* d..Matrix */
	    COPY_DIAGONAL(d, REAL);
	    break;
	case LGLSXP: /* [ln]..Matrix */
	    COPY_DIAGONAL(i, LOGICAL);
	    break;
	case INTSXP: /* i..Matrix */
	    COPY_DIAGONAL(i, INTEGER);
	    break;
	case CPLXSXP: /* z..Matrix */
	    COPY_DIAGONAL(z, COMPLEX);
	    break;
	default:
	    ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_force_symmetric");
	    break;
	}

#undef COPY_DIAGONAL
	
	SET_SLOT(to, Matrix_xSym, x_to);
	UNPROTECT(1); /* x_to */
    }
    
    UNPROTECT(2); /* x_from, to */
    return to;
}

#define PM_IS_DI(_RES_, _X_, _N_, _UPLO_, _METHOD_)			\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_packed_is_diagonal(REAL(_X_), _N_, _UPLO_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_packed_is_diagonal(LOGICAL(_X_), _N_, _UPLO_); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_packed_is_diagonal(INTEGER(_X_), _N_, _UPLO_); \
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_packed_is_diagonal(COMPLEX(_X_), _N_, _UPLO_); \
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", TYPEOF(_X_), _METHOD_);	\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

/* isTriangular(x, upper) */
SEXP packedMatrix_is_triangular(SEXP obj, SEXP upper)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(obj, "packedMatrix_is_triangular");

    int need_upper = asLogical(upper);

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
#define IF_DIAGONAL							\
    Rboolean res = FALSE;						\
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),			\
	dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));			\
    int n = INTEGER(dim)[0];						\
    PM_IS_DI(res, x, n, ul, "packedMatrix_is_triangular");		\
    UNPROTECT(2); /* dim, x */						\
    if (res)
    
    if (ivalid < 3) {
	/* .tpMatrix: fast if 'upper', 'uplo' agree; else need diagonal */
	if (need_upper == NA_LOGICAL)
	    RETURN_TRUE_OF_KIND((ul == 'U') ? "U" : "L");
	else if ((need_upper) ? ul == 'U' : ul != 'U')
	    return ScalarLogical(1);
	else {
	    IF_DIAGONAL {
		return ScalarLogical(1);
	    }
	}
    } else {
	/* .spMatrix: triangular iff diagonal */
	IF_DIAGONAL {
	    if (need_upper == NA_LOGICAL)
		RETURN_TRUE_OF_KIND("U");
	    else
		return ScalarLogical(1);
	}
    }
    
#undef IF_DIAGONAL

    return ScalarLogical(0);
}

/* isSymmetric(x, tol = 0, checkDN) */
/* FIXME: not checking for real diagonal in complex case */
SEXP packedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0) {
	ERROR_INVALID_CLASS(obj, "packedMatrix_is_symmetric");
	return R_NilValue;
    } else if (ivalid < 3) {
	/* .tpMatrix: symmetric iff diagonal */
	if (asLogical(checkDN) != 0) {
	    SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	    int s = DimNames_is_symmetric(dimnames);
	    UNPROTECT(1); /* dimnames */
	    if (!s)
		return ScalarLogical(0);
	}
	Rboolean res = FALSE;
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	    dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	    uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int n = INTEGER(dim)[0];
	char ul = *CHAR(STRING_ELT(uplo, 0));
	PM_IS_DI(res, x, n, ul, "packedMatrix_is_symmetric");
	UNPROTECT(3); /* uplo, dim, x */
	return ScalarLogical(res);
    } else {
	/* .spMatrix: symmetric by definition */
	return ScalarLogical(1);
    }
}

/* isDiagonal(x) */
SEXP packedMatrix_is_diagonal(SEXP obj)
{
    /* _Not_ checking class of 'obj' */
    Rboolean res = FALSE;
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int n = INTEGER(dim)[0];
    char ul = *CHAR(STRING_ELT(uplo, 0));
    PM_IS_DI(res, x, n, ul, "packedMatrix_is_diagonal");
    UNPROTECT(3); /* uplo, dim, x */
    return ScalarLogical(res);
}

#undef PM_IS_DI

/* t(x), typically preserving class */
SEXP packedMatrix_transpose(SEXP from)
{
    static const char *valid[] = {
	/* 0 */ "pCholesky", "pBunchKaufman",
	/* 2 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 5 */ "dppMatrix",
	/* 6 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(from, "packedMatrix_transpose");
    if (ivalid == 1)
	ivalid = 2; /* pBunchKaufman->dtpMatrix */

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid[ivalid]));
    
    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)); 
    int n = INTEGER(dim)[0];
    if (n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */

    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    if (ivalid < 5)
	set_reversed_DimNames(to, dimnames);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */
    
    SEXP uplo_from = PROTECT(GET_SLOT(from, Matrix_uploSym));
    char ulf = *CHAR(STRING_ELT(uplo_from, 0));
    UNPROTECT(1); /* uplo_from */

    if (ulf == 'U') {
	SEXP uplo_to = PROTECT(mkString("L"));
	SET_SLOT(to, Matrix_uploSym, uplo_to);
	UNPROTECT(1); /* uplo_to */
    }
    
    if (ivalid < 5) {
	/* .tpMatrix */
	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (di != 'N')
	    SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */
    } else {
	/* .spMatrix */
	SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
	if (LENGTH(factors) > 0)
	    SET_SLOT(to, Matrix_factorSym, factors);
	UNPROTECT(1); /* factors */
    }
    
    SEXP x_from = PROTECT(GET_SLOT(from, Matrix_xSym)),
	x_to = PROTECT(packed_transpose(x_from, n, ulf));
    SET_SLOT(to, Matrix_xSym, x_to);
    
    UNPROTECT(3); /* x_to, x_from, to */
    return to;
}

/* diag(x, names) */
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
	error(_("'names' must be TRUE or FALSE"));

    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    
    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */

    char di = '\0';
    if (HAS_SLOT(obj, Matrix_diagSym)) {
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
    }

    SEXPTYPE tx;
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	res = PROTECT(allocVector(tx = TYPEOF(x), n));

#define PM_D_G(_CTYPE_, _PTR_, _ONE_)					\
    do {								\
	_CTYPE_ *pres = _PTR_(res);					\
	int j;								\
	if (di == 'U') {						\
	    for (j = 0; j < n; ++j)					\
		*(pres++) = _ONE_;					\
	} else {							\
	    _CTYPE_ *px = _PTR_(x);					\
	    if (ul == 'U')						\
		for (j = 0; j < n; px += (++j)+1)			\
		    *(pres++) = *px;					\
	    else							\
		for (j = 0; j < n; px += n-(j++))			\
		    *(pres++) = *px;					\
	}								\
    } while (0)

    switch (tx) {
    case REALSXP: /* d..Matrix */
	PM_D_G(double, REAL, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_D_G(int, LOGICAL, 1);
	break;
    case INTSXP: /* i..Matrix */
	PM_D_G(int, INTEGER, 1);
	break;
    case CPLXSXP: /* z..Matrix */
	PM_D_G(Rcomplex, COMPLEX, Matrix_zone);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_diag_get");
	break;
    }
    
#undef PM_D_G
    
    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	if (isNull(cn)) {
	    if (di == '\0' && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (di == '\0' || (!isNull(rn) &&
			       (rn == cn || equal_string_vectors(rn, cn, n))))
		setAttrib(res, R_NamesSymbol, cn);
	}
	UNPROTECT(1); /* dn */
    }
    
    UNPROTECT(2); /* res, x */
    return res;
}

/* diag(x) <- value */
SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(obj, "packedMatrix_diag_set");

    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    
    SEXPTYPE tv = TYPEOF(val);
    if (tv < LGLSXP || tv > REALSXP)
	/* Upper bound can become CPLXSXP once we have proper zMatrix */
	error(_("replacement diagonal has incompatible type \"%s\""),
	      type2char(tv));
    
    R_xlen_t nv = XLENGTH(val);
    if (nv != 1 && nv != n)
	error(_("replacement diagonal has wrong length"));
    
    SEXP x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pid);
    SEXPTYPE tx = TYPEOF(x);
    
    /* Allocate and coerce as necessary */
    SEXP res;
    if (tv <= tx) {
	PROTECT(val = coerceVector(val, tv = tx));
	PROTECT(res = NEW_OBJECT_OF_CLASS(valid[ivalid]));
	REPROTECT(x = duplicate(x), pid);
    } else { /* tv > tx */
	/* dMatrix result is only possibility until we have proper [iz]Matrix */
	PROTECT(val = coerceVector(val, tv = REALSXP));
	char cl[] = "d.pMatrix";
	cl[1] = valid[ivalid][1];
	PROTECT(res = NEW_OBJECT_OF_CLASS(cl));
	REPROTECT(x = coerceVector(x, tx = tv), pid);
    }

    if (n > 0)
	SET_SLOT(res, Matrix_DimSym, dim);

    SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
    SET_SLOT(res, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */
    
    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    if (ul != 'U')
	SET_SLOT(res, Matrix_uploSym, uplo);
    UNPROTECT(1); /* uplo */
    
#define PM_D_S(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px = _PTR_(x), *pval = _PTR_(val);			\
	int j;								\
	if (nv == 1) {							\
	    if (ul == 'U')						\
		for (j = 0; j < n; px += (++j)+1)			\
		    *px = *pval;					\
	    else							\
		for (j = 0; j < n; px += n-(j++))			\
		    *px = *pval;					\
	} else {							\
	    if (ul == 'U')						\
		for (j = 0; j < n; px += (++j)+1)			\
		    *px = *(pval++);					\
	    else							\
		for (j = 0; j < n; px += n-(j++))			\
		    *px = *(pval++);					\
	}								\
    } while (0)

    switch (tx) {
    case REALSXP:
	PM_D_S(double, REAL);
	break;
    case LGLSXP:
	PM_D_S(int, LOGICAL);
	break;
    case INTSXP:
	PM_D_S(int, INTEGER);
	break;
    case CPLXSXP:
	PM_D_S(Rcomplex, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_diag_set");
	break;
    }
    
#undef PM_D_S
    
    SET_SLOT(res, Matrix_xSym, x);
    
    UNPROTECT(4); /* x, res, val, dim */
    return res;
}

/* symmpart(x) */
SEXP packedMatrix_symmpart(SEXP from)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(from, "packedMatrix_symmpart");

    const char *clf = valid[ivalid];
    if (clf[0] == 'd' && clf[1] == 's')
	return from;

    char clt[] = ".spMatrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
    
    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    if (n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */

    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] != 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */
	
    SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    if (ul != 'U')
	SET_SLOT(to, Matrix_uploSym, uplo);
    UNPROTECT(1); /* uplo */
    
    SEXP x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
    REPROTECT(x = (clf[0] == clt[0]) ? duplicate(x) : coerceVector(x, REALSXP),
	      pid);
    if (clf[0] == 'n')
	na2one(x);
    
    if (clf[1] != 's') {

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
	
	int i, j;

#define PM_SYMMPART_TP(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_) \
	do {								\
	    _CTYPE_ *px = _PTR_(x);					\
	    if (ul == 'U') {						\
		for (j = 0; j < n; ++j) {				\
		    for (i = 0; i < j; ++i, ++px)			\
			_ASSIGN_OFFDIAG_;				\
		    ++px;						\
		}							\
		if (di != 'N') {					\
		    px = _PTR_(x);					\
		    for (j = 0; j < n; px += (++j)+1)			\
			_ASSIGN_ONDIAG_;				\
		}							\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    ++px;						\
		    for (i = j+1; i < n; ++i, ++px)			\
			_ASSIGN_OFFDIAG_;				\
		}							\
		if (di != 'N') {					\
		    px = _PTR_(x);					\
		    for (j = 0; j < n; px += n-(j++))			\
			_ASSIGN_ONDIAG_;				\
		}							\
	    }								\
	} while (0)
	
	if (clt[0] != 'z')
	    PM_SYMMPART_TP(double, REAL,
			   *px *= 0.5,
			   *px  = 1.0);
	else
	    PM_SYMMPART_TP(Rcomplex, COMPLEX,
			   do { (*px).r *= 0.5; (*px).i *= 0.5; } while (0),
			   do { (*px).r  = 1.0; (*px).i  = 0.0; } while (0));
	
#undef PM_SYMMPART_TP
	
    } else { /* clf[1] == 's' */
	
	if (clt[0] == 'z')
	    /* Symmetric part of Hermitian matrix is real part */
	    zeroIm(x);

    }

    SET_SLOT(to, Matrix_xSym, x);
    
    UNPROTECT(2); /* x, to */
    return to;
}

/* skewpart(x) */
SEXP packedMatrix_skewpart(SEXP from)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(from, "packedMatrix_skewpart");
    const char *clf = valid[ivalid];

    char clt[] = "...Matrix";
    clt[0] = (clf[0] != 'z') ? 'd' : 'z';
    clt[1] = (clf[1] != 's') ? 'g' : 's';
    clt[2] = (clf[1] != 's') ? 'e' : ((clf[0] != 'z') ? 'C' : 'p');
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));

    SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    if (n > 0)
	SET_SLOT(to, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */
    
    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] != 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */

    SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    if (clf[1] == 's' && ul != 'U')
	SET_SLOT(to, Matrix_uploSym, uplo);
    UNPROTECT(1); /* uplo */
    
    SEXP x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);

    if (clf[1] != 's') {
	
	if ((double) n * n > R_XLEN_T_MAX)
	    error(_("attempt to allocate vector of length exceeding "
		    "R_XLEN_T_MAX"));
	
	SEXP y;
	int i, j;
	R_xlen_t upos = 0, lpos = 0;
	
#define PM_SKEWPART(_CTYPE_, _PTR_, _ASSIGN_OFFDIAG_, _ASSIGN_ONDIAG_)	\
	do {								\
	    _CTYPE_ *px = _PTR_(x), *py = _PTR_(y);			\
	    if (ul == 'U') {						\
		for (j = 0; j < n; ++j) {				\
		    lpos = j;						\
		    for (i = 0; i < j; ++i) {				\
			_ASSIGN_OFFDIAG_(upos, lpos);			\
			++px; ++upos; lpos += n;			\
		    }							\
		    _ASSIGN_ONDIAG_(upos);				\
		    ++px; upos += n-j;					\
		}							\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    upos = lpos;					\
		    _ASSIGN_ONDIAG_(lpos);				\
		    for (i = j+1; i < n; ++i) {				\
			++px; upos += n; ++lpos;			\
			_ASSIGN_OFFDIAG_(lpos, upos);			\
		    }							\
		    ++px; lpos += j+2;					\
		}							\
	    }								\
	} while (0)	    
	
	if (clf[0] != 'z') {
	    PROTECT(y = allocVector(REALSXP, (R_xlen_t) n * n));
	    REPROTECT(x = coerceVector(x, REALSXP), pid);
	    if (clf[0] == 'n')
		na2one(x);
	    
#define ASSIGN_OFFDIAG_DTP(_UPOS_, _LPOS_)			\
	    do {						\
		py[_UPOS_] = 0.5 * *px;				\
		py[_LPOS_] = -py[_UPOS_];			\
	    } while (0)

#define ASSIGN_ONDIAG_DTP(_DPOS_) py[_DPOS_] = 0.0

	    PM_SKEWPART(double, REAL,
			ASSIGN_OFFDIAG_DTP, ASSIGN_ONDIAG_DTP);
	    
#undef ASSIGN_OFFDIAG_DTP
#undef ASSIGN_ONDIAG_DTP
	    
	} else { /* clf[0] == 'z' */

	    PROTECT(y = allocVector(CPLXSXP, (R_xlen_t) n * n));
	    
#define ASSIGN_OFFDIAG_ZTP(_UPOS_, _LPOS_)			\
	    do {						\
		py[_UPOS_].r = 0.5 * (*px).r;			\
		py[_UPOS_].i = 0.5 * (*px).i;			\
		py[_LPOS_].r = -py[upos].r;			\
		py[_LPOS_].i = -py[upos].i;			\
	    } while (0)

#define ASSIGN_ONDIAG_ZTP(_DPOS_) py[_DPOS_].r = py[_DPOS_].i = 0.0
		
	    PM_SKEWPART(Rcomplex, COMPLEX,
			ASSIGN_OFFDIAG_ZTP, ASSIGN_ONDIAG_ZTP);

#undef ASSIGN_OFFDIAG_ZTP
#undef ASSIGN_ONDIAG_ZTP

	}

#undef PM_SKEWPART

	SET_SLOT(to, Matrix_xSym, y);
	UNPROTECT(1); /* y */
	
    } else { /* clf[1] == 's' */
	
	if (clf[0] != 'z') {
	    /* Skew-symmetric part of symmetric matrix is zero matrix */
	    R_xlen_t n1a = (R_xlen_t) n + 1;
	    SEXP p = PROTECT(allocVector(INTSXP, n1a));
	    int *pp = INTEGER(p);
	    Memzero(pp, n1a);
	    SET_SLOT(to, Matrix_pSym, p);
	    UNPROTECT(1); /* p */
	} else {
	    /* Skew-symmetric part of Hermitian matrix is imaginary part */
	    REPROTECT(x = duplicate(x), pid);
	    zeroRe(x);
	    SET_SLOT(to, Matrix_xSym, x);
	}
	
    }
    
    UNPROTECT(2); /* x, to */
    return to;
}

#define PM_XIJ_TR_UP_NUN(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ <= _J_							\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : _ZERO_)

#define PM_XIJ_TR_UP_UNT(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ < _J_							\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : (_I_ == _J_ ? _ONE_ : _ZERO_))

#define PM_XIJ_TR_LO_NUN(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : _ZERO_)

#define PM_XIJ_TR_LO_UNT(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ > _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : (_I_ == _J_ ? _ONE_ : _ZERO_))

#define PM_XIJ_SY_UP(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ <= _J_				       			\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))	       			\
     : *(_X_ + PM_AR21_UP(_J_, _I_)))

#define PM_XIJ_SY_LO(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : *(_X_ + PM_AR21_LO(_J_, _I_, _N2_)))

#define PM_SUB1_LOOP(_XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	int i, j;							\
	R_xlen_t k;							\
	if (TYPEOF(index) == INTSXP) {					\
	    int pos, *pindex = INTEGER(index);				\
	    for (k = 0; k < nindex; ++k) {				\
		pos = *(pindex++);					\
		if (pos == NA_INTEGER || pos > nn)			\
		    *(pres++) = _NA_;					\
		else {							\
		    pos -= 1; /* 1-index -> 0-index */			\
		    i = pos % n;					\
		    j = pos / n;					\
		    *(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);	\
		}							\
	    }								\
	} else {							\
	    double pos, *pindex = REAL(index);				\
	    R_xlen_t truncpos;						\
	    for (k = 0; k < nindex; ++k) {				\
		pos = *(pindex++);					\
		if (!R_FINITE(pos) || (truncpos = (R_xlen_t) pos) > nn) \
		    *(pres++) = _NA_;					\
		else {							\
		    truncpos -= 1; /* 1-index -> 0-index */		\
		    i = truncpos % n;					\
		    j = truncpos / n;					\
		    *(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);	\
		}							\
	    }								\
	}								\
    } while (0)

#define PM_SUB1_END(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)		\
    do {								\
	_CTYPE_ *pres = _PTR_(res), *px = _PTR_(x);			\
	if (di == '\0') { /* symmetric */				\
	    if (ul == 'U')						\
		PM_SUB1_LOOP(PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
	    else							\
		PM_SUB1_LOOP(PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
	} else if (di == 'N') { /* non-unit triangular */		\
	    if (ul == 'U')						\
		PM_SUB1_LOOP(PM_XIJ_TR_UP_NUN, _NA_, _ZERO_, _ONE_);	\
	    else							\
		PM_SUB1_LOOP(PM_XIJ_TR_LO_NUN, _NA_, _ZERO_, _ONE_);	\
	} else { /* unit triangular */					\
	    if (ul == 'U')						\
		PM_SUB1_LOOP(PM_XIJ_TR_UP_UNT, _NA_, _ZERO_, _ONE_);	\
	    else							\
		PM_SUB1_LOOP(PM_XIJ_TR_LO_UNT, _NA_, _ZERO_, _ONE_);	\
	}    								\
    } while (0)

#define PM_SUB1								\
    do {								\
	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));		\
	char ul = *CHAR(STRING_ELT(uplo, 0));				\
	UNPROTECT(1); /* uplo */					\
									\
	char di = '\0';							\
	if (HAS_SLOT(obj, Matrix_diagSym)) {				\
	    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));		\
	    di = *CHAR(STRING_ELT(diag, 0));				\
	    UNPROTECT(1); /* diag */					\
	}								\
									\
	SEXPTYPE tx;							\
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),			\
	    res = PROTECT(allocVector(tx = TYPEOF(x), nindex));		\
									\
	switch (tx) {							\
	case REALSXP: /* d..Matrix */					\
	    PM_SUB1_END(double, REAL, NA_REAL, 0.0, 1.0);		\
	    break;							\
	case LGLSXP: /* [ln]..Matrix */					\
	    PM_SUB1_END(int, LOGICAL, NA_LOGICAL, 0, 1);		\
	    break;							\
	case INTSXP: /* i..Matrix */					\
	    PM_SUB1_END(int, INTEGER, NA_INTEGER, 0, 1);		\
	    break;							\
	case CPLXSXP: /* z..Matrix */					\
	{								\
	    Rcomplex na, zero, one;					\
	    na.r = NA_REAL; zero.r = 0.0; one.r = 1.0;			\
	    na.i = NA_REAL; zero.i = 0.0; one.i = 0.0;			\
	    PM_SUB1_END(Rcomplex, COMPLEX, na, zero, one);		\
	    break;							\
	}								\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_sub1");	\
	    break;							\
	}								\
									\
	UNPROTECT(2); /* res, x */					\
	return res;							\
    } while (0)
	  
/* 'x[i]' where 'i' is an integer or double vector with elements 
   greater than or equal to 1
*/
SEXP packedMatrix_sub1(SEXP obj, SEXP index)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("indexing n-by-n packedMatrix is not supported "
		"for n*n exceeding R_XLEN_T_MAX"));

    R_xlen_t n2 = (R_xlen_t) n * 2, nn = (R_xlen_t) n * n,
	nindex = XLENGTH(index);
    PM_SUB1;
}

#undef PM_SUB1_LOOP

#define PM_SUB1_LOOP(_XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	int i, j, k, *pi = INTEGER(index), *pj = pi + nindex;		\
	for (k = 0; k < nindex; ++k) {					\
	    i = *(pi++);						\
	    j = *(pj++);						\
	    if (i == NA_INTEGER || j == NA_INTEGER)			\
		*(pres++) = _NA_;					\
	    else {							\
		i -= 1; /* 1-index -> 0-index */			\
		j -= 1;							\
		*(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);		\
	    }								\
	}								\
    } while (0)

/* 'x[i]' where 'i' is a 2-column integer matrix supplying 
   integers in 'c(1:n, NA)' _only_
*/
SEXP packedMatrix_sub1_mat(SEXP obj, SEXP index)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    
    PROTECT(dim = getAttrib(index, R_DimSymbol));
    int nindex = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("indexing n-by-n packedMatrix is not supported "
		"for n*n exceeding R_XLEN_T_MAX"));

    R_xlen_t n2 = (R_xlen_t) n * 2;
    PM_SUB1;
}

#undef PM_SUB1
#undef PM_SUB1_LOOP

#define PM_SUB2_LOOP(_FOR_, _XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
        int i, j, ki, kj;						\
	for (kj = 0; kj < nj; ++kj) {					\
	    if (mj)							\
		j = kj;							\
	    else {							\
		j = pj[kj];						\
		if (j != NA_INTEGER)					\
		    j -= 1;						\
		else {							\
		    _FOR_ {						\
			*(px1++) = _NA_;				\
		    }							\
		    if (do_cn)						\
			SET_STRING_ELT(cn1, kj, NA_STRING);		\
		    continue;						\
		}							\
	    }								\
	    _FOR_ {							\
		if (mi)							\
		    i = ki;						\
		else {							\
		    i = pi[ki];						\
		    if (i != NA_INTEGER)				\
			i -= 1;						\
		    else {						\
			*(px1++) = _NA_;				\
			continue;					\
		    }							\
		}							\
		*(px1++) = _XIJ_(px0, i, j, n2, _ZERO_, _ONE_);		\
	    }								\
	    if (do_cn)							\
		SET_STRING_ELT(cn1, kj, STRING_ELT(cn0, j));		\
	}								\
	if (do_rn) {							\
	    for (ki = 0; ki < ni; ++ki) {				\
		if (mi)							\
		    i = ki;						\
		else {							\
		    i = pi[ki];						\
		    if (i != NA_INTEGER)				\
			i -= 1;						\
		    else {						\
			SET_STRING_ELT(rn1, ki, NA_STRING);		\
			continue;					\
		    }							\
		}							\
		SET_STRING_ELT(rn1, ki, STRING_ELT(rn0, i));		\
	    }								\
	}								\
    } while (0)

#define PM_SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	if (di == '\0') { /* symmetric */				\
	    if (ul == 'U') {						\
		if (do_sp)						\
		    PM_SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),		\
				 PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
		else							\
		    PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				 PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		if (do_sp)						\
		    PM_SUB2_LOOP(for (ki = kj; ki < ni; ++ki),		\
				 PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
		else							\
		    PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				 PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else if (di == 'N') { /* non-unit triangular */		\
	    if (ul == 'U')						\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_UP_NUN, _NA_, _ZERO_, _ONE_);	\
	    else							\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_LO_NUN, _NA_, _ZERO_, _ONE_);	\
	} else { /* unit triangular */					\
	    if (ul == 'U')						\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_UP_UNT, _NA_, _ZERO_, _ONE_);	\
	    else							\
		PM_SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_LO_UNT, _NA_, _ZERO_, _ONE_);	\
	}								\
    } while (0)

/* 'x[i, ]', 'x[, j]', and 'x[i, j]' where 'i' and 'j' are integer vectors 
   supplying integers in 'c(1:n, NA)' _only_ ... NULL indicates missingness
*/
SEXP packedMatrix_sub2(SEXP obj, SEXP index1, SEXP index2, SEXP drop)
{
    static const char *valid[] = {
	/* 0 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/* 3 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(obj, "packedMatrix_sub2");
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    
    Rboolean mi = isNull(index1), mj = isNull(index2);
    int *pi = NULL, *pj = NULL;
    R_xlen_t ni, nj, n2 = (R_xlen_t) n * 2;
    
    if ((double) n * n > R_XLEN_T_MAX)
	error(_("indexing of n-by-n packedMatrix with n*n "
		"exceeding R_XLEN_T_MAX is not supported"));
    if (mi)
	ni = n;
    else {
	ni = XLENGTH(index1);
	if (ni > INT_MAX)
	    error(_("dimensions cannot exceed 2^31-1"));
	pi = INTEGER(index1);
    }
    if (mj)
	nj = n;
    else {
	nj = XLENGTH(index2);
	if (nj > INT_MAX)
	    error(_("dimensions cannot exceed 2^31-1"));
	pj = INTEGER(index2);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    
    char di = '\0';
    if (ivalid < 3) {
	SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
	di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
    }

    /* Initialize result of same "kind" but general "shape",
       except for symmetric indexing of symmetric matrix, 
       where shape is also retained */
    SEXP res;
    PROTECT_INDEX pid;
    Rboolean do_sp = ivalid >= 3 && !mi && !mj && ni == nj &&
	memcmp(pi, pj, ni * sizeof(int)) == 0;
    if (do_sp) {
	PROTECT_WITH_INDEX(res = NEW_OBJECT_OF_CLASS(valid[ivalid]), &pid);
	if (ul != 'U')
	    SET_SLOT(res, Matrix_uploSym, uplo);
    } else {
	char cl[] = ".geMatrix";
	cl[0] = valid[ivalid][0];
	PROTECT_WITH_INDEX(res = NEW_OBJECT_OF_CLASS(cl), &pid);
    }
    
    /* Set 'Dim' slot */
    PROTECT(dim = GET_SLOT(res, Matrix_DimSym));
    int *pdim = INTEGER(dim);
    pdim[0] = (int) ni;
    pdim[1] = (int) nj;
    UNPROTECT(1); /* dim */
    
    /* Set 'Dimnames' slot and 'names(Dimnames)' */
    SEXP
	dn0 = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	dn1 = PROTECT(GET_SLOT(res, Matrix_DimNamesSym)),
	ndn0 = PROTECT(getAttrib(dn0, R_NamesSymbol)),
	rn0, rn1, cn0, cn1;
    rn0 = rn1 = cn0 = cn1 = NULL;
    if (di == '\0') { /* symmetric */
    	int J;
	if (isNull(rn0 = cn0 = VECTOR_ELT(dn0, J = 1)))
	    rn0 = cn0 = VECTOR_ELT(dn0, J = 0);
	if (!isNull(ndn0)) {
	    SEXP ndn1 = PROTECT(allocVector(STRSXP, 2)),
		s = STRING_ELT(ndn0, J);
	    SET_STRING_ELT(ndn1, 0, s);
	    SET_STRING_ELT(ndn1, 1, s);
	    setAttrib(dn1, R_NamesSymbol, ndn1);
	    UNPROTECT(1); /* ndn1 */
	}
    } else { /* triangular */
	rn0 = VECTOR_ELT(dn0, 0);
	cn0 = VECTOR_ELT(dn0, 1);
	if (!isNull(ndn0))
	    setAttrib(dn1, R_NamesSymbol, ndn0);
    }
    UNPROTECT(1); /* ndn0 */
    
    Rboolean
	has_rn = !isNull(rn0) && ni > 0, do_rn = FALSE,
	has_cn = !isNull(cn0) && nj > 0, do_cn = FALSE;
    if (has_rn) {
	if (mi)
	    SET_VECTOR_ELT(dn1, 0, rn0);
	else {
	    PROTECT(rn1 = allocVector(STRSXP, ni));
	    do_rn = TRUE;
	}
    }
    if (has_cn && !do_sp) {
	if (mj)
	    SET_VECTOR_ELT(dn1, 1, cn0);
	else {
	    PROTECT(cn1 = allocVector(STRSXP, nj));
	    do_cn = TRUE;
	}
    }

    /* Set 'x' slot */
    SEXPTYPE tx;
    R_xlen_t nx = (do_sp ? ni + (ni * (ni - 1)) / 2 : ni * nj);
    SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	x1 = PROTECT(allocVector(tx = TYPEOF(x0), nx));
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	PM_SUB2(double, REAL, NA_REAL, 0.0, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_SUB2(int, LOGICAL, NA_LOGICAL, 0, 1);
	break;
    case INTSXP: /* i..Matrix */
	PM_SUB2(int, INTEGER, NA_INTEGER, 0, 1);
	break;
    case CPLXSXP: /* z..Matrix */
    {
	Rcomplex na, zero, one;
	na.r = NA_REAL; zero.r = 0.0; one.r = 1.0;
	na.i = NA_REAL; zero.i = 0.0; one.i = 0.0;
	PM_SUB2(Rcomplex, COMPLEX, na, zero, one);
	break;
    }
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "packedMatrix_sub2");
	break;
    }

    SET_SLOT(res, Matrix_xSym, x1);
    UNPROTECT(2); /* x1, x0 */
    if (do_cn) {
	SET_VECTOR_ELT(dn1, 1, cn1);
	UNPROTECT(1); /* cn1 */
    }
    if (do_rn) {
	SET_VECTOR_ELT(dn1, 0, rn1);
	UNPROTECT(1); /* rn1 */
    }

    /* Drop dimensions in this special case */
    if (asLogical(drop) != 0 && (ni == 1 || nj == 1)) {
	REPROTECT(res = GET_SLOT(res, Matrix_xSym), pid);
	if (has_rn && nj == 1 && ni != 1)
	    setAttrib(res, R_NamesSymbol, VECTOR_ELT(dn1, 0));
	else if (has_cn && ni == 1 && nj != 1)
	    setAttrib(res, R_NamesSymbol, VECTOR_ELT(dn1, 1));
    }
    
    UNPROTECT(4); /* dn0, dn1, res, uplo */
    return res;
}

#undef PM_SUB2
#undef PM_SUB2_LOOP
