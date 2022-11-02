#include <limits.h>
#include <R_ext/Lapack.h>
#include "Mutils.h"

/* More for 'Dimnames' ============================================== */

Rboolean DimNames_is_trivial(SEXP dn)
{
    if (!(isNull(VECTOR_ELT(dn, 0)) &&
	  isNull(VECTOR_ELT(dn, 1))))
	return FALSE;
    Rboolean res = TRUE;
    SEXP ndn = PROTECT(getAttrib(dn, R_NamesSymbol));
    if (!isNull(ndn))
	res = FALSE;
    UNPROTECT(1);
    return res;
}

Rboolean DimNames_is_symmetric(SEXP dn)
{
    /* NB: Assuming here that we have the 'Dimnames' slot 
       of a _valid_ Matrix, so that the elements are either 
       NULL or character vectors

       Keep synchronized with symmetricMatrix_validate() above,
       (which must do slightly more)! */

    SEXP rn, cn;
    int n;
    if (!isNull(rn = VECTOR_ELT(dn, 0)) &&
	!isNull(cn = VECTOR_ELT(dn, 1)) &&
	rn != cn &&
	((n = LENGTH(rn)) != LENGTH(cn) || !equal_string_vectors(rn, cn, n)))
	return FALSE;
    Rboolean res = TRUE;
    SEXP ndn = PROTECT(getAttrib(dn, R_NamesSymbol));
    const char *ndn0, *ndn1;
    if (!isNull(ndn) &&
	*(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
	*(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
	strcmp(ndn0, ndn1) != 0)
	res = FALSE;
    UNPROTECT(1);
    return res;
}

SEXP R_DimNames_is_symmetric(SEXP dn)
{
    return ScalarLogical(DimNames_is_symmetric(dn));
}

/**
 * @brief Produce symmetric `Dimnames` from possibly asymmetric ones.
 * 
 * Roughly `dest[1:2] <- rep(src[j], 2)`, where `j` is either 1 or 2
 * depending on `J`.  If `J` is 0 or 1, then `j = J+1`.  If `J` is -1, 
 * then `j = 1` if and only if `src[[2]]` is `NULL` and `src[[1]]`
 * is not.  For speed, it is assumed that `dest` is newly allocated,
 * i.e., that it is `list(NULL, NULL)`.
 *
 * @param dest,src Lists of length 2, typically the `Dimnames` slots
 *     of two square `Matrix` of equal size.
 * @param J An integer, one of -1, 0, and 1.
 */
void symmDN(SEXP dest, SEXP src, int J /* -1|0|1 */)
{
    SEXP s;
    if (J < 0) {
	if (!isNull(s = VECTOR_ELT(src, J = 1)) ||
	    !isNull(s = VECTOR_ELT(src, J = 0))) {
	    SET_VECTOR_ELT(dest, 0, s);
	    SET_VECTOR_ELT(dest, 1, s);
	} else {
	    J = 1;
	}
    } else {
	if (!isNull(s = VECTOR_ELT(src, J))) {
	    SET_VECTOR_ELT(dest, 0, s);
	    SET_VECTOR_ELT(dest, 1, s);
	}
    }
    /* names(dimnames(.)) */
    PROTECT(s = getAttrib(src, R_NamesSymbol));
    if (!isNull(s)) {
	SEXP destnms = PROTECT(allocVector(STRSXP, 2));
	if (*CHAR(s = STRING_ELT(s, J)) != '\0') {
	    SET_STRING_ELT(destnms, 0, s);
	    SET_STRING_ELT(destnms, 1, s);
	}
	setAttrib(dest, R_NamesSymbol, destnms);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return;
}

/**
 * @brief Reverse (or "transpose") `Dimnames`.
 * 
 * Roughly `dest[1:2] <- src[2:1]`.  For speed, it is assumed that 
 * `dest` is newly allocated, i.e., that it is `list(NULL, NULL)`.
 *
 * @param dest,src Lists of length 2, typically the `Dimnames` slots
 *     of two square `Matrix` of equal size.
 */
void revDN(SEXP dest, SEXP src) {
    SEXP s;
    if (!isNull(s = VECTOR_ELT(src, 0)))
	SET_VECTOR_ELT(dest, 1, s);
    if (!isNull(s = VECTOR_ELT(src, 1)))
	SET_VECTOR_ELT(dest, 0, s);
    PROTECT(s = getAttrib(src, R_NamesSymbol));
    if (!isNull(s)) {
	SEXP srcnms = s, destnms = PROTECT(allocVector(STRSXP, 2));
	if (*CHAR(s = STRING_ELT(srcnms, 0)) != '\0')
	    SET_STRING_ELT(destnms, 1, s);
	if (*CHAR(s = STRING_ELT(srcnms, 1)) != '\0')
	    SET_STRING_ELT(destnms, 0, s);
	setAttrib(dest, R_NamesSymbol, destnms);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return;
}

SEXP R_symmDN(SEXP dn)
{
    /* Be fast (do nothing!) when dimnames = list(NULL, NULL) */
    if (DimNames_is_trivial(dn))
	return dn;
    SEXP newdn = PROTECT(allocVector(VECSXP, 2));
    symmDN(newdn, dn, -1);
    UNPROTECT(1);
    return newdn;
}

SEXP R_revDN(SEXP dn)
{
    /* Be fast (do nothing!) when dimnames = list(NULL, NULL) */
    if (DimNames_is_trivial(dn))
	return dn;
    SEXP newdn = PROTECT(allocVector(VECSXP, 2));
    revDN(newdn, dn);
    UNPROTECT(1);
    return newdn;
}

SEXP get_symmetrized_DimNames(SEXP obj, int J) {
    SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
    if (DimNames_is_trivial(dn)) {
	UNPROTECT(1);
	return dn;
    }
    SEXP newdn = PROTECT(allocVector(VECSXP, 2));
    symmDN(newdn, dn, J);
    UNPROTECT(2);
    return newdn;
}

SEXP get_reversed_DimNames(SEXP obj) {
    SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
    if (DimNames_is_trivial(dn)) {
	UNPROTECT(1);
	return dn;
    }
    SEXP newdn = PROTECT(allocVector(VECSXP, 2));
    revDN(newdn, dn);
    UNPROTECT(2);
    return newdn;
}

void set_symmetrized_DimNames(SEXP obj, SEXP dn, int J) {
    if (!DimNames_is_trivial(dn)) {
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symmDN(newdn, dn, J);
	SET_SLOT(obj, Matrix_DimNamesSym, newdn);
	UNPROTECT(1);
    }
    return;
}

void set_reversed_DimNames(SEXP obj, SEXP dn) {
    if (!DimNames_is_trivial(dn)) {
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	revDN(newdn, dn);
	SET_SLOT(obj, Matrix_DimNamesSym, newdn);
	UNPROTECT(1);
    }
    return;
}

void set_DimNames(SEXP obj, SEXP dn)
{
    if (!DimNames_is_trivial(dn)) {
	SEXP s, newdn = PROTECT(allocVector(VECSXP, 2));
	if (!isNull(s = VECTOR_ELT(dn, 0)))
	    SET_VECTOR_ELT(newdn, 0, s);
	if (!isNull(s = VECTOR_ELT(dn, 1)))
	    SET_VECTOR_ELT(newdn, 1, s);
	PROTECT(s = getAttrib(dn, R_NamesSymbol));
	if (!isNull(s))
	    setAttrib(newdn, R_NamesSymbol, s);
	SET_SLOT(obj, Matrix_DimNamesSym, newdn);
	UNPROTECT(2);
    }
    return;
}


/* More for 'factors' =============================================== */

SEXP get_factor(SEXP obj, const char *nm)
{
    SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym));
    if (LENGTH(factors) > 0) {
	SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
	R_xlen_t i = strmatch(nm, valid);
	UNPROTECT(1);
	if (i >= 0) {
	    factors = VECTOR_ELT(factors, i);
	    UNPROTECT(1);
	    return factors;
	}
    }
    UNPROTECT(1);
    return R_NilValue;
}

void set_factor(SEXP obj, const char *nm, SEXP val)
{
    PROTECT(val);
    SEXP factors;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(factors = GET_SLOT(obj, Matrix_factorSym), &pid);
    if (LENGTH(factors) > 0) {
	SEXP valid = PROTECT(getAttrib(factors, R_NamesSymbol));
	R_xlen_t i = strmatch(nm, valid);
	UNPROTECT(1);
	if (i >= 0) {
	    SET_VECTOR_ELT(factors, i, val);
	    UNPROTECT(2);
	    return;
	}
    }
    REPROTECT(factors = append_to_named_list(factors, nm, val), pid);
    SET_SLOT(obj, Matrix_factorSym, factors);
    UNPROTECT(2);
    return;
}

/** 
 * @brief Subassign by name to the `factors` slot of a `compMatrix`.
 * 
 * Like `obj\@factors[[nm]] <- val`, but modifying `obj` (rather than a copy)
 * _even if_ `obj` is referenced elsewhere, supporting "automagic" caching of
 * factorizations by R functions taking `compMatrix` as an argument.
 * _Use with care!_
 * 
 * @param obj A `compMatrix`.
 * @param nm A length-1 `STRSXP` giving a factor name.
 * @param val A `SEXP`, usually a `MatrixFactorization`.
 * @param warn A length-1 `LGLSXP`. Warn if `obj` has no `factors` slot 
 *     (in which case `obj` is untouched)?
 *
 * @return `val`.
 */
SEXP R_set_factor(SEXP obj, SEXP nm, SEXP val, SEXP warn)
{
    if (TYPEOF(nm) != STRSXP || LENGTH(nm) < 1 ||
	(nm = STRING_ELT(nm, 0)) == NA_STRING)
	error(_("invalid factor name"));
    else if (HAS_SLOT(obj, Matrix_factorSym))
	set_factor(obj, CHAR(nm), val);
    else if (asLogical(warn) != 0)
	warning(_("attempt to set factor on Matrix "
		  "without 'factors' slot"));
    return val;
}

/** 
 * @brief Empty the 'factors' slot of a 'compMatrix'.
 * 
 * Like `obj\@factors <- list()`, but modifying `obj` (rather than a copy)
 * _even if_ `obj` is referenced elsewhere, supporting "automagic" clearing
 * of the `factors` slot by R functions taking `compMatrix` as an argument. 
 * _Use with care!_
 * 
 * @param obj A `compMatrix`.
 * @param warn A length-1 LGLSXP. Warn if `obj` has no `factors` slot 
 *     (in which case `obj` is untouched)?
 *
 * @return `TRUE` if `obj` has a nonempty `factors` slot, `FALSE` otherwise.
 */
SEXP R_empty_factors(SEXP obj, SEXP warn)
{
    /* If there is a nonempty 'factors' slot, then replace it with list() */
    if (HAS_SLOT(obj, Matrix_factorSym)) {
	SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym));
	if (LENGTH(factors) > 0) {
	    PROTECT(factors = allocVector(VECSXP, 0));
	    SET_SLOT(obj, Matrix_factorSym, factors);
	    UNPROTECT(2);
	    return ScalarLogical(1); /* slot was reset */
	}
	UNPROTECT(1);
    } else if (asLogical(warn) != 0)
	warning(_("attempt to discard factors from Matrix "
		  "without 'factors' slot"));
    return ScalarLogical(0); /* no-op */
}


/* For coercions ==================================================== */

#define PACK(_PREFIX_, _CTYPE_, _ONE_)					\
void _PREFIX_ ## dense_pack(_CTYPE_ *dest, const _CTYPE_ *src, int n,   \
			    char uplo, char diag)			\
{									\
    int i, j;								\
    R_xlen_t dpos = 0, spos = 0;					\
    if (uplo == 'U') {							\
	for (j = 0; j < n; spos += n-(++j))				\
	    for (i = 0; i <= j; ++i)					\
		dest[dpos++] = src[spos++];				\
	if (diag != 'N') {						\
	    dpos = 0;							\
	    for (j = 0; j < n; dpos += (++j)+1)				\
		dest[dpos] = _ONE_;					\
	}								\
    } else {								\
    	for (j = 0; j < n; spos += (++j))				\
	    for (i = j; i < n; ++i)					\
		dest[dpos++] = src[spos++];				\
	if (diag != 'N') {						\
	    dpos = 0;							\
	    for (j = 0; j < n; dpos += n-(j++))				\
		dest[dpos] = _ONE_;					\
	}								\
    }									\
    return;								\
}

/**
 * @brief Pack a square `unpackedMatrix`.
 * 
 * Copies the upper or lower triangular part of `src` to `dest`,
 * where it is stored contiguously ("packed"). Optionally resets 
 * the diagonal elements to 1.
 * 
 * @param dest,src Pointers to the first elements of length-`(n*(n+1))/2`
 *     and length-`n*n` (resp.) arrays, usually the "data" of the `x`
 *     slot of an `n`-by-`n` `packedMatrix` and `unpackedMatrix` (resp.).
 * @param n Size of matrix being packed.
 * @param uplo,diag `char` specifying whether the "nontrivial part" 
 *     is upper (`'U'`) or lower (`'L'`) and whether to "force" a
 *     unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_pack() */
PACK(d, double, 1.0)
/* idense_pack() */
PACK(i, int, 1)
/* zdense_pack() */
PACK(z, Rcomplex, Matrix_zone)

#undef PACK

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

#define UNPACKED_MAKE_SYMMETRIC_SET_EQ(_X_, _DEST_, _SRC_)	\
    _X_[_DEST_] = _X_[_SRC_]

#define UNPACKED_MAKE_SYMMETRIC_SET_CJ(_X_, _DEST_, _SRC_)		\
    do {								\
	_X_[_DEST_].r =  _X_[_SRC_].r;					\
	_X_[_DEST_].i = -_X_[_SRC_].i;					\
    } while (0)

#define UNPACKED_MAKE_SYMMETRIC(_PREFIX_, _CTYPE_, _SET_)		\
void _PREFIX_ ## dense_unpacked_make_symmetric(_CTYPE_ *x, int n, char uplo) \
{									\
    int i, j, n1s = n - 1;						\
    R_xlen_t upos = n, lpos = 1;					\
    if (uplo == 'U') {							\
	for (j = 0; j < n; upos = (lpos += (++j)+1) + n1s)		\
	    for (i = j+1; i < n; ++i, upos += n, ++lpos)		\
		_SET_(x, lpos, upos);					\
    } else {								\
	for (j = 0; j < n; upos = (lpos += (++j)+1) + n1s)		\
	    for (i = j+1; i < n; ++i, upos += n, ++lpos)		\
		_SET_(x, upos, lpos);					\
    }									\
    return;								\
}

/**
 * @brief Make symmetric a square `unpackedMatrix`.
 * 
 * "Symmetrizes" the elements of an `n`-by-`n` `unpackedMatrix`.
 * 
 * @param x A pointer to the first element of a length-`n*n` array,
 *     usually the "data" of the `x` slot of an `n`-by-`n` `unpackedMatrix`.
 * @param n Size of matrix being symmetrized.
 * @param uplo A `char` specifying whether to copy the upper triangle
 *     to the lower triangle (`'U'`) or to do the reverse (`'L'`).
 */
/* ddense_unpacked_make_symmetric() */
UNPACKED_MAKE_SYMMETRIC(d, double, UNPACKED_MAKE_SYMMETRIC_SET_EQ)
/* idense_unpacked_make_symmetric() */
UNPACKED_MAKE_SYMMETRIC(i, int, UNPACKED_MAKE_SYMMETRIC_SET_EQ)
/* zdense_unpacked_make_symmetric() */
UNPACKED_MAKE_SYMMETRIC(z, Rcomplex, UNPACKED_MAKE_SYMMETRIC_SET_CJ)

#undef UNPACKED_MAKE_SYMMETRIC
#undef UNPACKED_MAKE_SYMMETRIC_SET_CJ
#undef UNPACKED_MAKE_SYMMETRIC_SET_EQ

#define UNPACKED_MAKE_TRIANGULAR(_PREFIX_, _CTYPE_, _ZERO_, _ONE_)	\
void _PREFIX_ ## dense_unpacked_make_triangular(_CTYPE_ *x, int m, int n, \
						char uplo, char diag)	\
{									\
    int i, j, r = (m < n) ? m : n;					\
    R_xlen_t pos = 0;							\
    if (uplo == 'U') {							\
	for (j = 0; j < r; pos += (++j)+1)				\
	    for (i = j+1; i < m; ++i)					\
		x[++pos] = _ZERO_;					\
    } else {								\
	for (j = 0; j < r; pos += m-(j++))				\
	    for (i = 0; i < j; ++i)					\
		x[pos++] = _ZERO_;					\
	for (j = r; j < n; ++j)						\
	    for (i = 0; i < m; ++i)					\
		x[pos++] = _ZERO_;					\
    }									\
    if (diag != 'N') {							\
	pos = 0;							\
	R_xlen_t m1a = (R_xlen_t) m + 1;				\
	for (j = 0; j < r; ++j, pos += m1a)				\
	    x[pos] = _ONE_;						\
    }									\
    return;								\
}

/**
 * @brief Make triangular an `unpackedMatrix`.
 * 
 * "Triangularizes" the elements of an `m`-by-`n` `unpackedMatrix`,
 * which need not be square (though all `triangularMatrix` _are_).
 * 
 * @param x A pointer to the first element of a length-`m*n` array,
 *     usually the "data" of the `x` slot of an `m`-by-`n` `unpackedMatrix`.
 * @param m,n Dimensions of matrix being triangularized.
 * @param uplo,diag `char` constants specifying whether the matrix
 *     should be upper (`'U'`) or lower (`'L'`) triangular and
 *     whether it should have a unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_unpacked_make_triangular() */
UNPACKED_MAKE_TRIANGULAR(d, double, 0.0, 1.0)
/* idense_unpacked_make_triangular() */
UNPACKED_MAKE_TRIANGULAR(i, int, 0, 1)
/* zdense_unpacked_make_triangular() */
UNPACKED_MAKE_TRIANGULAR(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef UNPACKED_MAKE_TRIANGULAR

#define UNPACKED_MAKE_BANDED(_PREFIX_, _CTYPE_, _ZERO_, _ONE_)		\
void _PREFIX_ ## dense_unpacked_make_banded(_CTYPE_ *x,			\
					    int m, int n, int a, int b,	\
					    char diag)			\
{									\
    if (m == 0 || n == 0)						\
	return;								\
    if (a > b || a >= n || b <= -m) {					\
	Memzero(x, (R_xlen_t) m * n);					\
	return;								\
    }									\
    if (a <= -m) a = 1-m;						\
    if (b >=  n) b = n-1;						\
    									\
    int i, j, i0, i1,							\
	j0 = (a < 0) ? 0 : a,						\
	j1 = (b < n-m) ? m+b : n;					\
    									\
    if (j0 > 0) {							\
	R_xlen_t dx;							\
	Memzero(x, dx = (R_xlen_t) m * j0);				\
	x += dx;							\
    }									\
    for (j = j0; j < j1; ++j, x += m) {					\
	i0 = j - b;							\
	i1 = j - a + 1;							\
	for (i = 0; i < i0; ++i)					\
	    *(x + i) = _ZERO_;						\
	for (i = i1; i < m; ++i)					\
	    *(x + i) = _ZERO_;						\
    }									\
    if (j1 < n)								\
	Memzero(x, (R_xlen_t) m * (n - j1));				\
    if (diag != 'N' && a <= 0 && b >= 0) {				\
	x -= m * (R_xlen_t) j;						\
	R_xlen_t m1a = (R_xlen_t) m + 1;				\
	for (j = 0; j < n; ++j, x += m1a)				\
	    *x = _ONE_;							\
    }									\
    return;								\
}
UNPACKED_MAKE_BANDED(d, double, 0.0, 1.0)
UNPACKED_MAKE_BANDED(i, int, 0, 1)
UNPACKED_MAKE_BANDED(z, Rcomplex, Matrix_zzero, Matrix_zone)
#undef UNPACKED_MAKE_BANDED

#define PACKED_MAKE_BANDED(_PREFIX_, _CTYPE_, _ZERO_, _ONE_)		\
void _PREFIX_ ## dense_packed_make_banded(_CTYPE_ *x,		        \
					  int n, int a, int b,		\
					  char uplo, char diag)		\
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
PACKED_MAKE_BANDED(d, double, 0.0, 1.0)
PACKED_MAKE_BANDED(i, int, 0, 1)
PACKED_MAKE_BANDED(z, Rcomplex, Matrix_zzero, Matrix_zone)
#undef PACKED_MAKE_BANDED

#define UNPACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_, _ONE_)		\
void _PREFIX_ ## dense_unpacked_copy_diagonal(_CTYPE_ *dest,	        \
					      const _CTYPE_ *src,	\
					      int n, R_xlen_t len,	\
    					      char uplo, char diag)	\
{									\
    int j;								\
    R_xlen_t n1a = (R_xlen_t) n + 1;					\
    if (diag != 'N') {							\
	for (j = 0; j < n; ++j, dest += n1a)				\
	    *dest = _ONE_;						\
    } else if (len == n) {						\
	/* copying from diagonalMatrix */				\
	for (j = 0; j < n; ++j, dest += n1a, ++src)			\
	    *dest = *src;						\
    } else if (len == (n * n1a) / 2) {					\
	/* copying from packedMatrix */					\
	if (uplo == 'U') {						\
	    for (j = 0; j < n; dest += n1a, src += (++j)+1)		\
		*dest = *src;						\
	} else {							\
	    for (j = 0; j < n; dest += n1a, src += n-(j++))		\
		*dest = *src;						\
	}								\
    } else if (len == (R_xlen_t) n * n) {				\
	/* copying from square unpackedMatrix */			\
	for (j = 0; j < n; ++j, dest += n1a, src += n1a)		\
	    *dest = *src;						\
    } else {								\
	error(_("incompatible 'n' and 'len' to '*_copy_diagonal()'"));	\
    }									\
    return;								\
}

/**
 * Copy a length-`n` diagonal to a length-`n*n` array.
 * 
 * @param dest A pointer to the first element of a length-`n*n` array,
 *     usually the "data" of the `x` slot of an `n`-by-`n` `unpackedMatrix`.
 * @param src A pointer to the first element of a length-`n`, 
 *     length-`(n*(n+1))/2`, or length-`n*n` array, usually the "data"
 *     of the `x` slot of an `n`-by-`n` `diagonalMatrix`, `packedMatrix`,
 *     or `unpackedMatrix`, respectively.
 * @param n Size of matrix being copied from and to.
 * @param len Length of `src` array.
 * @param uplo,diag `char` constants specifying whether `src` stores an 
 *     upper (`'U'`) or lower (`'L'`) triangle when `len == (n*(n+1))/2` and
 *     whether the matrix should have a unit diagonal (`'U'`) or not (`'N'`).
 */
/* ddense_unpacked_copy_diagonal() */
UNPACKED_COPY_DIAGONAL(d, double, 1.0)
/* idense_unpacked_copy_diagonal() */
UNPACKED_COPY_DIAGONAL(i, int, 1)
/* zdense_unpacked_copy_diagonal() */
UNPACKED_COPY_DIAGONAL(z, Rcomplex, Matrix_zone)

#undef UNPACKED_COPY_DIAGONAL

#define PACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_, _ONE_)			\
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
PACKED_COPY_DIAGONAL(d, double, 1.0)
/* idense_packed_copy_diagonal() */
PACKED_COPY_DIAGONAL(i, int, 1)
/* zdense_packed_copy_diagonal() */
PACKED_COPY_DIAGONAL(z, Rcomplex, Matrix_zone)

#undef PACKED_COPY_DIAGONAL

#define UNPACKED_IS_SYMMETRIC(_PREFIX_, _CTYPE_,			\
			      _U_IS_NA_, _L_IS_NOT_NA_, _L_IS_NOT_EQUAL_) \
Rboolean _PREFIX_ ## dense_unpacked_is_symmetric(const _CTYPE_ *x, int n) \
{									\
    int i, j;								\
    R_xlen_t upos = 0, lpos = 0;					\
    for (j = 0; j < n; upos = (lpos += (++j)+1)) {			\
	for (i = j+1; i < n; ++i) {					\
	    upos += n; ++lpos;						\
	    if (_U_IS_NA_) {						\
		if (_L_IS_NOT_NA_)					\
		    return FALSE;					\
	    } else {							\
		if (_L_IS_NOT_EQUAL_)					\
		    return FALSE;					\
	    }								\
	}								\
    }									\
    return TRUE;							\
}

/* ddense_unpacked_is_symmetric() */
UNPACKED_IS_SYMMETRIC(d, double,
		      ISNAN(x[upos]),
		      !ISNAN(x[lpos]),
		      ISNAN(x[lpos]) || x[lpos] != x[upos])
/* ldense_unpacked_is_symmetric() */
UNPACKED_IS_SYMMETRIC(l, int,
		      x[upos] == NA_LOGICAL,
		      x[lpos] != NA_LOGICAL,
		      (x[upos] == 0) ? (x[lpos] != 0) : (x[lpos] == 0))
/* ndense_unpacked_is_symmetric() ... macro still works, if we cheat */
UNPACKED_IS_SYMMETRIC(n, int,
		      x[upos] == 0,
		      x[lpos] != 0,
		      x[lpos] == 0)
/* idense_unpacked_is_symmetric() */
UNPACKED_IS_SYMMETRIC(i, int,
		      x[upos] == NA_INTEGER,
		      x[lpos] != NA_INTEGER,
		      x[lpos] != x[upos])
/* zdense_unpacked_is_symmetric() */
UNPACKED_IS_SYMMETRIC(z, Rcomplex,
		      ISNAN(x[upos].r) || ISNAN(x[upos].i) ,
		      !(ISNAN(x[lpos].r) || ISNAN(x[lpos].i)),
		      ISNAN(x[lpos].r) || ISNAN(x[lpos].i) ||
		      x[upos].r != x[lpos].r || x[upos].i != -x[lpos].i)

#undef UNPACKED_IS_SYMMETRIC

#define UNPACKED_IS_TRIANGULAR(_PREFIX_, _CTYPE_,			\
			       _L_IS_NOT_ZERO_,				\
			       _U_IS_NOT_ZERO_)				\
Rboolean _PREFIX_ ## dense_unpacked_is_triangular(const _CTYPE_ *x, int n, \
						  char uplo)		\
{									\
    int i, j;								\
    if (uplo == 'U') {							\
	for (j = 0; j < n; x += (++j)+1)				\
	    for (i = j+1; i < n; ++i)					\
		if (_L_IS_NOT_ZERO_)					\
		    return FALSE;					\
    } else {								\
	for (j = 0; j < n; x += n-(j++))				\
	    for (i = 0; i < j; ++i)					\
		if (_U_IS_NOT_ZERO_)					\
		    return FALSE;					\
    }									\
    return TRUE;							\
}

/* ddense_unpacked_is_triangular() */
UNPACKED_IS_TRIANGULAR(d, double,
		       ISNAN(*(++x)) || *x != 0.0,
		       ISNAN(*x) || *(x++) != 0.0)
/* idense_unpacked_is_triangular() */
UNPACKED_IS_TRIANGULAR(i, int,
		       *(++x) != 0,
		       *(x++) != 0)
/* zdense_unpacked_is_triangular() */
UNPACKED_IS_TRIANGULAR(z, Rcomplex,
		       ISNAN((*(++x)).r) || (*x).r != 0.0 ||
		       ISNAN((*x).i) || (*x).i != 0.0,
		       ISNAN((*x).r) || (*x).r != 0.0 ||
		       ISNAN((*x).i) || (*(x++)).i != 0.0)

#undef UNPACKED_IS_TRIANGULAR

#define UNPACKED_IS_DIAGONAL(_PREFIX_, _CTYPE_, _OD_IS_NOT_ZERO_)	\
Rboolean _PREFIX_ ## dense_unpacked_is_diagonal(const _CTYPE_ *x, int n) \
{									\
    int i, j;								\
    for (j = 0; j < n; ++j) {						\
	for (i = 0; i < j; ++i)						\
	    if (_OD_IS_NOT_ZERO_)					\
		return FALSE;						\
	++x; /* skip over diagonal */					\
	for (i = j+1; i < n; ++i)					\
	    if (_OD_IS_NOT_ZERO_)					\
		return FALSE;						\
    }									\
    return TRUE;							\
}

/* ddense_unpacked_is_diagonal() */
UNPACKED_IS_DIAGONAL(d, double, ISNAN(*x) || *(x++) != 0.0)
/* idense_unpacked_is_diagonal() */
UNPACKED_IS_DIAGONAL(i, int, *(x++) != 0)
/* zdense_unpacked_is_diagonal() */
UNPACKED_IS_DIAGONAL(z, Rcomplex,
		     ISNAN((*x).r) || (*(x)).r != 0.0 ||
		     ISNAN((*x).i) || (*(x++)).i != 0.0)

#undef UNPACKED_IS_DIAGONAL

#define PACKED_IS_DIAGONAL(_PREFIX_, _CTYPE_,				\
			   _U_IS_NOT_ZERO_,				\
			   _L_IS_NOT_ZERO_)				\
Rboolean _PREFIX_ ## dense_packed_is_diagonal(const _CTYPE_ *x, int n,  \
					      char uplo)		\
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
PACKED_IS_DIAGONAL(d, double,
		   ISNAN(*x) || *(x++) != 0.0,
		   ISNAN(*(++x)) || *x != 0.0)
/* idense_packed_is_diagonal() */
PACKED_IS_DIAGONAL(i, int,
		   *(x++) != 0,
		   *(++x) != 0)
/* zdense_packed_is_diagonal() */
PACKED_IS_DIAGONAL(z, Rcomplex,
		   ISNAN((*x).r) || (*x).r != 0.0 ||
		   ISNAN((*x).i) || (*(x++)).i != 0.0,
		   ISNAN((*(++x)).r) || (*x).r != 0.0 ||
		   ISNAN((*x).i) || (*x).i != 0.0)

#undef PACKED_IS_DIAGONAL

#define PACKED_TRANSPOSE(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_transpose(_CTYPE_ *dest, const _CTYPE_ *src, \
					int n, char uplo)		\
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
PACKED_TRANSPOSE(d, double)
/* idense_packed_transpose() */
PACKED_TRANSPOSE(i, int)
/* zdense_packed_transpose() */
PACKED_TRANSPOSE(z, Rcomplex)
#undef PACKED_TRANSPOSE

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

SEXP unpacked_force(SEXP x, int n, char uplo, char diag)
{
    SEXPTYPE tx = TYPEOF(x);
    if (tx < LGLSXP || tx > CPLXSXP)
	ERROR_INVALID_TYPE("'x'", tx, "unpacked_force");
    R_xlen_t nx = XLENGTH(x);
    SEXP y = PROTECT(allocVector(tx, nx));
    
    if (diag == '\0') {

#define FORCE_SYMMETRIC(_PREFIX_, _CTYPE_, _PTR_)			\
	do {								\
	    _CTYPE_ *px = _PTR_(x), *py = _PTR_(y);			\
	    Memcpy(py, px, nx);						\
	    _PREFIX_ ## dense_unpacked_make_symmetric(py, n, uplo);	\
	} while (0)
	
	switch (tx) {
	case REALSXP:
	    FORCE_SYMMETRIC(d, double, REAL);
	    break;
	case LGLSXP:
	    FORCE_SYMMETRIC(i, int, LOGICAL);
	    break;
	case INTSXP:
	    FORCE_SYMMETRIC(i, int, INTEGER);
	    break;
	case CPLXSXP:
	    FORCE_SYMMETRIC(z, Rcomplex, COMPLEX);
	    break;
	default:
	    break;
	}

#undef FORCE_SYMMETRIC
	
    } else {

#define FORCE_TRIANGULAR(_PREFIX_, _CTYPE_, _PTR_, _ONE_)		\
	do {								\
	    _CTYPE_ *px = _PTR_(x), *py = _PTR_(y);			\
	    Memcpy(py, px, nx);						\
	    _PREFIX_ ## dense_unpacked_make_triangular(py, n, n, uplo, diag); \
	    if (diag != 'N') {						\
		R_xlen_t n1a = (R_xlen_t) n + 1;			\
		for (int j = 0; j < n; ++j, py += n1a)			\
		    *py = _ONE_;					\
	    }								\
	} while (0)
	
	switch (tx) {
	case REALSXP:
	    FORCE_TRIANGULAR(d, double, REAL, 1.0);
	    break;
	case LGLSXP:
	    FORCE_TRIANGULAR(i, int, LOGICAL, 1);
	    break;
	case INTSXP:
	    FORCE_TRIANGULAR(i, int, INTEGER, 1);
	    break;
	case CPLXSXP:
	    FORCE_TRIANGULAR(z, Rcomplex, COMPLEX, Matrix_zone);
	    break;
	default:
	    break;
	}

#undef FORCE_TRIANGULAR
	
    }
    
    UNPROTECT(1);
    return y;
}

char type2kind(SEXPTYPE type)
{
    switch (type) {
    case LGLSXP:
	return 'l';
    case INTSXP:
#ifdef HAVE_PROPER_IMATRIX
	return 'i';
#endif
    case REALSXP:
	return 'd';
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	return 'z';
#endif
    default:
	error(_("unexpected type \"%s\" in 'type2kind()'"), type2char(type)); 
	return '\0';
    }
}

SEXPTYPE kind2type(char kind)
{
    switch (kind) {
    case 'n':
    case 'l':
	return LGLSXP;
#ifdef HAVE_PROPER_IMATRIX
    case 'i':
	return INTSXP;
#endif
    case 'd':
	return REALSXP;
#ifdef HAVE_PROPER_ZMATRIX
    case 'z':
	return CPLXSXP;
#endif
    default:
	error(_("unexpected kind \"%c\" in 'kind2type()'"), kind);
	return NILSXP;
    }
}

size_t kind2size(char kind)
{
    switch (kind) {
    case 'n':
    case 'l':
#ifdef HAVE_PROPER_IMATRIX
    case 'i':
#endif
	return sizeof(int);
    case 'd':
	return sizeof(double);
#ifdef HAVE_PROPER_ZMATRIX
    case 'z':
	return sizeof(Rcomplex);
#endif
    default:
	error(_("unexpected kind \"%c\" in 'kind2size()'"), kind);
	return 0;
    }
}

/* TODO: compare with macros in ./Mutils.h */

#define VALID_NONVIRTUAL						\
/*  0 */ "dgCMatrix", "dgRMatrix", "dgTMatrix", "dgeMatrix",		\
/*  4 */ "dsCMatrix", "dsRMatrix", "dsTMatrix", "dspMatrix", "dsyMatrix", \
/*  9 */ "dtCMatrix", "dtRMatrix", "dtTMatrix", "dtpMatrix", "dtrMatrix", \
/* 14 */ "ddiMatrix", "dsparseVector",					\
/* 16 */ "lgCMatrix", "lgRMatrix", "lgTMatrix", "lgeMatrix",		\
/* 20 */ "lsCMatrix", "lsRMatrix", "lsTMatrix", "lspMatrix", "lsyMatrix", \
/* 25 */ "ltCMatrix", "ltRMatrix", "ltTMatrix", "ltpMatrix", "ltrMatrix", \
/* 30 */ "ldiMatrix", "lsparseVector",					\
/* 32 */ "ngCMatrix", "ngRMatrix", "ngTMatrix", "ngeMatrix",		\
/* 36 */ "nsCMatrix", "nsRMatrix", "nsTMatrix", "nspMatrix", "nsyMatrix", \
/* 41 */ "ntCMatrix", "ntRMatrix", "ntTMatrix", "ntpMatrix", "ntrMatrix", \
/* 46 */ /* "ndiMatrix", */ "nsparseVector",				\
/* 47 */ "igCMatrix", "igRMatrix", "igTMatrix", "igeMatrix",		\
/* 51 */ "isCMatrix", "isRMatrix", "isTMatrix", "ispMatrix", "isyMatrix", \
/* 56 */ "itCMatrix", "itRMatrix", "itTMatrix", "itpMatrix", "itrMatrix", \
/* 61 */ "idiMatrix", "isparseVector",					\
/* 63 */ "zgCMatrix", "zgRMatrix", "zgTMatrix", "zgeMatrix",		\
/* 67 */ "zsCMatrix", "zsRMatrix", "zsTMatrix", "zspMatrix", "zsyMatrix", \
/* 72 */ "ztCMatrix", "ztRMatrix", "ztTMatrix", "ztpMatrix", "ztrMatrix", \
/* 77 */ "zdiMatrix", "zsparseVector",					\
/* 79 */ "indMatrix", "pMatrix"

char Matrix_kind(SEXP obj, int i2d)
{
    if (IS_S4_OBJECT(obj)) {
	static const char *valid[] = { VALID_NONVIRTUAL, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
	    error(_("\"kind\" not yet defined for objects of class \"%s\""),
		  class_P(obj));
	if (ivalid >= 79)
	    return 'n'; /* indMatrix, pMatrix */
	return valid[ivalid][0];
    } else {
	switch (TYPEOF(obj)) {
	case LGLSXP:
	    return 'l';
	case INTSXP:
	    return (i2d) ? 'd' : 'i';
	case REALSXP:
	    return 'd';
	case CPLXSXP:
	    return 'z';
	default:
	    error(_("\"kind\" not yet defined for objects of type \"%s\""),
		  type2char(TYPEOF(obj)));
	    return '\0';
	}
    }
}

SEXP R_Matrix_kind(SEXP obj, SEXP i2d)
{
    char k = Matrix_kind(obj, asLogical(i2d));
    char s[] = { k, '\0' };
    return mkString(s);
}

char Matrix_shape(SEXP obj)
{
    if (!IS_S4_OBJECT(obj))
	error(_("\"shape\" not yet defined for objects of type \"%s\""),
		  type2char(TYPEOF(obj)));
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	error(_("\"shape\" not yet defined for objects of class \"%s\""),
	      class_P(obj));
    if (ivalid >= 79 || valid[ivalid][3] != 'M')
	return 'g'; /* indMatrix, pMatrix, .sparseVector */
    return valid[ivalid][1];
}

SEXP R_Matrix_shape(SEXP obj)
{
    char k = Matrix_shape(obj);
    char s[] = { k, '\0' };
    return mkString(s);
}

/* TODO: compare with macros in ./Mutils.h */

#define VALID_CRTSPARSE				\
"dgCMatrix", "dsCMatrix", "dtCMatrix",		\
"lgCMatrix", "lsCMatrix", "ltCMatrix",		\
"ngCMatrix", "nsCMatrix", "ntCMatrix",		\
"dgRMatrix", "dsRMatrix", "dtRMatrix",		\
"lgRMatrix", "lsRMatrix", "ltRMatrix",		\
"ngRMatrix", "nsRMatrix", "ntRMatrix",		\
"dgTMatrix", "dsTMatrix", "dtTMatrix",		\
"lgTMatrix", "lsTMatrix", "ltTMatrix",		\
"ngTMatrix", "nsTMatrix", "ntTMatrix"

char Matrix_repr(SEXP obj)
{
    if (!IS_S4_OBJECT(obj))
	error(_("\"repr\" not yet defined for objects of type \"%s\""),
	      type2char(TYPEOF(obj)));
    static const char *valid[] = { VALID_CRTSPARSE, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	return '\0'; /* useful to _not_ throw an error, for inheritance tests */
    return valid[ivalid][2];
}

SEXP R_Matrix_repr(SEXP obj)
{
    char k = Matrix_repr(obj);
    if (k == '\0')
	return mkString("");
    char s[] = { k, '\0' };
    return mkString(s);
}

SEXP R_matrix_as_dense(SEXP from, SEXP code, SEXP uplo, SEXP diag)
{
    const char *zzz;
    if (TYPEOF(code) != STRSXP || LENGTH(code) < 1 ||
	(code = STRING_ELT(code, 0)) == NA_STRING ||
	(zzz = CHAR(code))[0] == '\0' ||
	(zzz             )[1] == '\0' ||
	!((zzz[1] == 'g' && (zzz[2] == 'e'                 )) ||
	  (zzz[1] == 't' && (zzz[2] == 'r' || zzz[2] == 'p')) ||
	  (zzz[1] == 's' && (zzz[2] == 'y' || zzz[2] == 'p'))))
	error(_("invalid 'code' to 'R_matrix_as_dense()'"));
    char ul = 'U', di = 'N';
    if (zzz[1] != 'g') {
	if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
	    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
	    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
	    error(_("invalid 'uplo' to 'R_matrix_as_dense()'"));
	if (zzz[1] == 't') {
	    if (TYPEOF(diag) != STRSXP || LENGTH(diag) < 1 ||
		(diag = STRING_ELT(diag, 0)) == NA_STRING ||
		((di = *CHAR(diag)) != 'N' && di != 'U'))
		error(_("invalid 'diag' to 'R_matrix_as_dense()'"));
	}
    }
    return matrix_as_dense(from, zzz, ul, di, 0, 0);
}

SEXP matrix_as_dense(SEXP from, const char *code, char uplo, char diag,
		     int new, int transpose_if_vector)
{
    /* NB: also handing vectors 'from' _without_ a 'dim' attribute */
    SEXPTYPE tf = TYPEOF(from);
    switch (tf) {
    case LGLSXP:
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
#endif
    case REALSXP:
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
#endif
	break;
#ifndef HAVE_PROPER_IMATRIX
    case INTSXP:
	if (!inherits(from, "factor"))
	    break;
#endif
    default:
	if (OBJECT(from))
	    ERROR_INVALID_CLASS(from, "matrix_as_dense");
	else
	    ERROR_INVALID_TYPE("object", tf, "matrix_as_dense");
	break;
    }

    char cl[] = "...Matrix";
    cl[0] = (code[0] == '.') ? type2kind(tf) : code[0];
    cl[1] = code[1];
    cl[2] = code[2];
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    
    SEXP dim, dimnames;
    int *pdim, m, n, doDN, isM = (int) isMatrix(from), nprotect = 1;
    R_xlen_t len = XLENGTH(from);
    
    if (isM) {
	
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(dim = getAttrib(from, R_DimSymbol), &pid);
	pdim = INTEGER(dim);
	if ((m = pdim[0]) != (n = pdim[1]) || n > 0) {
	    if (new > 1)
		REPROTECT(dim = duplicate(dim), pid);
	    SET_SLOT(to, Matrix_DimSym, dim);
	}
	UNPROTECT(1); /* dim */
	
	PROTECT(dimnames = getAttrib(from, R_DimNamesSymbol));
	++nprotect;
	doDN = !isNull(dimnames);
	
    } else {
	
	if (len > INT_MAX)
	    error(_("dimensions cannot exceed 2^31-1"));
	PROTECT(dim = allocVector(INTSXP, 2));
	pdim = INTEGER(dim);
	if (transpose_if_vector) {
	    pdim[0] = m = 1;
	    pdim[1] = n = (int) len;
	} else {
	    pdim[0] = m = (int) len;
	    pdim[1] = n = 1;
	}
	SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP nms = PROTECT(getAttrib(from, R_NamesSymbol));
	++nprotect;
	doDN = !isNull(nms);
	if (doDN) {
	    PROTECT(dimnames = allocVector(VECSXP, 2));
	    ++nprotect;
	    SET_VECTOR_ELT(dimnames, transpose_if_vector ? 1 : 0, nms);
	}
	
    }
    
    if (cl[1] != 'g' && m != n)
	error(_("attempt to construct triangular or symmetric "
		"denseMatrix from non-square matrix"));

    if (doDN) {
	if (cl[1] == 's')
	    set_symmetrized_DimNames(to, dimnames, -1);
	else if (isM && new > 1)
	    set_DimNames(to, dimnames);
	else
	    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    }
    
    if (cl[1] != 'g' && uplo != 'U') {
	SEXP val = PROTECT(mkString("L"));
	SET_SLOT(to, Matrix_uploSym, val);
	UNPROTECT(1); /* val */
    }

    if (cl[1] == 't' && diag != 'N') {
	SEXP val = PROTECT(mkString("U"));
	SET_SLOT(to, Matrix_diagSym, val);
	UNPROTECT(1); /* val */
    }
    
    SEXPTYPE tt = kind2type(cl[0]);
    if (tf != tt) {
	PROTECT(from = coerceVector(from, tt));
	++nprotect;
    }

    SEXP x;
    
    if (cl[2] != 'p') {
	
	if (tf != tt || new < 0 || (new == 0 && !MAYBE_REFERENCED(from)))
	    x = from;
	else {
	    PROTECT(x = allocVector(tt, len));
	    ++nprotect;
	    switch (tt) {
	    case LGLSXP:
		Memcpy(LOGICAL(x), LOGICAL(from), len);
		break;
	    case INTSXP:
		Memcpy(INTEGER(x), INTEGER(from), len);
		break;
	    case REALSXP:
		Memcpy(REAL(x), REAL(from), len);
		break;
	    case CPLXSXP:
		Memcpy(COMPLEX(x), COMPLEX(from), len);
		break;
	    default:
		break;
	    }
	}
	if (!isNull(ATTRIB(x))) {
	    SET_ATTRIB(x, R_NilValue);
	    if (OBJECT(x))
		SET_OBJECT(x, 0);
	}
	
    } else {

	PROTECT(x = allocVector(tt, PM_LENGTH(n)));
	++nprotect;
	
#define PACK(_PREFIX_, _PTR_)						\
	_PREFIX_ ## dense_pack(_PTR_(x), _PTR_(from), n, uplo, diag)
	
	switch (tt) {
	case LGLSXP:
	    PACK(i, LOGICAL);
	    break;
#ifdef HAVE_PROPER_IMATRIX
	case INTSXP:
	    PACK(i, INTEGER);
	    break;
#endif
	case REALSXP:
	    PACK(d, REAL);
	    break;
#ifdef HAVE_PROPER_ZMATRIX
	case CPLXSXP:
	    PACK(z, COMPLEX);
	    break;
#endif
	default:
	    break;
	}

#undef PACK
	
    }
    
    SET_SLOT(to, Matrix_xSym, x);

    UNPROTECT(nprotect);
    return to;
}

SEXP R_dense_as_general(SEXP from, SEXP kind)
{
    char k;
    if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	(kind = STRING_ELT(kind, 0)) == NA_STRING ||
	(k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_dense_as_general()'"));
    return dense_as_general(from, k, 0, 0);
}

/** @brief Coerce `denseMatrix` (and others) to `.geMatrix`.
 *
 *  This utility supports the many `*_{prod,crossprod,tcrossprod,...}`
 *  functions that should work with both classed and unclassed matrices.
 *  It is used in many places for `.geMatrix` ("generalized") dispatch.
 *
 * @param from A `denseMatrix`, a `diagonalMatrix`, a numeric or logical 
 *     `matrix`, or a numeric or logical vector.
 * @param kind A `char` flag, one of `'.'`, `'d'`, `'l'`, and `'n'`, 
 *     indicating the "kind" of `.geMatrix` desired.  A dot `'.'` means 
 *     to preserve the "kind" of `from`.
 * @param new An `int` flag allowing the user to control allocation.
 *     If 0, then usual copy-on-modify rules are employed.  
 *     If less than 0, then the `x` slot of the result is the result 
 *     of modifying in place the `x` slot of `from` (or `from` itself 
 *     if it is a `matrix` or vector), unless the two have different 
 *     lengths or types.
 *     If greater than 0, then the `x` slot of the result is always 
 *     newly allocated.
 *     If greater than 1, then the `Dim` and `Dimnames` slots are also 
 *     always newly allocated (but never the _elements_ of `Dimnames`).
 * @param transpose_if_vector Should length-`n` vectors without a `dim` 
 *     attribute become 1-by-`n` (rather than `n`-by-1) matrices?
 *
 * @return A `.geMatrix`.
 */
SEXP dense_as_general(SEXP from, char kind, int new, int transpose_if_vector)
{
    /* NB: diagonalMatrix is no longer a subclass of denseMatrix 
           but .diMatrix are nonetheless dealt with here ... 
    */
    
    static const char *valid[] = {
	"dgeMatrix", "lgeMatrix", "ngeMatrix",
	"dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix", "ddiMatrix",
	"ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix", "ldiMatrix",
	"ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix", /* "ndiMatrix", */
	""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0) {
	/* dispatch to method for base matrix and base vector */
	char zzz[] = ".ge";
	zzz[0] = kind;
	return matrix_as_dense(from, zzz, '\0', '\0', new, transpose_if_vector);
    }
    
    /* Now handling just denseMatrix and diagonalMatrix ...
       We want to be fast if 'from' is already a .geMatrix,
       and _especially_ fast if it is already a .geMatrix
       of the right "kind" ... 
    */
    
    const char *clf = valid[ivalid];
    if (kind == '.')
	kind = clf[0];
    Rboolean ge = (clf[1] == 'g'), ge0 = ge && kind == clf[0];
    if (ge0 && new <= 0)
	return from;
    
    SEXP to;
    if (ge0)
	PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
    else {
	char clt[] = ".geMatrix";
	clt[0] = kind;
	PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
    }

    SEXP dim;
    PROTECT_INDEX pidA;
    PROTECT_WITH_INDEX(dim = GET_SLOT(from, Matrix_DimSym), &pidA);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (m != n || n > 0) {
	if (new > 1)
	    REPROTECT(dim = duplicate(dim), pidA);
	SET_SLOT(to, Matrix_DimSym, dim);
    }
    UNPROTECT(1); /* dim */
    
    SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else if (new > 1)
	set_DimNames(to, dimnames);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    UNPROTECT(1); /* dimnames */
    
    SEXP x0;
    PROTECT_INDEX pidB;
    PROTECT_WITH_INDEX(x0 = GET_SLOT(from, Matrix_xSym), &pidB);
    if (ge0) {
	REPROTECT(x0 = duplicate(x0), pidB);
	SET_SLOT(to, Matrix_xSym, x0);
	UNPROTECT(2); /* x0, to */
	return to;
    }
    
    SEXPTYPE tf = TYPEOF(x0), tt = kind2type(kind);
    int do_na2one = clf[0] == 'n' && kind != 'n';
    if (ge) {
	if (new == 0 && do_na2one) {
	    /* Try to avoid an allocation ... */
	    R_xlen_t ix, nx = XLENGTH(x0);
	    int *px0 = LOGICAL(x0);
	    for (ix = 0; ix < nx; ++ix)
		if (*(px0++) == NA_LOGICAL)
		    break;
	    do_na2one = (ix < nx);
	}
	if (tf != tt)
	    REPROTECT(x0 = coerceVector(x0, tt), pidB);
	else if (new > 0 || (new == 0 && do_na2one))
	    REPROTECT(x0 = duplicate(x0), pidB);
	if (do_na2one)
	    na2one(x0);
	SET_SLOT(to, Matrix_xSym, x0);
	UNPROTECT(2); /* x0, to */
	return to;
    }
    
    /* Now handling 'from' inheriting from .(tr|sy|tp|sp|di)Matrix ... */

    char ul = 'U', di = 'N';
    if (clf[1] != 'd') {
	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */
    }
    if (clf[1] != 's') {
	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	di = *CHAR(STRING_ELT(diag, 0));
	UNPROTECT(1); /* diag */
    } else if (kind == clf[0]) {
	SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorSym));
	if (LENGTH(factors) > 0)
	    SET_SLOT(to, Matrix_factorSym, factors);
	UNPROTECT(1); /* factors */
    }

    SEXP x1;
    int packed = clf[1] == 'd' || clf[2] == 'p';
    if (!packed) {
	/* (tr|sy)->ge */
	if (tf != tt)
	    REPROTECT(x1 = coerceVector(x0, tt), pidB);
	else if (new >= 0)
	    REPROTECT(x1 = duplicate(x0), pidB);
	else
	    x1 = x0;
    } else {
	/* (tp|sp|di)->ge */
	if ((double) n * n > R_XLEN_T_MAX)
	    error(_("attempt to allocate vector of length exceeding "
		    "R_XLEN_T_MAX"));
	if (tf != tt)
	    REPROTECT(x0 = coerceVector(x0, tt), pidB);
	PROTECT(x1 = allocVector(tt, (R_xlen_t) n * n));
    }

#define AS_GE(_PREFIX_, _CTYPE_, _PTR_)					\
    do {								\
	_CTYPE_ *px1 = _PTR_(x1);					\
	if (clf[1] == 'd') {						\
	    /* di->ge */						\
	    Memzero(px1, (R_xlen_t) n * n);				\
	    _PREFIX_ ## dense_unpacked_copy_diagonal(			\
		px1, _PTR_(x0), n, n, ul /* unused */, di);		\
	} else {							\
	    /* (tr|sy|tp|sp)->ge */					\
	    if (clf[2] == 'p')						\
		_PREFIX_ ## dense_unpack(px1, _PTR_(x0), n, ul, di);	\
	    if (clf[1] == 't')						\
		_PREFIX_ ## dense_unpacked_make_triangular(px1, n, n, ul, di); \
	    else							\
		_PREFIX_ ## dense_unpacked_make_symmetric(px1, n, ul);	\
	}								\
    } while (0)

    switch (kind) {
    case 'n':
    case 'l':
	AS_GE(i, int, LOGICAL);
	break;
    case 'i':
	AS_GE(i, int, INTEGER);
	break;
    case 'd': 
	AS_GE(d, double, REAL);
	break;
    case 'z':
	AS_GE(z, Rcomplex, COMPLEX);
	break;
    default:
	break;
    }

#undef AS_GE

    if (do_na2one)
	na2one(x1);
    SET_SLOT(to, Matrix_xSym, x1);

    if (packed)
	UNPROTECT(3); /* x1, x0, to */
    else
	UNPROTECT(2); /* x1,     to */
    return to;
}


/* "General" purpose ================================================ */

SEXP R_index_triangle(SEXP n_, SEXP upper_, SEXP diag_, SEXP packed_)
{
    int n = asInteger(n_), packed = asLogical(packed_);
    double nn = (double) n * n, nx = (packed) ? nn : 0.5 * (nn + n);
    if (nx > R_XLEN_T_MAX)
	error(_("cannot index a vector of length exceeding R_XLEN_T_MAX"));
    SEXP r;
    int i, j, upper = asLogical(upper_), diag = asLogical(diag_);
    double nr = (diag) ? 0.5 * (nn + n) : 0.5 * (nn - n);
    if (nx > INT_MAX) {
	
	PROTECT(r = allocVector(REALSXP, (R_xlen_t) nr));
	double k = 1.0, *pr = REAL(r);

#define DO_INDEX(_ONE_, _NR_)				\
	do {						\
	    if (packed) {				\
		if (diag) {				\
		    while (k <= _NR_) {			\
			*(pr++) = k;			\
			k += _ONE_;			\
		    }					\
		} else if (upper) {			\
		    for (j = 0; j < n; ++j) {		\
			for (i = 0; i < j; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
			k += _ONE_;			\
		    }					\
		} else {				\
		    for (j = 0; j < n; ++j) {		\
			k += 1.0;			\
			for (i = j+1; i < n; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
		    }					\
		}					\
	    } else if (diag) {				\
		if (upper) {				\
		    for (j = 0; j < n; ++j) {		\
			for (i = 0; i <= j; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
			k += n-j-1;			\
		    }					\
		} else {				\
		    for (j = 0; j < n; ++j) {		\
			k += j;				\
			for (i = j; i < n; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
		    }					\
		}					\
	    } else {					\
		if (upper) {				\
		    for (j = 0; j < n; ++j) {		\
			for (i = 0; i < j; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
			k += n-j;			\
		    }					\
		} else {				\
		    for (j = 0; j < n; ++j) {		\
			k += j+1;			\
			for (i = j+1; i < n; ++i) {	\
			    *(pr++) = k;		\
			    k += _ONE_;			\
			}				\
		    }					\
		}					\
	    }						\
	} while (0)

	DO_INDEX(1.0, nr);
	
    } else {
	
	PROTECT(r = allocVector(INTSXP, (R_xlen_t) nr));
	int k = 1, nr_ = (int) nr, *pr = INTEGER(r);

	DO_INDEX(1, nr_);

#undef DO_INDEX
	
    }
    
    UNPROTECT(1);
    return r;
}

SEXP R_index_diagonal(SEXP n_, SEXP upper_, SEXP packed_)
{
    int n = asInteger(n_), packed = asLogical(packed_);
    double nn = (double) n * n, nx = (packed) ? nn : 0.5 * (nn + n);
    if (nx > R_XLEN_T_MAX)
	error(_("cannot index a vector of length exceeding R_XLEN_T_MAX"));
    SEXP r;
    int j, upper = (packed) ? asLogical(upper_) : NA_LOGICAL;
    if (nx > INT_MAX) {
	
	PROTECT(r = allocVector(REALSXP, n));
	double k = 1.0, *pr = REAL(r);

#define DO_INDEX				\
	do {					\
	    if (!packed) {			\
		for (j = 0; j < n; ++j) {	\
		    *(pr++) = k;		\
		    k += n+1;			\
		}				\
	    } else if (upper) {			\
		for (j = 0; j < n; ++j) {	\
		    *(pr++) = k;		\
		    k += j+2;			\
		}				\
	    } else {				\
		for (j = 0; j < n; ++j) {	\
		    *(pr++) = k;		\
		    k += n-j;			\
		}				\
	    }					\
	} while (0)

	DO_INDEX;
	
    } else {

	PROTECT(r = allocVector(INTSXP, n));
	int k = 1, *pr = INTEGER(r);
	DO_INDEX;

#undef DO_INDEX
	
    }
    
    UNPROTECT(1);
    return r;
}

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax)
{
    int do_countNA = asLogical(countNA);
    R_xlen_t n = XLENGTH(x), nnz = 0;
    double n_ = asReal(nnzmax);
    if (!ISNAN(n_) && n_ >= 0.0 && n_ < (double) n)
	n = (R_xlen_t) n_;

#define DO_NNZ(_CTYPE_, _PTR_, _NA_, _NZ_, _STRICTLY_NZ_)	\
    do {							\
	_CTYPE_ *px = _PTR_(x);					\
	if (do_countNA == NA_LOGICAL) {				\
	    while (n-- > 0) {					\
		if (_NA_(*px))					\
		    return ScalarInteger(NA_INTEGER);		\
		if (_NZ_(*px))					\
		    ++nnz;					\
		++px;						\
	    }							\
	} else if (do_countNA != 0) {				\
	    while (n-- > 0) {					\
		if (_NZ_(*px))					\
		    ++nnz;					\
		++px;						\
	    }							\
	} else {						\
	    while (n-- > 0) {					\
		if (_STRICTLY_NZ_(*px))				\
		    ++nnz;					\
		++px;						\
	    }							\
	}							\
    } while (0)

    switch (TYPEOF(x)) {
    case LGLSXP:
	DO_NNZ(int, LOGICAL,
	       ISNA_LOGICAL, ISNZ_LOGICAL, STRICTLY_ISNZ_LOGICAL);
	break;
    case INTSXP:
    	DO_NNZ(int, INTEGER,
	       ISNA_INTEGER, ISNZ_INTEGER, STRICTLY_ISNZ_INTEGER);
	break;
    case REALSXP:
    	DO_NNZ(double, REAL,
	       ISNA_REAL, ISNZ_REAL, STRICTLY_ISNZ_REAL);
	break;
    case CPLXSXP:
    	DO_NNZ(Rcomplex, COMPLEX,
	       ISNA_COMPLEX, ISNZ_COMPLEX, STRICTLY_ISNZ_COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x'", TYPEOF(x), "R_nnz");
    }

#undef DO_NNZ

    return (nnz <= INT_MAX)
	? ScalarInteger((int) nnz) : ScalarReal((double) nnz);
}

void conjugate(SEXP x)
{
    Rcomplex *px = COMPLEX(x);
    R_xlen_t nx = XLENGTH(x);
    while (nx--) {
	(*px).i = -(*px).i;
	++px;
    }
    return;
}

void zeroRe(SEXP x)
{
    Rcomplex *px = COMPLEX(x);
    R_xlen_t nx = XLENGTH(x);
    while (nx--) {
	if (!ISNAN((*px).r))
	    (*px).r = 0.0;
	++px;
    }
    return;
}

void zeroIm(SEXP x)
{
    Rcomplex *px = COMPLEX(x);
    R_xlen_t nx = XLENGTH(x);
    while (nx--) {
	if (!ISNAN((*px).i))
	    (*px).i = 0.0;
	++px;
    }
    return;
}

void na2one(SEXP x)
{
    R_xlen_t i, n = XLENGTH(x);
    switch (TYPEOF(x)) {
    case LGLSXP:
    {
	int *px = LOGICAL(x);
	for (i = 0; i < n; ++i, ++px)
	    if (*px == NA_LOGICAL)
		*px = 1;
	break;
    }
    case INTSXP:
    {
	int *px = INTEGER(x);
	for (i = 0; i < n; ++i, ++px)
	    if (*px == NA_INTEGER)
		*px = 1;
	break;
    }
    case REALSXP:
    {
	double *px = REAL(x);
	for (i = 0; i < n; ++i, ++px)
	    if (ISNAN(*px))
		*px = 1.0;
	break;
    }
    case CPLXSXP:
    {
	Rcomplex *px = COMPLEX(x);
	for (i = 0; i < n; ++i, ++px)
	    if (ISNAN((*px).r) || ISNAN((*px).i))
		*px = Matrix_zone;
	break;
    }
    default:
	ERROR_INVALID_TYPE("'x'", TYPEOF(x), "na2one");
	break;
    }
    return;
}

SEXP v2spV(SEXP from)
{
    SEXPTYPE tx = TYPEOF(from);
    SEXP to = NULL, length = NULL, i = NULL, x = NULL;
    R_xlen_t n_ = XLENGTH(from);

#define V2SPV(_KIND_, _NZ_,						\
	      _CTYPE1_, _SEXPTYPE1_, _PTR1_,				\
	      _CTYPE2_, _SEXPTYPE2_, _PTR2_)				\
    do {								\
	PROTECT(to = NEW_OBJECT_OF_CLASS(#_KIND_ "sparseVector"));	\
	_CTYPE1_ *py = _PTR1_(from);					\
	for (k = 0; k < n; ++k)						\
	    if (_NZ_(py[k]))						\
		++nnz;							\
	PROTECT(i = allocVector(_SEXPTYPE2_, nnz));			\
	PROTECT(x = allocVector(_SEXPTYPE1_, nnz));			\
	_CTYPE2_ *pi = _PTR2_(i);					\
	_CTYPE1_ *px = _PTR1_(x);					\
	for (k = 0; k < n; ++k) {					\
	    if (_NZ_(py[k])) {						\
		*(pi++) = (_CTYPE2_) (k + 1);				\
		*(px++) = py[k];					\
	    }								\
	}								\
    } while (0)

#define V2SPV_CASES(_CTYPE2_, _SEXPTYPE2_, _PTR2_)			\
    do {								\
	switch (tx) {							\
	case LGLSXP:							\
	    V2SPV(l, ISNZ_LOGICAL, int, LGLSXP, LOGICAL,		\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	case INTSXP:							\
	    V2SPV(i, ISNZ_INTEGER, int, INTSXP, INTEGER,		\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	case REALSXP:							\
	    V2SPV(d, ISNZ_REAL, double, REALSXP, REAL,			\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	case CPLXSXP:							\
	    V2SPV(z, ISNZ_COMPLEX, Rcomplex, CPLXSXP, COMPLEX,		\
		  _CTYPE2_, _SEXPTYPE2_, _PTR2_);			\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("object", tx, "v2spV");			\
	    break;							\
	}								\
    } while (0)
    
    if (n_ <= INT_MAX) {
	int k, n = (int) n_, nnz = 0;
	PROTECT(length = ScalarInteger(n));
	V2SPV_CASES(int, INTSXP, INTEGER);
    } else {
	R_xlen_t k, n = n_, nnz = 0;
	PROTECT(length = ScalarReal((double) n));
	V2SPV_CASES(double, REALSXP, REAL);
    }

#undef V2SPV_CASES
#undef V2SPV

    SET_SLOT(to, Matrix_lengthSym, length);
    SET_SLOT(to, Matrix_iSym, i);
    SET_SLOT(to, Matrix_xSym, x);
    
    UNPROTECT(4); /* x, i, length, to */
    return to;
}

/* That both 's1' and 's2' are STRSXP of length at least 'n' must be 
   checked by the caller ... see, e.g., symmetricMatrix_validate() above
*/
Rboolean equal_string_vectors(SEXP s1, SEXP s2, int n)
{
    /* Only check the first 'n' elements, even if 's1' or 's2' is longer ...
    
       Note that 'R_compute_identical()' in src/main/identical.c
       is careful to distinguish between NA_STRING and "NA" in STRSXP, 
       but we need not be here ...
       
       MJ: Why not?
    */
	
    for (int i = 0; i < n; ++i)
	if (strcmp(CHAR(STRING_ELT(s1, i)), CHAR(STRING_ELT(s2, i))) != 0)
	    return FALSE;
    return TRUE;
}

/* That 'valid' is a STRSXP must be checked by the caller ...
   NOTE: ./Mutils.h has
         int Matrix_check_class_(char *x, const char **valid);
	 int Matrix_check_class(SEXP x, const char **valid);
	 ... and this is yet another variant ...
*/
R_xlen_t strmatch(const char *x, SEXP valid)
{
    R_xlen_t n = XLENGTH(valid);
    for (R_xlen_t i = 0; i < n; ++i)
	if (strcmp(x, CHAR(STRING_ELT(valid, i))) == 0)
	    return i;
    return -1;
}

SEXP append_to_named_list(SEXP x, const char *nm, SEXP val)
{
    PROTECT(val);
    R_xlen_t n = XLENGTH(x);
    SEXP y = PROTECT(allocVector(VECSXP, n + 1)),
	ny = PROTECT(allocVector(STRSXP, n + 1)),
	nval = PROTECT(mkChar(nm));
    if (n > 0) {
	SEXP nx = PROTECT(getAttrib(x, R_NamesSymbol));
	R_xlen_t i;
	for (i = 0; i < n; ++i) {
	    SET_VECTOR_ELT( y, i, VECTOR_ELT( x, i));
	    SET_STRING_ELT(ny, i, STRING_ELT(nx, i));
	}
	UNPROTECT(1);
    }
    SET_VECTOR_ELT( y, n,  val);
    SET_STRING_ELT(ny, n, nval);
    setAttrib(y, R_NamesSymbol, ny);
    UNPROTECT(4);
    return y;
}


/* ================================================================== */
/* ================================================================== */
    
/* La_norm_type() and La_rcond_type() have been in src/include/R_ext/Lapack.h
   and later in src/modules/lapack/Lapack.c but have still not been available 
   to package writers ...
*/
char La_norm_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* aliases */
    else if (typup == 'E')
	typup = 'F';
    else if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error(_("argument type[1]='%s' must be one of 'M','1','O','I','F', or 'E'"),
	      typstr);
    return typup;
}

char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error(_("argument type[1]='%s' must be one of '1','O', or 'I'"),
	      typstr);
    return typup; /* 'O' or 'I' */
}

SEXP as_det_obj(double mod, int log, int sign)
{
    SEXP det = PROTECT(allocVector(VECSXP, 2)),
	nms = PROTECT(allocVector(STRSXP, 2)),
	val = PROTECT(ScalarReal(mod));

    setAttrib(det, R_NamesSymbol, nms);
    SET_STRING_ELT(nms, 0, mkChar("modulus"));
    SET_STRING_ELT(nms, 1, mkChar("sign"));
    setAttrib(val, install("logarithm"), ScalarLogical(log));
    SET_VECTOR_ELT(det, 0, val);
    SET_VECTOR_ELT(det, 1, ScalarInteger(sign));
    setAttrib(det, R_ClassSymbol, mkString("det"));
    UNPROTECT(3);
    return det;
}

/* MJ: no longer needed ... prefer more general (un)?packedMatrix_diag_[gs]et() */
#if 0

/**
 * Copy the diagonal elements of the packed denseMatrix x to dest
 *
 * @param dest vector of length ncol(x)
 * @param x (pointer to) a "d?pMatrix" object
 * @param n number of columns in the matrix represented by x
 *
 * @return dest
 */
void d_packed_getDiag(double *dest, SEXP x, int n)
{
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

#define END_packed_getDiag						\
    int j, pos = 0;							\
									\
    if (*uplo_P(x) == 'U') {						\
	for(pos= 0, j=0; j < n; pos += 1+(++j))	 dest[j] = xx[pos];	\
    } else {								\
	for(pos= 0, j=0; j < n; pos += (n - j), j++) dest[j] = xx[pos]; \
    }									\
    return

    END_packed_getDiag;
}

void l_packed_getDiag(int *dest, SEXP x, int n)
{
    int *xx = LOGICAL(GET_SLOT(x, Matrix_xSym));

    END_packed_getDiag;
}

#undef END_packed_getDiag

/** diag(x) <- D  for   x a  <dspMatrix>  or dppMatrix, ..etc
 */
SEXP d_packed_setDiag(double *diag, int l_d, SEXP x, int n)
{
#define SET_packed_setDiag				\
    SEXP ret = PROTECT(duplicate(x)),			\
	r_x = GET_SLOT(ret, Matrix_xSym);		\
    Rboolean d_full = (l_d == n);			\
    if (!d_full && l_d != 1)				\
	error(_("replacement diagonal has wrong length"))

#define END_packed_setDiag						\
    int j, pos = 0;							\
									\
    if (*uplo_P(x) == 'U') {						\
	if(d_full)							\
	    for(pos= 0, j=0; j < n; pos += 1+(++j))	 xx[pos] = diag[j]; \
	else /* l_d == 1 */						\
	    for(pos= 0, j=0; j < n; pos += 1+(++j))	 xx[pos] = *diag; \
    } else {								\
	if(d_full)							\
	    for(pos= 0, j=0; j < n; pos += (n - j), j++) xx[pos] = diag[j]; \
	else /* l_d == 1 */						\
	    for(pos= 0, j=0; j < n; pos += (n - j), j++) xx[pos] = *diag; \
    }									\
    UNPROTECT(1);							\
    return ret

    SET_packed_setDiag; double *xx = REAL(r_x);
    END_packed_setDiag;
}

SEXP l_packed_setDiag(int *diag, int l_d, SEXP x, int n)
{
    SET_packed_setDiag; int *xx = LOGICAL(r_x);
    END_packed_setDiag;
}

#define tr_END_packed_setDiag						\
    if (*diag_P(x) == 'U') { /* uni-triangular */			\
	/* after setting, typically is not uni-triangular anymore: */	\
	SEXP ch_N = PROTECT(mkChar("N"));				\
	SET_STRING_ELT(GET_SLOT(ret, Matrix_diagSym), 0, ch_N);		\
	UNPROTECT(1);							\
    }									\
    END_packed_setDiag

SEXP tr_d_packed_setDiag(double *diag, int l_d, SEXP x, int n)
{
    SET_packed_setDiag; double *xx = REAL(r_x);
    tr_END_packed_setDiag;
}

SEXP tr_l_packed_setDiag(int *diag, int l_d, SEXP x, int n)
{
    SET_packed_setDiag; int *xx = LOGICAL(r_x);
    tr_END_packed_setDiag;
}

#undef SET_packed_setDiag
#undef END_packed_setDiag
#undef tr_END_packed_setDiag

void tr_d_packed_getDiag(double *dest, SEXP x, int n)
{
    if (*diag_P(x) == 'U') {
	for (int j = 0; j < n; j++) dest[j] = 1.;
    } else {
	d_packed_getDiag(dest, x, n);
    }
    return;
}

void tr_l_packed_getDiag(   int *dest, SEXP x, int n)
{
    if (*diag_P(x) == 'U')
	for (int j = 0; j < n; j++) dest[j] = 1;
    else
	l_packed_getDiag(dest, x, n);
    return;
}

/* These two *_addDiag() were unused and not replaced */
SEXP d_packed_addDiag(double *diag, int l_d, SEXP x, int n)
{
    SEXP ret = PROTECT(duplicate(x)),
	r_x = GET_SLOT(ret, Matrix_xSym);
    double *xx = REAL(r_x);
    int j, pos = 0;

    if (*uplo_P(x) == 'U') {
	for(pos= 0, j=0; j < n; pos += 1+(++j))	     xx[pos] += diag[j];
    } else {
	for(pos= 0, j=0; j < n; pos += (n - j), j++) xx[pos] += diag[j];
    }
    UNPROTECT(1);
    return ret;
}

SEXP tr_d_packed_addDiag(double *diag, int l_d, SEXP x, int n)
{
    SEXP ret = PROTECT(d_packed_addDiag(diag, l_d, x, n));
    if (*diag_P(x) == 'U') { /* uni-triangular */
	SEXP ch_N = PROTECT(mkChar("N"));
	SET_STRING_ELT(GET_SLOT(ret, Matrix_diagSym), 0, ch_N);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return ret;
}

#endif /* MJ */

#if 0 /* unused */

double get_double_by_name(SEXP obj, char *nm)
{
    SEXP nms = PROTECT(getAttrib(obj, R_NamesSymbol));
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error(_("object must be a named, numeric vector"));
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    UNPROTECT(1);
	    return REAL(obj)[i];
	}
    }
    UNPROTECT(1);
    return R_NaReal;
}

SEXP set_double_by_name(SEXP obj, double val, char *nm)
{
    SEXP nms = PROTECT(getAttrib(obj, R_NamesSymbol));
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error(_("object must be a named, numeric vector"));
    // case 1:  replace existing entry named  <nm>
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    REAL(obj)[i] = val;
	    UNPROTECT(1);
	    return obj;
	}
    }
    // case 2:  no such name --> add new entry with that name at end of vec
    {
	SEXP nx = PROTECT(allocVector(REALSXP, len + 1)),
	    nnms = allocVector(STRSXP, len + 1);

	setAttrib(nx, R_NamesSymbol, nnms);
	for (i = 0; i < len; i++) {
	    REAL(nx)[i] = REAL(obj)[i];
	    SET_STRING_ELT(nnms, i, duplicate(STRING_ELT(nms, i)));
	}
	REAL(nx)[len] = val;
	SET_STRING_ELT(nnms, len, mkChar(nm));
	UNPROTECT(2);
	return nx;
    }
}

/* useful for all the ..CMatrix classes (and ..R by [0] <-> [1]) */
SEXP CMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}

SEXP new_dgeMatrix(int nrow, int ncol)
{
    SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
	 ad = PROTECT(allocVector(INTSXP, 2));

    INTEGER(ad)[0] = nrow;
    INTEGER(ad)[1] = ncol;
    SET_SLOT(ans, Matrix_DimSym, ad);
    SET_SLOT(ans, Matrix_DimNamesSym, allocVector(VECSXP, 2));
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, (R_xlen_t) nrow * ncol);

    UNPROTECT(2);
    return ans;
}

#endif /* unused */

SEXP Matrix_expand_pointers(SEXP pP)
{
    int n = length(pP) - 1;
    int *p = INTEGER(pP);
    SEXP ans = PROTECT(allocVector(INTSXP, p[n]));

    expand_cmprPt(n, p, INTEGER(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param ij: 2-column integer matrix
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd(SEXP ij, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
    SEXP ans;
    int *ij_di = NULL, n, nprot=1;
    Rboolean check_bounds = asLogical(chk_bnds), one_ind = asLogical(orig_1);

    if(TYPEOF(di) != INTSXP) {di = PROTECT(coerceVector(di, INTSXP)); nprot++; }
    if(TYPEOF(ij) != INTSXP) {ij = PROTECT(coerceVector(ij, INTSXP)); nprot++; }
    if(!isMatrix(ij) ||
       (ij_di = INTEGER(getAttrib(ij, R_DimSymbol)))[1] != 2)
	error(_("Argument ij must be 2-column integer matrix"));
    n = ij_di[0];
    int *Di = INTEGER(di), *IJ = INTEGER(ij),
	*j_ = IJ+n;/* pointer offset! */

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];
#define do_ii_FILL(_i_, _j_)						\
	int i;								\
	if(check_bounds) {						\
	    for(i=0; i < n; i++) {					\
		if(_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER)	\
		    ii[i] = NA_INTEGER;					\
		else {							\
		    register int i_i, j_i;				\
	            if(one_ind) { i_i = _i_[i]-1; j_i = _j_[i]-1; }	\
	            else        { i_i = _i_[i]  ; j_i = _j_[i]  ; }	\
		    if(i_i < 0 || i_i >= Di[0])				\
			error(_("subscript 'i' out of bounds in M[ij]")); \
		    if(j_i < 0 || j_i >= Di[1])				\
			error(_("subscript 'j' out of bounds in M[ij]")); \
		    ii[i] = i_i + j_i * nr;				\
		}							\
	    }								\
	} else {							\
	    for(i=0; i < n; i++)					\
		ii[i] = (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER)	\
		    ? NA_INTEGER					\
 	            : (one_ind ? ((_i_[i]-1) + (_j_[i]-1)*nr)		\
	                       :   _i_[i]    +  _j_[i]   *nr);		\
	}

	do_ii_FILL(IJ, j_);
    } else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(IJ, j_);
    }
    UNPROTECT(nprot);
    return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param i: integer vector
 * @param j: integer vector of same length as 'i'
 * @param orig_1: logical: if TRUE, "1-origin" otherwise "0-origin"
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
    SEXP ans;
    int n = LENGTH(i), nprot = 1;
    Rboolean check_bounds = asLogical(chk_bnds), one_ind = asLogical(orig_1);

    if(TYPEOF(di)!= INTSXP) {di = PROTECT(coerceVector(di,INTSXP)); nprot++; }
    if(TYPEOF(i) != INTSXP) { i = PROTECT(coerceVector(i, INTSXP)); nprot++; }
    if(TYPEOF(j) != INTSXP) { j = PROTECT(coerceVector(j, INTSXP)); nprot++; }
    if(LENGTH(j) != n)
	error(_("i and j must be integer vectors of the same length"));
    int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];

	do_ii_FILL(i_, j_);
    } else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(i_, j_);
    }
    UNPROTECT(nprot);
    return ans;
}
#undef do_ii_FILL

// Almost "Cut n Paste" from ...R../src/main/array.c  do_matrix() :
// used in ../R/Matrix.R as
//
// .External(Mmatrix,
//	     data, nrow, ncol, byrow, dimnames,
//	     missing(nrow), missing(ncol))
SEXP Mmatrix(SEXP args)
{
    SEXP vals, ans, snr, snc, dimnames;
    int nr = 1, nc = 1, byrow, miss_nr, miss_nc;
    R_xlen_t lendat;

    args = CDR(args); /* skip 'name' */
    vals = CAR(args); args = CDR(args);
    /* Supposedly as.vector() gave a vector type, but we check */
    switch(TYPEOF(vals)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
	case STRSXP:
	case RAWSXP:
	case EXPRSXP:
	case VECSXP:
	    break;
	default:
	    error(_("'data' must be of a vector type"));
    }
    lendat = XLENGTH(vals);
    snr = CAR(args); args = CDR(args);
    snc = CAR(args); args = CDR(args);
    byrow = asLogical(CAR(args)); args = CDR(args);
    if (byrow == NA_INTEGER)
	error(_("invalid '%s' argument"), "byrow");
    dimnames = CAR(args);
    args = CDR(args);
    miss_nr = asLogical(CAR(args)); args = CDR(args);
    miss_nc = asLogical(CAR(args));

    if (!miss_nr) {
	if (!isNumeric(snr)) error(_("non-numeric matrix extent"));
	nr = asInteger(snr);
	if (nr == NA_INTEGER)
	    error(_("invalid 'nrow' value (too large or NA)"));
	if (nr < 0)
	    error(_("invalid 'nrow' value (< 0)"));
    }
    if (!miss_nc) {
	if (!isNumeric(snc)) error(_("non-numeric matrix extent"));
	nc = asInteger(snc);
	if (nc == NA_INTEGER)
	    error(_("invalid 'ncol' value (too large or NA)"));
	if (nc < 0)
	    error(_("invalid 'ncol' value (< 0)"));
    }
    if (miss_nr && miss_nc) {
	if (lendat > INT_MAX) error("data is too long");
	nr = (int) lendat;
    } else if (miss_nr) {
	if (lendat > (double) nc * INT_MAX) error("data is too long");
	nr = (int) ceil((double) lendat / (double) nc);
    } else if (miss_nc) {
	if (lendat > (double) nr * INT_MAX) error("data is too long");
	nc = (int) ceil((double) lendat / (double) nr);
    }

    if(lendat > 0) {
	R_xlen_t nrc = (R_xlen_t) nr * nc;
	if (lendat > 1 && nrc % lendat != 0) {
	    if (((lendat > nr) && (lendat / nr) * nr != lendat) ||
		((lendat < nr) && (nr / lendat) * lendat != nr))
		warning(_("data length [%d] is not a sub-multiple "
			  "or multiple of the number of rows [%d]"),
			lendat, nr);
	    else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
		     ((lendat < nc) && (nc / lendat) * lendat != nc))
		warning(_("data length [%d] is not a sub-multiple "
			  "or multiple of the number of columns [%d]"),
			lendat, nc);
	}
	else if ((lendat > 1) && (nrc == 0)){
	    warning(_("data length exceeds size of matrix"));
	}
    }

#ifndef LONG_VECTOR_SUPPORT
   if ((double)nr * (double)nc > INT_MAX)
	error(_("too many elements specified"));
#endif

    PROTECT(ans = allocMatrix(TYPEOF(vals), nr, nc));
    if(lendat) {
	if (isVector(vals))
	    copyMatrix(ans, vals, byrow);
	else
	    copyListMatrix(ans, vals, byrow);
    } else if (isVector(vals)) { /* fill with NAs */
	R_xlen_t N = (R_xlen_t) nr * nc, i;
	switch(TYPEOF(vals)) {
	case STRSXP:
	    for (i = 0; i < N; i++)
		SET_STRING_ELT(ans, i, NA_STRING);
	    break;
	case LGLSXP:
	    for (i = 0; i < N; i++)
		LOGICAL(ans)[i] = NA_LOGICAL;
	    break;
	case INTSXP:
	    for (i = 0; i < N; i++)
		INTEGER(ans)[i] = NA_INTEGER;
	    break;
	case REALSXP:
	    for (i = 0; i < N; i++)
		REAL(ans)[i] = NA_REAL;
	    break;
	case CPLXSXP:
	{
	    Rcomplex zna = { NA_REAL, 0.0 };
	    for (i = 0; i < N; i++)
		COMPLEX(ans)[i] = zna;
	    break;
	}
	case RAWSXP:
	    // FIXME:  N may overflow size_t !!
	    memset(RAW(ans), 0, N);
	    break;
	default:
	    /* don't fill with anything */
	    ;
	}
    }
    if(!isNull(dimnames)&& length(dimnames) > 0)
	ans = dimnamesgets(ans, dimnames);
    UNPROTECT(1);
    return ans;
}

/**
 * From the two 'x' slots of two dense matrices a and b,
 * compute the 'x' slot of rbind(a, b)
 *
 * Currently, an auxiliary only for setMethod rbind2(<denseMatrix>, <denseMatrix>)
 * in ../R/bind2.R
 *
 * @param a
 * @param b
 *
 * @return
 */SEXP R_rbind2_vector(SEXP a, SEXP b) {
    int *d_a = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*d_b = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	n1 = d_a[0], m = d_a[1],
	n2 = d_b[0];
    if(d_b[1] != m)
	error(_("the number of columns differ in R_rbind2_vector: %d != %d"),
	      m, d_b[1]);
    SEXP
	a_x = GET_SLOT(a, Matrix_xSym),
	b_x = GET_SLOT(b, Matrix_xSym);
    int nprot = 1;
    // Care: can have "ddenseMatrix" "ldenseMatrix" or "ndenseMatrix"
    if(TYPEOF(a_x) != TYPEOF(b_x)) { // choose the "common type"
	// Now know: either LGLSXP or REALSXP. FIXME for iMatrix, zMatrix,..
	if(TYPEOF(a_x) != REALSXP) {
	    a_x = PROTECT(duplicate(coerceVector(a_x, REALSXP))); nprot++;
	} else if(TYPEOF(b_x) != REALSXP) {
	    b_x = PROTECT(duplicate(coerceVector(b_x, REALSXP))); nprot++;
	}
    }

    SEXP ans = PROTECT(allocVector(TYPEOF(a_x), m * (n1 + n2)));
    int ii = 0;
    switch(TYPEOF(a_x)) {
    case LGLSXP: {
	int
	    *r = LOGICAL(ans),
	    *ax= LOGICAL(a_x),
	    *bx= LOGICAL(b_x);

#define COPY_a_AND_b_j					\
	for(int j=0; j < m; j++) {			\
	    Memcpy(r+ii, ax+ j*n1, n1); ii += n1;	\
	    Memcpy(r+ii, bx+ j*n2, n2); ii += n2;	\
	} ; break

	COPY_a_AND_b_j;
    }
    case REALSXP: {
	double
	    *r = REAL(ans),
	    *ax= REAL(a_x),
	    *bx= REAL(b_x);

	COPY_a_AND_b_j;
    }
    } // switch
    UNPROTECT(nprot);
    return ans;
}

#define TRUE_  ScalarLogical(1)
#define FALSE_ ScalarLogical(0)

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// all0     <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
// allFalse <- function(x) !any(x) && !any(is.na(x)) ## ~= all0
SEXP R_all0(SEXP x) {
    if (!isVectorAtomic(x)) {
	if(length(x) == 0) return TRUE_;
	// Typically S4.  TODO: Call the R code above, instead!
	error(_("Argument must be numeric-like atomic vector"));
    }
    R_xlen_t i, n = XLENGTH(x);
    if(n == 0) return TRUE_;

    switch(TYPEOF(x)) {
    case LGLSXP: {
	int *xx = LOGICAL(x);
	for(i=0; i < n; i++)
	    if(xx[i] == NA_LOGICAL || xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    case INTSXP: {
	int *xx = INTEGER(x);
	for(i=0; i < n; i++)
	    if(xx[i] == NA_INTEGER || xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    case REALSXP: {
	double *xx = REAL(x);
	for(i=0; i < n; i++)
	    if(ISNAN(xx[i]) || xx[i] != 0.) return FALSE_;
	return TRUE_;
    }
    case RAWSXP: {
	unsigned char *xx = RAW(x);
	for(i=0; i < n; i++)
	    if(xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    }
    error(_("Argument must be numeric-like atomic vector"));
    return R_NilValue; // -Wall
}

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// any0 <- function(x) isTRUE(any(x == 0)) ## ~= anyFalse
// anyFalse <- function(x) isTRUE(any(!x)) ## ~= any0
SEXP R_any0(SEXP x) {
    if (!isVectorAtomic(x)) {
	if(length(x) == 0) return FALSE_;
	// Typically S4.  TODO: Call the R code above, instead!
	error(_("Argument must be numeric-like atomic vector"));
    }
    R_xlen_t i, n = XLENGTH(x);
    if(n == 0) return FALSE_;

    switch(TYPEOF(x)) {
    case LGLSXP: {
	int *xx = LOGICAL(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    case INTSXP: {
	int *xx = INTEGER(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    case REALSXP: {
	double *xx = REAL(x);
	for(i=0; i < n; i++) if(xx[i] == 0.) return TRUE_;
	return FALSE_;
    }
    case RAWSXP: {
	unsigned char *xx = RAW(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    }
    error(_("Argument must be numeric-like atomic vector"));
    return R_NilValue; // -Wall
}

#undef TRUE_
#undef FALSE_

/**
 * A safe NEW_OBJECT(MAKE_CLASS(what)), where the caller must protect the
 * return value of this function
 *
 * @param an R character string specifying the name of a known S4 class
 */
SEXP NEW_OBJECT_OF_CLASS(const char* what)
{
    SEXP class = PROTECT(MAKE_CLASS(what)), obj = NEW_OBJECT(class);
    UNPROTECT(1);
    return obj;
}
