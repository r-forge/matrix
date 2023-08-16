#include "packedMatrix.h"

#define UNPACK(_PREFIX_, _CTYPE_, _ONE_) \
void _PREFIX_ ## dense_unpack(_CTYPE_ *dest, const _CTYPE_ *src, int n, \
                              char uplo, char diag) \
{ \
	int i, j; \
	R_xlen_t dpos = 0, spos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < n; dpos += n-(++j)) \
			for (i = 0; i <= j; ++i) \
				dest[dpos++] = src[spos++]; \
	} else { \
		for (j = 0; j < n; dpos += (++j)) \
			for (i = j; i <  n; ++i) \
				dest[dpos++] = src[spos++]; \
	} \
	if (diag != 'N') { \
		dpos = 0; \
		R_xlen_t n1a = (R_xlen_t) n + 1; \
		for (j = 0; j < n; ++j, dpos += n1a) \
			dest[dpos] = _ONE_; \
	} \
	return; \
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

#define MAKE_BANDED(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## dense_packed_make_banded(_CTYPE_ *x, int n, \
                                          int a, int b, \
                                          char uplo, char diag) \
{ \
	if (n == 0) \
		return; \
	if (a > b || a >= n || b <= -n) { \
		Matrix_memset(x, 0, PM_LENGTH(n), sizeof(_CTYPE_)); \
		return; \
	} \
	if (uplo == 'U') { \
		if (a <   0) a = 0; \
		if (b >=  n) b = n-1; \
	} else { \
		if (b >   0) b = 0; \
		if (a <= -n) a = 1-n; \
	} \
 \
	int i, j, i0, i1, \
		j0 = (a < 0) ? 0 : a, \
		j1 = (b < 0) ? n+b : n; \
 \
	if (uplo == 'U') { \
		if (j0 > 0) { \
			R_xlen_t dx; \
			Matrix_memset(x, 0, dx = PM_LENGTH(j0), \
			              sizeof(_CTYPE_)); \
			x += dx; \
		} \
		for (j = j0; j < j1; x += (++j)) { \
			i0 = j - b; \
			i1 = j - a + 1; \
			for (i = 0; i < i0; ++i) \
				*(x + i) = _ZERO_; \
			for (i = i1; i <= j; ++i) \
				*(x + i) = _ZERO_; \
		} \
		if (j1 < n) \
			Matrix_memset(x, 0, PM_LENGTH(n) - PM_LENGTH(j1), \
			              sizeof(_CTYPE_)); \
		if (diag != 'N' && a == 0) { \
			x -= PM_LENGTH(j); \
			for (j = 0; j < n; x += (++j)+1) \
				*x = _ONE_; \
		} \
	} else { \
		if (j0 > 0) { \
			R_xlen_t dx; \
			Matrix_memset(x, 0, dx = PM_LENGTH(n) - PM_LENGTH(j0), \
			              sizeof(_CTYPE_)); \
			x += dx; \
		} \
		for (j = j0; j < j1; x += n-(j++)) { \
			i0 = j - b; \
			i1 = j - a + 1; \
			for (i = j; i < i0; ++i) \
				*(x + i - j) = _ZERO_; \
			for (i = i1; i < n; ++i) \
				*(x + i - j) = _ZERO_; \
		} \
		if (j1 < n) \
			Matrix_memset(x, 0, PM_LENGTH(n - j1), \
			              sizeof(_CTYPE_)); \
		if (diag != 'N' && b == 0) { \
			x -= PM_LENGTH(n) - PM_LENGTH(j); \
			for (j = 0; j < n; x += n-(j++)) \
				*x = _ONE_; \
		} \
	} \
	return; \
}
MAKE_BANDED(d, double, 0.0, 1.0)
MAKE_BANDED(i, int, 0, 1)
MAKE_BANDED(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef MAKE_BANDED

#define COPY_DIAGONAL(_PREFIX_, _CTYPE_, _ONE_) \
void _PREFIX_ ## dense_packed_copy_diagonal(_CTYPE_ *dest, \
                                            const _CTYPE_ *src, \
                                            int n, R_xlen_t len, \
                                            char uplo_dest, \
                                            char uplo_src, \
                                            char diag) \
{ \
	int j; \
	if (diag != 'N') { \
		if (uplo_dest != 'L') { \
			for (j = 0; j < n; dest += (++j)+1) \
				*dest = _ONE_; \
		} else { \
			for (j = 0; j < n; dest += n-(j++)) \
				*dest = _ONE_; \
		} \
	} else if (len == n) { \
		/* copying from diagonalMatrix */ \
		if (uplo_dest != 'L') { \
			for (j = 0; j < n; dest += (++j)+1, ++src) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n-(j++), ++src) \
				*dest = *src; \
		} \
	} else if (len == PM_LENGTH(n)) { \
		/* copying from packedMatrix */ \
		if (uplo_dest != 'L') { \
			if (uplo_src != 'L') { \
				for (j = 0; j < n; src += (++j)+1, dest += j+1) \
					*dest = *src; \
			} else { \
				for (j = 0; j < n; src += n-j, dest += (++j)+1) \
					*dest = *src; \
			} \
		} else { \
			if (uplo_src != 'L') { \
				for (j = 0; j < n; dest += n-(j++), src += j+1) \
					*dest = *src; \
			} else { \
				for (j = 0; j < n; dest += n-j, src += n-(j++)) \
					*dest = *src; \
			} \
		} \
	} else if (len == (R_xlen_t) n * n) { \
		/* copying from square unpackedMatrix */ \
		R_xlen_t n1a = (R_xlen_t) n + 1; \
		if (uplo_dest != 'L') { \
			for (j = 0; j < n; dest += (++j)+1, src += n1a) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n-(j++), src += n1a) \
				*dest = *src; \
		} \
	} else { \
		error(_("incompatible '%s' and '%s' in %s()"), "n", "len", __func__); \
	} \
	return; \
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

#define TRANSPOSE(_PREFIX_, _CTYPE_) \
static void _PREFIX_ ## dense_packed_transpose(_CTYPE_ *dest, \
                                               const _CTYPE_ *src, \
                                               int n, char uplo) \
{ \
	int i, j; \
	if (uplo == 'U') { \
		for (j = 0; j < n; ++j) \
			for (i = j; i < n; ++i) \
				*(dest++) = *(src + PM_AR21_UP(j, i)); \
	} else { \
		R_xlen_t n2 = (R_xlen_t) n * 2; \
		for (j = 0; j < n; ++j) \
			for (i = 0; i <= j; ++i) \
				*(dest++) = *(src + PM_AR21_LO(j, i, n2)); \
	} \
	return; \
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
		ERROR_INVALID_TYPE(x, __func__);
	R_xlen_t nx = XLENGTH(x);
	SEXP y = PROTECT(allocVector(tx, nx));

#define TRANSPOSE(_PREFIX_, _PTR_) \
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
