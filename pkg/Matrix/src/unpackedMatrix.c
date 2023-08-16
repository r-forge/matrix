#include "unpackedMatrix.h"

#define PACK(_PREFIX_, _CTYPE_, _ONE_)								  \
void _PREFIX_ ## dense_pack(_CTYPE_ *dest, const _CTYPE_ *src, int n, \
                            char uplo, char diag) \
{ \
	int i, j; \
	R_xlen_t dpos = 0, spos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < n; spos += n-(++j)) \
			for (i = 0; i <= j; ++i) \
				dest[dpos++] = src[spos++]; \
		if (diag != 'N') { \
			dpos = 0; \
			for (j = 0; j < n; dpos += (++j)+1) \
				dest[dpos] = _ONE_; \
		} \
	} else { \
		for (j = 0; j < n; spos += (++j)) \
			for (i = j; i < n; ++i) \
				dest[dpos++] = src[spos++]; \
		if (diag != 'N') { \
			dpos = 0; \
			for (j = 0; j < n; dpos += n-(j++)) \
				dest[dpos] = _ONE_; \
		} \
	} \
	return; \
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

#define MAKE_TRIANGULAR(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## dense_unpacked_make_triangular(_CTYPE_ *x, int m, int n, \
                                                char uplo, char diag) \
{ \
	int i, j, r = (m < n) ? m : n; \
	R_xlen_t pos = 0; \
	if (uplo == 'U') { \
		for (j = 0; j < r; pos += (++j)+1) \
			for (i = j+1; i < m; ++i) \
				x[++pos] = _ZERO_; \
	} else { \
		for (j = 0; j < r; pos += m-(j++)) \
			for (i = 0; i < j; ++i) \
				x[pos++] = _ZERO_; \
		for (j = r; j < n; ++j) \
			for (i = 0; i < m; ++i) \
				x[pos++] = _ZERO_; \
	} \
	if (diag != 'N') { \
		pos = 0; \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		for (j = 0; j < r; ++j, pos += m1a) \
			x[pos] = _ONE_; \
	} \
	return; \
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
MAKE_TRIANGULAR(d, double, 0.0, 1.0)
/* idense_unpacked_make_triangular() */
MAKE_TRIANGULAR(i, int, 0, 1)
/* zdense_unpacked_make_triangular() */
MAKE_TRIANGULAR(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef MAKE_TRIANGULAR

#define MAKE_SYMMETRIC_SET_EQ(_X_, _DEST_, _SRC_) \
	_X_[_DEST_] = _X_[_SRC_]

#define MAKE_SYMMETRIC_SET_CJ(_X_, _DEST_, _SRC_) \
	do { \
		_X_[_DEST_].r =  _X_[_SRC_].r; \
		_X_[_DEST_].i = -_X_[_SRC_].i; \
	} while (0)

#define MAKE_SYMMETRIC(_PREFIX_, _CTYPE_, _SET_) \
void _PREFIX_ ## dense_unpacked_make_symmetric(_CTYPE_ *x, int n, \
                                               char uplo) \
{ \
	int i, j, n1s = n - 1; \
	R_xlen_t upos = n, lpos = 1; \
	if (uplo == 'U') { \
		for (j = 0; j < n; upos = (lpos += (++j)+1) + n1s) \
			for (i = j+1; i < n; ++i, upos += n, ++lpos) \
				_SET_(x, lpos, upos); \
	} else { \
		for (j = 0; j < n; upos = (lpos += (++j)+1) + n1s) \
			for (i = j+1; i < n; ++i, upos += n, ++lpos) \
				_SET_(x, upos, lpos); \
	} \
	return; \
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
MAKE_SYMMETRIC(d, double, MAKE_SYMMETRIC_SET_EQ)
/* idense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(i, int, MAKE_SYMMETRIC_SET_EQ)
/* zdense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(z, Rcomplex, MAKE_SYMMETRIC_SET_CJ)

#undef MAKE_SYMMETRIC
#undef MAKE_SYMMETRIC_SET_CJ
#undef MAKE_SYMMETRIC_SET_EQ

#define MAKE_BANDED(_PREFIX_, _CTYPE_, _ZERO_, _ONE_) \
void _PREFIX_ ## dense_unpacked_make_banded(_CTYPE_ *x, int m, int n, \
                                            int a, int b, char diag) \
{ \
	if (m == 0 || n == 0) \
		return; \
	if (a > b || a >= n || b <= -m) { \
		Matrix_memset(x, 0, (R_xlen_t) m * n, sizeof(_CTYPE_)); \
		return; \
	} \
	if (a <= -m) a = 1-m; \
	if (b >=  n) b = n-1; \
 \
	int i, j, i0, i1, \
		j0 = (a < 0) ? 0 : a, \
		j1 = (b < n-m) ? m+b : n; \
 \
	if (j0 > 0) { \
		R_xlen_t dx; \
		Matrix_memset(x, 0, dx = (R_xlen_t) m * j0, sizeof(_CTYPE_)); \
		x += dx; \
	} \
	for (j = j0; j < j1; ++j, x += m) { \
		i0 = j - b; \
		i1 = j - a + 1; \
		for (i = 0; i < i0; ++i) \
			*(x + i) = _ZERO_; \
		for (i = i1; i < m; ++i) \
			*(x + i) = _ZERO_; \
	} \
	if (j1 < n) \
		Matrix_memset(x, 0, (R_xlen_t) m * (n - j1), sizeof(_CTYPE_)); \
	if (diag != 'N' && a <= 0 && b >= 0) { \
		x -= m * (R_xlen_t) j; \
		R_xlen_t m1a = (R_xlen_t) m + 1; \
		for (j = 0; j < n; ++j, x += m1a) \
			*x = _ONE_; \
	} \
	return; \
}
MAKE_BANDED(d, double, 0.0, 1.0)
MAKE_BANDED(i, int, 0, 1)
MAKE_BANDED(z, Rcomplex, Matrix_zzero, Matrix_zone)

#undef MAKE_BANDED

#define COPY_DIAGONAL(_PREFIX_, _CTYPE_, _ONE_) \
void _PREFIX_ ## dense_unpacked_copy_diagonal(_CTYPE_ *dest, \
                                              const _CTYPE_ *src, \
                                              int n, R_xlen_t len, \
                                              char uplo, char diag) \
{ \
	int j; \
	R_xlen_t n1a = (R_xlen_t) n + 1; \
	if (diag != 'N') { \
		for (j = 0; j < n; ++j, dest += n1a) \
			*dest = _ONE_; \
	} else if (len == n) { \
		/* copying from diagonalMatrix */ \
		for (j = 0; j < n; ++j, dest += n1a, ++src) \
			*dest = *src; \
	} else if (len == (n * n1a) / 2) { \
		/* copying from packedMatrix */ \
		if (uplo == 'U') { \
			for (j = 0; j < n; dest += n1a, src += (++j)+1) \
				*dest = *src; \
		} else { \
			for (j = 0; j < n; dest += n1a, src += n-(j++)) \
				*dest = *src; \
		} \
	} else if (len == (R_xlen_t) n * n) { \
		/* copying from square unpackedMatrix */ \
		for (j = 0; j < n; ++j, dest += n1a, src += n1a) \
			*dest = *src; \
	} else { \
		error(_("incompatible '%s' and '%s' in %s()"), "n", "len", __func__); \
	} \
	return; \
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
COPY_DIAGONAL(d, double, 1.0)
/* idense_unpacked_copy_diagonal() */
COPY_DIAGONAL(i, int, 1)
/* zdense_unpacked_copy_diagonal() */
COPY_DIAGONAL(z, Rcomplex, Matrix_zone)

#undef COPY_DIAGONAL

SEXP unpacked_force(SEXP x, int n, char uplo, char diag)
{
	SEXPTYPE tx = TYPEOF(x);
	if (tx < LGLSXP || tx > CPLXSXP)
		ERROR_INVALID_TYPE(x, __func__);
	R_xlen_t nx = XLENGTH(x);
	SEXP y = PROTECT(allocVector(tx, nx));

	if (diag == '\0') {

#define FORCE_SYMMETRIC(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *py = _PTR_(y); \
		Matrix_memcpy(py, px, nx, sizeof(_CTYPE_)); \
		_PREFIX_ ## dense_unpacked_make_symmetric(py, n, uplo); \
	} while (0)

	switch (tx) {
	case LGLSXP:
		FORCE_SYMMETRIC(i, int, LOGICAL);
		break;
	case INTSXP:
		FORCE_SYMMETRIC(i, int, INTEGER);
		break;
	case REALSXP:
		FORCE_SYMMETRIC(d, double, REAL);
		break;
	case CPLXSXP:
		FORCE_SYMMETRIC(z, Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

#undef FORCE_SYMMETRIC

	} else {

#define FORCE_TRIANGULAR(_PREFIX_, _CTYPE_, _PTR_, _ONE_) \
		do { \
			_CTYPE_ *px = _PTR_(x), *py = _PTR_(y); \
			Matrix_memcpy(py, px, nx, sizeof(_CTYPE_)); \
			_PREFIX_ ## dense_unpacked_make_triangular( \
				py, n, n, uplo, diag); \
			if (diag != 'N') { \
				R_xlen_t n1a = (R_xlen_t) n + 1; \
				for (int j = 0; j < n; ++j, py += n1a) \
					*py = _ONE_; \
			} \
		} while (0)

		switch (tx) {
		case LGLSXP:
			FORCE_TRIANGULAR(i, int, LOGICAL, 1);
			break;
		case INTSXP:
			FORCE_TRIANGULAR(i, int, INTEGER, 1);
			break;
		case REALSXP:
			FORCE_TRIANGULAR(d, double, REAL, 1.0);
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
