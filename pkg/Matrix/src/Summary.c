/* C implementation of methods for sum, prod */

#include "Mdefines.h"
#include "M5.h"

#define LONGDOUBLE_AS_DOUBLE(x) \
(((x) > DBL_MAX) ? R_PosInf : (((x) < -DBL_MAX) ? R_NegInf : (double) (x)))

#define TRY_INCREMENT(x, y, label) \
do { \
	if ((x >= 0) \
		? ( y <= INT_FAST64_MAX - x) \
		: (-y <= x - INT_FAST64_MIN)) { \
		x += y; \
		y = 0; \
		count = 0; \
	} else \
		goto label; \
} while (0)

SEXP dense_sum(SEXP obj, const char *class, int narm)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym),
		ans = R_NilValue;
	int i, j, packed = class[2] == 'p',
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';

#define SUM \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || nu == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				SUM_KERNEL(for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(for (i = 0; i < j; ++i)); \
				px += 1; \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				px += 1; \
				SUM_KERNEL(for (i = j + 1; i < n; ++i)); \
			} \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	{
		int *px = LOGICAL(x);
		int_fast64_t s = (un) ? n : 0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != 0) \
					s += (sy && i != j) ? 2 : 1; \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		if (s <= INT_MAX)
			ans = Rf_ScalarInteger((int) s);
		else
			ans = Rf_ScalarReal((double) s);
		break;
	}
	case 'l':
	case 'i':
	{
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		int_fast64_t s = (un) ? n : 0, t = 0;
		unsigned int count = 0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) { \
					unsigned int d = (sy && i != j) ? 2 : 1; \
					if (count > UINT_MAX - d) \
						TRY_INCREMENT(s, t, over); \
					t += (sy && i != j) ? 2LL * *px : *px; \
					count += d; \
				} \
				else if (!narm) \
					return Rf_ScalarInteger(NA_INTEGER); \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		TRY_INCREMENT(s, t, over);

		if (s > INT_MIN && s <= INT_MAX)
			ans = Rf_ScalarInteger((int) s);
		else
			ans = Rf_ScalarReal((double) s);
		break;
over:
		px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		long double lr = (un) ? (long double) n : 0.0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) \
					lr += (sy && i != j) ? 2.0L * *px : *px; \
				else if (!narm) \
					return Rf_ScalarInteger(NA_INTEGER); \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
		break;
	}
	case 'd':
	{
		double *px = REAL(x);
		long double lr = (un) ? (long double) n : 0.0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && ISNAN(*px))) \
					lr += (sy && i != j) ? 2.0L * *px : *px; \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
		break;
	}
	case 'z':
	{
		Rcomplex *px = COMPLEX(x), tmp;
		long double lr = (un) ? (long double) n : 0.0;
		long double li = 0.0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && (ISNAN((*px).r) || (!he && ISNAN((*px).i))))) { \
					lr += (sy && i != j) ? 2.0L * (*px).r : (*px).r; \
					if (!he) \
					li += (sy && i != j) ? 2.0L * (*px).i : (*px).i; \
				} \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = Rf_ScalarComplex(tmp);
		break;
	}
	default:
		break;
	}

#undef SUM

	return ans;
}

SEXP sparse_sum(SEXP obj, const char *class, int narm)
{
	PROTECT(obj = sparse_aggregate(obj, class));

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = PROTECT((class[0] == 'n') ? R_NilValue : GET_SLOT(obj, Matrix_xSym)),
		ans = R_NilValue;
	int sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';

	if (class[2] != 'T') {

		if (class[2] == 'R')
			SWAP(m, n, int, );

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		pp++;

		UNPROTECT(4); /* i, p, x, obj */

		switch (class[0]) {
		case 'n':
		{
			int_fast64_t s = (int_fast64_t) pp[n - 1] + ((un) ? n : 0);
			if (sy) {
				int up = (class[2] == 'C') == (ul == 'U');
				s *= 2;
				for (j = 0, k = 0; j < n; ++j) {
					kend = pp[j];
					if (k < kend && pi[(up) ? kend - 1 : k] == j)
						--s;
					k = kend;
				}
			}
			if (s <= INT_MAX)
				ans = Rf_ScalarInteger((int) s);
			else
				ans = Rf_ScalarReal((double) s);
			break;
		}
		case 'l':
		case 'i':
		{
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			int_fast64_t s = (un) ? n : 0, t = 0;
			unsigned int count = 0;
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				while (k < kend) {
					if (px[k] != NA_INTEGER) {
						unsigned int d = (sy && pi[k] != j) ? 2 : 1;
						if (count > UINT_MAX - d)
							TRY_INCREMENT(s, t, overC);
						t += (sy && pi[k] != j) ? 2LL * px[k] : px[k];
						count += d;
					}
					else if (!narm)
						return Rf_ScalarInteger(NA_INTEGER);
					++k;
				}
			}
			TRY_INCREMENT(s, t, overC);
			if (s > INT_MIN && s <= INT_MAX)
				ans = Rf_ScalarInteger((int) s);
			else
				ans = Rf_ScalarReal((double) s);
			break;
overC:
			;
			long double lr = (long double) s + (long double) t;
			for (; j < n; ++j) {
				kend = pp[j];
				while (k < kend) {
					if (px[k] != NA_INTEGER)
						lr += (sy && pi[k] != j) ? 2.0L * px[k] : px[k];
					else if (!narm)
						return Rf_ScalarInteger(NA_INTEGER);
					++k;
				}
			}
			ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
			break;
		}
		case 'd':
		{
			double *px = REAL(x);
			long double lr = (un) ? (long double) n : 0.0;
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				while (k < kend) {
					if (!(narm && ISNAN(px[k])))
						lr += (sy && pi[k] != j) ? 2.0L * px[k] : px[k];
					++k;
				}
			}
			ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x), tmp;
			long double lr = (un) ? (long double) n : 0.0;
			long double li = 0.0;
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				while (k < kend) {
					if (!(narm && (ISNAN(px[k].r) || (!he && ISNAN(px[k].i))))) {
						lr += (sy && pi[k] != j) ? 2.0L * px[k].r : px[k].r;
						if (!he)
						li += (sy && pi[k] != j) ? 2.0L * px[k].i : px[k].i;
					}
					++k;
				}
			}
			tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
			tmp.i = LONGDOUBLE_AS_DOUBLE(li);
			ans = Rf_ScalarComplex(tmp);
			break;
		}
		default:
			break;
		}

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

		UNPROTECT(4); /* j, i, x, obj */

		switch (class[0]) {
		case 'n':
		{
			int_fast64_t s = (int_fast64_t) XLENGTH(i) + ((un) ? n : 0);
			if (sy) {
				s *= 2;
				for (k = 0; k < kend; ++k)
					if (pi[k] == pj[k])
						--s;
			}
			if (s <= INT_MAX)
				ans = Rf_ScalarInteger((int) s);
			else
				ans = Rf_ScalarReal((double) s);
			break;
		}
		case 'l':
		case 'i':
		{
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			int_fast64_t s = (un) ? n : 0, t = 0;
			unsigned int count = 0;
			for (k = 0; k < kend; ++k)
				if (px[k] != NA_INTEGER) {
					unsigned int d = (sy && pi[k] != pj[k]) ? 2 : 1;
					if (count > UINT_MAX - d)
						TRY_INCREMENT(s, t, overT);
					t += (sy && pi[k] != pj[k]) ? 2LL * px[k] : px[k];
					count += d;
				}
				else if (!narm)
					return Rf_ScalarInteger(NA_INTEGER);
			TRY_INCREMENT(s, t, overT);
			if (s > INT_MIN && s <= INT_MAX)
				ans = Rf_ScalarInteger((int) s);
			else
				ans = Rf_ScalarReal((double) s);
			break;
overT:
			;
			long double lr = (long double) s + (long double) t;
			for (; k < kend; ++k) {
				if (px[k] != NA_INTEGER)
					lr += (sy && pi[k] != pj[k]) ? 2.0L * px[k] : px[k];
				else if (!narm)
					return Rf_ScalarInteger(NA_INTEGER);
			}
			ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
			break;
		}
		case 'd':
		{
			double *px = REAL(x);
			long double lr = (un) ? (long double) n : 0.0;
			for (k = 0; k < kend; ++k)
				if (!(narm && ISNAN(px[k])))
					lr += (sy && pi[k] != pj[k]) ? 2.0L * px[k] : px[k];
			ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x), tmp;
			long double lr = (un) ? (long double) n : 0.0;
			long double li = 0.0;
			for (k = 0; k < kend; ++k)
				if (!(narm && (ISNAN(px[k].r) || (!he && ISNAN(px[k].i))))) {
					lr += (sy && pi[k] != pj[k]) ? 2.0L * px[k].r : px[k].r;
					if (!he)
					li += (sy && pi[k] != pj[k]) ? 2.0L * px[k].i : px[k].i;
				}
			tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
			tmp.i = LONGDOUBLE_AS_DOUBLE(li);
			ans = Rf_ScalarComplex(tmp);
			break;
		}
		default:
			break;
		}

	}

	return ans;
}

SEXP R_dense_sum(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return dense_sum(s_obj, class, narm);
}

SEXP R_sparse_sum(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return sparse_sum(s_obj, class, narm);
}

SEXP dense_prod(SEXP obj, const char *class, int narm)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p',
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';
	long double lr = 1.0, li = 0.0;

#define PROD \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				PROD_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				PROD_KERNEL(for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				PROD_KERNEL(for (i = j; i < n; ++i)); \
			} \
		} else if (nu == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				PROD_KERNEL(for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				if (!packed) \
					px += j; \
				PROD_KERNEL(for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				PROD_KERNEL(for (i = 0; i < j; ++i)); \
				px += 1; \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				if (!packed) \
					px += j; \
				px += 1; \
				PROD_KERNEL(for (i = j + 1; i < n; ++i)); \
			} \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	{
		int *px = LOGICAL(x);
		if (class[1] == 't') {
			if (n > 1 || (n == 1 && !un && *px == 0))
				return Rf_ScalarReal(0.0);
			break;
		}

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px == 0) \
					return Rf_ScalarReal(0.0); \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	case 'l':
	case 'i':
	{
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) \
					lr *= (sy && i != j) ? (long double) *px * *px : *px; \
				else if (!narm) \
					lr *= NA_REAL; \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	case 'd':
	{
		double *px = REAL(x);

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && ISNAN(*px))) \
					lr *= (sy && i != j) ? (long double) *px * *px : *px; \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	case 'z':
	{
		Rcomplex *px = COMPLEX(x);
		long double lr0, li0;

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && (ISNAN((*px).r) || (!(he && i == j) && ISNAN((*px).i))))) { \
					if (he) { \
						lr0 = (*px).r; \
						if (i != j) { \
						li0 = (*px).i; \
						lr0 = lr0 * lr0 + li0 * li0; \
						} \
						lr *= lr0; \
						li *= lr0; \
					} else { \
						lr0 = lr; li0 = li; \
						lr = lr0 * (*px).r - li0 * (*px).i; \
						li = li0 * (*px).r + lr0 * (*px).i; \
						if (i != j) { \
						lr0 = lr; li0 = li; \
						lr = lr0 * (*px).r - li0 * (*px).i; \
						li = li0 * (*px).r + lr0 * (*px).i; \
						} \
					} \
				} \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	default:
		break;
	}

#undef PROD

	SEXP ans;
	if (class[0] != 'z')
		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
	else {
		Rcomplex tmp;
		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = Rf_ScalarComplex(tmp);
	}
	return ans;
}

SEXP sparse_prod(SEXP obj, const char *class, int narm)
{
	PROTECT(obj = sparse_aggregate(obj, class));

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = PROTECT((class[0] == 'n') ? R_NilValue : GET_SLOT(obj, Matrix_xSym));
	int sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';
	long double lr = 1.0, li = 0.0;

	int_fast64_t mn = (int_fast64_t) m * n,
		q = (sy) ? (mn + n) / 2 : ((un) ? mn - n : mn);

	if (class[2] != 'T') {

		if (class[2] == 'R')
			SWAP(m, n, int, );

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		pp++;

		UNPROTECT(4); /* i, p, x, obj */

		int seen0 = 0, i0, i1,
			up = !sy || (class[2] == 'C') == (ul == 'U'),
			lo = !sy || (class[2] == 'C') == (ul != 'U');

		switch (class[0]) {
		case 'n':
			if (pp[n - 1] < q)
				lr = 0.0;
			break;
		case 'l':
		case 'i':
		{
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				if (!seen0) {
					i0 = (up) ? 0 : j;
					i1 = (lo) ? m : j + 1;
				}
				while (k < kend) {
					if (!seen0) {
						if (pi[k] == i0 || (un && i0 == j && pi[k] == ++i0))
							++i0;
						else {
							lr *= 0.0;
							seen0 = 1;
						}
					}
					if (px[k] != NA_INTEGER)
						lr *= (sy && pi[k] != j)
							? (long double) px[k] * px[k] : px[k];
					else if (!narm)
						lr *= NA_REAL;
					++k;
				}
				if (!(seen0 || i1 == i0 || (un && i0 == j && i1 == ++i0))) {
					lr *= 0.0;
					seen0 = 1;
				}
			}
			break;
		}
		case 'd':
		{
			double *px = REAL(x);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				if (!seen0) {
					i0 = (up) ? 0 : j;
					i1 = (lo) ? m : j + 1;
				}
				while (k < kend) {
					if (!seen0) {
						if (pi[k] == i0 || (un && i0 == j && pi[k] == ++i0))
							++i0;
						else {
							lr *= 0.0;
							seen0 = 1;
						}
					}
					if (!(narm && ISNAN(px[k])))
						lr *= (sy && pi[k] != j)
							? (long double) px[k] * px[k] : px[k];
					++k;
				}
				if (!(seen0 || i1 == i0 || (un && i0 == j && i1 == ++i0))) {
					lr *= 0.0;
					seen0 = 1;
				}
			}
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x);
			long double lr0, li0;
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				if (!seen0) {
					i0 = (up) ? 0 : j;
					i1 = (lo) ? m : j + 1;
				}
				while (k < kend) {
					if (!seen0) {
						if (pi[k] == i0 || (un && i0 == j && pi[k] == ++i0))
							++i0;
						else {
							lr *= 0.0;
							li *= 0.0;
							seen0 = 1;
						}
					}
					if (!(narm && (ISNAN(px[k].r) || (!(he && pi[k] == j) && ISNAN(px[k].i))))) {
						if (he) {
							lr0 = px[k].r;
							if (pi[k] != j) {
							li0 = px[k].i;
							lr0 = lr0 * lr0 + li0 * li0;
							}
							lr *= lr0;
							li *= lr0;
						} else {
							lr0 = lr; li0 = li;
							lr = lr0 * px[k].r - li0 * px[k].i;
							li = li0 * px[k].r + lr0 * px[k].i;
							if (pi[k] != j) {
							lr0 = lr; li0 = li;
							lr = lr0 * px[k].r - li0 * px[k].i;
							li = li0 * px[k].r + lr0 * px[k].i;
							}
						}
					}
					++k;
				}
				if (!(seen0 || i1 == i0 || (un && i0 == j && i1 == ++i0))) {
					lr *= 0.0;
					li *= 0.0;
					seen0 = 1;
				}
			}
			break;
		}
		default:
			break;
		}

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);
		if (XLENGTH(i) < q)
			lr = 0.0;

		UNPROTECT(4); /* j, i, x, obj */

		switch (class[0]) {
		case 'l':
		case 'i':
		{
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			for (k = 0; k < kend; ++k)
				if (px[k] != NA_INTEGER)
					lr *= (sy && pi[k] != pj[k])
						? (long double) px[k] * px[k] : px[k];
				else if (!narm)
					lr *= NA_REAL;
			break;
		}
		case 'd':
		{
			double *px = REAL(x);
			for (k = 0; k < kend; ++k)
				if (!(narm && ISNAN(px[k])))
					lr *= (sy && pi[k] != pj[k])
						? (long double) px[k] * px[k] : px[k];
			break;
		}
		case 'z':
		{
			Rcomplex *px = COMPLEX(x);
			long double lr0, li0;
			for (k = 0; k < kend; ++k)
				if (!(narm && (ISNAN(px[k].r) || (!(he && pi[k] == pj[k]) && ISNAN(px[k].i))))) {
					if (he) {
						lr0 = px[k].r;
						if (pi[k] != pj[k]) {
						li0 = px[k].i;
						lr0 = lr0 * lr0 + li0 * li0;
						}
						lr *= lr0;
						li *= lr0;
					} else {
						lr0 = lr; li0 = li;
						lr = lr0 * px[k].r - li0 * px[k].i;
						li = li0 * px[k].r + lr0 * px[k].i;
						if (pi[k] != pj[k]) {
						lr0 = lr; li0 = li;
						lr = lr0 * px[k].r - li0 * px[k].i;
						li = li0 * px[k].r + lr0 * px[k].i;
						}
					}
				}
			break;
		}
		default:
			break;
		}

	}

	SEXP ans;
	if (class[0] != 'z')
		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
	else {
		Rcomplex tmp;
		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = Rf_ScalarComplex(tmp);
	}
	return ans;
}

SEXP R_dense_prod(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return dense_prod(s_obj, class, narm);
}

SEXP R_sparse_prod(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return sparse_prod(s_obj, class, narm);
}
