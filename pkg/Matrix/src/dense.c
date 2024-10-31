#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "dense.h"

#define nCAST(x) (x != 0)
#define lCAST(x) (x)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void dense_colsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char nu, int narm, int mean,
                  SEXP ans)
{
	int i, j, count = -1, packed = class[2] == 'p';

#define SUM(c0, c1) \
	do { \
		c0##TYPE *px0 = c0##PTR(  x); \
		c1##TYPE *px1 = c1##PTR(ans); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##ZERO; \
				SUM_KERNEL(c0, c1, for (i = 0; i < m; ++i)); \
				px1 += 1; \
			} \
		} else if (nu == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##ZERO; \
				SUM_KERNEL(c0, c1, for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px0 += n - j - 1; \
				px1 += 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##ZERO; \
				if (!packed) \
					px0 += j; \
				SUM_KERNEL(c0, c1, for (i = j; i < n; ++i)); \
				px1 += 1; \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##UNIT; \
				SUM_KERNEL(c0, c1, for (i = 0; i < j; ++i)); \
				px0 += 1; \
				if (!packed) \
					px0 += n - j - 1; \
				px1 += 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##UNIT; \
				if (!packed) \
					px0 += j; \
				px0 += 1; \
				SUM_KERNEL(c0, c1, for (i = j + 1; i < n; ++i)); \
				px1 += 1; \
			} \
		} \
	} while (0)

#define SUM_KERNEL(c0, c1, __for__) \
	do { \
		if (mean) \
			count = m; \
		__for__ { \
			if (c0##NOT_NA(*px0)) \
				c1##INCREMENT_IDEN(*px1, c0##CAST(*px0)); \
			else if (!narm) \
				*px1 = c1##NA; \
			else if (mean) \
				--count; \
			px0 += 1; \
		} \
		if (mean) \
			c1##DIVIDE(*px1, count); \
	} while (0)

	switch (class[0]) {
	case 'n': if (mean) SUM(n, d); else SUM(n, i); break;
	case 'l': if (mean) SUM(l, d); else SUM(l, i); break;
	case 'i':                           SUM(i, d); break;
	case 'd':                           SUM(d, d); break;
	case 'z':                           SUM(z, z); break;
	default:                                       break;
	}

#undef SUM
#undef SUM_KERNEL

	return;
}

static
void dense_rowsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char nu, int narm, int mean,
                  SEXP ans)
{
	int i, j, *count = NULL, packed = XLENGTH(x) != (int_fast64_t) m * n,
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';

	if (mean && narm) {
		Matrix_Calloc(count, m, int);
		for (i = 0; i < m; ++i)
			count[i] = n;
	}

#define SUM(c0, c1) \
	do { \
		c0##TYPE *px0 = c0##PTR(  x), tmp0; \
		c1##TYPE *px1 = c1##PTR(ans), tmp1 = (un) ? c0##UNIT : c0##ZERO; \
		for (i = 0; i < m; ++i) \
			px1[i] = tmp1; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(c0, c1, for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || nu == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(c0, c1, for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px0 += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px0 += j; \
				SUM_KERNEL(c0, c1, for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(c0, c1, for (i = 0; i < j; ++i)); \
				px0 += 1; \
				if (!packed) \
					px0 += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px0 += j; \
				px0 += 1; \
				SUM_KERNEL(c0, c1, for (i = j + 1; i < n; ++i)); \
			} \
		} \
		if (mean) { \
			if (!narm) \
				for (i = 0; i < m; ++i) \
					c1##DIVIDE(px1[i], n); \
			else \
				for (i = 0; i < m; ++i) \
					c1##DIVIDE(px1[i], count[i]); \
		} \
	} while (0)

#define SUM_KERNEL(c0, c1, __for__) \
	do { \
		__for__ { \
			if (he && i == j) \
			c0##ASSIGN_PROJ_REAL(tmp0, *px0); \
			else \
			c0##ASSIGN_IDEN     (tmp0, *px0); \
			if (c0##NOT_NA(tmp0)) { \
				c1##INCREMENT_IDEN(px1[i], c0##CAST(tmp0)); \
				if (sy && i != j) { \
				if (he) \
				c1##INCREMENT_CONJ(px1[j], c0##CAST(tmp0)); \
				else \
				c1##INCREMENT_IDEN(px1[j], c0##CAST(tmp0)); \
				} \
			} else if (!narm) { \
				px1[i] = c1##NA; \
				if (sy && i != j) \
				px1[j] = c1##NA; \
			} else if (mean) { \
				--count[i]; \
				if (sy && i != j) \
				--count[j]; \
			} \
			px0 += 1; \
		} \
	} while (0)

	switch (class[0]) {
	case 'n': if (mean) SUM(n, d); else SUM(n, i); break;
	case 'l': if (mean) SUM(l, d); else SUM(l, i); break;
	case 'i':                           SUM(i, d); break;
	case 'd':                           SUM(d, d); break;
	case 'z':                           SUM(z, z); break;
	default:                                       break;
	}

#undef SUM
#undef SUM_KERNEL

	if (mean && narm)
		Matrix_Free(count, m);
	return;
}

SEXP dense_marginsum(SEXP obj, const char *class, int mg, int narm, int mean)
{
	narm = narm && class[0] != 'n';

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		r = (mg == 0) ? m : n;

	SEXP ans = PROTECT(Rf_allocVector(SUM_TYPEOF(class[0]), r)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));

	SEXP dimnames = DIMNAMES(obj, -(class[1] == 's')),
		marnames = VECTOR_ELT(dimnames, mg);
	if (marnames != R_NilValue) {
		PROTECT(marnames);
		Rf_setAttrib(ans, R_NamesSymbol, marnames);
		UNPROTECT(1); /* marnames */
	}

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	if (mg == 0 || class[1] == 's')
		dense_rowsum(x, class, m, n, ul, ct, nu, narm, mean, ans);
	else
		dense_colsum(x, class, m, n, ul, ct, nu, narm, mean, ans);

	UNPROTECT(2); /* x, ans */
	return ans;
}

/* (row|col)(Sums|Means)(<denseMatrix>, na.rm=) */
SEXP R_dense_marginsum(SEXP s_obj, SEXP s_margin, SEXP s_narm, SEXP s_mean)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int mg;
	VALID_MARGIN(s_margin, mg);

	int narm, mean;
	VALID_LOGIC2(s_narm, narm);
	VALID_LOGIC2(s_mean, mean);

	return dense_marginsum(s_obj, class, mg, narm, mean);
}

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

#define TRY_INCREMENT(x, y) \
		do { \
			if ((x >= 0) \
				? ( y <= INT_FAST64_MAX - x) \
				: (-y <= x - INT_FAST64_MIN)) { \
				x += y; \
				y = 0; \
				count = 0; \
			} else \
				goto over; \
		} while (0)

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) { \
					unsigned int d = (sy && i != j) ? 2 : 1; \
					if (count > UINT_MAX - d) \
						TRY_INCREMENT(s, t); \
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

		TRY_INCREMENT(s, t);

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

/* sum(<denseMatrix>, na.rm=) */
SEXP R_dense_sum(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return dense_sum(s_obj, class, narm);
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

/* prod(<denseMatrix>, na.rm=) */
SEXP R_dense_prod(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return dense_prod(s_obj, class, narm);
}
