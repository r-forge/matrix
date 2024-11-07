/* C implementation of methods for colSums, colMeans, rowSums, rowMeans */

#include "Mdefines.h"
#include "M5.h"

/* defined in ./aggregate.c : */
SEXP sparse_aggregate(SEXP, const char *);

#define SUM_TYPEOF(c) \
(c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

#define MAP(i) (map) ? map[i] : i

#define lCAST(x) (x)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#undef nCAST
#define nCAST(x) (x != 0)

static
void dense_colsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char nu,
                  int narm, int mean,
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
                  int m, int n, char ul, char ct, char nu,
                  int narm, int mean,
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

#undef nCAST
#define nCAST(x) (1)

static
void Csparse_colsum(SEXP obj, const char *class,
                    int m, int n, char ul, char ct, char nu,
                    int narm, int mean,
                    SEXP ans)
{
	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp0 = INTEGER(p0), j, k, kend, nnz1, count = -1,
		un = class[1] == 't' && nu != 'N';
	pp0++;

	SEXP x1;
	if (TYPEOF(ans) != OBJSXP) {
		nnz1 = n;
		x1 = ans;
	} else {
		if (un)
			nnz1 = n;
		else {
			nnz1 = 0;
			for (j = 0; j < n; ++j)
				if (pp0[j - 1] < pp0[j])
					++nnz1;
		}

		SEXP j1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		int *pj1 = INTEGER(j1);
		SET_SLOT(ans, Matrix_iSym, j1);
		if (un)
			for (j = 0; j < n; ++j)
				*(pj1++) = j + 1;
		else
			for (j = 0; j < n; ++j)
				if (pp0[j - 1] < pp0[j])
					*(pj1++) = j + 1;

		PROTECT(x1 = Rf_allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(ans, Matrix_xSym, x1);

		UNPROTECT(2); /* x1, j1 */
	}
	PROTECT(x1);

	int full = nnz1 == n;

#define SUM(c0, c1) \
	do { \
		c0##IF_NPATTERN( \
		SEXP x0 = GET_SLOT(obj, Matrix_xSym); \
		c0##TYPE *px0 = c0##PTR(x0); \
		); \
		c1##TYPE *px1 = c1##PTR(x1), tmp1 = (un) ? c1##UNIT : c1##ZERO; \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			if (full || k < kend) { \
				*px1 = tmp1; \
				if (mean) \
					count = m; \
				while (k < kend) { \
					if (c0##NOT_NA(*px0)) \
						c1##INCREMENT_IDEN(*px1, c0##CAST(*px0)); \
					else if (!narm) \
						*px1 = c1##NA; \
					else if (mean) \
						--count; \
					c0##IF_NPATTERN( \
					++px0; \
					); \
					++k; \
				} \
				if (mean) \
					c1##DIVIDE(*px1, count); \
				++px1; \
			} \
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

	UNPROTECT(2); /* x1, p0 */
	return;
}

static
void Csparse_rowsum(SEXP obj, const char *class,
                    int m, int n, char ul, char ct, char nu,
                    int narm, int mean,
                    SEXP ans, SEXP iSym)
{
	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(obj,        iSym));
	int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), i, j, k, kend, nnz1,
		*count = NULL, *map = NULL,
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';
	pp0++;

	SEXP x1;
	if (TYPEOF(ans) != OBJSXP) {
		nnz1 = m;
		x1 = ans;
		if (narm && mean) {
			Matrix_Calloc(count, m, int);
			for (i = 0; i < m; ++i)
				count[i] = n;
		}
	} else {
		if (un)
			nnz1 = m;
		else {
			nnz1 = 0;
			Matrix_Calloc(map, m, int);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					i = pi0[k];
					++map[i];
					if (sy && i != j)
					++map[j];
					++k;
				}
			}
			for (i = 0; i < m; ++i)
				map[i] = (map[i]) ? nnz1++ : -1;
		}

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		SET_SLOT(ans, Matrix_iSym, i1);
		if (narm && mean) {
			count = INTEGER(i1);
			for (i = 0; i < nnz1; ++i)
				count[i] = n;
		}

		PROTECT(x1 = Rf_allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(ans, Matrix_xSym, x1);

		UNPROTECT(2); /* x1, i1 */
	}
	PROTECT(x1);

#define SUM(c0, c1) \
	do { \
		c0##IF_NPATTERN( \
		SEXP x0 = GET_SLOT(obj, Matrix_xSym); \
		c0##TYPE *px0 = c0##PTR(x0); \
		); \
		c0##TYPE                     tmp0; \
		c1##TYPE *px1 = c1##PTR(x1), tmp1 = (un) ? c1##UNIT : c1##ZERO; \
		for (i = 0; i < nnz1; ++i) \
			px1[i] = tmp1; \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			while (k < kend) { \
				i = pi0[k]; \
				if (he && i == j) \
				c0##ASSIGN_PROJ_REAL(tmp0, c0##IFELSE_NPATTERN(px0[k], c0##UNIT)); \
				else \
				c0##ASSIGN_IDEN     (tmp0, c0##IFELSE_NPATTERN(px0[k], c0##UNIT)); \
				if (c0##NOT_NA(tmp0)) { \
					c1##INCREMENT_IDEN(px1[MAP(i)], c0##CAST(tmp0)); \
					if (sy && i != j) { \
					if (he) \
					c1##INCREMENT_CONJ(px1[MAP(j)], c0##CAST(tmp0)); \
					else \
					c1##INCREMENT_IDEN(px1[MAP(j)], c0##CAST(tmp0)); \
					} \
				} \
				else if (!narm) { \
					px1[MAP(i)] = c1##NA; \
					if (sy && i != j) \
					px1[MAP(j)] = c1##NA; \
				} \
				else if (mean) { \
					--count[MAP(i)]; \
					if (sy && i != j) \
					--count[MAP(j)]; \
				} \
				++k; \
			} \
		} \
		if (mean) { \
			if (!narm) \
				for (i = 0; i < nnz1; ++i) \
					c1##DIVIDE(px1[i], n); \
			else \
				for (i = 0; i < nnz1; ++i) \
					c1##DIVIDE(px1[i], count[i]); \
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

	if (TYPEOF(ans) != OBJSXP) {
		if (count)
			Matrix_Free(count, m);
	} else {
		SEXP i1 = GET_SLOT(ans, Matrix_iSym);
		int *pi1 = INTEGER(i1);
		if (map) {
			for (i = 0; i < m; ++i)
				if (map[i] >= 0)
					pi1[map[i]] = i + 1;
			Matrix_Free(map, m);
		}
		else
			for (i = 0; i < m; ++i)
				pi1[i] = i + 1;
	}

	UNPROTECT(3); /* x1, i0, p0 */
	return;
}

static
void Tsparse_colsum(SEXP obj, const char *class,
                    int m, int n, char ul, char ct, char nu,
                    int narm, int mean,
                    SEXP ans, SEXP iSym, SEXP jSym)
{
	if (narm && mean)
		obj = sparse_aggregate(obj, class);
	PROTECT(obj);

	SEXP i0 = PROTECT(GET_SLOT(obj, iSym)),
		j0 = PROTECT(GET_SLOT(obj, jSym));
	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), i, j, nnz1,
		*count = NULL, *map = NULL,
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';;
	R_xlen_t k, kend = XLENGTH(i0);

	SEXP x1;
	if (TYPEOF(ans) != OBJSXP) {
		nnz1 = n;
		x1 = ans;
		if (narm && mean) {
			Matrix_Calloc(count, n, int);
			for (j = 0; j < n; ++j)
				count[j] = m;
		}
	} else {
		if (un)
			nnz1 = n;
		else {
			nnz1 = 0;
			Matrix_Calloc(map, n, int);
			for (k = 0; k < kend; ++k) {
				i = pi0[k];
				j = pj0[k];
				++map[j];
				if (sy && i != j)
				++map[i];
			}
			for (j = 0; j < n; ++j)
				map[j] = (map[j]) ? nnz1++ : -1;
		}

		SEXP j1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		SET_SLOT(ans, Matrix_iSym, j1);
		if (narm && mean) {
			count = INTEGER(j1);
			for (j = 0; j < nnz1; ++j)
				count[j] = m;
		}

		PROTECT(x1 = Rf_allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(ans, Matrix_xSym, x1);

		UNPROTECT(2); /* x1, j1 */
	}
	PROTECT(x1);

#define SUM(c0, c1) \
	do { \
		c0##IF_NPATTERN( \
		SEXP x0 = GET_SLOT(obj, Matrix_xSym); \
		c0##TYPE *px0 = c0##PTR(x0); \
		); \
		c0##TYPE                     tmp0; \
		c1##TYPE *px1 = c1##PTR(x1), tmp1 = (un) ? c1##UNIT : c1##ZERO; \
		for (j = 0; j < nnz1; ++j) \
			px1[j] = tmp1; \
		for (k = 0; k < kend; ++k) { \
			i = pi0[k]; \
			j = pj0[k]; \
			if (he && i == j) \
			c0##ASSIGN_PROJ_REAL(tmp0, c0##IFELSE_NPATTERN(px0[k], c0##UNIT)); \
			else \
			c0##ASSIGN_IDEN     (tmp0, c0##IFELSE_NPATTERN(px0[k], c0##UNIT)); \
			if (c0##NOT_NA(tmp0)) { \
				c1##INCREMENT_IDEN(px1[MAP(j)], c0##CAST(tmp0)); \
				if (sy && i != j) { \
				if (he) \
				c1##INCREMENT_CONJ(px1[MAP(i)], c0##CAST(tmp0)); \
				else \
				c1##INCREMENT_IDEN(px1[MAP(i)], c0##CAST(tmp0)); \
				} \
			} \
			else if (!narm) { \
				px1[MAP(j)] = c1##NA; \
				if (sy && i != j) \
				px1[MAP(i)] = c1##NA; \
			} \
			else if (mean) { \
				--count[MAP(j)]; \
				if (sy && i != j) \
				--count[MAP(i)]; \
			} \
		} \
		if (mean) { \
			if (!narm) \
				for (j = 0; j < nnz1; ++j) \
					c1##DIVIDE(px1[j], m); \
			else \
				for (j = 0; j < nnz1; ++j) \
					c1##DIVIDE(px1[j], count[j]); \
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

	if (TYPEOF(ans) != OBJSXP) {
		if (count)
			Matrix_Free(count, n);
	} else {
		SEXP j1 = GET_SLOT(ans, Matrix_iSym);
		int *pj1 = INTEGER(j1);
		if (map) {
			for (j = 0; j < n; ++j)
				if (map[j] >= 0)
					pj1[map[j]] = j + 1;
			Matrix_Free(map, n);
		}
		else
			for (j = 0; j < n; ++j)
				pj1[j] = j + 1;
	}

	UNPROTECT(4); /* x1, j0, i0, obj */
	return;
}

SEXP dense_marginsum(SEXP obj, const char *class, int mg,
                     int narm, int mean)
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

SEXP sparse_marginsum(SEXP obj, const char *class, int mg,
                      int narm, int mean, int sparse)
{
	narm = narm && class[0] != 'n';

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		r = (mg == 0) ? m : n;

	SEXP ans;
	SEXPTYPE type = SUM_TYPEOF(class[0]);
	if (sparse) {
		char cl[] = ".sparseVector";
		cl[0] = typeToKind(type);
		PROTECT(ans = newObject(cl));

		SEXP length = GET_SLOT(ans, Matrix_lengthSym);
		INTEGER(length)[0] = r;
	} else {
		PROTECT(ans = Rf_allocVector(type, r));

		SEXP dimnames = DIMNAMES(obj, -(class[1] == 's')),
			marnames = VECTOR_ELT(dimnames, mg);
		if (marnames != R_NilValue) {
			PROTECT(marnames);
			Rf_setAttrib(ans, R_NamesSymbol, marnames);
			UNPROTECT(1); /* marnames */
		}
	}

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	if (mg == 0)
		switch (class[2]) {
		case 'C':
			Csparse_rowsum(obj, class, m, n, ul, ct, nu, narm, mean,
			               ans, Matrix_iSym);
			break;
		case 'R':
			if (class[1] == 's')
			Csparse_rowsum(obj, class, n, m, ul, ct, nu, narm, mean,
			               ans, Matrix_jSym);
			else
			Csparse_colsum(obj, class, n, m, ul, ct, nu, narm, mean,
			               ans);
			break;
		case 'T':
			Tsparse_colsum(obj, class, n, m, ul, ct, nu, narm, mean,
			               ans, Matrix_jSym, Matrix_iSym);
			break;
		default:
			break;
		}
	else
		switch (class[2]) {
		case 'C':
			if (class[1] == 's')
			Csparse_rowsum(obj, class, m, n, ul, ct, nu, narm, mean,
			               ans, Matrix_iSym);
			else
			Csparse_colsum(obj, class, m, n, ul, ct, nu, narm, mean,
			               ans);
			break;
		case 'R':
			Csparse_rowsum(obj, class, n, m, ul, ct, nu, narm, mean,
			               ans, Matrix_jSym);
			break;
		case 'T':
			Tsparse_colsum(obj, class, m, n, ul, ct, nu, narm, mean,
			               ans, Matrix_iSym, Matrix_jSym);
			break;
		default:
			break;
		}

	UNPROTECT(1); /* ans */
	return ans;
}

SEXP R_sparse_marginsum(SEXP s_obj, SEXP s_margin,
                        SEXP s_narm, SEXP s_mean, SEXP s_sparse)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int mg;
	VALID_MARGIN(s_margin, mg);

	int narm, mean, sparse;
	VALID_LOGIC2(s_narm  , narm  );
	VALID_LOGIC2(s_mean  , mean  );
	VALID_LOGIC2(s_sparse, sparse);

	return sparse_marginsum(s_obj, class, mg, narm, mean, sparse);
}

SEXP R_dense_marginsum(SEXP s_obj, SEXP s_margin,
                       SEXP s_narm, SEXP s_mean)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int mg;
	VALID_MARGIN(s_margin, mg);

	int narm, mean;
	VALID_LOGIC2(s_narm, narm);
	VALID_LOGIC2(s_mean, mean);

	return dense_marginsum(s_obj, class, mg, narm, mean);
}
