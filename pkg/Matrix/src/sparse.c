#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "sparse.h"

SEXP sparse_aggregate(SEXP from, const char *class)
{
	if (class[2] != 'T')
		return from;

	SEXP i0 = GET_SLOT(from, Matrix_iSym);
	if (XLENGTH(i0) < 2)
		return from;
	if (XLENGTH(i0) > INT_MAX)
		Rf_error(_("number of triplets to be aggregated exceeds %s"),
		         "2^31-1");
	PROTECT(i0);
	SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym)), i1, j1, to;
	int *pi0 = INTEGER(i0), *pi1 = NULL,
		*pj0 = INTEGER(j0), *pj1 = NULL,
		*pdim = DIM(from), m = pdim[0], n = pdim[1],
		nnz = (int) XLENGTH(i0), nnz_ = nnz, *iwork = NULL;
	size_t liwork = (size_t) ((int_fast64_t) n + 1 + n + m + nnz),
		lwork = (size_t) nnz;
	Matrix_Calloc(iwork, liwork, int);

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = NULL, *px1 = NULL, *work = NULL; \
		c##IF_NPATTERN( \
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
		px0 = c##PTR(x0); \
		Matrix_Calloc(work, lwork, c##TYPE); \
		); \
		c##tspaggr(pj1, pi1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
		if (nnz != nnz_) { \
		PROTECT(to = newObject(class)); \
		PROTECT(i1 = Rf_allocVector(INTSXP, nnz)), \
		PROTECT(j1 = Rf_allocVector(INTSXP, nnz)); \
		pi1 = INTEGER(i1); \
		pj1 = INTEGER(j1); \
		SET_SLOT(to, Matrix_iSym, i1); \
		SET_SLOT(to, Matrix_jSym, j1); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
		px1 = c##PTR(x1); \
		SET_SLOT(to, Matrix_xSym, x1); \
		); \
		c##tspaggr(pj1, pi1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
		c##IF_NPATTERN( \
		UNPROTECT(1); /* x1 */ \
		); \
		UNPROTECT(3); /* j1, i1, to */ \
		} \
		c##IF_NPATTERN( \
		Matrix_Free(work, lwork); \
		UNPROTECT(1); /* x0 */ \
		); \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

	Matrix_Free(iwork, liwork);
	UNPROTECT(2); /* j0, i0 */

	if (nnz == nnz_)
		return from;

	PROTECT(to);
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);
	UNPROTECT(1); /* to */
	return to;
}

/* aggregateT(<[CRT]sparseMatrix>) */
SEXP R_sparse_aggregate(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_aggregate(s_from, class);
}

SEXP sparse_drop0(SEXP from, const char *class, double tol)
{
	from = sparse_aggregate(from, class);
	if (class[0] == 'n')
		return from;
	PROTECT(from);

	char ct = '\0';
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(from);

	SEXP to;
	int strict = ISNAN(tol) || tol <= 0.0;

#define NZ(c, x) ((strict) ? c##NOT_ZERO(x) : c##NOT_ZERO_TOL(x, tol))

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			n = (int) (XLENGTH(p0) - 1), j, k, kend,
			nnz0 = pp0[n], nnz1 = 0;
		pp0++;

		if (ct != 'C') {

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0); \
				for (k = 0; k < nnz0; ++k) \
					if (NZ(c, px0[k])) \
						++nnz1; \
			} while (0)

			SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		} else {

			Rcomplex *px0 = COMPLEX(x0);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((pi0[k] == j) ? NZ(d, px0[k].r) : NZ(z, px0[k]))
						++nnz1;
					++k;
				}
			}

		}

		if (nnz1 == nnz0) {
			UNPROTECT(4); /* x0, i0, p0, from */
			return from;
		}

		PROTECT(to = newObject(class));

		SEXP p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0))),
			i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
			x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz1));
		int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		SET_SLOT(to, Matrix_xSym, x1);

		if (ct != 'C') {

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (j = 0, k = 0; j < n; ++j) { \
					pp1[j] = pp1[j - 1]; \
					kend = pp0[j]; \
					while (k < kend) { \
						if (NZ(c, px0[k])) { \
							pi1[pp1[j]  ] = pi0[k]; \
							px1[pp1[j]++] = px0[k]; \
						} \
						++k; \
					} \
				} \
			} while (0)

			SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		} else {

			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
			for (j = 0, k = 0; j < n; ++j) {
				pp1[j] = pp1[j - 1];
				kend = pp0[j];
				while (k < kend) {
					if ((pi0[k] == j) ? NZ(d, px0[k].r) : NZ(z, px0[k])) {
						pi1[pp1[j]  ] = pi0[k];
						px1[pp1[j]++] = px0[k];
					}
					++k;
				}
			}

		}

		UNPROTECT(7); /* x1, i1, p1, to, x0, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(x0), nnz1 = 0;

		if (ct != 'C') {

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0); \
				for (k = 0; k < nnz0; ++k) \
					if (NZ(c, px0[k])) \
						++nnz1; \
			} while (0)

			SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		} else {

			Rcomplex *px0 = COMPLEX(x0);
			for (k = 0; k < nnz0; ++k)
				if ((pi0[k] == pj0[k]) ? NZ(d, px0[k].r) : NZ(z, px0[k]))
					++nnz1;

		}

		if (nnz1 == nnz0) {
			UNPROTECT(4); /* x0, j0, i0, from */
			return from;
		}

		PROTECT(to = newObject(class));

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
			j1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
			x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);

		if (ct != 'C') {

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) \
					if (NZ(c, px0[k])) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						*(px1++) = px0[k]; \
					} \
			} while (0)

			SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		} else {

			Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
			for (k = 0; k < nnz0; ++k)
				if ((pi0[k] == pj0[k]) ? NZ(d, px0[k].r) : NZ(z, px0[k])) {
					*(pi1++) = pi0[k];
					*(pj1++) = pj0[k];
					*(px1++) = px0[k];
				}

		}

		UNPROTECT(7); /* x1, j1, i1, to, x0, j0, i0 */

	}

	PROTECT(to);
	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && ct != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	UNPROTECT(2); /* to, from */
	return to;
}

/* drop0(<[CRT]sparseMatrix>, tol) */
SEXP R_sparse_drop0(SEXP s_from, SEXP s_tol)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	double tol;
	if (TYPEOF(s_tol) != REALSXP || LENGTH(s_tol) < 1 ||
	    ISNAN(tol = REAL(s_tol)[0]))
		Rf_error(_("'%s' is not a number"), "tol");

	return sparse_drop0(s_from, class, tol);
}

SEXP sparse_diag_U2N(SEXP from, const char *class)
{
	/* defined in ./diag.c : */
	SEXP R_sparse_diag_set(SEXP, SEXP);

	if (class[1] != 't' || DIAG(from) == 'N')
		return from;
	SEXP value = PROTECT(Rf_ScalarLogical(1)),
		to = R_sparse_diag_set(from, value);
	UNPROTECT(1); /* value */
	return to;
}

/* diagU2N(<[CRT]sparseMatrix>) */
SEXP R_sparse_diag_U2N(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_diag_U2N(s_from, class);
}

SEXP sparse_diag_N2U(SEXP from, const char *class)
{
	/* defined in ./band.c : */
	SEXP sparse_band(SEXP, const char *, int, int);
	
	if (class[1] != 't' || DIAG(from) != 'N')
		return from;
	SEXP to;
	int n = DIM(from)[1];
	if (n == 0) {
		PROTECT(to = newObject(class));
		if (UPLO(from) != 'U')
			SET_UPLO(to);
	}
	else if (UPLO(from) == 'U')
		PROTECT(to = sparse_band(from, class,  1,  n));
	else
		PROTECT(to = sparse_band(from, class, -n, -1));
	SET_DIAG(to);
	UNPROTECT(1); /* to */
	return from;
}

/* diagN2U(<[CRT]sparseMatrix>) */
SEXP R_sparse_diag_N2U(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_diag_N2U(s_from, class);
}

#define MAP(i) (map) ? map[i] : i

#define nCAST(x) (1)
#define lCAST(x) (x)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void Csparse_colsum(SEXP obj, const char *class,
                    int m, int n, char ul, char ct, char nu, int narm, int mean,
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
                    int m, int n, char ul, char ct, char nu, int narm, int mean,
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
                    int m, int n, char ul, char ct, char nu, int narm, int mean,
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

SEXP sparse_marginsum(SEXP obj, const char *class, int mg, int narm, int mean,
                      int sparse)
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
			Csparse_rowsum(obj, class, m, n, ul, ct, nu, narm, mean, ans,
			               Matrix_iSym);
			break;
		case 'R':
			if (class[1] == 's')
			Csparse_rowsum(obj, class, n, m, ul, ct, nu, narm, mean, ans,
			               Matrix_jSym);
			else
			Csparse_colsum(obj, class, n, m, ul, ct, nu, narm, mean, ans);
			break;
		case 'T':
			Tsparse_colsum(obj, class, n, m, ul, ct, nu, narm, mean, ans,
			               Matrix_jSym, Matrix_iSym);
			break;
		default:
			break;
		}
	else
		switch (class[2]) {
		case 'C':
			if (class[1] == 's')
			Csparse_rowsum(obj, class, m, n, ul, ct, nu, narm, mean, ans,
			               Matrix_iSym);
			else
			Csparse_colsum(obj, class, m, n, ul, ct, nu, narm, mean, ans);
			break;
		case 'R':
			Csparse_rowsum(obj, class, n, m, ul, ct, nu, narm, mean, ans,
			               Matrix_jSym);
			break;
		case 'T':
			Tsparse_colsum(obj, class, m, n, ul, ct, nu, narm, mean, ans,
			               Matrix_iSym, Matrix_jSym);
			break;
		default:
			break;
		}

	UNPROTECT(1); /* ans */
	return ans;
}

/* (row|col)(Sums|Means)(<[CRT]sparseMatrix>, na.rm=) */
SEXP R_sparse_marginsum(SEXP s_obj, SEXP s_margin, SEXP s_narm, SEXP s_mean,
                        SEXP s_sparse)
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

#undef TRY_INCREMENT

	return ans;
}

/* sum(<[CRT]sparseMatrix>, na.rm=) */
SEXP R_sparse_sum(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return sparse_sum(s_obj, class, narm);
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

/* prod(<[CRT]sparseMatrix>, na.rm=) */
SEXP R_sparse_prod(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return sparse_prod(s_obj, class, narm);
}
