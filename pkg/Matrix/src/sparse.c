#include <math.h> /* fabs, hypot */
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
		error(_("number of triplets to be aggregated exceeds %s"),
		      "2^31-1");
	PROTECT(i0);
	int *pi0 = INTEGER(i0), nnz = (int) XLENGTH(i0), nnz_ = nnz;

	SEXP to;

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym)), i1, j1;
	int *pj0 = INTEGER(j0), *pi1 = NULL, *pj1 = NULL, *iwork = NULL;
	size_t liwork = (size_t) ((int_fast64_t) n + 1 + n + m + nnz),
		lwork = (size_t) nnz;
	Matrix_Calloc(iwork, liwork, int);

#define AGGR(c) \
	do { \
		c##TYPE *px0 = NULL, *px1 = NULL, *work = NULL; \
		c##IF_NPATTERN( \
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
		px0 = c##PTR(x0); \
		Matrix_Calloc(work, lwork, c##TYPE); \
		); \
		c##spaggr(pj1, pi1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
		if (nnz != nnz_) { \
		PROTECT(to = newObject(class)); \
		PROTECT(i1 = allocVector(INTSXP, nnz)), \
		PROTECT(j1 = allocVector(INTSXP, nnz)); \
		pi1 = INTEGER(i1); \
		pj1 = INTEGER(j1); \
		SET_SLOT(to, Matrix_iSym, i1); \
		SET_SLOT(to, Matrix_jSym, j1); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(allocVector(c##TYPESXP, nnz)); \
		px1 = c##PTR(x1); \
		SET_SLOT(to, Matrix_xSym, x1); \
		); \
		c##spaggr(pj1, pi1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
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

	SWITCH5(class[0], AGGR);

#undef AGGR

	Matrix_Free(iwork, liwork);
	UNPROTECT(2); /* j0, i0 */

	if (nnz == nnz_)
		return from;

	PROTECT(to);

	dim = GET_SLOT(to, Matrix_DimSym);
	pdim = INTEGER(dim);
	pdim[0] = m;
	pdim[1] = n;

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		char ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */
	}

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
	PROTECT(from = sparse_aggregate(from, class));

	if (class[0] == 'n') {
		UNPROTECT(1); /* from */
		return from;
	}

	SEXP trans = R_NilValue;
	char ct = '\0';
	if (class[1] == 's' && class[0] == 'z') {
		trans = GET_SLOT(from, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}
	PROTECT(trans);

	SEXP to;
	int strict = ISNAN(tol) || tol <= 0.0;

#define NZ(c, x) ((strict) ? c##NOT_ZERO(x) : c##NOT_ZERO_TOL(x, tol))

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			n = (int) (XLENGTH(p0) - 1), j, k, kend, nnz0 = pp0[n], nnz1 = 0;
		pp0++;

		if (ct != 'C') {

#define COUNT(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0); \
				for (k = 0; k < nnz0; ++k) \
					if (NZ(c, px0[k])) \
						++nnz1; \
			} while (0)

			SWITCH4(class[0], COUNT);

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
			UNPROTECT(5); /* x0, i0, p0, trans, from */
			return from;
		}

		PROTECT(to = newObject(class));

		SEXP p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0))),
			i1 = PROTECT(allocVector(INTSXP, nnz1)),
			x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
		int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		SET_SLOT(to, Matrix_xSym, x1);

		if (ct != 'C') {

#define DROP0(c) \
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

			SWITCH4(class[0], DROP0);

#undef DROP0

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

			SWITCH4(class[0], COUNT);

#undef COUNT

		} else {

			Rcomplex *px0 = COMPLEX(x0);
			for (k = 0; k < nnz0; ++k)
				if ((pi0[k] == pj0[k]) ? NZ(d, px0[k].r) : NZ(z, px0[k]))
					++nnz1;

		}

		if (nnz1 == nnz0) {
			UNPROTECT(5); /* x0, j0, i0, trans, from */
			return from;
		}

		PROTECT(to = newObject(class));

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1)),
			x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);
		SET_SLOT(to, Matrix_xSym, x1);

		if (ct != 'C') {

#define DROP0(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) \
					if (NZ(c, px0[k])) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						*(px1++) = px0[k]; \
					} \
			} while (0)

			SWITCH4(class[0], DROP0);

#undef DROP0

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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 's' && class[0] == 'z') {
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	UNPROTECT(3); /* to, trans, from */
	return to;
}

/* drop0(<[CRT]sparseMatrix>, tol) */
SEXP R_sparse_drop0(SEXP s_from, SEXP s_tol)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	double tol;
	if (TYPEOF(s_tol) != REALSXP || LENGTH(s_tol) < 1 ||
	    ISNAN(tol = REAL(s_tol)[0]))
		error(_("'%s' is not a number"), "tol");

	return sparse_drop0(s_from, class, tol);
}

SEXP sparse_diag_U2N(SEXP from, const char *class)
{
	if (class[1] != 't')
		return from;

	SEXP diag = GET_SLOT(from, Matrix_diagSym);
	char di = CHAR(STRING_ELT(diag, 0))[0];
	if (di == 'N')
		return from;

	SEXP value = PROTECT(ScalarLogical(1));
	from = R_sparse_diag_set(from, value);
	UNPROTECT(1); /* value */

	return from;
}

/* diagU2N(<[CRT]sparseMatrix>) */
SEXP R_sparse_diag_U2N(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_diag_U2N(s_from, class);
}

SEXP sparse_diag_N2U(SEXP from, const char *class)
{
	if (class[1] != 't')
		return from;

	SEXP diag = GET_SLOT(from, Matrix_diagSym);
	char di = CHAR(STRING_ELT(diag, 0))[0];
	if (di != 'N')
		return from;
	PROTECT(diag = mkString("U"));

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int n = INTEGER(dim)[0];

	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	if (n == 0) {
		PROTECT(from = newObject(class));
		if (ul != 'U')
			SET_SLOT(from, Matrix_DimSym, uplo);
	}
	else if (ul == 'U')
		PROTECT(from = sparse_band(from, class,  1,  n));
	else
		PROTECT(from = sparse_band(from, class, -n, -1));
	SET_SLOT(from, Matrix_diagSym, diag);

	UNPROTECT(3); /* from, uplo, diag */
	return from;
}

/* diagN2U(<[CRT]sparseMatrix>) */
SEXP R_sparse_diag_N2U(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_diag_N2U(s_from, class);
}

SEXP sparse_band(SEXP from, const char *class, int a, int b)
{
	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
	if ((m == 0 || n == 0 || (a <= 1 - m && b >= n - 1)) &&
	    (m != n || n > 1 || class[1] == 't'))
		return from;

	int ge, sy, tr;
	tr = class[1] == 't' || (m == n && (a >= 0 || b <= 0));
	sy = !tr && class[1] == 's' && a == -b;
	ge = !tr && !sy;

	char ul0 = '\0', ct0 = '\0', di0 = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul0 = CHAR(STRING_ELT(uplo, 0))[0];
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct0 = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		/* Be fast if band contains entire triangle */
		if ((ul0 == 'U') ? (a <= 0 && b >= n - 1) : (a <= 1 - m && b >= 0))
			return from;
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di0 = CHAR(STRING_ELT(diag, 0))[0];
	}

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (ge) ? 'g' : ((sy) ? 's' : 't');
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

	dim = GET_SLOT(to, Matrix_DimSym);
	pdim = INTEGER(dim);
	pdim[0] = m;
	pdim[1] = n;

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's' && !sy)
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul1 = (tr && class[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ul0;
	if (ul1 != '\0' && ul1 != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (ct0 != '\0' && ct0 != 'C' && sy) {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (di0 != '\0' && di0 != 'N' && tr && a <= 0 && b >= 0) {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	int d;

	if (class[2] != 'T') {

		if (class[2] == 'R') {
			int tmp;
			tmp = m; m =  n; n =  tmp;
			tmp = a; a = -b; b = -tmp;
		}

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0)));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			*pp1 = INTEGER(p1), *iwork = NULL,
			j, k, kend, nnz0 = pp0[n], nnz1 = 0;
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] != 's' || sy) {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
			if (nnz0 == nnz1) {
				SET_SLOT(to, iSym, i0);
				if (class[0] != 'n') {
					SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
					SET_SLOT(to, Matrix_xSym, x0);
					UNPROTECT(1); /* x0 */
				}
				UNPROTECT(4); /* p1, i0, p0, to */
				return to;
			}
		} else {
			memset(pp1, 0, sizeof(int) * (size_t) n);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++pp1[j];
					if (d != 0 && -d >= a && -d <= b)
						++pp1[pi0[k]];
					++k;
				}
			}
			Matrix_Calloc(iwork, n, int);
			for (j = 0; j < n; ++j) {
				iwork[j] = nnz1;
				nnz1 += pp1[j];
				pp1[j] = nnz1;
			}
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define BAND(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz1)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] != 's' || sy) \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							*(pi1++) = pi0[k]; \
							c##IF_NPATTERN( \
							*(px1++) = px0[k]; \
							); \
						} \
						++k; \
					} \
				} \
			else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							c##IF_NPATTERN( \
							if (ct0 == 'C' && d == 0) \
								c##ASSIGN_PROJ_REAL(px1[iwork[j]], px0[k]); \
							else \
								c##ASSIGN_IDEN     (px1[iwork[j]], px0[k]); \
							); \
							pi1[iwork[j]++] = pi0[k]; \
						} \
						if (d != 0 && -d >= a && -d <= b) { \
							c##IF_NPATTERN( \
							if (ct0 == 'C') \
								c##ASSIGN_CONJ(px1[iwork[pi0[k]]], px0[k]); \
							else \
								c##ASSIGN_IDEN(px1[iwork[pi0[k]]], px0[k]); \
							); \
							pi1[iwork[pi0[k]]++] = j; \
						} \
						++k; \
					} \
				} \
				Matrix_Free(iwork, n); \
			} \
		} while (0)

		SWITCH5(class[0], BAND);

#undef BAND

		UNPROTECT(4); /* i1, p1, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;

		if (class[1] != 's' || sy) {
			for (k = 0; k < nnz0; ++k)
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
			if (nnz0 == nnz1) {
				SET_SLOT(to, Matrix_iSym, i0);
				SET_SLOT(to, Matrix_jSym, j0);
				if (class[0] != 'n') {
					SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
					SET_SLOT(to, Matrix_xSym, x0);
					UNPROTECT(1); /* x0 */
				}
				UNPROTECT(3); /* j0, i0, to */
				return to;
			}
		} else
			for (k = 0; k < nnz0; ++k) {
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
				if (d != 0 && -d >= a && -d <= b)
					++nnz1;
			}

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#define BAND(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz1)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] != 's' || sy) \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						c##IF_NPATTERN( \
						*(px1++) = px0[k]; \
						); \
					} \
				} \
			else \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						c##IF_NPATTERN( \
						if (ct0 == 'C' && d == 0) \
							c##ASSIGN_PROJ_REAL(*px1, px0[k]); \
						else \
							c##ASSIGN_IDEN     (*px1, px0[k]); \
						px1++; \
						); \
					} \
					if (d != 0 && -d >= a && -d <= b) { \
						*(pi1++) = pj0[k]; \
						*(pj1++) = pi0[k]; \
						c##IF_NPATTERN( \
						if (ct0 == 'C') \
							c##ASSIGN_CONJ(*px1, px0[k]); \
						else \
							c##ASSIGN_IDEN(*px1, px0[k]); \
						px1++; \
						); \
					} \
				} \
		} while (0)

		SWITCH5(class[0], BAND);

#undef BAND

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

/* band(<[CRT]sparseMatrix>, k1, k2), tri[ul](<[CRT]sparseMatrix>, k) */
/* NB: argument validation more or less copied from R_dense_band() */
SEXP R_sparse_band(SEXP s_from, SEXP s_a, SEXP s_b)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	SEXP dim = GET_SLOT(s_from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	int a, b;
	if (s_a == R_NilValue)
		a = -m;
	else if ((a = asInteger(s_a)) == NA_INTEGER || a < -m || a > n)
		error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		      "k1", a, "-Dim[1]", -m, "Dim[2]", n);
	if (s_b == R_NilValue)
		b = n;
	else if ((b = asInteger(s_b)) == NA_INTEGER || b < -m || b > n)
		error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		      "k2", b, "-Dim[1]", -m, "Dim[2]", n);
	else if (b < a)
		error(_("'%s' (%d) must be less than or equal to '%s' (%d)"),
		      "k1", a, "k2", b);

	return sparse_band(s_from, class, a, b);
}

SEXP sparse_diag_get(SEXP obj, const char *class, int names)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP ans;

	if (di != '\0' && di != 'N') {

		PROTECT(ans = allocUnit  (kindToType(class[0]), r));

	} else if (class[2] != 'T') {

		PROTECT(ans = allocVector(kindToType(class[0]), r));

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend,
			upper = (class[2] == 'C') == (ul == 'U');
		pp++;

#define DIAG(c) \
		do { \
			c##TYPE *pa = c##PTR(ans); \
			c##IF_NPATTERN( \
			SEXP x = GET_SLOT(obj, Matrix_xSym); \
			c##TYPE *px = c##PTR(x); \
			); \
			if (class[1] == 'g') \
				for (j = 0, k = 0; j < r; ++j) { \
					pa[j] = c##ZERO; \
					kend = pp[j]; \
					while (k < kend) { \
						if (pi[k] != j) \
							++k; \
						else { \
							pa[j] = c##IFELSE_NPATTERN(px[k], c##UNIT); \
							k = kend; \
						} \
					} \
				} \
			else \
				for (j = 0, k = 0; j < r; ++j) { \
					kend = pp[j]; \
					pa[j] = (k < kend && pi[(upper) ? kend - 1 : k] == j) \
						? c##IFELSE_NPATTERN(px[(upper) ? kend - 1 : k], c##UNIT) \
						: c##ZERO; \
					k = kend; \
				} \
		} while (0)

		SWITCH5(class[0], DIAG);

#undef DIAG

		UNPROTECT(2); /* i, p */

		if (class[0] == 'z' && class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), NULL, (size_t) r);

	} else {

		PROTECT(ans = allocVector(kindToType(class[0]), r));

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

#define DIAG(c) \
		do { \
			c##TYPE *pa = c##PTR(ans); \
			memset(pa, 0, sizeof(c##TYPE) * (size_t) r); \
			c##IF_NPATTERN( \
			SEXP x = GET_SLOT(obj, Matrix_xSym); \
			c##TYPE *px = c##PTR(x); \
			); \
			for (k = 0; k < kend; ++k) \
				if (pi[k] == pj[k]) \
					c##INCREMENT_IDEN(pa[pi[k]], px[k]); \
		} while (0)

		SWITCH5(class[0], DIAG);

#undef DIAG

		UNPROTECT(2); /* j0, i0 */

		if (class[0] == 'z' && class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), NULL, (size_t) r);

	}

	if (names) {
		/* NB: The logic here must be adjusted once the validity method
		   for 'symmetricMatrix' enforces symmetric 'Dimnames'
		*/
		SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (cn == R_NilValue) {
			if (class[1] == 's' && rn != R_NilValue)
				setAttrib(ans, R_NamesSymbol, rn);
		} else {
			if (class[1] == 's')
				setAttrib(ans, R_NamesSymbol, cn);
			else if (rn != R_NilValue &&
			         (rn == cn || equalString(rn, cn, r)))
				setAttrib(ans, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(1); /* ans */
	return ans;
}

/* diag(<[CRT]sparseMatrix>, names=) */
SEXP R_sparse_diag_get(SEXP s_obj, SEXP s_names)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int names;
	VALID_LOGIC2(s_names, names);

	return sparse_diag_get(s_obj, class, names);
}

SEXP sparse_diag_set(SEXP from, const char *class, SEXP value)
{
	SEXP to = PROTECT(newObject(class));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	int mode = 0;

	if (class[2] != 'T') {

		if (class[2] == 'R') {
			int tmp;
			tmp = m; m = n; n = tmp;
		}

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0)));
		int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0),
			j, k, kend, nd0 = 0, nd1 = 0,
			upper = (class[2] == 'C') == (ul == 'U');
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] == 'g') {
			for (j = 0, k = 0; j < r; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] < j)
						++k;
					else {
						if (pi0[k] == j)
							++nd0;
						k = kend;
					}
				}
				pp1[j] = kend - nd0;
			}
			for (j = r; j < n; ++j)
				pp1[j] = pp0[j] - nd0;
		}
		else if (class[1] == 's' || di == 'N')
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[(upper) ? kend - 1 : k] == j)
					++nd0;
				k = kend;
				pp1[j] = kend - nd0;
			}
		else
			for (j = 0; j < n; ++j)
				pp1[j] = pp0[j];

#define COUNT(c) \
		do { \
			c##TYPE *pv = c##PTR(value); \
			if (LENGTH(value) == r) { \
				mode = 1; \
				for (j = 0; j < r; ++j) { \
					if (c##NOT_ZERO(pv[j])) \
						++nd1; \
					pp1[j] += nd1; \
				} \
				for (j = r; j < n; ++j) \
					pp1[j] += nd1; \
			} else if (c##NOT_ZERO(pv[0])) { \
				mode = 2; \
				for (j = 0; j < r; ++j) \
					pp1[j] += ++nd1; \
				for (j = r; j < n; ++j) \
					pp1[j] +=   nd1; \
			} \
		} while (0)

		SWITCH4(class[0], COUNT);

#undef COUNT

		if (nd1 - nd0 > INT_MAX - pp0[n - 1])
			error(_("%s cannot exceed %s"), "p[length(p)]", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define DIAG(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, pp1[n - 1])); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			c##TYPE *pv = c##PTR(value); \
			for (j = 0, k = 0; j < r; ++j) { \
				kend = pp0[j]; \
				while (k < kend && pi0[k] < j) { \
					*(pi1++) = pi0[k]; \
					c##IF_NPATTERN( \
					*(px1++) = px0[k]; \
					); \
					++k; \
				} \
				if (k < kend && pi0[k] == j) \
					++k; \
				if (mode > 1 || (mode > 0 && c##NOT_ZERO(pv[j]))) { \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = pv[(mode > 1) ? 0 : j]; \
					); \
				} \
				while (k < kend) { \
					*(pi1++) = pi0[k]; \
					c##IF_NPATTERN( \
					*(px1++) = px0[k]; \
					); \
					++k; \
				} \
			} \
			for (j = r; j < n; ++j) { \
				kend = pp0[j]; \
				while (k < kend) { \
					*(pi1++) = pi0[k]; \
					c##IF_NPATTERN( \
					*(px1++) = px0[k]; \
					); \
					++k; \
				} \
			} \
		} while (0)

		SWITCH5(class[0], DIAG);

#undef DIAG

		UNPROTECT(4); /* i1, p1, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), j;
		R_xlen_t k, nd0 = 0, nd1 = 0, nnz0 = XLENGTH(i0), nnz1 = nnz0;

		if (di == '\0' || di == 'N')
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					++nd0;

#define COUNT(c) \
		do { \
			c##TYPE *pv = c##PTR(value); \
			if (LENGTH(value) == r) { \
				mode = 1; \
				for (j = 0; j < r; ++j) \
					if (c##NOT_ZERO(pv[j])) \
						++nd1; \
			} else if (c##NOT_ZERO(pv[0])) { \
				mode = 2; \
				nd1 = r; \
			} \
		} while (0)

		SWITCH4(class[0], COUNT);

#undef COUNT

		if (nd1 - nd0 > R_XLEN_T_MAX - nnz0)
			error(_("%s cannot exceed %s"), "length(i)", "R_XLEN_T_MAX");
		nnz1 += nd1 - nd0;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#define DIAG(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz1)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			c##TYPE *pv = c##PTR(value); \
			for (k = 0; k < nnz0; ++k) { \
				if (pi0[k] != pj0[k]) { \
					*(pi1++) = pi0[k]; \
					*(pj1++) = pj0[k]; \
					c##IF_NPATTERN( \
					*(px1++) = px0[k]; \
					); \
				} \
			} \
			if (mode > 1) \
			for (j = 0; j < r; ++j) { \
				*(pi1++) = *(pj1++) = j; \
				c##IF_NPATTERN( \
				*(px1++) = pv[0]; \
				); \
			} \
			else if (mode > 0) \
			for (j = 0; j < r; ++j) { \
				if (c##NOT_ZERO(pv[j])) { \
					*(pi1++) = *(pj1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = pv[j]; \
					); \
				} \
			} \
		} while (0)

		SWITCH5(class[0], DIAG);

#undef DIAG

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

/* diag(<[CRT]sparseMatrix>) <- value */
SEXP R_sparse_diag_set(SEXP s_from, SEXP s_value)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);
	SEXPTYPE tx = kindToType(class[0]), tv = TYPEOF(s_value);

	switch (tv) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		error(_("replacement diagonal has incompatible type \"%s\""),
		      type2char(tv));
		break;
	}

	SEXP dim = GET_SLOT(s_from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	if (XLENGTH(s_value) != 1 && XLENGTH(s_value) != r)
		error(_("replacement diagonal has wrong length"));

	if (tv <= tx) {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		if (class[1] == 's' && class[0] == 'z' && tv == tx) {
			SEXP trans = GET_SLOT(s_from, Matrix_transSym);
			int ct = CHAR(STRING_ELT(trans, 0))[0];
			if (ct == 'C') {
				PROTECT(s_from = sparse_as_general(s_from, class));
				class = Matrix_class(s_from, valid_sparse, 6, __func__);
				UNPROTECT(1); /* s_from */
			}
		}
		PROTECT(s_from);
		PROTECT(s_value = coerceVector(s_value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
#ifndef MATRIX_ENABLE_IMATRIX
		if (tv == INTSXP) {
		PROTECT(s_from = sparse_as_kind(s_from, class, 'd'));
		PROTECT(s_value = coerceVector(s_value, REALSXP));
		} else {
#endif
		PROTECT(s_from = sparse_as_kind(s_from, class, typeToKind(tv)));
		PROTECT(s_value);
#ifndef MATRIX_ENABLE_IMATRIX
		}
#endif
		class = Matrix_class(s_from, valid_sparse, 6, __func__);
	}

	s_from = sparse_diag_set(s_from, class, s_value);
	UNPROTECT(2); /* s_value, s_from */
	return s_from;
}

SEXP sparse_transpose(SEXP from, const char *class, char op_ct, int lazy)
{
	if (class[0] != 'z')
		op_ct = '\0';

	SEXP to;
	if (class[2] == 'T' || !lazy)
		PROTECT(to = newObject(class));
	else {
		char cl[] = "...Matrix";
		cl[0] = class[0];
		cl[1] = class[1];
		cl[2] = (class[2] == 'C') ? 'R' : 'C';
		PROTECT(to = newObject(cl));
	}

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	dim = GET_SLOT(to, Matrix_DimSym);
	pdim = INTEGER(dim);
	pdim[0] = n;
	pdim[1] = m;

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's' || class[1] == 'p')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_reversed_DimNames(to, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul == 'U') {
			PROTECT(uplo = mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		char ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	if (class[2] != 'T') {

		if (class[2] == 'R') {
			int tmp;
			tmp = m; m = n; n = tmp;
		}

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			jSym = (class[2] == 'C') ? Matrix_jSym : Matrix_iSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), nnz = INTEGER(p0)[n];

		if (lazy) {

		SET_SLOT(to, Matrix_pSym, p0);
		SET_SLOT(to,        jSym, i0);

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (op_ct != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(allocVector(CPLXSXP, nnz));
				zvconj(COMPLEX(x1), COMPLEX(x0), (size_t) nnz);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
			UNPROTECT(1); /* x0 */
		}

		} else {

		SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) m + 1)),
			i1 = PROTECT(allocVector(INTSXP, nnz));
		int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1), *iwork = NULL;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		Matrix_Calloc(iwork, m, int);

#define TRANS(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz)); \
			px0 = c##PTR(x0); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			c##sptrans(pp1, pi1, px1, pp0, pi0, px0, m, n, op_ct, iwork); \
		} while (0)

		SWITCH5(class[0], TRANS);

#undef TRANS

		Matrix_Free(iwork, m);
		UNPROTECT(2); /* j1, p1 */

		}

		UNPROTECT(2); /* i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_iSym, j0);
		SET_SLOT(to, Matrix_jSym, i0);

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (op_ct != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(allocVector(CPLXSXP, XLENGTH(x0)));
				zvconj(COMPLEX(x1), COMPLEX(x0), (size_t) XLENGTH(x0));
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
			UNPROTECT(1); /* x0 */
		}

		UNPROTECT(2); /* j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

/* t(<[CRT]sparseMatrix>) */
SEXP R_sparse_transpose(SEXP s_from, SEXP s_trans, SEXP s_lazy)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	int lazy;
	VALID_LOGIC2(s_lazy, lazy);

	return sparse_transpose(s_from, class, ct, lazy);
}

SEXP sparse_force_symmetric(SEXP from, const char *class, char op_ul, char op_ct)
{
	if (class[0] != 'z')
		op_ct = '\0';

	char
		ul0 = '\0', ul1 = 'U',
		ct0 = '\0', ct1 = (class[0] == 'z') ? 'C' : '\0',
		di0 = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul0 = ul1 = CHAR(STRING_ELT(uplo, 0))[0];
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct0 = ct1 = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di0 = CHAR(STRING_ELT(diag, 0))[0];
	}
	if (op_ul != '\0')
		ul1 = op_ul;
	if (op_ct != '\0')
		ct1 = op_ct;

	if (class[1] == 's' && ul0 == ul1 && ct0 == ct1)
		return from;

	char cl[] = ".s.Matrix";
	cl[0] = class[0];
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to symmetrize a non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (ul1 != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (ct1 != 'C' && ct1 != '\0') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	int upper = (class[0] != 'R') == (ul1 == 'U');

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			j, k, kend, nnz0 = pp0[n], nnz1 = 0;
		pp0++;

		if ((class[1] == 's' && ul0 == ul1) ||
		    (class[1] == 't' && ul0 == ul1 && di0 == 'N')) {

		SET_SLOT(to, Matrix_pSym, p0);
		SET_SLOT(to,        iSym, i0);

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (class[1] != 's' || ct0 != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(allocVector(CPLXSXP, nnz0));
				Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
				for (j = 0, k = 0; j < n; ++j) {
					kend = pp0[j];
					while (k < kend) {
						if (pi0[k] == j)
						zASSIGN_PROJ_REAL(px1[k], px0[k]);
						else
						zASSIGN_IDEN     (px1[k], px0[k]);
						++k;
					}
				}
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
			UNPROTECT(1); /* x0 */
		}

		} else {

		SEXP p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0)));
		int *pp1 = INTEGER(p1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] == 'g')
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((upper) ? pi0[k] <= j : pi0[k] >= j)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		else if (class[1] == 's')
			nnz1 = nnz0;
		else if (di0 == 'N')
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[(upper) ? k : kend - 1] == j)
					++nnz1;
				pp1[j] = nnz1;
				k = kend;
			}
		else if (ul0 == ul1) {
			for (j = 0; j < n; ++j)
				pp1[j] = ++nnz1 + pp0[j];
			nnz1 = nnz0 + n;
		}
		else
			for (j = 0; j < n; ++j)
				pp1[j] = ++nnz1;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define FORCE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz1)); \
			px0 = c##PTR(x0); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] == 'g') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((upper) ? pi0[k] <= j : pi0[k] >= j) { \
							*(pi1++) = pi0[k]; \
							c##IF_NPATTERN( \
							*(px1++) = px0[k]; \
							); \
						} \
						++k; \
					} \
				} \
			else if (class[1] == 's') { \
				int *iwork = NULL; \
				Matrix_Calloc(iwork, n, int); \
				c##sptrans(pp1, pi1, px1, pp0, pi0, px0, n, n, ct0, iwork); \
				Matrix_Free(iwork, n); \
				if (ct0 == 'C') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp1[j]; \
					if (k < kend && pi1[(upper) ? kend - 1 : k] == j) \
						c##SET_PROJ_REAL(px1[(upper) ? kend - 1 : k]); \
					k = kend; \
				} \
			} \
			else if (di0 == 'N') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					if (k < kend && pi0[(upper) ? k : kend - 1] == j) { \
						*(pi1++) = j; \
						c##IF_NPATTERN( \
						*(px1++) = px0[(upper) ? k : kend - 1]; \
						); \
					} \
					k = kend; \
				} \
			else if (ul0 == ul1) \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					if (!upper) { \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
					} \
					while (k < kend) { \
						*(pi1++) = pi0[k]; \
						c##IF_NPATTERN( \
						*(px1++) = px0[k]; \
						); \
						++k; \
					} \
					if (upper) { \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
					} \
				} \
			else \
				for (j = 0; j < n; ++j) { \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
				} \
		} while (0)

		SWITCH5(class[0], FORCE);

#undef FORCE

		UNPROTECT(2); /* i1, p1 */

		}

		UNPROTECT(2); /* i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;

		if ((class[1] == 's') ||
		    (class[1] == 't' && ul0 == ul1 && di0 == 'N')) {

		if (class[1] != 's' || ul0 == ul1) {
			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);
		} else {
			SET_SLOT(to, Matrix_iSym, j0);
			SET_SLOT(to, Matrix_jSym, i0);
		}

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (class[1] != 's' || ct0 != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(allocVector(CPLXSXP, XLENGTH(x0)));
				Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
				for (k = 0; k < nnz0; ++k)
					if (pi0[k] == pj0[k])
						zASSIGN_PROJ_REAL(px1[k], px0[k]);
					else if (ul0 == ul1)
						zASSIGN_IDEN     (px1[k], px0[k]);
					else
						zASSIGN_CONJ     (px1[k], px0[k]);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
			UNPROTECT(1); /* x0 */
		}

		} else {

		if (class[1] == 'g' || di0 == 'N') {
			for (k = 0; k < nnz0; ++k)
				if ((upper) ? pi0[k] <= pj0[k] : pi0[k] >= pj0[k])
					++nnz1;
		}
		else
			nnz1 = (ul0 == ul1) ? n + nnz0 : n;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1), j;
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#define FORCE(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz1)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] == 'g' || di0 == 'N') \
				for (k = 0; k < nnz0; ++k) { \
					if ((upper) ? pi0[k] <= pj0[k] : pi0[k] >= pj0[k]) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						c##IF_NPATTERN( \
						*(px1++) = px0[k]; \
						); \
					} \
				} \
			else { \
				if (ul0 == ul1) { \
					memcpy(pi1, pi0, sizeof(    int) * (size_t) nnz0); \
					memcpy(pj1, pj0, sizeof(    int) * (size_t) nnz0); \
					pi1 += nnz0; \
					pj1 += nnz0; \
					c##IF_NPATTERN( \
					memcpy(px1, px0, sizeof(c##TYPE) * (size_t) nnz0); \
					px1 += nnz0; \
					); \
				} \
				for (j = 0; j < n; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
				} \
			} \
		} while (0)

		SWITCH5(class[0], FORCE);

#undef FORCE

		UNPROTECT(2); /* j1, i1 */

		}

		UNPROTECT(2); /* j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

/* forceSymmetric(<[CRT]sparseMatrix>, uplo, trans) */
SEXP R_sparse_force_symmetric(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	char ul = '\0', ct = '\0';
	if (s_uplo != R_NilValue)
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return sparse_force_symmetric(s_from, class, ul, ct);
}

SEXP sparse_symmpart(SEXP from, const char *class, char op_ul, char op_ct)
{
	/* defined in ./coerce.c : */
	SEXP sparse_as_kind(SEXP, const char *, char);
	PROTECT(from = sparse_as_kind(from, class, ','));

	if (class[1] != 'g')
		op_ul = '\0';
	if (class[0] != 'z')
		op_ct = '\0';

	char ct = '\0';
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}

	if (class[1] == 's' && op_ct == ct) {
		UNPROTECT(1); /* from */
		return from;
	}

	char cl[] = ".s.Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error((op_ct == 'C')
		      ? _("attempt to get Hermitian part of non-square matrix")
		      : _("attempt to get symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	char ul = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	if (op_ul != '\0' && op_ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (op_ct != '\0' && op_ct != 'C') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	int upper = (class[0] != 'R') == (op_ul == 'U' || ul == 'U');

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			j, k, kend, nnz0 = pp0[n], nnz1;
		pp0++;

		if (class[1] == 'g') {

			int *iwork = NULL;
			size_t liwork = (size_t) ((int_fast64_t) n + n + 1 + nnz0);
			Matrix_Calloc(iwork, liwork, int);

			int *pp0_ = iwork + n + 1, *pi0_ = iwork + n + n + 1;
			nsptrans(pp0_ - 1, pi0_, NULL, pp0 - 1, pi0, NULL, n, n, 'T', iwork);
			memcpy(iwork, pp0 - 1, sizeof(int) * (size_t) n);

			SEXP p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0)));
			int *pp1 = INTEGER(p1), k_, kend_;
			*(pp1++) = 0;
			SET_SLOT(to, Matrix_pSym, p1);

			for (j = 0, k = 0, k_ = 0; j < n; ++j) {
				kend  = pp0 [j];
				kend_ = pp0_[j];
				pp1[j] = pp1[j - 1];
				while (k < kend) {
					if (pi0[k] > j)
						k = kend;
					else {
						while (k_ < kend_ && pi0_[k_] < pi0[k]) {
							++pp1[j];
							++k_;
						}
						++pp1[j];
						if (k_ < kend_ && pi0_[k_] == pi0[k])
							++k_;
						++k;
					}
				}
				while (k_ < kend_) {
					if (pi0_[k_] > j)
						k_ = kend_;
					else {
						++pp1[j];
						++k_;
					}
				}
			}
			nnz1 = pp1[n - 1];

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), l;
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				if (op_ct == 'C') \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (upper && pi0[k] > j) \
							k = kend; \
						else if (!upper && pi0[k] < j) \
							++k; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								*pi1 = pi0_[k_]; \
								c##ASSIGN_CONJ(*px1, px0[l]); \
								c##MULTIPLY(*px1, 0.5); \
								++k_; ++pi1; ++px1; \
							} \
							l = iwork[j]++; \
							*pi1 = pi0[k]; \
							if (pi0[k] == j) { \
								c##ASSIGN_PROJ_REAL(*px1, px0[k]); \
								++k_; \
							} else { \
								c##ASSIGN_IDEN(*px1, px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##INCREMENT_CONJ(*px1, px0[l]); \
									++k_; \
								} \
								c##MULTIPLY(*px1, 0.5); \
							} \
							++k; ++pi1; ++px1; \
						} \
					} \
					while (k_ < kend_) { \
						if (upper && pi0_[k_] > j) \
							k_ = kend_; \
						else if (!upper && pi0_[k_] < j) \
							++k_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							*pi1 = pi0_[k_]; \
							c##ASSIGN_CONJ(*px1, px0[l]); \
							c##MULTIPLY(*px1, 0.5); \
							++k_; ++pi1; ++px1; \
						} \
					} \
				} \
				else \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (upper && pi0[k] > j) \
							k = kend; \
						else if (!upper && pi0[k] < j) \
							++k; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								*pi1 = pi0_[k_]; \
								c##ASSIGN_IDEN(*px1, px0[l]); \
								c##MULTIPLY(*px1, 0.5); \
								++k_; ++pi1; ++px1; \
							} \
							l = iwork[j]++; \
							*pi1 = pi0[k]; \
							if (pi0[k] == j) { \
								c##ASSIGN_IDEN(*px1, px0[k]); \
								++k_; \
							} else { \
								c##ASSIGN_IDEN(*px1, px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##INCREMENT_IDEN(*px1, px0[l]); \
									++k_; \
								} \
								c##MULTIPLY(*px1, 0.5); \
							} \
							++k; ++pi1; ++px1; \
						} \
					} \
					while (k_ < kend_) { \
						if (upper && pi0_[k_] > j) \
							k_ = kend_; \
						else if (!upper && pi0_[k_] < j) \
							++k_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							*pi1 = pi0_[k_]; \
							c##ASSIGN_IDEN(*px1, px0[l]); \
							c##MULTIPLY(*px1, 0.5); \
							++k_; ++pi1; ++px1; \
						} \
					} \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			Matrix_Free(iwork, liwork);
			UNPROTECT(3); /* x1, i1, p1 */

		} else if (class[1] == 's') {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			/* Symmetric part of Hermitian matrix is real part */
			/* Hermitian part of symmetric matrix is real part */
			SEXP x1 = PROTECT(allocVector(CPLXSXP, nnz0));
			zvreal(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else if (di == 'N') {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			SEXP x1 = PROTECT(allocVector(TYPEOF(x0), nnz0));
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] != j) { \
							c##ASSIGN_IDEN(*px1, px0[k]); \
							c##MULTIPLY(*px1, 0.5); \
						} \
						else if (op_ct == 'C') \
							c##ASSIGN_PROJ_REAL(*px1, px0[k]); \
						else \
							c##ASSIGN_IDEN(*px1, px0[k]); \
						++k; ++px1; \
					} \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			UNPROTECT(1); /* x1 */

		} else {

			nnz1 = nnz0 + n;

			SEXP p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0))),
				i1 = PROTECT(allocVector(INTSXP, nnz1)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1);
			*(pp1++) = 0;
			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					if (!upper) { \
						*pi1 = j; \
						*px1 = c##UNIT; \
						++pi1; ++px1; \
					} \
					while (k < kend) { \
						*pi1 = pi0[k]; \
						c##ASSIGN_IDEN(*px1, px0[k]); \
						c##MULTIPLY(*px1, 0.5); \
						++k; ++pi1; ++px1; \
					} \
					if (upper) { \
						*pi1 = j; \
						*px1 = c##UNIT; \
						++pi1; ++px1; \
					} \
					pp1[j] = kend + j + 1; \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			UNPROTECT(3); /* x1, i1, p1 */

		}

		UNPROTECT(3); /* x0, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1;

		if (class[1] == 'g') {

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz0)),
				j1 = PROTECT(allocVector(INTSXP, nnz0)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz0));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				if (op_ct == 'C') \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						if ((upper) ? *pi0 < *pj0 : *pi0 > *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						c##ASSIGN_IDEN(*px1, *px0); \
						} else { \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
						c##ASSIGN_CONJ(*px1, *px0); \
						} \
						c##MULTIPLY(*px1, 0.5); \
					} else { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						c##ASSIGN_PROJ_REAL(*px1, *px0); \
					} \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
				else \
				for (k = 0; k < nnz0; ++k) { \
					if ((upper) ? *pi0 <= *pj0 : *pi0 >= *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
					} else { \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
					} \
					c##ASSIGN_IDEN(*px1, *px0); \
					if (*pi0 != *pj0) \
					c##MULTIPLY(*px1, 0.5); \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			UNPROTECT(3); /* x1, j1, i1 */

		} else if (class[1] == 's') {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			/* Symmetric part of Hermitian matrix is real part */
			/* Hermitian part of symmetric matrix is real part */
			SEXP x1 = PROTECT(allocVector(CPLXSXP, nnz0));
			zvreal(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else if (di == 'N') {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			SEXP x1 = PROTECT(allocVector(TYPEOF(x0), nnz0));
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						c##ASSIGN_IDEN(*px1, *px0); \
						c##MULTIPLY(*px1, 0.5); \
					} \
					else if (op_ct == 'C') \
						c##ASSIGN_PROJ_REAL(*px1, *px0); \
					else \
						c##ASSIGN_IDEN(*px1, *px0); \
					++pi0; ++pj0; ++px0; ++px1; \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			UNPROTECT(1); /* x1 */

		} else {

			nnz1 = nnz0 + n;

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
				j1 = PROTECT(allocVector(INTSXP, nnz1)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1), j;
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) { \
					*pi1 = *pi0; \
					*pj1 = *pj0; \
					c##ASSIGN_IDEN(*px1, *px0); \
					c##MULTIPLY(*px1, 0.5); \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
				for (j = 0; j < n; ++j) { \
					*pi1 = *pj1 = j; \
					*px1 = c##UNIT; \
					++pi1; ++pj1; ++px1; \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			UNPROTECT(3); /* x1, j1, i1 */

		}

		UNPROTECT(3); /* x0, j0, i1 */

	}

	UNPROTECT(2); /* to, from */
	return to;
}

/* symmpart(<[CRT]sparseMatrix>, uplo, trans) */
SEXP R_sparse_symmpart(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	char ul, ct;
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return sparse_symmpart(s_from, class, ul, ct);
}

SEXP sparse_skewpart(SEXP from, const char *class, char op_ct)
{
	/* defined in ./coerce.c : */
	SEXP sparse_as_kind(SEXP, const char *, char);
	PROTECT(from = sparse_as_kind(from, class, ','));

	if (class[0] != 'z')
		op_ct = '\0';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 's') ? 's' : 'g';
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error((op_ct == 'C')
		      ? _("attempt to get skew-Hermitian part of non-square matrix")
		      : _("attempt to get skew-symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	char ul = '\0', ct = '\0';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U' && class[1] == 's')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	if (class[1] == 's' && op_ct == ct) {
		if (class[2] != 'T') {
			SEXP p1 = PROTECT(allocZero(INTSXP, (R_xlen_t) n + 1));
			SET_SLOT(to, Matrix_pSym, p1);
			UNPROTECT(1); /* p1 */
		}
		UNPROTECT(2); /* to, from */
		return to;
	}

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			j, k, kend, nnz0 = pp0[n], nnz1;
		pp0++;

		if (class[1] == 's') {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			/* Skew-symmetric part of Hermitian matrix is imaginary part */
			/* Skew-Hermitian part of symmetric matrix is imaginary part */
			SEXP x1 = PROTECT(allocVector(CPLXSXP, nnz0));
			zvimag(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else {

			int *iwork = NULL;
			size_t liwork = (size_t) ((int_fast64_t) n + n + 1 + nnz0);
			Matrix_Calloc(iwork, liwork, int);

			int *pp0_ = iwork + n + 1, *pi0_ = iwork + n + n + 1;
			nsptrans(pp0_ - 1, pi0_, NULL, pp0 - 1, pi0, NULL, n, n, 'T', iwork);
			memcpy(iwork, pp0 - 1, sizeof(int) * (size_t) n);

			SEXP p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0)));
			int *pp1 = INTEGER(p1), k_, kend_;
			*(pp1++) = 0;
			SET_SLOT(to, Matrix_pSym, p1);

			for (j = 0, k = 0, k_ = 0; j < n; ++j) {
				kend  = pp0 [j];
				kend_ = pp0_[j];
				pp1[j] = 0;
				while (k < kend) {
					if (pi0[k] >= j)
						k = kend;
					else {
						while (k_ < kend_ && pi0_[k_] < pi0[k]) {
							++pp1[j];
							++pp1[pi0_[k_]];
							++k_;
						}
						++pp1[j];
						++pp1[pi0[k]];
						if (k_ < kend_ && pi0_[k_] == pi0[k])
							++k_;
						++k;
					}
				}
				while (k_ < kend_) {
					if (pi0_[k_] >= j)
						k_ = kend_;
					else {
						++pp1[j];
						++pp1[pi0_[k_]];
						++k_;
					}
				}
			}
			for (j = 0; j < n; ++j)
				pp1[j] += pp1[j - 1];
			nnz1 = pp1[n - 1];
			for (j = n - 1; j >= 0; --j)
				pp1[j] = pp1[j - 1];

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), l;
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				if (op_ct == 'C') \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (pi0[k] > j) \
							k = kend; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								c##ASSIGN_CONJ(px1[pp1[       j]], px0[l]); \
								c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
								c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
								c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
								pi1[pp1[       j]++] = pi0_[k_]; \
								pi1[pp1[pi0_[k_]]++] =        j; \
								++k_; \
							} \
							l = iwork[j]++; \
							if (pi0[k] == j) \
								++k_; \
							else { \
								c##ASSIGN_IDEN(px1[pp1[     j]], px0[k]); \
								c##ASSIGN_CONJ(px1[pp1[pi0[k]]], px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##DECREMENT_CONJ(px1[pp1[     j]], px0[l]); \
									c##DECREMENT_IDEN(px1[pp1[pi0[k]]], px0[l]); \
									++k_; \
								} \
								c##MULTIPLY   (px1[pp1[     j]],    0.5); \
								c##MULTIPLY   (px1[pp1[pi0[k]]],   -0.5); \
								pi1[pp1[     j]++] = pi0[k]; \
								pi1[pp1[pi0[k]]++] =      j; \
							} \
							++k; \
						} \
					} \
					while (k_ < kend_) { \
						if (pi0_[k_] >= j) \
							k_ = kend_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							c##ASSIGN_CONJ(px1[pp1[       j]], px0[l]); \
							c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
							c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
							c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
							pi1[pp1[       j]++] = pi0_[k_]; \
							pi1[pp1[pi0_[k_]]++] =        j; \
							++k_; \
						} \
					} \
				} \
				else \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (pi0[k] > j) \
							k = kend; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								c##ASSIGN_IDEN(px1[pp1[       j]], px0[l]); \
								c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
								c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
								c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
								pi1[pp1[       j]++] = pi0_[k_]; \
								pi1[pp1[pi0_[k_]]++] =        j; \
								++k_; \
							} \
							l = iwork[j]++; \
							if (pi0[k] == j) \
								++k_; \
							else { \
								c##ASSIGN_IDEN(px1[pp1[     j]], px0[k]); \
								c##ASSIGN_IDEN(px1[pp1[pi0[k]]], px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##DECREMENT_IDEN(px1[pp1[     j]], px0[l]); \
									c##DECREMENT_IDEN(px1[pp1[pi0[k]]], px0[l]); \
									++k_; \
								} \
								c##MULTIPLY   (px1[pp1[     j]],    0.5); \
								c##MULTIPLY   (px1[pp1[pi0[k]]],   -0.5); \
								pi1[pp1[     j]++] = pi0[k]; \
								pi1[pp1[pi0[k]]++] =      j; \
							} \
							++k; \
						} \
					} \
					while (k_ < kend_) { \
						if (pi0_[k_] >= j) \
							k_ = kend_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							c##ASSIGN_IDEN(px1[pp1[       j]], px0[l]); \
							c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
							c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
							c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
							pi1[pp1[       j]++] = pi0_[k_]; \
							pi1[pp1[pi0_[k_]]++] =        j; \
							++k_; \
						} \
					} \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			Matrix_Free(iwork, liwork);
			UNPROTECT(3); /* x1, i1, p1 */

		}

		UNPROTECT(3); /* x0, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1;

		if (class[1] == 's') {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			/* Skew-symmetric part of Hermitian matrix is imaginary part */
			/* Skew-Hermitian part of symmetric matrix is imaginary part */
			SEXP x1 = PROTECT(allocVector(CPLXSXP, nnz0));
			zvimag(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else {

			nnz1 = nnz0;
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					--nnz1;
			nnz1 *= 2;

			SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
				j1 = PROTECT(allocVector(INTSXP, nnz1)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);

#define SPART(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						c##ASSIGN_IDEN(*px1, *px0); \
						c##MULTIPLY(*px1,  0.5); \
						++pi1; ++pj1; ++px1; \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
						if (op_ct == 'C') \
						c##ASSIGN_CONJ(*px1, *px0); \
						else \
						c##ASSIGN_IDEN(*px1, *px0); \
						c##MULTIPLY(*px1, -0.5); \
						++pi1; ++pj1; ++px1; \
					} \
					++pi0; ++pj0; ++px0; \
				} \
			} while (0)

			SWITCH2(cl[0], SPART);

#undef SPART

			UNPROTECT(3); /* x1, j1, i1 */

		}

		UNPROTECT(3); /* x0, j0, i0 */

	}

	UNPROTECT(2); /* to, from */
	return to;
}

/* skewpart(<[CRT]sparseMatrix>, trans) */
SEXP R_sparse_skewpart(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	return sparse_skewpart(s_from, class, ct);
}

int sparse_is_symmetric(SEXP obj, const char *class,
                        char op_ct, int exact, int checkDN)
{
	if (class[0] != 'z') {
		op_ct = '\0';
		if (class[0] != 'd')
			exact = 1;
	}

	char ct = '\0';
	if (class[1] == 's') {
		if (class[0] != 'z')
			return 1;
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (op_ct == ct)
			return 1;
		checkDN = 0;
	}

	if (checkDN) {
		SEXP dimnames = GET_SLOT(obj, Matrix_DimNamesSym);
		if (!DimNames_is_symmetric(dimnames))
			return 0;
	}

	char di = '\0';
	if (class[1] == 't') {
		if (!sparse_is_diagonal(obj, class))
			return 0;
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N' || op_ct != 'C')
			return 1;
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n == 0 || (n == 1 && op_ct != 'C'))
		return 1;
	if (!exact && class[0] != 'g')
		return NA_LOGICAL; /* do inexact numerical test in R */

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		if (exact)
		obj = sparse_as_Csparse(obj, class);
		else {
		/* Testing for symmetric nonzero pattern, hence: */
		char cl[] = "n.TMatrix";
		cl[1] = class[1];
		obj = sparse_as_Csparse(obj, cl);
		}
	}
	PROTECT(obj);

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p_ = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i_ = PROTECT(GET_SLOT(obj,        iSym));
	int *pp = INTEGER(p_), *pi = INTEGER(i_), i, j, k, kend, *iwork = NULL;
	if (class[1] == 'g') {
	Matrix_Calloc(iwork, n, int);
	memcpy(iwork, pp, sizeof(int) * (size_t) n);
	}
	pp++;

	int ans = 0;

	if (class[1] == 'g' && !(exact && op_ct == 'C')) {

#define ISS(c) \
	do { \
		c##IF_NPATTERN( \
		SEXP x = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px = c##PTR(x); \
		); \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				i = pi[k]; \
				if (i < j) { \
					if (iwork[i] == pp[i] || pi[iwork[i]] != j) \
						goto done; \
					c##IF_NPATTERN( \
					if (c##NOT_IDEN(px[k], px[iwork[i]])) \
						goto done; \
					); \
					++iwork[i]; \
					++iwork[j]; \
					++k; \
				} else { \
					if (pi[k] == j) \
						++iwork[j]; \
					k = kend; \
				} \
			} \
		} \
	} while (0)

	SWITCH5((exact) ? class[0] : 'n', ISS);

#undef ISS

	} else {

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	Rcomplex *px = COMPLEX(x);
	if (class[1] == 'g')
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			while (k < kend) {
				i = pi[k];
				if (i < j) {
					if (iwork[i] == pp[i] || pi[iwork[i]] != j) \
						goto done;
					if (zNOT_CONJ(px[k], px[iwork[i]]))
						goto done;
					++iwork[i];
					++iwork[j];
					++k;
				} else {
					if (i == j) {
						if (zNOT_ZERO_IMAG(px[k]))
							goto done;
						++iwork[j];
					}
					k = kend;
				}
			}
		}
	else if (class[1] == 's')
		/* Testing if Hermitian matrix is symmetric */
		/*      or if symmetric matrix is Hermitian */
		/* <=====> if matrix is real                */
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			while (k < kend) {
				if (zNOT_ZERO_IMAG(px[k]) && (pi[k] != j || op_ct == 'C'))
				    goto done;
				++k;
			}
		}
	else
		/* Testing if non-unit triangular matrix is Hermitian */
		/* <=====> if matrix is real and diagonal             */
		/* NB: diagonal nonzero pattern is already known      */
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && zNOT_ZERO_IMAG(px[k]))
				goto done;
			k = kend;
		}

	}

	if (class[1] == 'g')
	for (j = 0; j < n; ++j)
		if (iwork[j] != pp[j])
			goto done;

	ans = (exact) ? 1 : NA_LOGICAL;

done:
	if (class[1] == 'g')
	Matrix_Free(iwork, n);
	UNPROTECT(3); /* i_, p_, obj */
	return ans;
}

/* isSymmetric(<[CRT]sparseMatrix>, tol, tol1, trans, checkDN)
   NB: requires symmetric nonzero pattern
*/
SEXP R_sparse_is_symmetric(SEXP s_obj,
                           SEXP s_trans, SEXP s_exact, SEXP s_checkDN)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	int exact, checkDN;
	VALID_LOGIC2(s_exact  , exact  );
	VALID_LOGIC2(s_checkDN, checkDN);

	int ans_ = sparse_is_symmetric(s_obj, class, ct, exact, checkDN);
	SEXP ans = ScalarLogical(ans_);
	return ans;
}

int sparse_is_triangular(SEXP obj, const char *class, char op_ul)
{
	if (class[1] == 't') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (op_ul == '\0' || op_ul == ul)
			return (   ul == 'U') ? 1 : -1;
		else if (sparse_is_diagonal(obj, class))
			return (op_ul == 'U') ? 1 : -1;
		else
			return 0;
	}

	if (class[1] == 's') {
		if (!sparse_is_diagonal(obj, class))
			return 0;
		else if (op_ul != '\0')
			return (op_ul == 'U') ? 1 : -1;
		else {
			SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
			char ul = CHAR(STRING_ELT(uplo, 0))[0];
			return (   ul == 'U') ? 1 : -1;
		}
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return ((op_ul == '\0') ? class[2] != 'R' : op_ul == 'U') ? 1 : -1;

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		pp++;
		UNPROTECT(2); /* i, p */

#define IST_CURL(t) \
		do { \
			for (j = 0, k = 0; j < n; ++j) { \
				kend = pp[j]; \
				if (k < kend && pi[kend - 1] > j) \
					break; \
				k = kend; \
			} \
			if (j == n) \
				return t; \
		} while (0)

#define IST_CLRU(t) \
		do { \
			for (j = 0, k = 0; j < n; ++j) { \
				kend = pp[j]; \
				if (k < kend && pi[k] < j) \
					break; \
				k = kend; \
			} \
			if (j == n) \
				return t; \
		} while (0)

		if (op_ul == '\0' || op_ul == 'U') {
			if (class[2] == 'C')
			IST_CURL(1);
			else
			IST_CLRU(1);
		}
		if (op_ul == '\0' || op_ul != 'U') {
			if (class[2] == 'C')
			IST_CLRU(-1);
			else
			IST_CURL(-1);
		}

#undef IST_CURL
#undef IST_CLRU

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);
		UNPROTECT(2); /* i, j */

		if (op_ul == '\0' || op_ul == 'U') {
			for (k = 0; k < kend; ++k)
				if (pi[k] > pj[k])
					break;
			if (k == kend)
				return  1;
		}
		if (op_ul == '\0' || op_ul != 'U') {
			for (k = 0; k < kend; ++k)
				if (pi[k] < pj[k])
					break;
			if (k == kend)
				return -1;
		}

	}

	return 0;
}

/* isTriangular(<[CRT]sparseMatrix>, upper)
   NB: requires triangular nonzero pattern
*/
SEXP R_sparse_is_triangular(SEXP s_obj, SEXP s_upper)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int upper;
	VALID_LOGIC3(s_upper, upper);

	int ans_ = sparse_is_triangular(s_obj, class,
		(upper == NA_LOGICAL) ? '\0' : ((upper != 0) ? 'U' : 'L'));
	SEXP ans = allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = ans_ != 0;
	if (upper == NA_LOGICAL && ans_ != 0) {
		PROTECT(ans);
		static
		SEXP kindSym = NULL;
		SEXP kindVal = PROTECT(mkString((ans_ > 0) ? "U" : "L"));
		if (!kindSym)
			kindSym = install("kind");
		setAttrib(ans, kindSym, kindVal);
		UNPROTECT(2); /* kindVal, ans */
	}
	return ans;
}

int sparse_is_diagonal(SEXP obj, const char *class)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return 1;

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		pp++;
		UNPROTECT(2); /* i, p */

		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && (kend - k > 1 || pi[k] != j))
				return 0;
			k = kend;
		}
		return 1;

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);
		UNPROTECT(2); /* j, i */

		for (k = 0; k < kend; ++k)
			if (pi[k] != pj[k])
				return 0;
		return 1;

	}
}

/* isDiagonal(<[CRT]sparseMatrix>)
   NB: requires diagonal nonzero pattern
*/
SEXP R_sparse_is_diagonal(SEXP s_obj)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);
	return ScalarLogical(sparse_is_diagonal(s_obj, class));
}

#define MAP(i) (map) ? map[i] : i

#define nCAST(x) (1)
#define lCAST(x) (x != 0)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void Csparse_colsum(SEXP obj, const char *class,
                    int m, int n, char ul, char ct, char di, int narm, int mean,
                    SEXP ans)
{
	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp0 = INTEGER(p0), j, k, kend, nnz1, count = -1,
		un = class[1] == 't' && di != 'N';
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

		SEXP j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pj1 = INTEGER(j1);
		SET_SLOT(ans, Matrix_iSym, j1);
		if (un)
			for (j = 0; j < n; ++j)
				*(pj1++) = j + 1;
		else
			for (j = 0; j < n; ++j)
				if (pp0[j - 1] < pp0[j])
					*(pj1++) = j + 1;

		PROTECT(x1 = allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(ans, Matrix_xSym, x1);

		UNPROTECT(2); /* x1, j1 */
	}
	PROTECT(x1);

	int full = nnz1 == n;

#define SUM(c, d) \
	do { \
		c##IF_NPATTERN( \
		SEXP x0 = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px0 = c##PTR(x0); \
		); \
		d##TYPE *px1 = d##PTR(x1), tmp1 = (un) ? d##UNIT : d##ZERO; \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			if (full || k < kend) { \
				*px1 = tmp1; \
				if (mean) \
					count = m; \
				while (k < kend) { \
					if (c##NOT_NA(*px0)) \
						d##INCREMENT_IDEN(*px1, c##CAST(*px0)); \
					else if (!narm) \
						*px1 = d##NA; \
					else if (mean) \
						--count; \
					c##IF_NPATTERN( \
					++px0; \
					); \
					++k; \
				} \
				if (mean) \
					d##DIVIDE(*px1, count); \
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
                    int m, int n, char ul, char ct, char di, int narm, int mean,
                    SEXP ans, SEXP iSym)
{
	SEXP p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(obj,        iSym));
	int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), i, j, k, kend, nnz1,
		*count = NULL, *map = NULL,
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && di != 'N';
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

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1));
		SET_SLOT(ans, Matrix_iSym, i1);
		if (narm && mean) {
			count = INTEGER(i1);
			for (i = 0; i < nnz1; ++i)
				count[i] = n;
		}

		PROTECT(x1 = allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(ans, Matrix_xSym, x1);

		UNPROTECT(2); /* x1, i1 */
	}
	PROTECT(x1);

#define SUM(c, d) \
	do { \
		c##IF_NPATTERN( \
		SEXP x0 = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px0 = c##PTR(x0); \
		); \
		c##TYPE                    tmp0; \
		d##TYPE *px1 = d##PTR(x1), tmp1 = (un) ? d##UNIT : d##ZERO; \
		for (i = 0; i < nnz1; ++i) \
			px1[i] = tmp1; \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			while (k < kend) { \
				i = pi0[k]; \
				if (he && i == j) \
				c##ASSIGN_PROJ_REAL(tmp0, c##IFELSE_NPATTERN(px0[k], c##UNIT)); \
				else \
				c##ASSIGN_IDEN     (tmp0, c##IFELSE_NPATTERN(px0[k], c##UNIT)); \
				if (c##NOT_NA(tmp0)) { \
					d##INCREMENT_IDEN(px1[MAP(i)], c##CAST(tmp0)); \
					if (sy && i != j) { \
					if (he) \
					d##INCREMENT_CONJ(px1[MAP(j)], c##CAST(tmp0)); \
					else \
					d##INCREMENT_IDEN(px1[MAP(j)], c##CAST(tmp0)); \
					} \
				} \
				else if (!narm) { \
					px1[MAP(i)] = d##NA; \
					if (sy && i != j) \
					px1[MAP(j)] = d##NA; \
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
					d##DIVIDE(px1[i], n); \
			else \
				for (i = 0; i < nnz1; ++i) \
					d##DIVIDE(px1[i], count[i]); \
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
                    int m, int n, char ul, char ct, char di, int narm, int mean,
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
		un = class[1] == 't' && di != 'N';;
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

		SEXP j1 = PROTECT(allocVector(INTSXP, nnz1));
		SET_SLOT(ans, Matrix_iSym, j1);
		if (narm && mean) {
			count = INTEGER(j1);
			for (j = 0; j < nnz1; ++j)
				count[j] = m;
		}

		PROTECT(x1 = allocVector(SUM_TYPEOF(class[0]), nnz1));
		SET_SLOT(ans, Matrix_xSym, x1);

		UNPROTECT(2); /* x1, j1 */
	}
	PROTECT(x1);

#define SUM(c, d) \
	do { \
		c##IF_NPATTERN( \
		SEXP x0 = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px0 = c##PTR(x0); \
		); \
		c##TYPE                    tmp0; \
		d##TYPE *px1 = d##PTR(x1), tmp1 = (un) ? d##UNIT : d##ZERO; \
		for (j = 0; j < nnz1; ++j) \
			px1[j] = tmp1; \
		for (k = 0; k < kend; ++k) { \
			i = pi0[k]; \
			j = pj0[k]; \
			if (he && i == j) \
			c##ASSIGN_PROJ_REAL(tmp0, c##IFELSE_NPATTERN(px0[k], c##UNIT)); \
			else \
			c##ASSIGN_IDEN     (tmp0, c##IFELSE_NPATTERN(px0[k], c##UNIT)); \
			if (c##NOT_NA(tmp0)) { \
				d##INCREMENT_IDEN(px1[MAP(j)], c##CAST(tmp0)); \
				if (sy && i != j) { \
				if (he) \
				d##INCREMENT_CONJ(px1[MAP(i)], c##CAST(tmp0)); \
				else \
				d##INCREMENT_IDEN(px1[MAP(i)], c##CAST(tmp0)); \
				} \
			} \
			else if (!narm) { \
				px1[MAP(j)] = d##NA; \
				if (sy && i != j) \
				px1[MAP(i)] = d##NA; \
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
					d##DIVIDE(px1[j], m); \
			else \
				for (j = 0; j < nnz1; ++j) \
					d##DIVIDE(px1[j], count[j]); \
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

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (mg == 0) ? m : n;

	SEXP ans;
	SEXPTYPE type = SUM_TYPEOF(class[0]);
	if (sparse) {
		char cl[] = ".sparseVector";
		cl[0] = typeToKind(type);
		PROTECT(ans = newObject(cl));

		SEXP length = GET_SLOT(ans, Matrix_lengthSym);
		INTEGER(length)[0] = r;
	} else {
		PROTECT(ans = allocVector(type, r));

		SEXP dimnames = (class[1] == 's')
			? get_symmetrized_DimNames(obj, -1)
			: GET_SLOT(obj, Matrix_DimNamesSym),
			marnames = VECTOR_ELT(dimnames, mg);
		if (marnames != R_NilValue) {
			PROTECT(marnames);
			setAttrib(ans, R_NamesSymbol, marnames);
			UNPROTECT(1); /* marnames */
		}
	}

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	if (mg == 0)
		switch (class[2]) {
		case 'C':
			Csparse_rowsum(obj, class, m, n, ul, ct, di, narm, mean, ans,
			               Matrix_iSym);
			break;
		case 'R':
			if (class[1] == 's')
			Csparse_rowsum(obj, class, n, m, ul, ct, di, narm, mean, ans,
			               Matrix_jSym);
			else
			Csparse_colsum(obj, class, n, m, ul, ct, di, narm, mean, ans);
			break;
		case 'T':
			Tsparse_colsum(obj, class, n, m, ul, ct, di, narm, mean, ans,
			               Matrix_jSym, Matrix_iSym);
			break;
		default:
			break;
		}
	else
		switch (class[2]) {
		case 'C':
			if (class[1] == 's')
			Csparse_rowsum(obj, class, m, n, ul, ct, di, narm, mean, ans,
			               Matrix_iSym);
			else
			Csparse_colsum(obj, class, m, n, ul, ct, di, narm, mean, ans);
			break;
		case 'R':
			Csparse_rowsum(obj, class, n, m, ul, ct, di, narm, mean, ans,
			               Matrix_jSym);
			break;
		case 'T':
			Tsparse_colsum(obj, class, m, n, ul, ct, di, narm, mean, ans,
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

#undef SUM_CASES

#define TRY_INCREMENT(_LABEL_) \
	do { \
		if ((s >= 0) \
			? ( t <= INT_FAST64_MAX - s) \
			: (-t <= s - INT_FAST64_MIN)) { \
			s += t; \
			t = 0; \
			count = 0; \
		} else { \
			over = 1; \
			goto _LABEL_; \
		} \
	} while (0)

SEXP sparse_sum(SEXP obj, const char *class, int narm)
{
	PROTECT(obj = sparse_aggregate(obj, class));

	SEXP ans;

	if (!narm && (class[0] == 'l' || class[0] == 'i')) {
		SEXP x = GET_SLOT(obj, Matrix_xSym);
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		R_xlen_t nx = XLENGTH(x);
		while (nx--) {
			if (*px == NA_INTEGER) {
				ans = allocVector(INTSXP, 1);
				INTEGER(ans)[0] = NA_INTEGER;
				UNPROTECT(1); /* obj */
				return ans;
			}
			++px;
		}
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char di = 'N';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	int symmetric = class[1] == 's';

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j_, k, kend,
			n_ = (class[2] == 'C') ? n : m;

		if (class[0] == 'n') {
			int_fast64_t nnz = pp[n_ - 1];
			if (di != 'N')
				nnz += n;
			if (symmetric) {
				SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
				char ul = CHAR(STRING_ELT(uplo, 0))[0];

				nnz *= 2;
				for (j_ = 0, k = 0; j_ < n_; ++j_) {
					kend = pp[j_];
					if (k < kend && pi[(ul == 'U') ? kend - 1 : k] == j_)
						--nnz;
					k = kend;
				}
			}
			if (nnz <= INT_MAX) {
				ans = allocVector(INTSXP, 1);
				INTEGER(ans)[0] = (int) nnz;
			} else {
				ans = allocVector(REALSXP, 1);
				REAL(ans)[0] = (double) nnz;
			}
			UNPROTECT(3); /* i, p, obj */
			return ans;
		}

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* i, p */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr = (di == 'N') ? 0.0L : n, zi = 0.0L;
			for (j_ = 0, k = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				while (k < kend) {
					if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
						zr += (symmetric && pi[k] != j_)
							? 2.0L * px[k].r : px[k].r;
						zi += (symmetric && pi[k] != j_)
							? 2.0L * px[k].i : px[k].i;
					}
					++k;
				}
			}
			ans = allocVector(CPLXSXP, 1);
			COMPLEX(ans)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
			COMPLEX(ans)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			long double zr = (di == 'N') ? 0.0L : n;
			for (j_ = 0, k = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				while (k < kend) {
					if (!(narm && ISNAN(px[k])))
						zr += (symmetric && pi[k] != j_)
							? 2.0L * px[k] : px[k];
					++k;
				}
			}
			ans = allocVector(REALSXP, 1);
			REAL(ans)[0] = LONGDOUBLE_AS_DOUBLE(zr);
		} else {
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			int_fast64_t s = (di == 'N') ? 0LL : n, t = 0LL;
			unsigned int count = 0;
			int over = 0;
			for (j_ = 0, k = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				while (k < kend) {
					if (!narm || px[k] != NA_INTEGER) {
						unsigned int d = (symmetric && pi[k] != j_) ? 2 : 1;
						if (count > UINT_MAX - d)
							TRY_INCREMENT(ifoverC);
						t += (d == 2) ? 2LL * px[k] : px[k];
						count += d;
					}
					++k;
				}
			}
			TRY_INCREMENT(ifoverC);
		ifoverC:
			if (over) {
				long double zr = (long double) s + (long double) t;
				for (; j_ < n_; ++j_) {
					kend = pp[j_];
					while (k < kend) {
						if (!narm || px[k] != NA_INTEGER)
							zr += (symmetric && pi[k] != j_)
								? 2.0L * px[k] : px[k];
						++k;
					}
				}
				ans = allocVector(REALSXP, 1);
				REAL(ans)[0] = LONGDOUBLE_AS_DOUBLE(zr);
			} else if (s > INT_MIN && s <= INT_MAX) {
				ans = allocVector(INTSXP, 1);
				INTEGER(ans)[0] = (int) s;
			} else {
				ans = allocVector(REALSXP, 1);
				REAL(ans)[0] = (double) s;
			}
		}

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

		if (class[0] == 'n') {
			int_fast64_t nnz = (int_fast64_t) kend;
			if (di != 'N')
				nnz += n;
			if (symmetric) {
				nnz *= 2;
				for (k = 0; k < kend; ++k)
					if (pi[k] == pj[k])
						--nnz;
			}
			if (nnz <= INT_MAX) {
				ans = allocVector(INTSXP, 1);
				INTEGER(ans)[0] = (int) nnz;
			} else {
				ans = allocVector(REALSXP, 1);
				REAL(ans)[0] = (double) nnz;
			}
			UNPROTECT(3); /* j, i, obj */
			return ans;
		}

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* j, i */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr = (di == 'N') ? 0.0L : n, zi = 0.0L;
			for (k = 0; k < kend; ++k)
				if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
					zr += (symmetric && pi[k] != pj[k])
						? 2.0L * px[k].r : px[k].r;
					zi += (symmetric && pi[k] != pj[k])
						? 2.0L * px[k].i : px[k].i;
				}
			ans = allocVector(CPLXSXP, 1);
			COMPLEX(ans)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
			COMPLEX(ans)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			long double zr = (di == 'N') ? 0.0L : n;
			for (k = 0; k < kend; ++k)
				if (!(narm && ISNAN(px[k])))
					zr += (symmetric && pi[k] != pj[k])
						? 2.0L * px[k] : px[k];
			ans = allocVector(REALSXP, 1);
			REAL(ans)[0] = LONGDOUBLE_AS_DOUBLE(zr);
		} else {
			int *px = (class[0] == 'i') ? INTEGER(x) : LOGICAL(x);
			int_fast64_t s = (di == 'N') ? 0LL : n, t = 0LL;
			unsigned int count = 0;
			int over = 0;
			for (k = 0; k < kend; ++k) {
				if (!narm || px[k] != NA_INTEGER) {
					unsigned int d = (symmetric && pi[k] != pj[k]) ? 2 : 1;
					if (count > UINT_MAX - d)
						TRY_INCREMENT(ifoverT);
					t += (d == 2) ? 2LL * px[k] : px[k];
					count += d;
				}
			}
			TRY_INCREMENT(ifoverT);
		ifoverT:
			if (over) {
				long double zr = (long double) s + (long double) t;
				for (; k < kend; ++k)
					if (!(narm && px[k] == NA_INTEGER))
						zr += (symmetric && pi[k] != pj[k])
							? 2.0L * px[k] : px[k];
				ans = allocVector(REALSXP, 1);
				REAL(ans)[0] = LONGDOUBLE_AS_DOUBLE(zr);
			} else if (s > INT_MIN && s <= INT_MAX) {
				ans = allocVector(INTSXP, 1);
				INTEGER(ans)[0] = (int) s;
			} else {
				ans = allocVector(REALSXP, 1);
				REAL(ans)[0] = (double) s;
			}
		}

	}

	UNPROTECT(1); /* obj */
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

	SEXP ans = PROTECT(allocVector((class[0] == 'z') ? CPLXSXP : REALSXP, 1));

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}

	int symmetric = (class[1] != 's')
		? 0 : (((class[2] == 'C') == (ul == 'U')) ? 1 : -1);
	long double zr = 1.0L, zi = 0.0L;

	int_fast64_t mn = (int_fast64_t) m * n,
		nnz, nnzmax = (symmetric) ? (mn + n) / 2 : mn;

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p) + 1, *pi = INTEGER(i), i_, j_, k, kend,
			m_ = (class[2] == 'C') ? m : n, n_ = (class[2] == 'C') ? n : m,
			seen0 = 0;

		nnz = pp[n_ - 1];
		if (di != 'N')
			nnz += n;
		if (class[0] == 'n') {
			REAL(ans)[0] = (nnz == nnzmax) ? 1.0 : 0.0;
			UNPROTECT(4); /* i, p, ans, obj */
			return ans;
		}

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* i, p */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr0, zi0;
			for (j_ = 0, k = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				if (seen0 || kend - k == m_) {
					while (k < kend) {
						if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							if (symmetric && pi[k] != j_) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							}
						}
						++k;
					}
				} else {
					int i0 = (symmetric >= 0) ? 0 : j_,
						i1 = (symmetric <= 0) ? m_ : j_ + 1;
					for (i_ = i0; i_ < i1; ++i_) {
						if (seen0 || (k < kend && pi[k] == i_)) {
						if (k >= kend)
							break;
						if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							if (symmetric && pi[k] != j_) {
							zr0 = zr; zi0 = zi;
							zr = zr0 * px[k].r - zi0 * px[k].i;
							zi = zr0 * px[k].i + zi0 * px[k].r;
							}
						}
						++k;
						} else if (di == 'N' || i_ != j_) {
						zr *= 0.0L;
						zi *= 0.0L;
						seen0 = 1;
						}
					}
				}
			}
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			for (j_ = 0, k = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				if (seen0 || kend - k == m_) {
					while (k < kend) {
						if (!(narm && ISNAN(px[k])))
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						++k;
					}
				} else {
					int i0 = (symmetric >= 0) ? 0 : j_,
						i1 = (symmetric <= 0) ? m_ : j_ + 1;
					for (i_ = i0; i_ < i1; ++i_) {
						if (seen0 || (k < kend && pi[k] == i_)) {
						if (k >= kend)
							break;
						if (!(narm && ISNAN(px[k])))
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						++k;
						} else if (di == 'N' || i_ != j_) {
						zr *= 0.0L;
						seen0 = 1;
						}
					}
				}
			}
		} else {
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			for (j_ = 0, k = 0; j_ < n_; ++j_) {
				kend = pp[j_];
				if (seen0 || kend - k == m_) {
					while (k < kend) {
						if (px[k] != NA_INTEGER)
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						else if (!narm)
							zr *= NA_REAL;
						++k;
					}
				} else {
					int i0 = (symmetric >= 0) ? 0 : j_,
						i1 = (symmetric <= 0) ? m_ : j_ + 1;
					for (i_ = i0; i_ < i1; ++i_) {
						if (seen0 || (k < kend && pi[k] == i_)) {
						if (k >= kend)
							break;
						if (px[k] != NA_INTEGER)
							zr *= (symmetric && pi[k] != j_)
								? (long double) px[k] * px[k] : px[k];
						else if (!narm)
							zr *= NA_REAL;
						++k;
						} else if (di == 'N' || i_ != j_) {
						zr *= 0.0L;
						seen0 = 1;
						}
					}
				}
			}
		}

	} else {

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

		nnz = (int_fast64_t) kend;
		if (di != 'N')
			nnz += n;
		if (class[0] == 'n') {
			REAL(ans)[0] = (nnz == nnzmax) ? 1.0 : 0.0;
			UNPROTECT(4); /* j, i, ans, obj */
			return ans;
		}
		if (nnz < nnzmax)
			zr = 0.0;

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		UNPROTECT(2); /* j, i */

		if (class[0] == 'z') {
			Rcomplex *px = COMPLEX(x);
			long double zr0, zi0;
			for (k = 0; k < kend; ++k)
				if (!(narm && (ISNAN(px[k].r) || ISNAN(px[k].i)))) {
					zr0 = zr; zi0 = zi;
					zr = zr0 * px[k].r - zi0 * px[k].i;
					zi = zr0 * px[k].i + zi0 * px[k].r;
					if (symmetric && pi[k] != pj[k]) {
					zr0 = zr; zi0 = zi;
					zr = zr0 * px[k].r - zi0 * px[k].i;
					zi = zr0 * px[k].i + zi0 * px[k].r;
					}
				}
		} else if (class[0] == 'd') {
			double *px = REAL(x);
			for (k = 0; k < kend; ++k)
				if (!(narm && ISNAN(px[k])))
					zr *= (symmetric && pi[k] != pj[k])
						? (long double) px[k] * px[k] : px[k];
		} else {
			int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
			for (k = 0; k < kend; ++k)
				if (px[k] != NA_INTEGER)
					zr *= (symmetric && pi[k] != pj[k])
						? (long double) px[k] * px[k] : px[k];
				else if (!narm)
					zr *= NA_REAL;
		}

	}

	if (class[0] == 'z') {
		COMPLEX(ans)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
		COMPLEX(ans)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
	} else
		   REAL(ans)[0]   = LONGDOUBLE_AS_DOUBLE(zr);
	UNPROTECT(2); /* ans, obj */
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
