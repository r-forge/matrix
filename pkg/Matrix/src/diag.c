/* C implementation of methods for diag, diag<- */

#include <math.h> /* sqrt */
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

SEXP dense_diag_get(SEXP obj, const char *class, int names)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n, packed = class[2] == 'p';

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP ans;

	if (nu != '\0' && nu != 'N') {

		PROTECT(ans = allocUnit(kindToType(class[0]), r));

	} else {

		PROTECT(ans = Rf_allocVector(kindToType(class[0]), r));

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		size_t m_ = (size_t) m, n_ = (size_t) n, r_ = (size_t) r;

#define TEMPLATE(c) \
		do { \
			c##TYPE *pa = c##PTR(ans), *px = c##PTR(x); \
			if (!packed) \
				c##NAME(copy2)(r_, pa, 1, px, m_ + 1); \
			else if (ul == 'U') \
				c##NAME(copy1)(n_, pa, 1, 0, 0, px, 2 , 1, 0); \
			else \
				c##NAME(copy1)(n_, pa, 1, 0, 0, px, n_, 1, 1); \
		} while (0)

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		if (class[0] == 'n')
			naToUnit(ans);
		if (class[0] == 'z' && class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), NULL, r_);

	}

	if (names) {
		/* NB: Logic here must be adjusted once the validity method */
		/*     for 'symmetricMatrix' enforces symmetric 'Dimnames'  */
		SEXP dn = PROTECT(DIMNAMES(obj, 0)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (cn == R_NilValue) {
			if (class[1] == 's' && rn != R_NilValue)
				Rf_setAttrib(ans, R_NamesSymbol, rn);
		} else {
			if (class[1] == 's')
				Rf_setAttrib(ans, R_NamesSymbol, cn);
			else if (rn != R_NilValue &&
			         (rn == cn || equalString(rn, cn, r)))
				Rf_setAttrib(ans, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(1); /* ans */
	return ans;
}

SEXP sparse_diag_get(SEXP obj, const char *class, int names)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP ans;

	if (nu != '\0' && nu != 'N') {

		PROTECT(ans = allocUnit(kindToType(class[0]), r));

	} else if (class[2] != 'T') {

		PROTECT(ans = Rf_allocVector(kindToType(class[0]), r));

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj,        iSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend,
			up = (class[2] == 'C') == (ul == 'U');
		pp++;

#define TEMPLATE(c) \
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
					pa[j] = (k < kend && pi[(up) ? kend - 1 : k] == j) \
						? c##IFELSE_NPATTERN(px[(up) ? kend - 1 : k], c##UNIT) \
						: c##ZERO; \
					k = kend; \
				} \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(2); /* i, p */

		if (class[0] == 'z' && class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), NULL, (size_t) r);

	} else {

		PROTECT(ans = Rf_allocVector(kindToType(class[0]), r));

		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t k, kend = XLENGTH(i);

#define TEMPLATE(c) \
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

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(2); /* j0, i0 */

		if (class[0] == 'z' && class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), NULL, (size_t) r);

	}

	if (names) {
		/* NB: Logic here must be adjusted once the validity method */
		/*     for 'symmetricMatrix' enforces symmetric 'Dimnames'  */
		SEXP dn = PROTECT(DIMNAMES(obj, 0)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (cn == R_NilValue) {
			if (class[1] == 's' && rn != R_NilValue)
				Rf_setAttrib(ans, R_NamesSymbol, rn);
		} else {
			if (class[1] == 's')
				Rf_setAttrib(ans, R_NamesSymbol, cn);
			else if (rn != R_NilValue &&
			         (rn == cn || equalString(rn, cn, r)))
				Rf_setAttrib(ans, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(1); /* ans */
	return ans;
}

SEXP R_dense_diag_get(SEXP s_obj, SEXP s_names)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int names;
	VALID_LOGIC2(s_names, names);

	return dense_diag_get(s_obj, class, names);
}

SEXP R_sparse_diag_get(SEXP s_obj, SEXP s_names)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int names;
	VALID_LOGIC2(s_names, names);

	return sparse_diag_get(s_obj, class, names);
}

SEXP dense_diag_set(SEXP from, const char *class, SEXP value, int new)
{
	SEXP to = PROTECT(newObject(class));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n, packed = class[2] == 'p';
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	char ul = '\0', ct = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (new) {
		x = duplicateVector(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
	}

	size_t m_ = (size_t) m, n_ = (size_t) n, r_ = (size_t) r,
		d_ = (LENGTH(value) == r) ? 1 : 0;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x), *pv = c##PTR(value); \
		if (!packed) \
			c##NAME(copy2)(r_, px, m_ + 1, pv, d_); \
		else if (ul == 'U') \
			c##NAME(copy1)(n_, px, 2 , 1, 0, pv, d_, 0, 0); \
		else \
			c##NAME(copy1)(n_, px, n_, 1, 1, pv, d_, 0, 0); \
	} while (0)

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(2); /* x, to */
	return to;
}

SEXP sparse_diag_set(SEXP from, const char *class, SEXP value)
{
	SEXP to = PROTECT(newObject(class));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && (nu = DIAG(from)) != 'N')
		;

	int mode = 0;

	if (class[2] != 'T') {

		if (class[2] == 'R')
			SWAP(m, n, int, );

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0)));
		int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0),
			j, k, kend, nd0 = 0, nd1 = 0,
			up = (class[2] == 'C') == (ul == 'U');
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
		else if (class[1] == 's' || nu == 'N')
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[(up) ? kend - 1 : k] == j)
					++nd0;
				k = kend;
				pp1[j] = kend - nd0;
			}
		else
			for (j = 0; j < n; ++j)
				pp1[j] = pp0[j];

#define TEMPLATE(c) \
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

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		if (nd1 - nd0 > INT_MAX - pp0[n - 1])
			Rf_error(_("%s cannot exceed %s"), "p[length(p)]", "2^31-1");

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, pp1[n - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define TEMPLATE(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, pp1[n - 1])); \
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

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(4); /* i1, p1, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), j;
		R_xlen_t k, nd0 = 0, nd1 = 0, nnz0 = XLENGTH(i0), nnz1 = nnz0;

		if (nu == '\0' || nu == 'N')
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					++nd0;

#define TEMPLATE(c) \
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

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		if (nd1 - nd0 > R_XLEN_T_MAX - nnz0)
			Rf_error(_("%s cannot exceed %s"), "length(i)", "R_XLEN_T_MAX");
		nnz1 += nd1 - nd0;

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
			j1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#define TEMPLATE(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz1)); \
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

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

SEXP R_dense_diag_set(SEXP s_from, SEXP s_value)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);
	SEXPTYPE tx = kindToType(class[0]), tv = TYPEOF(s_value);

	switch (tv) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		Rf_error(_("replacement diagonal has incompatible type \"%s\""),
		         Rf_type2char(tv));
		break;
	}

	int *pdim = DIM(s_from), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;
	if (XLENGTH(s_value) != 1 && XLENGTH(s_value) != r)
		Rf_error(_("replacement diagonal has wrong length"));

	int new = 1;
	if (tv <= tx) {
		/* defined in ./coerce.c : */
		SEXP dense_as_general(SEXP, const char *, int);
		if (class[1] == 's' && class[0] == 'z' && tv == tx &&
		    TRANS(s_from) == 'C') {
			PROTECT(s_from = dense_as_general(s_from, class, 1));
			class = Matrix_class(s_from, valid_dense, 6, __func__);
			UNPROTECT(1);
			new = 0;
		}
		PROTECT(s_from);
		PROTECT(s_value = Rf_coerceVector(s_value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
		PROTECT(s_from = dense_as_kind(s_from, class, typeToKind(tv), 0));
		PROTECT(s_value);
		class = Matrix_class(s_from, valid_dense, 6, __func__);
		new = 0;
	}

	s_from = dense_diag_set(s_from, class, s_value, new);
	UNPROTECT(2);
	return s_from;
}

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
		Rf_error(_("replacement diagonal has incompatible type \"%s\""),
		         Rf_type2char(tv));
		break;
	}

	int *pdim = DIM(s_from), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;
	if (XLENGTH(s_value) != 1 && XLENGTH(s_value) != r)
		Rf_error(_("replacement diagonal has wrong length"));

	if (tv <= tx) {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		if (class[1] == 's' && class[0] == 'z' && tv == tx &&
		    TRANS(s_from) == 'C') {
			PROTECT(s_from = sparse_as_general(s_from, class));
			class = Matrix_class(s_from, valid_sparse, 6, __func__);
			UNPROTECT(1);
		}
		PROTECT(s_from);
		PROTECT(s_value = Rf_coerceVector(s_value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
		PROTECT(s_from = sparse_as_kind(s_from, class, typeToKind(tv)));
		PROTECT(s_value);
		class = Matrix_class(s_from, valid_sparse, 6, __func__);
	}

	s_from = sparse_diag_set(s_from, class, s_value);
	UNPROTECT(2);
	return s_from;
}

SEXP sparse_diag_U2N(SEXP from, const char *class)
{
	if (class[1] != 't' || DIAG(from) == 'N')
		return from;
	SEXP value = PROTECT(Rf_ScalarLogical(1)),
		to = R_sparse_diag_set(from, value);
	UNPROTECT(1); /* value */
	return to;
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

SEXP R_sparse_diag_U2N(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_diag_U2N(s_from, class);
}

SEXP R_sparse_diag_N2U(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_diag_N2U(s_from, class);
}

SEXP denseCholesky_diag_get(SEXP s_trf, SEXP s_root)
{
	SEXP x = PROTECT(GET_SLOT(s_trf, Matrix_xSym));
	int n = DIM(s_trf)[1];
	char ul = UPLO(s_trf);
	SEXP ans = Rf_allocVector(TYPEOF(x), n);
	if (n > 0) {
	int j, root = Rf_asLogical(s_root),
		packed = XLENGTH(x) != (int_fast64_t) n * n;
	R_xlen_t dx = (!packed) ? (R_xlen_t) n + 1 : (ul == 'U') ? 1 : (R_xlen_t) n + 1,
		ddx = (!packed) ? 0 : (ul == 'U') ? 1 : -1;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *pa = COMPLEX(ans), *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			(*pa).r = (*px).r;
			(*pa).i = 0.0;
			if (!root)
				(*pa).r *= (*pa).r;
			pa += 1;
			px += (dx += ddx);
		}
	} else {
		double *pa = REAL(ans), *px = REAL(x);
		for (j = 0; j < n; ++j) {
			*pa = *px;
			if (!root)
				*pa *= *pa;
			pa += 1;
			px += (dx += ddx);
		}
	}
	}
	UNPROTECT(1);
	return ans;
}

SEXP sparseCholesky_diag_get(SEXP s_trf, SEXP s_root)
{
	cholmod_factor *L = M2CHF(s_trf, 1);
	int n = (int) L->n;
	SEXP ans = Rf_allocVector((L->xtype == CHOLMOD_COMPLEX) ? CPLXSXP : REALSXP, n);
	if (n > 0) {
	int j, root = Rf_asLogical(s_root);
	if (L->is_super) {
		int k, nc,
			nsuper = (int) L->nsuper,
			*psuper = (int *) L->super,
			*ppi = (int *) L->pi,
			*ppx = (int *) L->px;
		R_xlen_t nr1a;
		if (L->xtype == CHOLMOD_COMPLEX) {
			Rcomplex *pa = COMPLEX(ans), *px = (Rcomplex *) L->x, *py;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k + 1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k + 1] - ppi[k]) + 1;
				py = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					(*pa).r = (*py).r;
					(*pa).i = 0.0;
					if (!root)
						(*pa).r *= (*pa).r;
					pa += 1;
					py += nr1a;
				}
			}
		} else {
			double *pa = REAL(ans), *px = (double *) L->x, *py;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k + 1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k + 1] - ppi[k]) + 1;
				py = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					*pa = *py;
					if (!root)
						*pa *= *pa;
					pa += 1;
					py += nr1a;
				}
			}
		}
	} else {
		int *pp = (int *) L->p;
		if (L->xtype == CHOLMOD_COMPLEX) {
			Rcomplex *pa = COMPLEX(ans), *px = (Rcomplex *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j) {
					pa[j].r = px[pp[j]].r;
					pa[j].i = 0.0;
					if (!root)
						pa[j].r *= pa[j].r;
				}
			} else {
				for (j = 0; j < n; ++j) {
					pa[j].r = px[pp[j]].r;
					pa[j].i = 0.0;
					if (root)
						pa[j].r = (pa[j].r >= 0.0) ? sqrt(pa[j].r) : R_NaN;
				}
			}
		} else {
			double *pa = REAL(ans), *px = (double *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j) {
					pa[j] = px[pp[j]];
					if (!root)
						pa[j] *= pa[j];
				}
			} else {
				for (j = 0; j < n; ++j) {
					pa[j] = px[pp[j]];
					if (root)
						pa[j] = (pa[j] >= 0.0) ? sqrt(pa[j]) : R_NaN;
				}
			}
		}
	}
	}
	return ans;
}
