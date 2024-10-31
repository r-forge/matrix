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
