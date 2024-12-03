/* C implementation of methods for forceSymmetric */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

SEXP dense_force_symmetric(SEXP from, const char *class,
                           char op_ul, char op_ct)
{
	if (class[0] != 'z')
		op_ct = '\0';

	char
		ul0 = '\0', ul1 = 'U',
		ct0 = '\0', ct1 = (class[0] == 'z') ? 'C' : '\0',
		nu0 = '\0';
	if (class[1] != 'g')
		ul0 = ul1 = UPLO(from);
	if (class[1] == 's' && class[0] == 'z')
		ct0 = ct1 = TRANS(from);
	if (class[1] == 't')
		nu0 = DIAG(from);
	if (op_ul != '\0')
		ul1 = op_ul;
	if (op_ct != '\0')
		ct1 = op_ct;

	if (class[1] == 's' && ul0 == ul1 && ct0 == ct1)
		return from;

	int packed = class[2] == 'p';

	char cl[] = ".s.Matrix";
	cl[0] = class[0];
	cl[2] = (packed) ? 'p' : 'y';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("matrix is not square"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));
	if (ul1 != 'U')
		SET_UPLO(to);
	if (ct1 != 'C' && ct1 != '\0')
		SET_TRANS(to);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));

	if ((class[1] == 'g') ||
	    (class[1] == 's' && ul0 == ul1 && ct0 != 'C') ||
	    (class[1] == 't' && ul0 == ul1 && nu0 == 'N'))
		SET_SLOT(to, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));
		size_t n_ = (size_t) n;

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			if (!packed) \
				c##NAME(force2)(px1,  px0, n_, ul0, ct0, nu0); \
			else if (ul0 == ul1) \
				c##NAME(force1)(px1,  px0, n_, ul0, ct0, nu0); \
			else { \
				c##NAME(trans1)(px1,  px0, n_, ul0, ct0); \
				c##NAME(force1)(px1, NULL, n_, ul1, ct0, nu0); \
			} \
		} while (0)

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		SET_SLOT(to, Matrix_xSym, x1);
		UNPROTECT(1); /* x1 */
	}

	UNPROTECT(2); /* x0, to */
	return to;
}

SEXP sparse_force_symmetric(SEXP from, const char *class,
                            char op_ul, char op_ct)
{
	if (class[0] != 'z')
		op_ct = '\0';

	char
		ul0 = '\0', ul1 = 'U',
		ct0 = '\0', ct1 = (class[0] == 'z') ? 'C' : '\0',
		nu0 = '\0';
	if (class[1] != 'g')
		ul0 = ul1 = UPLO(from);
	if (class[1] == 's' && class[0] == 'z')
		ct0 = ct1 = TRANS(from);
	if (class[1] == 't')
		nu0 = DIAG(from);
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

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("matrix is not square"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));
	if (ul1 != 'U')
		SET_UPLO(to);
	if (ct1 != 'C' && ct1 != '\0')
		SET_TRANS(to);

	int up = (class[0] != 'R') == (ul1 == 'U');

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			j, k, kend, nnz0 = pp0[n], nnz1 = 0;
		pp0++;

		if ((class[1] == 's' && ul0 == ul1) ||
		    (class[1] == 't' && ul0 == ul1 && nu0 == 'N')) {

		SET_SLOT(to, Matrix_pSym, p0);
		SET_SLOT(to,        iSym, i0);

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (class[1] != 's' || ct0 != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, nnz0));
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

		SEXP p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0)));
		int *pp1 = INTEGER(p1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] == 'g')
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((up) ? pi0[k] <= j : pi0[k] >= j)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
		else if (class[1] == 's')
			nnz1 = nnz0;
		else if (nu0 == 'N')
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				if (k < kend && pi0[(up) ? k : kend - 1] == j)
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

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz1)); \
			px0 = c##PTR(x0); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] == 'g') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((up) ? pi0[k] <= j : pi0[k] >= j) { \
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
				c##csptrans(pp1, pi1, px1, pp0, pi0, px0, n, n, ct0, iwork); \
				Matrix_Free(iwork, n); \
				if (ct0 == 'C') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp1[j]; \
					if (k < kend && pi1[(up) ? kend - 1 : k] == j) \
						c##SET_PROJ_REAL(px1[(up) ? kend - 1 : k]); \
					k = kend; \
				} \
			} \
			else if (nu0 == 'N') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					if (k < kend && pi0[(up) ? k : kend - 1] == j) { \
						*(pi1++) = j; \
						c##IF_NPATTERN( \
						*(px1++) = px0[(up) ? k : kend - 1]; \
						); \
					} \
					k = kend; \
				} \
			else if (ul0 == ul1) \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					if (!up) { \
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
					if (up) { \
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

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(2); /* i1, p1 */

		}

		UNPROTECT(2); /* i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;

		if ((class[1] == 's') ||
		    (class[1] == 't' && ul0 == ul1 && nu0 == 'N')) {

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
				SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, XLENGTH(x0)));
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

		if (class[1] == 'g' || nu0 == 'N') {
			for (k = 0; k < nnz0; ++k)
				if ((up) ? pi0[k] <= pj0[k] : pi0[k] >= pj0[k])
					++nnz1;
		}
		else
			nnz1 = (ul0 == ul1) ? n + nnz0 : n;

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
			j1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1), j;
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
			if (class[1] == 'g' || nu0 == 'N') \
				for (k = 0; k < nnz0; ++k) { \
					if ((up) ? pi0[k] <= pj0[k] : pi0[k] >= pj0[k]) { \
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

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(2); /* j1, i1 */

		}

		UNPROTECT(2); /* j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

SEXP R_dense_force_symmetric(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ul = '\0', ct = '\0';
	if (s_uplo != R_NilValue)
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return dense_force_symmetric(s_from, class, ul, ct);
}

SEXP R_sparse_force_symmetric(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	char ul = '\0', ct = '\0';
	if (s_uplo != R_NilValue)
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return sparse_force_symmetric(s_from, class, ul, ct);
}
