/* C implementation of methods for forceCanonical */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

/* defined in ./aggregate.c : */
SEXP sparse_aggregate(SEXP, const char *);

SEXP dense_force_canonical(SEXP from, const char *class, int check)
{
	if (class[1] == 'g' && class[0] != 'n')
		return from;
	SEXP x = GET_SLOT(from, Matrix_xSym);
	int naToUnitIfPattern = (check) ? 0 : 1;
	if (class[0] == 'n' && !naToUnitIfPattern) {
		int *px = LOGICAL(x);
		R_xlen_t nx = XLENGTH(x);
		while (nx-- > 0) if (*(px++) == NA_LOGICAL) break;
		naToUnitIfPattern = nx >= 0;
		if (class[1] == 'g' && !naToUnitIfPattern)
			return from;
	}
	PROTECT(x);
	int *pdim = DIM(from), m = pdim[0], n = pdim[1],
		packed = class[2] == 'p', canonical = !naToUnitIfPattern;
	char ul = '\0', ct = (class[1] == 'p') ? 'C' : '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(from);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(from);
	if (class[1] == 't')
		nu = DIAG(from);
	if (canonical) {
#define TEMPLATE(c) \
		do { \
			c##TYPE *px = c##PTR(x); \
			canonical = (!packed) \
				? !c##NAME(test2)(px, (size_t) n, ul, ct, nu) \
				: !c##NAME(test1)(px, (size_t) n, ul, ct, nu); \
		} while (0)
		SWITCH5(class[0], TEMPLATE);
#undef TEMPLATE
		if (canonical) {
			UNPROTECT(1); /* x */
			return from;
		}
	}
	SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x), *py = c##PTR(y); \
		if (class[1] == 'g') \
			memcpy(py, px, sizeof(c##TYPE) * (size_t) m * (size_t) n); \
		else if (!packed) \
			c##NAME(force2)(py, px, (size_t) n, ul, ct, nu); \
		else \
			c##NAME(force1)(py, px, (size_t) n, ul, ct, nu); \
	} while (0)
	SWITCH4(class[0], TEMPLATE);
#undef TEMPLATE
	if (class[0] == 'n' && naToUnitIfPattern) {
		int *py = LOGICAL(y);
		R_xlen_t ny = XLENGTH(y);
		while (ny-- > 0) { *py = *py != 0; ++py; }
	}
	SEXP to = PROTECT(newObject(class));
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && ul != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && ct != 'C' && class[0] == 'z')
		SET_TRANS(to);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);
	if (class[1] == 'o')
		COPY_SLOT(to, from, Matrix_sdSym);
	SET_SLOT(to, Matrix_xSym, y);
	UNPROTECT(3); /* to, y, x */
	return to;
}

SEXP sparse_force_canonical(SEXP from, const char *class, int check)
{
	SEXP to = PROTECT(sparse_aggregate(from, class));
	switch (class[1]) {
	case 'g':
		break;
	case 's':
	case 'p':
		if (class[0] == 'z' && TRANS(from) == 'C') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		Rcomplex *px = COMPLEX(x);
		int n = DIM(from)[1], canonical = (check) ? 1 : 0;
		char ul = UPLO(from);
		if (class[2] != 'T') {
			SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
				p = PROTECT(GET_SLOT(from, Matrix_pSym)),
				i = PROTECT(GET_SLOT(from, iSym));
			int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k_, k, kend,
				up = (class[2] == 'C') == (ul == 'U');
			j = 0; k = 0;
			if (canonical) {
				for (; j < n; ++j) {
					kend = pp[j];
					if (k < kend && pi[k_ = (up) ? kend - 1 : k] == j &&
					    (ISNAN(px[k_].i) || px[k_].i != 0.0))
						break;
					k = kend;
				}
				canonical = j == n;
			}
			if (!canonical) {
				SEXP y = PROTECT(duplicateVector(x));
				Rcomplex *py = COMPLEX(y);
				for (; j < n; ++j) {
					kend = pp[j];
					if (k < kend && pi[k_ = (up) ? kend - 1 : k] == j)
						py[k_].r = 0.0;
				}
				PROTECT(to = newObject(class));
				SET_DIM(to, n, n);
				SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
				if (ul != 'U')
					SET_UPLO(to);
				SET_SLOT(to, Matrix_pSym, p);
				SET_SLOT(to, Matrix_iSym, i);
				SET_SLOT(to, Matrix_xSym, y);
				UNPROTECT(2); /* to, y */
			}
			UNPROTECT(2); /* i, p */
		} else {
			SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym)),
				j = PROTECT(GET_SLOT(from, Matrix_jSym));
			int *pi = INTEGER(i), *pj = INTEGER(j);
			R_xlen_t k = 0, kend = XLENGTH(i);
			if (canonical) {
				for (; k < kend; ++k)
					if (pi[k] == pj[k] &&
					    (ISNAN(px[k].i) || px[k].i != 0.0))
						break;
				canonical = k == kend;
			}
			if (!canonical) {
				SEXP y = PROTECT(duplicateVector(x));
				Rcomplex *py = COMPLEX(y);
				for (; k < kend; ++k)
					if (pi[k] == pj[k])
						py[k].i = 0.0;
				PROTECT(to = newObject(class));
				SET_DIM(to, n, n);
				SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
				if (ul != 'U')
					SET_UPLO(to);
				SET_SLOT(to, Matrix_iSym, i);
				SET_SLOT(to, Matrix_jSym, j);
				SET_SLOT(to, Matrix_xSym, y);
				UNPROTECT(2); /* to, y */
			}
			UNPROTECT(2); /* j, i */
		}
		UNPROTECT(1); /* x */
		}
	case 'o':
		break;
	case 't':
		if (DIAG(from) != 'N') {
		/* defined in ./diag.c : */
		SEXP sparse_diag_set(SEXP, const char *, SEXP);
		SEXP value = R_NilValue;
#define TEMPLATE(c) \
		do { \
			value = Rf_allocVector(c##TYPESXP, 1); \
			c##PTR(value)[0] = c##UNIT; \
		} while (0)
		SWITCH4(class[0], TEMPLATE);
#undef TEMPLATE
		PROTECT(value);
		to = sparse_diag_set(from, class, value);
		UNPROTECT(1); /* value */
		}
		break;
	default:
		Rf_error("should never happen ...");
		to = R_NilValue;
		break;
	}
	UNPROTECT(1); /* to */
	return to;
}

SEXP R_dense_force_canonical(SEXP s_from, SEXP s_check)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);

	int check;
	VALID_LOGIC2(s_check, check);

	return dense_force_canonical(s_from, class, check);
}

SEXP R_sparse_force_canonical(SEXP s_from, SEXP s_check)
{
	const char *class = Matrix_class(s_from, valid_sparse, 0, __func__);

	int check;
	VALID_LOGIC2(s_check, check);

	return sparse_force_canonical(s_from, class, check);
}
