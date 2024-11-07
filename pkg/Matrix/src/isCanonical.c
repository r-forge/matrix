/* C implementation of methods for isCanonical */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

int dense_is_canonical(SEXP obj, const char *class)
{
	if (class[1] == 'g' && class[0] != 'n')
		return 1;
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (class[0] == 'n') {
		int *px = LOGICAL(x);
		R_xlen_t nx = XLENGTH(x);
		while (nx-- > 0) if (*(px++) == NA_LOGICAL) return 0;
		if (class[1] == 'g')
			return 1;
	}
	PROTECT(x);
	int n = DIM(obj)[1], packed = class[2] == 'p';
	char ul = UPLO(obj), ct = (class[1] == 'p') ? 'C' : '\0', nu = '\0';
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);
	UNPROTECT(1); /* x */
#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x); \
		return (!packed) \
			? !c##NAME(test2)(px, (size_t) n, ul, ct, nu) \
			: !c##NAME(test1)(px, (size_t) n, ul, ct, nu); \
	} while (0)
	SWITCH5(class[0], TEMPLATE);
#undef TEMPLATE
	Rf_error("should never happen ...");
	return 0;
}

int sparse_is_canonical(SEXP obj, const char *class)
{
	switch (class[1]) {
	case 'g':
		return 1;
	case 's':
		if (class[0] == 'z' && TRANS(obj) == 'C') {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		Rcomplex *px = COMPLEX(x);
		int n = DIM(obj)[1];
		char ul = UPLO(obj);
		if (class[2] != 'T') {
			SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
				p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
				i = PROTECT(GET_SLOT(obj, iSym));
			int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k_, k, kend,
				up = (class[2] == 'C') == (ul == 'U');
			UNPROTECT(3); /* i, p, x */
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				if (k < kend && pi[k_ = (up) ? kend - 1 : k] == j &&
				    (ISNAN(px[k_].i) || px[k_].i != 0.0))
					return 0;
				k = kend;
			}
		} else {
			SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
				j = PROTECT(GET_SLOT(obj, Matrix_jSym));
			int *pi = INTEGER(i), *pj = INTEGER(j);
			R_xlen_t k, kend = XLENGTH(i);
			UNPROTECT(3); /* j, i, x */
			for (k = 0; k < kend; ++k)
				if (pi[k] == pj[k] &&
				    (ISNAN(px[k].i) || px[k].i != 0.0))
					return 0;			
		}
		}
		return 1;
	case 't':
		return DIAG(obj) == 'N';
	default:
		Rf_error("should never happen ...");
		return 0;
	}
}

SEXP R_dense_is_canonical(SEXP s_obj)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);
	return Rf_ScalarLogical(dense_is_canonical(s_obj, class));
}

SEXP R_sparse_is_canonical(SEXP s_obj)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);
	return Rf_ScalarLogical(sparse_is_canonical(s_obj, class));
}
