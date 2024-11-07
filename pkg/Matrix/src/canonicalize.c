#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

SEXP dense_canonicalize(SEXP from, const char *class)
{
	int *pdim = DIM(from), m = pdim[0], n = pdim[1],
		packed = class[2] == 'p';

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(from);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(from);
	if (class[1] == 'p' && class[0] == 'z')
		ct = 'C';
	if (class[1] == 't')
		nu = DIAG(from);

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	int patternWithNA = 0;
	if (class[0] == 'n') {
		int *px = LOGICAL(x);
		R_xlen_t nx = XLENGTH(x);
		while (nx-- > 0) if (*(px++) == NA_LOGICAL) break;
		patternWithNA = nx >= 0;
	}
	int canonical = !patternWithNA;
	if (canonical && class[1] != 'g') {
#define TEMPLATE(c) \
		do { \
			c##TYPE *px = c##PTR(x); \
			canonical = (!packed) \
				? !c##NAME(test2)(px, (size_t) n, ul, ct, nu) \
				: !c##NAME(test1)(px, (size_t) n, ul, ct, nu); \
		} while (0)
		SWITCH5(class[0], TEMPLATE);
#undef TEMPLATE
	}

	SEXP to = from;
	if (!canonical) {
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
	if (patternWithNA) {
		int *py = LOGICAL(y);
		R_xlen_t ny = XLENGTH(y);
		while (ny-- > 0) { *py = *py != 0; ++py; }
	}
	PROTECT(to = newObject(class));
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
	UNPROTECT(2); /* to, y */
	}
	UNPROTECT(1); /* x */ 
	return to;
}

SEXP sparse_canonicalize(SEXP from, const char *class)
{
	SEXP to = from;
	switch (class[1]) {
	case 'g':
		break;
	case 's':
	case 'p':
		if (class[0] == 'z' && TRANS(from) == 'C') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		Rcomplex *px = COMPLEX(x);
		int n = DIM(from)[1];
		char ul = UPLO(from);
		if (class[2] != 'T') {
			SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
				p = PROTECT(GET_SLOT(from, Matrix_pSym)),
				i = PROTECT(GET_SLOT(from, iSym));
			int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k_, k, kend,
				up = (class[2] == 'C') == (ul == 'U');
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j];
				if (k < kend && pi[k_ = (up) ? kend - 1 : k] == j &&
				    (ISNAN(px[k_].i) || px[k_].i != 0.0))
					break;
				k = kend;
			}
			if (j < n) {
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
			R_xlen_t k, kend = XLENGTH(i);
			for (k = 0; k < kend; ++k) {
				if (pi[k] == pj[k] &&
				    (ISNAN(px[k].i) || px[k].i != 0.0))
					break;
			}
			if (k < kend) {
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
		break;
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
	}
	return to;
}

SEXP R_dense_canonicalize(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);
	return dense_canonicalize(s_from, class);
}

SEXP R_sparse_canonicalize(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 0, __func__);
	return sparse_canonicalize(s_from, class);
}
