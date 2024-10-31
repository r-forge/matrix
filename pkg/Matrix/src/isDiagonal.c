/* C implementation of methods for isDiagonal */

#include "Mdefines.h"
#include "M5.h"

int dense_is_diagonal(SEXP obj, const char *class)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		return 0;
	if (n <= 1)
		return 1;

	char ul = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p';

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
				px += 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
				px += 1; \
				if (!packed) \
				px += n - j - 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
				px += j; \
				px += 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} \
		return 1; \
	} while (0)

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	return 1;
}

int sparse_is_diagonal(SEXP obj, const char *class)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
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

SEXP R_dense_is_diagonal(SEXP s_obj)
{
	if (TYPEOF(s_obj) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_obj = matrix_as_dense(s_obj, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_obj);
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int ans_ = dense_is_diagonal(s_obj, class);
	SEXP ans = Rf_ScalarLogical(ans_ != 0);
	UNPROTECT(1);
	return ans;
}

/* NB: requires diagonal nonzero pattern */
SEXP R_sparse_is_diagonal(SEXP s_obj)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);
	return Rf_ScalarLogical(sparse_is_diagonal(s_obj, class));
}
