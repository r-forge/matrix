/* C implementation of methods for isTriangular */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

/* defined in ./isDiagonal.c :*/
int dense_is_diagonal(SEXP, const char *);
int sparse_is_diagonal(SEXP, const char *);

int dense_is_triangular(SEXP obj, const char *class, char op_ul)
{
	if (class[1] == 't') {
		char ul = UPLO(obj);
		if (op_ul == '\0' || op_ul == ul)
			return (   ul == 'U') ? 1 : -1;
		else if (dense_is_diagonal(obj, class))
			return (op_ul == 'U') ? 1 : -1;
		else
			return 0;
	}

	if (class[1] == 's') {
		if (!dense_is_diagonal(obj, class))
			return 0;
		else if (op_ul != '\0')
			return (op_ul == 'U') ? 1 : -1;
		else {
			char ul = UPLO(obj);
			return (   ul == 'U') ? 1 : -1;
		}
	}

	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		return 0;
	if (n <= 1)
		return (op_ul == '\0' || op_ul == 'U') ? 1 : -1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x); \
		if ((op_ul == '\0' || op_ul == 'U') && \
		    !c##NAME(test2)(px, (size_t) n, 'U', '\0', 'N')) \
			return  1; \
		if ((op_ul == '\0' || op_ul != 'U') && \
		    !c##NAME(test2)(px, (size_t) n, 'L', '\0', 'N')) \
			return -1; \
	} while (0)

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	return 0;
}

int sparse_is_triangular(SEXP obj, const char *class, char op_ul)
{
	if (class[1] == 't') {
		char ul = UPLO(obj);
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
			char ul = UPLO(obj);
			return (   ul == 'U') ? 1 : -1;
		}
	}

	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
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

#define CURL(t) \
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

#define CLRU(t) \
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
			CURL(1);
			else
			CLRU(1);
		}
		if (op_ul == '\0' || op_ul != 'U') {
			if (class[2] == 'C')
			CLRU(-1);
			else
			CURL(-1);
		}

#undef CURL
#undef CLRU

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

SEXP R_dense_is_triangular(SEXP s_obj, SEXP s_upper)
{
	if (TYPEOF(s_obj) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_obj = matrix_as_dense(s_obj, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_obj);
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int up;
	VALID_LOGIC3(s_upper, up);

	int ans_ = dense_is_triangular(s_obj, class,
		(up == NA_LOGICAL) ? '\0' : ((up != 0) ? 'U' : 'L'));
	SEXP ans = Rf_allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = ans_ != 0;
	if (up == NA_LOGICAL && ans_ != 0) {
		SEXP kind;
		PROTECT(ans);
		PROTECT(kind = Rf_mkString((ans_ > 0) ? "U" : "L"));
		Rf_setAttrib(ans, Matrix_kindSym, kind);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

/* NB: requires triangular nonzero pattern */
SEXP R_sparse_is_triangular(SEXP s_obj, SEXP s_upper)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	int up;
	VALID_LOGIC3(s_upper, up);

	int ans_ = sparse_is_triangular(s_obj, class,
		(up == NA_LOGICAL) ? '\0' : ((up != 0) ? 'U' : 'L'));
	SEXP ans = Rf_allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = ans_ != 0;
	if (up == NA_LOGICAL && ans_ != 0) {
		SEXP kind;
		PROTECT(ans);
		PROTECT(kind = Rf_mkString((ans_ > 0) ? "U" : "L"));
		Rf_setAttrib(ans, Matrix_kindSym, kind);
		UNPROTECT(2);
	}
	return ans;
}
