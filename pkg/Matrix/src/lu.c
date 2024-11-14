#include "Lapack-etc.h"
#include "cs-etc.h"
#include "Mdefines.h"

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP dense_as_general(SEXP, const char *, int);
SEXP sparse_as_kind(SEXP, const char *, char);
SEXP sparse_as_general(SEXP, const char *);
SEXP sparse_as_Csparse(SEXP, const char *);

SEXP dense_lu(SEXP obj, const char *class, int warn)
{
	char cl[] = ".denseLU";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	SEXP ans = PROTECT(newObject(cl));
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	SET_DIM(ans, m, n);
	SET_DIMNAMES(ans, -(class[1] == 's'), DIMNAMES(obj, 0));
	if (r > 0) {
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(obj, &pid);
	if (class[0] != 'z' && class[0] != 'd') {
		REPROTECT(obj = dense_as_kind(obj, class, ',', 1), pid);
		class = Matrix_class(obj, valid_dense, 6, __func__);
	}
	if (class[1] != 'g')
		REPROTECT(obj = dense_as_general(obj, class, 1), pid);
	SEXP perm = PROTECT(Rf_allocVector(INTSXP, r)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
	int *pperm = INTEGER(perm), info;
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
	memcpy(py, px, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	F77_CALL(zgetrf)(&m, &n, py, &m, pperm, &info);
	ERROR_LAPACK_2(zgetrf, info, warn, U);
	} else {
	double *px = REAL(x), *py = REAL(y);
	memcpy(py, px, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dgetrf)(&m, &n, py, &m, pperm, &info);
	ERROR_LAPACK_2(dgetrf, info, warn, U);
	}
	SET_SLOT(ans, Matrix_permSym, perm);
	SET_SLOT(ans, Matrix_xSym, y);
	UNPROTECT(4); /* y, x, perm, obj */
	}
	UNPROTECT(1); /* ans */
	return ans;
}

/* copy and paste in ./qr.c : */
#define DO_FREE(_T_, _S_, _N_, _P_) \
do { \
	if (!(_T_)) \
		_T_ = Matrix_cs_spfree(_T_); \
	if (!(_S_)) \
		_S_ = Matrix_cs_sfree (_S_); \
	if (!(_N_)) \
		_N_ = Matrix_cs_nfree (_N_); \
	if (!(_P_)) \
		_P_ = Matrix_cs_free  (_P_); \
} while (0)

/* copy and paste in ./qr.c : */
#define DO_SORT(_A_, _T_) \
do { \
	Matrix_cs_dropzeros(_A_); \
	_T_ = Matrix_cs_transpose(_A_, 1); \
	if (!_T_) \
		goto oom; \
	_A_ = Matrix_cs_spfree(_A_); \
	_A_ = Matrix_cs_transpose(_T_, 1); \
	if (!_A_) \
		goto oom; \
	_T_ = Matrix_cs_spfree(_T_); \
} while (0)

SEXP sparse_lu(SEXP obj, const char *class, int warn, int order,
               double tol)
{
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(obj, &pid);
	if (class[0] != 'z' && class[0] != 'd') {
		REPROTECT(obj = sparse_as_kind(obj, class, ','), pid);
		class = Matrix_class(obj, valid_sparse, 6, __func__);
	}
	if (class[1] != 'g') {
		REPROTECT(obj = sparse_as_general(obj, class), pid);
		class = Matrix_class(obj, valid_sparse, 6, __func__);
	}
	if (class[2] != 'C') {
		REPROTECT(obj = sparse_as_Csparse(obj, class), pid);
		class = Matrix_class(obj, valid_sparse, 6, __func__);
	}

	Matrix_cs *A = M2CXS(obj, 1);
	CXSPARSE_XTYPE_SET(A->xtype);

	if (A->m < A->n)
		Rf_error(_("sparse LU factorization of m-by-n matrix requires m == n"));

	Matrix_cs  *T = NULL;
	Matrix_css *S = NULL;
	Matrix_csn *N = NULL;
	int        *P = NULL;

	if (!(S = Matrix_cs_sqr(order, A, 0)) ||
	    !(N = Matrix_cs_lu(A, S, tol)) ||
	    !(P = Matrix_cs_pinv(N->pinv, A->m)))
		goto oom;
	DO_SORT(N->L, T);
	DO_SORT(N->U, T);

	char cl[] = ".sparseLU";
	cl[0] = (A->xtype == CXSPARSE_COMPLEX) ? 'z' : 'd';
	SEXP ans = PROTECT(newObject(cl));

	SET_DIM(ans, A->m, A->n);
	SET_DIMNAMES(ans, 0, DIMNAMES(obj, 0));

	SEXP L = PROTECT(CXS2M(N->L, 1, 't')),
		U = PROTECT(CXS2M(N->U, 1, 't'));
	SET_UPLO(L);
	SET_SLOT(ans, Matrix_LSym, L);
	SET_SLOT(ans, Matrix_USym, U);
	UNPROTECT(2); /* U, L */

	SEXP p = PROTECT(Rf_allocVector(INTSXP, A->m));
	memcpy(INTEGER(p), P, sizeof(int) * (size_t) A->m);
	SET_SLOT(ans, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order > 0) {
	SEXP q = PROTECT(Rf_allocVector(INTSXP, A->n));
	memcpy(INTEGER(q), S->q, sizeof(int) * (size_t) A->n);
	SET_SLOT(ans, Matrix_qSym, q);
	UNPROTECT(1); /* q */
	}

	DO_FREE(T, S, N, P);
	UNPROTECT(2); /* ans, obj */
	return ans;

oom:
	DO_FREE(T, S, N, P);
	if (warn > 1)
		Rf_error  (_("sparse LU factorization failed: out of memory or near-singular"));
	else if (warn > 0)
		Rf_warning(_("sparse LU factorization failed: out of memory or near-singular"));
	UNPROTECT(1); /* obj */
	return R_NilValue;
}

SEXP R_dense_lu(SEXP s_obj, SEXP s_warn)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);
	int cache =
		(class[1] == 'g' || class[1] == 's') &&
		(class[0] == 'z' || class[0] == 'd');
	const char *nm = "denseLU";
	SEXP ans = (cache) ? get_factor(s_obj, nm) : R_NilValue;
	if (ans == R_NilValue) {
		int warn = Rf_asLogical(s_warn);
		ans = dense_lu(s_obj, class, warn);
		if (cache) {
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
		}
	}
	return ans;
}

SEXP R_sparse_lu(SEXP s_obj, SEXP s_warn, SEXP s_order, SEXP s_tol)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);
	double tol = Rf_asReal(s_tol);
	if (ISNAN(tol))
		Rf_error(_("'%s' is not a number"), "tol");
	int order = Rf_asInteger(s_order);
	if (order == NA_INTEGER)
		order = (tol == 1.0) ? 2 : 1;
	else if (order < 0 || order > 3)
		order = 0;
	int cache =
		(class[1] == 'g' || class[1] == 's') &&
		(class[0] == 'z' || class[0] == 'd');
	const char *nm = (order == 0) ? "sparseLU-" : "sparseLU+";
	SEXP ans = (cache) ? get_factor(s_obj, nm) : R_NilValue;
	if (ans == R_NilValue) {
		int warn = Rf_asInteger(s_warn);
		ans = sparse_lu(s_obj, class, warn, order, tol);
		if (cache) {
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
		}
	}
	return ans;
}
