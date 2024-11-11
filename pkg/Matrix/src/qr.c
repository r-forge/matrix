#include "cs-etc.h"
#include "Mdefines.h"

/* defined in ./coerce.c : */
SEXP sparse_as_kind(SEXP, const char *, char);
SEXP sparse_as_general(SEXP, const char *);
SEXP sparse_as_Csparse(SEXP, const char *);

/* copy and paste in ./lu.c : */
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

/* copy and paste in ./lu.c : */
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

SEXP sparse_qr(SEXP obj, const char *class, int warn, int order)
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
		Rf_error(_("sparse QR factorization of m-by-n matrix requires m >= n"));

	Matrix_cs  *T = NULL;
	Matrix_css *S = NULL;
	Matrix_csn *N = NULL;
	int        *P = NULL;

	if (!(S = Matrix_cs_sqr(order, A, 1)) ||
	    !(N = Matrix_cs_qr(A, S)) ||
		!(P = Matrix_cs_pinv(S->pinv, S->m2)))
		goto oom;
	DO_SORT(N->L, T);
	DO_SORT(N->U, T);

	char cl[] = ".sparseQR";
	cl[0] = (A->xtype == CXSPARSE_COMPLEX) ? 'z' : 'd';
	SEXP ans = PROTECT(newObject(cl));

	SET_DIM(ans, A->m, A->n);
	SET_DIMNAMES(ans, 0, DIMNAMES(obj, 0));

	SEXP V = PROTECT(CXS2M(N->L, 1, 'g')),
		R = PROTECT(CXS2M(N->U, 1, 'g'));
	SET_SLOT(ans, Matrix_VSym, V);
	SET_SLOT(ans, Matrix_RSym, R);
	UNPROTECT(2); /* R, V */

	SEXP beta = PROTECT(Rf_allocVector(REALSXP, A->n));
	memcpy(REAL(beta), N->B, sizeof(double) * (size_t) A->n);
	SET_SLOT(ans, Matrix_betaSym, beta);
	UNPROTECT(1); /* beta */

	SEXP p = PROTECT(Rf_allocVector(INTSXP, S->m2));
	memcpy(INTEGER(p), P, sizeof(int) * (size_t) S->m2);
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
		Rf_error  (_("sparse QR factorization failed: out of memory"));
	else if (warn > 0)
		Rf_warning(_("sparse QR factorization failed: out of memory"));
	UNPROTECT(1); /* obj */
	return R_NilValue;
}

SEXP R_sparse_qr(SEXP s_obj, SEXP s_warn, SEXP s_order)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);
	int order = Rf_asInteger(s_order);
	if (order < 0 || order > 3)
		order = 0;
	int cache =
		(class[1] == 'g' || class[1] == 's') &&
		(class[0] == 'z' || class[0] == 'd');
	const char *nm = (order == 0) ? "sparseQR-" : "sparseQR+";
	SEXP ans = (cache) ? get_factor(s_obj, nm) : R_NilValue;
	if (ans == R_NilValue) {
		int warn = Rf_asLogical(s_warn);
		ans = sparse_qr(s_obj, class, warn, order);
		if (cache) {
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
		}
	}
	return ans;
}
