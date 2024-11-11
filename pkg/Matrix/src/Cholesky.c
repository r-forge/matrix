#include "Lapack-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP sparse_as_kind(SEXP, const char *, char);
SEXP sparse_as_general(SEXP, const char *);
SEXP sparse_as_Csparse(SEXP, const char *);

/* defined in ./forceCanonical.c : */
SEXP dense_force_canonical(SEXP, const char *, int);

/* defined in ./forceSymmetric.c : */
SEXP sparse_force_symmetric(SEXP, const char *, char, char);

/* defined in ./t.c : */
SEXP sparse_transpose(SEXP, const char *, char, int);

SEXP dense_cholesky(SEXP obj, const char *class, int warn, int pivot,
                    double tol, char ul)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("matrix is not square"));
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(obj, &pid);
	if (class[0] != 'z' && class[0] != 'd') {
		REPROTECT(obj = dense_as_kind(obj, class, ',', 1), pid);
		class = Matrix_class(obj, valid_dense, 6, __func__);
	}
	if (class[1] == 't' && DIAG(obj) != 'N')
		REPROTECT(obj = dense_force_canonical(obj, class, 0), pid);
	ul = (class[1] != 'g') ? UPLO(obj) : (ul == '\0') ? 'U' : ul;
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseCholesky";
	cl[0] = class[0];
	SEXP ans = PROTECT(newObject(cl));
	SET_DIM(ans, n, n);
	SET_DIMNAMES(ans, -(class[1] == 's'), DIMNAMES(obj, 0));
	if (ul != 'U')
		SET_UPLO(ans);
	if (n > 0) {
	SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
	int info;
	if (class[2] != 'p') {
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
	memset(py, 0, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	F77_CALL(zlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
	if (!pivot) {
	F77_CALL(zpotrf)(&ul, &n, py, &n, &info FCONE);
	ERROR_LAPACK_3(zpotrf, info, warn);
	} else {
	SEXP perm = PROTECT(Rf_allocVector(INTSXP, n));
	int *pperm = INTEGER(perm), rank;
	double *work = (double *) R_alloc((size_t) n * 2, sizeof(double));
	F77_CALL(zpstrf)(&ul, &n, py, &n, pperm, &rank, &tol, work,
	                 &info FCONE);
	ERROR_LAPACK_4(zpstrf, info, warn, rank);
	if (info > 0) {
		int j, d = n - rank;
		py += (R_xlen_t) rank * n + rank;
		for (j = rank; j < n; ++j) {
			memset(py, 0, sizeof(Rcomplex) * (size_t) d);
			py += n;
		}
	}
	SET_SLOT(ans, Matrix_permSym, perm);
	UNPROTECT(1); /* perm */
	}
	} else {
	double *px = REAL(x), *py = REAL(y);
	memset(py, 0, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
	if (!pivot) {
	F77_CALL(dpotrf)(&ul, &n, py, &n, &info FCONE);
	ERROR_LAPACK_3(dpotrf, info, warn);
	} else {
	SEXP perm = PROTECT(Rf_allocVector(INTSXP, n));
	int *pperm = INTEGER(perm), rank;
	double *work = (double *) R_alloc((size_t) n * 2, sizeof(double));
	F77_CALL(dpstrf)(&ul, &n, py, &n, pperm, &rank, &tol, work,
	                 &info FCONE);
	ERROR_LAPACK_4(dpstrf, info, warn, rank);
	if (info > 0) {
		int j, d = n - rank;
		py += (R_xlen_t) rank * n + rank;
		for (j = rank; j < n; ++j) {
			memset(py, 0, sizeof(double) * (size_t) d);
			py += n;
		}
	}
	SET_SLOT(ans, Matrix_permSym, perm);
	UNPROTECT(1); /* perm */
	}
	}
	} else {
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
	memcpy(py, px, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	F77_CALL(zpptrf)(&ul, &n, py, &info FCONE);
	ERROR_LAPACK_3(zpptrf, info, warn);
	} else {
	double *px = REAL(x), *py = REAL(y);
	memcpy(py, px, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dpptrf)(&ul, &n, py, &info FCONE);
	ERROR_LAPACK_3(dpptrf, info, warn);
	}
	}
	SET_SLOT(ans, Matrix_xSym, y);
	UNPROTECT(1); /* y */
	}
	UNPROTECT(3); /* trf, x, obj */
	return ans;
}

SEXP sparse_cholesky(SEXP obj, const char *class, int warn, int order,
                     int *ll, int *super, Rcomplex beta, int force,
                     char ul, SEXP trf)
{
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(obj, &pid);
	if (class[0] != 'z' && class[0] != 'd') {
		REPROTECT(obj = sparse_as_kind(obj, class, ','), pid);
		class = Matrix_class(obj, valid_sparse, 6, __func__);
	}
	if (class[1] == 's' && (class[0] != 'z' || TRANS(obj) == 'C'))
		;
	else {
		if (force)
		REPROTECT(obj = sparse_force_symmetric(obj, class, ul, 'C'), pid);
		else
		REPROTECT(obj = sparse_as_general(obj, class), pid);
		class = Matrix_class(obj, valid_sparse, 6, __func__);
	}
	if (class[2] != 'C') {
		if (class[2] == 'R' && class[1] == 's')
		REPROTECT(obj = sparse_transpose(obj, class, 'C', 1), pid);
		else
		REPROTECT(obj = sparse_as_Csparse(obj, class), pid);
		class = Matrix_class(obj, valid_sparse, 6, __func__);
	}

	cholmod_sparse *A =                              M2CHS(obj, 1);
	cholmod_factor *L = (trf == R_NilValue) ? NULL : M2CHF(trf, 1);
	double b[2]; b[0] = beta.r; b[1] = beta.i;

	A->stype = (class[1] != 's') ? 0 : (UPLO(obj) == 'U') ? 1 : -1;

	if (L)
		L = cholmod_copy_factor(L, &c);
	else {
		if (order == 0) {
			c.nmethods = 1;
			c.method[0].ordering = CHOLMOD_NATURAL;
			c.postorder = 0;
		}
		c.supernodal = (!super || *super == NA_LOGICAL) ? CHOLMOD_AUTO :
			((*super != 0) ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL);
		L = cholmod_analyze(A, &c);
	}

	if (super)
		*super = L->is_super != 0;
	if (ll)
		*ll    = L->is_super != 0 || L->is_ll != 0;

	c.final_asis = 0;
	c.final_ll = ((ll) ? *ll : L->is_ll) != 0;
	c.final_super = ((super) ? *super : L->is_super) != 0;
	c.final_pack = 1;
	c.final_monotonic = 1;

	cholmod_factorize_p(A, b, NULL, 0, L, &c);
	cholmod_defaults(&c);

	SEXP ans = PROTECT(CHF2M(L, 1));
	cholmod_free_factor(&L, &c);
	if (TYPEOF(ans) == CHARSXP) {
		if (warn > 1)
			Rf_error  ("%s", CHAR(ans));
		else if (warn > 0)
			Rf_warning("%s", CHAR(ans));
		ans = R_NilValue;
	} else {
		SEXP adimnames = PROTECT(Rf_allocVector(VECSXP, 2)),
			odimnames = PROTECT(DIMNAMES(obj, 0));
		symDN(adimnames, odimnames, (class[1] != 's') ? 0 : -1);
		SET_DIMNAMES(ans, 0, adimnames);
		UNPROTECT(2); /* odimnames, adimnames */
	}
	UNPROTECT(2); /* ans, obj */
	return ans;
}

SEXP R_dense_cholesky(SEXP s_obj, SEXP s_warn, SEXP s_pivot,
                      SEXP s_tol, SEXP s_uplo)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);
	int pivot = Rf_asLogical(s_pivot);
	double tol = Rf_asReal(s_tol);
	if (ISNAN(tol))
		Rf_error(_("'%s' is not a number"), "tol");
	char ul = '\0';
	if (s_uplo != R_NilValue)
		VALID_UPLO(s_uplo, ul);
	int cache = (class[1] == 's') &&
		((class[0] == 'z' && TRANS(s_obj) == 'C') || class[0] == 'd');
	char nm[] = "denseCholesky..";
	if (class[1] == 's' && class[2] == 'p')
		nm[13] = nm[14] = '-';
	else {
		nm[13] = '+';
		nm[14] = (!pivot) ? '-' : '+';
	}
	SEXP ans = (cache) ? get_factor(s_obj, nm) : R_NilValue;
	if (ans == R_NilValue) {
		int warn = Rf_asInteger(s_warn);
		ans = dense_cholesky(s_obj, class, warn, pivot, tol, ul);
		if (cache) {
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
		}
	}
	return ans;
}

SEXP R_sparse_cholesky(SEXP s_obj, SEXP s_warn, SEXP s_order,
                       SEXP s_ll, SEXP s_super, SEXP s_beta,
                       SEXP s_force, SEXP s_uplo)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);
	int warn = Rf_asInteger(s_warn), order = Rf_asInteger(s_order),
		ll = Rf_asLogical(s_ll), super = Rf_asLogical(s_super),
		force = Rf_asLogical(s_force);
	Rcomplex beta = Rf_asComplex(s_beta);
	if (order < 0 || order > 1)
		order = 0;
	if (!R_FINITE(beta.r) || !R_FINITE(beta.i))
		Rf_error(_("'%s' is not a number or not finite"), "beta");
	char ul = '\0';
	if (s_uplo != R_NilValue)
		VALID_UPLO(s_uplo, ul);
	int cache = (class[1] == 's') &&
		((class[0] == 'z' && TRANS(s_obj) == 'C') || class[0] == 'd');
	char nm[] = "..........Cholesky.";
	nm[18] = (order == 0) ? '-' : '+';
	SEXP ans = R_NilValue;
	if (cache) {
	if (super == NA_LOGICAL || super == 0) {
		memcpy(nm, "simplicial", 10);
		ans = get_factor(s_obj, nm);
		if (ans != R_NilValue) super = 0;
	}
	if (ans == R_NilValue && (super == NA_LOGICAL || super != 0)) {
		memcpy(nm, "supernodal", 10);
		ans = get_factor(s_obj, nm);
		if (ans != R_NilValue) super = 1;
	}
	}
	if (beta.r != 0.0 || beta.i != 0.0 || ans == R_NilValue ||
	    (super == 0 && ll != 0)) {
		if (beta.r != 0.0 || beta.i != 0.0) {
			PROTECT(ans);
			ans = sparse_cholesky(s_obj, class, warn,
			                      order, &ll, &super, beta,
			                      force, ul, ans);
			UNPROTECT(1);
		} else {
			if (ans == R_NilValue) {
			int zz = 0;
			ans = sparse_cholesky(s_obj, class, warn,
			                      order, &zz, &super, beta,
			                      force, ul, ans);
			if (cache && ans != R_NilValue) {
			memcpy(nm, (super == 0) ? "simplicial" : "supernodal", 10);
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
			}
			}
			if (ans != R_NilValue && super == 0 && ll != 0) {
			PROTECT(ans);
			cholmod_factor *L = M2CHF(ans, 1);
			L = cholmod_copy_factor(L, &c);
			cholmod_change_factor(L->xtype, 1, 0, 1, 1, L, &c);
			SEXP tmp = PROTECT(CHF2M(L, 1));
			cholmod_free_factor(&L, &c);
			if (TYPEOF(tmp) == CHARSXP) {
				if (warn > 1)
					Rf_error  ("%s", CHAR(ans));
				else if (warn > 0)
					Rf_warning("%s", CHAR(ans));
				ans = R_NilValue;
			} else {
				SET_DIMNAMES(tmp, 0, DIMNAMES(ans, 0));
				ans = tmp;
			}
			UNPROTECT(2); /* tmp, ans */
			}
		}
	}
	return ans;
}
