#include "Lapack-etc.h"
#include "Mdefines.h"

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);

/* defined in ./forceCanonical.c : */
SEXP dense_force_canonical(SEXP, const char *, int);

SEXP dense_bunchkaufman(SEXP obj, const char *class, int warn,
                        char ul, char ct)
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
	ct = (class[1] == 's' && class[0] == 'z') ? TRANS(obj) : ct;
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseBunchKaufman";
	cl[0] = class[0];
	SEXP ans = PROTECT(newObject(cl));
	SET_DIM(ans, n, n);
	SET_DIMNAMES(ans, -(class[1] == 's'), DIMNAMES(obj, 0));
	if (ul != 'U')
		SET_UPLO(ans);
	if (ct != 'C' && class[0] == 'z')
		SET_TRANS(ans);
	if (n > 0) {
	SEXP perm = PROTECT(Rf_allocVector(INTSXP, n)),
		y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
	int *pperm = INTEGER(perm), info;
	if (class[2] != 'p') {
	int lwork = -1;
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y), tmp, *work;
	memset(py, 0, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	F77_CALL(zlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
	if (ct == 'C') {
	F77_CALL(zhetrf)(&ul, &n, py, &n, pperm, &tmp, &lwork, &info FCONE);
	lwork = (int) tmp.r;
	work = (Rcomplex *) R_alloc((size_t) lwork, sizeof(Rcomplex));
	F77_CALL(zhetrf)(&ul, &n, py, &n, pperm, work, &lwork, &info FCONE);
	ERROR_LAPACK_2(zhetrf, info, warn, D);
	} else {
	F77_CALL(zsytrf)(&ul, &n, py, &n, pperm, &tmp, &lwork, &info FCONE);
	lwork = (int) tmp.r;
	work = (Rcomplex *) R_alloc((size_t) lwork, sizeof(Rcomplex));
	F77_CALL(zsytrf)(&ul, &n, py, &n, pperm, work, &lwork, &info FCONE);
	ERROR_LAPACK_2(zsytrf, info, warn, D);
	}
	} else {
	double *px = REAL(x), *py = REAL(y), tmp, *work;
	memset(py, 0, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
	F77_CALL(dsytrf)(&ul, &n, py, &n, pperm, &tmp, &lwork, &info FCONE);
	lwork = (int) tmp;
	work = (double *) R_alloc((size_t) lwork, sizeof(double));
	F77_CALL(dsytrf)(&ul, &n, py, &n, pperm, work, &lwork, &info FCONE);
	ERROR_LAPACK_2(dsytrf, info, warn, D);
	}
	} else {
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
	memcpy(py, px, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	if (ct == 'C') {
	F77_CALL(zhptrf)(&ul, &n, py, pperm, &info FCONE);
	ERROR_LAPACK_2(zhptrf, info, warn, D);
	} else {
	F77_CALL(zsptrf)(&ul, &n, py, pperm, &info FCONE);
	ERROR_LAPACK_2(zsptrf, info, warn, D);
	}
	} else {
	double *px = REAL(x), *py = REAL(y);
	memcpy(py, px, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dsptrf)(&ul, &n, py, pperm, &info FCONE);
	ERROR_LAPACK_2(dsptrf, info, warn, D);
	}
	}
	SET_SLOT(ans, Matrix_permSym, perm);
	SET_SLOT(ans, Matrix_xSym, y);
	UNPROTECT(2); /* y, perm */
	}
	UNPROTECT(3); /* ans, x, obj */
	return ans;
}

SEXP R_dense_bunchkaufman(SEXP s_obj, SEXP s_warn,
                          SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);
	char ul = '\0', ct = '\0';
	if (s_uplo != R_NilValue)
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);
	int cache = (class[1] == 's') &&
		(class[0] == 'z' || class[0] == 'd');
	const char *nm = (class[1] == 's' && class[2] == 'p')
		? "denseBunchKaufman-"
		: "denseBunchKaufman+";
	SEXP ans = (cache) ? get_factor(s_obj, nm) : R_NilValue;
	if (ans == R_NilValue) {
		int warn = Rf_asInteger(s_warn);
		ans = dense_bunchkaufman(s_obj, class, warn, ul, ct);
		if (cache) {
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
		}
	}
	return ans;
}
