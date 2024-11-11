#include "Lapack-etc.h"
#include "Mdefines.h"

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP dense_as_general(SEXP, const char *, int);

SEXP dense_schur(SEXP obj, const char *class, int warn, int vectors)
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
	if ((class[1] == 's' && class[0] == 'z' && TRANS(obj) != 'C') ||
	    (class[1] == 't')) {
		REPROTECT(obj = dense_as_general(obj, class, 1), pid);
		class = Matrix_class(obj, valid_dense, 6, __func__);
	}
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseSchur";
	cl[0] = class[0];
	SEXP ans = PROTECT(newObject(cl));
	SET_DIM(ans, n, n);
	SET_DIMNAMES(ans, -(class[1] == 's'), DIMNAMES(obj, 0));
	if (n > 0) {
	SEXP y, v, w;
	int info;
	if (class[1] == 'g') {
	PROTECT(y = Rf_allocVector(TYPEOF(x), XLENGTH(x)));
	PROTECT(v = Rf_allocVector(TYPEOF(x), (vectors) ? XLENGTH(x) : 0));
	PROTECT(w = Rf_allocVector(TYPEOF(x), n));
	int lwork = -1, zero = 0;
	double *rwork = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y),
		*pv = COMPLEX(v), *pw = COMPLEX(w), tmp, *work = &tmp;
	memcpy(py, px, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	F77_CALL(zgees)((vectors) ? "V" : "N", "N", NULL,
					&n, py, &n, &zero, pw,        pv, &n,
					work, &lwork, rwork, NULL, &info FCONE FCONE);
	lwork = (int) tmp.r;
	work = (Rcomplex *) R_alloc((size_t) lwork, sizeof(Rcomplex));
	F77_CALL(zgees)((vectors) ? "V" : "N", "N", NULL,
					&n, py, &n, &zero, pw,        pv, &n,
					work, &lwork, rwork, NULL, &info FCONE FCONE);
	ERROR_LAPACK_5(zgees, info, warn);
	} else {
	double *px = REAL(x), *py = REAL(y),
		*pv = REAL(v), *pw = REAL(w), tmp, *work = &tmp;
	memcpy(py, px, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dgees)((vectors) ? "V" : "N", "N", NULL,
					&n, py, &n, &zero, pw, rwork, pv, &n,
					work, &lwork,        NULL, &info FCONE FCONE);
	lwork = (int) tmp;
	work = (double *) R_alloc((size_t) lwork, sizeof(double));
	F77_CALL(dgees)((vectors) ? "V" : "N", "N", NULL,
					&n, py, &n, &zero, pw, rwork, pv, &n,
					work, &lwork,        NULL, &info FCONE FCONE);
	ERROR_LAPACK_5(dgees, info, warn);
	for (int j = 0; j < n; ++j) {
		if (fabs(rwork[j]) > 10.0 * DBL_EPSILON * fabs(pw[j])) {
			SEXP w_ = Rf_allocVector(CPLXSXP, n);
			Rcomplex *pw_ = COMPLEX(w_);
			for (j = 0; j < n; ++j) {
				pw_[j].r =    pw[j];
				pw_[j].i = rwork[j];
			}
			UNPROTECT(1); /* w */
			PROTECT(w = w_);
			break;
		}
	}
	}
	} else if (class[2] != 'p') {
	PROTECT(y = Rf_allocVector(TYPEOF(x), 0));
	PROTECT(v = Rf_allocVector(TYPEOF(x), XLENGTH(x)));
	PROTECT(w = Rf_allocVector(REALSXP, n));
	char ul = UPLO(obj);
	int lwork = -1;
	double *pw = REAL(w);
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *pv = COMPLEX(v), tmp, *work = &tmp;
	double *rwork = (double *) R_alloc((size_t) n * 3, sizeof(double));
	memcpy(pv, px, sizeof(Rcomplex) * (size_t) XLENGTH(v));
	F77_CALL(zheev)((vectors) ? "V" : "N", &ul,
					&n, pv, &n, pw, work, &lwork, rwork, &info FCONE FCONE);
	lwork = (int) tmp.r;
	work = (Rcomplex *) R_alloc((size_t) lwork, sizeof(Rcomplex));
	F77_CALL(zheev)((vectors) ? "V" : "N", &ul,
					&n, pv, &n, pw, work, &lwork, rwork, &info FCONE FCONE);
	ERROR_LAPACK_5(zheev, info, warn);
	} else {
	double *px = REAL(x), *pv = REAL(v), tmp, *work = &tmp;
	memcpy(pv, px, sizeof(double) * (size_t) XLENGTH(v));
	F77_CALL(dsyev)((vectors) ? "V" : "N", &ul,
					&n, pv, &n, pw, work, &lwork,        &info FCONE FCONE);
	lwork = (int) tmp;
	work = (double *) R_alloc((size_t) lwork, sizeof(double));
	F77_CALL(dsyev)((vectors) ? "V" : "N", &ul,
					&n, pv, &n, pw, work, &lwork,        &info FCONE FCONE);
	ERROR_LAPACK_5(dsyev, info, warn);
	}
	} else {
	PROTECT(y = Rf_allocVector(TYPEOF(x), XLENGTH(x)));
	PROTECT(v = Rf_allocVector(TYPEOF(x), (vectors) ? (R_xlen_t) n * n : 0));
	PROTECT(w = Rf_allocVector(REALSXP, n));
	char ul = UPLO(obj);
	double *pw = REAL(w);
	if (TYPEOF(x) == CPLXSXP) {
	Rcomplex *px = COMPLEX(x), *py = COMPLEX(y), *pv = COMPLEX(v),
		*work = (Rcomplex *) R_alloc((size_t) n * 2, sizeof(Rcomplex));
	double *rwork = (double *) R_alloc((size_t) n * 3, sizeof(double));
	memcpy(py, px, sizeof(Rcomplex) * (size_t) XLENGTH(y));
	F77_CALL(zhpev)((vectors) ? "V" : "N", &ul,
					&n, py, pw, pv, &n, work, rwork, &info FCONE FCONE);
	ERROR_LAPACK_5(zhpev, info, warn);
	} else {
	double *px = REAL(x), *py = REAL(y), *pv = REAL(v),
		*work = (double *) R_alloc((size_t) n * 3, sizeof(double));
	memcpy(py, px, sizeof(double) * (size_t) XLENGTH(y));
	F77_CALL(dspev)((vectors) ? "V" : "N", &ul,
					&n, py, pw, pv, &n, work,        &info FCONE FCONE);
	ERROR_LAPACK_5(dspev, info, warn);
	}
	}
	if (class[1] != 's')
	SET_SLOT(ans, Matrix_xSym, y);
	if (vectors)
	SET_SLOT(ans, Matrix_vectorsSym, v);
	SET_SLOT(ans, Matrix_valuesSym, w);
	UNPROTECT(3); /* w, v, y */
	}
	UNPROTECT(3); /* ans, x, obj */
	return ans;
}

SEXP R_dense_schur(SEXP s_obj, SEXP s_warn, SEXP s_vectors)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);
	int vectors = Rf_asLogical(s_vectors);
	int cache = vectors &&
		(class[1] == 'g' || class[1] == 's') &&
		(class[0] == 'z' || class[0] == 'd');
	const char *nm = "denseSchur";
	SEXP ans = (cache) ? get_factor(s_obj, nm) : R_NilValue;
	if (ans == R_NilValue) {
		int warn = Rf_asInteger(s_warn);
		ans = dense_schur(s_obj, class, warn, vectors);
		if (cache) {
			PROTECT(ans);
			set_factor(s_obj, nm, ans);
			UNPROTECT(1);
		}
	}
	return ans;
}
