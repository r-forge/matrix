#include <math.h> /* fabs */
#include "Lapack-etc.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "M5.h"

/* defined in ./attrib.c : */
SEXP get_factor(SEXP, const char *);
void set_factor(SEXP, const char *, SEXP);

static
SEXP geMatrix_scf_(SEXP obj, int warn, int vectors)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseSchur";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP scf = PROTECT(newObject(cl));
	SET_DIM(scf, n, n);
	SET_DIMNAMES(scf, 0, DIMNAMES(obj, 0));
	if (n > 0) {
		SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x))),
			v = PROTECT(Rf_allocVector(TYPEOF(x), (vectors) ? XLENGTH(x) : 0)),
			w = PROTECT(Rf_allocVector(TYPEOF(x), n));
		int info, lwork = -1, zero = 0;
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
		int j;
		for (j = 0; j < n; ++j) {
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
		SET_SLOT(scf, Matrix_xSym, y);
		SET_SLOT(scf, Matrix_vectorsSym, v);
		SET_SLOT(scf, Matrix_valuesSym, w);
		UNPROTECT(3); /* w, v, y */
	}
	UNPROTECT(2); /* scf, x */
	return scf;
}

static
SEXP syMatrix_scf_(SEXP obj, int warn, int vectors)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) == CPLXSXP && TRANS(obj) != 'C') {
		/* defined in ./coerce.c : */
		SEXP dense_as_general(SEXP, const char *, int);
		obj = dense_as_general(obj, "zsyMatrix", 1);
		PROTECT(obj);
		obj = geMatrix_scf_(obj, warn, vectors);
		UNPROTECT(1);
		return obj;
	}
	PROTECT(x);
	char cl[] = ".denseSchur";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP scf = PROTECT(newObject(cl));
	int n = DIM(obj)[1];
	SET_DIM(scf, n, n);
	SET_DIMNAMES(scf, -1, DIMNAMES(obj, 0));
	if (n > 0) {
		SEXP v = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x))),
			w = PROTECT(Rf_allocVector(REALSXP, n));
		char ul = UPLO(obj);
		int info, lwork = -1;
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
		if (vectors)
		SET_SLOT(scf, Matrix_vectorsSym, v);
		SET_SLOT(scf, Matrix_valuesSym, w);
		UNPROTECT(2); /* w, v */
	}
	UNPROTECT(2); /* scf, x */
	return scf;
}

static
SEXP spMatrix_scf_(SEXP obj, int warn, int vectors)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) == CPLXSXP && TRANS(obj) != 'C') {
		/* defined in ./coerce.c : */
		SEXP dense_as_general(SEXP, const char *, int);
		obj = dense_as_general(obj, "zspMatrix", 1);
		PROTECT(obj);
		obj = geMatrix_scf_(obj, warn, vectors);
		UNPROTECT(1);
		return obj;
	}
	PROTECT(x);
	char cl[] = ".denseSchur";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP scf = PROTECT(newObject(cl));
	int n = DIM(obj)[1];
	SET_DIM(scf, n, n);
	SET_DIMNAMES(scf, -1, DIMNAMES(obj, 0));
	if (n > 0) {
		SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x))),
			v = PROTECT(Rf_allocVector(TYPEOF(x), (vectors) ? (R_xlen_t) n * n : 0)),
			w = PROTECT(Rf_allocVector(REALSXP, n));
		char ul = UPLO(obj);
		int info;
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
		if (vectors)
		SET_SLOT(scf, Matrix_vectorsSym, v);
		SET_SLOT(scf, Matrix_valuesSym, w);
		UNPROTECT(3); /* w, v, y */
	}
	UNPROTECT(2); /* scf, x */
	return scf;
}

static
SEXP geMatrix_trf_(SEXP obj, int warn)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseLU";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	SET_DIM(trf, m, n);
	SET_DIMNAMES(trf, 0, DIMNAMES(obj, 0));
	if (r > 0) {
		SEXP perm = PROTECT(Rf_allocVector(INTSXP, r)),
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
		SET_SLOT(trf, Matrix_permSym, perm);
		SET_SLOT(trf, Matrix_xSym, y);
		UNPROTECT(2); /* y, perm */
	}
	UNPROTECT(2); /* trf, x */
	return trf;
}

static
SEXP syMatrix_trf_(SEXP obj, int warn)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseBunchKaufman";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	int n = DIM(obj)[1];
	SET_DIM(trf, n, n);
	SET_DIMNAMES(trf, -1, DIMNAMES(obj, 0));
	char ul = UPLO(obj),
		ct = (TYPEOF(x) == CPLXSXP) ? TRANS(obj) : 'C';
	if (ul != 'U')
		SET_UPLO(trf);
	if (ct != 'C')
		SET_TRANS(trf);
	if (n > 0) {
		SEXP perm = PROTECT(Rf_allocVector(INTSXP, n)),
			y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
		int *pperm = INTEGER(perm), info, lwork = -1;
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
		SET_SLOT(trf, Matrix_permSym, perm);
		SET_SLOT(trf, Matrix_xSym, y);
		UNPROTECT(2); /* y, perm */
	}
	UNPROTECT(2); /* trf, x */
	return trf;
}

static
SEXP spMatrix_trf_(SEXP obj, int warn)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseBunchKaufman";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	int n = DIM(obj)[1];
	SET_DIM(trf, n, n);
	SET_DIMNAMES(trf, -1, DIMNAMES(obj, 0));
	char ul = UPLO(obj),
		ct = (TYPEOF(x) == CPLXSXP) ? TRANS(obj) : 'C';
	if (ul != 'U')
		SET_UPLO(trf);
	if (ct != 'C')
		SET_TRANS(trf);
	if (n > 0) {
		SEXP perm = PROTECT(Rf_allocVector(INTSXP, n)),
			y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
		int *pperm = INTEGER(perm), info;
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
		SET_SLOT(trf, Matrix_permSym, perm);
		SET_SLOT(trf, Matrix_xSym, y);
		UNPROTECT(2); /* y, perm */
	}
	UNPROTECT(2); /* trf, x */
	return trf;
}

static
SEXP poMatrix_trf_(SEXP obj, int warn, int pivot, double tol)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseCholesky";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	int n = DIM(obj)[1];
	SET_DIM(trf, n, n);
	SET_DIMNAMES(trf, -1, DIMNAMES(obj, 0));
	char ul = UPLO(obj);
	if (ul != 'U')
		SET_UPLO(trf);
	if (n > 0) {
		SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
		int info;
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
		F77_CALL(zpstrf)(&ul, &n, py, &n, pperm, &rank, &tol, work, &info FCONE);
		ERROR_LAPACK_4(zpstrf, info, warn, rank);
		if (info > 0) {
			int j, d = n - rank;
			py += (R_xlen_t) rank * n + rank;
			for (j = rank; j < n; ++j) {
				memset(py, 0, sizeof(Rcomplex) * (size_t) d);
				py += n;
			}
		}
		SET_SLOT(trf, Matrix_permSym, perm);
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
		F77_CALL(dpstrf)(&ul, &n, py, &n, pperm, &rank, &tol, work, &info FCONE);
		ERROR_LAPACK_4(dpstrf, info, warn, rank);
		if (info > 0) {
			int j, d = n - rank;
			py += (R_xlen_t) rank * n + rank;
			for (j = rank; j < n; ++j) {
				memset(py, 0, sizeof(double) * (size_t) d);
				py += n;
			}
		}
		SET_SLOT(trf, Matrix_permSym, perm);
		UNPROTECT(1); /* perm */
		}
		}
		SET_SLOT(trf, Matrix_xSym, y);
		UNPROTECT(1); /* y */
	}
	UNPROTECT(2); /* trf, x */
	return trf;
}

static
SEXP ppMatrix_trf_(SEXP obj, int warn)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	char cl[] = ".denseCholesky";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	int n = DIM(obj)[1];
	SET_DIM(trf, n, n);
	SET_DIMNAMES(trf, -1, DIMNAMES(obj, 0));
	char ul = UPLO(obj);
	if (ul != 'U')
		SET_UPLO(trf);
	if (n > 0) {
		SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(x)));
		int info;
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
		SET_SLOT(trf, Matrix_xSym, y);
		UNPROTECT(1); /* y */
	}
	UNPROTECT(2); /* trf, x */
	return trf;
}

SEXP geMatrix_scf(SEXP s_obj, SEXP s_warn, SEXP s_vectors)
{
	int vectors = Rf_asLogical(s_vectors);
	const char *nm = "denseSchur";
	SEXP scf = (vectors) ? get_factor(s_obj, nm) : R_NilValue;
	if (scf == R_NilValue) {
		scf = geMatrix_scf_(s_obj, Rf_asInteger(s_warn), vectors);
		if (vectors) {
		PROTECT(scf);
		set_factor(s_obj, nm, scf);
		UNPROTECT(1);
		}
	}
	return scf;
}

SEXP syMatrix_scf(SEXP s_obj, SEXP s_warn, SEXP s_vectors)
{
	int vectors = Rf_asLogical(s_vectors);
	const char *nm = "denseSchur";
	SEXP scf = (vectors) ? get_factor(s_obj, nm) : R_NilValue;
	if (scf == R_NilValue) {
		scf = syMatrix_scf_(s_obj, Rf_asInteger(s_warn), vectors);
		if (vectors) {
		PROTECT(scf);
		set_factor(s_obj, nm, scf);
		UNPROTECT(1);
		}
	}
	return scf;
}

SEXP spMatrix_scf(SEXP s_obj, SEXP s_warn, SEXP s_vectors)
{
	int vectors = Rf_asLogical(s_vectors);
	const char *nm = "denseSchur";
	SEXP scf = (vectors) ? get_factor(s_obj, nm) : R_NilValue;
	if (scf == R_NilValue) {
		scf = spMatrix_scf_(s_obj, Rf_asInteger(s_warn), vectors);
		if (vectors) {
		PROTECT(scf);
		set_factor(s_obj, nm, scf);
		UNPROTECT(1);
		}
	}
	return scf;
}

SEXP geMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseLU";
	SEXP trf = get_factor(s_obj, nm);
	if (trf == R_NilValue) {
		PROTECT(trf = geMatrix_trf_(s_obj, Rf_asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP syMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseBunchKaufman+";
	SEXP trf = get_factor(s_obj, nm);
	if (trf == R_NilValue) {
		PROTECT(trf = syMatrix_trf_(s_obj, Rf_asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP spMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseBunchKaufman-";
	SEXP trf = get_factor(s_obj, nm);
	if (trf == R_NilValue) {
		PROTECT(trf = spMatrix_trf_(s_obj, Rf_asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP poMatrix_trf(SEXP s_obj, SEXP s_warn, SEXP s_pivot, SEXP s_tol)
{
	int pivot = Rf_asLogical(s_pivot);
	const char *nm = (pivot) ? "denseCholesky++" : "denseCholesky+-";
	SEXP trf = get_factor(s_obj, nm);
	if (trf == R_NilValue) {
		double tol = Rf_asReal(s_tol);
		PROTECT(trf = poMatrix_trf_(s_obj, Rf_asInteger(s_warn), pivot, tol));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP ppMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseCholesky--";
	SEXP trf = get_factor(s_obj, nm);
	if (trf == R_NilValue) {
		PROTECT(trf = ppMatrix_trf_(s_obj, Rf_asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

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

static
SEXP gCMatrix_orf_(SEXP obj, int warn, int order)
{
	Matrix_cs *A = M2CXS(obj, 1);
	CXSPARSE_XTYPE_SET(A->xtype);

	if (A->m < A->n)
		Rf_error(_("QR factorization of m-by-n %s requires m >= n"),
		         ".gCMatrix");

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
	SEXP orf = PROTECT(newObject(cl));

	SET_DIM(orf, A->m, A->n);
	SET_DIMNAMES(orf, 0, DIMNAMES(obj, 0));

	SEXP V = PROTECT(CXS2M(N->L, 1, 'g')),
		R = PROTECT(CXS2M(N->U, 1, 'g'));
	SET_SLOT(orf, Matrix_VSym, V);
	SET_SLOT(orf, Matrix_RSym, R);
	UNPROTECT(2); /* R, V */

	SEXP beta = PROTECT(Rf_allocVector(REALSXP, A->n));
	memcpy(REAL(beta), N->B, sizeof(double) * (size_t) A->n);
	SET_SLOT(orf, Matrix_betaSym, beta);
	UNPROTECT(1); /* beta */

	SEXP p = PROTECT(Rf_allocVector(INTSXP, S->m2));
	memcpy(INTEGER(p), P, sizeof(int) * (size_t) S->m2);
	SET_SLOT(orf, Matrix_pSym, p);
	UNPROTECT(1); /* p */

	if (order > 0) {
	SEXP q = PROTECT(Rf_allocVector(INTSXP, A->n));
	memcpy(INTEGER(q), S->q, sizeof(int) * (size_t) A->n);
	SET_SLOT(orf, Matrix_qSym, q);
	UNPROTECT(1); /* q */
	}

	DO_FREE(T, S, N, P);
	UNPROTECT(1); /* orf */
	return orf;

oom:
	DO_FREE(T, S, N, P);
	if (warn > 1)
		Rf_error  (_("QR factorization of %s failed: out of memory"),
		           ".gCMatrix");
	else if (warn > 0)
		Rf_warning(_("QR factorization of %s failed: out of memory"),
		           ".gCMatrix");
	return R_NilValue;
}

static
SEXP gCMatrix_trf_(SEXP obj, int warn, int order, double tol)
{
	Matrix_cs *A = M2CXS(obj, 1);
	CXSPARSE_XTYPE_SET(A->xtype);

	if (A->m != A->n)
		Rf_error(_("LU factorization of m-by-n %s requires m == n"),
		         ".gCMatrix");

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
	SEXP trf = PROTECT(newObject(cl));

	SET_DIM(trf, A->m, A->n);
	SET_DIMNAMES(trf, 0, DIMNAMES(obj, 0));

	SEXP L = PROTECT(CXS2M(N->L, 1, 't')),
		U = PROTECT(CXS2M(N->U, 1, 't'));
	SET_UPLO(L);
	SET_SLOT(trf, Matrix_LSym, L);
	SET_SLOT(trf, Matrix_USym, U);
	UNPROTECT(2); /* U, L */

	SEXP p = PROTECT(Rf_allocVector(INTSXP, A->m));
	memcpy(INTEGER(p), P, sizeof(int) * (size_t) A->m);
	SET_SLOT(trf, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order > 0) {
	SEXP q = PROTECT(Rf_allocVector(INTSXP, A->n));
	memcpy(INTEGER(q), S->q, sizeof(int) * (size_t) A->n);
	SET_SLOT(trf, Matrix_qSym, q);
	UNPROTECT(1); /* q */
	}

	DO_FREE(T, S, N, P);
	UNPROTECT(1); /* trf */
	return trf;

oom:
	DO_FREE(T, S, N, P);
	if (warn > 1)
		Rf_error  (_("LU factorization of %s failed: out of memory or near-singular"),
		           ".gCMatrix");
	else if (warn > 0)
		Rf_warning(_("LU factorization of %s failed: out of memory or near-singular"),
		           ".gCMatrix");
	return R_NilValue;
}

#undef DO_FREE
#undef DO_SORT

static
SEXP pCMatrix_trf_(SEXP obj, SEXP trf,
                   int warn, int order, int *ll, int *super, Rcomplex beta)
{
	cholmod_sparse *A =                              M2CHS(obj, 1);
	cholmod_factor *L = (trf == R_NilValue) ? NULL : M2CHF(trf, 1);
	double betaRI[2]; betaRI[0] = beta.r; betaRI[1] = beta.i;

	A->stype = (UPLO(obj) == 'U') ? 1 : -1;

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

	cholmod_factorize_p(A, betaRI, NULL, 0, L, &c);
	cholmod_defaults(&c);

#define PCTRF_FINISH(_OBJ_, _WARN_) \
	do { \
		SEXP dimnames = PROTECT(DIMNAMES(_OBJ_, 0)); \
		PROTECT(trf = CHF2M(L, 1)); \
		cholmod_free_factor(&L, &c); \
		if (TYPEOF(trf) == CHARSXP) { \
			if (_WARN_ > 1) \
				Rf_error  ("%s", CHAR(trf)); \
			else if (_WARN_ > 0) \
				Rf_warning("%s", CHAR(trf)); \
			UNPROTECT(2); \
			return R_NilValue; \
		} \
		SET_DIMNAMES(trf, -1, dimnames); \
		UNPROTECT(2); \
	} while (0)

	PCTRF_FINISH(obj, warn);
	return trf;
}

SEXP gCMatrix_orf(SEXP s_obj, SEXP s_warn, SEXP s_order)
{
	int order = Rf_asInteger(s_order);
	if (order < 0 || order > 3)
		order = 0;
	const char *nm = (order > 0) ? "sparseQR+" : "sparseQR-";
	SEXP orf = get_factor(s_obj, nm);
	if (orf == R_NilValue) {
		orf = gCMatrix_orf_(s_obj, Rf_asInteger(s_warn), order);
		if (orf != R_NilValue) {
		PROTECT(orf);
		set_factor(s_obj, nm, orf);
		UNPROTECT(1);
		}
	}
	return orf;
}

SEXP gCMatrix_trf(SEXP s_obj, SEXP s_warn, SEXP s_order, SEXP s_tol)
{
	double tol = Rf_asReal(s_tol);
	if (ISNAN(tol))
		Rf_error(_("'%s' is not a number"), "tol");
	int order = Rf_asInteger(s_order);
	if (order == NA_INTEGER)
		order = (tol == 1.0) ? 2 : 1;
	else if (order < 0 || order > 3)
		order = 0;
	const char *nm = (order > 0) ? "sparseLU+" : "sparseLU-";
	SEXP trf = get_factor(s_obj, nm);
	if (trf == R_NilValue) {
		trf = gCMatrix_trf_(s_obj, Rf_asInteger(s_warn), order, tol);
		if (trf != R_NilValue) {
		PROTECT(trf);
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
		}
	}
	return trf;
}

SEXP pCMatrix_trf(SEXP s_obj, SEXP s_warn, SEXP s_order,
                  SEXP s_ll, SEXP s_super, SEXP s_beta)
{
	int warn = Rf_asInteger(s_warn), order = Rf_asInteger(s_order),
		ll = Rf_asLogical(s_ll), super = Rf_asLogical(s_super);
	Rcomplex beta = Rf_asComplex(s_beta);
	if (order < 0 || order > 1)
		order = 0;
	if (!R_FINITE(beta.r) || !R_FINITE(beta.i))
		Rf_error(_("'%s' is not a number or not finite"), "beta");
	SEXP trf = R_NilValue;
	char nm[] = "..........Cholesky.";
	nm[18] = (order > 0) ? '+' : '-';
	if (super == NA_LOGICAL || super == 0) {
		memcpy(nm, "simplicial", 10);
		trf = get_factor(s_obj, nm);
		if (trf != R_NilValue) super = 0;
	}
	if (trf == R_NilValue && (super == NA_LOGICAL || super != 0)) {
		memcpy(nm, "supernodal", 10);
		trf = get_factor(s_obj, nm);
		if (trf != R_NilValue) super = 1;
	}
	if (beta.r != 0.0 || beta.i != 0.0 || trf == R_NilValue ||
	    (super == 0 && ll != 0)) {
		if (beta.r != 0.0 || beta.i != 0.0) {
			PROTECT(trf);
			trf = pCMatrix_trf_(s_obj, trf, warn, order, &ll, &super, beta);
			UNPROTECT(1);
		} else {
			if (trf == R_NilValue) {
			int zz_ = 0;
			trf = pCMatrix_trf_(s_obj, trf, warn, order, &zz_, &super, beta);
			if (trf != R_NilValue) {
			memcpy(nm, (super == 0) ? "simplicial" : "supernodal", 10);
			PROTECT(trf);
			set_factor(s_obj, nm, trf);
			UNPROTECT(1);
			}
			}
			if (trf != R_NilValue && super == 0 && ll != 0) {
			PROTECT(trf);
			cholmod_factor *L = M2CHF(trf, 1);
			L = cholmod_copy_factor(L, &c);
			cholmod_change_factor(L->xtype, 1, 0, 1, 1, L, &c);
			PCTRF_FINISH(s_obj, warn);
			UNPROTECT(1);
			}
		}
	}
	return trf;
}
