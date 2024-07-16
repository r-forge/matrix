#include <math.h> /* fabs */
#include "Lapack-etc.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "factor.h"

/* defined in ./attrib.c : */
SEXP get_factor(SEXP, const char *);
void set_factor(SEXP, const char *, SEXP);

static
SEXP geMatrix_scf_(SEXP obj, int warn, int vectors)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int *pdim = INTEGER(dim), n = pdim[1];
	if (pdim[0] != n)
		error(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");
	char cl[] = ".denseSchur";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP scf = PROTECT(newObject(cl));
	SET_SLOT(scf, Matrix_DimSym, dim);
	SET_SLOT(scf, Matrix_DimNamesSym, dimnames);
	if (n > 0) {
		SEXP y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x))),
			v = PROTECT(allocVector(TYPEOF(x), (vectors) ? XLENGTH(x) : 0)),
			w = PROTECT(allocVector(TYPEOF(x), n));
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
				SEXP w_ = allocVector(CPLXSXP, n);
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
	UNPROTECT(4); /* scf, x, dimnames, dim */
	return scf;
}

static
SEXP syMatrix_scf_(SEXP obj, int warn, int vectors)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) == CPLXSXP) {
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		char ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C') {
			/* defined in ./coerce.c : */
			SEXP dense_as_general(SEXP, const char *, int);
			PROTECT(obj = dense_as_general(obj, "zsyMatrix", 1));
			obj = geMatrix_scf_(obj, warn, vectors);
			UNPROTECT(1);
			return obj;
		}
	}
	PROTECT(x);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	int n = INTEGER(dim)[1];
	char cl[] = ".denseSchur";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP scf = PROTECT(newObject(cl));
	SET_SLOT(scf, Matrix_DimSym, dim);
	set_symmetrized_DimNames(scf, dimnames, -1);
	if (n > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
			v = PROTECT(allocVector(TYPEOF(x), XLENGTH(x))),
			w = PROTECT(allocVector(REALSXP, n));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
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
		UNPROTECT(3); /* w, v, uplo */
	}
	UNPROTECT(4); /* scf, dimnames, dim, x */
	return scf;
}

static
SEXP spMatrix_scf_(SEXP obj, int warn, int vectors)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) == CPLXSXP) {
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		char ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C') {
			/* defined in ./coerce.c : */
			SEXP dense_as_general(SEXP, const char *, int);
			PROTECT(obj = dense_as_general(obj, "zspMatrix", 1));
			obj = geMatrix_scf_(obj, warn, vectors);
			UNPROTECT(1);
			return obj;
		}
	}
	PROTECT(x);
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	int n = INTEGER(dim)[1];
	char cl[] = ".denseSchur";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP scf = PROTECT(newObject(cl));
	SET_SLOT(scf, Matrix_DimSym, dim);
	set_symmetrized_DimNames(scf, dimnames, -1);
	if (n > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x))),
			v = PROTECT(allocVector(TYPEOF(x), (vectors) ? (R_xlen_t) n * n : 0)),
			w = PROTECT(allocVector(REALSXP, n));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
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
		UNPROTECT(4); /* w, v, y, uplo */
	}
	UNPROTECT(4); /* scf, dimnames, dim, x */
	return scf;
}

static
SEXP geMatrix_trf_(SEXP obj, int warn)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	char cl[] = ".denseLU";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	SET_SLOT(trf, Matrix_DimSym, dim);
	SET_SLOT(trf, Matrix_DimNamesSym, dimnames);
	if (r > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, r)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
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
	UNPROTECT(4); /* trf, x, dimnames, dim */
	return trf;
}

static
SEXP syMatrix_trf_(SEXP obj, int warn)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int n = INTEGER(dim)[1];
	char ul = CHAR(STRING_ELT(uplo, 0))[0], ct = 'C';
	char cl[] = ".denseBunchKaufman";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	SET_SLOT(trf, Matrix_DimSym, dim);
	set_symmetrized_DimNames(trf, dimnames, -1);
	SET_SLOT(trf, Matrix_uploSym, uplo);
	if (TYPEOF(x) == CPLXSXP) {
		SEXP trans = PROTECT(GET_SLOT(obj, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(trf, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (n > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, n)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
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
	UNPROTECT(5); /* trf, x, uplo, dimnames, dim */
	return trf;
}

static
SEXP spMatrix_trf_(SEXP obj, int warn)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int n = INTEGER(dim)[1];
	char ul = CHAR(STRING_ELT(uplo, 0))[0], ct = 'C';
	char cl[] = ".denseBunchKaufman";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	SET_SLOT(trf, Matrix_DimSym, dim);
	set_symmetrized_DimNames(trf, dimnames, -1);
	SET_SLOT(trf, Matrix_uploSym, uplo);
	if (TYPEOF(x) == CPLXSXP) {
		SEXP trans = PROTECT(GET_SLOT(obj, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(trf, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (n > 0) {
		SEXP perm = PROTECT(allocVector(INTSXP, n)),
			y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
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
	UNPROTECT(5); /* trf, x, uplo, dimnames, dim */
	return trf;
}

static
SEXP poMatrix_trf_(SEXP obj, int warn, int pivot, double tol)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int n = INTEGER(dim)[1];
	char ul = CHAR(STRING_ELT(uplo, 0))[0];
	char cl[] = ".denseCholesky";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	SET_SLOT(trf, Matrix_DimSym, dim);
	set_symmetrized_DimNames(trf, dimnames, -1);
	SET_SLOT(trf, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
		int info;
		if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), *py = COMPLEX(y);
		memset(py, 0, sizeof(Rcomplex) * (size_t) XLENGTH(y));
		F77_CALL(zlacpy)(&ul, &n, &n, px, &n, py, &n FCONE);
		if (!pivot) {
		F77_CALL(zpotrf)(&ul, &n, py, &n, &info FCONE);
		ERROR_LAPACK_3(zpotrf, info, warn);
		} else {
		SEXP perm = PROTECT(allocVector(INTSXP, n));
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
		SEXP perm = PROTECT(allocVector(INTSXP, n));
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
	UNPROTECT(5); /* trf, x, uplo, dimnames, dim */
	return trf;
}

static
SEXP ppMatrix_trf_(SEXP obj, int warn)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int n = INTEGER(dim)[1];
	char ul = CHAR(STRING_ELT(uplo, 0))[0];
	char cl[] = ".denseCholesky";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));
	SET_SLOT(trf, Matrix_DimSym, dim);
	set_symmetrized_DimNames(trf, dimnames, -1);
	SET_SLOT(trf, Matrix_uploSym, uplo);
	if (n > 0) {
		SEXP y = PROTECT(allocVector(TYPEOF(x), XLENGTH(x)));
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
	UNPROTECT(5); /* trf, x, uplo, dimnames, dim */
	return trf;
}

SEXP geMatrix_scf(SEXP s_obj, SEXP s_warn, SEXP s_vectors)
{
	int vectors = asLogical(s_vectors);
	const char *nm = "denseSchur";
	SEXP scf = (vectors) ? get_factor(s_obj, nm) : R_NilValue;
	if (isNull(scf)) {
		scf = geMatrix_scf_(s_obj, asInteger(s_warn), vectors);
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
	int vectors = asLogical(s_vectors);
	const char *nm = "denseSchur";
	SEXP scf = (vectors) ? get_factor(s_obj, nm) : R_NilValue;
	if (isNull(scf)) {
		scf = syMatrix_scf_(s_obj, asInteger(s_warn), vectors);
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
	int vectors = asLogical(s_vectors);
	const char *nm = "denseSchur";
	SEXP scf = (vectors) ? get_factor(s_obj, nm) : R_NilValue;
	if (isNull(scf)) {
		scf = spMatrix_scf_(s_obj, asInteger(s_warn), vectors);
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
	if (isNull(trf)) {
		PROTECT(trf = geMatrix_trf_(s_obj, asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP syMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseBunchKaufman+";
	SEXP trf = get_factor(s_obj, nm);
	if (isNull(trf)) {
		PROTECT(trf = syMatrix_trf_(s_obj, asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP spMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseBunchKaufman-";
	SEXP trf = get_factor(s_obj, nm);
	if (isNull(trf)) {
		PROTECT(trf = spMatrix_trf_(s_obj, asInteger(s_warn)));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP poMatrix_trf(SEXP s_obj, SEXP s_warn, SEXP s_pivot, SEXP s_tol)
{
	int pivot = asLogical(s_pivot);
	const char *nm = (pivot) ? "denseCholesky++" : "denseCholesky+-";
	SEXP trf = get_factor(s_obj, nm);
	if (isNull(trf)) {
		double tol = asReal(s_tol);
		PROTECT(trf = poMatrix_trf_(s_obj, asInteger(s_warn), pivot, tol));
		set_factor(s_obj, nm, trf);
		UNPROTECT(1);
	}
	return trf;
}

SEXP ppMatrix_trf(SEXP s_obj, SEXP s_warn)
{
	const char *nm = "denseCholesky--";
	SEXP trf = get_factor(s_obj, nm);
	if (isNull(trf)) {
		PROTECT(trf = ppMatrix_trf_(s_obj, asInteger(s_warn)));
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
		error(_("QR factorization of m-by-n %s requires m >= n"),
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
	cl[0] = (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX) ? 'z' : 'd';
	SEXP orf = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(orf, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(orf, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP V = PROTECT(CXS2M(N->L, 1, 'g')),
		R = PROTECT(CXS2M(N->U, 1, 'g'));
	SET_SLOT(orf, Matrix_VSym, V);
	SET_SLOT(orf, Matrix_RSym, R);
	UNPROTECT(2); /* R, V */

	SEXP beta = PROTECT(allocVector(REALSXP, A->n));
	memcpy(REAL(beta), N->B, sizeof(double) * (size_t) A->n);
	SET_SLOT(orf, Matrix_betaSym, beta);
	UNPROTECT(1); /* beta */

	SEXP p = PROTECT(allocVector(INTSXP, S->m2));
	memcpy(INTEGER(p), P, sizeof(int) * (size_t) S->m2);
	SET_SLOT(orf, Matrix_pSym, p);
	UNPROTECT(1); /* p */

	if (order > 0) {
	SEXP q = PROTECT(allocVector(INTSXP, A->n));
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
		error  (_("QR factorization of %s failed: out of memory"),
		        ".gCMatrix");
	else if (warn > 0)
		warning(_("QR factorization of %s failed: out of memory"),
		        ".gCMatrix");
	return R_NilValue;
}

static
SEXP gCMatrix_trf_(SEXP obj, int warn, int order, double tol)
{
	Matrix_cs *A = M2CXS(obj, 1);
	CXSPARSE_XTYPE_SET(A->xtype);

	if (A->m != A->n)
		error(_("LU factorization of m-by-n %s requires m == n"),
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
	cl[0] = (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX) ? 'z' : 'd';
	SEXP trf = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(trf, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(trf, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP L = PROTECT(CXS2M(N->L, 1, 't')),
		U = PROTECT(CXS2M(N->U, 1, 't')),
		uplo = PROTECT(mkString("L"));
	SET_SLOT(L, Matrix_uploSym, uplo);
	SET_SLOT(trf, Matrix_LSym, L);
	SET_SLOT(trf, Matrix_USym, U);
	UNPROTECT(3); /* uplo, U, L */

	SEXP p = PROTECT(allocVector(INTSXP, A->m));
	memcpy(INTEGER(p), P, sizeof(int) * (size_t) A->m);
	SET_SLOT(trf, Matrix_pSym, p);
	UNPROTECT(1); /* p */
	if (order > 0) {
	SEXP q = PROTECT(allocVector(INTSXP, A->n));
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
		error  (_("LU factorization of %s failed: out of memory or near-singular"),
		        ".gCMatrix");
	else if (warn > 0)
		warning(_("LU factorization of %s failed: out of memory or near-singu"),
		        ".gCMatrix");
	return R_NilValue;
}

#undef DO_FREE
#undef DO_SORT

static
SEXP pCMatrix_trf_(SEXP obj, SEXP trf,
                   int warn, int order, int *ll, int *super, Rcomplex beta)
{
	cholmod_sparse *A =                        M2CHS(obj, 1);
	cholmod_factor *L = (isNull(trf)) ? NULL : M2CHF(trf, 1);
	double betaRI[2]; betaRI[0] = beta.r; betaRI[1] = beta.i;

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];
	A->stype = (ul == 'U') ? 1 : -1;

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
		SEXP dimnames = PROTECT(GET_SLOT(_OBJ_, Matrix_DimNamesSym)); \
		PROTECT(trf = CHF2M(L, 1)); \
		cholmod_free_factor(&L, &c); \
		if (TYPEOF(trf) == CHARSXP) { \
			if (_WARN_ > 1) \
				error  ("%s", CHAR(trf)); \
			else if (_WARN_ > 0) \
				warning("%s", CHAR(trf)); \
			UNPROTECT(2); \
			return R_NilValue; \
		} \
		set_symmetrized_DimNames(trf, dimnames, -1); \
		UNPROTECT(2); \
	} while (0)

	PCTRF_FINISH(obj, warn);
	return trf;
}

SEXP gCMatrix_orf(SEXP s_obj, SEXP s_warn, SEXP s_order)
{
	int order = asInteger(s_order);
	if (order < 0 || order > 3)
		order = 0;
	const char *nm = (order > 0) ? "sparseQR~" : "sparseQR";
	SEXP orf = get_factor(s_obj, nm);
	if (isNull(orf)) {
		orf = gCMatrix_orf_(s_obj, asInteger(s_warn), order);
		if (!isNull(orf)) {
		PROTECT(orf);
		set_factor(s_obj, nm, orf);
		UNPROTECT(1);
		}
	}
	return orf;
}

SEXP gCMatrix_trf(SEXP s_obj, SEXP s_warn, SEXP s_order, SEXP s_tol)
{
	double tol = asReal(s_tol);
	if (ISNAN(tol))
		error(_("'%s' is not a number"), "tol");
	int order = asInteger(s_order);
	if (order == NA_INTEGER)
		order = (tol == 1.0) ? 2 : 1;
	else if (order < 0 || order > 3)
		order = 0;
	const char *nm = (order > 0) ? "sparseLU~" : "sparseLU";
	SEXP trf = get_factor(s_obj, nm);
	if (isNull(trf)) {
		trf = gCMatrix_trf_(s_obj, asInteger(s_warn), order, tol);
		if (!isNull(trf)) {
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
	int warn = asInteger(s_warn), order = asInteger(s_order),
		ll = asLogical(s_ll), super = asLogical(s_super);
	Rcomplex beta = asComplex(s_beta);
	if (order < 0 || order > 1)
		order = 0;
	if (!R_FINITE(beta.r) || !R_FINITE(beta.i))
		error(_("'%s' is not a number or not finite"), "beta");
	SEXP trf = R_NilValue;
	char nm[] = "..........Cholesky.";
	nm[18] = (order > 0) ? '+' : '-';
	if (super == NA_LOGICAL || super == 0) {
		memcpy(nm, "simplicial", 10);
		trf = get_factor(s_obj, nm);
		if (!isNull(trf)) super = 0;
	}
	if (isNull(trf) && (super == NA_LOGICAL || super != 0)) {
		memcpy(nm, "supernodal", 10);
		trf = get_factor(s_obj, nm);
		if (!isNull(trf)) super = 1;
	}
	if (beta.r != 0.0 || beta.i != 0.0 || isNull(trf) ||
	    (super == 0 && ll != 0)) {
		if (beta.r != 0.0 || beta.i != 0.0) {
			PROTECT(trf);
			trf = pCMatrix_trf_(s_obj, trf, warn, order, &ll, &super, beta);
			UNPROTECT(1);
		} else {
			if (isNull(trf)) {
			int zz_ = 0;
			trf = pCMatrix_trf_(s_obj, trf, warn, order, &zz_, &super, beta);
			if (!isNull(trf)) {
			memcpy(nm, (super == 0) ? "simplicial" : "supernodal", 10);
			PROTECT(trf);
			set_factor(s_obj, nm, trf);
			UNPROTECT(1);
			}
			}
			if (!isNull(trf) && super == 0 && ll != 0) {
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

SEXP sparseCholesky_update(SEXP s_trf, SEXP s_obj, SEXP s_beta)
{
	/* defined in ./objects.c : */
	char Matrix_shape(SEXP);

	Rcomplex beta = asComplex(s_beta);
	if (!R_FINITE(beta.r) || !R_FINITE(beta.i))
		error(_("'%s' is not a number or not finite"), "beta");

	cholmod_sparse *A = M2CHS(s_obj, 1);
	cholmod_factor *L = M2CHF(s_trf, 1);
	double betaRI[2]; betaRI[0] = beta.r; betaRI[1] = beta.i;

	if (Matrix_shape(s_obj) == 's') {
		SEXP uplo = GET_SLOT(s_obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		A->stype = (ul == 'U') ? 1 : -1;
	}

	L = cholmod_copy_factor(L, &c);

	c.final_asis = 0;
	c.final_ll = L->is_ll;
	c.final_super = L->is_super;
	c.final_pack = 1;
	c.final_monotonic = 1;

	cholmod_factorize_p(A, betaRI, NULL, 0, L, &c);
	cholmod_defaults(&c);

#define UPDOWN_FINISH \
	do { \
		SEXP dimnames = PROTECT(GET_SLOT(s_trf, Matrix_DimNamesSym)); \
		PROTECT(s_trf = CHF2M(L, 1)); \
		cholmod_free_factor(&L, &c); \
		if (TYPEOF(s_trf) == CHARSXP) \
			error("%s", CHAR(s_trf)); \
		SET_SLOT(s_trf, Matrix_DimNamesSym, dimnames); \
		UNPROTECT(2); \
	} while (0)

	UPDOWN_FINISH;
	return s_trf;
}

SEXP sparseCholesky_updown(SEXP s_trf, SEXP s_obj, SEXP s_update)
{
	/* defined in ./objects.c : */
	char Matrix_shape(SEXP);

	cholmod_sparse *A = M2CHS(s_obj, 1);
	cholmod_factor *L = M2CHF(s_trf, 1);

	if (Matrix_shape(s_obj) == 's') {
		SEXP uplo = GET_SLOT(s_obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		A->stype = (ul == 'U') ? 1 : -1;
	}

	L = cholmod_copy_factor(L, &c);
	cholmod_updown(asLogical(s_update) != 0, A, L, &c);

	UPDOWN_FINISH;
	return s_trf;
}

SEXP sparseCholesky_diag_get(SEXP s_trf, SEXP s_square)
{
	cholmod_factor *L = M2CHF(s_trf, 1);
	int n = (int) L->n, square = asLogical(s_square);
	SEXP y = allocVector(REALSXP, n);
	double *py = REAL(y);
	if (L->is_super) {
		int k, j, nc,
			nsuper = (int) L->nsuper,
			*psuper = (int *) L->super,
			*ppi = (int *) L->pi,
			*ppx = (int *) L->px;
		double *px = (double *) L->x, *px_;
		R_xlen_t nr1a;
		for (k = 0; k < nsuper; ++k) {
			nc = psuper[k+1] - psuper[k];
			nr1a = (R_xlen_t) (ppi[k+1] - ppi[k]) + 1;
			px_ = px + ppx[k];
			for (j = 0; j < nc; ++j) {
				*py = *px_;
				if (square)
					*py *= *py;
				++py;
				px_ += nr1a;
			}
		}
	} else {
		square = square && L->is_ll;
		int j, *pp = (int *) L->p;
		double *px = (double *) L->x;
		for (j = 0; j < n; ++j) {
			*py = px[pp[j]];
			if (square)
				*py *= *py;
			++py;
		}
	}
	return y;
}

SEXP denseBunchKaufman_expand(SEXP s_trf)
{
	SEXP dim = PROTECT(GET_SLOT(s_trf, Matrix_DimSym));
	int n = INTEGER(dim)[1];

	SEXP x = PROTECT(GET_SLOT(s_trf, Matrix_xSym));
	int packed = XLENGTH(x) != (int_fast64_t) n * n;

	SEXP P_ = PROTECT(newObject("pMatrix"));
	char cl[] = "..CMatrix";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	cl[1] = 't';
	SEXP T_ = PROTECT(newObject(cl));
	cl[1] = 's';
	SEXP D_ = PROTECT(newObject(cl));

	if (n > 0) {
		SET_SLOT(P_, Matrix_DimSym, dim);
		SET_SLOT(T_, Matrix_DimSym, dim);
		SET_SLOT(D_, Matrix_DimSym, dim);
	}

	SEXP uplo = PROTECT(GET_SLOT(s_trf, Matrix_uploSym));
	char ul = CHAR(STRING_ELT(uplo, 0))[0];
	if (ul != 'U') {
		SET_SLOT(T_, Matrix_uploSym, uplo);
		SET_SLOT(D_, Matrix_uploSym, uplo);
	}
	UNPROTECT(1); /* uplo */

	if (TYPEOF(x) == CPLXSXP) {
		SEXP trans = PROTECT(GET_SLOT(s_trf, Matrix_transSym));
		char ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(D_, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	SEXP diag = PROTECT(mkString("U"));
	SET_SLOT(T_, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	int i, j, s;
	R_xlen_t n1a = (R_xlen_t) n + 1;

	SEXP pivot = PROTECT(GET_SLOT(s_trf, Matrix_permSym)),
		D_p = PROTECT(allocVector(INTSXP, n1a));
	int *ppivot = INTEGER(pivot), *D_pp = INTEGER(D_p),
		b = n, dp = (ul == 'U') ? 1 : 2;
	D_pp[0] = 0;
	j = 0;
	while (j < n) {
		if (ppivot[j] > 0) {
			D_pp[j+1] = D_pp[j] + 1;
			j += 1;
		} else {
			D_pp[j+1] = D_pp[j] + dp;
			D_pp[j+2] = D_pp[j] + 3;
			j += 2;
			b -= 1;
		}
	}
	SET_SLOT(D_, Matrix_pSym, D_p);
	UNPROTECT(1); /* D_p */

	SEXP P, P_perm, T, T_p, T_i, T_x,
		D_i = PROTECT(allocVector(INTSXP, D_pp[n])),
		D_x = PROTECT(allocVector(TYPEOF(x), D_pp[n]));
	int *P_pperm, *T_pp, *T_pi, *D_pi = INTEGER(D_i);

	R_xlen_t len = (R_xlen_t) 2 * b + 1, k = (ul == 'U') ? len - 2 : 0;
	SEXP ans = PROTECT(allocVector(VECSXP, len));

#define EXPAND(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *T_px, *D_px = _PTR_(D_x), *px = _PTR_(x); \
		 \
		j = 0; \
		while (b--) { \
			s = (ppivot[j] > 0) ? 1 : 2; \
			dp = (ul == 'U') ? j : n - j - s; \
			 \
			PROTECT(P = duplicate(P_)); \
			PROTECT(P_perm = allocVector(INTSXP, n)); \
			PROTECT(T = duplicate(T_)); \
			PROTECT(T_p = allocVector(INTSXP, n1a)); \
			PROTECT(T_i = allocVector(INTSXP, (R_xlen_t) s * dp)); \
			PROTECT(T_x = allocVector(TYPEOF(x), (R_xlen_t) s * dp)); \
			 \
			P_pperm = INTEGER(P_perm); \
			T_pp = INTEGER(T_p); \
			T_pi = INTEGER(T_i); \
			T_px = _PTR_(T_x); \
			T_pp[0] = 0; \
			 \
			for (i = 0; i < j; ++i) { \
				T_pp[i+1] = 0; \
				P_pperm[i] = i + 1; \
			} \
			for (i = j; i < j+s; ++i) { \
				T_pp[i+1] = T_pp[i] + dp; \
				P_pperm[i] = i + 1; \
			} \
			for (i = j+s; i < n; ++i) { \
				T_pp[i+1] = T_pp[i]; \
				P_pperm[i] = i + 1; \
			} \
			 \
			if (s == 1) { \
				P_pperm[j] = ppivot[j]; \
				P_pperm[ppivot[j]-1] = j + 1; \
			} else if (ul == 'U') { \
				P_pperm[j] = -ppivot[j]; \
				P_pperm[-ppivot[j]-1] = j + 1; \
			} else { \
				P_pperm[j+1] = -ppivot[j]; \
				P_pperm[-ppivot[j]-1] = j + 2; \
			} \
			 \
			if (ul == 'U') { \
				for (i = 0; i < j; ++i) { \
					*(T_pi++) = i; \
					*(T_px++) = *(px++); \
				} \
				*(D_pi++) = j; \
				*(D_px++) = *(px++); \
				++j; \
				if (!packed) \
					px += n - j; \
				if (s == 2) { \
					for (i = 0; i < j-1; ++i) { \
						*(T_pi++) = i; \
						*(T_px++) = *(px++); \
					} \
					*(D_pi++) = j - 1; \
					*(D_pi++) = j; \
					*(D_px++) = *(px++); \
					*(D_px++) = *(px++); \
					++j; \
					if (!packed) \
						px += n - j; \
				} \
			} else { \
				if (s == 2) { \
					*(D_pi++) = j; \
					*(D_pi++) = j + 1; \
					*(D_px++) = *(px++); \
					*(D_px++) = *(px++); \
					for (i = j+2; i < n; ++i) { \
						*(T_pi++) = i; \
						*(T_px++) = *(px++); \
					} \
					++j; \
					if (!packed) \
						px += j; \
				} \
				*(D_pi++) = j; \
				*(D_px++) = *(px++); \
				for (i = j+1; i < n; ++i) { \
					*(T_pi++) = i; \
					*(T_px++) = *(px++); \
				} \
				++j; \
				if (!packed) \
					px += j; \
			} \
			 \
			SET_SLOT(P, Matrix_permSym, P_perm); \
			SET_SLOT(T, Matrix_pSym, T_p); \
			SET_SLOT(T, Matrix_iSym, T_i); \
			SET_SLOT(T, Matrix_xSym, T_x); \
			 \
			if (ul == 'U') { \
				SET_VECTOR_ELT(ans, k-1, P); \
				SET_VECTOR_ELT(ans, k  , T); \
				k -= 2; \
			} else { \
				SET_VECTOR_ELT(ans, k  , P); \
				SET_VECTOR_ELT(ans, k+1, T); \
				k += 2; \
			} \
			UNPROTECT(6); /* T_x, T_i, T_p, T, P_perm, P */ \
		} \
	} while (0)

	if (TYPEOF(x) == CPLXSXP)
		EXPAND(Rcomplex, COMPLEX);
	else
		EXPAND(double, REAL);

	SET_SLOT(D_, Matrix_iSym, D_i);
	SET_SLOT(D_, Matrix_xSym, D_x);
	SET_VECTOR_ELT(ans, len-1, D_);

	UNPROTECT(9); /* ans, D_x, D_i, pivot, D_, T_, P_, x, dim */
	return ans;
}
