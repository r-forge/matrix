#include "Lapack-etc.h"
#include "Mdefines.h"
#include "kappa.h"

static
char La_norm_type(SEXP s)
{
#define ARGNAME "type"
	if (TYPEOF(s) != STRSXP)
		error(_("argument '%s' is not of type \"%s\""),
		      ARGNAME, "character");
	if (LENGTH(s) == 0)
		error(_("argument '%s' has length %d"),
		      ARGNAME, 0);
	const char *type = CHAR(STRING_ELT(s, 0));
	if (type[0] == '\0' || type[1] != '\0')
		error(_("argument '%s' (\"%s\") does not have string length %d"),
		      ARGNAME, type, 1);
	char t = '\0';
	switch (type[0]) {
	case 'M':
	case 'm':
		t = 'M';
		break;
	case 'O':
	case 'o':
	case '1':
		t = 'O';
		break;
	case 'I':
	case 'i':
		t = 'I';
		break;
	case 'F':
	case 'f':
	case 'E':
	case 'e':
		t = 'F';
		break;
	default:
		error(_("argument '%s' (\"%s\") is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\""),
		      ARGNAME, type, "M", "O", "1", "I", "F", "E");
		break;
	}
	return t;
#undef ARGNAME
}

static
char La_rcond_type(SEXP s)
{
#define ARGNAME "norm"
	if (TYPEOF(s) != STRSXP)
		error(_("argument '%s' is not of type \"%s\""),
		      ARGNAME, "character");
	if (LENGTH(s) == 0)
		error(_("argument '%s' has length %d"),
		      ARGNAME, 0);
	const char *type = CHAR(STRING_ELT(s, 0));
	if (type[0] == '\0' || type[1] != '\0')
		error(_("argument '%s' (\"%s\") does not have string length %d"),
		      ARGNAME, type, 1);
	char t = '\0';
	switch (type[0]) {
	case 'O':
	case 'o':
	case '1':
		t = 'O';
		break;
	case 'I':
	case 'i':
		t = 'I';
		break;
	default:
		error(_("argument '%s' (\"%s\") is not \"%s\", \"%s\", or \"%s\""),
		      ARGNAME, type, "O", "1", "I");
		break;
	}
	return t;
#undef ARGNAME
}

SEXP geMatrix_norm(SEXP obj, SEXP type)
{
	char t = La_norm_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m == 0 || n == 0)
		return ScalarReal(0.0);

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (t == 'I')
		work = (double *) R_alloc((size_t) m, sizeof(double));
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlange)(&t, &m, &n, COMPLEX(x), &m, work FCONE);
	else
	norm =
	F77_CALL(dlange)(&t, &m, &n,    REAL(x), &m, work FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP syMatrix_norm(SEXP obj, SEXP type)
{
	char t = La_norm_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[1];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = *CHAR(STRING_ELT(uplo, 0));

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (t == 'O' || t == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	char ct = *CHAR(STRING_ELT(trans, 0));
	if (ct == 'C')
	norm =
	F77_CALL(zlanhe)(&t, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	else
	norm =
	F77_CALL(zlansy)(&t, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	}
	else
	norm =
	F77_CALL(dlansy)(&t, &ul, &n,    REAL(x), &n, work FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP spMatrix_norm(SEXP obj, SEXP type)
{
	char t = La_norm_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[1];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (t == 'O' || t == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	char ct = *CHAR(STRING_ELT(trans, 0));
	if (ct == 'C')
	norm =
	F77_CALL(zlanhp)(&t, &ul, &n, COMPLEX(x), work FCONE FCONE);
	else
	norm =
	F77_CALL(zlansp)(&t, &ul, &n, COMPLEX(x), work FCONE FCONE);
	}
	else
	norm =
	F77_CALL(dlansp)(&t, &ul, &n,    REAL(x), work FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP trMatrix_norm(SEXP obj, SEXP type)
{
	char t = La_norm_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (t == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlantr)(&t, &ul, &di, &n, &n, COMPLEX(x), &n, work FCONE FCONE FCONE);
	else
	norm =
	F77_CALL(dlantr)(&t, &ul, &di, &n, &n,    REAL(x), &n, work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP tpMatrix_norm(SEXP obj, SEXP type)
{
	char t = La_norm_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (t == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlantp)(&t, &ul, &di, &n, COMPLEX(x), work FCONE FCONE FCONE);
	else
	norm =
	F77_CALL(dlantp)(&t, &ul, &di, &n,    REAL(x), work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP geMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n)
		error(_("%s(%s) is undefined: '%s' is not square"), "rcond", "x", "x");
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	double *rwork = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	norm =
	F77_CALL(zlange)(&t, &n, &n, COMPLEX(x), &n, work FCONE);
	F77_CALL(zgecon)(&t, &n,     COMPLEX(y), &n, &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlange)(&t, &n, &n,    REAL(x), &n, work FCONE);
	F77_CALL(dgecon)(&t, &n,        REAL(y), &n, &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP syMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	char ct = *CHAR(STRING_ELT(trans, 0));
	if (ct == 'C') {
	norm =
	F77_CALL(zlanhe)(&t, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	F77_CALL(zhecon)(    &ul, &n, COMPLEX(y), &n, INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	} else {
	norm =
	F77_CALL(zlansy)(&t, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	F77_CALL(zsycon)(    &ul, &n, COMPLEX(y), &n, INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	}
	} else {
	double * work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansy)(&t, &ul, &n,    REAL(x), &n, work FCONE FCONE);
	F77_CALL(dsycon)(    &ul, &n,    REAL(y), &n, INTEGER(pivot), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(3); /* x, y, pivot */

	return ScalarReal(rcond);
}

SEXP spMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	char ct = *CHAR(STRING_ELT(trans, 0));
	if (ct == 'C') {
	norm =
	F77_CALL(zlanhp)(&t, &ul, &n, COMPLEX(x), work FCONE FCONE);
	F77_CALL(zhpcon)(    &ul, &n, COMPLEX(y), INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	} else {
	norm =
	F77_CALL(zlansp)(&t, &ul, &n, COMPLEX(x), work FCONE FCONE);
	F77_CALL(zspcon)(    &ul, &n, COMPLEX(y), INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	}
	} else {
	double * work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansp)(&t, &ul, &n,    REAL(x), work FCONE FCONE);
	F77_CALL(dspcon)(    &ul, &n,    REAL(y), INTEGER(pivot), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(3); /* x, y, pivot */

	return ScalarReal(rcond);
}

SEXP poMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	norm =
	F77_CALL(zlansy)(&t, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	F77_CALL(zpocon)(    &ul, &n, COMPLEX(y), &n, &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansy)(&t, &ul, &n,    REAL(x), &n, work FCONE FCONE);
	F77_CALL(dpocon)(    &ul, &n, REAL(y), &n, &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP ppMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	norm =
	F77_CALL(zlansp)(&t, &ul, &n, COMPLEX(x), work FCONE FCONE);
	F77_CALL(zppcon)(    &ul, &n, COMPLEX(y), &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansp)(&t, &ul, &n,    REAL(x), work FCONE FCONE);
	F77_CALL(dppcon)(    &ul, &n,    REAL(y), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP trMatrix_rcond(SEXP obj, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	F77_CALL(ztrcon)(&t, &ul, &di, &n, COMPLEX(x), &n, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE FCONE FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	F77_CALL(dtrcon)(&t, &ul, &di, &n,    REAL(x), &n, &rcond,
	                 (double   *) work, iwork, &info FCONE FCONE FCONE);
	}
	UNPROTECT(1); /* x */

	return ScalarReal(rcond);
}

SEXP tpMatrix_rcond(SEXP obj, SEXP type)
{
	char t = La_rcond_type(type);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char di = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	F77_CALL(ztpcon)(&t, &ul, &di, &n, COMPLEX(x), &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE FCONE FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	F77_CALL(dtpcon)(&t, &ul, &di, &n,    REAL(x), &rcond,
	                 (double   *) work, iwork, &info FCONE FCONE FCONE);
	}
	UNPROTECT(1); /* x */

	return ScalarReal(rcond);
}
