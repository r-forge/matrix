#include "Lapack-etc.h"
#include "Mutils.h"
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
	char type_ = '\0';
	switch (type[0]) {
	case 'M':
	case 'm':
		type_ = 'M';
		break;
	case 'O':
	case 'o':
	case '1':
		type_ = 'O';
		break;
	case 'I':
	case 'i':
		type_ = 'I';
		break;
	case 'F':
	case 'f':
	case 'E':
	case 'e':
		type_ = 'F';
		break;
	default:
		error(_("argument '%s' (\"%s\") is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\""),
		      ARGNAME, type, "M", "O", "1", "I", "F", "E");
		break;
	}
	return type_;
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
	char type_ = '\0';
	switch (type[0]) {
	case 'O':
	case 'o':
	case '1':
		type_ = 'O';
		break;
	case 'I':
	case 'i':
		type_ = 'I';
		break;
	default:
		error(_("argument '%s' (\"%s\") is not \"%s\", \"%s\", or \"%s\""),
		      ARGNAME, type, "O", "1", "I");
		break;
	}
	return type_;
#undef ARGNAME
}

SEXP dgeMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m == 0 || n == 0)
		return ScalarReal(0.0);

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'I')
		work = (double *) R_alloc((size_t) m, sizeof(double));
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlange)(&type_, &m, &n, COMPLEX(x), &m,
	                 work FCONE);
	else
#endif
	norm =
	F77_CALL(dlange)(&type_, &m, &n,    REAL(x), &m,
	                 work FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dsyMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[1];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'O' || type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlansy)(&type_, &uplo_, &n, COMPLEX(x), &n,
	                 work FCONE FCONE);
	else
#endif
	norm =
	F77_CALL(dlansy)(&type_, &uplo_, &n,    REAL(x), &n,
	                 work FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dspMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[1];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'O' || type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlansp)(&type_, &uplo_, &n, COMPLEX(x),
	                 work FCONE FCONE);
	else
#endif
	norm =
	F77_CALL(dlansp)(&type_, &uplo_, &n,    REAL(x),
	                 work FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dtrMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char diag_ = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlantr)(&type_, &uplo_, &diag_, &n, &n, COMPLEX(x), &n,
	                 work FCONE FCONE FCONE);
	else
#endif
	norm =
	F77_CALL(dlantr)(&type_, &uplo_, &diag_, &n, &n,    REAL(x), &n,
	                 work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dtpMatrix_norm(SEXP obj, SEXP type)
{
	char type_ = La_norm_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return ScalarReal(0.0);

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char diag_ = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type_ == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlantp)(&type_, &uplo_, &diag_, &n, COMPLEX(x),
	                 work FCONE FCONE FCONE);
	else
#endif
	norm =
	F77_CALL(dlantp)(&type_, &uplo_, &diag_, &n,    REAL(x),
	                 work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return ScalarReal(norm);
}

SEXP dgeMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
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
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double *rwork = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	norm =
	F77_CALL(zlange)(&type_, &n, &n, COMPLEX(x), &n,
	                 work FCONE);
	F77_CALL(zgecon)(&type_, &n,     COMPLEX(y), &n, &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
#endif
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlange)(&type_, &n, &n,    REAL(x), &n,
	                 work FCONE);
	F77_CALL(dgecon)(&type_, &n,        REAL(y), &n, &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP dsyMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond;
	int info;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	norm =
	F77_CALL(zlansy)(&type_, &uplo_, &n, COMPLEX(x), &n,
	                 work FCONE FCONE);
	F77_CALL(zsycon)(        &uplo_, &n, COMPLEX(y), &n,
	                 INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	} else {
#endif
	double * work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansy)(&type_, &uplo_, &n,    REAL(x), &n,
	                 work FCONE FCONE);
	F77_CALL(dsycon)(        &uplo_, &n,    REAL(y), &n,
	                 INTEGER(pivot), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(3); /* x, y, pivot */

	return ScalarReal(rcond);
}

SEXP dspMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond;
	int info;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	norm =
	F77_CALL(zlansp)(&type_, &uplo_, &n, COMPLEX(x),
	                 work FCONE FCONE);
	F77_CALL(zspcon)(        &uplo_, &n, COMPLEX(y),
	                 INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	} else {
#endif
	double * work = (double *) R_alloc((size_t) 2 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansp)(&type_, &uplo_, &n,    REAL(x),
	                 work FCONE FCONE);
	F77_CALL(dspcon)(        &uplo_, &n,    REAL(y),
	                 INTEGER(pivot), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(3); /* x, y, pivot */

	return ScalarReal(rcond);
}

SEXP dpoMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	norm =
	F77_CALL(zlansy)(&type_, &uplo_, &n, COMPLEX(x), &n,
	                 work FCONE FCONE);
	F77_CALL(zpocon)(        &uplo_, &n, COMPLEX(y), &n, &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
#endif
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansy)(&type_, &uplo_, &n, REAL(x), &n,
	                 work FCONE FCONE);
	F77_CALL(dpocon)(        &uplo_, &n, REAL(y), &n, &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP dppMatrix_rcond(SEXP obj, SEXP trf, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	norm =
	F77_CALL(zlansp)(&type_, &uplo_, &n, COMPLEX(x),
	                 work FCONE FCONE);
	F77_CALL(zppcon)(        &uplo_, &n, COMPLEX(y), &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
#endif
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	norm =
	F77_CALL(dlansp)(&type_, &uplo_, &n,    REAL(x),
	                 work FCONE FCONE);
	F77_CALL(dppcon)(        &uplo_, &n,    REAL(y), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(2); /* x, y */

	return ScalarReal(rcond);
}

SEXP dtrMatrix_rcond(SEXP obj, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char diag_ = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double rcond;
	int info;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	F77_CALL(ztrcon)(&type_, &uplo_, &diag_, &n, COMPLEX(x), &n, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE FCONE FCONE);
	} else {
#endif
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	F77_CALL(dtrcon)(&type_, &uplo_, &diag_, &n,    REAL(x), &n, &rcond,
	                 (double   *) work, iwork, &info FCONE FCONE FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(1); /* x */

	return ScalarReal(rcond);
}

SEXP dtpMatrix_rcond(SEXP obj, SEXP type)
{
	char type_ = La_rcond_type(type);
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	if (n == 0)
		return(ScalarReal(R_PosInf));

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	char uplo_ = CHAR(STRING_ELT(uplo, 0))[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char diag_ = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double rcond;
	int info;
#ifdef MATRIX_ENABLE_ZMATRIX
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) 4 * n, sizeof(double));
	double *rwork = (double *) R_alloc((size_t)     n, sizeof(double));
	F77_CALL(ztpcon)(&type_, &uplo_, &diag_, &n, COMPLEX(x), &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE FCONE FCONE);
	} else {
#endif
	double * work = (double *) R_alloc((size_t) 3 * n, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t)     n, sizeof(int   ));
	F77_CALL(dtpcon)(&type_, &uplo_, &diag_, &n,    REAL(x), &rcond,
	                 (double   *) work, iwork, &info FCONE FCONE FCONE);
#ifdef MATRIX_ENABLE_ZMATRIX
	}
#endif
	UNPROTECT(1); /* x */

	return ScalarReal(rcond);
}
