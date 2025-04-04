/* C implementation of methods for rcond */

#include "Lapack-etc.h"
#include "Mdefines.h"

static
char La_rcond_type(SEXP s)
{
#define ARGNAME "norm"
	if (TYPEOF(s) != STRSXP)
		Rf_error(_("argument '%s' is not of type \"%s\""),
		         ARGNAME, "character");
	if (LENGTH(s) == 0)
		Rf_error(_("argument '%s' has length %d"),
		         ARGNAME, 0);
	s = STRING_ELT(s, 0);
	if (CHAR(s)[0] == '\0' || CHAR(s)[1] != '\0')
		Rf_error(_("argument '%s' (\"%s\") does not have string length %d"),
		         ARGNAME, CHAR(s), 1);
	char type = '\0';
	switch (CHAR(s)[0]) {
	case 'O':
	case 'o':
	case '1':
		type = 'O';
		break;
	case 'I':
	case 'i':
		type = 'I';
		break;
	default:
		Rf_error(_("argument '%s' (\"%s\") is not \"%s\", \"%s\", or \"%s\""),
		         ARGNAME, CHAR(s), "O", "1", "I");
		break;
	}
	return type;
#undef ARGNAME
}

SEXP geMatrix_rcond(SEXP s_obj, SEXP trf, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int *pdim = DIM(s_obj), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("non-square matrix arguments not supported by '%s'"),
		         __func__);
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	double *rwork = (double *) R_alloc((size_t) n * 2, sizeof(double));
	norm =
	F77_CALL(zlange)(&type, &n, &n, COMPLEX(x), &n, work FCONE);
	F77_CALL(zgecon)(&type, &n,     COMPLEX(y), &n, &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	norm =
	F77_CALL(dlange)(&type, &n, &n,    REAL(x), &n, work FCONE);
	F77_CALL(dgecon)(&type, &n,        REAL(y), &n, &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(2); /* x, y */

	return Rf_ScalarReal(rcond);
}

SEXP syMatrix_rcond(SEXP s_obj, SEXP trf, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int n = DIM(s_obj)[0];
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));
	char ul = UPLO(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	if (TRANS(s_obj) == 'C') {
	norm =
	F77_CALL(zlanhe)(&type, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	F77_CALL(zhecon)(       &ul, &n, COMPLEX(y), &n, INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	} else {
	norm =
	F77_CALL(zlansy)(&type, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	F77_CALL(zsycon)(       &ul, &n, COMPLEX(y), &n, INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	}
	} else {
	double * work = (double *) R_alloc((size_t) n * 2, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	norm =
	F77_CALL(dlansy)(&type, &ul, &n,    REAL(x), &n, work FCONE FCONE);
	F77_CALL(dsycon)(       &ul, &n,    REAL(y), &n, INTEGER(pivot), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(3); /* x, y, pivot */

	return Rf_ScalarReal(rcond);
}

SEXP spMatrix_rcond(SEXP s_obj, SEXP trf, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));
	char ul = UPLO(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym)),
		pivot = PROTECT(GET_SLOT(trf, Matrix_permSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	if (TRANS(s_obj) == 'C') {
	norm =
	F77_CALL(zlanhp)(&type, &ul, &n, COMPLEX(x), work FCONE FCONE);
	F77_CALL(zhpcon)(       &ul, &n, COMPLEX(y), INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	} else {
	norm =
	F77_CALL(zlansp)(&type, &ul, &n, COMPLEX(x), work FCONE FCONE);
	F77_CALL(zspcon)(       &ul, &n, COMPLEX(y), INTEGER(pivot), &norm, &rcond,
	                 (Rcomplex *) work,        &info FCONE);
	}
	} else {
	double * work = (double *) R_alloc((size_t) n * 2, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	norm =
	F77_CALL(dlansp)(&type, &ul, &n,    REAL(x), work FCONE FCONE);
	F77_CALL(dspcon)(       &ul, &n,    REAL(y), INTEGER(pivot), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(3); /* x, y, pivot */

	return Rf_ScalarReal(rcond);
}

SEXP poMatrix_rcond(SEXP s_obj, SEXP trf, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));
	char ul = UPLO(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	double *rwork = (double *) R_alloc((size_t) n    , sizeof(double));
	norm =
	F77_CALL(zlansy)(&type, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	F77_CALL(zpocon)(       &ul, &n, COMPLEX(y), &n, &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) n * 3, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	norm =
	F77_CALL(dlansy)(&type, &ul, &n,    REAL(x), &n, work FCONE FCONE);
	F77_CALL(dpocon)(       &ul, &n, REAL(y), &n, &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(2); /* x, y */

	return Rf_ScalarReal(rcond);
}

SEXP ppMatrix_rcond(SEXP s_obj, SEXP trf, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));
	char ul = UPLO(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym)),
		y = PROTECT(GET_SLOT(trf, Matrix_xSym));
	double norm, rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	double *rwork = (double *) R_alloc((size_t) n    , sizeof(double));
	norm =
	F77_CALL(zlansp)(&type, &ul, &n, COMPLEX(x), work FCONE FCONE);
	F77_CALL(zppcon)(       &ul, &n, COMPLEX(y), &norm, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	norm =
	F77_CALL(dlansp)(&type, &ul, &n,    REAL(x), work FCONE FCONE);
	F77_CALL(dppcon)(       &ul, &n,    REAL(y), &norm, &rcond,
	                 (double   *) work, iwork, &info FCONE);
	}
	UNPROTECT(2); /* x, y */

	return Rf_ScalarReal(rcond);
}

SEXP trMatrix_rcond(SEXP s_obj, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));
	char ul = UPLO(s_obj), nu = DIAG(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	double *rwork = (double *) R_alloc((size_t) n    , sizeof(double));
	F77_CALL(ztrcon)(&type, &ul, &nu, &n, COMPLEX(x), &n, &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE FCONE FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) n * 3, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	F77_CALL(dtrcon)(&type, &ul, &nu, &n,    REAL(x), &n, &rcond,
	                 (double   *) work, iwork, &info FCONE FCONE FCONE);
	}
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(rcond);
}

SEXP tpMatrix_rcond(SEXP s_obj, SEXP s_type)
{
	char type = La_rcond_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return(Rf_ScalarReal(R_PosInf));
	char ul = UPLO(s_obj), nu = DIAG(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double rcond;
	int info;
	if (TYPEOF(x) == CPLXSXP) {
	double * work = (double *) R_alloc((size_t) n * 4, sizeof(double));
	double *rwork = (double *) R_alloc((size_t) n    , sizeof(double));
	F77_CALL(ztpcon)(&type, &ul, &nu, &n, COMPLEX(x), &rcond,
	                 (Rcomplex *) work, rwork, &info FCONE FCONE FCONE);
	} else {
	double * work = (double *) R_alloc((size_t) n * 3, sizeof(double));
	int    *iwork = (int    *) R_alloc((size_t) n    , sizeof(int   ));
	F77_CALL(dtpcon)(&type, &ul, &nu, &n,    REAL(x), &rcond,
	                 (double   *) work, iwork, &info FCONE FCONE FCONE);
	}
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(rcond);
}
