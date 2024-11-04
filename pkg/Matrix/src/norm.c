/* C implementation of methods for norm */

#include "Lapack-etc.h"
#include "Mdefines.h"

static
char La_norm_type(SEXP s)
{
#define ARGNAME "type"
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
	case 'M':
	case 'm':
		type = 'M';
		break;
	case 'O':
	case 'o':
	case '1':
		type = 'O';
		break;
	case 'I':
	case 'i':
		type = 'I';
		break;
	case 'F':
	case 'f':
	case 'E':
	case 'e':
		type = 'F';
		break;
	default:
		Rf_error(_("argument '%s' (\"%s\") is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\""),
		         ARGNAME, CHAR(s), "M", "O", "1", "I", "F", "E");
		break;
	}
	return type;
#undef ARGNAME
}

SEXP geMatrix_norm(SEXP s_obj, SEXP s_type)
{
	char type = La_norm_type(s_type);

	int *pdim = DIM(s_obj), m = pdim[0], n = pdim[1];
	if (m == 0 || n == 0)
		return Rf_ScalarReal(0.0);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type == 'I')
		work = (double *) R_alloc((size_t) m, sizeof(double));
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlange)(&type, &m, &n, COMPLEX(x), &m, work FCONE);
	else
	norm =
	F77_CALL(dlange)(&type, &m, &n,    REAL(x), &m, work FCONE);
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(norm);
}

SEXP syMatrix_norm(SEXP s_obj, SEXP s_type)
{
	char type = La_norm_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return Rf_ScalarReal(0.0);
	char ul = UPLO(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type == 'O' || type == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	if (TRANS(s_obj) == 'C')
	norm =
	F77_CALL(zlanhe)(&type, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	else
	norm =
	F77_CALL(zlansy)(&type, &ul, &n, COMPLEX(x), &n, work FCONE FCONE);
	}
	else
	norm =
	F77_CALL(dlansy)(&type, &ul, &n,    REAL(x), &n, work FCONE FCONE);
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(norm);
}

SEXP spMatrix_norm(SEXP s_obj, SEXP s_type)
{
	char type = La_norm_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return Rf_ScalarReal(0.0);
	char ul = UPLO(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type == 'O' || type == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP) {
	if (TRANS(s_obj) == 'C')
	norm =
	F77_CALL(zlanhp)(&type, &ul, &n, COMPLEX(x), work FCONE FCONE);
	else
	norm =
	F77_CALL(zlansp)(&type, &ul, &n, COMPLEX(x), work FCONE FCONE);
	}
	else
	norm =
	F77_CALL(dlansp)(&type, &ul, &n,    REAL(x), work FCONE FCONE);
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(norm);
}

SEXP trMatrix_norm(SEXP s_obj, SEXP s_type)
{
	char type = La_norm_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return Rf_ScalarReal(0.0);
	char ul = UPLO(s_obj), nu = DIAG(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlantr)(&type, &ul, &nu, &n, &n, COMPLEX(x), &n, work FCONE FCONE FCONE);
	else
	norm =
	F77_CALL(dlantr)(&type, &ul, &nu, &n, &n,    REAL(x), &n, work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(norm);
}

SEXP tpMatrix_norm(SEXP s_obj, SEXP s_type)
{
	char type = La_norm_type(s_type);

	int n = DIM(s_obj)[1];
	if (n == 0)
		return Rf_ScalarReal(0.0);
	char ul = UPLO(s_obj), nu = DIAG(s_obj);

	SEXP x = PROTECT(GET_SLOT(s_obj, Matrix_xSym));
	double norm, *work = NULL;
	if (type == 'I')
		work = (double *) R_alloc((size_t) n, sizeof(double));
	if (TYPEOF(x) == CPLXSXP)
	norm =
	F77_CALL(zlantp)(&type, &ul, &nu, &n, COMPLEX(x), work FCONE FCONE FCONE);
	else
	norm =
	F77_CALL(dlantp)(&type, &ul, &nu, &n,    REAL(x), work FCONE FCONE FCONE);
	UNPROTECT(1); /* x */

	return Rf_ScalarReal(norm);
}
