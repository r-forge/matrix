#include <ctype.h> /* toupper */
#include "kappa.h"

/* La_norm_type() and La_rcond_type() have been in src/include/R_ext/Lapack.h
   and later in src/modules/lapack/Lapack.c but have still not been available
   to package writers ...
*/

static char La_norm_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* aliases */
    else if (typup == 'E')
	typup = 'F';
    else if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error(_("argument type[1]='%s' must be one of 'M','1','O','I','F', or 'E'"),
	      typstr);
    return typup;
}

static char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error(_("argument type[1]='%s' must be one of '1','O', or 'I'"),
	      typstr);
    return typup; /* 'O' or 'I' */
}

static double get_norm_dge(SEXP obj, const char *typstr)
{
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    R_xlen_t i, nx = XLENGTH(x);
    double *px = REAL(x);
    for (i = 0; i < nx; ++i) {
	if (ISNAN(px[i])) {
	    UNPROTECT(1);
	    return NA_REAL;
	}
    }
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim);
    double norm, *work = NULL;

    if (typstr[0] == 'I')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlange)(typstr, pdim, pdim + 1, px, pdim,
			    work FCONE);
    
    UNPROTECT(2);
    return norm;
}

SEXP dgeMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dge(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dgeMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim);
    if (pdim[0] != pdim[1] || pdim[0] < 1)
	error(_("'rcond' requires a square, nonempty matrix"));

    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_rcond_type(CHAR(type));

    SEXP trf = PROTECT(dgeMatrix_trf_(obj, 0)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));
    int info;
    double *px = REAL(x), norm = get_norm_dge(obj, typstr), rcond;
    
    F77_CALL(dgecon)(typstr, pdim, px, pdim, &norm, &rcond,
		     (double *) R_alloc((size_t) 4 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);

    UNPROTECT(4);
    return ScalarReal(rcond);
}

static double get_norm_dtr(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    if (typstr[0] == 'I')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlantr)(typstr, ul, di, pdim, pdim + 1, px, pdim,
			    work FCONE FCONE FCONE);

    UNPROTECT(4);
    return norm;
}

SEXP dtrMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dtr(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dtrMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));

    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_rcond_type(CHAR(type));
    
    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    F77_CALL(dtrcon)(typstr, ul, di, pdim, px, pdim, &rcond,
		     (double *) R_alloc((size_t) 3 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE FCONE FCONE);

    UNPROTECT(5);
    return ScalarReal(rcond);
}

static double get_norm_dtp(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));
    
    if (typstr[0] == 'I')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlantp)(typstr, ul, di, pdim, px,
			    work FCONE FCONE FCONE);

    UNPROTECT(4);
    return norm;
}

SEXP dtpMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dtp(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dtpMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	diag = PROTECT(GET_SLOT(obj, Matrix_diagSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_rcond_type(CHAR(type));
    
    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0)), *di = CHAR(STRING_ELT(diag, 0));

    F77_CALL(dtpcon)(typstr, ul, di, pdim, px, &rcond,
		     (double *) R_alloc((size_t) 3 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE FCONE FCONE);

    UNPROTECT(5);
    return ScalarReal(rcond);
}

static double get_norm_dsy(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    if (typstr[0] == 'I' || typstr[0] == 'O')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlansy)(typstr, ul, pdim, px, pdim, work FCONE FCONE);

    UNPROTECT(3);
    return norm;
}

SEXP dsyMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dsy(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dsyMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dsyMatrix_trf_(obj, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));
    
    int *pdim = INTEGER(dim), *pperm = INTEGER(perm), info;
    double *px = REAL(x), norm = get_norm_dsy(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dsycon)(ul, pdim, px, pdim, pperm, &norm, &rcond,
		     (double *) R_alloc((size_t) 2 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);
    
    UNPROTECT(5);
    return ScalarReal(rcond);
}

static double get_norm_dsp(SEXP obj, const char *typstr)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *pdim = INTEGER(dim);
    double *px = REAL(x), norm, *work = NULL;
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    if (typstr[0] == 'I' || typstr[0] == 'O')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlansp)(typstr, ul, pdim, px, work FCONE FCONE);

    UNPROTECT(3);
    return norm;
}

SEXP dspMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dsp(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dspMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dspMatrix_trf_(obj, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	perm = PROTECT(GET_SLOT(trf, Matrix_permSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));
    
    int *pdim = INTEGER(dim), *pperm = INTEGER(perm), info;
    double *px = REAL(x), norm = get_norm_dsp(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));

    F77_CALL(dspcon)(ul, pdim, px, pperm, &norm, &rcond, 
		     (double *) R_alloc((size_t) 2 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);
    
    UNPROTECT(5);
    return ScalarReal(rcond);
}

SEXP dpoMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dpoMatrix_trf_(obj, 2, 0, -1.0)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));

    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), norm = get_norm_dsy(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));

    F77_CALL(dpocon)(ul, pdim, px, pdim, &norm, &rcond,
		     (double *) R_alloc((size_t) 3 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);

    UNPROTECT(4);
    return ScalarReal(rcond);
}

SEXP dppMatrix_rcond(SEXP obj)
{
    SEXP trf = PROTECT(dppMatrix_trf_(obj, 2)),
	dim = PROTECT(GET_SLOT(trf, Matrix_DimSym)),
	uplo = PROTECT(GET_SLOT(trf, Matrix_uploSym)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));

    int *pdim = INTEGER(dim), info;
    double *px = REAL(x), norm = get_norm_dsp(obj, "O"), rcond;
    const char *ul = CHAR(STRING_ELT(uplo, 0));

    F77_CALL(dppcon)(ul, pdim, px, &norm, &rcond,
		     (double *) R_alloc((size_t) 3 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);

    UNPROTECT(4);
    return ScalarReal(rcond);
}
