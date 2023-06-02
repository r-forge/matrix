#include "dspMatrix.h"

double get_norm_dsp(SEXP obj, const char *typstr)
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

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0

SEXP dspMatrix_getDiag(SEXP x)

{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(REALSXP, n));

    d_packed_getDiag(REAL(val), x, n);
    UNPROTECT(1);
    return val;
}

SEXP lspMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(LGLSXP, n));

    l_packed_getDiag(LOGICAL(val), x, n);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return d_packed_setDiag(REAL(d), LENGTH(d), x, n);
}

SEXP lspMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return l_packed_setDiag(INTEGER(d), LENGTH(d), x, n);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0

SEXP dspMatrix_as_dsyMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dsyMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_unpack(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n*n)),
		  REAL(GET_SLOT(from, Matrix_xSym)),
		  n,
		  *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
		  NUN);
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
