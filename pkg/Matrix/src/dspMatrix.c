/* MJ: no longer needed ... nothing below */
#if 0
#include "dspMatrix.h"
#endif /* MJ */

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
