/* double (precision) Triangular Packed Matrices
 * Note: this means *square* {n x n} matrices
*/

/* MJ: no longer needed ... nothing below */
#if 0
#include "dtpMatrix.h"
#endif /* MJ */

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0

// also applicable to dspMatrix , dppMatrix :
SEXP dtpMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(REALSXP, n));

    tr_d_packed_getDiag(REAL(val), x, n);
    UNPROTECT(1);
    return val;
}

// also applicable to lspMatrix :
SEXP ltpMatrix_getDiag(SEXP x)
{
    int n = *INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP val = PROTECT(allocVector(LGLSXP, n));

    tr_l_packed_getDiag(LOGICAL(val), x, n);
    UNPROTECT(1);
    return val;
}

SEXP dtpMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return tr_d_packed_setDiag(REAL(d), LENGTH(d), x, n);
}

SEXP ltpMatrix_setDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return tr_l_packed_setDiag(INTEGER(d), LENGTH(d), x, n);
}

/* was unused, not replaced: */
SEXP dtpMatrix_addDiag(SEXP x, SEXP d)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    return tr_d_packed_addDiag(REAL(d), LENGTH(d), x, n);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0

SEXP dtpMatrix_as_dtrMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_unpack(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n*n)),
		  REAL(GET_SLOT(from, Matrix_xSym)),
		  n,
		  *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
		  *CHAR(STRING_ELT(diag, 0)) == 'N' ? NUN : UNT);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

