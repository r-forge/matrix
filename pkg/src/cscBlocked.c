#include "cscBlocked.h"

SEXP cscBlocked_validate(SEXP x)
{
    SEXP pp = GET_SLOT(x, Matrix_pSym),
	ip = GET_SLOT(x, Matrix_iSym),
	xp = GET_SLOT(x, Matrix_xSym),
	dp = getAttrib(xp, R_DimSymbol);
    int *pv = INTEGER(pp),
	*iv = INTEGER(ip),
	*dim = INTEGER(dp),
	ncol = length(pp) - 1;
    int nnz = pv[ncol];

    if (!isReal(xp))
	return ScalarString(mkChar("slot x should be a real array"));
    if (length(dp) != 3)
	return ScalarString(mkChar("slot x should be a 3-dimensional array"));
    if (length(ip) != nnz)
	return ScalarString(mkChar("length of slot i does not matck last element of slot p"));
    if (dim[2] != nnz)
	return ScalarString(mkChar("third dimension of slot x does not match number of nonzeros"));
    return ScalarLogical(1);
}

    
