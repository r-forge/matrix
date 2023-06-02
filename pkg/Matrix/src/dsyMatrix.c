/* MJ: no longer needed ... nothing below */
#if 0
#include "dsyMatrix.h"
#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0

// this is very close to lsyMatrix_as_lsp*() in ./ldense.c  -- keep synced !
SEXP dsyMatrix_as_dspMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dspMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_pack(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL( GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	NUN);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    SET_SLOT(val, Matrix_factorSym,
	     duplicate(GET_SLOT(from, Matrix_factorSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0

SEXP dsyMatrix_as_matrix(SEXP from, SEXP keep_dimnames)
{
    int n = INTEGER(GET_SLOT(from, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocMatrix(REALSXP, n, n));
    R_xlen_t nsqr = n; nsqr *= n;

    ddense_unpacked_make_symmetric(Memcpy(REAL(val),
					  REAL(GET_SLOT(from, Matrix_xSym)),
					  nsqr),
				   from);
    if(asLogical(keep_dimnames))
	setAttrib(val, R_DimNamesSymbol, get_symmetrized_DimNames(from, -1));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */


