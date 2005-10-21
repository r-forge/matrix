#include "ldense.h"

/* this is very close to dspMatrix_as_dsy* () in ./dspMatrix.c : */
SEXP lspMatrix_as_lsyMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("lsyMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    packed_to_full_int(LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, n*n)),
		       LOGICAL( GET_SLOT(from, Matrix_xSym)), n,
		       *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW);
    UNPROTECT(1);
    return val;
}

/* this is very close to dsyMatrix_as_dsp* () in ./dsyMatrix.c : */
SEXP lsyMatrix_as_lspMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("lspMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    full_to_packed_int(
	LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, (n*(n+1))/2)),
	LOGICAL( GET_SLOT(from, Matrix_xSym)), n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW, NUN);
    UNPROTECT(1);
    return val;
}

/* this is very close to dtpMatrix_as_dtr* () in ./dtpMatrix.c : */
SEXP ltpMatrix_as_ltrMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("ltrMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    packed_to_full_int(LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, n*n)),
		       LOGICAL(GET_SLOT(from, Matrix_xSym)), n,
		       *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW);
    UNPROTECT(1);
    return val;
}

/* this is very close to dtrMatrix_as_dtp* () in ./dtrMatrix.c : */
SEXP ltrMatrix_as_ltpMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("ltpMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    full_to_packed_int(
	LOGICAL(ALLOC_SLOT(val, Matrix_xSym, LGLSXP, (n*(n+1))/2)),
	LOGICAL(GET_SLOT(from, Matrix_xSym)), n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	*CHAR(STRING_ELT(diag, 0)) == 'U' ? UNT : NUN);
    UNPROTECT(1);
    return val;
}
