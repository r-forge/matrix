#include "LU.h"

/* FIXME - add the permutation matrix P to the result */
SEXP LU_expand(SEXP x)
{
    SEXP L = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))),
	U = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))),
	val = PROTECT(allocVector(VECSXP, 2)),
	nms = PROTECT(allocVector(STRSXP, 2)),
	lux = GET_SLOT(x, Matrix_xSym),
	dd = GET_SLOT(x, Matrix_DimSym);

    SET_STRING_ELT(nms, 0, mkChar("L"));
    SET_STRING_ELT(nms, 1, mkChar("U"));
    setAttrib(val, R_NamesSymbol, nms);
    SET_VECTOR_ELT(val, 0, L);
    SET_VECTOR_ELT(val, 1, U);
    SET_SLOT(L, Matrix_xSym, duplicate(lux));
    SET_SLOT(L, Matrix_DimSym, dd);
    SET_SLOT(L, Matrix_uploSym, mkString("L"));
    SET_SLOT(L, Matrix_diagSym, mkString("U"));
    SET_SLOT(L, Matrix_rcondSym, allocVector(REALSXP, 0));
    SET_SLOT(L, Matrix_factorSym, allocVector(VECSXP, 0));
    make_array_triangular(REAL(GET_SLOT(L, Matrix_xSym)), L);
    SET_SLOT(U, Matrix_xSym, duplicate(lux));
    SET_SLOT(U, Matrix_DimSym, dd);
    SET_SLOT(U, Matrix_uploSym, mkString("U"));
    SET_SLOT(U, Matrix_diagSym, mkString("N"));
    SET_SLOT(U, Matrix_rcondSym, allocVector(REALSXP, 0));
    SET_SLOT(U, Matrix_factorSym, allocVector(VECSXP, 0));
    make_array_triangular(REAL(GET_SLOT(U, Matrix_xSym)), U);
    UNPROTECT(4);
    return val;
}
