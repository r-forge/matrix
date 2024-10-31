#include <Rinternals.h>
#include "version.h"

SEXP R_Matrix_version(void)
{
	SEXP ans, nms;
	PROTECT(ans = Rf_allocVector(INTSXP, 3));
	INTEGER(ans)[0] = MATRIX_PACKAGE_VERSION;
	INTEGER(ans)[1] = MATRIX_ABI_VERSION;
	INTEGER(ans)[2] = MATRIX_SUITESPARSE_VERSION;
	PROTECT(nms = Rf_allocVector(STRSXP, 3));
	SET_STRING_ELT(nms, 0, Rf_mkChar("package"));
	SET_STRING_ELT(nms, 1, Rf_mkChar("abi"));
	SET_STRING_ELT(nms, 2, Rf_mkChar("suitesparse"));
	Rf_setAttrib(ans, R_NamesSymbol, nms);
	UNPROTECT(2);
	return ans;
}
