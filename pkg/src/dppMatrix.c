#include "dppMatrix.h"

SEXP dppMatrix_validate(SEXP obj)
{
/*     int i, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0]; */
/*     double *x = REAL(GET_SLOT(obj, Matrix_xSym)); */

    /* quick but nondefinitive check on positive definiteness */
/*     for (i = 0; i < n; i++) */
/* 	if (x[i * np1] < 0) */
/* 	    return mkString(_("dppMatrix is not positive definite")); */
    return ScalarLogical(1);
}

SEXP dppMatrix_chol(SEXP x)
{
    SEXP val = get_factors(x, "pCholesky"),
	dimP = GET_SLOT(x, Matrix_DimSym),
	uploP = GET_SLOT(x, Matrix_uploSym);
    char *uplo = CHAR(STRING_ELT(uploP, 0));
    int *dims = INTEGER(dimP), info;

    if (val != R_NilValue) return val;
    dims = INTEGER(dimP);
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("pCholesky")));
    SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    F77_CALL(dpptrf)(uplo, dims, REAL(GET_SLOT(val, Matrix_xSym)), &info);
    UNPROTECT(1);
    return set_factors(x, val, "pCholesky");
}

static
double set_rcond(SEXP obj, char *typstr)
{
    char typnm[] = {'O', '\0'};	/* always use the one norm */
    SEXP rcv = GET_SLOT(obj, Matrix_rcondSym);
    double rcond = get_double_by_name(rcv, typnm);

    if (R_IsNA(rcond)) {
        SEXP Chol = dppMatrix_chol(obj);
	int *dims = INTEGER(GET_SLOT(Chol, Matrix_DimSym)), info;
	double anorm = get_norm_sp(obj, typnm);

	F77_CALL(dppcon)(CHAR(asChar(GET_SLOT(Chol, Matrix_uploSym))),
			 dims, REAL(GET_SLOT(Chol, Matrix_xSym)),
			 &anorm, &rcond,
			 (double *) R_alloc(3*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, Matrix_rcondSym,
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP dppMatrix_rcond(SEXP obj, SEXP type)
{
    return ScalarReal(set_rcond(obj, CHAR(asChar(type))));
}

SEXP dppMatrix_solve(SEXP x)
{
    SEXP Chol = dppMatrix_chol(x);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dppMatrix")));
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), info;

    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(Chol, Matrix_uploSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(Chol, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(Chol, Matrix_DimSym)));
    F77_CALL(dpptri)(CHAR(asChar(GET_SLOT(val, Matrix_uploSym))),
		     dims, REAL(GET_SLOT(val, Matrix_xSym)), &info);
    UNPROTECT(1);
    return val;
}

SEXP dppMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP Chol = dppMatrix_chol(a),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = (cl ? INTEGER(GET_SLOT(b, Matrix_DimSym)) :
		  INTEGER(getAttrib(b, R_DimSymbol)));
    int n = bdims[0], nrhs = bdims[1], info;
    int sz = n * nrhs;

    if (!cl && !(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    Memcpy(INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2)), bdims, 2);
    F77_CALL(dpptrs)
	(CHAR(asChar(GET_SLOT(Chol, Matrix_uploSym))), &n, &nrhs,
	 REAL(GET_SLOT(Chol, Matrix_xSym)),
	 Memcpy(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, sz)),
		REAL(cl ? GET_SLOT(b, Matrix_xSym) : b), sz), &n, &info);
    UNPROTECT(1);
    return val;
}
