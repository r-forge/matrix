#include "dspMatrix.h"

SEXP dspMatrix_validate(SEXP obj)
{
    SEXP val;
    int *Dim = INTEGER(GET_SLOT(obj, Matrix_DimSym));

    if (isString(val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym),
					   "LU", "uplo"))) return val;
    if (Dim[0] != Dim[1])
	return mkString(_("Symmetric matrix must be square"));
    return ScalarLogical(1);
}

double get_norm_sp(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    double *work = (double *) NULL;

    typnm[0] = norm_type(typstr);
    if (*typnm == 'I' || *typnm == 'O') {
        work = (double *) R_alloc(dims[0], sizeof(double));
    }
    return F77_CALL(dlansp)(typnm,
			    CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))),
			    dims, REAL(GET_SLOT(obj, Matrix_xSym)),
			    work);
}

SEXP dspMatrix_norm(SEXP obj, SEXP type)
{
    return ScalarReal(get_norm_sp(obj, CHAR(asChar(type))));
}

static
double set_rcond_sp(SEXP obj, char *typstr)
{
    char typnm[] = {'\0', '\0'};
    SEXP rcv = GET_SLOT(obj, Matrix_rcondSym);
    double rcond;

    typnm[0] = rcond_type(typstr);
    rcond = get_double_by_name(rcv, typnm);

    if (R_IsNA(rcond)) {
	SEXP trf = dspMatrix_trf(obj);
	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)), info;
	double anorm = get_norm_sp(obj, "O");

	F77_CALL(dspcon)(CHAR(asChar(GET_SLOT(trf, Matrix_uploSym))),
			 dims, REAL(GET_SLOT(trf, Matrix_xSym)),
			 INTEGER(GET_SLOT(trf, Matrix_permSym)),
			 &anorm, &rcond,
			 (double *) R_alloc(2*dims[0], sizeof(double)),
			 (int *) R_alloc(dims[0], sizeof(int)), &info);
	SET_SLOT(obj, Matrix_rcondSym,
		 set_double_by_name(rcv, rcond, typnm));
    }
    return rcond;
}

SEXP dspMatrix_rcond(SEXP obj, SEXP type)
{
    return ScalarReal(set_rcond_sp(obj, CHAR(asChar(type))));
}

SEXP dspMatrix_solve(SEXP a)
{
    SEXP trf = dspMatrix_trf(a);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dspMatrix")));
    int *dims = INTEGER(GET_SLOT(trf, Matrix_DimSym)), info;

    SET_SLOT(val, Matrix_uploSym, duplicate(GET_SLOT(trf, Matrix_uploSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(trf, Matrix_xSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(trf, Matrix_DimSym)));
    SET_SLOT(val, Matrix_rcondSym, duplicate(GET_SLOT(a, Matrix_rcondSym)));
    F77_CALL(dsptri)(CHAR(asChar(GET_SLOT(val, Matrix_uploSym))),
		     dims, REAL(GET_SLOT(val, Matrix_xSym)),
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     (double *) R_alloc((long) dims[0], sizeof(double)),
		     &info);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_dgeMatrix_solve(SEXP a, SEXP b)
{
    SEXP trf = dspMatrix_trf(a),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	info;

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    SET_SLOT(val, Matrix_DimSym, duplicate(GET_SLOT(b, Matrix_DimSym)));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(b, Matrix_xSym)));
    F77_CALL(dsptrs)(CHAR(asChar(GET_SLOT(trf, Matrix_uploSym))),
		     adims, bdims + 1,
		     REAL(GET_SLOT(trf, Matrix_xSym)),
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     REAL(GET_SLOT(val, Matrix_xSym)),
		     bdims, &info);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP trf = dspMatrix_trf(a),
	val = PROTECT(duplicate(b));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol)),
	info;

    if (!(isReal(b) && isMatrix(b)))
	error(_("Argument b must be a numeric matrix"));
    if (*adims != *bdims || bdims[1] < 1 || *adims < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    F77_CALL(dsptrs)(CHAR(asChar(GET_SLOT(trf, Matrix_uploSym))),
		     adims, bdims + 1,
		     REAL(GET_SLOT(trf, Matrix_xSym)),
		     INTEGER(GET_SLOT(trf, Matrix_permSym)),
		     REAL(val), bdims, &info);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_as_dsyMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dsyMatrix"))),
	uplo = GET_SLOT(from, Matrix_uploSym),
	dimP = GET_SLOT(from, Matrix_DimSym),
	dmnP = GET_SLOT(from, Matrix_DimNamesSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_rcondSym,
	     duplicate(GET_SLOT(from, Matrix_rcondSym)));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_DimNamesSym, duplicate(dmnP));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    packed_to_full(REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n*n)),
		   REAL(GET_SLOT(from, Matrix_xSym)), n,
		   *CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_dgeMatrix_mm(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym));
    int i, ione = 1, n = adims[0], nrhs = bdims[1];
    SEXP val = PROTECT(duplicate(b));
    char *uplo = CHAR(STRING_ELT(GET_SLOT(a, Matrix_uploSym), 0));
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)),
	*vx = REAL(GET_SLOT(val, Matrix_xSym)), one = 1., zero = 0.;

    if (bdims[0] != n)
	error(_("Matrices are not conformable for multiplication"));
    if (nrhs < 1 || n < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    for (i = 0; i < nrhs; i++)
	F77_CALL(dspmv)(uplo, &n, &one, ax, vx + i * n, &ione,
			&zero, vx + i * n, &ione);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_matrix_mm(SEXP a, SEXP b)
{
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(getAttrib(b, R_DimSymbol));
    int i, ione = 1, n = adims[0], nrhs = bdims[1];
    SEXP val = PROTECT(duplicate(b));
    char *uplo = CHAR(STRING_ELT(GET_SLOT(a, Matrix_uploSym), 0));
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)),
	*vx = REAL(val), one = 1., zero = 0.;

    if (bdims[0] != n)
	error(_("Matrices are not conformable for multiplication"));
    if (nrhs < 1 || n < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    for (i = 0; i < nrhs; i++)
	F77_CALL(dspmv)(uplo, &n, &one, ax, vx + i * n, &ione,
			&zero, vx + i * n, &ione);
    UNPROTECT(1);
    return val;
}

SEXP dspMatrix_trf(SEXP x)
{
    SEXP val = get_factors(x, "pBunchKaufman"),
	dimP = GET_SLOT(x, Matrix_DimSym),
	uploP = GET_SLOT(x, Matrix_uploSym);
    int *dims = INTEGER(dimP), *perm, info;
    int n = dims[0];
    char *uplo = CHAR(STRING_ELT(uploP, 0));

    if (val != R_NilValue) return val;
    dims = INTEGER(dimP);
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("pBunchKaufman")));
    SET_SLOT(val, Matrix_uploSym, duplicate(uploP));
    SET_SLOT(val, Matrix_diagSym, mkString("N"));
    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_xSym, duplicate(GET_SLOT(x, Matrix_xSym)));
    perm = INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, n));
    F77_CALL(dsptrf)(uplo, dims, REAL(GET_SLOT(val, Matrix_xSym)), perm, &info);
    if (info) error(_("Lapack routine %s returned error code %d"), "dsptrf", info);
    UNPROTECT(1);
    return set_factors(x, val, "pBunchKaufman");
}

