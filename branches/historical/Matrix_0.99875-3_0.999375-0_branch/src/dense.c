#include "dense.h"
#include "Mutils.h"
#include "chm_common.h"

/**
 * Perform a left cyclic shift of columns j to k in the upper triangular
 * matrix x, then restore it to upper triangular form with Givens rotations.
 * The algorithm is based on the Fortran routine DCHEX from Linpack.
 *
 * The lower triangle of x is not modified.
 *
 * @param x Matrix stored in column-major order
 * @param ldx leading dimension of x
 * @param j column number (0-based) that will be shifted to position k
 * @param k last column number (0-based) to be shifted
 * @param cosines cosines of the Givens rotations
 * @param sines sines of the Givens rotations
 *
 * @return 0 for success
 */
static
int left_cyclic(double x[], int ldx, int j, int k,
		double cosines[], double sines[])
{
    double *lastcol;
    int i, jj;

    if (j >= k)
	error(_("incorrect left cyclic shift, j (%d) >= k (%d)"), j, k);
    if (j < 0)
	error(_("incorrect left cyclic shift, j (%d) < 0"), j, k);
    if (ldx < k)
	error(_("incorrect left cyclic shift, k (%d) > ldx (%d)"), k, ldx);
    lastcol = (double*) R_alloc(k+1, sizeof(double));
				/* keep a copy of column j */
    for(i = 0; i <= j; i++) lastcol[i] = x[i + j*ldx];
				/* For safety, zero the rest */
    for(i = j+1; i <= k; i++) lastcol[i] = 0.;
    for(jj = j+1; jj <= k; jj++) { /* columns to be shifted */
	int diagind = jj*(ldx+1), ind = (jj-j) - 1;
	double tmp = x[diagind], cc, ss;
				/* Calculate the Givens rotation. */
				/* This modified the super-diagonal element */
	F77_CALL(drotg)(x + diagind-1, &tmp, cosines + ind, sines + ind);
	cc = cosines[ind]; ss = sines[ind];
				/* Copy column jj+1 to column jj. */
	for(i = 0; i < jj; i++) x[i + (jj-1)*ldx] = x[i+jj*ldx];
				/* Apply rotation to columns up to k */
	for(i = jj; i < k; i++) {
	    tmp = cc*x[(jj-1)+i*ldx] + ss*x[jj+i*ldx];
	    x[jj+i*ldx] = cc*x[jj+i*ldx] - ss*x[(jj-1)+i*ldx];
	    x[(jj-1)+i*ldx] = tmp;
	}
				/* Apply rotation to lastcol */
	lastcol[jj] = -ss*lastcol[jj-1]; lastcol[jj-1] *= cc;
    }
				/* Copy lastcol to column k */
    for(i = 0; i <= k; i++) x[i+k*ldx] = lastcol[i];
    return 0;
}

static
SEXP getGivens(double x[], int ldx, int jmin, int rank)
{
    int shiftlen = (rank - jmin) - 1;
    SEXP ans = PROTECT(allocVector(VECSXP, 4)), nms, cosines, sines;

    SET_VECTOR_ELT(ans, 0, ScalarInteger(jmin));
    SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
    SET_VECTOR_ELT(ans, 2, cosines = allocVector(REALSXP, shiftlen));
    SET_VECTOR_ELT(ans, 3, sines = allocVector(REALSXP, shiftlen));
    setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("jmin"));
    SET_STRING_ELT(nms, 1, mkChar("rank"));
    SET_STRING_ELT(nms, 2, mkChar("cosines"));
    SET_STRING_ELT(nms, 3, mkChar("sines"));
    if (left_cyclic(x, ldx, jmin, rank - 1, REAL(cosines), REAL(sines)))
	error(_("Unknown error in getGivens"));
    UNPROTECT(1);
    return ans;
}

SEXP checkGivens(SEXP X, SEXP jmin, SEXP rank)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 2)),
	Xcp = PROTECT(duplicate(X));
    int  *Xdims;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
    SET_VECTOR_ELT(ans, 1, getGivens(REAL(Xcp), Xdims[0],
				     asInteger(jmin), asInteger(rank)));
    SET_VECTOR_ELT(ans, 0, Xcp);
    UNPROTECT(2);
    return ans;
}

SEXP lsq_dense_Chol(SEXP X, SEXP y)
{
    SEXP ans;
    int info, n, p, k, *Xdims, *ydims;
    double *xpx, d_one = 1., d_zero = 0.;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
    n = Xdims[0];
    p = Xdims[1];
    if (!(isReal(y) & isMatrix(y)))
	error(_("y must be a numeric (double precision) matrix"));
    ydims = INTEGER(coerceVector(getAttrib(y, R_DimSymbol), INTSXP));
    if (ydims[0] != n)
	error(_(
	    "number of rows in y (%d) does not match number of rows in X (%d)"),
	    ydims[0], n);
    k = ydims[1];
    if (k < 1 || p < 1) return allocMatrix(REALSXP, p, k);
    ans = PROTECT(allocMatrix(REALSXP, p, k));
    F77_CALL(dgemm)("T", "N", &p, &k, &n, &d_one, REAL(X), &n, REAL(y), &n,
		    &d_zero, REAL(ans), &p);
    xpx = (double *) R_alloc(p * p, sizeof(double));
    F77_CALL(dsyrk)("U", "T", &p, &n, &d_one, REAL(X), &n, &d_zero,
		    xpx, &p);
    F77_CALL(dposv)("U", &p, &k, xpx, &p, REAL(ans), &p, &info);
    if (info) error(_("Lapack routine dposv returned error code %d"), info);
    UNPROTECT(1);
    return ans;
}


SEXP lsq_dense_QR(SEXP X, SEXP y)
{
    SEXP ans;
    int info, n, p, k, *Xdims, *ydims, lwork;
    double *work, tmp, *xvals;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
    n = Xdims[0];
    p = Xdims[1];
    if (!(isReal(y) & isMatrix(y)))
	error(_("y must be a numeric (double precision) matrix"));
    ydims = INTEGER(coerceVector(getAttrib(y, R_DimSymbol), INTSXP));
    if (ydims[0] != n)
	error(_(
	    "number of rows in y (%d) does not match number of rows in X (%d)"),
	    ydims[0], n);
    k = ydims[1];
    if (k < 1 || p < 1) return allocMatrix(REALSXP, p, k);
    xvals = (double *) R_alloc(n * p, sizeof(double));
    Memcpy(xvals, REAL(X), n * p);
    ans = PROTECT(duplicate(y));
    lwork = -1;
    F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
		    &tmp, &lwork, &info);
    if (info)
	error(_("First call to Lapack routine dgels returned error code %d"),
	      info);
    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
		    work, &lwork, &info);
    if (info)
	error(_("Second call to Lapack routine dgels returned error code %d"),
	      info);
    UNPROTECT(1);
    return ans;
}

SEXP lapack_qr(SEXP Xin, SEXP tl)
{
    SEXP ans, Givens, Gcpy, nms, pivot, qraux, X;
    int i, n, nGivens = 0, p, trsz, *Xdims, rank;
    double rcond = 0., tol = asReal(tl), *work;

    if (!(isReal(Xin) & isMatrix(Xin)))
	error(_("X must be a real (numeric) matrix"));
    if (tol < 0.) error(_("tol, given as %g, must be non-negative"), tol);
    if (tol > 1.) error(_("tol, given as %g, must be <= 1"), tol);
    ans = PROTECT(allocVector(VECSXP,5));
    SET_VECTOR_ELT(ans, 0, X = duplicate(Xin));
    Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
    n = Xdims[0]; p = Xdims[1];
    SET_VECTOR_ELT(ans, 2, qraux = allocVector(REALSXP, (n < p) ? n : p));
    SET_VECTOR_ELT(ans, 3, pivot = allocVector(INTSXP, p));
    for (i = 0; i < p; i++) INTEGER(pivot)[i] = i + 1;
    trsz = (n < p) ? n : p;	/* size of triangular part of decomposition */
    rank = trsz;
    Givens = PROTECT(allocVector(VECSXP, rank - 1));
    setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 5));
    SET_STRING_ELT(nms, 0, mkChar("qr"));
    SET_STRING_ELT(nms, 1, mkChar("rank"));
    SET_STRING_ELT(nms, 2, mkChar("qraux"));
    SET_STRING_ELT(nms, 3, mkChar("pivot"));
    SET_STRING_ELT(nms, 4, mkChar("Givens"));
    if (n > 0 && p > 0) {
	int  info, *iwork, lwork;
	double *xpt = REAL(X), tmp;

	lwork = -1;
	F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), &tmp, &lwork, &info);
	if (info)
	    error(_("First call to dgeqrf returned error code %d"), info);
	lwork = (int) tmp;
	work = (double *) R_alloc((lwork < 3*trsz) ? 3*trsz : lwork,
				  sizeof(double));
	F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), work, &lwork, &info);
	if (info)
	    error(_("Second call to dgeqrf returned error code %d"), info);
	iwork = (int *) R_alloc(trsz, sizeof(int));
	F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			 work, iwork, &info);
	if (info)
	    error(_("Lapack routine dtrcon returned error code %d"), info);
	while (rcond < tol) {	/* check diagonal elements */
	    double minabs = (xpt[0] < 0.) ? -xpt[0]: xpt[0];
	    int jmin = 0;
	    for (i = 1; i < rank; i++) {
		double el = xpt[i*(n+1)];
		el = (el < 0.) ? -el: el;
		if (el < minabs) {
		    jmin = i;
		    minabs = el;
		}
	    }
	    if (jmin < (rank - 1)) {
		SET_VECTOR_ELT(Givens, nGivens, getGivens(xpt, n, jmin, rank));
		nGivens++;
	    }
	    rank--;
	    F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			     work, iwork, &info);
	    if (info)
		error(_("Lapack routine dtrcon returned error code %d"), info);
	}
    }
    SET_VECTOR_ELT(ans, 4, Gcpy = allocVector(VECSXP, nGivens));
    for (i = 0; i < nGivens; i++)
	SET_VECTOR_ELT(Gcpy, i, VECTOR_ELT(Givens, i));
    SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
    setAttrib(ans, install("useLAPACK"), ScalarLogical(1));
    setAttrib(ans, install("rcond"), ScalarReal(rcond));
    UNPROTECT(2);
    return ans;
}

SEXP dense_to_Csparse(SEXP x)
{
    CHM_DN chxd = AS_CHM_DN(PROTECT(mMatrix_as_geMatrix(x)));
    R_CheckStack();

    /* cholmod_dense_to_sparse() in CHOLMOD/Core/ below does only work for
       "REAL" 'xtypes', i.e. *not* for "nMatrix".
       ===> need "_x" in above call.

       Also it cannot keep symmetric / triangular, hence the
       as_geMatrix() above.  Note that this is already a *waste* for
       symmetric matrices; However, we could conceivably use an
       enhanced cholmod_dense_to_sparse(), with an extra boolean
       argument for symmetry.
    */
    CHM_SP chxs = cholmod_dense_to_sparse(chxd, 1, &c);
    int Rkind = (chxd->xtype == CHOLMOD_REAL) ? Real_KIND(x) : 0;

    UNPROTECT(1);
    /* chm_sparse_to_SEXP() *could* deal with symmetric
     * if chxs had such an stype; and we should be able to use uplo below */
    return chm_sparse_to_SEXP(chxs, 1, 0/*TODO: uplo_P(x) if x has an uplo slot*/,
			      Rkind, "",
			      isMatrix(x) ? getAttrib(x, R_DimNamesSymbol)
			      : GET_SLOT(x, Matrix_DimNamesSym));
}


/* FIXME: generalize this to  dense_band() : */

SEXP ddense_band(SEXP x, SEXP k1P, SEXP k2P)
/* Always returns a full matrix with entries outside the band zeroed
 * Class of the value can be dtrMatrix or dgeMatrix
 */
{
    SEXP aa, ans = PROTECT(dup_mMatrix_as_dgeMatrix(x));
    int *adims = INTEGER(GET_SLOT(ans, Matrix_DimSym)),
	i, j, k1 = asInteger(k1P), k2 = asInteger(k2P);
    int m = adims[0], n = adims[1], sqr = (adims[0] == adims[1]),
	tru = (k1 >= 0), trl = (k2 <= 0);
    double *ax = REAL(GET_SLOT(ans, Matrix_xSym));

    if (k1 > k2)
	error(_("Lower band %d > upper band %d"), k1, k2);
    for (j = 0; j < n; j++) {
	int i1 = j - k2, i2 = j + 1 - k1;
	for (i = 0; i < i1; i++) ax[i + j * m] = 0.;
	for (i = i2; i < m; i++) ax[i + j * m] = 0.;
    }
    if (!sqr || (!tru && !trl)) { /* return the dgeMatrix */
	UNPROTECT(1);
	return ans;
    }
    /* Copy ans to a dtrMatrix object (must be square) */
    /* Because slots of ans are freshly allocated and ans will not be
     * used, we use the slots themselves and don't duplicate */
    aa = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix")));
    SET_SLOT(aa, Matrix_xSym, GET_SLOT(ans, Matrix_xSym));
    SET_SLOT(aa, Matrix_DimSym, GET_SLOT(ans, Matrix_DimSym));
    SET_SLOT(aa, Matrix_DimNamesSym, GET_SLOT(ans, Matrix_DimNamesSym));
    SET_SLOT(aa, Matrix_diagSym, mkString("N"));
    SET_SLOT(aa, Matrix_uploSym, mkString(tru ? "U" : "L"));
    UNPROTECT(2);
    return aa;
}
