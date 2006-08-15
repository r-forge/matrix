#include "dgCMatrix.h"

#include "chm_common.h"

/* FIXME -- we "forget" about dimnames almost everywhere : */

SEXP dgCMatrix_validate(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym),
	xslot = GET_SLOT(x, Matrix_xSym);
    int j,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow = dims[0],
	ncol = dims[1],
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    if (length(islot) != length(xslot))
	return mkString(_("lengths of slots i and x must match"));
    if (length(pslot) != ncol + 1)
	return mkString(_("slot p must have length ncol + 1"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slots i and x"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi))
	csc_sort_columns(ncol, xp, xi, REAL(xslot));

    return ScalarLogical(1);
}


SEXP compressed_to_dgTMatrix(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    SEXP indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix"))),
	indP = GET_SLOT(x, indSym),
	pP = GET_SLOT(x, Matrix_pSym);
    int npt = length(pP) - 1;

    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_xSym,  duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(ans, indSym, duplicate(indP));
    expand_cmprPt(npt, INTEGER(pP),
		  INTEGER(ALLOC_SLOT(ans, col ? Matrix_jSym : Matrix_iSym,
				     INTSXP, length(indP))));
    UNPROTECT(1);
    return ans;
}

SEXP compressed_non_0_ij(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    SEXP ans, indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP indP = GET_SLOT(x, indSym),
	pP = GET_SLOT(x, Matrix_pSym);
    int n_el = length(indP), i, *ij;

    ij = INTEGER(ans = PROTECT(allocMatrix(INTSXP, n_el, 2)));
    /* expand the compressed margin to 'i' or 'j' : */
    expand_cmprPt(length(pP) - 1, INTEGER(pP), &ij[col ? n_el : 0]);
    /* and copy the other one: */
    if (col)
	for(i = 0; i < n_el; i++)
	    ij[i] = INTEGER(indP)[i];
    else /* row compressed */
	for(i = 0; i < n_el; i++)
	    ij[i + n_el] = INTEGER(indP)[i];

    UNPROTECT(1);
    return ans;
}

SEXP csc_matrix_mm(SEXP a, SEXP b, SEXP classed, SEXP right)
{
    int cl = asLogical(classed), rt = asLogical(right);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*ai = INTEGER(GET_SLOT(a, Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(a, Matrix_pSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol)),
	*cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2)),
	chk, ione = 1, j, jj, k, m, n;
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)),
	*bx = REAL(cl ? GET_SLOT(b, Matrix_xSym) : b), *cx;

    if (rt) {
	m = bdims[0]; n = adims[1]; k = bdims[1]; chk = adims[0];
    } else {
	m = adims[0]; n = bdims[1]; k = adims[1]; chk = bdims[0];
    }
    if (chk != k)
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1 || k < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    cx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n));
    AZERO(cx, m * n); /* zero the accumulators */
    for (j = 0; j < n; j++) { /* across columns of c */
	if (rt) {
	    int kk, k2 = ap[j + 1];
	    for (kk = ap[j]; kk < k2; kk++) {
		F77_CALL(daxpy)(&m, &ax[kk], &bx[ai[kk]*m],
				&ione, &cx[j*m], &ione);
	    }
	} else {
	    double *ccol = cx + j * m,
		*bcol = bx + j * k;

	    for (jj = 0; jj < k; jj++) { /* across columns of a */
		int kk, k2 = ap[jj + 1];
		for (kk = ap[jj]; kk < k2; kk++) {
		    ccol[ai[kk]] += ax[kk] * bcol[jj];
		}
	    }
	}
    }
    cdims[0] = m; cdims[1] = n;
    UNPROTECT(1);
    return val;
}

SEXP dgCMatrix_lusol(SEXP x, SEXP y)
{
    SEXP ycp = PROTECT(duplicate(y));
    cs *xc = Matrix_as_cs(x);

    if (xc->m != xc->n || xc->m <= 0)
	error(_("dgCMatrix_lusol requires a square, non-empty matrix"));
    if (!isReal(ycp) || LENGTH(ycp) != xc->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    if (!cs_lusol(xc, REAL(ycp), /*order*/ 1, /*tol*/ 1e-7))
	error(_("cs_lusol failed"));
    Free(xc);
    UNPROTECT(1);
    return ycp;
}

SEXP dgCMatrix_qrsol(SEXP x, SEXP y)
{
    SEXP ycp = PROTECT(duplicate(y));
    cs *xc = Matrix_as_cs(x);

    if (xc->m < xc->n || xc->n <= 0)
	error(_("dgCMatrix_qrsol requires a 'tall' rectangular matrix"));
    if (!isReal(ycp) || LENGTH(ycp) != xc->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    if (!cs_qrsol(xc, REAL(ycp), /*order*/ 1))
	error(_("cs_qrsol failed"));
    Free(xc);
    UNPROTECT(1);
    return ycp;
}

/* Patterned directly on Tim Davis's cs_qr_mex.c file for MATLAB */
SEXP dgCMatrix_QR(SEXP Ap, SEXP order)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sparseQR")));
    cs *A = Matrix_as_cs(Ap), *D;
    css *S;
    csn *N;
    int ord = asLogical(order);
    int m = A->m, n = A->n, *p;

    if (m < n) error("A must have # rows >= # columns") ;
    S = cs_sqr (A, ord, 1);	/* symbolic QR ordering & analysis*/
    N = cs_qr (A, S);		/* numeric QR factorization */
    if (!N) error("cs_qr failed") ;
    cs_dropzeros (N->L);		    /* drop zeros from V and sort */
    D = cs_transpose (N->L, 1);
    cs_spfree (N->L);
    N->L = cs_transpose (D, 1);
    cs_spfree (D);
    cs_dropzeros (N->U);		    /* drop zeros from R and sort */
    D = cs_transpose (N->U, 1);
    cs_spfree (N->U) ;
    N->U = cs_transpose (D, 1);
    cs_spfree (D);
    m = N->L->m;				    /* m may be larger now */
    Rprintf("m = %d, n = %d\n", m, n);
    p = cs_pinv (S->Pinv, m);			    /* p = pinv' */
    SET_SLOT(ans, install("V"),			    /* return V */
	     Matrix_cs_to_SEXP(N->L, "dgCMatrix", 0));
    Memcpy(REAL(ALLOC_SLOT(ans, install("beta"),    /* return beta */
			   REALSXP, n)), N->B, n);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym,     /* return p */
			      INTSXP, m)), p, m);
    SET_SLOT(ans, install("R"),                     /* return R */
	     Matrix_cs_to_SEXP(N->U, "dgCMatrix", 0));
    Memcpy(INTEGER(ALLOC_SLOT(ans, install("q"),    /* return q */
			      INTSXP, n)), S->Q, n);
    cs_nfree(N);
    cs_sfree(S);
    cs_free(p);
    UNPROTECT(1);
    return ans;
}

