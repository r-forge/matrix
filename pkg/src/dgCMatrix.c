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
	return
	    mkString(_("last element of slot p must match length of slot i"));
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
    if (!cs_lusol(/*order*/ 1, xc, REAL(ycp), /*tol*/ 1e-7))
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
    if (!cs_qrsol(/*order*/ 1, xc, REAL(ycp)))
	error(_("cs_qrsol failed"));
    Free(xc);
    UNPROTECT(1);
    return ycp;
}

/* Modified version of Tim Davis's cs_qr_mex.c file for MATLAB */
SEXP dgCMatrix_QR(SEXP Ap, SEXP order)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sparseQR")));
    cs *A = Matrix_as_cs(Ap), *D;
    css *S;
    csn *N;
    int m = A->m, n = A->n, ord = asLogical(order) ? 3 : 0, *p;

    if (m < n) error("A must have # rows >= # columns") ;
    S = cs_sqr(ord, A, 1);	/* symbolic QR ordering & analysis*/
    if (!S) error("cs_sqr failed");
    N = cs_qr(A, S);		/* numeric QR factorization */
    if (!N) error("cs_qr failed") ;
    cs_dropzeros(N->L);		/* drop zeros from V and sort */
    D = cs_transpose(N->L, 1);
    cs_spfree(N->L);
    N->L = cs_transpose(D, 1);
    cs_spfree(D);
    cs_dropzeros(N->U);		/* drop zeros from R and sort */
    D = cs_transpose(N->U, 1);
    cs_spfree(N->U) ;
    N->U = cs_transpose(D, 1);
    cs_spfree(D);
    m = N->L->m;		/* m may be larger now */
    p = cs_pinv(S->pinv, m);	/* p = pinv' */
    SET_SLOT(ans, install("V"),
	     Matrix_cs_to_SEXP(N->L, "dgCMatrix", 0));
    Memcpy(REAL(ALLOC_SLOT(ans, install("beta"),
			   REALSXP, n)), N->B, n);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym,
			      INTSXP, m)), p, m);
    SET_SLOT(ans, install("R"),
	     Matrix_cs_to_SEXP(N->U, "dgCMatrix", 0));
    if (ord)
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("q"),
				  INTSXP, n)), S->q, n);
    else
	ALLOC_SLOT(ans, install("q"), INTSXP, 0);
    cs_nfree(N);
    cs_sfree(S);
    cs_free(p);
    UNPROTECT(1);
    return ans;
}

/* Modified version of Tim Davis's cs_qr_mex.c file for MATLAB */
SEXP dgCMatrix_LU(SEXP Ap, SEXP orderp, SEXP tolp)
{
    SEXP ans = get_factors(Ap, "LU");
    cs *A, *D;
    css *S;
    csn *N;
    int n, order = asInteger(orderp), *p;
    double tol = asReal(tolp);

    /* FIXME: dgCMatrix_LU should check ans for consistency in
     * permutation type with the requested value - Should have two
     * classes or two different names in the factors list for LU with
     * permuted columns or not. */

    if (ans != R_NilValue) return ans;
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sparseLU")));
    A = Matrix_as_cs(Ap);
    n = A->n;
    if (A->m != n)
	error("LU decomposition applies only to square matrices");
    if (order) {		/* not using natural order */
	order = (tol == 1) ? 2	/* amd(S'*S) w/dense rows or I */
	    : 1;		/* amd (A+A'), or natural */
    }
    S = cs_sqr (order, A, 0) ;	/* symbolic ordering, no QR bound */
    N = cs_lu (A, S, tol) ;	/* numeric factorization */
    if (!N) error ("cs_lu failed (singular, or out of memory)") ;
    cs_dropzeros (N->L) ;	/* drop zeros from L and sort it */
    D = cs_transpose (N->L, 1) ;
    cs_spfree (N->L) ;
    N->L = cs_transpose (D, 1) ;
    cs_spfree (D) ;
    cs_dropzeros (N->U) ;	/* drop zeros from U and sort it */
    D = cs_transpose (N->U, 1) ;
    cs_spfree (N->U) ;
    N->U = cs_transpose (D, 1) ;
    cs_spfree (D) ;
    p = cs_pinv (N->pinv, n) ;	/* p=pinv' */
    SET_SLOT(ans, install("L"),
	     Matrix_cs_to_SEXP(N->L, "dgCMatrix", 0));
    SET_SLOT(ans, install("U"),
	     Matrix_cs_to_SEXP(N->U, "dgCMatrix", 0));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym,
			      INTSXP, n)), p, n);
    if (order)
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("q"),
				  INTSXP, n)), S->q, n);
    cs_nfree(N);
    cs_sfree(S);
    cs_free(p);
    Free(A);
    UNPROTECT(1);
    return set_factors(Ap, ans, "LU");
}

SEXP dgCMatrix_matrix_solve(SEXP Ap, SEXP b)
{
    SEXP ans = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    SEXP lu = dgCMatrix_LU(Ap, ScalarLogical(1), ScalarReal(1));
    SEXP qslot = GET_SLOT(lu, install("q"));
    cs *L = Matrix_as_cs(GET_SLOT(lu, install("L"))),
	*U = Matrix_as_cs(GET_SLOT(lu, install("U")));
    int *bdims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    int j, n = bdims[0], nrhs = bdims[1];
    int *p = INTEGER(GET_SLOT(lu, Matrix_pSym)),
	*q = LENGTH(qslot) ? INTEGER(qslot) : (int *) NULL;
    double *ax = REAL(GET_SLOT(ans, Matrix_xSym)),
	*x = Calloc(n, double);

    if (U->n != n || nrhs < 1 || n < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    for (j = 0; j < nrhs; j++) {
	cs_pvec(p, ax + j * n, x, n);  /* x = b(p) */
	cs_lsolve(L, x);	       /* x = L\x */
	cs_usolve(U, x);	       /* x = U\x */
	if (q)			       /* b(q) = x */
	    cs_ipvec(q, x, ax + j * n, n);
	else
	    Memcpy(ax + j * n, x, n);
    }
    Free(L); Free(U); Free(x);
    UNPROTECT(1);
    return ans;
}

