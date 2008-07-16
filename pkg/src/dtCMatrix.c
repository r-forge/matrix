				/* Sparse triangular numeric matrices */
#include "dtCMatrix.h"
#include "cs_utils.h"

#define RETURN(_CH_)   UNPROTECT(1); return (_CH_);

/* This is used for *BOTH* triangular and symmetric Csparse: */
SEXP tCMatrix_validate(SEXP x)
{
    SEXP val = xCMatrix_validate(x);/* checks x slot */
    if(isString(val))
	return(val);
    else {
	SEXP
	    islot = GET_SLOT(x, Matrix_iSym),
	    pslot = GET_SLOT(x, Matrix_pSym);
	int uploT = (*uplo_P(x) == 'U'),
	    k, nnz = length(islot),
	    *xi = INTEGER(islot),
	    *xj = INTEGER(PROTECT(allocVector(INTSXP, nnz)));

	expand_cmprPt(length(pslot) - 1, INTEGER(pslot), xj);

	/* Maybe FIXME: ">" should be ">="	for diag = 'U' (uplo = 'U') */
	if(uploT) {
	    for (k = 0; k < nnz; k++)
		if(xi[k] > xj[k]) {
		    RETURN(mkString(_("uplo='U' must not have sparse entries below the diagonal")));
		}
	}
	else {
	    for (k = 0; k < nnz; k++)
		if(xi[k] < xj[k]) {
		    RETURN(mkString(_("uplo='L' must not have sparse entries above the diagonal")));
		}
	}

	RETURN(ScalarLogical(1));
    }
}

/* This is used for *BOTH* triangular and symmetric Rsparse: */
SEXP tRMatrix_validate(SEXP x)
{
    SEXP val = xRMatrix_validate(x);/* checks x slot */
    if(isString(val))
	return(val);
    else {
	SEXP
	    jslot = GET_SLOT(x, Matrix_jSym),
	    pslot = GET_SLOT(x, Matrix_pSym);
	int uploT = (*uplo_P(x) == 'U'),
	    k, nnz = length(jslot),
	    *xj = INTEGER(jslot),
	    *xi = INTEGER(PROTECT(allocVector(INTSXP, nnz)));

	expand_cmprPt(length(pslot) - 1, INTEGER(pslot), xi);

	/* Maybe FIXME: ">" should be ">="	for diag = 'U' (uplo = 'U') */
	if(uploT) {
	    for (k = 0; k < nnz; k++)
		if(xi[k] > xj[k]) {
		    RETURN(mkString(_("uplo='U' must not have sparse entries below the diagonal")));
		}
	}
	else {
	    for (k = 0; k < nnz; k++)
		if(xi[k] < xj[k]) {
		    RETURN(mkString(_("uplo='L' must not have sparse entries above the diagonal")));
		}
	}

	RETURN(ScalarLogical(1));
    }
}

#undef RETURN

SEXP dtCMatrix_solve(SEXP a)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    CSP A = AS_CSP(Csparse_diagU2N(a));
    int *bp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, (A->n) + 1)),
	bnz = 10 * A->n,	/* initial estimate of nnz in b */
	lo = uplo_P(a)[0] == 'L', top;
    /* These arrays must use Calloc because of possible Realloc */
    int *ti = Calloc(bnz, int), p, j, nz, pos = 0;
    double *tx = Calloc(bnz, double);
    cs *u = cs_spalloc(A->n, 1,1,1,0);	/* Sparse unit vector */
    double  *wrk = Alloca(A->n, double);
    int *xi = Alloca(2*A->n, int);	/* for cs_reach */
    R_CheckStack();

    slot_dup(ans, a, Matrix_DimSym);
    SET_DimNames(ans, a);
    slot_dup(ans, a, Matrix_uploSym);
    slot_dup(ans, a, Matrix_diagSym);
    /* initialize the "constant part" of the sparse unit vector */
    u->x[0] = 1.;
    u->p[0] = 0; u->p[1] = 1;
    bp[0] = 0;
    for (j = 0; j < A->n; j++) {
	u->i[0] = j;			/* u := j'th unit vector */
	/* (wrk[top:n],xi[top:n]) :=  A^{-1} u  : */
	top = cs_spsolve (A, u, 0, xi, wrk, 0, lo);
	nz = A->n - top;
	bp[j + 1] = nz + bp[j];
	if (bp[j + 1] > bnz) {
	    while (bp[j + 1] > bnz) bnz *= 2;
	    ti = Realloc(ti, bnz, int);
	    tx = Realloc(tx, bnz, double);
	}
	if (lo)
	    for(p = top; p < A->n; p++, pos++) {
		ti[pos] = xi[p];
		tx[pos] = wrk[xi[p]];
	    }
	else /* upper triagonal */
	    for(p = A->n - 1; p >= top; p--, pos++) {
		ti[pos] = xi[p];
		tx[pos] = wrk[xi[p]];
	    }
    }
    nz = bp[A->n];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)), ti, nz);
    Memcpy(   REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)), tx, nz);

    Free(ti); Free(tx); cs_spfree(u);
    UNPROTECT(1);
    return ans;
}

SEXP dtCMatrix_matrix_solve(SEXP a, SEXP b, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    CSP A = AS_CSP(Csparse_diagU2N(a));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol));
    int j, n = bdims[0], nrhs = bdims[1], lo = (*uplo_P(a) == 'L');
    double *bx;
    R_CheckStack();

    if (*adims != n || nrhs < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)), bdims, 2);
    /* FIXME: copy dimnames or Dimnames as well */
    bx = Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * nrhs)),
		REAL(cl ? GET_SLOT(b, Matrix_xSym):b), n * nrhs);
    for (j = 0; j < nrhs; j++)
	lo ? cs_lsolve(A, bx + n * j) : cs_usolve(A, bx + n * j);
    UNPROTECT(1);
    return ans;
}

SEXP dtCMatrix_sparse_solve(SEXP a, SEXP b)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    CSP A = AS_CSP(Csparse_diagU2N(a)), B = AS_CSP(Csparse_diagU2N(b));
    int *xp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, (B->n) + 1)),
	xnz = 10 * B->p[B->n],	/* initial estimate of nnz in x */
	lo = uplo_P(a)[0] == 'L', top;
    /* These arrays must use Calloc because of possible Realloc */
    int *ti = Calloc(xnz, int), p, j, nz, pos = 0;
    double *tx = Calloc(xnz, double);
    double  *wrk = Alloca(A->n, double);
    int *xi = Alloca(2*A->n, int);	/* for cs_reach */
    R_CheckStack();
    
    if (A->m != A->n || B->n < 1 || A->n < 1 || A->n != B->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    slot_dup(ans, b, Matrix_DimSym);
    SET_DimNames(ans, b);
    xp[0] = 0;
    for (j = 0; j < B->n; j++) {
	/* (wrk[top:n],xi[top:n]) :=  A^{-1} B[,j] */
	top = cs_spsolve (A, B, j, xi, wrk, 0, lo);
	nz = A->n - top;
	xp[j + 1] = nz + xp[j];
	if (xp[j + 1] > xnz) {
	    while (xp[j + 1] > xnz) xnz *= 2;
	    ti = Realloc(ti, xnz, int);
	    tx = Realloc(tx, xnz, double);
	}
	if (lo)
	    for(p = top; p < A->n; p++, pos++) {
		ti[pos] = xi[p];
		tx[pos] = wrk[xi[p]];
	    }
	else /* upper triagonal */
	    for(p = A->n - 1; p >= top; p--, pos++) {
		ti[pos] = xi[p];
		tx[pos] = wrk[xi[p]];
	    }
    }
    nz = xp[A->n];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)), ti, nz);
    Memcpy(   REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)), tx, nz);
    
    Free(ti); Free(tx);
    UNPROTECT(1);
    return ans;
}
