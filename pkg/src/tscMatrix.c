				/* Sparse triangular matrices */
#include "tscMatrix.h"

SEXP tsc_validate(SEXP x)
{
    return ScalarLogical(1);
}

SEXP tsc_transpose(SEXP x)
{
    SEXP
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("tscMatrix"))),
	islot = GET_SLOT(x, Matrix_iSym);
    int nnz = length(islot),
	*adims, *xdims = INTEGER(GET_SLOT(x, Matrix_DimSym));


    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    adims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    adims[0] = xdims[1]; adims[1] = xdims[0];
    if (toupper(CHAR(asChar(GET_SLOT(x, Matrix_uploSym)))[0]) == 'U')
	SET_SLOT(ans, Matrix_uploSym, ScalarString(mkChar("L")));
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, xdims[0] + 1));
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    csc_components_transpose(xdims[0], xdims[1], nnz,
			     INTEGER(GET_SLOT(x, Matrix_pSym)),
			     INTEGER(islot),
			     REAL(GET_SLOT(x, Matrix_xSym)),
			     INTEGER(GET_SLOT(ans, Matrix_pSym)),
			     INTEGER(GET_SLOT(ans, Matrix_iSym)),
			     REAL(GET_SLOT(ans, Matrix_xSym)));
    UNPROTECT(1);
    return ans;
}

SEXP tsc_to_triplet(SEXP x)
{
    SEXP ans;
    if (toupper(CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0))[0]) != 'U')
	ans = csc_to_triplet(x);
    else {			/* unit triangular matrix */
	SEXP islot = GET_SLOT(x, Matrix_iSym), 
	    pslot = GET_SLOT(x, Matrix_pSym);
	int *ai, *aj, j,
	    n = length(pslot) - 1,
	    nod = length(islot),
	    nout = n + nod,
	    *p = INTEGER(pslot);
	double *ax;
    
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("tripletMatrix")));
	SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
	SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nout));
	ai = INTEGER(GET_SLOT(ans, Matrix_iSym));
	Memcpy(ai, INTEGER(islot), nod);
	SET_SLOT(ans, Matrix_jSym, allocVector(INTSXP, nout));
	aj = INTEGER(GET_SLOT(ans, Matrix_jSym));
	SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nout));
	ax = REAL(GET_SLOT(ans, Matrix_xSym));
	Memcpy(ax, REAL(GET_SLOT(x, Matrix_xSym)), nod);
	for (j = 0; j < n; j++) {
	    int jj, npj = nod + j,  p2 = p[j+1];
	    aj[npj] = ai[npj] = j;
	    ax[npj] = 1.;
	    for (jj = p[j]; jj < p2; jj++) aj[jj] = j;
	}
	UNPROTECT(1);
    }
    return ans;
}

/** 
 * Derive the column pointer vector for the inverse of L from the parent array
 * 
 * @param n length of parent array
 * @param countDiag 0 for a unit triangular matrix with implicit diagonal, otherwise 1
 * @param pr parent vector describing the elimination tree
 * @param ap array of length n+1 to be filled with the column pointers
 * 
 * @return the number of non-zero entries (ap[n])
 */
int parent_inv_ap(int n, int countDiag, const int pr[], int ap[])
{
    int *sz = Calloc(n, int), j;

    for (j = n - 1; j >= 0; j--) {
	int parent = pr[j];
	sz[j] = (parent < 0) ?  countDiag : (1 + sz[parent]);
    }
    ap[0] = 0;
    for (j = 0; j < n; j++)
	ap[j+1] = ap[j] + sz[j];
    Free(sz);
    return ap[n];
}

/** 
 * Derive the row index array for the inverse of L from the parent array
 * 
 * @param n length of parent array
 * @param countDiag 0 for a unit triangular matrix with implicit diagonal, otherwise 1
 * @param pr parent vector describing the elimination tree
 * @param ai row index vector of length ap[n]
 */
void parent_inv_ai(int n, int countDiag, const int pr[], int ai[])
{
    int j, k, pos = 0;
    for (j = 0; j < n; j++) {
	if (countDiag) ai[pos++] = j;
	for (k = pr[j]; k >= 0; k = pr[k]) ai[pos++] = k;
    }
}
    
SEXP Parent_inverse(SEXP par, SEXP unitdiag)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("tscMatrix")));
    int *ap, *ai, *dims, *pr = INTEGER(par),
	countDiag = 1 - asLogical(unitdiag),
	j, n = length(par), nnz;
    double *ax;
    
    if (!isInteger(par)) error("par argument must be an integer vector");
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, n + 1));
    ap = INTEGER(GET_SLOT(ans, Matrix_pSym));
    nnz = parent_inv_ap(n, countDiag, pr, ap);
    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    ai = INTEGER(GET_SLOT(ans, Matrix_iSym));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    ax = REAL(GET_SLOT(ans, Matrix_xSym));
    for (j = 0; j < nnz; j++) ax[j] = 1.;
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    dims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    dims[0] = dims[1] = n;
    SET_SLOT(ans, Matrix_uploSym, ScalarString(mkChar("L")));
    SET_SLOT(ans, Matrix_diagSym,
	     (countDiag ? ScalarString(mkChar("N")) :
		 ScalarString(mkChar("U"))));
    parent_inv_ai(n, countDiag, pr, ai);
    UNPROTECT(1);
    return ans;
}
