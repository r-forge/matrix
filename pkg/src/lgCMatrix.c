#include "lgCMatrix.h"

SEXP lgCMatrix_validate(SEXP x)
{
   SEXP pslot = GET_SLOT(x, Matrix_pSym),
      islot = GET_SLOT(x, Matrix_iSym);
    int j,
	ncol = length(pslot) - 1,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow,
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    nrow = dims[0];
    if (length(pslot) <= 0)
	return mkString(_("slot p must have length > 0"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slot i"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi)) {
      csc_sort_columns(ncol, xp, xi, (double *) NULL);
    }
    return ScalarLogical(1);
}

/** 
 * C := op(A) %*% op(B) + beta ^ C for logical sparse column-oriented matrices
 * 
 * @param tra nonzero if A is to be transposed
 * @param trb nonzero if B is to be transposed
 * @param m number of rows in C
 * @param n number of columns in C
 * @param k number of columns in A if tra == 0, otherwise number of
 *          rows in A
 * @param ai vector of row indices of TRUE elements in A
 * @param ap column pointers for A
 * @param bi vector of row indices of TRUE elements in B
 * @param bp column pointers for B
 * @param beta if non-zero existing TRUE elements in C are retained
 * @param ciP SEXP whose INTEGER part is the column indices of TRUE
 * elements in C (not used if beta == 0).
 * @param cp column pointers for C
 *
 * @return SEXP whose INTEGER part is the column indices of TRUE
 * elements in the product.  Note that the contents of cp may be modified.
 */
SEXP Matrix_lgClgCmm(int tra, int trb, int m, int n, int k,
		     const int ai[], const int ap[],
		     const int bi[], const int bp[],
		     int beta, SEXP CIP, int cp[])
{
/*     int annz = ap[tra ? m : k], bnnz = bp[trb ? k : n]; */
    int cnnz = cp[n], extra = 0;
    int *ci, i, j, prot = 0;	/* prot is the number of PROTECTs to UNPROTECT */
    
    if (beta) {
	ci = INTEGER(CIP);
    } else {			/* blank the C matrix */
	for (j = 0; j <= n; j++) cp[j] = 0;
	cnnz = 0;
	ci = (int *) NULL;
    } 
    if (tra) {
	if (trb) {		/* t(A) %*% t(B) */
	} else {		/* t(A) %*% B */
	}
    } else {
	if (trb) {		/* A %*% t(B) */
	} else {		/* A %*% B */
	    for (j = 0; j < n; j++) { /* col index for B and C */
		int ii, ii2 = bp[j + 1];
		for (ii = bp[j]; ii < ii2; ii++) { /* index into bi */
		    int jj = bi[ii]; /* row index of B; col index of A */
		    int i, i2 = ap[jj + 1]; /* index into ai */
		    for (i = ap[jj]; i < i2; i++)
			if (check_csc_index(cp, ci, ai[i], j, -1) < 0) extra++;
		}
	    }
	    if (extra) {
		int ntot = cnnz + extra;
		int *Ti = Calloc(ntot, int),
		    *rwInd = Calloc(m, int), /* indicator of TRUE in column j */
		    pos = 0;

		cp[0] = 0;
		for (j = 0; j < n; j++) {
		    int ii, ii2 = bp[j + 1];
		    cp[j + 1] = cp[j];
		    AZERO(rwInd, m); /* initialize column j of C to FALSE */
		    for (ii = bp[j]; ii < ii2; ii++) { /* index into bi */
			int jj = bi[ii]; /* row index of B; col index of A */
			int i, i2 = ap[jj + 1]; /* index into ai */
			for (i = ap[jj]; i < i2; i++) rwInd[ai[i]] = 1;
		    }
		    for (i = 0; i < m; i++)
			if (rwInd[i]) {cp[j + 1]++; Ti[pos++] = i;}
		}
		PROTECT(CIP = allocVector(INTSXP, cp[n])); prot++;
		Memcpy(INTEGER(CIP), Ti, cp[n]);
		Free(Ti); Free(rwInd);
	    }
	}
    }
    UNPROTECT(prot);
    return CIP;
}
	    
SEXP lgCMatrix_lgCMatrix_mm(SEXP a, SEXP b)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lgCMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    int k = adims[1], m = adims[0], n = bdims[1];
    int *cp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1));
    
    if (bdims[0] != k)
	error(_("Matrices are not conformable for multiplication"));
    cdims[0] = m; cdims[1] = n;
    SET_SLOT(ans, Matrix_iSym,
	     Matrix_lgClgCmm(0, 0, m, n, k,
			     INTEGER(GET_SLOT(a, Matrix_iSym)),
			     INTEGER(GET_SLOT(a, Matrix_pSym)),
			     INTEGER(GET_SLOT(b, Matrix_iSym)),
			     INTEGER(GET_SLOT(b, Matrix_pSym)),
			     0, (SEXP) NULL, cp));
    UNPROTECT(1);
    return ans;
}

SEXP lgCMatrix_trans(SEXP x)
{
    SEXP xi = GET_SLOT(x, Matrix_iSym);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lgCMatrix")));
    int *adims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	*xdims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nz = length(xi);
    int *xj = Calloc(nz, int);
    
    adims[1] = xdims[0];
    adims[0] = xdims[1];
    triplet_to_col(adims[0], adims[1], nz, 
		   expand_cmprPt(xdims[1], INTEGER(GET_SLOT(x, Matrix_pSym)), xj),
		   INTEGER(GET_SLOT(x, Matrix_iSym)), (double *) NULL,
		   INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP,  adims[1] + 1)),
		   INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nz)),
		   (double *) NULL);
    Free(xj);
    UNPROTECT(1);
    return ans;
}
