#include "sscCrosstab.h"

SEXP sscCrosstab(SEXP flist, SEXP upper)
{
    int **fpt, *Ap, *Gp, *TTi, *Ti, *Tj, *dims, i,
	ncol = 0,
	nfac = length(flist),
	nfc2 = (nfac * (nfac - 1))/2, /* nfac choose 2 */
	nobs = length(VECTOR_ELT(flist, 0)),
	ntrpl, nz, pos, up = asLogical(upper);
    double *TTx, *Tx;
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("sscCrosstab")));

    if (!isNewList(flist) || nfac < 1)
	error("flist must be a non-empty list");
    SET_SLOT(val, Matrix_GpSym, allocVector(INTSXP, nfac + 1));
    Gp = INTEGER(GET_SLOT(val, Matrix_GpSym));
    fpt = (int **) R_alloc(nfac, sizeof(int *));
    for (i = 0; i < nfac; i++) {
	SEXP el = VECTOR_ELT(flist, i);
	if (!inherits(el, "factor"))
	    error("flist must be a non-empty list of factors");
	if (length(el) != nobs)
	    error("All elements of flist must have the same length");
	Gp[i] = ncol;
	ncol += length(getAttrib(el, R_LevelsSymbol));
	fpt[i] = INTEGER(el);
    }
    Gp[nfac] = ncol;
    SET_SLOT(val, Matrix_uploSym, ScalarString(mkChar(up ? "U" : "L")));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    dims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    dims[0] = dims[1] = ncol;
    ntrpl = nfc2 * nobs + ncol;
    Ti = Calloc(ntrpl, int); Tj = Calloc(ntrpl, int); TTi = Calloc(ntrpl, int);
    Tx = Calloc(ntrpl, double); TTx = Calloc(ntrpl, double);
				/* Generate the triplet form of the result */
    for (i = 0; i < ncol; i++) {
	Ti[i] = Tj[i] = i;	/* The diagonals - these will store counts */
	Tx[i] = 0.0;
    }
    pos = ncol;
    for (i = 0; i < nobs ; i++) {
	int j, jcol, k;
	for (j = 0; j < nfac; j++) {
	    jcol = Gp[j] + fpt[j][i] - 1;
	    Tx[jcol] += 1.;	/* increment diagonal count */
	    for (k = j + 1; k < nfac; k++) { /* off-diagonals */
		int irow = Gp[k] + fpt[k][i] - 1;
		if (up) {
		    Ti[pos] = jcol; Tj[pos] = irow;
		} else {
		    Tj[pos] = jcol; Ti[pos] = irow;
		}
		Tx[pos] = 1.;
		pos++;
	    }
	}
    }
    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    Ap = INTEGER(GET_SLOT(val, Matrix_pSym));
    triplet_to_col(ncol, ncol, ntrpl, Ti, Tj, Tx, Ap, TTi, TTx);
    nz = Ap[ncol];		/* non-zeros in Z'Z crosstab */
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nz));
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), TTi, nz);
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), TTx, nz);
    Free(Ti); Free(Tj); Free(Tx);
    Free(TTi); Free(TTx);
    
    UNPROTECT(1);
    return val;
}

SEXP sscCrosstab_L_LI_sizes(SEXP ctab, SEXP permexp)
{
    SEXP ans = PROTECT(allocVector(INTSXP, 4));
    int *Ai = INTEGER(GET_SLOT(ctab, Matrix_iSym)),
	*Ap = INTEGER(GET_SLOT(ctab, Matrix_pSym)),
	*aa = INTEGER(ans),
	*perm = INTEGER(permexp),
	n = INTEGER(GET_SLOT(ctab, Matrix_DimSym))[1],
	*Lp = Calloc(n + 1, int),
	*Parent = Calloc(n, int),
	*Lnz = Calloc(n, int),
	*Flag = Calloc(n, int);

    ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag,
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    aa[0] = Lp[n];
    ssclme_fill_LIp(n, Parent, Lp);
    aa[1] = Lp[n];
    ssc_symbolic_permute(n, 1, perm, Ap, Ai);
    ldl_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag,
		 (int *) NULL, (int *) NULL); /* P & Pinv */
    aa[2] = Lp[n];
    ssclme_fill_LIp(n, Parent, Lp);
    aa[3] = Lp[n];
    Free(Flag); Free(Lnz); Free(Parent); Free(Lp); 
    UNPROTECT(1);
    return ans;
}

static
void col_metis_order(int j0, int j1, int i2,
		     const int Tp[], const int Ti[], int ans[])
{
    int j, nz = 0;		/* count off-diagonal pairs */
    for (j = j0; j < j1; j++) {	/* columns of interest */
	int ii, nr = 0, p2 = Tp[j + 1];
	for (ii = Tp[j]; ii < p2; ii++) {
	    int i = Ti[ii];
	    if (j1 <= i && i < i2) nr++; /* verify row index */
	}
	nz += (nr * (nr - 1))/2; /* add number of pairs of rows */
    }
    if (nz > 0) {		/* Form an ssc Matrix */
	int j, n = i2 - j1,	/* number of rows */
	    nnz = n + nz, pos;
	int *Ap = Calloc(n + 1, int),
	    *Ai = Calloc(nnz, int),
	    *Tj = Calloc(nnz, int),
	    *TTi = Calloc(nnz, int),
	    *perm = Calloc(n, int),
	    *iperm = Calloc(n, int);

	for (j = 0; j < n; j++) { /* diagonals */
	    TTi[j] = Tj[j] = j;
	    Tx[j] = 1.;
	}
	pos = n;
	for (j = j0; j < j1; j++) { /* create the pairs */
	    int ii, nr = 0, p2 = Tp[j + 1];
	    for (ii = Tp[j]; ii < p2; ii++) {
		int r1 = Ti[ii], i1;
		if (j1 <= r1 && r1 < i2) {
		    for (i1 = ii + 1; i1 < p2; i1++) {
			int r2 = Ti[i1];
			if (r2 < i2) {
			    TTi[pos] = r2 - j1;
			    Tj[pos] = r1 - j1;
			    Tx[pos] = 1.;
			    pos++;
			}
		    }
		}
	    }
	}
	triplet_to_col(n, n, nnz, TTi, Tj, (double *) 0, Ap, Ai, (double *) 0);
	ssc_metis_order(n, Ap, Ai, perm, iperm);
	for (j = j1; j < i2; j++) ans[j] = j1 + iperm[j - j1];
	Free(TTi); Free(Tj); Free(Ai); Free(Ap);
	Free(perm); Free(iperm);
    }
}

SEXP sscCrosstab_groupedPerm(SEXP ctab)
{
    SEXP
	GpSlot = GET_SLOT(ctab, Matrix_GpSym),
	iSlot = GET_SLOT(ctab, Matrix_iSym),
	pSlot = GET_SLOT(ctab, Matrix_pSym);
    int *Ai = INTEGER(iSlot),
	*Ap = INTEGER(pSlot),
	*Gp = INTEGER(GpSlot),
	i,
	n = length(pSlot) - 1,	/* number of columns */
	nf = length(GpSlot) - 1, /* number of factors */
	up;
    SEXP ans = PROTECT(allocVector(INTSXP, n));

    up = toupper(*CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0))) != 'L';
    if (nf > 1 && up) {			/* transpose */
	int nz = length(iSlot);
	int *ai = Calloc(nz, int),
	    *ap = Calloc(n + 1, int);
	double *ax = Calloc(nz, double);

	csc_components_transpose(n, n, nz, Ap, Ai,
				 REAL(GET_SLOT(ctab, Matrix_xSym)),
				 ap, ai, ax);
	Ap = ap;
	Ai = ai;
	Free(ax);		/* don't need values, only positions */
    }
    for (i = 0; i < n; i++) {
	INTEGER(ans)[i] = i;    /* initialize permutation to identity */
    }
    for (i = 1; i < nf; i++) {
	col_metis_order(Gp[i - 1], Gp[i], Gp[i+1], Ap, Ai, INTEGER(ans));
    }
    if (nf > 1 && up) {Free(Ap); Free(Ai);}
    UNPROTECT(1);
    return ans;
}
