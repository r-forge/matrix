#include "sscCrosstab.h"

SEXP sscCrosstab(SEXP flist, SEXP upper)
{
    int
	**fpt,
	*Ap,
	*Gp,
	*TTi,
	*Ti,
	*Tj,
	*dims,
	i,
	ncol = 0,
	nfac = length(flist),
	nfc2 = (nfac * (nfac - 1))/2, /* nfac choose 2 */
	nobs = length(VECTOR_ELT(flist, 0)),
	ntrpl,
	nz,
	pos,
	up = asLogical(upper);
    double
	*TTx,
	*Tx;
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

extern void ssclme_fill_LIp(int n, const int Parent[], int LIp[]);

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
void make_icounts(int nod, int ni, const int i[], const int j[],
		  const char ind[], int icounts[])
{
    int k, *lastj = memset(Calloc(ni, int), 0, sizeof(int) * ni);

    memset(icounts, 0, sizeof(int) * ni);
    for (k = 0; k < nod; k++) {
	int ik = i[k], jk = j[k];
	if (ind[k] && jk != lastj[ik]) {
	    icounts[ik]++;
	    lastj[ik] = jk;
	}
    }
    Free(lastj);
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
	nl1 = Gp[1],		/* number of levels of first factor */
	*icounts = Calloc(nl1, int),
	nl2, *jcounts,
	j, jj,
	n = length(pSlot) - 1,	/* number of columns */
	nf = length(GpSlot) - 1, /* number of factors */
	nz = length(iSlot),	/* number of non-zeros */
	nod = nz - n,		/* number of off-diagonals */
	nuse = nod,
	*iv = Calloc(nod, int),
	*jv = Calloc(nod, int),
	p1, p2;
    char *active = (char *) memset(Calloc(n, char), 1, (size_t) n),
	*ind = (char *) memset(Calloc(nod, char), 1, (size_t) nod);

    SEXP ans = PROTECT(allocVector(INTSXP, n));
    int *perm = INTEGER(ans);

    if (toupper(*CHAR(STRING_ELT(GET_SLOT(ctab, Matrix_uploSym), 0))) != 'L')
	error("Lower triangle required in sscCrosstab object");
    if (nf < 2) error("At least two factors are required");
    nl2 = Gp[2] - Gp[1];
    jcounts = Calloc(nl2, int);

    p1 = 0; p2 = nl1;		/* copy and expand off-diagonal indices */
    for (j = 0; j < p2; j++) {
	int p3 = Ap[j + 1];
	for (jj = Ap[j]; jj < p3; jj++) {
	    int i = Ai[jj];
/* FIXME: iv and jv are named backwards.  Reconcile this. */
	    if (i > j) {
		jv[p1] = i;
		iv[p1] = j;
		p1++;
	    }
	}
    }
    if (p1 != nod) error("Logic error - counts of off-diagonals disagree");
    
    p1 = 0; p2 = Gp[2] - 1;	/* fill 1st level from LHS, 2nd from RHS */
    while (nuse > 0) {
	int k, mcount, mind;
	make_icounts(nod, nl1, iv, jv, ind, icounts);

	mind = -1;
	mcount = Gp[2];	/* must be bigger than Gp[2] - Gp[1] */
	for (k = 0; k < nl1; k++) { /* determine column with min count */
	    int ic = icounts[k];
	    if (active[k]) {	
		if (ic < 1) {	/* remove active columns with zero counts */
		    perm[p1++] = k;
		    active[k] = 0;
		} else if (ic < mcount) {
		    mcount = ic;
		    mind = k;
		}
	    }
	}
	if (mind < 0)
	    error("Logic error - ran out of columns before nuse == 0");
	
				/* Count rows for j  */
	memset(jcounts, 0, sizeof(int) * nl2);
	for (k = 0; k < nod; k++) {
	    if (ind[k] && icounts[iv[k]] == mcount) {
		jcounts[jv[k] - nl1]++;
	    }
	}
	
	mcount = -1; mind = -1;	/* determine j with max count */
	for (k = 0; k < nl2; k++) {
	    int jc = jcounts[k];
				/* use >= below so ties give last row */
	    if (active[k + nl1] && jc >= mcount) {
		mcount = jc;
		mind = k;
	    }
	}
	if (mind < 0)
	    error("Logic error - ran out of rows before nuse == 0");

	perm[p2--] = mind + nl1;
	for (k = 0; k < nod; k++) {
	    if (jv[k] - nl1 == mind) {
		ind[k] = 0;
		nuse--;
	    }
	}
    }

    Free(ind); Free(active); Free(jv); Free(iv); Free(jcounts); Free(icounts);
    UNPROTECT(1);
    return ans;
}
