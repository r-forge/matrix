#include "bCrosstab.h"

/** 
 * Calculate the zero-based index in a row-wise packed lower
 * triangular matrix.  This is used for the arrays of blocked sparse matrices.
 * 
 * @param i row number (0-based)
 * @param k column number (0-based)
 * 
 * @return The 0-based index of the (i,k) element of a row-wise packed lower
 * triangular matrix.
 */    
static R_INLINE
int Lind(int i, int k)
{
    return (i * (i + 1))/2 + k;
}

/** 
 * Permute an index vector
 * 
 * @param i vector of 0-based indices
 * @param nnz length of vector i
 * @param perm 0-based permutation vector of length max(i)
 */
static R_INLINE
void ind_permute(int i[], int nnz, const int perm[])
{
    int j;
    for (j = 0; j < nnz; j++) i[j] = perm[i[j]];
}

/** 
 * Force indices to be in the upper triangle of a matrix
 * 
 * @param i vector of 0-based row indices
 * @param j vector of 0-based column indices
 * @param nnz length of index vectors
 */
static R_INLINE
void make_upper_triangular(int i[], int j[], int nnz)
{
    int k;
    for (k = 0; k < nnz; k++) {
	if (i[k] > j[k]) {
	    int tmp = i[k];
	    i[k] = j[k];
	    j[k] = tmp;
	}
    }
}

/** 
 * Create a named list of length n
 * 
 * @param n length of list to return
 * @param names names of list elements
 * 
 * @return pointer to a named list of length n
 */
static
SEXP make_named_list(int n, char **names)
{
    SEXP ans = PROTECT(allocVector(VECSXP, n)),
	nms = PROTECT(allocVector(STRSXP, n));
    int i;

    for (i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

static R_INLINE
int check_csc_index(const int p[], const int i[], int row, int col)
{
    int k, k2 = p[col + 1];
				/* use a linear search for now */
				/* perhaps replace by bsearch later */
    for (k = p[col]; k < k2; k++) {
	if (i[k] == row) return k;
    }
    return -1;
}
    
/** 
 * Update a diagonal block
 * 
 * @param ctab pointer to a blocked crosstabulation object
 * @param j index of updating column
 * @param i index of diagonal block to be updated
 */
static void diag_update(SEXP ctab, int j, int i)
{
    SEXP db = VECTOR_ELT(ctab, Lind(i, i)),
	jb = VECTOR_ELT(ctab, Lind(i, j));
    SEXP dpp = GET_SLOT(db, Matrix_pSym),
	jpp = GET_SLOT(jb, Matrix_pSym);
    int *di = INTEGER(GET_SLOT(db, Matrix_iSym)),
	*dp = INTEGER(dpp),
	*ji = INTEGER(GET_SLOT(jb, Matrix_iSym)),
	*jp = INTEGER(jpp),
	dnc = length(dpp) - 1,
	jnc = length(jpp) - 1;
    int jj, extra;
				/* bound the number of extra elements */
    extra = 0;
    for (jj = 0; jj < jnc; jj++) {
	int k, kk, k2 = jp[jj + 1];
	for (k = jp[jj]; k < k2; k++) {
	    for (kk = k; kk < k2; kk++) {
		if (check_csc_index(dp, di, ji[k], ji[kk]) < 0) extra++;
	    }
	}
    }
    if (!extra) return;
    {
	int pos, nnz = dp[dnc];
	int ntot = nnz + extra;
	int *Ai = Calloc(ntot, int),
	    *Ti = Calloc(ntot, int),
	    *Tj = Calloc(ntot, int);
	double *Ax;

	Memcpy(Ti, di, nnz);	/* make a copy of the row indices */
	pos = 0;		/* fill in the column indices */
	for (jj = 0; jj < dnc; jj++) {
	    int j2 = dp[jj + 1];
	    while (pos < j2) {
		Tj[pos] = jj;
		pos++;
	    }
	}
				/* add the extra elements */
	for (jj = 0; jj < jnc; jj++) {
	    int k, kk, k2 = jp[jj + 1];
	    for (k = jp[jj]; k < k2; k++) {
		for (kk = k; kk < k2; kk++) {
		    if (check_csc_index(dp, di, ji[k], ji[kk]) < 0) {
			Ti[pos] = ji[k];
			Tj[pos] = ji[kk];
			pos++;
		    }
		}
	    }
	}
	triplet_to_col(dnc, dnc, ntot, Ti, Tj, (double *) NULL,
		       dp, Ai, (double *) NULL);
	nnz = dp[dnc];
	SET_SLOT(db, Matrix_iSym, allocVector(INTSXP, nnz));
	Memcpy(INTEGER(GET_SLOT(db, Matrix_iSym)), Ai, nnz);
	SET_SLOT(db, Matrix_xSym, allocVector(REALSXP, nnz));
	Ax = REAL(GET_SLOT(db, Matrix_xSym));
	for (j = 0; j < nnz; j++) Ax[j] = 1.;
	Free(Ai); Free(Ti); Free(Tj);
	return;
    }
}

/** 
 * Update a block
 * 
 * @param ctab pointer to a blocked crosstabulation object
 * @param j index of updating column
 * @param k column index of block to be updated 
 * @param i row index of block to be updated (j < k <= i)
 */
static void block_update(SEXP ctab, int j, int k, int i)
{
    SEXP tb = VECTOR_ELT(ctab, Lind(i, k)),
	ib = VECTOR_ELT(ctab, Lind(i, j)),
	kb = VECTOR_ELT(ctab, Lind(k, j));
    SEXP tpp = GET_SLOT(tb, Matrix_pSym),
	kpp = GET_SLOT(kb, Matrix_pSym);
    int *ti = INTEGER(GET_SLOT(tb, Matrix_iSym)),
	*tp = INTEGER(tpp),
	*ii = INTEGER(GET_SLOT(ib, Matrix_iSym)),
	*ip = INTEGER(GET_SLOT(ib, Matrix_pSym)),
	*ki = INTEGER(GET_SLOT(kb, Matrix_iSym)),
	*kp = INTEGER(kpp),
	tnc = length(tpp) - 1,
	knc = length(kpp) - 1;
    int jj, extra;

    if (k > i || j >= k)
	error("i,j,k values of %d,%d,%d do not satisfy j < k <= i",
	      i, j, k);
				/* bound the number of extra elements */
    extra = 0;
    for (jj = 0; jj < knc; jj++) {
	int i1, kk, i2 = ip[jj + 1], k2 = kp[jj + 1];
	for (kk = kp[jj]; kk < k2; kk++) {
	    for (i1 = ip[jj]; i1 < i2; i1++) {
		    if ((check_csc_index(tp, ti, ii[i1], ki[kk]) < 0) &&
				/* only update upper triangle of
				 * diagonal blocks */
			((k != i) || (ii[i1] <= ki[kk]))) extra++;
	    }
	}
    }
    if (!extra) return;
    {
	int pos, nnz = tp[tnc];
	int ntot = nnz + extra;
	int *Ai = Calloc(ntot, int),
	    *Ti = Calloc(ntot, int),
	    *Tj = Calloc(ntot, int),
	    *Dims = INTEGER(GET_SLOT(tb, Matrix_DimSym));
	double *Ax;

	Memcpy(Ti, ti, nnz);	/* make a copy of the row indices */
	for (pos = 0, jj = 0; jj < tnc; jj++) {	/* fill in the column indices */
	    int j2 = tp[jj + 1];
	    for (; pos < j2; pos++) Tj[pos] = jj;
	}
				/* add the extra elements */
	for (jj = 0; jj < knc; jj++) {
	    int i1, kk, i2 = ip[jj + 1], k2 = kp[jj + 1];
	    for (kk = kp[jj]; kk < k2; kk++) {
		for (i1 = ip[jj]; i1 < i2; i1++) {
		    if ((check_csc_index(tp, ti, ii[i1], ki[kk]) < 0) &&
			((k != i) || (ii[i1] <= ki[kk]))) { 
			Ti[pos] = ii[i1];
			Tj[pos] = ki[kk];
			pos++;
		    }
		}
	    }
	}
	triplet_to_col(Dims[0], Dims[1], ntot, Ti, Tj, (double *) NULL,
		       tp, Ai, (double *) NULL);
	nnz = tp[tnc];
	SET_SLOT(tb, Matrix_iSym, allocVector(INTSXP, nnz));
	Memcpy(INTEGER(GET_SLOT(tb, Matrix_iSym)), Ai, nnz);
	SET_SLOT(tb, Matrix_xSym, allocVector(REALSXP, nnz));
	Ax = REAL(GET_SLOT(tb, Matrix_xSym));
	for (j = 0; j < nnz; j++) Ax[j] = 1.;
	Free(Ai); Free(Ti); Free(Tj);
	return;
    }
}

/** 
 * Permute the levels of one of the grouping factors in a bCrosstab object
 * 
 * @param ctab Pointer to a bCrosstab object
 * @param nf number of factors in ctab
 * @param jj index (0-based) of the factor levels to permute
 * @param ncj number of columns in level jj
 * @param perm permutation (0-based) to apply
 * @param pperm inverse of the permutation
 */
void bCrosstab_permute(SEXP ctab, int nf, int jj,
		       const int perm[], const int iperm[])
{
    SEXP trip, ipt, jpt;
    int j, nnz;
    
    for (j = 0; j < nf; j++) {
	int ind = (j < jj ? Lind(jj, j) : Lind(j, jj));
	SEXP csc = VECTOR_ELT(ctab, ind);
	int *Dims = INTEGER(GET_SLOT(csc, Matrix_DimSym));

	trip = csc_to_triplet(csc);
	ipt = GET_SLOT(trip, Matrix_iSym);
	nnz = length(ipt);
	jpt = GET_SLOT(trip, Matrix_jSym);
	if (j <= jj) ind_permute(INTEGER(ipt), nnz, iperm);
	if (j >= jj) ind_permute(INTEGER(jpt), nnz, iperm);
	if (j == jj)
	    make_upper_triangular(INTEGER(ipt), INTEGER(jpt), nnz);
	triplet_to_col(Dims[0], Dims[1], nnz, INTEGER(ipt),
		       INTEGER(jpt),
		       REAL(GET_SLOT(trip, Matrix_xSym)),
		       INTEGER(GET_SLOT(csc, Matrix_pSym)),
		       INTEGER(GET_SLOT(csc, Matrix_iSym)),
		       REAL(GET_SLOT(csc, Matrix_xSym)));
    }
}

static void
factor_levels_permute(SEXP dest, SEXP src, const int perm[],
		      const int iperm[])
{
    SEXP dlev = getAttrib(dest, R_LevelsSymbol),
	slev = getAttrib(src, R_LevelsSymbol);
    int nlev = length(dlev), flen = length(dest);
    int *d = INTEGER(dest), *s = INTEGER(src), i;

    if (length(slev) != nlev)
	error("number of levels in src and dest must match");
    if (length(src) != flen)
	error("length of src and dest must match");
    for (i = 0; i < nlev; i++)
	SET_STRING_ELT(dlev, i, STRING_ELT(slev, perm[i]));
    for (i = 0; i < flen; i++)
	d[i] = 1 + iperm[s[i]-1];
}

SEXP bCrosstab_convert(SEXP bCtab)
{
    char *anms[] = {"flist", "L", "Linv", "perm", "Parent"};
    SEXP flist = VECTOR_ELT(bCtab, 0),
	ctab = PROTECT(duplicate(VECTOR_ELT(bCtab, 1))),
	ans = PROTECT(make_named_list(5, anms));
    SEXP fcp, perm, L, Linv, Parent;
    int ctbl = length(ctab), j, nf = length(flist);

    if (ctbl != (nf*(nf + 1))/2)
	error("length of ctab = %d is not permisable", ctbl);
    SET_VECTOR_ELT(ans, 0, duplicate(flist));
    SET_VECTOR_ELT(ans, 1, allocVector(VECSXP, ctbl));
    SET_VECTOR_ELT(ans, 2, allocVector(VECSXP, nf));
    SET_VECTOR_ELT(ans, 3, allocVector(VECSXP, nf));
    SET_VECTOR_ELT(ans, 4, allocVector(VECSXP, nf));
    fcp = VECTOR_ELT(ans, 0);
    L = VECTOR_ELT(ans, 1);
    setAttrib(L, R_NamesSymbol, duplicate(getAttrib(ctab, R_NamesSymbol)));
    Linv = VECTOR_ELT(ans, 2);
    setAttrib(Linv, R_NamesSymbol, duplicate(getAttrib(flist, R_NamesSymbol)));
    perm = VECTOR_ELT(ans, 3);
    setAttrib(perm, R_NamesSymbol, duplicate(getAttrib(flist, R_NamesSymbol)));    
    Parent = VECTOR_ELT(ans, 4);
    setAttrib(Parent, R_NamesSymbol, duplicate(getAttrib(flist, R_NamesSymbol)));
    for (j = 0; j < nf; j++) {
	int dind = Lind(j, j), i, k;
	SEXP ctd = VECTOR_ELT(ctab, dind); /* diagonal in crosstab */
	SEXP Dimslot = GET_SLOT(ctd, Matrix_DimSym),
	    Linvj, Ljj, cpp = GET_SLOT(ctd, Matrix_pSym),
	    cip = GET_SLOT(ctd, Matrix_iSym);
	int *Lp, *Linvp, *Perm, *cp = INTEGER(cpp),
	    *ci = INTEGER(cip), *parent,
	    ncj = length(cpp) - 1,
	    nnz = length(cip);
	double *dtmp;
				
	SET_VECTOR_ELT(Parent, j, allocVector(INTSXP, ncj));
	parent = INTEGER(VECTOR_ELT(Parent, j));
	SET_VECTOR_ELT(perm, j, allocVector(INTSXP, ncj));
	Perm = INTEGER(VECTOR_ELT(perm, j));
	SET_VECTOR_ELT(L, dind, NEW_OBJECT(MAKE_CLASS("cscMatrix")));
	Ljj = VECTOR_ELT(L, dind);
	SET_VECTOR_ELT(Linv, j, NEW_OBJECT(MAKE_CLASS("cscMatrix")));
	Linvj = VECTOR_ELT(Linv, j);
	SET_SLOT(Ljj, Matrix_DimSym, duplicate(Dimslot));
	SET_SLOT(Linvj, Matrix_DimSym, duplicate(Dimslot));
	SET_SLOT(Ljj, Matrix_factorization, allocVector(VECSXP, 0));
	SET_SLOT(Linvj, Matrix_factorization, allocVector(VECSXP, 0));
	SET_SLOT(Ljj, Matrix_pSym, allocVector(INTSXP, ncj + 1));
	SET_SLOT(Linvj, Matrix_pSym, allocVector(INTSXP, ncj + 1));
	Lp = INTEGER(GET_SLOT(Ljj, Matrix_pSym));
	Linvp = INTEGER(GET_SLOT(Linvj, Matrix_pSym));
	if (nnz > ncj) {	/* calculate fill-reducing permutation */
	    int *Li, *Lnz = Calloc(ncj, int), *tmp = Calloc(ncj, int), info;

	    ssc_metis_order(ncj, cp, ci, Perm, tmp);
				/* apply to the crosstabulation */
	    bCrosstab_permute(ctab, nf, j, Perm, tmp);
				/* apply to the factor */
	    factor_levels_permute(VECTOR_ELT(fcp, j),
				  VECTOR_ELT(flist, j), Perm, tmp);
				/* symbolic analysis to get Parent */
	    ldl_symbolic(ncj, cp, ci, Lp, parent, Lnz, tmp,
			 (int *) NULL, (int *) NULL);
				/* decompose the identity to get the row pointers */
	    dtmp = REAL(GET_SLOT(ctd, Matrix_xSym));
	    for (i = 0; i < nnz; i++) dtmp[i] = 0.; /* initialize */
				/* diagonal el is last element in the column */
	    for (i = 0; i < ncj; i++) dtmp[cp[i+1] - 1] = 1.;
	    nnz = Lp[ncj];
	    SET_SLOT(Ljj, Matrix_iSym, allocVector(INTSXP, nnz));
	    SET_SLOT(Ljj, Matrix_xSym, allocVector(REALSXP, nnz));
	    info = ldl_numeric(ncj, cp, ci, dtmp, Lp, parent, Lnz,
			       INTEGER(GET_SLOT(Ljj, Matrix_iSym)),
			       REAL(GET_SLOT(Ljj, Matrix_xSym)),
			       (double *) R_alloc(ncj, sizeof(double)),	/* D */
			       (double *) R_alloc(ncj, sizeof(double)),	/* Y */
			       (int *) R_alloc(ncj, sizeof(int)),	/* Pattern */
			       (int *) R_alloc(ncj, sizeof(int)),	/* Flag */
			       (int *) NULL, (int *) NULL);
	    if (info < ncj) error("identity matrix is not positive definite?");
	    Free(Lnz); Free(tmp);
	} else {
	    for (i = 0; i <= ncj; i++) {
		Lp[i] = 0;
		parent[i] = -1;
		Perm[i] = i;
	    }
	    SET_SLOT(Ljj, Matrix_iSym, allocVector(INTSXP, 0));
	    SET_SLOT(Ljj, Matrix_xSym, allocVector(REALSXP, 0));
	}
				/* Derive the diagonal block of Linv */
	nnz = parent_inv_ap(ncj, 0, parent, Linvp);
	SET_SLOT(Linvj, Matrix_iSym, allocVector(INTSXP, nnz));
	parent_inv_ai(ncj, 0, parent, INTEGER(GET_SLOT(Linvj, Matrix_iSym)));
	SET_SLOT(Linvj, Matrix_xSym, allocVector(REALSXP, nnz));
	dtmp = REAL(GET_SLOT(Linvj, Matrix_xSym));
	for (i = 0; i < nnz; i++) dtmp[i] = 1.;
/* FIXME: Update any blocks below the diagonal in this column if necessary*/
	for (k = j+1; k < nf; k++) { /* update remaining columns, if any */
	    for (i = k; i < nf; i++) block_update(ctab, j, k, i);
	}
	for (i= 0; i < j; i++) { /* copy blocks to the left */
	    SET_VECTOR_ELT(L, Lind(j,i), duplicate(VECTOR_ELT(ctab, Lind(j,i))));
	}
    }
    UNPROTECT(2);
    return ans;
}



	
