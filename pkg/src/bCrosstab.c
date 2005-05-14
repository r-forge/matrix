#include "bCrosstab.h"
/* TODO:
 * - Only do a fill-reducing permutation on the first non-nested factor
 * - Alternatively: change the algorithm for the fill-reducing
 *   permutation to a greedy or a picky algorithm.
 */


/** 
 * Replace the structure of C by the structure of CA^{-T}
 * 
 * @param anc number of column blocks in A
 * @param Parent parent array for column blocks of A
 * @param C a dgBCMatrix object to be updated
 */
static void
symbolic_right_unit_sm_trans(int anc, const int Parent[], SEXP C)
{
    SEXP cip = GET_SLOT(C, Matrix_iSym),
	cpp = GET_SLOT(C, Matrix_pSym);
    int *ci = INTEGER(cip),
	*cp = INTEGER(cpp),
	cnz = length(cip),
	i, j, nextra = 0;

    if ((length(cpp) - 1) != anc)
	error(_("No. of cols in A (%d) does not match no. of cols in C (%d)"),
	      anc, length(cpp) - 1);
    
    i = 1;			/* check for A being the identity */
    for (j = 0; j < anc; j++) {
	if (Parent[j] >= 0) {
	    i = 0;
	    break;
	}
    }
    if (i) return;		/* A is the identity */

    for (j = 0; j < anc; j++) { /* bound the number of extra triplets */
	int cj2 = cp[j + 1], ka, kc;
	for (ka = Parent[j]; ka >= 0; ka = Parent[ka]) {
	    for (kc = cp[j]; kc < cj2; kc++) {
		if (check_csc_index(cp, ci, ci[kc], ka, -1) < 0) nextra++;
	    }
	}
    }
    if (nextra) {
	int cnr, ntot = cnz + nextra, pos;
	int
	    *dims = INTEGER(getAttrib(GET_SLOT(C, Matrix_xSym), R_DimSymbol)),
	    *Ti = Memcpy((int *) Calloc(ntot, int), ci, cnz),
	    *Tj = expand_cmprPt(anc, cp, Calloc(ntot, int)),
	    *Ci = Calloc(ntot, int);

	for (j = 0, pos = cnz; j < anc; j++) {
	    int cj2 = cp[j + 1], ka, kc;
	    for (ka = Parent[j]; ka >= 0; ka = Parent[ka]) {
		for (kc = cp[j]; kc < cj2; kc++) {
		    if (check_csc_index(cp, ci, ci[kc], ka, -1) < 0) {
			Tj[pos] = ka;
			Ti[pos] = ci[kc];
			pos++;
		    }
		}
	    }
	}
	for (j = 0, cnr = 0; j < cnz; j++) { /* determine number of rows in C */
	    int rr = ci[j] + 1;
	    if (rr > cnr) cnr = rr;
	}
	triplet_to_col(cnr, anc, ntot, Ti, Tj, (double *) NULL,
		       INTEGER(cpp), Ci, (double *) NULL);
	cnz = cp[anc];
	SET_SLOT(C, Matrix_iSym, allocVector(INTSXP, cnz));
	SET_SLOT(C, Matrix_xSym, alloc3Darray(REALSXP, dims[0], dims[1], cnz));
	Free(Ti); Free(Tj); Free(Ci);
    }
}
    
/** 
 * Update a diagonal block of ZZpO in the blocked crosstabulation
 * 
 * @param db pointer to the diagonal block
 * @param odb pointer to the off-diagonal block
 */
static void
diag_update(SEXP db, SEXP odb)
{
    SEXP dbpP = GET_SLOT(db, Matrix_pSym),
	odbpP = GET_SLOT(odb, Matrix_pSym);

    SET_SLOT(db, Matrix_iSym,
	     Matrix_lgCsyrk(1, 0, length(dbpP) - 1, length(odbpP) - 1,
			    INTEGER(GET_SLOT(odb, Matrix_iSym)),
			    INTEGER(odbpP), 1,
			    GET_SLOT(db, Matrix_iSym),
			    INTEGER(dbpP)));
}

/** 
 * Update an off-diagonal block of L from the blocked crosstabulation
 * 
 * @param A upper block
 * @param B lower block
 * @param C product block
 * @param nrA number of rows in A
 */
static void
offdiag_update(SEXP A, SEXP B, SEXP C, int nrA)
{
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	CpP = GET_SLOT(C, Matrix_pSym);

    SET_SLOT(C, Matrix_iSym,
	     Matrix_lgClgCmm(0, 1, nrA, length(CpP) - 1,
			     length(ApP) - 1,
			     INTEGER(GET_SLOT(A, Matrix_iSym)),
			     INTEGER(ApP),
			     INTEGER(GET_SLOT(B, Matrix_iSym)),
			     INTEGER(GET_SLOT(B, Matrix_pSym)),
			     1,
			     GET_SLOT(C, Matrix_iSym),
			     INTEGER(CpP)));
}

	/* check dimensions of x slots and modify if necessary */
/** 
 * Check the dimensions of the x slot and modify if necessary
 * 
 * @param mm Pointer to a dgBCMatrix object
 */
static void
check_x_slot_dim(SEXP mm)
{
    SEXP xP = GET_SLOT(mm, Matrix_xSym);
    int *dims = INTEGER(getAttrib(xP, R_DimSymbol)),
	nnz = length(GET_SLOT(mm, Matrix_iSym));

    if (dims[2] != nnz)
	SET_SLOT(mm, Matrix_xSym,
		 alloc3Darray(REALSXP, dims[0], dims[1], nnz));
}

/** 
 * Permute the levels of one of the grouping factors in a bCrosstab object
 * 
 * @param ctab Pointer to a bCrosstab object
 * @param nf number of factors in ctab
 * @param jj index (0-based) of the factor levels to permute
 * @param nlev number of levels of the grouping factors
 * @param iperm inverse of the permutation
 */
static void
bCrosstab_permute(SEXP ctab, int nf, int jj,
		  const int nlev[], const int iperm[])
{
    int j;
    for (j = 0; j < nf; j++) {
	int ind = (j < jj ? Lind(jj, j) : Lind(j, jj)),
	    ncol = (j < jj ? nlev[j] : nlev[jj]),
	    nrow = (j < jj ? nlev[jj] : nlev[j]);
	SEXP cscb = VECTOR_ELT(ctab, ind),
	    cscbi = GET_SLOT(cscb, Matrix_iSym);
	int *cp = INTEGER(GET_SLOT(cscb, Matrix_pSym)),
	    nnz = length(cscbi);
/* 	double *cx = REAL(GET_SLOT(cscb, Matrix_xSym)); */
	int *mj = expand_cmprPt(ncol, cp, Calloc(nnz, int));
	int *mi = Memcpy(Calloc(nnz, int), INTEGER(cscbi), nnz);
/* 	double *mx = Memcpy(Calloc(nnz, double), cx, nnz); */

	if (j <= jj) int_permute(mi, nnz, iperm);
	if (j >= jj) int_permute(mj, nnz, iperm);
	if (j == jj) make_upper_triangular(mi, mj, nnz);
/* 	triplet_to_col(nrow, ncol, nnz, mi, mj, mx, cp, INTEGER(cscbi), cx); */
	triplet_to_col(nrow, ncol, nnz, mi, mj, (double *) NULL,
		       cp, INTEGER(cscbi), (double *) NULL);
	Free(mi); Free(mj);
/* 	Free(mx); */
    }
}

static void
symmetric_permute(SEXP A, int nlev, const int iperm[])
{
    SEXP AiP = GET_SLOT(A, Matrix_iSym);
    int *Ap = INTEGER(GET_SLOT(A, Matrix_pSym)),
	nnz = length(AiP);
/*     double *Ax = REAL(GET_SLOT(A, Matrix_xSym)); */
    int *mj = expand_cmprPt(nlev, Ap, Calloc(nnz, int));
    int *mi = Memcpy(Calloc(nnz, int), INTEGER(AiP), nnz);
/*     double *mx = Memcpy(Calloc(nnz, double), Ax, nnz); */

    int_permute(mi, nnz, iperm);
    int_permute(mj, nnz, iperm);
    make_upper_triangular(mi, mj, nnz);
/*     triplet_to_col(nlev, nlev, nnz, mi, mj, mx, Ap, INTEGER(AiP), Ax); */
    triplet_to_col(nlev, nlev, nnz, mi, mj, (double *) NULL,
		   Ap, INTEGER(AiP), (double *) NULL);
    Free(mi); Free(mj);
/*     Free(mx); */
}

/** 
 * Apply a permutation vector to the levels of a factor.
 *
 * The dest pointer is assumed to point to a copy of the src pointer's
 * contents.
 * 
 * @param dest pointer to the destination factor
 * @param src pointer to the source factor
 * @param perm permutation vector (0-based)
 * @param iperm inverse permutation vector (0-based)
 */
static void
factor_levels_permute(SEXP dest, SEXP src, const int perm[],
		      const int iperm[])
{
    SEXP dlev = getAttrib(dest, R_LevelsSymbol),
	slev = getAttrib(src, R_LevelsSymbol);
    int nlev = length(dlev), flen = length(dest);
    int *d = INTEGER(dest), *s = INTEGER(src), i;

    if (length(slev) != nlev)
	error(_("number of levels in src and dest must match"));
    if (length(src) != flen)
	error(_("length of src and dest must match"));
    for (i = 0; i < nlev; i++)
	SET_STRING_ELT(dlev, i, STRING_ELT(slev, perm[i]));
    for (i = 0; i < flen; i++)
	d[i] = 1 + iperm[s[i]-1];
}

/** 
 * Create and populate slots in an lmer object from the blocked crosstabulation.
 * 
 * @param val Pointer to an lmer object
 */
void
lmer_populate(SEXP val)
{
    SEXP D, L, Parent, ZZpO, flist = GET_SLOT(val, Matrix_flistSym),
	perm, Omega, ZtZ = GET_SLOT(val, Matrix_ZtZSym);
    SEXP fnms = getAttrib(flist, R_NamesSymbol);
    int j, k, nf = length(flist);
    int *nc = INTEGER(GET_SLOT(val, Matrix_ncSym)), *Gp,
	*nlev = Calloc(nf, int), npairs = (nf * (nf + 1))/2;
    char *statnms[] = {"factored", "inverted", ""},
	*devnms[] = {"ML", "REML", ""},
	*pnms[] = {"index", "block", ""};
	
    /* Allocate fixed-sized slots */
    SET_SLOT(val, Matrix_statusSym, Matrix_make_named(LGLSXP, statnms));
    SET_SLOT(val, Matrix_devianceSym, Matrix_make_named(REALSXP, devnms));
    SET_SLOT(val, Matrix_devCompSym, allocVector(REALSXP, 4));
    /* Allocate slots that are lists of length nf */
    ZZpO = ALLOC_SLOT(val, Matrix_ZZpOSym, VECSXP, nf);
    setAttrib(ZZpO, R_NamesSymbol, duplicate(fnms));
    D = ALLOC_SLOT(val, Matrix_DSym, VECSXP, nf);
    setAttrib(D, R_NamesSymbol, duplicate(fnms));
    perm = ALLOC_SLOT(val, Matrix_permSym, VECSXP, nf);
    setAttrib(perm, R_NamesSymbol, duplicate(fnms));    
    Parent = ALLOC_SLOT(val, Matrix_ParentSym, VECSXP, nf);
    setAttrib(Parent, R_NamesSymbol, duplicate(fnms));
    Omega = ALLOC_SLOT(val, Matrix_OmegaSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    
    /* Allocate peculiar length slots */
    L = ALLOC_SLOT(val, Matrix_LSym, VECSXP, npairs);
    Gp = INTEGER(ALLOC_SLOT(val, Matrix_GpSym, INTSXP, nf + 1));
    Gp[0] = 0;
    for (j = 0; j < nf; j++) {
	nlev[j] = length(getAttrib(VECTOR_ELT(flist, j), R_LevelsSymbol));
	Gp[j + 1] = Gp[j] + nc[j] * nlev[j];
	SET_VECTOR_ELT(D, j, alloc3Darray(REALSXP, nc[j], nc[j], nlev[j]));
	SET_VECTOR_ELT(Omega, j, allocMatrix(REALSXP, nc[j], nc[j]));
	SET_VECTOR_ELT(ZZpO, j, duplicate(VECTOR_ELT(ZtZ, Lind(j, j))));
	for (k = j; k < nf; k++)
	    SET_VECTOR_ELT(L, Lind(k, j),
			   duplicate(VECTOR_ELT(ZtZ, Lind(k, j))));
    }
    SET_SLOT(val, Matrix_XtXSym, allocMatrix(REALSXP, nc[nf], nc[nf]));
    AZERO(REAL(GET_SLOT(val, Matrix_XtXSym)), nc[nf] * nc[nf]);
    SET_SLOT(val, Matrix_RXXSym, allocMatrix(REALSXP, nc[nf], nc[nf]));
    AZERO(REAL(GET_SLOT(val, Matrix_RXXSym)), nc[nf] * nc[nf]);
    SET_SLOT(val, Matrix_ZtXSym, allocMatrix(REALSXP, Gp[nf], nc[nf]));
    SET_SLOT(val, Matrix_RZXSym, allocMatrix(REALSXP, Gp[nf], nc[nf]));
    for (j = 0; j < nf; j++) {
	int dind = Lind(j, j), i;
	SEXP ctd = VECTOR_ELT(ZZpO, j); /* diagonal in crosstab */
	SEXP Ljj = VECTOR_ELT(L, dind),
	    cpp = GET_SLOT(ctd, Matrix_pSym),
	    cip = GET_SLOT(ctd, Matrix_iSym), parent;
	int *Lp = INTEGER(GET_SLOT(Ljj, Matrix_pSym)), *Perm,
	    *cp = INTEGER(cpp),
	    *ci = INTEGER(cip),
	    ncj = length(cpp) - 1,
	    nnz = length(cip);
				
	SET_VECTOR_ELT(Parent, j, Matrix_make_named(VECSXP, pnms));
	parent = VECTOR_ELT(Parent, j);
	SET_VECTOR_ELT(parent, 0, allocVector(INTSXP, ncj));
	SET_VECTOR_ELT(parent, 1, allocVector(INTSXP, ncj));
	SET_VECTOR_ELT(perm, j, allocVector(INTSXP, ncj));
	Perm = INTEGER(VECTOR_ELT(perm, j));
	if (nnz > ncj) {	/* calculate fill-reducing permutation */
	    SEXP fac = VECTOR_ELT(flist, j);
	    SEXP fcp = PROTECT(duplicate(fac));
	    int *iPerm = Calloc(ncj, int);

	    ssc_metis_order(ncj, cp, ci, Perm, iPerm);
				/* apply to the crosstabulation, L, and ZZpO */
	    bCrosstab_permute(ZtZ, nf, j, nlev, iPerm);
	    bCrosstab_permute(L, nf, j, nlev, iPerm);
	    symmetric_permute(VECTOR_ELT(ZZpO, j), nlev[j], iPerm);
				/* apply to the factor */
	    factor_levels_permute(fac, fcp, Perm, iPerm);
				/* symbolic analysis to get Parent */
	    R_ldl_symbolic(ncj, cp, ci, Lp, INTEGER(VECTOR_ELT(parent, 0)), 
			 (int *) NULL, (int *) NULL);
	    for (i = 0; i < ncj; i++)
		INTEGER(VECTOR_ELT(parent, 1))[i] =
		    (INTEGER(VECTOR_ELT(parent, 0))[i] < 0) ? -1 : j;
	    nnz = Lp[ncj];
	    SET_SLOT(Ljj, Matrix_iSym, allocVector(INTSXP, nnz));
	    SET_SLOT(Ljj, Matrix_xSym,
		     alloc3Darray(REALSXP, nc[j], nc[j], nnz));
	    Free(iPerm); UNPROTECT(1);
	} else {
	    for (i = 0; i < ncj; i++) {
		Lp[i] = 0;
		INTEGER(VECTOR_ELT(parent,0))[i] = -1;
		INTEGER(VECTOR_ELT(parent,1))[i] = -1;
		Perm[i] = i;
	    }
	    Lp[ncj] = 0;
	    SET_SLOT(Ljj, Matrix_iSym, allocVector(INTSXP, 0));
	    SET_SLOT(Ljj, Matrix_xSym,
		     alloc3Darray(REALSXP, nc[j], nc[j], 0));
	}
	for (k = j+1; k < nf; k++) { /* Update other blocks in this column */
	    symbolic_right_unit_sm_trans(ncj, INTEGER(VECTOR_ELT(parent, 0)),
					 VECTOR_ELT(L, Lind(k,j)));
	}
	for (k = j+1; k < nf; k++) { /* Update remaining columns */
	    diag_update(VECTOR_ELT(ZZpO, k), VECTOR_ELT(L, Lind(k, j)));
	    for (i = k + 1; i < nf; i++)
		offdiag_update(VECTOR_ELT(L, Lind(k, j)),
			       VECTOR_ELT(L, Lind(i, j)),
			       VECTOR_ELT(L, Lind(i, k)),
			       nlev[k]);
	}
	check_x_slot_dim(VECTOR_ELT(ZZpO, j));
    }

    /* Convert blockwise Parent arrays to extended Parent arrays */
    for (j = 0; j < (nf - 1); j++) { /* Parent[nf] does not need conversion */
	SEXP Ljp1j = VECTOR_ELT(L, Lind(j + 1, j)),
	    LpP = GET_SLOT(Ljp1j, Matrix_pSym);
	int *Li = INTEGER(GET_SLOT(Ljp1j, Matrix_iSym)),
	    *Lp = INTEGER(LpP),
	    *block = INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, j), 1)),
	    *parent = INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, j), 0)),
	    i, nlev = length(LpP) - 1;
	for (i = 0; i < nlev; i++) {
	    if (block[i] < 0) {
		block[i] = j + 1;
		parent[i] = Li[Lp[i]];
	    }
	}
    }
    Free(nlev);
}
