#include "lmer.h"

/**
 * Check validity of an lmer object.
 *
 * @param x Pointer to an lmer object
 *
 * @return TRUE if the object is a valid lmer object, else a string
 * describing the nature of the violation.
 */
SEXP lmer_validate(SEXP x)
{
    SEXP
	/* ZZxP = GET_SLOT(x, Matrix_ZZxSym), */
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	XtXP = GET_SLOT(x, Matrix_XtXSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	RXXP = GET_SLOT(x, Matrix_RXXSym)
	/* , cnames = GET_SLOT(x, Matrix_cnamesSym) */
	;
    int *ZtXd = INTEGER(getAttrib(ZtXP, R_DimSymbol)),
	*XtXd = INTEGER(getAttrib(XtXP, R_DimSymbol));

    if (!(isReal(ZtXP) && isReal(XtXP) && isReal(RZXP) && isReal(RXXP) ))
	return mkString(_("Slots ZtX, XtX, RZX, and RXX must be real matrices"));
    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, R_DimSymbol))))
	return mkString(_("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, R_DimSymbol))))
	return mkString(_("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return mkString(_("Slots XtX must be a square matrix with same no. of cols as ZtX"));
    return ScalarLogical(1);
}

/**
 * Create the pairwise crosstabulation of the elements of flist.
 *
 * @param flist pointer to the factor list.
 * @param nobs number of observations.
 * @param nc number of columns in the model matrices.
 *
 * @return the pairwise crosstabulation in the form of the ZtZ array.
 * This version does not fill in the counts as they are not needed.
 */
static SEXP
internal_crosstab(SEXP flist, int nobs, const int nc[])
{
    int i, nf = length(flist);
    int npairs = (nf * (nf + 1))/2;
    SEXP val = PROTECT(allocVector(VECSXP, npairs));
    SEXP cscbCl = PROTECT(MAKE_CLASS("dgBCMatrix"));
    int *Ti = Calloc(nobs, int),
	*nlevs = Calloc(nf, int),
	**zb = Calloc(nf, int*); /* zero-based indices */

    for (i = 0; i < nf; i++) {	/* populate the zb vectors */
	SEXP fi = VECTOR_ELT(flist, i);
	int j;

	zb[i] = Calloc(nobs, int);
	nlevs[i] = length(getAttrib(fi, R_LevelsSymbol));
	for (j = 0; j < nobs; j++) zb[i][j] = INTEGER(fi)[j] - 1;
	for (j = 0; j <= i; j++) {
	    int *ijp, ind = Lind(i, j), nnz;
	    SEXP ZZij;

	    SET_VECTOR_ELT(val, ind, ZZij = NEW_OBJECT(cscbCl));
	    ijp = INTEGER(ALLOC_SLOT(ZZij, Matrix_pSym,
				     INTSXP, nlevs[j] + 1));
	    triplet_to_col(nlevs[i], nlevs[j], nobs,
			   zb[i], zb[j], (double *) NULL,
			   ijp, Ti, (double *) NULL);
	    nnz = ijp[nlevs[j]];
	    Memcpy(INTEGER(ALLOC_SLOT(ZZij, Matrix_iSym, INTSXP, nnz)),
		   Ti, nnz);
	}
    }

    for (i = 0; i < nf; i++) Free(zb[i]);
    Free(zb); Free(nlevs); Free(Ti);
    UNPROTECT(2);
    return val;
}

SEXP lmer_Crosstab(SEXP flist)
{
    SEXP val;
    int i, nf = length(flist), nobs;
    int *nc = Calloc(nf, int);

    if (!(nf > 0 && isNewList(flist)))
	error(_("flist must be a non-empty list"));
    nobs = length(VECTOR_ELT(flist, 0));
    if (nobs < 1) error(_("flist[[1]] must be a non-null factor"));
    for (i = 0; i < nf; i++) {
	SEXP fi = VECTOR_ELT(flist, i);
	if (!(isFactor(fi) && length(fi) == nobs))
	    error(_("flist[[%d]] must be a factor of length %d"),
		  i + 1, nobs);
	nc[i] = 1;
    }
    val = internal_crosstab(flist, nobs, nc);
    Free(nc);
    return val;
}

/** 
 * Allocate the x slot in an dgBCMatrix object
 * 
 * mm Pointer to a dgBCMatrix object
 * nr number of rows per block
 * nc number of columns per block
 */
#define ALLOC_X_SLOT(mm, nr, nc) \
    SET_SLOT(mm, Matrix_xSym, alloc3Darray(REALSXP, nr, nc, \
					   length(GET_SLOT(mm, Matrix_iSym))))

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
	int *mj = expand_cmprPt(ncol, cp, Calloc(nnz, int));
	int *mi = Memcpy(Calloc(nnz, int), INTEGER(cscbi), nnz);

	if (j <= jj) int_permute(mi, nnz, iperm);
	if (j >= jj) int_permute(mj, nnz, iperm);
	if (j == jj) make_upper_triangular(mi, mj, nnz);
	triplet_to_col(nrow, ncol, nnz, mi, mj, (double *) NULL,
		       cp, INTEGER(cscbi), (double *) NULL);
	Free(mi); Free(mj);
    }
}

/** 
 * Apply a permutation of the rows and columns to a sparse symmetric
 * matrix object.
 * 
 * @param A A sparse, symmetric matrix object stored in the upper
 * triangle
 * @param nlev order of A
 * @param iperm A 0-based permutation of length nlev
 */
static void
symmetric_permute(int Ap[], int Ai[], int n, const int iperm[])
{
    int nnz = Ap[n];
    int *mj = expand_cmprPt(n, Ap, Calloc(nnz, int));
    int *mi = Memcpy(Calloc(nnz, int), Ai, nnz);

    int_permute(mi, nnz, iperm);
    int_permute(mj, nnz, iperm);
    make_upper_triangular(mi, mj, nnz);
    triplet_to_col(n, n, nnz, mi, mj, (double *) NULL,
		   Ap, Ai, (double *) NULL);
    Free(mi); Free(mj);
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
    SEXP D, L, Parent, ZZpO, 
	flist = GET_SLOT(val, Matrix_flistSym),
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
	AZERO(REAL(VECTOR_ELT(D, j)), nc[j] * nc[j] * nlev[j]);
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
	    symmetric_permute(cp, ci, nlev[j], iPerm);
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
	}
	for (k = j+1; k < nf; k++) { /* Update other blocks in this column */
	    SEXP Lkj = VECTOR_ELT(L, Lind(k,j));
	    SET_SLOT(Lkj, Matrix_iSym,
		     lCholClgCsm(RGT, TRN, nlev[k], nlev[j],
				 INTEGER(VECTOR_ELT(parent, 0)),
				 GET_SLOT(Lkj, Matrix_iSym),
				 INTEGER(GET_SLOT(Lkj, Matrix_pSym))));
	}
	for (k = j + 1; k < nf; k++) { /* Update remaining columns */
	    SEXP db = VECTOR_ELT(ZZpO, k), Lkj = VECTOR_ELT(L, Lind(k, j));
	    int *Lkji = INTEGER(GET_SLOT(Lkj, Matrix_iSym)),
		*Lkjp = INTEGER(GET_SLOT(Lkj, Matrix_pSym));
	    SET_SLOT(db, Matrix_iSym,
		     Matrix_lgCsyrk(1, 0, nlev[k], nlev[j], Lkji, Lkjp,
				    1, GET_SLOT(db, Matrix_iSym),
				    INTEGER(GET_SLOT(db, Matrix_pSym))));
	    for (i = k + 1; i < nf; i++) {
		SEXP Lij = VECTOR_ELT(L, Lind(i, j)),
		    Lik = VECTOR_ELT(L, Lind(i, k));
		SET_SLOT(Lik, Matrix_iSym,
			 Matrix_lgClgCmm(0, 1, nlev[i], nlev[k], nlev[j],
					 INTEGER(GET_SLOT(Lij, Matrix_iSym)),
					 INTEGER(GET_SLOT(Lij, Matrix_pSym)),
					 Lkji, Lkjp,
					 1, GET_SLOT(Lik, Matrix_iSym),
					 INTEGER(GET_SLOT(Lik, Matrix_pSym))));
	    }
	}
    }
				
    for (j = 0; j < nf; j++) {	/* allocate x slots in dgBCMatrix objects */
	ALLOC_X_SLOT(VECTOR_ELT(ZZpO, j), nc[j], nc[j]);
	for (k = j; k < nf; k++) {
	    int indkj = Lind(k,j);
	    ALLOC_X_SLOT(VECTOR_ELT(L, indkj), nc[k], nc[j]);
	    ALLOC_X_SLOT(VECTOR_ELT(ZtZ, indkj), nc[k], nc[j]);
	}
    }
/* FIXME: Use these macros from Tim Davis instead */
#define EMPTY -1
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP(i) : (i))
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

/**
 * Update the arrays ZtZ, ZtX, and XtX in an lme object
 * according to a list of model matrices.
 *
 * @param x pointer to an lmer object
 * @param mmats pointer to a list of model matrices
 *
 * @return NULL
 */
SEXP lmer_update_mm(SEXP x, SEXP mmats)
{
    SEXP
	ZtZP = GET_SLOT(x, Matrix_ZtZSym),
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(ZtXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	nf = length(flist), nfp1 = nf + 1,
	i, ione = 1,
	nobs = nc[nfp1],
	pp1 = nc[nf];
    double
	*X,
	*XtX = REAL(GET_SLOT(x, Matrix_XtXSym)),
	*ZtX = REAL(ZtXP),
	one = 1.0, zero = 0.0;

    if (pp1 < 0) pp1 = -pp1;
    if (!isNewList(mmats) || length(mmats) != nfp1)
	error(_("mmats must be a list of %d model matrices"), nfp1);
    for (i = 0; i <= nf; i++) {
	SEXP mmat = VECTOR_ELT(mmats, i);
	int *mdims = INTEGER(getAttrib(mmat, R_DimSymbol));

	if (!isMatrix(mmat) || !isReal(mmat))
	    error(_("mmats[[%d]] is not a numeric matrix"), i + 1);
	if (nobs != mdims[0])
	    error(_("Expected %d rows in mmats[[%d]]. Got %d"),
		  nobs, i+1, mdims[0]);
	if ((nc[i] < 0 ? -nc[i] : nc[i]) != mdims[1])
	    error(_("Expected %d columns in mmats[[%d]]. Got %d"),
		  nc[i], i+1, mdims[1]);
    }
				/* Create XtX */
    X = REAL(VECTOR_ELT(mmats, nf));
    F77_CALL(dsyrk)("U", "T", &pp1, &nobs, &one, X, &nobs, &zero, XtX, &pp1);
				/* Zero an accumulator */
    AZERO(ZtX, pp1 * Gp[nf]);
    for (i = 0; i < nf; i++) {
	int *fac = INTEGER(VECTOR_ELT(flist, i)),
	    j, k, nci = nc[i], ZtXrows = Gp[i+1] - Gp[i];
	int ncisqr = nci * nci, nlev = ZtXrows/nci;
	double *Z = REAL(VECTOR_ELT(mmats, i)), *ZZx;

	for (k = 0; k < i; k++) {
	    SEXP ZZxM = VECTOR_ELT(ZtZP, Lind(i, k));
	    int *rowind = INTEGER(GET_SLOT(ZZxM, Matrix_iSym)),
		*colptr = INTEGER(GET_SLOT(ZZxM, Matrix_pSym));
	    int *f2 = INTEGER(VECTOR_ELT(flist, k)), nck = nc[k];
	    double *Zk = REAL(VECTOR_ELT(mmats, k));

	    ZZx = REAL(GET_SLOT(ZZxM, Matrix_xSym));
	    AZERO(ZZx, length(GET_SLOT(ZZxM, Matrix_xSym)));
	    for (j = 0; j < nobs; j++) {
		F77_CALL(dgemm)("T", "N", nc + i, nc + k, &ione, &one,
				Z + j, &nobs, Zk + j, &nobs, &one,
				ZZx + check_csc_index(colptr, rowind,
						      fac[j] - 1, f2[j] - 1, 0)
				* (nci * nck), &nci);
	    }
	}
	ZZx = REAL(GET_SLOT(VECTOR_ELT(ZtZP, Lind(i, i)), Matrix_xSym));
	AZERO(ZZx, nci * nci * nlev);
	if (nci == 1) {		/* single column in Z */
	    for (j = 0; j < nobs; j++) {
		int fj = fac[j] - 1; /* factor indices are 1-based */
		ZZx[fj] += Z[j] * Z[j];
		F77_CALL(daxpy)(&pp1, Z + j, X + j, &nobs, ZtX + fj, dims);
	    }
	} else {
	    for (j = 0; j < nobs; j++) {
		int fj = fac[j] - 1; /* factor indices are 1-based */

		F77_CALL(dsyr)("U", nc + i, &one, Z + j, &nobs,
			       ZZx + fj * ncisqr, nc + i);
		F77_CALL(dgemm)("T", "N", nc + i, &pp1, &ione,
				&one, Z + j, &nobs,
				X + j, &nobs, &one,
				ZtX + fj * nci, dims);
	    }
	}
	ZtX += ZtXrows;
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/**
 * Create an lmer object from a list of grouping factors and a list of model
 * matrices.  There is one more model matrix than grouping factor.  The last
 * model matrix is the fixed effects and the response.
 *
 * @param flist pointer to a list of grouping factors
 * @param mmats pointer to a list of model matrices
 *
 * @return pointer to an lmer object
 */
SEXP lmer_create(SEXP flist, SEXP mmats, SEXP method)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer")));
    SEXP ZtZ, cnames, fnms, nms;
    int *nc, i, nf = length(flist), nobs;

				/* Check validity of flist */
    if (!(nf > 0 && isNewList(flist)))
	error(_("flist must be a non-empty list"));
    nobs = length(VECTOR_ELT(flist, 0));
    if (nobs < 1) error(_("flist[[0]] must be a non-null factor"));
    for (i = 0; i < nf; i++) {
	SEXP fi = VECTOR_ELT(flist, i);
	if (!(isFactor(fi) && length(fi) == nobs))
	    error(_("flist[[%d]] must be a factor of length %d"),
		  i + 1, nobs);
    }
    SET_SLOT(val, Matrix_flistSym, duplicate(flist));
				/* Check mmats; allocate and populate nc */
    if (!(isNewList(mmats) && length(mmats) == (nf + 1)))
	error(_("mmats must be a list of length %d"), nf + 1);
    nc = INTEGER(ALLOC_SLOT(val, Matrix_ncSym, INTSXP, nf + 2));
    nc[nf + 1] = nobs;
    for (i = 0; i <= nf; i++) {
	SEXP mi = VECTOR_ELT(mmats, i);
	int *dims;

	if (!(isMatrix(mi) && isReal(mi)))
	    error(_("mmats[[%d]] must be a numeric matrix"), i + 1);
	dims = INTEGER(getAttrib(mi, R_DimSymbol));
	if (dims[0] != nobs)
	    error(_("mmats[[%d]] must have %d rows"), i + 1, nobs);
	if (dims[1] < 1)
	    error(_("mmats[[%d]] must have at least 1 column"), i + 1);
	nc[i] = dims[1];
    }   /* Arguments have now been checked for type, dimension, etc. */
				/* Create pairwise crosstabulation in ZtZ */
    SET_SLOT(val, Matrix_ZtZSym, internal_crosstab(flist, nobs, nc));
    SET_SLOT(val, Matrix_methodSym, duplicate(method));
    lmer_populate(val);
    ZtZ = GET_SLOT(val, Matrix_ZtZSym);
    /* FIXME: Check for possible reordering of the factors to maximize the
     * number of levels (columns?) in the leading nested sequence. */
    fnms = getAttrib(flist, R_NamesSymbol);
				/* Allocate and populate cnames */
    cnames = ALLOC_SLOT(val, Matrix_cnamesSym, VECSXP, nf + 1);
    setAttrib(cnames, R_NamesSymbol, allocVector(STRSXP, nf + 1));
    nms = getAttrib(cnames, R_NamesSymbol);
    for (i = 0; i <= nf; i++) {
	SEXP mi = VECTOR_ELT(mmats, i);
	SET_VECTOR_ELT(cnames, i,
		       duplicate(VECTOR_ELT(getAttrib(mi, R_DimNamesSymbol),
					    1)));
	SET_STRING_ELT(nms, i, (i < nf) ? duplicate(STRING_ELT(fnms, i)) :
		       mkChar(".fixed"));
    }
    lmer_update_mm(val, mmats);
    SET_SLOT(val, Matrix_bVarSym, duplicate(GET_SLOT(val, Matrix_DSym)));
    UNPROTECT(1);
    return val;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP lmer_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym);
    int	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)), i, nf = length(Omg);

    for (i = 0; i < nf; i++) {
	SEXP ZZxP = GET_SLOT(VECTOR_ELT(GET_SLOT(x, Matrix_ZtZSym), Lind(i, i)),
			     Matrix_xSym);
	int *dims = INTEGER(getAttrib(ZZxP, R_DimSymbol));
	int j, k, nzc = dims[0], nlev = dims[2];
	int nzcsqr = nzc * nzc, nzcp1 = nzc + 1;
	double *Omega = REAL(VECTOR_ELT(Omg, i)),
	    mi = 0.375 / ((double) nlev);

	AZERO(Omega, nzc * nzc);
	for (j = 0; j < nlev; j ++) {
	    for (k = 0; k < nzc; k++) {
		Omega[k * nzcp1] += REAL(ZZxP)[k * nzcp1 + j * nzcsqr] * mi;
	    }
	}
    }
    status[0] = status[1] = 0;
    return R_NilValue;
}

/**
 * Copy ZtZ to ZZpO and L.  Inflate diagonal blocks of ZZpO by Omega.
 * Update devComp[1].
 *
 * @param x pointer to an lmer object
 */
SEXP
lmer_inflate(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym),
	ZZpO = GET_SLOT(x, Matrix_ZZpOSym),
	ZtZ = GET_SLOT(x, Matrix_ZtZSym),
	LP = GET_SLOT(x, Matrix_LSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, k, nf = length(Omg);
    double *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym));

    for (i = 0; i < nf; i++) {
	SEXP ZZOel = VECTOR_ELT(ZZpO, i);
	SEXP ZZOm = GET_SLOT(ZZOel, Matrix_xSym);
	SEXP ZZel = VECTOR_ELT(ZtZ, Lind(i, i));
	int *Di = INTEGER(GET_SLOT(ZZOel, Matrix_iSym)),
	    *Dp = INTEGER(GET_SLOT(ZZOel, Matrix_pSym)),
	    *Si = INTEGER(GET_SLOT(ZZel, Matrix_iSym)),
	    *Sp = INTEGER(GET_SLOT(ZZel, Matrix_pSym)),
	    *dims = INTEGER(getAttrib(ZZOm, R_DimSymbol));
	int sz = dims[0] * dims[1];
	int ii, j, nci = nc[i], ncisqr = nci * nci;
	int nlev = (Gp[i + 1] - Gp[i])/nci;
	double *Omega = REAL(VECTOR_ELT(Omg, i)),
	    *ZZ = REAL(GET_SLOT(ZZel, Matrix_xSym)),
	    *tmp = Memcpy(Calloc(ncisqr, double), Omega, ncisqr);

	F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
	if (j)
	    error(_("Leading %d minor of Omega[[%d]] not positive definite"),
		  j, i + 1);
				/* update dcmp[1] */
	for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
	    dcmp[1] += nlev * 2. * log(tmp[j * (nci + 1)]);
	}
	Free(tmp);
	AZERO(REAL(ZZOm), dims[0] * dims[1] * dims[2]);
	for (j = 0; j < nlev; j++) { /* copy diagonal block and inflate */
	    double *ZZOkk = REAL(ZZOm) + check_csc_index(Dp, Di, j, j, 0) * sz;
	    int kk, k2 = Sp[j + 1];
	    for (kk = Sp[j]; kk < k2; kk++) {
		Memcpy(REAL(ZZOm) + check_csc_index(Dp, Di, Si[kk], j, 0) * sz,
		       ZZ + kk * sz, sz);
	    }
	    for (kk = 0; kk < nci; kk++) {
		for (ii = 0; ii <= kk; ii++) {
		    int ind = ii + kk * nci;
		    ZZOkk[ind] += Omega[ind];
		}
	    }
	}
	for (k = i + 1; k < nf; k++) {
	    int ind = Lind(k, i);
	    SEXP Lel = VECTOR_ELT(LP, ind),
		Lm = GET_SLOT(Lel, Matrix_xSym);
	    double *L = REAL(Lm);

	    dims = INTEGER(getAttrib(Lm, R_DimSymbol));
	    ZZel = VECTOR_ELT(ZtZ, ind);
	    ZZ = REAL(GET_SLOT(ZZel, Matrix_xSym));
	    Di = INTEGER(GET_SLOT(Lel, Matrix_iSym));
	    Dp = INTEGER(GET_SLOT(Lel, Matrix_pSym));
	    Si = INTEGER(GET_SLOT(ZZel, Matrix_iSym));
	    Sp = INTEGER(GET_SLOT(ZZel, Matrix_pSym));
	    sz = dims[0] * dims[1];

	    AZERO(L, sz * dims[2]); /* zero L  */
	    for (j = 0; j < nlev; j++) { /* copy src blocks to dest */
		int kk, k2 = Sp[j + 1];
		for (kk = Sp[j]; kk < k2; kk++) {
		    Memcpy(L + check_csc_index(Dp, Di, Si[kk], j, 0) * sz,
			   ZZ + kk * sz, sz);
		}
	    }
	}
    }
    return R_NilValue;
}

/**
 * Convert the extended parent pair (Parent, Block) to a parent array
 * for the jth diagonal block of size n.
 *
 * @param j index (0-based) of the diagonal outer block
 * @param n number of inner column blocks in the outer block
 * @param par array of length n to be filled with the parent array
 * @param ParP pointer to the extended parent structure
 *
 * @return par
 */
static R_INLINE
int *block_parent(int j, int n, int par[], SEXP ParP)
{
    SEXP Parj = VECTOR_ELT(ParP, j);
    int *Parent = INTEGER(VECTOR_ELT(Parj, 0)),
	*Block = INTEGER(VECTOR_ELT(Parj, 1)), i;
    for (i = 0; i < n; i++) par[i] = (Block[i] == j) ? Parent[i] : -1;
    return par;
}

/**
 * If status[["factored"]] is FALSE, create and factor Z'Z+Omega.  Also
 * create RZX and RXX, the deviance components, and the value of the
 * deviance for both ML and REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP lmer_factor(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));

    if (!status[0]) {
	SEXP DP = GET_SLOT(x, Matrix_DSym),
	    LP = GET_SLOT(x, Matrix_LSym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    ZZOP = GET_SLOT(x, Matrix_ZZpOSym),
	    Parent = GET_SLOT(x, Matrix_ParentSym);
	int *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    i, j, nf = length(DP);
	int nml = nc[nf + 1], nreml = nml + 1 - nc[nf], ncX = dims[1];
	double
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXP),
	    *ZtX = REAL(GET_SLOT(x, Matrix_ZtXSym)),
	    *XtX = REAL(GET_SLOT(x, Matrix_XtXSym)),
	    *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	    *deviance = REAL(GET_SLOT(x, Matrix_devianceSym)),
	    minus1 = -1., one = 1.;

	if (nc[nf] < 0) {	/* indicates that fixed effects should be skipped */
	    int ZXshift = dims[0] * (dims[1] - 1),
		XXshift = dims[1]*dims[1] - 1;
	    
	    RZX += ZXshift; ZtX += ZXshift; RXX += XXshift; XtX += XXshift;
	    ncX = 1;
	    nreml = nml + 1 + nc[nf];
	}
	dcmp[0] = dcmp[1] = dcmp[2] = dcmp[3] = 0.;
	Memcpy(RZX, ZtX, dims[0] * ncX);
	lmer_inflate(x);	/* initialize ZZpO and L */
	for (i = 0; i < nf; i++) {
	    SEXP ZZOiP = VECTOR_ELT(ZZOP, i);
	    SEXP DiP = VECTOR_ELT(DP, i);
	    SEXP LiP = VECTOR_ELT(LP, Lind(i, i));
	    int nlev = INTEGER(getAttrib(DiP, R_DimSymbol))[2];
	    int jj, nci = nc[i], ncisqr = nci * nci;
	    int *Pari = block_parent(i, nlev, Calloc(nlev, int), Parent);
	    double *D = REAL(DiP);

	    jj = cscb_ldl(ZZOiP, Pari, LiP, DiP);
	    if (jj != nlev) error(_("cscb_ldl returned %d < nlev = %d"), jj, nlev);
	    for (j = 0; j < nlev; j++) { /* accumulate dcmp[0] */
		double *Dj = D + j * ncisqr;
		for (jj = 0; jj < nci; jj++) /* accumulate determinant */
		    dcmp[0] += 2. * log(Dj[jj * (nci + 1)]);
	    }
	    /* Solve L_{i,i} %*% RZX_i := RZX_i */
	    cscb_trsm(LOW, NTR, UNT, 1., LiP,
		      Gp[i+1] - Gp[i], ncX, RZX + Gp[i], dims[0]);
	    /* Solve D_i^{T/2} %*% RZX_i := RZX_i */
	    for (jj = 0; jj < nlev; jj++) {
		F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &ncX,
				&one, D + jj * ncisqr, &nci,
				RZX + Gp[i] + jj * nci, &dims[0]);
	    }
	    for (j = i + 1; j < nf; j++) { /*  further blocks */
		SEXP Lji = VECTOR_ELT(LP, Lind(j, i));
		SEXP Lx = GET_SLOT(Lji, Matrix_xSym);
		double *L = REAL(Lx);
		int *xdims = INTEGER(getAttrib(Lx, R_DimSymbol)),
		    *Lp = INTEGER(GET_SLOT(Lji, Matrix_pSym));
		int ntot = xdims[0] * xdims[1];

		/* L_{j,i} := L_{j,i} %*% L_{i,i}^{-T} %*% D_i^{-1/2} */
		cscb_trcbsm(RGT, LOW, TRN, UNT, 1.0, LiP, Pari, Lji);
		for (jj = 0; jj < nlev; jj++) {
		    int k, k2 = Lp[jj + 1];
		    for (k = Lp[jj]; k < k2; k++)
			F77_CALL(dtrsm)("R", "U", "N", "N", xdims, xdims + 1,
					&one, D + jj * ncisqr, &nci,
					L + k * ntot, xdims);
		}
		/* RZX_j := RZX_j - (L_{j,i} %*% D_i^{T/2}) %*% RZX_i */
		/* At this point Lji contains L_{j,i} %*% D_i^{T/2} */
		cscb_mm(LFT, NTR, Gp[j + 1] - Gp[j], ncX, Gp[i+1] - Gp[i],
			-1.0, Lji, RZX + Gp[i], dims[0],
			1.0, RZX + Gp[j], dims[0]);
	    }
	    for (j = i + 1; j < nf; j++) { /* block pairs and final update */
		SEXP Lji = VECTOR_ELT(LP, Lind(j, i));
		SEXP Lx = GET_SLOT(Lji, Matrix_xSym);
		double *L = REAL(Lx);
		int *xdims = INTEGER(getAttrib(Lx, R_DimSymbol)),
		    *Lp = INTEGER(GET_SLOT(Lji, Matrix_pSym));
		int ntot = xdims[0] * xdims[1];


		/* ZZpO_{j,j} := ZZpO_{j,j} - L{j,i} %*% L_{j,i}^T */
		cscb_syrk(UPP, NTR, -1.0, Lji, 1.0, VECTOR_ELT(ZZOP, j));
		for (jj = j+1; jj < nf; jj++) {
		    /* L_{jj,j} := L_{jj,j} - L{jj,i} %*% L_{j,i}^T */
		    cscb_cscbm(NTR, TRN, -1.0, VECTOR_ELT(LP, Lind(jj, i)),
			    Lji, 1.0, VECTOR_ELT(LP, Lind(jj, j)));
		}
		/* L_{j,i} := L_{j,i} %*% D_i^{-T/2} */
		for (jj = 0; jj < nlev; jj++) {
		    int k, k2 = Lp[jj + 1];
		    for (k = Lp[jj]; k < k2; k++)
			F77_CALL(dtrsm)("R", "U", "T", "N", xdims, xdims + 1,
					&one, D + jj * ncisqr, &nci,
					L + k * ntot, xdims);
		}
	    }
	    Free(Pari);
	}
				/* downdate and factor XtX */
	Memcpy(RXX, XtX, ncX * ncX);
	F77_CALL(dsyrk)("U", "T", &ncX, &dims[0],
			&minus1, RZX, &dims[0], &one, RXX, &ncX);
	F77_CALL(dpotrf)("U", &ncX, RXX, &ncX, &j);
	if (j) {
	    warning(_("Leading minor of size %d of downdated X'X is indefinite"),
		    j);
	    dcmp[2] = dcmp[3] = deviance[0] = deviance[1] = NA_REAL;
	} else {
	    for (j = 0; j < (ncX - 1); j++) /* 2 logDet(RXX) */
		dcmp[2] += 2 * log(RXX[j * (dims[1] + 1)]);
	    dcmp[3] = 2. * log(RXX[ncX * ncX - 1]); /* 2 log(ryy) */
	    deviance[0] =	/* ML criterion */
		dcmp[0] - dcmp[1] + nml*(1.+dcmp[3]+log(2.*PI/nml));
	    deviance[1] = dcmp[0] - dcmp[1] + /* REML */
		dcmp[2] + nreml*(1.+dcmp[3]+log(2.*PI/nreml));
	}
	status[0] = 1; status[1] = 0; /* factored but not inverted */
    }
    return R_NilValue;
}

/**
 * Solve one of the matrix equations op(L)*X=alpha*B or
 * X*op(L)=alpha*B where L is a sparse, blocked, unit lower triangular matrix.
 *
 * @param side LFT or RGT for left or right
 * @param trans TRN or NTR for transpose or no transpose
 * @param nf number of grouping factors
 * @param Gp group pointers for the rows
 * @param n number of columns
 * @param alpha multiplier
 * @param L pointer to the L cscb object
 * @param B pointer to the matrix of right-hand sides
 * @param ldb leading dimension of array B as declared in the caller
 */
static void
internal_sm(enum CBLAS_SIDE side, enum CBLAS_TRANSPOSE trans, int nf,
	    const int Gp[], int n, double alpha, SEXP L, double B[], int ldb)
{
    int j, k;

    if (side == LFT) {
	if (trans == TRN) {
	    for (j = nf - 1; j >= 0; j--) {
		int nrj = Gp[j + 1] - Gp[j];

		cscb_trsm(LOW, TRN, UNT, alpha, VECTOR_ELT(L, Lind(j, j)),
			  nrj, n, B + Gp[j], ldb);
		for (k = 0; k < j; k++) {
		    cscb_mm(LFT, TRN, Gp[k + 1] - Gp[k], n, nrj,
			    -1., VECTOR_ELT(L, Lind(j, k)),
			    B + Gp[j], ldb, alpha, B + Gp[k], ldb);
		}
	    }
	} else error(_("Code for non-transpose case not yet written"));
    } else error(_("Code for right-side solutions not yet written"));
}

/** 
 * Determine the maximum number of nonzero elements in a column and
 * allocate storage for the tmp and ind arrays.
 * 
 * @param j level
 * @param Parent Parent list
 * 
 * @return Maximum number of nonzero elements in a column
 */
static void
alloc_tmp_ind(int nf, const int nc[], const int nlevs[], SEXP Parent,
	      double *tmp[], int *ind[])
{
    int j, maxnc;
    for (maxnc = -1, j = 0; j < nf; j++) {
	SEXP lst = VECTOR_ELT(Parent, j);
	SEXP blk = VECTOR_ELT(lst, 1), par = VECTOR_ELT(lst, 0);
	int *nfj = Calloc(nlevs[j], int), i, val;

	
	if (nc[j] > maxnc) maxnc = nc[j];
	for (val = -1, i = nlevs[j] - 1; i >= 0; i--) {
	    int thisnnz = (INTEGER(blk)[i] != j) ? 1 : nfj[INTEGER(par)[i]] + 1;
	    if (thisnnz > val) val = thisnnz;
	    nfj[i] = thisnnz;
	}
	ind[j] = Calloc(val, int);
	tmp[j] = Calloc(val * nc[j] * maxnc, double);
	Free(nfj);
    }
}

#define BLK(i,j) INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, i), 1))[j]
#define PAR(i,j) INTEGER(VECTOR_ELT(VECTOR_ELT(Parent, i), 0))[j]

/**
 * Fill the nnz array with the number of nonzero inner blocks in each
 * outer block of the jth inner column block of the ith outer block of
 * L^{-1}.  Also fill the ind array.
 *
 * @param i outer block index
 * @param j inner block index within the ith outer block
 * @param nf number of factors
 * @param Parent pointer to the extended parent pairs
 * @param nc
 * @param nnz array of length nf
 * @param tmp array of length nf of pointers to doubles
 * @param ind array of length nf of pointers to ints
 *
 */
static
void fill_ind(int i, int j, int nf, SEXP Parent, int nnz[], int *ind[])
{
    int blk, k, par;

    AZERO(nnz, nf);
    for (blk = BLK(i,j), par = PAR(i,j); blk >= 0;
	 k = BLK(blk,par), par = PAR(blk,par), blk = k) {
	ind[blk][nnz[blk]++] = par;
    }
}

static R_INLINE
int fsrch(int target, const int vals[], int nvals)
{
    int i;
    for (i = 0; i < nvals; i++) if (vals[i] == target) return i;
    error(_("fsrch: unable to find target %d in nvals %d "), target, nvals);
    return -1;			/* -Wall */
}

/**
 * If necessary, factor Z'Z+Omega, ZtX, and XtX then, if necessary,
 * replace the RZX and RXX slots by the corresponding parts of the
 * inverse of the Cholesky factor.  Replace the elements of the D slot
 * by the blockwise inverses and evaluate bVar.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL (x is updated in place)
 */
SEXP lmer_invert(SEXP x)
{
    int *status = LOGICAL(GET_SLOT(x, Matrix_statusSym));
    if (!status[0]) lmer_factor(x);
    if (!R_FINITE(REAL(GET_SLOT(x, Matrix_devianceSym))[0]))
	error(_("Unable to invert singular factor of downdated X'X"));
    if (!status[1]) {
	SEXP DP = GET_SLOT(x, Matrix_DSym),
	    LP = GET_SLOT(x, Matrix_LSym),
	    ParP = GET_SLOT(x, Matrix_ParentSym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    bVarP = GET_SLOT(x, Matrix_bVarSym);
	int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i, nf = length(DP);
	int **ind = Calloc(nf, int *),
	    *nlevs = Calloc(nf, int),
	    *nnz = Calloc(nf, int), ncX = dims[1];
	double **tmp = Calloc(nf, double *),
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXP),
	    minus1 = -1., one = 1., zero = 0.;

	if (nc[nf] < 0) {	/* indicates that fixed effects should be skipped */
	    RZX += dims[0] * (dims[1] - 1);
	    RXX += dims[1]*dims[1] - 1;
	    ncX = 1;
	}
	/* RXX := RXX^{-1} */
	F77_CALL(dtrtri)("U", "N", &ncX, RXX, &ncX, &i);
	if (i)
	    error(_("Leading minor of size %d of downdated X'X,is indefinite"),
		  i + 1);

	/* RZX := - RZX %*% RXX */
	F77_CALL(dtrmm)("R", "U", "N", "N", &dims[0], &ncX, &minus1,
			RXX, &ncX, RZX, &dims[0]);
	for(i = 0; i < nf; i++) {
	    int info, j, jj, nci = nc[i];
	    int ncisqr = nci * nci;
	    double *Di = REAL(VECTOR_ELT(DP, i)),
		*RZXi = RZX + Gp[i];

	    nlevs[i] = (Gp[i+1] - Gp[i])/nci;
	    /* D_i := D_i^{-1}; RZX_i := D_i %*% RZX_i */
	    if (nci == 1) {
		for (j = 0; j < nlevs[i]; j++) {
		    Di[j] = 1./Di[j];
		    for (jj = 0; jj < ncX; jj++)
			RZXi[j + jj * dims[0]] *= Di[j];
		}
	    } else {
		for (j = 0; j < nlevs[i]; j++) {
		    F77_CALL(dtrtri)("U", "N", &nci, Di + j * ncisqr, &nci, &info);
		    if (info)
			error(_("D[,,%d] for factor %d is singular"), j + 1, i + 1);
		    F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &ncX, &one,
				    Di + j * ncisqr, &nci, RZXi + j * nci, &dims[0]);
		}
	    }
	}

	/* RZX := L^{-T} %*% RZX */
	internal_sm(LFT, TRN, nf, Gp, ncX, 1.0, LP, RZX, dims[0]);

	alloc_tmp_ind(nf, nc, nlevs, ParP, tmp, ind);
	/* Create bVar arrays as crossprod of column blocks of D^{-T/2}%*%L^{-1} */
	for (i = 0; i < nf; i++) { /* ith column of outer blocks */
	    int j, k, kj, nci = nc[i];
	    int ncisqr = nci * nci;
	    double *Di = REAL(VECTOR_ELT(DP, i)),
		*bVi = REAL(VECTOR_ELT(bVarP, i));

	    AZERO(bVi, ncisqr * nlevs[i]);
	    for (j = 0; j < nlevs[i]; j++) {
		double *bVij = bVi + j * ncisqr, *Dij = Di + j * ncisqr;
		
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &one, Dij,
				&nci, &zero, bVij, &nci);
		/* count non-zero blocks; allocate and zero storage */
		fill_ind(i, j, nf, ParP, nnz, ind);

		for (k = i; k < nf; k++) { /* kth row of outer blocks */
		    SEXP Lki = VECTOR_ELT(LP, Lind(k, i));
		    int *Lkii = INTEGER(GET_SLOT(Lki, Matrix_iSym)),
			*Lkip = INTEGER(GET_SLOT(Lki, Matrix_pSym));
		    double *Lkix = REAL(GET_SLOT(Lki, Matrix_xSym));
		    int kk, sz = nc[i] * nc[k];
		    
		    AZERO(tmp[k], sz * nnz[k]);
		    /* initialize tmp from jth column of (k,i)th block */
		    /* - sign in sol'n incorporated in dtrmm call below */
		    for (kk = Lkip[j]; kk < Lkip[j + 1]; kk++)
			Memcpy(tmp[k] + fsrch(Lkii[kk], ind[k], nnz[k]) * sz,
			       Lkix + kk * sz, sz);
		    /* columns in ind[kk] for (k,kk)th block */
		    for (kk = i; kk <= k; kk++) {
			int szk = nc[k] * nc[kk];
			/* skip getting slots if not using them */
			if (!nnz[kk]) continue;
			Lki = VECTOR_ELT(LP, Lind(k, kk));
			Lkii = INTEGER(GET_SLOT(Lki, Matrix_iSym));
			Lkip = INTEGER(GET_SLOT(Lki, Matrix_pSym));
			Lkix = REAL(GET_SLOT(Lki, Matrix_xSym));
			for (kj = 0; kj < nnz[kk]; kj++) {
			    int col = ind[kk][kj], k1, szkk = nc[i] * nc[kk];
			    
			    for (k1 = Lkip[col]; k1 < Lkip[col + 1]; k1++) {
				if ((kk == k) && col >= Lkii[k1]) break;
				F77_CALL(dgemm)("N", "N", &nc[k], &nci, &nc[kk],
						&minus1, Lkix + k1 * szk,
						&nc[k], tmp[kk] + kj * szkk,
						&nc[kk], &one,
						tmp[k] +
						fsrch(Lkii[k1],ind[k],nnz[k])*sz,
						&nc[k]);
			    }
			}
		    }
		}
		for (k = 0; k < nf; k++) {
		    for (kj = 0; kj < nnz[k]; kj++) {
			F77_CALL(dtrmm)("L", "U", "T", "N", nc + k, &nci, &minus1,
					REAL(VECTOR_ELT(DP, k))+ind[k][kj]*nc[k]*nc[k],
					nc + k, tmp[k] + kj * nc[i] * nc[k],
					nc + k);
		    }
		    if (nnz[k] > 0) {
			kj = nc[k] * nnz[k];
			F77_CALL(dsyrk)("U", "T", &nci, &kj, &one, tmp[k], &kj,
					&one, bVij, &nci);
		    }
		}
	    }
	}
	for (i = 0; i < nf; i++) {
	    if (tmp[i]) Free(tmp[i]);
	    if (ind[i]) Free(ind[i]);
	}
	Free(tmp); Free(nlevs); Free(nnz); Free(ind);
	status[1] = 1;
    }
    return R_NilValue;
}

static double
internal_sigma(SEXP x, int REML)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1],
	nobs = INTEGER(GET_SLOT(x, Matrix_ncSym))
	[length(GET_SLOT(x, Matrix_OmegaSym)) + 1];

    lmer_invert(x);
    return 1./(REAL(RXXsl)[pp1*pp1 - 1] *
	       sqrt((double)(REML ? nobs + 1 - pp1 : nobs)));
}

/**
 * Extract the ML or REML conditional estimate of sigma
 *
 * @param x pointer to an lme object
 * @param REML logical scalar - TRUE if REML estimates are requested
 *
 * @return pointer to a numeric scalar
 */
SEXP lmer_sigma(SEXP x, SEXP REML)
{
    return ScalarReal(internal_sigma(x, asLogical(REML)));
}

/**
 * Calculate the length of the parameter vector (historically called "coef"
 * even though these are not coefficients).
 *
 * @param nf number of factors
 * @param nc number of columns in the model matrices for each factor
 *
 * @return total length of the coefficient vector
 */
static R_INLINE
int coef_length(int nf, const int nc[])
{
    int i, ans = 0;
    for (i = 0; i < nf; i++) ans += (nc[i] * (nc[i] + 1))/2;
    return ans;
}

/**
 * Extract parameters from the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 * 
 * The parameters can be returned in three forms:
 *   0 - nonlinearly constrained - elements of the relative precision matrix
 *   1 - unconstrained - from the LDL' decomposition - logarithms of
 *       the diagonal elements of D
 *   2 - box constrained - also from the LDL' decomposition - inverses
 *       of the diagonal elements of D
 *
 * @param x pointer to an lme object
 * @param pType pointer to an integer scalar indicating the form of the 
 *        parameters to be returned.
 *
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP lmer_coef(SEXP x, SEXP pType)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), ptyp = asInteger(pType), vind;
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double *vv = REAL(val);

    vind = 0;			/* index in vv */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    double dd = REAL(VECTOR_ELT(Omega, i))[0];
	    vv[vind++] = ptyp ? ((ptyp == 1) ? log(dd) : 1./dd) : dd;
	} else {
	    if (ptyp) {	/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
		double *tmp = Memcpy(Calloc(ncisq, double),
				     REAL(VECTOR_ELT(Omega, i)), ncisq);
		F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
		if (j)		/* should never happen */
		    error(_("DPOTRF returned error code %d on Omega[[%d]]"),
			  j, i+1);
		for (j = 0; j < nci; j++) {
		    double diagj = tmp[j * ncip1];
		    vv[vind++] = (ptyp == 1) ? (2. * log(diagj)) :
			1./(diagj * diagj);
		    for (k = j + 1; k < nci; k++) {
			tmp[j + k * nci] /= diagj;
		    }
		}
		for (j = 0; j < nci; j++) {
		    for (k = j + 1; k < nci; k++) {
			vv[vind++] = tmp[j + k * nci];
		    }
		}
		Free(tmp);
	    } else {		/* upper triangle of Omega[,,i] */
		int j, k, odind = vind + nci;
		double *omgi = REAL(VECTOR_ELT(Omega, i));

		for (j = 0; j < nci; j++) {
		    vv[vind++] = omgi[j * ncip1];
		    for (k = j + 1; k < nci; k++) {
			vv[odind++] = omgi[k*nci + j];
		    }
		}
		vind = odind;
	    }
	}
    }
    UNPROTECT(1);
    return val;
}

static
void internal_coefGets(SEXP x, const double cc[], int ptyp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	cind, i, nf = length(Omega);

    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    double dd = cc[cind++];
	    REAL(VECTOR_ELT(Omega, i))[0] =
		ptyp ? ((ptyp == 1) ? exp(dd) : 1./dd) : dd;
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		j, k,
		ncip1 = nci + 1,
		ncisq = nci * nci;
	    double
		*omgi = REAL(VECTOR_ELT(Omega, i));
	    if (ptyp) {
		double *tmp = Calloc(ncisq, double),
		    diagj, one = 1., zero = 0.;

		AZERO(omgi, ncisq);
		for (j = 0; j < nci; j++) {
		    double dd = cc[cind++];
		    tmp[j * ncip1] = diagj =
			(ptyp == 1) ? exp(dd/2.) : sqrt(1./dd);
		    for (k = j + 1; k < nci; k++) {
			tmp[k*nci + j] = cc[odind++] * diagj;
		    }
		}
		F77_CALL(dsyrk)("U", "T", &nci, &nci, &one,
				tmp, &nci, &zero, omgi, &nci);
		Free(tmp);
	    } else {
		for (j = 0; j < nci; j++) {
		    omgi[j * ncip1] = cc[cind++];
		    for (k = j + 1; k < nci; k++) {
			omgi[k*nci + j] = cc[odind++];
		    }
		}
	    }
	    cind = odind;
	}
    }
    status[0] = status[1] = 0;
}

/**
 * Assign the upper triangles of the Omega matrices according to a
 * vector of parameters.
 *
 * @param x pointer to an lme object
 * @param coef pointer to an numeric vector of appropriate length
 * @param pType pointer to an integer scalar 
 *
 * @return R_NilValue
 */
SEXP lmer_coefGets(SEXP x, SEXP coef, SEXP pType)
{
    int clen = coef_length(LENGTH(GET_SLOT(x, Matrix_flistSym)),
			   INTEGER(GET_SLOT(x, Matrix_ncSym)));   
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_coefGets(x, REAL(coef), asInteger(pType));
    return x;
}

static double*
internal_fixef(SEXP x, double beta[])
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int j, pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1];
    int p = pp1 - 1;
    double nryyinv;		/* negative ryy-inverse */

    lmer_invert(x);
    Memcpy(beta, REAL(RXXsl) + pp1 * (pp1 - 1), p);
    nryyinv = -REAL(RXXsl)[pp1*pp1 - 1];
    for (j = 0; j < p; j++) beta[j] /= nryyinv;
    return beta;
}

/**
 * Extract the conditional estimates of the fixed effects
 *
 * @param x Pointer to an mer object
 *
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP lmer_fixef(SEXP x)
{
    SEXP cnames = GET_SLOT(x, Matrix_cnamesSym);
    int *dims = INTEGER(getAttrib(GET_SLOT(x, Matrix_RZXSym), R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(cnames) - 1;
    int i, p = dims[1] - 1;
    SEXP ans, nms, cnms;
		
    if (nc[nf] < 0 || !p) return(allocVector(REALSXP, 0));
    ans = PROTECT(allocVector(REALSXP, p));
    internal_fixef(x, REAL(ans));
    nms = PROTECT(allocVector(STRSXP, p));
    cnms = VECTOR_ELT(cnames, nf);
    for (i = 0; i < p; i++) SET_STRING_ELT(nms, i, STRING_ELT(cnms, i));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Extract the conditional modes of the random effects.
 *
 * @param x Pointer to an lme object
 *
 * @return a list of matrices containing the conditional modes of the random effects
 */
SEXP lmer_ranef(SEXP x)
{
    SEXP RZXP = GET_SLOT(x, Matrix_RZXSym),
	cnames = GET_SLOT(x, Matrix_cnamesSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ii, jj,
	nf = length(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double
	*b = REAL(RZXP) + dims[0] * (dims[1] - 1),
	nryyinv;		/* negative ryy-inverse */

    lmer_invert(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
    nryyinv = -REAL(GET_SLOT(x, Matrix_RXXSym))[dims[1] * dims[1] - 1];
    for (i = 0; i < nf; i++) {
	SEXP nms, rnms = getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol);
	int nci = nc[i], mi = length(rnms);
	double *bi = b + Gp[i], *mm;

	SET_VECTOR_ELT(val, i, allocMatrix(REALSXP, mi, nci));
	setAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol, allocVector(VECSXP, 2));
	nms = getAttrib(VECTOR_ELT(val, i), R_DimNamesSymbol);
	SET_VECTOR_ELT(nms, 0, duplicate(rnms));
	SET_VECTOR_ELT(nms, 1, duplicate(VECTOR_ELT(cnames, i)));
	mm = REAL(VECTOR_ELT(val, i));
	for (jj = 0; jj < nci; jj++)
	    for(ii = 0; ii < mi; ii++)
		mm[ii + jj * mi] = bi[jj + ii * nci]/nryyinv;
    }
    UNPROTECT(1);
    return val;
}

/**
 * Fill in four symmetric matrices for each level, providing the
 * information to generate the gradient or the ECME step.  The four
 * matrices are
 *  1) -m_i\bOmega_i^{-1}
 *  2) \bB_i\bB_i\trans
 *  3) \tr\left[\der_{\bOmega_i}\bOmega\left(\bZ\trans\bZ+\bOmega\right)\inv\right]
 *  4) The term added to 3) to get \tr\left[\der_{\bOmega_i}\bOmega\vb\right]
 *
 * @param x pointer to an lme object
 * @param val pointer to a list of matrices of the correct sizes
 *
 * @return val
 */
/* static */
SEXP lmer_firstDer(SEXP x, SEXP val)
{
    SEXP bVarP = GET_SLOT(x, Matrix_bVarSym),
	OmegaP = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym);
    int *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	i, nf = length(OmegaP), p = dims[1] - 1;
    double *RZX = REAL(RZXP),
	*b = REAL(RZXP) + dims[0] * p;

    lmer_invert(x);
    /* FIXME: Why is this loop run backwards?  It appears it could run forwards. */
    for (i = nf - 1; i >= 0; i--) {
	SEXP bVPi = VECTOR_ELT(bVarP, i);
	int *ddims = INTEGER(getAttrib(bVPi, R_DimSymbol)), j, k;
	int nci = ddims[0];
	int ncisqr = nci * nci, RZXrows = Gp[i + 1] - Gp[i];
	int nlev = RZXrows/nci;
	double *RZXi = RZX + Gp[i], *bVi = REAL(bVPi),
	    *bi = b + Gp[i], *mm = REAL(VECTOR_ELT(val, i)),
	    *tmp = Memcpy(Calloc(ncisqr, double),
			  REAL(VECTOR_ELT(OmegaP, i)), ncisqr),
	    dlev = (double) nlev,
	    one = 1., zero = 0.;

 	if (nci == 1) {
	    int ione = 1;
 	    mm[0] = ((double) nlev)/tmp[0];
 	    mm[1] = F77_CALL(ddot)(&nlev, bi, &ione, bi, &ione);
	    mm[2] = 0.;
	    for (k = 0; k < nlev; k++) mm[2] += bVi[k];
	    mm[3] = 0.;
  	    for (j = 0; j < p; j++) {
  		mm[3] += F77_CALL(ddot)(&RZXrows, RZXi + j * dims[0], &ione,
					RZXi + j * dims[0], &ione);
  	    }
 	} else {
	    AZERO(mm, 4 * ncisqr);
	    F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
	    if (j)
		error(_("Omega[[%d]] is not positive definite"), i + 1);
	    F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
	    if (j)
		error(_("Omega[[%d]] is not positive definite"), i + 1);
	    F77_CALL(dsyrk)("U", "N", &nci, &nci, &dlev, tmp, &nci,
			    &zero, mm, &nci);
	    mm += ncisqr;	/* \bB_i term */
	    F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one, bi, &nci,
			    &zero, mm, &nci);
	    mm += ncisqr;     /* Sum of diagonal blocks of the inverse
			       * (Z'Z+Omega)^{-1} */
	    for (j = 0; j < ncisqr; j++) {
		for (k = 0; k < nlev; k++) mm[j] += bVi[j + k*ncisqr];
	    }
	    mm += ncisqr;	/* Extra term for \vb */
	    for (j = 0; j < p; j++) {
		F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one,
				RZXi + j * dims[0], &nci,
				&one, mm, &nci);
	    }
	}
	Free(tmp);
    }
    return val;
}

/**
 * Return a length nf list of arrays of dimension (nci, nci, 4).  The
 * values of these arrays are assigned in lmer_firstDer.
 *
 * @param nf number of factors
 * @param nc vector of number of columns per factor
 *
 * @return pointer to a list of REAL arrays
 */
static
SEXP EM_grad_array(int nf, const int nc[])
{
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    int i;

    for (i = 0; i < nf; i++) {
	SET_VECTOR_ELT(val, i, alloc3Darray(REALSXP, nc[i], nc[i], 4));
    }
    UNPROTECT(1);
    return val;
}

/**
 * Fill in the 4-dimensional vector of linear combinations of the
 * firstDer array according to whether ECME steps or the gradient are
 * needed and to whether or not REML is being used.
 *
 * @param cc coefficient vector to be filled in
 * @param EM non-zero for ECME steps, zero for gradient
 * @param REML non-zero for REML, zero for ML
 * @param ns ns[0] is p+1, ns[1] is n
 *
 * @return cc with the coefficients filled in
 */
static R_INLINE
double *EM_grad_lc(double *cc, int EM, int REML, int ns[])
{
    cc[0] = EM ? 0. : -1.;
    cc[1] = (double)(ns[1] - (REML ? ns[0] - 1 : 0));
    cc[2] = 1.;
    cc[3] = REML ? 1. : 0.;
    return cc;
}


/**
 * Print the verbose output in the ECME iterations
 *
 * @param x pointer to an ssclme object
 * @param iter iteration number
 * @param REML non-zero for REML, zero for ML
 * @param firstDer arrays for calculating ECME steps and the first derivative
 * @param val Pointer to a list of arrays to receive the calculated values
 */
static
void EMsteps_verbose_print(SEXP x, int iter, int REML, SEXP firstDer, SEXP val)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	pMat = VECTOR_ELT(val, 2);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Its = INTEGER(VECTOR_ELT(val, 0)),
	i, ifour = 4, ii, ione = 1, jj, nf = length(Omega),
	niter = INTEGER(getAttrib(pMat, R_DimSymbol))[0];
    double
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	*cc = EM_grad_lc(Calloc(4, double), 0, REML, nc + nf),
	*Devs = REAL(VECTOR_ELT(val, 1)),
	*pars = REAL(pMat) + iter,
	*grds = REAL(VECTOR_ELT(val, 3)) + iter,
	one = 1., zero = 0.;

    lmer_factor(x);
    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", Its[iter] = iter, Devs[iter] = dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *Grad = Calloc(ncisqr, double);

				/* diagonals */
	Rprintf(" (%#8g", *pars = Omgi[0]);
	pars += niter;
	for (jj = 1; jj < nci; jj++, pars += niter) {
	    Rprintf(" %#8g", *pars = Omgi[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++, pars += niter)
		Rprintf(" %#8g", *pars = Omgi[ii + jj * nci]);
				/* Evaluate and print the gradient */
	F77_CALL(dgemv)("N", &ncisqr, &ifour, &one,
			REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			cc, &ione, &zero, Grad, &ione);
	Rprintf(":%#8.3g", *grds = Grad[0]);
	grds += niter;
				/* diagonals */
	for (jj = 1; jj < nci; jj++, grds += niter) {
	    Rprintf(" %#8.3g", *grds = Grad[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++, grds += niter)
		Rprintf(" %#8.3g", *grds = Grad[ii + jj * nci]);
	Rprintf(")");
	Free(Grad);
    }
    Rprintf("\n");
    Free(cc);
}

/**
 * Perform ECME steps for the REML or ML criterion.
 *
 * @param x pointer to an ssclme object
 * @param nsteps pointer to an integer scalar - the number of ECME steps to perform
 * @param Verbp pointer to a logical scalar indicating verbose output
 *
 * @return R_NilValue if verb == FALSE, otherwise a list of iteration
 *numbers, deviances, parameters, and gradients.
 */
SEXP lmer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	flist = GET_SLOT(x, Matrix_flistSym),
	val = R_NilValue;
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))), "REML"),
	i, ifour = 4, info, ione = 1, iter,
	nEM = asInteger(nsteps),
	nf = length(Omega),
	verb = asLogical(Verbp);
    double
	*cc = EM_grad_lc(Calloc(4, double), 1, REML, nc + nf),
	zero = 0.0;
    SEXP firstDer = PROTECT(EM_grad_array(nf, nc));

    lmer_firstDer(x, firstDer);
    if (verb) {
	int nEMp1 = nEM + 1, npar = coef_length(nf, nc);
	val = PROTECT(allocVector(VECSXP, 4));
	SET_VECTOR_ELT(val, 0, allocVector(INTSXP, nEMp1));
	SET_VECTOR_ELT(val, 1, allocVector(REALSXP, nEMp1));
	SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, nEMp1, npar));
	SET_VECTOR_ELT(val, 3, allocMatrix(REALSXP, nEMp1, npar));
	EMsteps_verbose_print(x, 0, REML, firstDer, val);
    }
    for (iter = 0; iter < nEM; iter++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], ncisqr = nci * nci;
	    double *Omgi = REAL(VECTOR_ELT(Omega, i)),
		mult = 1./
		((double) length(getAttrib(VECTOR_ELT(flist, i),
				 R_LevelsSymbol)));

	    F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &mult,
			    REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			    cc, &ifour, &zero, Omgi, &ncisqr);
	    F77_CALL(dpotrf)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error(_("DPOTRF in ECME update gave code %d"), info);
	    F77_CALL(dpotri)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error(_("Matrix inverse in ECME update gave code %d"), info);
	}
	status[0] = status[1] = 0;
	lmer_firstDer(x, firstDer);
	if (verb) EMsteps_verbose_print(x, iter + 1, REML, firstDer, val);
    }
    lmer_factor(x);
    if (verb) UNPROTECT(1);
    UNPROTECT(1);
    return val;
}

/** 
 * Evaluate the gradient vector
 * 
 * @param x Pointer to an lmer object
 * @param pType Pointer to an integer indicator of the parameterization being used
 * 
 * @return pointer to a gradient vector
 */
SEXP lmer_gradient(SEXP x, SEXP pType)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	dind, i, ifour = 4, info, ione = 1, nf = length(Omega),
	odind, ptyp = asInteger(pType);
    SEXP
	firstDer = lmer_firstDer(x, PROTECT(EM_grad_array(nf, nc))),
	val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double
	*cc = EM_grad_lc(Calloc(4, double), 0,
			 !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
				 "REML"), nc + nf),
	one = 1.0, zero = 0.0;

    dind = 0;			/* index into val for diagonals */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *tmp = Calloc(ncisqr, double);

	F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &one,
			REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			cc, &ifour, &zero, tmp, &ncisqr);
	if (nci == 1) {
	    REAL(val)[dind++] =
		(ptyp?((ptyp == 1)?Omgi[0]: -Omgi[0] * Omgi[0]) : 1) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;

	    odind = dind + nci; /* index into val for off-diagonals */
	    if (ptyp) {
		double *chol = Memcpy(Calloc(ncisqr, double),
				      REAL(VECTOR_ELT(Omega, i)), ncisqr),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

		F77_CALL(dpotrf)("U", &nci, chol, &nci, &info);
		if (info)
		    error(_("Omega[[%d]] is not positive definite"), i + 1);
		/* tmp2 := chol %*% tmp using only upper triangle of tmp */
		F77_CALL(dsymm)("R", "U", &nci, &nci, &one, tmp, &nci,
				chol, &nci, &zero, tmp2, &nci);
		/* full symmetric product gives diagonals */
		F77_CALL(dtrmm)("R", "U", "T", "N", &nci, &nci, &one, chol, &nci,
				Memcpy(tmp, tmp2, ncisqr), &nci);
		/* overwrite upper triangle with gradients for positions in L' */
		for (ii = 1; ii < nci; ii++) {
		    for (j = 0; j < ii; j++) {
			tmp[j + ii*nci] = chol[j*ncip1] * tmp2[j + ii*nci];
			tmp[ii + j*nci] = 0.;
		    }
		}
		if (ptyp > 1)
		    for (ii = 0; ii < nci; ii++) {
			int ind = ii * ncip1;
			double sqrtd = chol[ind];
			tmp[ind] *= -(sqrtd*sqrtd);
		    }
	    }
	    for (j = 0; j < nci; j++) {
		REAL(val)[dind + j] = tmp[j * ncip1];
		for (ii = 0; ii < j; ii++) /* offdiagonals count twice */
		    REAL(val)[odind++] = 2. * tmp[ii + j * nci];
	    }
	    dind = odind;
	}
	Free(tmp);
    }
    UNPROTECT(2);
    Free(cc);
    return val;
}

/**
 * Fill in five symmetric matrices, providing the
 * information to generate the Hessian.

 * @param x pointer to an lme object
 * @param Valp ignored at present
 *
 * @return Valp an array consisting of five symmetric faces
 */
static
SEXP lmer_secondDer(SEXP x, SEXP Valp)
{
    SEXP
	D = GET_SLOT(x, Matrix_DSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	levels = GET_SLOT(x, R_LevelsSymbol),
	val;
    int *dRZX = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    SEXP
	firstDer = lmer_firstDer(x, PROTECT(EM_grad_array(nf, nc)));
    double
	*RZX = REAL(RZXP),
	*b = REAL(RZXP) + dRZX[0] * p,
	*bbface,		/* vec of second faces of firstDer elts */
	one = 1.,
	zero = 0.;

    Q = 0;			/* number of rows and columns in the result */
    for (i = 0; i < nf; i++) Q += nc[i] * nc[i];
    Qsqr = Q * Q;
    bbface = Calloc(Q, double);
    val = PROTECT(alloc3Darray(REALSXP, Q, Q, 5));
    AZERO(REAL(val), Qsqr * 5);

    pos = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncisqr = nci * nci;
	double *fDi = REAL(VECTOR_ELT(firstDer, i)),
	    mult = 1./((double) length(VECTOR_ELT(levels, i)));

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of firstDer on the diagonal
	 * of the third face of val */
	F77_CALL(dsyr)("U", &ncisqr, &mult, fDi + 2 * ncisqr, &ione,
		       REAL(val) + 2 * Qsqr + pos * Q, &Q);
	pos += ncisqr;
    }
				/* fifth face of val is outer product of bbface */
    F77_CALL(dsyr)("U", &Q, &one, bbface, &ione, REAL(val) + 4 * Qsqr, &Q);
				/* fourth face from \bb\trans\der\vb\der\bb */
    AZERO(REAL(val) + 3 * Qsqr, Qsqr); /* zero accumulator */
    RZXpos = 0;
    facepos = 0;
    for (i = 0; i < nf; i++) {
	int ii, jj, nci = nc[i], ncisqr = nci * nci, nctp = nci * p,
	    nlev = length(VECTOR_ELT(levels, i));
	int maxpq = (p > nci) ? p : nci;
	double
	    *Di = REAL(VECTOR_ELT(D, i)),
	    *arr = Calloc(ncisqr * maxpq, double), /* tmp 3Darray */
	    *face = REAL(val) + 3 * Qsqr,
	    *mat = Calloc(nci * maxpq, double); /* tmp matrix */

	for (j = 0; j < nlev; j++) {
	    F77_CALL(dgemm)("T", "T", &p, &nci, &nci,
			    &one, RZX + j * nci, dRZX, Di + j * ncisqr, &nci,
			    &zero, mat, &p);
	    F77_CALL(dgemm)("N", "N", &nctp, &nci, &ione,
			    &one, mat, &nctp, b + j * nci, &ione,
			    &zero, arr, &nctp);
	    F77_CALL(dsyrk)("U", "T", &ncisqr, &p, &one, arr, &p,
			    &one, face + facepos, &Q);
				/* Add the D_{i,j}^{-T/2} term */
	    Memcpy(mat, Di + j * ncisqr, ncisqr);
	    for (jj = 1; jj < nci; jj++) { /* transpose mat */
		for (ii = 0; ii < jj; ii++) {
		    mat[jj + ii * nci] = mat[ii + jj * nci];
		    mat[ii + jj * nci] = 0.;
		}
	    }
	    F77_CALL(dgemm)("N", "N", &ncisqr, &nci, &ione,
			    &one, mat, &ncisqr, b + j * nci, &ione,
			    &zero, arr, &ncisqr);
	    /* FIXME: Next call could be dsyr (it's rank one). */
	    F77_CALL(dsyrk)("U", "T", &ncisqr, &nci, &one, arr, &nci,
			    &one, face + facepos, &Q);

	}
	RZXpos += nci * nlev;
	facepos += ncisqr;
	Free(arr); Free(mat);
    }
    UNPROTECT(2);
    Free(bbface);
    return val;
}

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static double*
internal_symmetrize(double *a, const int nc)
{
    int i, j;

    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/**
 * Return the unscaled variances
 *
 * @param x pointer to an lmer object
 *
 * @return a list similar to the Omega list with the unscaled variances
 */
SEXP lmer_variances(SEXP x)
{
    SEXP Omg = PROTECT(duplicate(GET_SLOT(x, Matrix_OmegaSym)));
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omg);

    for (i = 0; i < nf; i++) {
	double *mm = REAL(VECTOR_ELT(Omg, i));
	int j, nci = nc[i];

	F77_CALL(dpotrf)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error(_("DPOTRF returned error code %d on Omega[%d]"),
		  j, i + 1);
	F77_CALL(dpotri)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error(_("DTRTRI returned error code %d on Omega[%d]"),
		  j, i + 1);
	internal_symmetrize(mm, nci);
    }
    UNPROTECT(1);
    return Omg;
}

/** 
 * Calculate the fitted values.
 * 
 * @param x pointer to an lmer object
 * @param mmats list of model matrices
 * @param useRand logical scalar indicating if the random
 * effects should be used
 * @param val array to hold the fitted values
 * 
 * @return pointer to a numeric array of fitted values
 */
static double*
internal_fitted(SEXP x, SEXP mmats, int useRand, double val[])
{
    SEXP flist = GET_SLOT(x, Matrix_flistSym),
	fixef = PROTECT(lmer_fixef(x));
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)), ione = 1,
	nf = LENGTH(flist), nobs = length(VECTOR_ELT(flist, 0)),
	p = LENGTH(fixef);
    double one = 1.0, zero = 0.0;

    if (p) F77_CALL(dgemm)("N", "N", &nobs, &ione, &p, &one,
			   REAL(VECTOR_ELT(mmats, nf)), &nobs,
			   REAL(fixef), &p, &zero, val, &nobs);
    else AZERO(val, nobs);
    if (useRand) {
	int i;
	SEXP b = PROTECT(lmer_ranef(x));
	for (i = 0; i < nf; i++) {
	    SEXP bi = VECTOR_ELT(b, i);
	    int mi = INTEGER(getAttrib(bi, R_DimSymbol))[0];
	    int *ff = INTEGER(VECTOR_ELT(flist, i)), j, nci = nc[i];
	    double *mm = REAL(VECTOR_ELT(mmats, i));

	    for (j = 0; j < nobs; ) {
		int nn = 1, lev = ff[j];
		/* check for adjacent rows with same factor level */
		while ((j + nn) < nobs && ff[j + nn] == lev) nn++; 
		F77_CALL(dgemm)("N", "T", &nn, &ione, &nci,
				&one, mm + j, &nobs,
				REAL(bi) + (lev - 1), &mi,
				&one, &val[j], &nobs);
		j += nn;
	    }
	}
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return val;
}

/** 
 * Return the fitted values.
 * 
 * @param x pointer to an lmer object
 * @param mmats list of model matrices
 * @param useRf pointer to a logical scalar indicating if the random
 * effects should be used
 * 
 * @return pointer to a numeric array of fitted values
 */

SEXP lmer_fitted(SEXP x, SEXP mmats, SEXP useRf)
{
    int nobs = LENGTH(VECTOR_ELT(GET_SLOT(x, Matrix_flistSym), 0));
    SEXP ans = PROTECT(allocVector(REALSXP, nobs));

    internal_fitted(x, mmats, asLogical(useRf), REAL(ans));
    UNPROTECT(1);
    return ans;
}

/** 
 * Copy an lmer object collapsing the fixed effects slots to the response only.
 * 
 * @param x pointer to an lmer object
 * 
 * @return a duplicate of x with the fixed effects slots collapsed to the response only
 */
SEXP lmer_collapse(SEXP x)
{
    SEXP 
        ans = PROTECT(NEW_OBJECT(MAKE_CLASS("lmer"))),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	Dim = getAttrib(GET_SLOT(x, Matrix_ZtXSym), R_DimSymbol);
    int 
        nf = length(Omega), 
        nz = INTEGER(Dim)[0];

    slot_dup(ans, x, Matrix_flistSym);
    slot_dup(ans, x, Matrix_permSym);
    slot_dup(ans, x, Matrix_ParentSym);
    slot_dup(ans, x, Matrix_DSym);
    slot_dup(ans, x, Matrix_bVarSym);
    slot_dup(ans, x, Matrix_LSym);
    slot_dup(ans, x, Matrix_ZZpOSym);
    slot_dup(ans, x, Matrix_OmegaSym);
    slot_dup(ans, x, Matrix_methodSym);
    slot_dup(ans, x, Matrix_ZtZSym);

    slot_dup(ans, x, Matrix_cnamesSym);
    slot_dup(ans, x, Matrix_devCompSym);
    slot_dup(ans, x, Matrix_devianceSym);
    slot_dup(ans, x, Matrix_ncSym);
    slot_dup(ans, x, Matrix_GpSym);
    slot_dup(ans, x, Matrix_statusSym);

    INTEGER(GET_SLOT(ans, Matrix_ncSym))[nf] = 1;
    SET_SLOT(ans, Matrix_XtXSym, allocMatrix(REALSXP, 1, 1));
    REAL(GET_SLOT(ans, Matrix_XtXSym))[0] = NA_REAL;
    SET_SLOT(ans, Matrix_RXXSym, allocMatrix(REALSXP, 1, 1));
    REAL(GET_SLOT(ans, Matrix_RXXSym))[0] = NA_REAL;
    SET_SLOT(ans, Matrix_ZtXSym, allocMatrix(REALSXP, nz, 1));
    SET_SLOT(ans, Matrix_RZXSym, allocMatrix(REALSXP, nz, 1));
    LOGICAL(GET_SLOT(ans, Matrix_statusSym))[0] = 0;
    UNPROTECT(1);
    return ans;
}

/* Gauss-Hermite Quadrature x positions */
static const double
    GHQ_x1[1] = {0},
    GHQ_w1[1] = {1},
    GHQ_x2[1] = {1},
    GHQ_w2[1] = {0.5},
    GHQ_x3[2] = {1.7320507779261, 0},
    GHQ_w3[2] = {0.166666666666667, 0.666666666666667},
    GHQ_x4[2] = {2.3344141783872, 0.74196377160456},
    GHQ_w4[2] = {0.0458758533899086, 0.454124131589555},
    GHQ_x5[3] = {2.85696996497785, 1.35562615677371, 0},
    GHQ_w5[3] = {0.0112574109895360, 0.222075915334214, 0.533333317311434},
    GHQ_x6[3] = {3.32425737665988, 1.88917584542184, 0.61670657963811},
    GHQ_w6[3] = {0.00255578432527774, 0.0886157433798025, 0.408828457274383},
    GHQ_x7[4] = {3.7504396535397, 2.36675937022918, 1.15440537498316, 0},
    GHQ_w7[4] = {0.000548268839501628, 0.0307571230436095, 0.240123171391455,
		 0.457142843409801},
    GHQ_x8[4] = {4.14454711519499, 2.80248581332504, 1.63651901442728,
		 0.539079802125417},
    GHQ_w8[4] = {0.000112614534992306, 0.0096352198313359, 0.117239904139746,
		 0.373012246473389},
    GHQ_x9[5] = {4.51274578616743, 3.20542894799789, 2.07684794313409,
		 1.02325564627686, 0},
    GHQ_w9[5] = {2.23458433364535e-05, 0.00278914123744297, 0.0499164052656755,
		 0.244097495561989, 0.406349194142045},
    GHQ_x10[5] = {4.85946274516615, 3.58182342225163, 2.48432579912153,
		  1.46598906930182, 0.484935699216176},
    GHQ_w10[5] = {4.31065250122166e-06, 0.000758070911538954, 0.0191115799266379,
		  0.135483698910192, 0.344642324578594},
    GHQ_x11[6] = {5.18800113558601, 3.93616653976536, 2.86512311160915,
		  1.87603498804787, 0.928868981484148, 0},
    GHQ_w11[6] = {8.12184954622583e-07, 0.000195671924393029, 0.0067202850336527,
		  0.066138744084179, 0.242240292596812, 0.36940835831095};

static const double
    *GHQ_x[12] = {(double *) NULL, GHQ_x1, GHQ_x2, GHQ_x3, GHQ_x4, GHQ_x5,
		  GHQ_x6, GHQ_x7, GHQ_x8, GHQ_x9, GHQ_x10, GHQ_x11},
    *GHQ_w[12] = {(double *) NULL, GHQ_w1, GHQ_w2, GHQ_w3, GHQ_w4, GHQ_w5,
		  GHQ_w6, GHQ_w7, GHQ_w8, GHQ_w9, GHQ_w10, GHQ_w11};

/** 
 * Compute certain components of the Laplace likelihood approximation 
 * 
 * @param x pointer to an lmer object
 * 
 * @return log likelihood
 */
SEXP glmer_Laplace_devComp(SEXP x) {
    SEXP 
        ranef = PROTECT(lmer_ranef(x)),
        bVar = GET_SLOT(x, Matrix_bVarSym),
        Omg = GET_SLOT(x, Matrix_OmegaSym);
    int 
        *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
        i, ione = 1, nf = length(Omg);
    double ans = 0, one = 1, zero = 0;

    for (i = 0; i < nf; i++) {
        int j, k, nci = nc[i];
        int ncip1 = nci + 1, ncisqr = nci * nci,
	    nlev = (Gp[i + 1] - Gp[i])/nci;
	int ntot = nlev * nci;
	double *bVi = REAL(VECTOR_ELT(bVar, i)),
	    *tmp = Memcpy(Calloc(ncisqr, double),
			  REAL(VECTOR_ELT(Omg, i)), ncisqr),
	    *tmp2 = Calloc(ntot, double);

        F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
        if (j)
            error(_("Leading %d minor of Omega[[%d]] not positive definite"),
                  j, i + 1);
        for (j = 0; j < nci; j++) { /* 0.5 * nlev * logDet(Omega_i) */
            ans += nlev * log(tmp[j * ncip1]); /* (2 * 0.5) since factoring */
        }
	F77_CALL(dgemm)("N", "T", &nlev, &nci, &nci, &one,
			REAL(VECTOR_ELT(ranef, i)), &nlev,
			tmp, &nci, &zero, tmp2, &nlev);
        ans -= 0.5 * F77_CALL(ddot)(&ntot, tmp2, &ione, tmp2, &ione);

        for (k = 0; k < nlev; k++) {
            Memcpy(tmp, bVi + k * ncisqr, ncisqr);
            F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
            if (j)
                error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
            for (j = 0; j < nci; j++) {
                ans += log(tmp[j * ncip1]);
            }
        }
        Free(tmp); Free(tmp2);
    }
    UNPROTECT(1);
    return ScalarReal(ans);
}
				 
static void
glmer_wt_lst(SEXP MLin, double *wts, double *adjst, int n, SEXP MLout)
{
    int i, j, nf = LENGTH(MLin);
    SEXP lastM;

    for (i = 0; i < nf; i++) {
	SEXP Min = VECTOR_ELT(MLin, i),
	    Mout = VECTOR_ELT(MLout, i);
	int *din, *dout, k, nc;

	if (!(isMatrix(Min) && isReal(Min)))
	    error(_("component %d of MLin is not a numeric matrix"), i + 1);
	din = INTEGER(getAttrib(Min, R_DimSymbol));
	nc = din[1];
	if (din[0] != n)
	    error(_("component %d of MLin has %d rows, expected %d"), i + 1,
		  din[0], n);
	if (!(isMatrix(Mout) && isReal(Mout)))
	    error(_("component %d of MLout is not a numeric matrix"), i + 1);
	dout = INTEGER(getAttrib(Mout, R_DimSymbol));
	if (dout[0] != n)
	    error(_("component %d of MLout has %d rows, expected %d"), i + 1,
		  dout[0], n);
	if (dout[1] != nc)
	    error(_("component %d of MLout has %d columns, expected %d"), i + 1,
		  dout[1], nc);
	for (k = 0; k < nc; k++) {
	    for (j = 0; j < n; j++) {
		REAL(Mout)[j + k * n] = REAL(Min)[j + k * n] * wts[j];
	    }
	}
    }
    lastM = VECTOR_ELT(MLout, nf - 1);
    j = INTEGER(getAttrib(lastM, R_DimSymbol))[1] - 1;
    for (i = 0; i < n; i++)
	REAL(lastM)[j*n + i] = adjst[i] * wts[i];
}    

static
SEXP find_and_check(SEXP rho, SEXP nm, SEXPTYPE mode, int len)
{    
    SEXP ans;
    if (R_NilValue == PROTECT(ans = findVarInFrame(rho, nm)))
	error(_("environment `rho' must contain an object `%s'"),
	      CHAR(PRINTNAME(nm)));
    if (TYPEOF(ans) != mode)
	error(_("object `%s' of incorrect type"), CHAR(PRINTNAME(nm)));
    if (len && LENGTH(ans) != len)
	error(_("object `%s' must be of length `%d'"),
	      CHAR(PRINTNAME(nm)), len);
    UNPROTECT(1);
    return ans;
}

static
SEXP eval_check_store(SEXP fcn, SEXP rho, SEXP vv)
{
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != TYPEOF(vv) || LENGTH(v) != LENGTH(vv))
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), TYPEOF(vv), LENGTH(vv));
    switch (TYPEOF(v)) {
    case LGLSXP:
	Memcpy(LOGICAL(vv), LOGICAL(v), LENGTH(vv));
	break;
    case INTSXP:
	Memcpy(INTEGER(vv), INTEGER(v), LENGTH(vv));
	break;
    case REALSXP:
	Memcpy(REAL(vv), REAL(v), LENGTH(vv));
	break;
    default:
	error(_("invalid type for eval_check_store"));
    }
    UNPROTECT(1);
    return vv;
}

static
SEXP eval_check(SEXP fcn, SEXP rho, SEXPTYPE mode, int len)
{
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != mode || LENGTH(v) != len)
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), mode, len);
    UNPROTECT(1);
    return v;
}

#define BHAT_NITER 20

SEXP lmer_devLaplace(SEXP pars, SEXP tolp, SEXP rho)
{
    SEXP dev_resids, eta, linkinv, mu, mu_eta, offset, rdobj,
	unwtd, variance, wts, wtd, x, y;
    int conv, i, ii, ione = 1, n, *nc, nf, p;
    SEXP dmu_deta, var, devr;
    double *etaold, *fitted, *off, *w, *z, ans,
	one = 1, tol = asReal(tolp), zero = 0;
	
    if (!isReal(pars) || LENGTH(pars) < 1)
	error(_("`%s' must be a nonempty, numeric vector"), "pars");
    if (!isEnvironment(rho))
	error(_("`rho' must be an environment"));
    y = find_and_check(rho, install("y"), REALSXP, 0);
    if ((n = LENGTH(y)) < 1)
	error(_("`%s' must be a nonempty, numeric vector"), "y");	
    etaold = Calloc(n, double);
    fitted = Calloc(n, double);
    off = Calloc(n, double);
    w = Calloc(n, double);
    z = Calloc(n, double);
    rdobj = find_and_check(rho, install("mer"), VECSXP, 0);
    nc = INTEGER(GET_SLOT(rdobj, Matrix_ncSym));
    nf = LENGTH(GET_SLOT(rdobj, Matrix_flistSym));
    p = LENGTH(pars) - coef_length(nf, nc);
    mu = find_and_check(rho, install("mu"), REALSXP, n);
    offset = find_and_check(rho, install("offset"), REALSXP, n);
    x = find_and_check(rho, install("x"), REALSXP, n * p);
    eta = find_and_check(rho, install("eta"), REALSXP, n);
    unwtd = find_and_check(rho, install("mmo"), VECSXP, nf + 1);
    wts = find_and_check(rho, install("wts"), REALSXP, n);
    wtd = find_and_check(rho, install("mmats"), VECSXP, nf + 1);

    dev_resids = find_and_check(rho, install("dev.resids"), LANGSXP, 0);
    linkinv = find_and_check(rho, install("linkinv"), LANGSXP, 0);
    mu_eta = find_and_check(rho, install("mu.eta"), LANGSXP, 0);
    variance = find_and_check(rho, install("variance"), LANGSXP, 0);

    F77_CALL(dgemv)("N", &n, &p, &one, REAL(x), &n,
		    REAL(pars), &ione, &zero, off, &ione);
    for (ii = 0; ii < n; ii++)
	REAL(eta)[ii] = (off[ii] += REAL(offset)[ii]);
    internal_coefGets(rdobj, REAL(pars) + p, 2);
    conv = 0;
    
    Memcpy(etaold, REAL(eta), n);
    for (i = 0; i < BHAT_NITER && !conv; i++) {
	double max_eta, max_diff;

	eval_check_store(linkinv, rho, mu);
	dmu_deta = PROTECT(eval_check(mu_eta, rho, REALSXP, n));
	var = PROTECT(eval_check(variance, rho, REALSXP, n));
	for (ii = 0; ii < n; ii++) {
	    w[ii] = REAL(wts)[ii] * REAL(dmu_deta)[ii]/sqrt(REAL(var)[ii]);
	    z[ii] = REAL(eta)[ii] - off[ii] +
		(REAL(y)[ii] - REAL(mu)[ii])/REAL(dmu_deta)[ii];
	}
	UNPROTECT(2);
	glmer_wt_lst(unwtd, w, z, n, wtd);
	lmer_update_mm(rdobj, wtd);
	internal_fitted(rdobj, unwtd, 1, fitted);
	max_eta = max_diff = -1.;
	for (ii = 0; ii < n; ii++) {
	    double abs_eta, abs_diff;
	    
	    REAL(eta)[ii] = off[ii] + fitted[ii];
	    abs_eta = fabs(REAL(eta)[ii]);
	    if (abs_eta > max_eta) max_eta = abs_eta;
	    abs_diff = fabs(REAL(eta)[ii] - etaold[ii]);
	    if (abs_diff > max_diff) max_diff = abs_diff;
	    etaold[ii] = REAL(eta)[ii];
	}
	if (max_diff < (0.1 + max_eta) * tol) {
	    conv = 1;
	    break;
	}
    }
    if (!conv) warning(_("iterations for bhat did not converge"));
    Free(etaold); Free(fitted); Free(off); Free(w); Free(z);
    devr = PROTECT(eval_check(dev_resids, rho, REALSXP, n));
    ans = 0;
    for (ii = 0; ii < n; ii++) ans += REAL(devr)[ii];
    UNPROTECT(1);
    return ScalarReal(ans - 2 * asReal(glmer_Laplace_devComp(rdobj)));
}

/**
 * Produce a weighted copy of the matrices in MLin in the storage
 * allocated to MLout
 *
 * @param MLin input matrix list
 * @param wts real vector of weights
 * @param adjst adjusted response
 * @param MLout On input a list of matrices of the same dimensions as MLin.
 *
 * @return MLout with its contents overwritten by a weighted copy of
 * MLin according to wts with adjst overwriting the response.
 */
SEXP glmer_weight_matrix_list(SEXP MLin, SEXP wts, SEXP adjst, SEXP MLout)
{
    int n, nf;
    
    if (!(isNewList(MLin) && isReal(wts) && isReal(adjst) && isNewList(MLout)))
	error(_("Incorrect argument type"));
    nf = LENGTH(MLin);
    if (LENGTH(MLout) != nf)
	error(_("Lengths of MLin (%d) and MLout (%d) must match"), nf,
	      LENGTH(MLout));
    n = LENGTH(wts);
    if (LENGTH(adjst) != n)
	error(_("Expected adjst to have length %d, got %d"), n, LENGTH(adjst));
    glmer_wt_lst(MLin, REAL(wts), REAL(adjst), n, MLout);
    return MLout;
}
