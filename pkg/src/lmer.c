#include "lmer.h"
/* TODO 
 * - Remove the Linv slot.
 * - Change the structure so that ZZpO is diagonal and the off-diagonals of L are used
 * - The egsingle example with ~year|childid+schoolid shows an unusual
 * drop in the deviance when switching from ECME to optim.  Is it real?
 * (Apparently so.) 
 */

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
 * Calculate the zero-based index in a packed lower triangular matrix.  This is
 * used for the arrays of blocked sparse matrices.
 * 
 * @param i column number (zero-based)
 * @param k row number (zero-based)
 * 
 * @return The index of the (k,i) element of a packed lower triangular matrix
 */    
static R_INLINE
int Lind(int i, int k)
{
    return (i * (i + 1))/2 + k;
}

/** 
 * Allocate a 3-dimensional array
 * 
 * @param TYP The R Type code (e.g. INTSXP)
 * @param nr number of rows
 * @param nc number of columns
 * @param nf number of faces
 * 
 * @return A 3-dimensional array of the indicated dimensions and type
 */
static
SEXP alloc3Darray(int TYP, int nr, int nc, int nf)
{
    SEXP val, dd = PROTECT(allocVector(INTSXP, 3));
    
    INTEGER(dd)[0] = nr; INTEGER(dd)[1] = nc; INTEGER(dd)[2] = nf;
    val = allocArray(TYP, dd);
    UNPROTECT(1);
    return val;
}

static R_INLINE
int match_mat_dims(const int xd[], const int yd[])
{
    return xd[0] == yd[0] && xd[1] == yd[1];
}
    
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
	return ScalarString(mkChar("Slots ZtX, XtX, RZX, and RXX must be real matrices"));
    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, R_DimSymbol))))
	return ScalarString(mkChar("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, R_DimSymbol))))
	return ScalarString(mkChar("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return ScalarString(mkChar("Slots XtX must be a square matrix with same no. of cols as ZtX"));
    return ScalarLogical(1);
}

static R_INLINE
int Tind(const int rowind[], const int colptr[], int i, int j)
{
    int k, k2 = colptr[j + 1];
    for (k = colptr[j]; k < k2; k++)
	if (rowind[k] == i) return k;
    error("row %d and column %d not defined in rowind and colptr",
	  i, j);
    return -1;			/* to keep -Wall happy */
}

/** 
 * Create the pairwise crosstabulation of the elements of flist.
 * 
 * @param flist pointer to the factor list.
 * @param nobs number of observations.
 * @param nc number of columns in the model matrices.
 * 
 * @return the pairwise crosstabulation in the form of the ZtZ array.
 */
static SEXP
lmer_crosstab(SEXP flist, int nobs, const int nc[])
{
    int i, nf = length(flist);
    int npairs = (nf * (nf + 1))/2;
    SEXP val = PROTECT(allocVector(VECSXP, npairs));
    SEXP cscbCl = MAKE_CLASS("cscBlocked");
    int *Ti = Calloc(nobs, int),
	*nlevs = Calloc(nf, int),
	**zb = Calloc(nf, int*); /* zero-based indices */
    double *ones = Calloc(nobs, double),
	*Tx = Calloc(nobs, double);
    
    for (i = 0; i < nobs; i++) ones[i] = 1.;
    for (i = 0; i < nf; i++) {	/* populate the zb vectors */
	SEXP fi = VECTOR_ELT(flist, i);
	int j;

	zb[i] = Calloc(nobs, int);
	nlevs[i] = length(getAttrib(fi, R_LevelsSymbol));
	for (j = 0; j < nobs; j++) zb[i][j] = INTEGER(fi)[j] - 1; 
	for (j = 0; j <= i; j++) {
	    int *ijp, ind = Lind(i, j), nnz;
	    SEXP ZZij;
	    
	    SET_VECTOR_ELT(val, ind, NEW_OBJECT(cscbCl));
	    ZZij = VECTOR_ELT(val, ind);
	    SET_SLOT(ZZij, Matrix_pSym, allocVector(INTSXP, nlevs[j] + 1));
	    ijp = INTEGER(GET_SLOT(ZZij, Matrix_pSym));
	    triplet_to_col(nlevs[i], nlevs[j], nobs, zb[i], zb[j], ones,
			   ijp, Ti, Tx);
	    nnz = ijp[nlevs[j]];
	    SET_SLOT(ZZij, Matrix_iSym, allocVector(INTSXP, nnz));
	    Memcpy(INTEGER(GET_SLOT(ZZij, Matrix_iSym)), Ti, nnz);
	    SET_SLOT(ZZij, Matrix_xSym, alloc3Darray(REALSXP, nc[i], nc[j], nnz));
	    /* The crosstab counts are copied into the first nnz elements */
	    /* of the x slot.  These aren't the correct array positions */
	    /* unless nc[i] == nc[j] == 1 but we don't use them. */
	    Memcpy(REAL(GET_SLOT(ZZij, Matrix_xSym)), Tx, nnz);
	}
    }
    
    for (i = 0; i < nf; i++) Free(zb[i]);
    Free(zb); Free(nlevs); Free(ones); Free(Ti); Free(Tx);
    UNPROTECT(1);
    return val;
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

    if (!isNewList(mmats) || length(mmats) != nfp1)
	error("mmats must be a list of %d model matrices", nfp1);
    for (i = 0; i <= nf; i++) {
	SEXP mmat = VECTOR_ELT(mmats, i);
	int *mdims = INTEGER(getAttrib(mmat, R_DimSymbol));
	
	if (!isMatrix(mmat) || !isReal(mmat))
	    error("element %d of mmats is not a numeric matrix", i + 1);
	if (nobs != mdims[0])
	    error("Expected %d rows in the %d'th model matrix. Got %d",
		  nobs, i+1, mdims[0]);
	if (nc[i] != mdims[1])
	    error("Expected %d columns in the %d'th model matrix. Got %d",
		  nc[i], i+1, mdims[1]);
    }
				/* Create XtX */
    X = REAL(VECTOR_ELT(mmats, nf));
    F77_CALL(dsyrk)("U", "T", &pp1, &nobs, &one, X, &nobs, &zero, XtX, nc + nf);
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
				ZZx + Tind(rowind, colptr, fac[j] - 1, f2[j] - 1)
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
 * Create an lmer object from a list grouping factors and a list of model
 * matrices.  There is one more model matrix than grouping factor.  The last
 * model matrix is the fixed effects and the response.
 * 
 * @param facs pointer to a list of grouping factors
 * @param ncv pointer to a list of model matrices
 * 
 * @return pointer to an lmer object
 */
SEXP lmer_create(SEXP flist, SEXP mmats)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("lmer")));
    SEXP ZtZ, cnames, fnms, nms;
    int *nc, i, nf = length(flist), nobs;
    
				/* Check validity of flist */
    if (!(nf > 0 && isNewList(flist)))
	error("flist must be a non-empty list");
    nobs = length(VECTOR_ELT(flist, 0));
    if (nobs < 1) error("flist[[0]] must be a non-null factor");
    for (i = 0; i < nf; i++) {
	SEXP fi = VECTOR_ELT(flist, i);
	if (!(isFactor(fi) && length(fi) == nobs))
	    error("flist[[%d]] must be a factor of length %d",
		  i + 1, nobs);
    }
    SET_SLOT(val, Matrix_flistSym, flist);
    if (!(isNewList(mmats) && length(mmats) == (nf + 1)))
	error("mmats must be a list of length %d", nf + 1);
    SET_SLOT(val, Matrix_ncSym, allocVector(INTSXP, nf + 2));
    nc = INTEGER(GET_SLOT(val, Matrix_ncSym));
    nc[nf + 1] = nobs;
    for (i = 0; i <= nf; i++) {
	SEXP mi = VECTOR_ELT(mmats, i);
	int *dims;

	if (!(isMatrix(mi) && isReal(mi)))
	    error("mmats[[%d]] must be a numeric matrix", i + 1);
	dims = INTEGER(getAttrib(mi, R_DimSymbol));
	if (dims[0] != nobs)
	    error("mmats[[%d]] must have %d rows", i + 1, nobs);
	if (dims[1] < 1)
	    error("mmats[[%d]] must have at least 1 column", i + 1);
	nc[i] = dims[1];
    }   /* Arguments have now been checked for type, dimension, etc. */
				/* Create pairwise crosstabulation in ZtZ */
    SET_SLOT(val, Matrix_ZtZSym, lmer_crosstab(flist, nobs, nc));
    lmer_populate(val);
    ZtZ = GET_SLOT(val, Matrix_ZtZSym);
    /* FIXME: Check for possible reordering of the factors to maximize the
     * number of levels (columns?) in the leading nested sequence. */
    fnms = getAttrib(flist, R_NamesSymbol);
				/* Allocate and populate nc and cnames */
    SET_SLOT(val, Matrix_cnamesSym, allocVector(VECSXP, nf + 1));
    cnames = GET_SLOT(val, Matrix_cnamesSym);
    setAttrib(cnames, R_NamesSymbol, allocVector(STRSXP, nf + 1));
    nms = getAttrib(cnames, R_NamesSymbol);
    for (i = 0; i <= nf; i++) {
	SEXP mi = VECTOR_ELT(mmats, i);
	SET_VECTOR_ELT(cnames, i,
		       duplicate(VECTOR_ELT(getAttrib(mi, R_DimNamesSymbol),
					    1)));
	if (i < nf)
	    SET_STRING_ELT(nms, i, duplicate(STRING_ELT(fnms, i)));
	else
	    SET_STRING_ELT(nms, nf, mkChar(".fixed"));
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
	ZtZ = GET_SLOT(x, Matrix_ZtZSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, k, nf = length(Omg);
    double *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym));

    for (i = 0; i < nf; i++) {
	int j, nci = nc[i], ncisqr = nci * nci;
	int nlev = (Gp[i + 1] - Gp[i])/nci;
	double *Omega = REAL(VECTOR_ELT(Omg, i));
	double *tmp = Memcpy(Calloc(ncisqr, double), Omega, ncisqr);

	F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j); /* update dcmp[1] */
	if (j)
	    error("Leading %d minor of Omega[[%d]] not positive definite",
		  j, i + 1);
	for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
	    dcmp[1] += nlev * 2. * log(tmp[j * (nci + 1)]);
	}
	Free(tmp);
	for (k = i; k < nf; k++) {
	    int ind = Lind(k, i);
	    SEXP ZZOel = VECTOR_ELT(ZZpO, ind);
	    SEXP ZZel = VECTOR_ELT(ZtZ, ind);
	    SEXP ZZm = GET_SLOT(ZZel, Matrix_xSym),
		ZZOm = GET_SLOT(ZZOel, Matrix_xSym);
	    int *Di = INTEGER(GET_SLOT(ZZOel, Matrix_iSym)),
		*Dp = INTEGER(GET_SLOT(ZZOel, Matrix_pSym)),
		*Si = INTEGER(GET_SLOT(ZZel, Matrix_iSym)),
		*Sp = INTEGER(GET_SLOT(ZZel, Matrix_pSym)),
		*dims = INTEGER(getAttrib(ZZm, R_DimSymbol));
	    int ii, jj, nblk = dims[0] * dims[1];
	    double *ZZ = REAL(ZZm), *ZZO = REAL(ZZOm);

				/* zero the whole of the ZZpO component */
	    AZERO(ZZO, INTEGER(getAttrib(ZZOm, R_DimSymbol))[2]);
	    for (j = 0; j < nlev; j++) { /* copy src blocks to dest */
		int kk, k2 = Sp[j + 1];
		for (kk = Sp[j]; kk < k2; kk++) {
		    Memcpy(ZZO + Tind(Di, Dp, Si[kk], j) * nblk,
			   ZZ + kk * nblk, nblk);
		}
	    }
	    
	    if (k == i) { /* inflate diagonal blocks */
		for (j = 0; j < nlev; j++) { 
		    double *ZZOkk = ZZO + Tind(Di, Dp, j, j) * nblk;
		    for (jj = 0; jj < nci; jj++) { 
			for (ii = 0; ii <= jj; ii++) {
			    int ind = ii + jj * nci;
			    ZZOkk[ind] += Omega[ind];
			}
		    }
		}
	    }
	}
    }
    return R_NilValue;
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
	int nml = nc[nf + 1], nreml = nml + 1 - nc[nf];
	double
	    *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXP),
	    *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	    *deviance = REAL(GET_SLOT(x, Matrix_devianceSym)),
	    minus1 = -1., one = 1.;


	dcmp[0] = dcmp[1] = dcmp[2] = dcmp[3] = 0.;
	Memcpy(RZX, REAL(GET_SLOT(x, Matrix_ZtXSym)), dims[0] * dims[1]);
	Memcpy(RXX, REAL(GET_SLOT(x, Matrix_XtXSym)), dims[1] * dims[1]);
	lmer_inflate(x);	/* initialize ZZpO */
	for (i = 0; i < nf; i++) {
	    int dind = Lind(i, i);
	    SEXP ZZOiP = VECTOR_ELT(ZZOP, dind);
	    SEXP DiP = VECTOR_ELT(DP, i);
	    SEXP LiP = VECTOR_ELT(LP, dind);
	    int nlev = INTEGER(getAttrib(DiP, R_DimSymbol))[2];
	    int jj, nci = nc[i], ncisqr = nci * nci;
	    int *Pari = INTEGER(VECTOR_ELT(Parent, i));
	    double *D = REAL(DiP);
	
	    jj = cscb_ldl(ZZOiP, Pari, LiP, DiP);
	    if (jj != nlev) error("cscb_ldl returned %d < nlev = %d", jj, nlev);
	    for (j = 0; j < nlev; j++) { /* accumulate dcmp[0] */
		double *Dj = D + j * ncisqr;
		for (jj = 0; jj < nci; jj++) /* accumulate determinant */
		    dcmp[0] += 2. * log(Dj[jj * (nci + 1)]);
	    }
	    /* Solve L_{i,i} %*% RZX_i := RZX_i */
	    cscb_trsm('L', 'N', 'U', 1., LiP, RZX + Gp[i],
		      Gp[i+1] - Gp[i], dims[1], dims[0]);
	    /* Solve D_i^{1/2} %*% RZX_i := RZX_i */
	    for (jj = 0; jj < nlev; jj++) {
		F77_CALL(dtrsm)("L", "U", "T", "N", &nci, dims + 1,
				&one, D + jj * ncisqr, &nci,
				RZX + Gp[i] + jj * nci, dims);
	    }
	    for (j = i + 1; j < nf; j++) { /*  further blocks */
		SEXP ZZOji = VECTOR_ELT(ZZOP, Lind(j, i));
		SEXP ZZOx = GET_SLOT(ZZOji, Matrix_xSym);
		double *ZZO = REAL(ZZOx);
		int *xdims = INTEGER(getAttrib(ZZOx, R_DimSymbol)),
		    *ZZp = INTEGER(GET_SLOT(ZZOji, Matrix_pSym));
		int ntot = xdims[0] * xdims[1];

		/* ZZpO_{j,i} := ZZpO_{j,i} %*% L_{i,i}^{-T} %*% D_i^{-1/2} */
		/* FIXME: Change off-diagonal blocks to L, not ZZO */
		cscb_trcbsm('R', 'L', 'T', 'U', 1.0, LiP, Pari, ZZOji);
		for (jj = 0; jj < nlev; jj++) {
		    int k, k2 = ZZp[jj + 1];
		    for (k = ZZp[jj]; k < k2; k++)
			F77_CALL(dtrsm)("R", "U", "N", "N", xdims, xdims + 1,
					&one, D + jj * ncisqr, &nci,
					ZZO + k * ntot, xdims);
		}
		/* RZX_j := RZX_j - ZZpO_{j,i} %*% RZX_i */
		cscb_mm('L', 'N', Gp[j + 1] - Gp[j], dims[1], Gp[i+1] - Gp[i],
			-1.0, ZZOji, RZX + Gp[i], dims[0],
			1.0, RZX + Gp[j], dims[0]);
	    }
	    for (j = i + 1; j < nf; j++) { /* block pairs and final update */
		SEXP ZZOji = VECTOR_ELT(ZZOP, Lind(j, i));
		SEXP ZZOx = GET_SLOT(ZZOji, Matrix_xSym);
		double *ZZO = REAL(ZZOx);
		int *xdims = INTEGER(getAttrib(ZZOx, R_DimSymbol)),
		    *ZZp = INTEGER(GET_SLOT(ZZOji, Matrix_pSym));
		int ntot = xdims[0] * xdims[1];

		
		/* ZZpO_{j,j} := ZZpO_{j,j} - ZZpO{j,i}%*%ZZpO_{j,i}^T */
		cscb_syrk('U', 'N', -1.0, ZZOji,
			  1.0, VECTOR_ELT(ZZOP, Lind(j, j)));
		for (jj = j+1; jj < nf; jj++) {
		    /* ZZpO_{jj,j} := ZZpO_{jj,j} - ZZpO{jj,i}%*%ZZpO_{j,i}^T */
		    cscb_cscbm('N', 'T', -1.0, VECTOR_ELT(ZZOP, Lind(jj, i)),
			    ZZOji, 1.0, VECTOR_ELT(ZZOP, Lind(jj, j)));
		}
		/* ZZpO_{j,i} := ZZpO_{j,i} %*% D_i^{-T/2} */
		for (jj = 0; jj < nlev; jj++) {
		    int k, k2 = ZZp[jj + 1];
		    for (k = ZZp[jj]; k < k2; k++)
			F77_CALL(dtrsm)("R", "U", "T", "N", xdims, xdims + 1,
					&one, D + jj * ncisqr, &nci,
					ZZO + k * ntot, xdims);
		}
	    }
	}
				/* downdate and factor XtX */
	F77_CALL(dsyrk)("U", "T", dims + 1, dims,
			&minus1, RZX, dims, &one, RXX, dims + 1);
	F77_CALL(dpotrf)("U", dims + 1, RXX, dims + 1, &j);
	if (j) {
	    warning("Leading minor of size %d of downdated X'X is indefinite",
		    j);
	    dcmp[2] = dcmp[3] = deviance[0] = deviance[1] = NA_REAL;
	} else {
	    for (j = 0; j < (dims[1] - 1); j++) /* 2 logDet(RXX) */
		dcmp[2] += 2 * log(RXX[j * (dims[1] + 1)]);
	    dcmp[3] = 2. * log(RXX[dims[1] * dims[1] - 1]); /* 2 log(ryy) */
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
 * @param side 'L' for left, 'R' for right
 * @param trans 'T' for transpose, otherwise no transpose
 * @param nf number of grouping factors
 * @param Gp group pointers for the rows
 * @param n number of columns
 * @param alpha multiplier
 * @param ZZpO pointer to the Z'Z+Omega sparse blocked matrix
 * @param L pointer to the L cscb object
 * @param mm pointer to the matrix of right-hand sides
 */
static void
lmer_sm(char side, char trans, int nf, const int Gp[], int n,
	double alpha, SEXP ZZpO, SEXP L, double B[], int ldb)
{
    int itr = (trans == 'T' || trans == 't'), j, k,
	lside = (side == 'L' || side == 'l');
    
    if (lside) {
	if (itr) {
	    for (j = nf - 1; j >= 0; j--) {
		int nrj = Gp[j + 1] - Gp[j];
		
		cscb_trsm('L', 'T', 'U', alpha, VECTOR_ELT(L, Lind(j, j)),
			   B + Gp[j], nrj, n, ldb);
		for (k = 0; k < j; k++) {
		    cscb_mm('L', 'T', Gp[k + 1] - Gp[k], n, nrj,
			    -1., VECTOR_ELT(ZZpO, Lind(j,k)),
			    B + Gp[j], ldb, alpha, B + Gp[k], ldb);
		}
	    }
	} else error("Code for non-transpose case not yet written");
    } else error("Code for right-side solutions not yet written");
}

/** 
 * Return the number of non-zero blocks below the diagonal in a jth column
 * block of a diagonal outer block of L^{-1}.
 * 
 * @param j column block index
 * @param Pari parent array
 * @param L pointer to the list of cscBlocked objects composing L
 * @param nnz vector to be filled in
 */
static R_INLINE
int nnzcb_diag(int j, const int Pari[])
{
    int k, val;
    for (k = Pari[j], val = 0; k > 0; k = Pari[k]) val++;
    return val;
}

static
int nnzcb_offdiag(int j, int nrb, const int Lp[], const int Li[], const int Park[])
{
    int k, k2 = Lp[j + 1], val;
    char *ind = memset(Calloc(nrb, char), 0, nrb);

    for (k = Park[j], val = 0; k > 0; k = Park[k]) val++;
    return val;
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
	error("Unable to invert singular factor of downdated X'X");
    if (!status[1]) {
	SEXP DP = GET_SLOT(x, Matrix_DSym),
	    LP = GET_SLOT(x, Matrix_LSym),
	    ParP = GET_SLOT(x, Matrix_ParentSym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    ZZpOP = GET_SLOT(x, Matrix_ZZpOSym),
	    bVarP = GET_SLOT(x, Matrix_bVarSym);
	int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	    *dims = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	    *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	    i, nn, nf = length(DP);
	double *RXX = REAL(GET_SLOT(x, Matrix_RXXSym)),
	    *RZX = REAL(RZXP),
	     minus1 = -1., one = 1., zero = 0.;

	/* RXX := RXX^{-1} */
	F77_CALL(dtrtri)("U", "N", dims + 1, RXX, dims + 1, &i);
	if (i)
	    error("Leading minor of size %d of downdated X'X,is indefinite",
		    i + 1);
	/* RZX := - RZX %*% RXX */
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims + 1, &minus1,
			RXX, dims + 1, RZX, dims);
	for(i = 0; i < nf; i++) {
	    int info, j, jj, nci = nc[i];
	    int ncisqr = nci * nci, nlev = (Gp[i+1] - Gp[i])/nci;
	    double *Di = REAL(VECTOR_ELT(DP, i)),
		*RZXi = RZX + Gp[i];
	    
	    /* D_i := D_i^{-1}; RZX_i := D_i %*% RZX_i */
	    if (nci == 1) {
		for (j = 0; j < nlev; j++) {
		    Di[j] = 1./Di[j];
		    for (jj = 0; jj < dims[1]; jj++)
			RZXi[j + jj * dims[0]] *= Di[j];
		}
	    } else {
		for (j = 0; j < nlev; j++) {
		    F77_CALL(dtrtri)("U", "N", &nci, Di + j * ncisqr, &nci, &info);
		    if (info)
			error("D[,,%d] for factor %d is singular", j + 1, i + 1);
		    F77_CALL(dtrmm)("L", "U", "N", "N", &nci, dims + 1, &one,
				    Di + j * ncisqr, &nci, RZXi + j * nci, dims);
		}
	    }
	}
	/* RZX := L^{-T} %*% RZX */
	lmer_sm('L', 'T', nf, Gp, dims[1], 1.0, GET_SLOT(x, Matrix_ZZpOSym),
		GET_SLOT(x, Matrix_LSym), RZX, dims[0]);

	/* Create bVar arrays as cprod of column blocks of D^{-T/2}%*%L^{-1} */
	for (i = 0; i < nf; i++) {
	    SEXP Lii = VECTOR_ELT(LP, Lind(i, i));
	    SEXP LxP = GET_SLOT(Lii, Matrix_xSym);
	    int *Pari = INTEGER(VECTOR_ELT(ParP, i)),
		*xdims = INTEGER(getAttrib(LxP, R_DimSymbol)),
		*Lp = INTEGER(GET_SLOT(Lii, Matrix_pSym)),
		*Li = INTEGER(GET_SLOT(Lii, Matrix_iSym)),
		ip1 = i + 1, j, k, nci = nc[i];
	    int ncisqr = nci * nci, nlev = (Gp[i+1] - Gp[i])/nci;
	    double *Di = REAL(VECTOR_ELT(DP, i)), *bVi = REAL(VECTOR_ELT(bVarP, i)),
		*Lx = REAL(GET_SLOT(Lii, Matrix_xSym));

	    for (j = 0; j < nlev; j++) {
		int nzj = nnzcb_diag(j, Pari);
		double *bVij = bVi + j * ncisqr, *Dij = Di + j * ncisqr;

		AZERO(bVij, ncisqr);
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &one, Dij, &nci, &zero, bVij, &nci);
		if (nzj) {
		    int *ind = Calloc(nzj, int), kk;
		    double *tmp = Calloc(nzj * ncisqr, double);
		    for (k = Pari[j], kk = 0; kk < nzj; kk++, k = Pari[k]) {
			int ix = Tind(Li, Lp, k, i), jj;
			if (ix < 0)
			    error("%d, %d block not found in L[%d,%d]",
				  k+1, ip1, ip1, ip1);
			ind[kk] = k;
			/* - sign in sol'n incorporated in dtrsm call below */
			Memcpy(tmp + kk * ncisqr, Lx + ix * ncisqr, ncisqr);
			for (jj = 0; jj < kk; jj++) {
			    int ijj = ind[jj];
			    int jx = Tind(Li, Lp, k, ijj);
			    if (jx < 0)
				error("%d, %d block not found in L[%d,%d]",
				      k+1, ijj+1, ip1, ip1);
			    F77_CALL(dgemm)("N", "N", &nci, &nci, &nci,
					    &minus1, tmp + jj * ncisqr, &nci,
					    Lx + jx * ncisqr, &nci,
					    &one, tmp + kk * ncisqr, &nci);
			}
		    }
		    for (k = 0; k < nzj; k++) {
			F77_CALL(dtrsm)("L", "U", "T", "U", &nci, &nci,
					&minus1, Di + ind[k] * ncisqr, &nci,
					tmp + k * ncisqr, &nci);
		    }
		    kk = nci * nzj;
		    F77_CALL(dsyrk)("U", "T", &kk, &nci, &one, tmp, &kk,
				    &one, bVij, &nci);
		    Free(tmp); Free(ind);
		}
	    }
	    
/* 	    for (k = i + 1; k < nn; k++) tpt[k] = Calloc(nc[k] * nci, double); */
/* 	    for (k = i + 1; k < nn; k++) { */
/* 		SEXP ZZki = VECTOR_ELT(ZZpOP, Lind(k, i)); */
/* 		double *Dk = REAL(VECTOR_ELT(DP, k)); */
/* 		int *rowind = INTEGER(GET_SLOT(ZZki, Matrix_iSym)); */
/* 		int nck = nc[k]; */
/* 		int sz = nci * nck; */
/* 		double *ZZkix = REAL(GET_SLOT(ZZki, Matrix_xSym)); */
		
/* 		for (j = 0; j < nlev; j++) { */
/* 		    Memcpy(tpt[k], ZZkix + j * sz, sz); */
/* 		    F77_CALL(dtrmm)("L", "U", "T", "N", &nck, &nci, */
/* 				    &one, Dk + rowind[j] * sz, &nck, tpt[k], &nck); */
/* 		    F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tpt[k], &nci, */
/* 				    &one, bVi + j * ncisqr, &nci); */
/* 		} */
/* 	    } */
/* 	    for (k = i + 1; k < nn; k++) Free(tpt[k]); */
/* 	    Free(tpt); */
/* 	    } */
	}
	status[1] = 1;
    }
    return R_NilValue;
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
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym);
    int pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1],
	nobs = INTEGER(GET_SLOT(x, Matrix_ncSym))
	[length(GET_SLOT(x, Matrix_OmegaSym)) + 1];

    lmer_invert(x);
    return ScalarReal(1./(REAL(RXXsl)[pp1*pp1 - 1] *
			  sqrt((double)(asLogical(REML) ?
					nobs + 1 - pp1 : nobs))));
}

/** 
 * Extract the upper triangles of the Omega matrices.  These aren't
 * "coefficients" but the extractor is called coef for historical
 * reasons.  Within each group these values are in the order of the
 * diagonal entries first then the strict upper triangle in row
 * order.
 * 
 * @param x pointer to an lme object
 * @param Unc pointer to a logical scalar indicating if the parameters
 * are in the unconstrained form.
 * 
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP lmer_coef(SEXP x, SEXP Unc)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), unc = asLogical(Unc), vind;
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double *vv = REAL(val);

    vind = 0;			/* index in vv */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    vv[vind++] = (unc ?
			  log(REAL(VECTOR_ELT(Omega, i))[0]) :
			  REAL(VECTOR_ELT(Omega, i))[0]);
	} else {
	    if (unc) {		/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
		double *tmp = Memcpy(Calloc(ncisq, double), 
				     REAL(VECTOR_ELT(Omega, i)), ncisq);
		F77_CALL(dpotrf)("U", &nci, tmp, &nci, &j);
		if (j)		/* should never happen */
		    error("DPOTRF returned error code %d on Omega[[%d]]",
			  j, i+1);
		for (j = 0; j < nci; j++) {
		    double diagj = tmp[j * ncip1];
		    vv[vind++] = 2. * log(diagj);
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


/** 
 * Assign the upper triangles of the Omega matrices.
 * (Called coef for historical reasons.)
 * 
 * @param x pointer to an lme object
 * @param coef pointer to an numeric vector of appropriate length
 * @param Unc pointer to a logical scalar indicating if the parameters
 * are in the unconstrained form.
 * 
 * @return R_NilValue
 */
SEXP lmer_coefGets(SEXP x, SEXP coef, SEXP Unc)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	cind, i, nf = length(Omega),
	unc = asLogical(Unc);
    double *cc = REAL(coef);

    if (length(coef) != coef_length(nf, nc) || !isReal(coef))
	error("coef must be a numeric vector of length %d",
	      coef_length(nf, nc));
    cind = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	if (nci == 1) {
	    REAL(VECTOR_ELT(Omega, i))[0] = (unc ?
					     exp(cc[cind++]) :
					     cc[cind++]);
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		j, k,
		ncip1 = nci + 1,
		ncisq = nci * nci;
	    double
		*omgi = REAL(VECTOR_ELT(Omega, i));
	    if (unc) {
		double
		    *tmp = Calloc(ncisq, double),
		    diagj, one = 1., zero = 0.;

		AZERO(omgi, ncisq);
		for (j = 0; j < nci; j++) {
		    tmp[j * ncip1] = diagj = exp(cc[cind++]/2.);
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
    return x;
}

/** 
 * Extract the conditional estimates of the fixed effects
 * 
 * @param x Pointer to an lme object
 * 
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP lmer_fixef(SEXP x)
{
    SEXP RXXsl = GET_SLOT(x, Matrix_RXXSym),
	cnames = GET_SLOT(x, Matrix_cnamesSym);
    int j, pp1 = INTEGER(getAttrib(RXXsl, R_DimSymbol))[1];
    SEXP val = PROTECT(allocVector(REALSXP, pp1));
    double
	*beta = REAL(val),
	nryyinv;		/* negative ryy-inverse */

    lmer_invert(x);
    Memcpy(beta, REAL(RXXsl) + pp1 * (pp1 - 1), pp1);
    nryyinv = -REAL(RXXsl)[pp1*pp1 - 1];
    for (j = 0; j < pp1; j++) beta[j] /= nryyinv;
    setAttrib(val, R_NamesSymbol,
	      duplicate(VECTOR_ELT(cnames, length(cnames) - 1)));
    UNPROTECT(1);
    return val;
}

/** 
 * Extract the conditional modes of the random effects.
 * 
 * @param x Pointer to an lme object
 * 
 * @return a vector containing the conditional modes of the random effects
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
static
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
		error("Omega[[%d]] is not positive definite", i + 1);
	    F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
	    if (j)
		error("Omega[[%d]] is not positive definite", i + 1);
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
static
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
			       
    /* FIXME: Check with MM for format. */
    lmer_factor(x);
    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", Its[iter] = iter, Devs[iter] = dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *Grad = Calloc(ncisqr, double);
	
				/* diagonals */
	for (jj = 0; jj < nci; jj++, pars += niter) {
	    Rprintf(" %#8g", *pars = Omgi[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++, pars += niter)
		Rprintf(" %#8g", *pars = Omgi[ii + jj * nci]);
				/* Evaluate and print the gradient */
	F77_CALL(dgemv)("N", &ncisqr, &ifour, &one,
			REAL(VECTOR_ELT(firstDer, i)), &ncisqr,
			cc, &ione, &zero, Grad, &ione);
	Rprintf(":");
				/* diagonals */
	for (jj = 0; jj < nci; jj++, grds += niter) {
	    Rprintf(" %#8g", *grds = Grad[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++, grds += niter)
		Rprintf(" %#8g", *grds = Grad[ii + jj * nci]);
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
 * @param REMLp pointer to a logical scalar indicating if REML is to be used
 * @param verb pointer to a logical scalar indicating verbose output
 * 
 * @return R_NilValue if verb == FALSE, otherwise a list of iteration
 *numbers, deviances, parameters, and gradients.
 */
SEXP lmer_ECMEsteps(SEXP x, SEXP nsteps, SEXP REMLp, SEXP Verbp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	flist = GET_SLOT(x, Matrix_flistSym),
	val = R_NilValue;
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*status = LOGICAL(GET_SLOT(x, Matrix_statusSym)),
	REML = asLogical(REMLp),
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
		error("DPOTRF in ECME update gave code %d", info);
	    F77_CALL(dpotri)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error("Matrix inverse in ECME update gave code %d", info);
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

SEXP lmer_gradient(SEXP x, SEXP REMLp, SEXP Uncp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	dind, i, ifour = 4, info, ione = 1, nf = length(Omega),
	odind, unc = asLogical(Uncp);
    SEXP
	firstDer = lmer_firstDer(x, PROTECT(EM_grad_array(nf, nc))),
	val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double
	*cc = EM_grad_lc(Calloc(4, double), 0,
			 asInteger(REMLp), nc + nf),
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
	    REAL(val)[dind++] = (unc ? Omgi[0] : 1.) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;
	    
	    odind = dind + nci; /* index into val for off-diagonals */
	    if (unc) {
		double *chol = Memcpy(Calloc(ncisqr, double),
				      REAL(VECTOR_ELT(Omega, i)), ncisqr),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

		F77_CALL(dpotrf)("U", &nci, chol, &nci, &info);
		if (info)
		    error("Omega[[%d]] is not positive definite", i + 1);
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
		Free(chol); Free(tmp2);
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
 * @param val ignored at present
 * 
 * @return val an array consisting of five symmetric faces
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
	    error("DPOTRF returned error code %d on Omega[%d]",
		  j, i + 1);
	F77_CALL(dpotri)("U", &nci, mm, &nci, &j);
	if (j)			/* shouldn't happen */
	    error("DTRTRI returned error code %d on Omega[%d]",
		  j, i + 1);
	nlme_symmetrize(mm, nci);
    }
    UNPROTECT(1);
    return Omg;
}

