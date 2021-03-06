#include "lmer.h"
#include <float.h>

/* TODO: */
/* Consider the steps in reimplementing AGQ.  First you need to find
   bhat, then evaluate the posterior variances, then step out
   according to the posterior variance, evaluate the integrand
   relative to the step. */
/* Because the Gauss-Hermite quadrature is formed as a sum, it is
 * necessary to divide the contributions to the deviance according to
 * the levels of the random effects.  This means that it is only
 * practical to use AGQ when the response vector can be split into
 * sections that are conditionally independent. As far as I can see
 * this will mean a single grouping factor only. */

/* Some of these utility functions should be moved to the R sources */

/* Internally used utilities */

/**
 * Create the RZX section of R^{-1} and store it in the RZXinv slot
 *
 * @param x - an mer object
 */
static void
internal_mer_RZXinv(SEXP x)
{
    SEXP RZXP = GET_SLOT(x, lme4_RZXSym);
    int *dims = INTEGER(GET_SLOT(RZXP, lme4_DimSym));
    int i, j, p = dims[1], q = dims[0];
    double m1[2] = {-1, 0},
	*RZXinv = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym), lme4_xSym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    CHM_DN RZX = AS_CHM_DN(RZXP), tmp1;
    int *Perm = (int *)(L->Perm), *iperm = Alloca(q, int);
    R_CheckStack();
				/* create the inverse permutation */
    for (j = 0; j < q; j++) iperm[Perm[j]] = j;
				/* solve system in L' */
    tmp1 = M_cholmod_solve(CHOLMOD_Lt, L, RZX, &c);
    /* copy columns of tmp1 to RZXinv applying the inverse permutation */
    for (j = 0; j < p; j++) {
	double *dest = RZXinv + j * q, *src = ((double*)(tmp1->x)) + j * q;
	for (i = 0; i < q; i++) dest[i] = src[iperm[i]];
    }
    M_cholmod_free_dense(&tmp1, &c);
    F77_CALL(dtrsm)("R", "U", "N", "N", &q, &p, m1,
		    REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
		    &p, RZXinv, &q);
}


/**
 * Evaluate the elements of the bVar slot (marginal
 * variance-covariance matrices of the random effects)
 *
 * @param x - an mer object
 */
static void
internal_mer_bVar(SEXP x)
{
    SEXP bVar = GET_SLOT(x, lme4_bVarSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym)), i, nf = LENGTH(bVar);
    CHM_SP rhs, sm1, sm2;
    CHM_DN dm1;
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    int *Perm = (int *)(L->Perm), *iperm = Alloca(q, int);
    R_CheckStack();

    for (i = 0; i < q; i++) iperm[Perm[i]] = i; /* inverse permutation */
    for (i = 0; i < nf; i++) {
	int j, nci = nc[i];
	int ncisqr = nci * nci, nlev = (Gp[i + 1] - Gp[i])/nci;
	double *bVi = REAL(VECTOR_ELT(bVar, i));

	rhs = M_cholmod_allocate_sparse((size_t) q, (size_t) nci,
					(size_t) nci, 1 /*sorted*/,
					1 /*packed*/, 0 /*stype*/,
					CHOLMOD_REAL, &c);
	for (j = 0; j <= nci; j++) ((int *)(rhs->p))[j] = j;
	for (j = 0; j < nci; j++) ((double *)(rhs->x))[j] = 1;
	for (j = 0; j < nlev; j++) {
	    int base = Gp[i] + j * nci, k;

	    for (k = 0; k < nci; k++) ((int*)(rhs->i))[k] = iperm[base + k];
	    sm1 = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
	    sm2 = M_cholmod_transpose(sm1, 1 /*values*/, &c);
	    M_cholmod_free_sparse(&sm1, &c);
	    sm1 = M_cholmod_aat(sm2, (int*)NULL, (size_t)0, 1 /*mode*/, &c);
	    dm1 = M_cholmod_sparse_to_dense(sm1, &c);
	    M_cholmod_free_sparse(&sm1, &c); M_cholmod_free_sparse(&sm2, &c);
	    Memcpy(bVi + j * ncisqr, (double*)(dm1->x), ncisqr);
	    M_cholmod_free_dense(&dm1, &c);
	}
    }
}

/**
 * Calculate the deviance for a linear mixed model at arbitrary
 * parameter values.  REML is not used - I don't think it applies.
 *
 * @param x a mixed-effects model representation
 * @param sigma standard deviation of per-observation noise
 * @param beta fixed-effects parameter vector
 *
 * @return deviance
 */
static
double lmm_deviance(SEXP x, double sigma, const double beta[])
{
    SEXP rXy = GET_SLOT(x, lme4_rXySym);
    int i, ione = 1, p = LENGTH(rXy);
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	*betacp = Memcpy((double*)alloca(p * sizeof(double)), beta, p),
	*rXyp = REAL(rXy),
	sprss; /* scaled penalized rss */
    R_CheckStack();

    mer_factor(x);
    F77_CALL(dtrmv)("U", "N", "N", &p,
		    REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
		    &p, betacp, &ione);
    sprss = exp(dcmp[3])/(sigma * sigma);
    for (i = 0; i < p; i++) {
	double ri = (rXyp[i] - betacp[i])/sigma;
	sprss += ri * ri;
    }

    return dcmp[4] - dcmp[5] + dcmp[0] * log(2.*PI*sigma*sigma) + sprss;
}

/**
 * Check the ZtZ matrix to see if it is a simple design from a nested
 * sequence of grouping factors.
 *
 * @param nf number of factors
 * @param nc[] number of columns per factor
 * @param Gp[] group pointers
 * @param p[] column pointers for the lower triangle of ZtZ
 *
 * @return 1 for a simple nested sequence, 0 otherwise.
 */
static int
internal_mer_isNested(int nf, const int nc[], const int Gp[], const int p[])
{
    int **cnz = Alloca(nf, int*), ans = 1, i, j, k, nct;

    for (i = 0, nct = 0; i < nf; i++) { /* total number of columns */
	nct += nc[i];
	cnz[i] = Alloca(nc[i], int);
    }
    R_CheckStack();
    for (i = 0; i < nf; i++) {	/* target number of nonzeros per column */
	for (j = 0; j < nc[i]; j++) cnz[i][j] = nct - j;
	nct -= nc[i];
    }
    for (i = 0; i < nf && ans; i++) { /* check for consistent nonzeros*/
	int nlev = (Gp[i + 1] - Gp[i])/nc[i];
	for (j = 0; j < nlev && ans; j++) {
	    for (k = 0; k < nc[i] && ans; k++) {
		int jj =  Gp[i] + j * nc[i] + k; /* column in ZtZ */
		if ((p[jj + 1] - p[jj]) != cnz[i][k]) ans = 0;
	    }
	}
    }
    return ans;
}

/**
 * Extract the coefficients
 *
 * @param x pointer to an mer object
 * @param ptyp parameter type to extract
 * @param ans vector to hold the extracted values
 *
 * @return ans
 */
static double *
internal_mer_coef(SEXP x, int ptyp, double ans[])
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, nf = length(Omega), vind;

    vind = 0;			/* index in ans */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    double dd = REAL(GET_SLOT(VECTOR_ELT(Omega, i), lme4_xSym))[0];
	    ans[vind++] = ptyp ? ((ptyp == 1) ? log(dd) : 1./dd) : dd;
	} else {
	    if (ptyp) {	/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
	        double *tmp = Memcpy(Calloc(ncisq, double),
				     REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
						   lme4_xSym)), ncisq);
		for (j = 0; j < nci; j++) {
		    double diagj = tmp[j * ncip1];
		    ans[vind++] = (ptyp == 1) ? (2. * log(diagj)) :
			1./(diagj * diagj);
		    for (k = j + 1; k < nci; k++) {
			tmp[j + k * nci] /= diagj;
		    }
		}
		for (j = 0; j < nci; j++) {
		    for (k = j + 1; k < nci; k++) {
			ans[vind++] = tmp[j + k * nci];
		    }
		}
		Free(tmp);
	    } else {		/* upper triangle of Omega[,,i] */
		int j, k, odind = vind + nci;
		double *omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), lme4_xSym));

		for (j = 0; j < nci; j++) {
		    ans[vind++] = omgi[j * ncip1];
		    for (k = j + 1; k < nci; k++) {
			ans[odind++] = omgi[k*nci + j];
		    }
		}
		vind = odind;
	    }
	}
    }
    return ans;
}

/**
 * Evaluate current estimate of sigma from an mer object
 *
 * @param x pointer to an mer object
 * @param REML indicator of whether to use REML.
 *           < 0  -> determine REML or ML from x@status
 *           == 0 -> use ML unconditionally
 *           > 0  -> use REML unconditionally
 *
 * @return
 */
static double
internal_mer_sigma(SEXP x, int REML)
{
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym));

    if (REML < 0)		/* get REML from x */
	REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1];
    mer_factor(x);
    return exp(dcmp[3]/2)/sqrt(dcmp[0] - (REML ? dcmp[1] : 0));
}

/**
 * Update the derived quantities (ZtZ, ZtX, XtX, Zty, Xty
 * and dcmp[2] = y'y when Z, X, y, wts or wrkres has been changed.
 *
 * @param x pointer to an mer object
 * @param perm permutation from the Cholesky factor
 */

static void
internal_mer_update_ZXy(SEXP x, int *perm)
{
    SEXP Xp = GET_SLOT(x, lme4_XSym), ZtZ = GET_SLOT(x, lme4_ZtZSym),
	Ztyp = GET_SLOT(x, lme4_ZtySym);
    SEXP ZtZx = GET_SLOT(ZtZ, lme4_xSym),
	ZtZp = GET_SLOT(ZtZ, lme4_pSym), ZtZi = GET_SLOT(ZtZ, lme4_iSym);
    int *dims = INTEGER(getAttrib(Xp, R_DimSymbol)), i, ione = 1, j;
    int n = dims[0], nnz, p = dims[1], q = LENGTH(Ztyp);
    CHM_SP ts1, ts2, Zt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym));
    double *XtX = REAL(GET_SLOT(GET_SLOT(x, lme4_XtXSym), lme4_xSym)),
	*Xty = REAL(GET_SLOT(x, lme4_XtySym)),
	*ZtX = REAL(GET_SLOT(GET_SLOT(x, lme4_ZtXSym), lme4_xSym)),
	*Zty = REAL(Ztyp),
	*wts = REAL(GET_SLOT(x, lme4_wtsSym)),
	one[] = {1, 0}, zero[] = {0,0};
    CHM_DN td1, Xd = AS_CHM_DN(Xp),
	wkrd = AS_CHM_DN(GET_SLOT(x, lme4_wrkresSym));
    CHM_DN Xcp = M_cholmod_copy_dense(Xd, &c),
	wkrcp = M_cholmod_copy_dense(wkrd, &c);
    CHM_SP Ztcp = M_cholmod_copy_sparse(Zt, &c);
    int *Zp = (int*)Ztcp->p;
    double *X = (double*)(Xcp->x), *Ztx = (double*)(Ztcp->x),
	*wtres = (double*)(wkrcp->x);
    R_CheckStack();
				/* Apply weights */
    for (i = 0; i < n; i++)
	wtres[i] = ((double*)(wkrd->x))[i] * wts[i];
    for (j = 0; j < p; j++)
	for (i = 0; i < n; i++)
	    X[i + j * n] = ((double*)(Xd->x))[i + j * n] * wts[i];
    for (j = 0; j < n; j++)
	for (i = Zp[j]; i < Zp[j + 1]; i++)
	    Ztx[i] = ((double*)(Zt->x))[i] * wts[j];

				/* y'y */
    REAL(GET_SLOT(x, lme4_devCompSym))[2] =
	F77_CALL(ddot)(&n, wtres, &ione, wtres, &ione);
				/* ZtZ */
    ts1 = M_cholmod_aat(Ztcp, (int *) NULL, (size_t) 0, 1/* mode */, &c);
    /* cholmod_aat returns stype == 0; copy to set stype == 1 */
    ts2 = M_cholmod_copy(ts1, 1/* stype */, 1/* mode */, &c);
    nnz = M_cholmod_nnz(ts2, &c);
    if (((int)(ts2->ncol) + 1) != LENGTH(ZtZp))
	error(_("Order of Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZp) - 1, (int)(ts2->ncol));
    /* double transpose to sort the columns */
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = M_cholmod_transpose(ts2, 1, &c);
    M_cholmod_free_sparse(&ts2, &c);
    ts2 = M_cholmod_transpose(ts1, 1, &c);
    Memcpy(INTEGER(ZtZp), (int*)(ts2->p), LENGTH(ZtZp));
    if (nnz != LENGTH(ZtZx))
	error(_("Number of nonzeros in Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZx), nnz);
    Memcpy(INTEGER(ZtZi), (int*)(ts2->i), nnz);
    Memcpy(REAL(ZtZx), (double*)(ts2->x), nnz);
    M_cholmod_free_sparse(&ts1, &c); M_cholmod_free_sparse(&ts2, &c);
				/* PZ'X into ZtX */
    td1 = M_cholmod_allocate_dense(q, p, q, CHOLMOD_REAL, &c);
    if (!M_cholmod_sdmult(Ztcp, 0, one, zero, Xcp, td1, &c))
	error(_("Error return from cholmod_sdmult"));
    for (j = 0; j < p; j++) { 	/* apply the permutation to each column */
	double *dcol = ZtX + j * q,
	    *scol = (double*)(td1->x) + j * q;
	for (i = 0; i < q; i++) dcol[i] = scol[perm[i]];
    }
    M_cholmod_free_dense(&td1, &c);
				/* PZ'y into Zty */
    td1 = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!M_cholmod_sdmult(Ztcp, 0, one, zero, wkrcp, td1, &c))
	error(_("Error return from cholmod_sdmult"));
    for (i = 0; i < q; i++) Zty[i] = ((double *)(td1->x))[perm[i]];
    M_cholmod_free_dense(&td1, &c); M_cholmod_free_sparse(&Ztcp, &c);
				/* XtX and Xty */
    AZERO(XtX, p * p);
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, XtX, &p);
    F77_CALL(dgemv)("T", &n, &p, one, X, &n, wtres, &ione, zero, Xty, &ione);
    M_cholmod_free_dense(&Xcp, &c); M_cholmod_free_dense(&wkrcp, &c);
    flag_not_factored(x);
}

/**
 * Update the fixef slot on a factored mer object.
 *
 * @param x Pointer to an mer object
 *
 * @return fixed effects vector
 */
static double*
internal_mer_fixef(SEXP x)
{
    SEXP fixef = GET_SLOT(x, lme4_fixefSym);
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));
    if (!status[0]) {
	error("Applying internal_mer_fixef without factoring");
	return (double*)NULL;	/* -Wall */
    }
    if (status[0] < 2) {
	int ione = 1, p = LENGTH(fixef);
	Memcpy(REAL(fixef), REAL(GET_SLOT(x, lme4_rXySym)), p);
	F77_CALL(dtrsv)("U", "N", "N", &p,
			REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym),
				      lme4_xSym)),
			&p, REAL(fixef), &ione);
    }
    return REAL(fixef);
}

/**
 * Perform ECME steps for the REML or ML criterion.
 *
 * @param x pointer to an mer object
 * @param nsteps pointer to an integer scalar - the number of ECME
 * steps to perform
 * @param Verbp pointer to a logical scalar indicating verbose output
 *
 * @return R_NilValue
 */
SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp)
{
    int nstp = asInteger(nsteps);
    if (nstp > 0) internal_ECMEsteps(x, nstp, asLogical(Verbp));
    return R_NilValue;
}

/**
 * Fill in five symmetric matrices, providing the information to
 * generate the Hessian.
 *
 * @param x pointer to an mer object
 *
 * @return an array consisting of five symmetric faces
 */
SEXP mer_Hessian(SEXP x)
{
    SEXP
	D = GET_SLOT(x, lme4_DSym),
	Omega = GET_SLOT(x, lme4_OmegaSym),
	RZXP = GET_SLOT(x, lme4_RZXSym),
	gradComp = GET_SLOT(x, lme4_gradCompSym),
	val;
    int *dRZX = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    double
	*RZX = REAL(RZXP),
	*b = REAL(RZXP) + dRZX[0] * p,
	*valp, *bbface,		/* vec of second faces of gradComp elts */
	one = 1.,
	zero = 0.;

    mer_gradComp(x);
    Q = 0;			/* number of rows and columns in the result */
    for (i = 0; i < nf; i++) Q += nc[i] * nc[i];
    Qsqr = Q * Q;
    bbface = Calloc(Q, double);
    val = PROTECT(alloc3Darray(REALSXP, Q, Q, 5));
    valp= REAL(val);
    AZERO(valp, Qsqr * 5);

    pos = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int ncisqr = nci * nci;
	double *fDi = REAL(VECTOR_ELT(gradComp, i)),
	    mult = 1./((double)(Gp[i + 1] - Gp[i])/nci);

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of gradComp on the diagonal
	 * of the third face of val */
	F77_CALL(dsyr)("U", &ncisqr, &mult, fDi + 2 * ncisqr, &ione,
		       valp + 2 * Qsqr + pos * Q, &Q);
	pos += ncisqr;
    }
				/* fifth face is outer product of bbface */
    F77_CALL(dsyr)("U", &Q, &one, bbface, &ione, valp + 4 * Qsqr, &Q);
				/* fourth face from \bb\trans\der\vb\der\bb */
    AZERO(valp + 3 * Qsqr, Qsqr); /* zero accumulator */
    RZXpos = 0;
    facepos = 0;
    for (i = 0; i < nf; i++) {
	int ii, jj, nci = nc[i];
	int ncisqr = nci * nci, nctp = nci * p,
	    nlev = (Gp[i + 1] - Gp[i])/nci;
	int maxpq = (p > nci) ? p : nci;
	double
	    *Di = REAL(VECTOR_ELT(D, i)),
	    *arr = Calloc(ncisqr * maxpq, double), /* tmp 3Darray */
	    *face = valp + 3 * Qsqr,
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
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 *
 * @param x pointer to an mer object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbosep, SEXP deviancep)
{
    SEXP ans, Omega = GET_SLOT(x, lme4_OmegaSym),
	Omegacp = PROTECT(duplicate(Omega));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1],
	i, j, n = LENGTH(GET_SLOT(x, lme4_ySym)),
	nf = LENGTH(Omega), nsamp = asInteger(nsampp),
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym)),
	saveb = asLogical(savebp),
	trans = asLogical(transp),
	verbose = asLogical(verbosep),
	deviance = asLogical(deviancep);
    double
	*RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	*bhat = REAL(GET_SLOT(x, lme4_ranefSym)),
	*betahat = REAL(GET_SLOT(x, lme4_fixefSym)),
	*dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	*ansp, df = n - (REML ? p : 0);
    int nrbase = p + 1 + coef_length(nf, nc); /* rows always included */
    int nrtot = nrbase + deviance + (saveb ? q : 0);
    double *bnew = Alloca(q, double), *betanew = Alloca(p, double);
    R_CheckStack();

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    GetRNGstate();
    if (verbose) Rprintf("%12s %14s\n", "sigma", "deviance");

    for (i = 0; i < nsamp; i++) {
	double *col = ansp + i * nrtot, sigma, dev;
				/* simulate and store new value of sigma */
	sigma = exp(dcmp[3]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
				/* simulate new fixed and random effects */
	internal_betab_update(p, q, sigma, L, RZX, RXX, betahat, bhat,
			      betanew, bnew);
				/* Store beta */
	for (j = 0; j < p; j++) col[j] = betanew[j];
				/* Optionally store b */
	if (saveb) for (j = 0; j < q; j++)
	    col[nrbase + deviance + j] = bnew[j];
				/* Update and store Omega */
	internal_Omega_update(Omega, bnew, sigma, nf, nc, Gp,
				     col + p + 1, trans);
	internal_mer_refactor(x);
	mer_secondary(x);

	dev = lmm_deviance(x, sigma, betanew);
	if (deviance) col[nrbase] = dev; /* store deviance */
	if (verbose) Rprintf("%12.6g %14.8g\n", sigma, dev);
    }
    PutRNGstate();
				/* Restore original Omega */
    SET_SLOT(x, lme4_OmegaSym, Omegacp);
    internal_mer_refactor(x);

    UNPROTECT(2);
    return ans;
}

/**
 * Create an mer object from a list of grouping factors and a list of model
 * matrices.
 *
 * @param fl named list of grouping factors
 * @param ZZt transpose of Z as a sparse matrix
 * @param Xp model matrix for the fixed effects
 * @param yp response vector
 * @param REMLp logical scalar indicating if REML is to be used
 * @param ncp integer vector of the number of random effects per level
 *        of each grouping factors
 *
 * @return pointer to an mer object
 */
SEXP mer_create(SEXP fl, SEXP ZZt, SEXP Xp, SEXP yp, SEXP REMLp,
		SEXP ncp, SEXP cnames)
{
    SEXP Omega, bVar, gradComp, fnms = getAttrib(fl, R_NamesSymbol),
	stat, val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer"))), xnms;
    CHM_SP ts1, ts2, Zt;
    cholmod_factor *F;
    int *nc = INTEGER(ncp), *Gp, *xdims, REML = asInteger(REMLp),
	i, nested, nf = LENGTH(fl), nobs = LENGTH(yp), p, q;
    char *devCmpnms[] = {"n", "p", "yty", "logryy2", "logDetL2",
			 "logDetOmega", "logDetRXX", "scale", ""};
    char *devnms[] = {"ML", "REML", ""};
    char *statusnms[] =
	{"stage", "REML", "glmm", ""};
    double *dcmp, *wts, *wrkres, *y;
				/* Check arguments to be duplicated */
    if (!isReal(yp)) error(_("yp must be a real vector"));
    SET_SLOT(val, lme4_ySym, duplicate(yp));
    if (!isMatrix(Xp) || !isReal(Xp))
	error(_("Xp must be a real matrix"));
    xdims = INTEGER(getAttrib(Xp, R_DimSymbol));
    if (xdims[0] != nobs) error(_("Xp must have %d rows"), nobs);
    p = xdims[1];
    xnms = VECTOR_ELT(getAttrib(Xp, R_DimNamesSymbol), 1);
    SET_SLOT(val, lme4_XSym, duplicate(Xp));
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);
    }
    SET_SLOT(val, lme4_flistSym, duplicate(fl));
    if (!isNewList(cnames) || LENGTH(cnames) != nf + 1)
	error(_("cnames must be a list of length %d"), nf + 1);
    SET_SLOT(val, lme4_cnamesSym, duplicate(cnames));
    if (!isInteger(ncp) || LENGTH(ncp) != nf)
	error(_("ncp must be an integer vector of length %d"), nf);
    SET_SLOT(val, lme4_ncSym, duplicate(ncp));
    SET_SLOT(val, lme4_ZtSym, duplicate(ZZt));
    Zt = AS_CHM_SP(GET_SLOT(val, lme4_ZtSym));
    R_CheckStack();
    q = Zt->nrow;
				/* allocate other slots */
    SET_SLOT(val, lme4_devianceSym, internal_make_named(REALSXP, devnms));
    SET_SLOT(val, lme4_devCompSym, internal_make_named(REALSXP, devCmpnms));
    SET_SLOT(val, lme4_statusSym, internal_make_named(INTSXP, statusnms));
    stat = GET_SLOT(val, lme4_statusSym);
    AZERO(INTEGER(stat), LENGTH(stat));
    INTEGER(stat)[1] = REML;
    dcmp = REAL(GET_SLOT(val, lme4_devCompSym));
    AZERO(dcmp, 8);		/* cosmetic */
    dcmp[0] = (double) nobs;
    dcmp[1] = (double) p;
    dcmp[7] = 1.;		/* initialize to a positive value */
				/* allocate and populate list slots */
    Omega = ALLOC_SLOT(val, lme4_OmegaSym, VECSXP, nf);
    bVar = ALLOC_SLOT(val, lme4_bVarSym, VECSXP, nf);
    gradComp = ALLOC_SLOT(val, lme4_gradCompSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    setAttrib(bVar, R_NamesSymbol, duplicate(fnms));
    setAttrib(gradComp, R_NamesSymbol, duplicate(fnms));
    Gp = INTEGER(ALLOC_SLOT(val, lme4_GpSym, INTSXP, nf + 1));
    Gp[0] = 0;
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int nlev = LENGTH(getAttrib(VECTOR_ELT(fl, i), R_LevelsSymbol));
	SET_VECTOR_ELT(Omega, i,
		       alloc_dpoMatrix(nci, "U",
					 VECTOR_ELT(cnames, i),
					 VECTOR_ELT(cnames, i)));
	SET_VECTOR_ELT(bVar, i, alloc3Darray(REALSXP, nci, nci, nlev));
	SET_VECTOR_ELT(gradComp, i, alloc3Darray(REALSXP, nci, nci, 4));
	Gp[i + 1] = Gp[i] + nlev * nci;
    }
				/* analyze Zt and ZtZ */
    ts1 = M_cholmod_aat(Zt, (int*)NULL/* fset */,(size_t)0,
		      CHOLMOD_PATTERN, &c);
    ts2 = M_cholmod_copy(ts1, -1/*lower triangle*/, CHOLMOD_PATTERN, &c);
    SET_SLOT(val, lme4_ZtZSym,
	     alloc_dsCMatrix(q, M_cholmod_nnz(ts2, &c), "U", R_NilValue,
			       R_NilValue));
    i = c.supernodal;
    c.supernodal = CHOLMOD_SUPERNODAL; /* force a supernodal decomposition */
    nested = internal_mer_isNested(nf, nc, Gp, (int*)(ts2->p));
    if (nested) {		/* require identity permutation */
	int nmethods = c.nmethods, ord0 = c.method[0].ordering;
	c.nmethods = 1;
	c.method[0].ordering = CHOLMOD_NATURAL;
	c.postorder = FALSE;
	F = M_cholmod_analyze(Zt, &c);
	c.nmethods = nmethods; c.method[0].ordering = ord0;
	c.postorder = TRUE;
    } else F = M_cholmod_analyze(Zt, &c);
    c.supernodal = i;		/* restore previous setting */
    M_cholmod_free_sparse(&ts1, &c); M_cholmod_free_sparse(&ts2, &c);
				/* create ZtX, RZX, XtX, RXX */
    SET_SLOT(val, lme4_ZtXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, lme4_RZXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, lme4_XtXSym, alloc_dpoMatrix(p, "U", xnms, xnms));
    SET_SLOT(val, lme4_RXXSym, alloc_dtrMatrix(p, "U", "N", xnms, xnms));
    SET_SLOT(val, lme4_ZtySym, allocVector(REALSXP, q));
    SET_SLOT(val, lme4_rZySym, allocVector(REALSXP, q));
    SET_SLOT(val, lme4_XtySym, allocVector(REALSXP, p));
    SET_SLOT(val, lme4_rXySym, allocVector(REALSXP, p));
				/* create weights and working residuals */
    wts = REAL(ALLOC_SLOT(val, lme4_wtsSym, REALSXP, nobs));
    wrkres = REAL(ALLOC_SLOT(val, lme4_wrkresSym, REALSXP, nobs));
    y = REAL(GET_SLOT(val, lme4_ySym));
    for (i = 0; i < nobs; i++) {wts[i] = 1.; wrkres[i] = y[i];}
    internal_mer_update_ZXy(val, (int*)(F->Perm));
				/* secondary slots */
    SET_SLOT(val, lme4_ranefSym, allocVector(REALSXP, q));
    SET_SLOT(val, lme4_fixefSym, allocVector(REALSXP, p));
    SET_SLOT(val, lme4_RZXinvSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
				/* initialize */
    mer_initial(val);
    /* The next calls are simply to set up the L slot.  At present the
     * factor F is a symbolic factor.  We need to force it to become
     * numeric before allocating the L slot in the object. */
    internal_mer_Zfactor(val, F);
    /* One side-effect of this call is to set the status as
     * factored.  We need to reset it */
    INTEGER(GET_SLOT(val, lme4_statusSym))[0] = 0;
    /* Create the dCHMfactor object and store it in the L slot.  This
     * also frees the storage. */
    SET_SLOT(val, lme4_LSym, M_chm_factor_to_SEXP(F, 1));
    /* OK, done now. */
    UNPROTECT(1);
    return val;
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
 * @param x pointer to an mer object
 * @param pType pointer to an integer scalar indicating the form of the
 *        parameters to be returned.
 *
 * @return numeric vector of the values in the upper triangles of the
 * Omega matrices
 */
SEXP mer_coef(SEXP x, SEXP pType)
{
    int	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	nf = LENGTH(GET_SLOT(x, lme4_OmegaSym));
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));

    internal_mer_coef(x, asInteger(pType), REAL(val));
    UNPROTECT(1);
    return val;
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
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType)
{
    int clen = coef_length(LENGTH(GET_SLOT(x, lme4_flistSym)),
			   INTEGER(GET_SLOT(x, lme4_ncSym)));
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_mer_coefGets(x, REAL(coef), asInteger(pType));
    return x;
}

/**
 * The naive way of calculating the trace of the hat matrix
 *
 * @param x pointer to an mer object
 *
 * @return trace of the hat matrix
 */

SEXP mer_hat_trace(SEXP x)
{
    SEXP Zt = GET_SLOT(x, lme4_ZtSym);
    int *Zti = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	*Ztp = INTEGER(GET_SLOT(Zt, lme4_pSym)), i, ione = 1, j,
	n = INTEGER(GET_SLOT(Zt, lme4_DimSym))[1],
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));
    double *RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	*Ztx = REAL(GET_SLOT(Zt, lme4_xSym)),
	*wrk = Calloc(q, double), m1 = -1, one = 1, tr;
    double *Xcp = Alloca(n * p, double);
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    CHM_DN zrow = N_AS_CHM_DN(wrk, q, 1);
    R_CheckStack();

    mer_factor(x);
    Memcpy(Xcp, REAL(GET_SLOT(x, lme4_XSym)), n * p);

    /* Accumulate F-norm of L^{-1}Zt and downdate rows of Xcp */
/* FIXME: Does this handle a non-trivial permutation properly? */
    for (j = 0, tr = 0; j < n; j++) { /* j'th column of Zt */
	CHM_DN sol; double *sx;
	for (i = 0; i < q; i++) wrk[i] = 0;
	for (i = Ztp[j]; i < Ztp[j + 1]; i++) wrk[Zti[i]] = Ztx[i];
	sol = M_cholmod_solve(CHOLMOD_L, L, zrow, &c);
	sx = (double*)(sol->x);
	for (i = 0; i < q; i++) tr += sx[i] * sx[i];
				/* downdate jth row of Xcp */
 	F77_CALL(dgemv)("T", &q, &p, &m1, RZX, &q, sx, &ione,
 			&one, Xcp + j, &n);
	M_cholmod_free_dense(&sol, &c);
    }
    F77_CALL(dtrsm)("R", "U", "N", "N", &n, &p, &one, RXX, &p, Xcp, &n);
    for (i = 0; i < n * p; i++) tr += Xcp[i] * Xcp[i];

    return ScalarReal(tr);
}

/**
 * Faster calculation of the trace of the hat matrix (due to Jialiang Li)
 *
 * @param x pointer to an mer object
 *
 * @return trace of the hat matrix
 */

SEXP mer_hat_trace2(SEXP x)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	ncp = GET_SLOT(x, lme4_ncSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(ncp),
	nf = LENGTH(ncp), i, j, k,
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    double one = 1, tr = p + q,
	*RZXicp = Alloca(q * p, double);
    R_CheckStack();
				/* factor and evaluate RZXinv */
    mer_factor(x);
    internal_mer_RZXinv(x);
    Memcpy(RZXicp, REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym),
				 lme4_xSym)), q * p);
    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int ncisqr = nci * nci, nlev = (Gp[i + 1] - Gp[i])/nci;
	double *deli = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				     lme4_xSym));
	CHM_SP sol, Prhs,
	    rhs = M_cholmod_allocate_sparse(q, nci, ncisqr, TRUE, TRUE,
					    0, CHOLMOD_REAL, &c);
	int *rhi = (int *)(rhs->i), *rhp = (int *)(rhs->p);
	double *rhx = (double *)(rhs->x);

	rhp[0] = 0;		/* Establish column pointers and value of rhs */
	for (j = 0; j < nci; j++) {
	    rhp[j+1] = rhp[j] + nci;
	    for (k = 0; k < nci; k++) { /* transpose of deli */
		rhx[j * nci + k] = deli[k * nci + j];
		rhi[j * nci + k] = Gp[i] + k; /* initial row numbers */
	    }
	    /* zero the upper triangle (just in case) */
	    for (k = 0; k < j; k++) rhx[j * nci + k] = 0;
	}
	for (j = 0; j < nlev; j++) {
	    /* Evaluate nci rows of Delta RZXinv */
	    F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &p, &one, deli, &nci,
			    RZXicp + Gp[i] + j * nci, &q);
	    /* Solve for and accumulate nci columns of L^{-1} P Delta' */
	    Prhs = M_cholmod_spsolve(CHOLMOD_P, L, rhs, &c);
	    sol = M_cholmod_spsolve(CHOLMOD_L, L, Prhs, &c);
	    M_cholmod_free_sparse(&Prhs, &c);
	    for (k = 0; k < ((int*)(sol->p))[nci]; k++) {
		double sxk = ((double*)(sol->x))[k];
		tr -= sxk * sxk;
	    }
	    M_cholmod_free_sparse(&sol, &c);
	    /* Update rhs for the next set of rows */
	    for (k = 0; k < ncisqr; k++) rhi[k] += nci;
	}
	M_cholmod_free_sparse(&rhs, &c);
    }
    for (i = 0; i < q * p; i++) tr -= RZXicp[i] * RZXicp[i];
    return ScalarReal(tr);
}

/**
 * Return L as a dtCMatrix object
 *
 * @param x pointer to an mer object
 *
 * @return L as an dtCMatrix object
 */
SEXP mer_dtCMatrix(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2)), nz, q;
    CHM_SP Lm;
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym)), Lcp;
    R_CheckStack();

    dims[0] = dims[1] = q = (int)(L->n);
    Lcp = M_cholmod_copy_factor(L, &c); /* next call changes Lcp */
    Lm = M_cholmod_factor_to_sparse(Lcp, &c); M_cholmod_free_factor(&Lcp, &c);
    SET_SLOT(ans, lme4_uploSym, mkString("L"));
    SET_SLOT(ans, lme4_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, q + 1)),
	   (int*) Lm->p, q + 1);
    nz = ((int*)(Lm->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz)),
	   (int*) Lm->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz)),
	   (double*) Lm->x, nz);
    M_cholmod_free_sparse(&Lm, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Return L^{-1} as a dtCMatrix object
 *
 * @param x pointer to an mer object
 *
 * @return L^{-1} as an dtCMatrix object
 */
SEXP mer_dtCMatrix_inv(SEXP x)
{
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    cholmod_sparse
	*b = M_cholmod_allocate_sparse(L->n, L->n, L->n, 1, 1,
				       0, CHOLMOD_REAL, &c),
	*Linv;
    double *bx = (double *)(b->x);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *bi = (int *) (b->i), *bp = (int *) (b->p),
	*dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2)),
	j, nz, q;

    dims[0] = dims[1] = q = (int)(L->n);
    for (j = 0; j < q; j++) {
	bp[j] = bi[j] = j;
	bx[j] = 1;
    }
    bp[q] = q;
    Linv = M_cholmod_spsolve(CHOLMOD_L, L, b, &c);
    M_cholmod_free_sparse(&b, &c);
    SET_SLOT(ans, lme4_uploSym, mkString("L"));
    SET_SLOT(ans, lme4_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, q + 1)),
	   (int *) Linv->p, q + 1);
    nz = ((int *)(Linv->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz)),
	   (int *) Linv->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz)),
	   (double *) Linv->x, nz);
    M_cholmod_free_sparse(&Linv, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Create and factor Z'Z+Omega if it has not already.
 * Also create RZX and RXX, the deviance components,
 * and the value of the deviance for both ML and  REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP mer_factor(SEXP x)
{
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));
    if (!status[0]) {
	SEXP rXyP = GET_SLOT(x, lme4_rXySym),
	    rZyP = GET_SLOT(x, lme4_rZySym);
	int i, info, ione = 1, p = LENGTH(rXyP), q = LENGTH(rZyP);
	double *RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	    *RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	    *rXy = REAL(rXyP), *rZy = REAL(rZyP),
	    *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	    *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	    one[2] = {1, 0}, m1[2] = {-1, 0};
	double nml = dcmp[0], nreml = dcmp[0] - dcmp[1];
	CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
	R_CheckStack();

	/* Inflate Z'Z to Z'Z+Omega and factor to form L. Form RZX and
	 * rZy. Update stage flag, dcmp[4] and dcmp[5]. */
	internal_mer_Zfactor(x, L);
				/* downdate XtX and factor */
	if ((info = internal_mer_Xfactor(x))) /* unable to factor downdated XtX */
	    error(_("Leading minor of order %d in downdated X'X is not positive definite"),
		  info);
	for (dcmp[6] = 0, i = 0; i < p; i++) /* 2 * logDet(RXX) */
	    dcmp[6] += 2. * log(RXX[i * (p + 1)]);
				/* solve for rXy  and ryy^2 */
	Memcpy(rXy, REAL(GET_SLOT(x, lme4_XtySym)), p);
	F77_CALL(dgemv)("T", &q, &p, m1, RZX, &q, rZy, &ione, one, rXy, &ione);
	F77_CALL(dtrsv)("U", "T", "N", &p, RXX, &p, rXy, &ione);
	dcmp[3] = log(dcmp[2] /* dcmp[3] = log(ryy^2); dcmp[2] = y'y; */
		      - F77_CALL(ddot)(&p, rXy, &ione, rXy, &ione)
		      - F77_CALL(ddot)(&q, rZy, &ione, rZy, &ione));
				/* evaluate ML and REML deviance */
	dev[0] = dcmp[4] - dcmp[5] +
	    nml*(1.+dcmp[3]+log(2.*PI/nml));
	dev[1] = dcmp[4] - dcmp[5] + dcmp[6] +
	    nreml*(1.+dcmp[3]+log(2.*PI/nreml));
	if (dcmp[7] >= 0) dcmp[7] = internal_mer_sigma(x, -1);
	status[0] = 1;
    }
    return R_NilValue;
}

/**
 * Return the fitted values as an SEXP
 *
 * @param x pointer to an mer object
 * @param useFe pointer to a logical scalar indicating if the fixed
 * effects should be used
 * @param useRe pointer to a logical scalar indicating if the random
 * effects should be used
 *
 * @return pointer to a numeric array of fitted values
 */

SEXP mer_fitted(SEXP x)
{
    int n = LENGTH(GET_SLOT(x, lme4_ySym));
    SEXP ans = PROTECT(allocVector(REALSXP, n));

    internal_mer_fitted(x, (double*) NULL, REAL(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the conditional estimates of the fixed effects
 *
 * @param x Pointer to an mer object
 *
 * @return a numeric vector containing the conditional estimates of
 * the fixed effects
 */
SEXP mer_fixef(SEXP x)
{
    int nf = LENGTH(GET_SLOT(x, lme4_OmegaSym));
    SEXP ans;

    mer_secondary(x);
    ans = PROTECT(duplicate(GET_SLOT(x, lme4_fixefSym)));
    setAttrib(ans, R_NamesSymbol,
	      duplicate(VECTOR_ELT(GET_SLOT(x, lme4_cnamesSym), nf)));
    UNPROTECT(1);
    return ans;
}

/**
 * Fill in the gradComp and bVar slots.  Each component in the gradComp slot
 * consists of four symmetric matrices used to generate the gradient or the ECME
 * step.  They are
 *  1) -m_i\bOmega_i^{-1}
 *  2) \bB_i\bB_i\trans
 *  3) \tr\left[\der_{\bOmega_i}\bOmega\left(\bZ\trans\bZ+\bOmega\right)\inv\right]
 *  4) The term added to 3) to get \tr\left[\der_{\bOmega_i}\bOmega\vb\right]
 *
 * @param x pointer to an mer object
 * @param val pointer to a list of matrices of the correct sizes
 *
 * @return NULL
 */
SEXP mer_gradComp(SEXP x)
{
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));

    if (status[0] < 3) {
	SEXP bVarP = GET_SLOT(x, lme4_bVarSym),
	    OmegaP = GET_SLOT(x, lme4_OmegaSym),
	    gradComp = GET_SLOT(x, lme4_gradCompSym),
	    ranefP = GET_SLOT(x, lme4_ranefSym);
	int q = LENGTH(ranefP), p = LENGTH(GET_SLOT(x, lme4_rXySym));
	int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	    *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	    i, j, k, nf = length(OmegaP);
	double *b = REAL(GET_SLOT(x, lme4_ranefSym)),
	    *RZXinv = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym),
				    lme4_xSym)),
	    alpha;
/* 	CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym)); */
	R_CheckStack();

	mer_secondary(x);
	alpha = 1./internal_mer_sigma(x, -1);
	alpha = alpha * alpha;

	internal_mer_RZXinv(x);
	internal_mer_bVar(x);
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], RZXrows = Gp[i + 1] - Gp[i];
	    int ncisq = nci * nci, nlev = RZXrows/nci;
	    double *bVi = REAL(VECTOR_ELT(bVarP, i)),
		*bi = b + Gp[i], *mm = REAL(VECTOR_ELT(gradComp, i)),
		*tmp = Memcpy(Calloc(ncisq, double),
			      REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(OmegaP, i)),
					    lme4_xSym)), ncisq),
		*RZXi = RZXinv + Gp[i],
		dlev = (double) nlev,
		one[] = {1,0}, zero[] = {0,0};

	    if (nci == 1) {
		int ione = 1;
		mm[0] = ((double) nlev)/(tmp[0] * tmp[0]);
		mm[1] = alpha * F77_CALL(ddot)(&nlev, bi, &ione, bi, &ione);
		mm[2] = 0.;
		for (k = 0; k < nlev; k++) mm[2] += bVi[k];
		mm[3] = 0.;
		for (j = 0; j < p; j++) {
		    mm[3] += F77_CALL(ddot)(&RZXrows, RZXi + j * q, &ione,
					    RZXi + j * q, &ione);
		}
	    } else {
		AZERO(mm, 4 * ncisq);
		F77_CALL(dtrtri)("U", "N", &nci, tmp, &nci, &j);
		if (j)
		    error(_("Omega[[%d]] is not positive definite"), i + 1);
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &dlev, tmp, &nci,
				zero, mm, &nci);
		mm += ncisq;	/* \bB_i term */
		F77_CALL(dsyrk)("U", "N", &nci, &nlev, &alpha, bi, &nci,
				zero, mm, &nci);
		mm += ncisq;     /* Sum of diagonal blocks of the inverse
				  * (Z'Z+Omega)^{-1} */
		for (j = 0; j < ncisq; j++) {
		    for (k = 0; k < nlev; k++) mm[j] += bVi[j + k*ncisq];
		}
		mm += ncisq;	/* Extra term for \vb */
		for (j = 0; j < p; j++) {
		    F77_CALL(dsyrk)("U", "N", &nci, &nlev, one,
				    RZXi + j * q, &nci,
				    one, mm, &nci);
		}
	    }
	    Free(tmp);
	}
	status[0] = 3;
    }
    return R_NilValue;
}

/**
 * Evaluate the gradient vector
 *
 * @param x Pointer to an lmer object
 * @param pType Pointer to an integer indicator of the
 * parameterization being used
 *
 * @return pointer to a gradient vector
 */
SEXP mer_gradient(SEXP x, SEXP pType)
{
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym);
    SEXP gradComp = GET_SLOT(x, lme4_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	dind, i, ifour = 4, ione = 1, nf = length(Omega),
	odind, ptyp = asInteger(pType);
    int REML = INTEGER(GET_SLOT(x, lme4_statusSym))[1];
    SEXP val = PROTECT(allocVector(REALSXP, coef_length(nf, nc)));
    double cc[] = {-1, 1, 1, REML ? 1 : 0},
      	*valp = REAL(val),
	one = 1.0, zero = 0.0;

    mer_gradComp(x);
    dind = 0;			/* index into val for diagonals */
    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omega, i);
	int nci = nc[i], ncisqr = nci * nci;
	double *tmp = Calloc(ncisqr, double);

	F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &one,
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			cc, &ifour, &zero, tmp, &ncisqr);
	if (nci == 1) {
	    double omg = *REAL(GET_SLOT(Omgi, lme4_xSym));
	    valp[dind++] =
		(ptyp?((ptyp == 1)? omg : -omg * omg) : 1) * tmp[0];
	} else {
	    int ii, j, ncip1 = nci + 1;

	    odind = dind + nci; /* index into val for off-diagonals */
	    if (ptyp) {
		double *chol = REAL(GET_SLOT(M_dpoMatrix_chol(Omgi), lme4_xSym)),
		    *tmp2 = Calloc(ncisqr, double);

		/* Overwrite the gradient with respect to positions in
		 * Omega[[i]] by the gradient with respect to the
		 * unconstrained parameters.*/

		/* tmp2 := chol %*% tmp using only upper triangle of tmp */
		F77_CALL(dsymm)("R", "U", &nci, &nci, &one, tmp, &nci,
				chol, &nci, &zero, tmp2, &nci);
		/* full symmetric product gives diagonals */
		F77_CALL(dtrmm)("R", "U", "T", "N", &nci, &nci, &one, chol,
				&nci, Memcpy(tmp, tmp2, ncisqr), &nci);
		/* overwrite upper triangle with gradients for L' */
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
		Free(tmp2);
	    }
	    for (j = 0; j < nci; j++) {
		valp[dind + j] = tmp[j * ncip1];
		for (ii = 0; ii < j; ii++) /* offdiagonals count twice */
		    valp[odind++] = 2. * tmp[ii + j * nci];
	    }
	    dind = odind;
	}
	Free(tmp);
    }
    UNPROTECT(1);
    return val;
}

/**
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an mer object
 *
 * @return NULL
 */
SEXP mer_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, lme4_OmegaSym),
	ZtZ = GET_SLOT(x, lme4_ZtZSym);
    int	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	*p = INTEGER(GET_SLOT(ZtZ, lme4_pSym)),
	i, nf = length(Omg);
    double *xx = REAL(GET_SLOT(ZtZ, lme4_xSym));

    for (i = 0; i < nf; i++) {
	SEXP Omgi = VECTOR_ELT(Omg, i);
	double *omgi = REAL(GET_SLOT(Omgi, lme4_xSym));
	int bb = Gp[i], j, k, nci = nc[i];
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - bb)/nci;

	AZERO(omgi, nci * nci);
	for (j = 0; j < nlev; j++) {
	    int base = bb + j * nci;
	    for (k = 0; k < nci; k++)
		/* add the last element in the column */
		omgi[k * ncip1] += xx[p[base + k + 1] - 1];
	}
	for (k = 0; k < nci; k++) omgi[k * ncip1] *= 0.375/nlev;
	SET_SLOT(Omgi, lme4_factorSym, allocVector(VECSXP, 0));
	M_dpoMatrix_chol(Omgi);
    }
    flag_not_factored(x);
    return R_NilValue;
}

/**
 * Externally callable check on nesting
 *
 * @param x Pointer to an mer object
 *
 * @return a scalar logical value indicating if ZtZ corresponds to a
 * simple nested structure.
 */
SEXP mer_isNested(SEXP x)
{
    CHM_SP ZtZ = AS_CHM_SP(GET_SLOT(x, lme4_ZtZSym));
    CHM_SP ZtZl = M_cholmod_transpose(ZtZ, (int) ZtZ->xtype, &c);
    SEXP ncp = GET_SLOT(x, lme4_ncSym);
    int ans = internal_mer_isNested(LENGTH(ncp), INTEGER(ncp),
				    INTEGER(GET_SLOT(x, lme4_GpSym)),
				    ZtZl->p);
    M_cholmod_free_sparse(&ZtZl, &c);
    return ScalarLogical(ans);
}

/**
 * Extract the conditional modes of the random effects.
 *
 * @param x Pointer to an mer object
 *
 * @return a list of matrices containing the conditional modes of the
 * random effects
 */
SEXP mer_ranef(SEXP x)
{
    SEXP cnames = GET_SLOT(x, lme4_cnamesSym),
	flist = GET_SLOT(x, lme4_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, ii, jj,
	nf = length(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, lme4_ranefSym));

    mer_secondary(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
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
		mm[ii + jj * mi] = bi[jj + ii * nci];
    }
    UNPROTECT(1);
    return val;
}

/**
 * Update the secondary slots - fixef and ranef
 *
 * @param x pointer to a mer object
 *
 */
SEXP mer_secondary(SEXP x)
{
    int *status = INTEGER(GET_SLOT(x, lme4_statusSym));

    if (status[0] < 2) {
	mer_factor(x);
	internal_mer_fixef(x);
	internal_mer_ranef(x);
	status[0] = 2;
    }
    return R_NilValue;
}

/**
 * Extract the posterior variances of the random effects
 *
 * @param x pointer to a mer object
 *
 * @return pointer to a list of arrays
 */
SEXP mer_postVar(SEXP x)
{
    SEXP ans;
    double *RZXi = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXinvSym), lme4_xSym)),
	*dcmp = REAL(GET_SLOT(x, lme4_devCompSym)), one = 1, sc;
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ione = 1, nf,
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym));

    sc = dcmp[7] * dcmp[7];
    mer_factor(x);
    internal_mer_RZXinv(x);
    internal_mer_bVar(x);
    ans = PROTECT(duplicate(GET_SLOT(x, lme4_bVarSym)));
    nf = LENGTH(ans);
    for (i = 0; i < nf; i++) {
	SEXP ansi = VECTOR_ELT(ans, i);
	int *dims = INTEGER(getAttrib(ansi, R_DimSymbol));
	int j, nci = dims[1], nlev = dims[2], ntot = LENGTH(ansi);
	int ncisqr = nci * nci;
	double *vv = REAL(ansi);

	if (dims[0] != nci)
	    error(_("rows and columns of element %d of bVar do not match"),
		  i + 1);
	for (j = 0; j < nlev; j++)
	    F77_CALL(dsyrk)("U", "N", &nci, &p,
			    &one, RZXi + Gp[i] + j * nci, &q,
			    &one, vv + j * ncisqr, &nci);
	if (sc != 1) F77_CALL(dscal)(&ntot, &sc, vv, &ione);
	if (nci > 1) {
	    for (j = 0; j < nlev; j++)
		internal_symmetrize(vv + j * ncisqr, nci);
	}
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the ML or REML conditional estimate of sigma
 *
 * @param x pointer to an mer object
 * @param REML logical scalar - TRUE if REML estimates are requested
 *
 * @return pointer to a numeric scalar
 */
SEXP mer_sigma(SEXP x, SEXP REML)
{
    return ScalarReal(
	internal_mer_sigma(x,
			   (REML == R_NilValue) ? -1 :
			   (asLogical(REML))));
}

/**
 * Simulate a set of linear predictors from the random effects part of
 * an mer object
 *
 * @param x Pointer to an mer object
 * @param np Pointer to an integer giving the number of values to simulate
 *
 * @return a matrix of simulated linear predictors
 */
SEXP mer_simulate(SEXP x, SEXP nsimP)
{
    int *nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ii, j, nsim = asInteger(nsimP),
	nf = LENGTH(GET_SLOT(x, lme4_OmegaSym)),
	n = LENGTH(GET_SLOT(x, lme4_ySym)),
	q = LENGTH(GET_SLOT(x, lme4_ZtySym));
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, nsim)),
	Omega = GET_SLOT(x, lme4_OmegaSym);
    double *dcmp = REAL(GET_SLOT(x, lme4_devCompSym)),
	one[] = {1,0}, zero[] = {0,0};
    double scale = (dcmp[7] < 0) ? -dcmp[7] : dcmp[7];
    CHM_DN cha = AS_CHM_DN(ans),
	chb = M_cholmod_allocate_dense(q, nsim, q, CHOLMOD_REAL, &c);
    CHM_SP Zt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym));
    R_CheckStack();

    GetRNGstate();
    for (ii = 0; ii < nsim; ii++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], relen = Gp[i + 1] - Gp[i];
	    int nlev = relen/nci;
	    double *bi = (double *)(chb->x) + ii * q + Gp[i],
		*Rmat = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				      lme4_xSym));

	    for (j = 0; j < relen; j++) bi[j] = norm_rand();
	    F77_CALL(dtrsm)("L", "U", "N", "N", &nci, &nlev, &scale,
			    Rmat, &nci, bi, &nci);
	}
    }
    PutRNGstate();

    if (!M_cholmod_sdmult(Zt, 1, one, zero, chb, cha, &c))
	error(_("cholmod_sdmult failed"));
    M_cholmod_free_dense(&chb, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Externally callable version of internal_mer_update_ZXy
 *
 * @param x pointer to an mer object
 *
 * @return NULL
 */
SEXP mer_update_ZXy(SEXP x)
{
    internal_mer_update_ZXy(x,
			    INTEGER(GET_SLOT(GET_SLOT(x, lme4_LSym),
					     lme4_permSym)));
    return R_NilValue;
}

/**
 * Update the y slot (and slots derived from it) in an mer object
 *
 * @param x pointer to an mer object
 * @param ynew pointer to a numeric vector of length n
 *
 * @return NULL
 */
SEXP mer_update_y(SEXP x, SEXP ynew)
{
    SEXP y = GET_SLOT(x, lme4_ySym),
	Xty = GET_SLOT(x, lme4_XtySym),
	Zty = GET_SLOT(x, lme4_ZtySym);
    double one[] = {1,0}, zero[] = {0,0};
    int *perm, i, ione = 1, n = LENGTH(y), p = LENGTH(Xty), q = LENGTH(Zty);
    CHM_SP Zt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym));
    CHM_DN td1, yd = AS_CHM_DN(GET_SLOT(x, lme4_ySym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    R_CheckStack();

    perm = (int*)(L->Perm);
    if (!isReal(ynew) || LENGTH(ynew) != n)
	error(_("ynew must be a numeric vector of length %d"), n);
    Memcpy(REAL(y), REAL(ynew), n);
    				/* y'y */
    REAL(GET_SLOT(x, lme4_devCompSym))[2] =
	F77_CALL(ddot)(&n, REAL(y), &ione, REAL(y), &ione);
				/* PZ'y into Zty */
    td1 = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!M_cholmod_sdmult(Zt, 0, one, zero, yd, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (i = 0; i < q; i++) REAL(Zty)[i] = ((double *)(td1->x))[perm[i]];
    M_cholmod_free_dense(&td1, &c);
    				/* Xty */
    F77_CALL(dgemv)("T", &n, &p, one, REAL(GET_SLOT(x, lme4_XSym)),
		    &n, REAL(y), &ione, zero, REAL(Xty), &ione);
    flag_not_factored(x);
    return R_NilValue;
}

/**
 * Check validity of an mer object.
 *
 * @param x Pointer to an mer object
 *
 * @return TRUE if the object is a valid lmer object, else a string
 * describing the nature of the violation.
 */
SEXP mer_validate(SEXP x)
{
    SEXP
	ZtXP = GET_SLOT(x, lme4_ZtXSym),
	XtXP = GET_SLOT(x, lme4_XtXSym),
	RZXP = GET_SLOT(x, lme4_RZXSym),
	RXXP = GET_SLOT(x, lme4_RXXSym)
	/* , cnames = GET_SLOT(x, lme4_cnamesSym) */
	;
    int *ZtXd = INTEGER(getAttrib(ZtXP, lme4_DimSym)),
	*XtXd = INTEGER(getAttrib(XtXP, lme4_DimSym));

    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, lme4_DimSym))))
	return mkString(_("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, lme4_DimSym))))
	return mkString(_("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return mkString(_("Slot XtX must be a square matrix with same ncol as ZtX"));
    return ScalarLogical(1);
}

static
SEXP SEXP_Zt(int n, int ii, SEXP fi, SEXP tmmat)
{
    int *dims = INTEGER(getAttrib(tmmat, R_DimSymbol)), *fac =INTEGER(fi), j, k;
    int m = dims[0], nlev = LENGTH(getAttrib(fi, R_LevelsSymbol));
    SEXP ans = PROTECT(alloc_dgCMatrix(m * nlev, n, m * n, R_NilValue, R_NilValue));
    int *i = INTEGER(GET_SLOT(ans, lme4_iSym)), *p = INTEGER(GET_SLOT(ans, lme4_pSym));

    if (!isFactor(fi) || LENGTH(fi) != n)
	error(_("fl[[%d]] must be a factor of length %d"), ii + 1, n);
    if (!isMatrix(tmmat) || !isReal(tmmat))
	error(_("Ztl[[%d]] must be real matrix"), ii + 1);
    if (dims[1] != n)
	error(_("Ztl[[%d]] must have %d columns"), ii + 1, n);

    p[n] = m * n;
    for (j = 0; j < n; j++) {
	p[j] = m * j;
	for (k = 0; k < m; k++) i[j * m + k] = (fac[j] - 1) * m + k;
    }
    Memcpy(REAL(GET_SLOT(ans, lme4_xSym)), REAL(tmmat), m * n);
    UNPROTECT(1);
    return ans;
}

/**
 * Create a list of sparse Zt matrices from a factor list and a list
 * of dense, skinny model matrices
 *
 * @param fl list of factors
 * @param Ztl list of transposes of model matrices
 *
 * @return a list of sparse (full) Zt matrices
 */
SEXP Ztl_sparse(SEXP fl, SEXP Ztl)
{
    int i, nf = LENGTH(fl), nobs = LENGTH(VECTOR_ELT(fl, 0));
    SEXP ans = PROTECT(allocVector(VECSXP, nf));

    setAttrib(ans, R_NamesSymbol, duplicate(getAttrib(fl, R_NamesSymbol)));
    for (i = 0; i < nf; i++)
	SET_VECTOR_ELT(ans, i, SEXP_Zt(nobs, i, VECTOR_ELT(fl, i), VECTOR_ELT(Ztl, i)));
    UNPROTECT(1);
    return ans;
}

/**
 * Create a new sparse Zt matrix by carrying over elements for the
 * same level of f
 *
 * @param f factor determining the carryover (e.g. student)
 * @param Zt sparse model matrix for another factor (e.g. teacher)
 * @param tvar numeric vector of time values
 * @param discount numeric vector of discounting fractions
 *
 * @return modified model matrix
 */
SEXP Zt_carryOver(SEXP fp, SEXP Zt, SEXP tvar, SEXP discount)
{
    int *cct, *f = INTEGER(fp), *p;
    int cmax, j, jj, k, last, n = LENGTH(fp), nlev, nnz, ntot, q;
    int *ii, *ij, *oi, *oj, dl = LENGTH(discount), ip, op;
    double *ix, *ox, *disc, *tv;
    CHM_SP ans, chsz = AS_CHM_SP(Zt);
    CHM_TR ant, chtz = M_cholmod_sparse_to_triplet(chsz, &c);
    R_CheckStack();

    p = (int*)(chsz->p);
    q = p[1] - p[0];
    if (!isReal(discount))
	error(_("discount must be a numeric vector"));
    if (!isReal(tvar))
	error(_("tvar must be a numeric vector"));
    if (LENGTH(tvar) != n)
	error(_("tvar must have length %d"), n);
    tv = REAL(tvar);
    disc = REAL(discount);
    if (!isFactor(fp)) error(_("f must be a factor"));
    nlev = LENGTH(getAttrib(fp, R_LevelsSymbol));
    cct = Alloca(nlev, int);
    R_CheckStack();

    if (chtz->ncol != n) error(_("ncol(Zt) must match length(fp)"));
    for (j = 0; j < n; j++)	/* check consistency of p */
	if (p[j+1] - p[j] != q)
	    error(_("nonzeros per column in Zt must be constant"));
				/* create column counts */
    for (last = -1, j = 0; j < n; j++) {
	int ll = f[j] - 1;
	if (last > ll) error(_("f is not in increasing order"));
	if (ll == last) cct[ll]++;
	else {cct[ll] = 1; last = ll;};
    }
    if (nlev != f[n - 1]) error(_("number of levels of f is not consistent"));
				/*  max column count and total nonzeros*/
    for (cmax = -1, ntot = 0, k = 0; k < nlev; k++) {
	cmax = (cct[k] > cmax) ? cct[k] : cmax;
	ntot += (cct[k] * (cct[k] + 1))/2;
    }
    nnz = ntot * q;
    ant = M_cholmod_allocate_triplet(chtz->nrow, chtz->ncol, (size_t)(nnz),
				     0 /*stype*/, CHOLMOD_REAL, &c);
    ip = 0; ii = (int*)(chtz->i); ij = (int*)(chtz->j); ix = (double*)(chtz->x);
    op = 0; oi = (int*)(ant->i); oj = (int*)(ant->j); ox = (double*)(ant->x);
    for (k = 0; k < nlev; k++) {
	for (j = 0; j < cct[k]; j++) {
	    for (jj = 0; jj < cct[k] - j; jj++) {
		int dj = (int)(tv[ip + jj] - tv[ip]);
		if (dj > dl)
		    error(_("diff(tvar) (= %d) > length(discount) (= %d)"),
			  dj, dl);
		oi[op] = ii[ip]; oj[op] = ij[ip] + jj;
		ox[op] = ix[ip] * disc[dj];
		op++;
	    }
	    ip++;
	}
    }
    ant->nnz = nnz;
    ans = M_cholmod_triplet_to_sparse(ant, nnz, &c);
    M_cholmod_free_triplet(&chtz, &c);
    M_cholmod_free_triplet(&ant, &c);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "", R_NilValue);
}

