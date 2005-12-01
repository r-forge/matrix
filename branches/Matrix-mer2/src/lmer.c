#include "lmer.h"
#include <float.h>

/* These should be moved to the R sources */

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
internal_symmetrize(double *a, int nc)
{
    int i, j;

    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/** 
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and df degrees of freedom.
 * 
 * @param df degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param ans array of size p * p to hold the result
 * 
 * @return ans
 */
static double*
std_rWishart_factor(double df, int p, double ans[])
{
    int i, j, pp1 = p + 1;

    if (df < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");
    for (j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(df - (double) j));
	for (i = 0; i < j; i++) ans[i + j * p] = norm_rand();
    }
    return ans;
}

/* These should be moved to other source files such as chm_common.c */
static cholmod_dense
*numeric_as_chm_dense(double *v, int n)
{
    cholmod_dense *ans = Calloc(1, cholmod_dense);
    
    ans->d = ans->nzmax = ans->nrow = n;
    ans->ncol = 1;
    ans->x = (void *) v;
    ans->xtype = CHOLMOD_REAL;
    ans->dtype = CHOLMOD_DOUBLE;
    return ans;
}

static
SEXP alloc_dgeMatrix(int m, int n, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, m * n);
    dn = GET_SLOT(ans, Matrix_DimNamesSym);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

static
SEXP alloc_dpoMatrix(int n, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * n);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    dn = GET_SLOT(ans, Matrix_DimNamesSym);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

static
SEXP alloc_dtrMatrix(int n, char *uplo, char *diag, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, n * n);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    dn = GET_SLOT(ans, Matrix_DimNamesSym);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

static
SEXP alloc_dsCMatrix(int n, int nz, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz);
    ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz);
    ALLOC_SLOT(ans, Matrix_pSym, INTSXP, n + 1);
    SET_SLOT(ans, Matrix_uploSym, mkString(uplo));
    dn = GET_SLOT(ans, Matrix_DimNamesSym);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

/* Internally used utilities */

static double *
cholmod_ata(cholmod_sparse *A, double val[], const char uplo[],
	    double alpha, double beta);
{
    int nnz = cholmod_nnz(A, &c);
    int *ai = (int*)(A->i), *ap = (int*)(A->p),
	*ind = Calloc(nnz, int), i, j, k, nr;
    double *ax = (double*)(A->x),
	*tmp = AZERO(Calloc(nnz * A->ncol, double), nnz * A->ncol);

    nr = ap[1];
    Memcpy(ind, ai, nr);
    Memcpy(tmp, ax, nr);
    for (j = 1; j < A->ncol; j++) {
	for (i = ap[j]; i < ap[j + 1]; i++) {
	    int ii = ai[i];
	    for (k = 0; k < nr; k++) {
		if (ii == ind[k]) {
		    tmp[k + j * nnz] = ax[i];
		    ii = -1;
		    break;
		}
	    }
	    if (ii >= 0) {	/* did not find the row index */
		ind[nr] = ii;
		tmp[nr + j * nnz] = ax[i];
		nr++;
	    }
	}
    }
    F77_CALL(dsyrk)(uplo, "T", &(A->ncol), &nr, &alpha, tmp, &nnz,
		    &beta, val, &(A->ncol));
    Free(ind); Free(tmp);
    return val;
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
 * Find a variable of a given name in a given environment and check
 * that its length and mode are correct.
 * 
 * @param rho Environment in which to find the variable
 * @param nm Name of the variable to find
 * @param mode Desired mode
 * @param len Desired length
 * 
 * @return 
 */
static
SEXP find_and_check(SEXP rho, SEXP nm, SEXPTYPE mode, int len)
{    
    SEXP ans;
    if (R_NilValue == PROTECT(ans = findVarInFrame(rho, nm)))
	error(_("environment `rho' must contain an object `%s'"),
	      CHAR(PRINTNAME(nm)));
    if (TYPEOF(ans) != mode)
	error(_("object `%s' of incorrect type"),
	      CHAR(PRINTNAME(nm)));
    if (len && LENGTH(ans) != len)
	error(_("object `%s' must be of length `%d'"),
	      CHAR(PRINTNAME(nm)), len);
    UNPROTECT(1);
    return ans;
}

/** 
 * Evaluate an expression in an environment, check that the length and
 * mode are as expected and store the result.
 * 
 * @param fcn expression to evaluate
 * @param rho environment in which to evaluate it
 * @param vv position to store the result
 * 
 * @return vv with new contents
 */
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

/** 
 * Evaluate an expression in an environment, check that the length and
 * mode are as expected and return the result.
 * 
 * @param fcn expression to evaluate
 * @param rho environment in which to evaluate it
 * @param mode desired mode
 * @param len desired length
 * 
 * @return evaluated expression
 */
static SEXP
eval_check(SEXP fcn, SEXP rho, SEXPTYPE mode, int len) {
    SEXP v = PROTECT(eval(fcn, rho));
    if (TYPEOF(v) != mode || LENGTH(v) != len)
	error(_("fcn produced mode %d, length %d - wanted mode %d, length %d"),
	      TYPEOF(v), LENGTH(v), mode, len);
    UNPROTECT(1);
    return v;
}

typedef struct glmer_struct
{
    SEXP cv;         /* control values */
    SEXP mer;	     /* mixed-effects representation */
    SEXP rho;        /* environment in which to evaluate the calls */
    SEXP eta;        /* linear predictor */
    SEXP mu;         /* mean vector */
    SEXP unwtd;      /* unweighted model matrices */
    SEXP wtd;        /* weighted model matrices */
    SEXP linkinv;    /* expression for inverse link evaluation */
    SEXP mu_eta;     /* expression for dmu/deta evaluation */
    SEXP LMEopt;     /* expression for LME optimization */
    SEXP dev_resids; /* expression for deviance residuals */
    SEXP var;        /* expression for variance evaluation */
    double *off;     /* offset for random effects only */
    double *offset;  /* offset for GLM */
    double *wts;     /* prior weights for GLM */
    double *w;	     /* vector of weights */
    double *X;	     /* copy of fixed-effects model matrix */
    double *y;       /* copy of response vector */
    double *z;       /* working residual */
    double *Zt;      /* copy of x slot in transpose of random-effects model matrix */
    double *etaold;  /* previous value of eta for evaluating  */
    int n;	     /* length of the response vector */
    int p;	     /* length of fixed effects vector */
    int nf;	     /* number of grouping factors */
    int npar;        /* total length of the parameter */
    int niterEM;     /* default number of ECME iterations */
    int EMverbose;   /* logical indicator */
    int maxiter;     /* maximum number of IRLS iterations */
    double tol;      /* convergence tolerance for IRLS iterations */
} glmer_struct, *GlmerStruct;

/** 
 * Increment the fitted values from an mer object using the (possibly
 * externally stored) contents of the model matrices.
 * 
 * @param x pointer to an mer object
 * @param X contents of the X model matrix or (double *) NULL
 * @param Zt values of non-zeros in the Zt model matrix or (double *) NULL
 * @param val array to hold the fitted values
 * 
 * @return pointer to a numeric array of fitted values
 */
static double *
internal_mer_fitted(SEXP x, double X[], double Ztx[], double val[])
{
    int ione = 1, n = LENGTH(GET_SLOT(x, Matrix_ySym));
    double one[] = {1,0};
    

    mer_secondary(x); 		/* has no effect if fixef is not 0 length */
    if (X) {
	SEXP fixef = GET_SLOT(x, Matrix_fixefSym);
	int p = LENGTH(fixef);

	F77_CALL(dgemv)("N", &n, &p, one, X, &n, REAL(fixef),
			&ione, one, val, &ione);
    }
    if (Ztx) {
	SEXP ranef = GET_SLOT(x, Matrix_ranefSym);
	cholmod_sparse *Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
	cholmod_dense *chv = numeric_as_chm_dense(val, n),
	    *chb = numeric_as_chm_dense(REAL(ranef), LENGTH(ranef));
	cholmod_sparse *Ztcp = cholmod_copy_sparse(Zt, &c);

	free(Zt);
	if (Ztx != (double *)(Zt->x)) Memcpy((double *)(Ztcp->x), Ztx, n);
	if (!cholmod_sdmult(Ztcp, 1, one, one, chb, chv, &c))
	    error(_("Error return from sdmult"));
	Free(chv); Free(chb); cholmod_free_sparse(&Ztcp, &c);
    }
    return val;
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
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, nf = length(Omega), vind;

    vind = 0;			/* index in ans */
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1;
	if (nci == 1) {
	    double dd = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym))[0];
	    ans[vind++] = ptyp ? ((ptyp == 1) ? log(dd) : 1./dd) : dd;
	} else {
	    if (ptyp) {	/* L log(D) L' factor of Omega[,,i] */
		int j, k, ncisq = nci * nci;
	        double *tmp = Memcpy(Calloc(ncisq, double),
				     REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)),
						   Matrix_xSym)), ncisq);
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
		double *omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym));

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
 * Set the coefficient vector and perform a factorization
 * 
 * @param x pointer to an mer object
 * @param cc vector of coefficients to assign
 * @param ptyp indicator of the parametrization being used
 */
static
void internal_mer_coefGets(SEXP x, const double cc[], int ptyp)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	cind, i, nf = length(Omega);

    cind = 0;
    for (i = 0; i < nf; i++) {
	SEXP Omegai = VECTOR_ELT(Omega, i);
	int j, k, nci = nc[i], ncisq = nc[i] * nc[i];
	double *choli, *omgi = REAL(GET_SLOT(Omegai, Matrix_xSym));

	if (nci == 1) {
	    double dd = cc[cind++];
	    *omgi = ptyp ? ((ptyp == 1) ? exp(dd) : 1./dd) : dd;
	} else {
	    int odind = cind + nci, /* off-diagonal index */
		ncip1 = nci + 1;

	    if (ptyp) {
		/* FIXME: Write this as an LDL decomposition */
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
	choli = REAL(GET_SLOT(dpoMatrix_chol(Omegai), Matrix_xSym));
	Memcpy(choli, omgi, ncisq);
	F77_CALL(dpotrf)("U", &nci, choli, &nci, &j);
	/* Yes, you really do need to do that decomposition.
	   The contents of choli before the decomposition are
	   from the previous value of Omegai. */
	if (j)
	    error(_("Omega[[%d]] is not positive definite"), i + 1);
    }
    mer_factor(x);
}

/** 
 * Evaluate current estimate of sigma from an mer object
 * 
 * @param x pointer to an mer object
 * @param REML indicator of whether to use REML.
 *           < 0  -> determine REML or ML from x@method
 *           == 0 -> use ML unconditionally
 *           > 0  -> use REML unconditionally
 * 
 * @return 
 */
static double
internal_mer_sigma(SEXP x, int REML)
{
    double *dcmp = REAL(GET_SLOT(x, Matrix_devCompSym));

    if (REML < 0)		/* get REML from x */
	REML = !strcmp(CHAR(asChar(GET_SLOT(x,
					    Matrix_methodSym))),
		       "REML");
/*     mer_factor(x); */ 
    return exp(dcmp[3]/2)/sqrt(dcmp[0] - (REML ? dcmp[1] : 0));
}

/** 
 * Update the relative precision matrices by sampling from a Wishart
 * distribution with scale factor determined by the current sample of
 * random effects.
 * 
 * @param Omega pointer to the list of relative precision matrices
 * @param b current sample from the random effects
 * @param sigma current value of sigma
 * @param nf number of grouping factors
 * @param nc number columns per grouping factor
 * @param Gp integer vector pointing to the beginning of each outer
 * group of columns
 * @param vals vector in which to store values
 * @param trans logical value indicating if variance components are
 * stored in the transformed scale.
 */
static void
internal_Omega_update(SEXP Omega, const double b[], double sigma, int nf,
		  const int nc[], const int Gp[], double *vals, int trans)
{
    int i, j, k, info;
    double one = 1, zero = 0;

    for (i = 0; i < nf; i++) {
	int nci = nc[i];
	int nlev = (Gp[i + 1] - Gp[i])/nci, ncip1 = nci + 1,
	    ncisqr = nci * nci;
	double
	    *scal = Calloc(ncisqr, double), /* factor of scale matrix */
	    *tmp = Calloc(ncisqr, double),
	    *var = Calloc(ncisqr, double), /* simulated variance-covariance */
	    *wfac = Calloc(ncisqr, double); /* factor of Wishart variate */

	/* generate and factor the scale matrix */
	AZERO(scal, ncisqr);
	F77_CALL(dsyrk)("U", "N", &nci, &nlev, &one, b + Gp[i], &nci,
			&zero, scal, &nci);
	F77_CALL(dpotrf)("U", &nci, scal, &nci, &info);
	if (info)
	    error(_("Singular random effects varcov at level %d"), i + 1);

	/* generate a random factor from a standard Wishart distribution */
	AZERO(wfac, ncisqr);
	std_rWishart_factor((double) (nlev - nci + 1), nci, wfac);

	/* form the variance-covariance matrix and store elements */
	Memcpy(tmp, scal, ncisqr);
	F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &nci,
			&one, wfac, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci,
			&zero, var, &nci);
	for (j = 0; j < nci; j++) {
	    *vals++ = (trans ? log(var[j * ncip1]) : var[j * ncip1]);
	}
	for (j = 1; j < nci; j++) {
	    for (k = 0; k < j; k++) {
		*vals++ = (trans ? atanh(var[k + j * nci]/
				     sqrt(var[j * ncip1] * var[k * ncip1]))
		       : var[k + j * nci]);
	    }
	}
	/* calculate and store the relative precision matrix */
	Memcpy(tmp, wfac, ncisqr);
	F77_CALL(dtrsm)("R", "U", "T", "N", &nci, &nci,
			&sigma, scal, &nci, tmp, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, tmp, &nci, &zero,
			REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym)),
			&nci);
	dpoMatrix_chol(VECTOR_ELT(Omega, i));
	Free(scal); Free(tmp); Free(wfac); Free(var); 
    }
}

static void
internal_weight_ZXy(SEXP mer, const double X[], const double Zt[],
		    const double w[], const double z[])
{
    SEXP mZt = GET_SLOT(mer, Matrix_ZtSym);
    int *Zp = INTEGER(GET_SLOT(mZt, Matrix_pSym));
    int i, j,
	n = LENGTH(GET_SLOT(mer, Matrix_ySym)), 
	p = LENGTH(GET_SLOT(mer, Matrix_rXySym));
    double *mX = REAL(GET_SLOT(mer, Matrix_XSym)),
	*my = REAL(GET_SLOT(mer, Matrix_ySym)),
	*mZx = REAL(GET_SLOT(mZt, Matrix_xSym));

    for (i = 0; i < n; i++) my[i] = z[i] * w[i];
    if (X) 
	for (j = 0; j < p; j++)
	    for (i = 0; i < n; i++) mX[i + j * n] = X[i + j * n] * w[i];
    if (Zt)
	for (j = 0; j < n; j++) {
	    int i2 = Zp[j + 1];
	    for (i = Zp[j]; i < i2; i++) mZx[i] = Zt[i] * w[j];
	}
}

/** 
 * Evaluate new weights and working residuals.
 * 
 * @param GS a GlmerStruct object
 */
static void
internal_weights(GlmerStruct GS) {
    SEXP dmu_deta, var;
    int i;
				/* reweight mer */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				  REALSXP, GS->n));
    var = PROTECT(eval_check(GS->var, GS->rho,
			     REALSXP, GS->n));
    for (i = 0; i < GS->n; i++) {
	GS->w[i] = GS->wts[i] *
	    REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	GS->z[i] = REAL(GS->eta)[i] - GS->offset[i] +
	    (GS->y[i] - REAL(GS->mu)[i])/REAL(dmu_deta)[i];
    }
    UNPROTECT(2);
}

/** 
 * Update eta, evaluate the convergence criterion, then copy eta to
 * etaold
 * 
 * @param GS a GlmerStruct object
 * @param etaold previous values of the linear predictors
 * 
 * @return convergence criterion
 */
static double
conv_crit(double etaold[], double eta[], int n) {
    double max_abs_eta = -1, max_abs_diff = -1;
    int i;

    for (i = 0; i < n; i++) {
	double abs_eta, abs_diff;

	abs_eta = fabs(eta[i]);
	if (abs_eta > max_abs_eta) max_abs_eta = abs_eta;
	abs_diff = fabs(eta[i] - etaold[i]);
	if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
	etaold[i] = eta[i];
    }
    return max_abs_diff / (0.1 + max_abs_eta);
}

/** 
 * Evaluate the quadratic form (x-mn)'A'A(x-mn) from the multivariate
 * normal kernel.
 * 
 * @param n dimension of random variate
 * @param mn mean vector
 * @param a upper Cholesky factor of the precision matrix
 * @param lda leading dimension of A
 * @param x vector at which to evaluate the kernel
 * 
 * @return value of the normal kernel
 */
static double
normal_kernel(int n, const double mn[],
	      const double a[], int lda, const double x[])
{
    int i, ione = 1;
    double *tmp = Calloc(n, double), ans;
    
    for (i = 0; i < n; i++) tmp[i] = x[i] - mn[i];
    F77_CALL(dtrmv)("U", "N", "N", &n, a, &lda, tmp, &ione);
    for (i = 0, ans = 0; i < n; i++) ans += tmp[i] * tmp[i];
    Free(tmp);
    return ans;
}  

/** 
 * Evaluate the conditional deviance for a given set of fixed effects.
 * 
 * @param GS Pointer to a GlmerStruct
 * @param fixed value of the fixed effects
 * 
 * @return conditional deviance
 */
static double
fixed_effects_deviance(GlmerStruct GS, const double fixed[])
{
    SEXP devs;
    int i, ione = 1;
    double ans, one = 1, zero = 0;

    F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
		    GS->X, &(GS->n), fixed,
		    &ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
    vecIncrement(REAL(GS->eta), GS->off, GS->n);
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
}

/** 
 * Establish off, the effective offset for the fixed effects, and
 * iterate to determine the conditional modes.
 * 
 * @param GS a GlmerStruct object
 * @param fixed vector of fixed effects
 * @param varc vector of parameters for the variance-covariance
 * 
 * @return An indicator of whether the iterations converged
 */
static int
internal_bhat(GlmerStruct GS, const double fixed[], const double varc[])
{
    int i, ione = 1;
    double crit, one = 1;

    if (varc)	  /* skip this step if varc == (double*) NULL */	
	internal_mer_coefGets(GS->mer, varc, 2);

    Memcpy(GS->off, GS->offset, GS->n);
    F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
		    GS->X, &(GS->n), fixed,
		    &ione, &one, GS->off, &ione);
    Memcpy(REAL(GS->eta), GS->off, GS->n);
    Memcpy(GS->etaold, REAL(GS->eta), GS->n);

    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	internal_weights(GS);
	internal_weight_ZXy(GS->mer, (double *) NULL, GS->Zt, GS->w, GS->z);
	mer_update_ZXy(GS->mer);
	mer_factor(GS->mer);
	mer_secondary(GS->mer);
	Memcpy(REAL(GS->eta), GS->off, GS->n);
	internal_mer_fitted(GS->mer, (double *) NULL, GS->Zt, REAL(GS->eta));
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
    }
    return (crit > GS->tol) ? 0 : 1;
}

/**
 * Fill in the 4-dimensional vector of linear combinations of the
 * gradComp array according to whether ECME steps or the gradient are
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
 */
static void
EMsteps_verbose_print(SEXP x, int iter, int REML)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	gradComp = GET_SLOT(x, Matrix_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ifour = 4, ii, ione = 1, jj, nf = LENGTH(Omega);
    double
	*cc = EM_grad_lc(Calloc(4, double), 0, REML, nc + nf),
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	one = 1., zero = 0.;

    mer_factor(x);
    if (iter == 0) Rprintf("  EM iterations\n");
    Rprintf("%3d %.3f", iter, dev[REML ? 1 : 0]);
    for (i = 0; i < nf; i++) {
	int nci = nc[i], ncip1 = nci + 1, ncisqr = nci * nci;
	double
	    *Omgi = REAL(VECTOR_ELT(Omega, i)),
	    *Grad = Calloc(ncisqr, double);

				/* diagonals */
	Rprintf(" (%#8g", Omgi[0]);
	for (jj = 1; jj < nci; jj++) {
	    Rprintf(" %#8g", Omgi[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++)
		Rprintf(" %#8g", Omgi[ii + jj * nci]);
				/* Evaluate and print the gradient */
	F77_CALL(dgemv)("N", &ncisqr, &ifour, &one,
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			cc, &ione, &zero, Grad, &ione);
	Rprintf(":%#8.3g", Grad[0]);
				
	for (jj = 1; jj < nci; jj++) { /* diagonals */
	    Rprintf(" %#8.3g", Grad[jj * ncip1]);
	}
	for (jj = 1; jj < nci; jj++) /* offdiagonals */
	    for (ii = 0; ii < jj; ii++)
		Rprintf(" %#8.3g", Grad[ii + jj * nci]);
	Rprintf(")");
	Free(Grad);
    }
    Rprintf("\n");
    Free(cc);
}

/** 
 * Perform a number of ECME steps
 * 
 * @param x pointer to an mer object
 * @param nEM number of iteration to perform
 * @param verb indicator of verbose output
 */
static void
internal_ECMEsteps(SEXP x, int nEM, int verb)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym),
	flist = GET_SLOT(x, Matrix_flistSym),
	gradComp = GET_SLOT(x, Matrix_gradCompSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),
		       "REML"),
	i, ifour = 4, info, ione = 1, iter,
	nf = LENGTH(Omega);
    double
	*cc = EM_grad_lc(Calloc(4, double), 1, REML, nc + nf),
	zero = 0.0;

/* FIXME: This is currently a stub. */
    return;
    mer_gradComp(x);
    if (verb)
	EMsteps_verbose_print(x, 0, REML);
    for (iter = 0; iter < nEM; iter++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], ncisqr = nci * nci;
	    double *Omgi = REAL(VECTOR_ELT(Omega, i)),
		mult = 1./
		((double) length(getAttrib(VECTOR_ELT(flist, i),
				 R_LevelsSymbol)));

	    F77_CALL(dgemm)("N", "N", &ncisqr, &ione, &ifour, &mult,
			    REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
			    cc, &ifour, &zero, Omgi, &ncisqr);
	    F77_CALL(dpotrf)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error(_("DPOTRF in ECME update gave code %d"),
		      info);
	    F77_CALL(dpotri)("U", &nci, Omgi, &nci, &info);
	    if (info)
		error(_("Matrix inverse in ECME update gave code %d"), info);
	}
	mer_gradComp(x);
	if (verb)
	    EMsteps_verbose_print(x, iter + 1, REML);
    }
    mer_factor(x);
    Free(cc);
    UNPROTECT(1);
}

static double chm_log_abs_det(cholmod_factor *F)
{
    double ans = 0;

    if (F->is_super) {
	int i;
	for (i = 0; i < F->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(F->pi))[i + 1] - ((int *)(F->pi))[i],
		nc = ((int *)(F->super))[i + 1] - ((int *)(F->super))[i];
	    double *x = (double *)(F->x) + ((int *)(F->px))[i];

	    for (j = 0; j < nc; j++) ans += log(fabs(x[j * nrp1]));
	}
    } else
	error(_("code for simplicial factorization not yet written"));
    return ans;
}

static double Omega_log_det(SEXP Omega, int nf, int *nc, int *Gp)
{
    double ans = 0;
    int i;

    for (i = 0; i < nf; i++) {
	int j, nci = nc[i], ncip1 = nc[i] + 1, nlev = (Gp[i + 1] - Gp[i])/nc[i];
	double *omgi = REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				     Matrix_xSym));

	for (j = 0; j < nci; j++) ans += 2. * nlev * log(fabs(omgi[j * ncip1]));
    }
    return ans;
}

/** 
 * Determine the deviance components associated with each of the
 * levels of a grouping factor at the conditional modes or a value
 * offset from the conditional modes by delb.
 * 
 * @param GS pointer to a GlmerStruct
 * @param b pointer to the conditional modes of the random effects 
 * @param nlev number of levels of the kth grouping factor
 * @param nc number of columns in the model matrix for the kth
 * grouping factor
 * @param k index (0-based) of the grouping factor
 * @param delb vector of length nc giving the changes in the
 * orthonormalized random effects 
 * @param OmgFac Cholesky factor of the inverse of the penalty matrix
 * for this grouping factor 
 * @param bVfac 3-dimensional array holding the factors of the
 * conditional variance-covariance matrix of the random effects 
 * @param devcmp array to hold the deviance components
 * 
 * @return devcmp
 */
static double*
rel_dev_1(GlmerStruct GS, SEXP b, int nlev, int nc, int k,
	  const double delb[], const double OmgFac[],
	  const double bVfac[], double devcmp[])
{
    SEXP devs, flist = GET_SLOT(GS->mer, Matrix_flistSym);
    int *fv = INTEGER(VECTOR_ELT(flist, k)), i, ione = 1,
	j, ntot = nlev * nc;
    double *bb = REAL(VECTOR_ELT(b, k)), *bcp = (double *) NULL;

    AZERO(devcmp, nlev);
    if (delb) {
	double sumsq = 0;
				/* copy the contents of b */
	bcp = Memcpy(Calloc(ntot, double), bb, ntot);
	if (nc == 1) {
	    sumsq = delb[0] * delb[0];
	    for (i = 0; i < nlev; i++) bb[i] += delb[0] * bVfac[i];
	} else {
	    int ncsq = nc * nc;
	    double *tmp = Calloc(nc, double);
	    for (i = 0; i < nlev; i++) {
		Memcpy(tmp, delb, nc);
		F77_CALL(dtrmv)("U", "N", "N", &nc, &(bVfac[i * ncsq]),
				&nc, tmp, &ione);
		for (j = 0; j < nc; j++) bb[i + j * nc] = tmp[j];
	    }
				/* sum of squares of delb */
	    for (j = 0; j < nc; j++) sumsq += delb[j] * delb[j];
	}
	for (i = 0; i < nlev; i++) devcmp[i] = -sumsq;
    }
    Memcpy(REAL(GS->eta), GS->off, GS->n);
/*     fitted_ranef(flist, GS->unwtd, b, &nc, REAL(GS->eta)); */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0; i < GS->n; i++)
	devcmp[fv[i] - 1] += REAL(devs)[i];
    UNPROTECT(1);
    if (nc == 1) {
	for (i = 0; i < nlev; i++) {
	    double tmp = *OmgFac * bb[i];
	    devcmp[i] += tmp * tmp;
	}
    } else {
	double *tmp = Calloc(nc, double);
	int ione = 1;
	
	for (i = 0; i < nlev; i++) {
	    for (j = 0; j < nc; j++) tmp[j] = bb[i + j * nlev];
	    F77_CALL(dtrmv)("U", "N", "N", &nc, OmgFac, &nc,
			    tmp, &ione);
	    for (j = 0; j < nc; j++) 
		devcmp[i] += tmp[j] * tmp[j];
	}
    }
    if (delb) {
	Memcpy(bb, bcp, ntot);
	Free(bcp);
    }
    return devcmp;
}

/** 
 * Create a copy of ZtZ with the diagonal blocks inflated according to Omega
 * 
 * @param zz cholmod_sparse version of ZtZ
 * @param nf number of factors
 * @param Omega Omega list
 * @param nc number of columns in model matrices
 * @param Gp group pointers
 * 
 * @return a freshly allocated cholmod_sparse version of the sum
 */
static cholmod_sparse *
ZZ_inflate(cholmod_sparse *zz, int nf, SEXP Omega, int *nc, int *Gp)
{
    cholmod_sparse *Omg, *ans;
    int *omp, *nnz = Calloc(nf + 1, int), i;
    double one = 1;

    for (nnz[0] = 0, i = 0; i < nf; i++)
	nnz[i + 1] = nnz[i] + (Gp[i + 1] - Gp[i])*(nc[i] + 1)/2;
    Omg = cholmod_allocate_sparse(zz->nrow, zz->ncol, (size_t) nnz[nf],
				  TRUE, TRUE, 1, CHOLMOD_REAL, &c);
    omp = (int *) Omg->p;
    for (i = 0; i < nf; i++) {
	int bb = Gp[i], j, jj, k, nci = nc[i];
	int nlev = (Gp[i + 1] - bb)/nci;
	double *Omgi = REAL(GET_SLOT(VECTOR_ELT(Omega, i), Matrix_xSym));

	for (j = 0; j < nlev; j++) { /* column of result */
	    int col0 = bb + j * nci; /* absolute column number */

	    for (jj = 0; jj < nci; jj++) { /* column of Omega_i */
		int coljj = col0 + jj;

		omp[coljj + 1] = omp[coljj] + jj + 1;
		for (k = 0; k <= jj; k++) { /* row of Omega_i */
		    int ind = omp[coljj];
		    ((int *)Omg->i)[ind + k] = col0 + k;
		    ((double *)Omg->x)[ind + k] = Omgi[jj * nci + k];
		}
	    }
	}
    }
    ans = cholmod_add(zz, Omg, &one, &one, TRUE, TRUE, &c);

    Free(nnz); cholmod_free_sparse(&Omg, &c);
    return ans;
}

/** 
 * Evaluate the conditional deviance for a given set of fixed effects.
 * 
 * @param GS Pointer to a GlmerStruct
 * @param fixed value of the fixed effects
 * 
 * @return conditional deviance
 */
static double
random_effects_deviance(GlmerStruct GS, SEXP b)
{
    SEXP devs;
    int i;
    double ans;

    Memcpy(REAL(GS->eta), GS->off, GS->n);
/*     fitted_ranef(GET_SLOT(GS->mer, Matrix_flistSym), GS->unwtd, b, */
/* 		 INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)), REAL(GS->eta)); */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n));
    for (i = 0, ans = 0; i < GS->n; i++) ans += REAL(devs)[i];
    UNPROTECT(1);
    return ans;
}

static void
internal_glmer_ranef_update(GlmerStruct GS, SEXP b)
{
    SEXP bhat, bprop = PROTECT(duplicate(b)), 
	bVar = GET_SLOT(GS->mer, Matrix_bVarSym),
	Omega = GET_SLOT(GS->mer, Matrix_OmegaSym);
    int i, ione = 1, j, k;
    double devr, one = 1;

/*     bhat = PROTECT(lmer_ranef(GS->mer)); */
				/* subtract deviance at b */
    devr = -random_effects_deviance(GS, b);
    for (i = 0; i < GS->nf; i++) {
	SEXP Bi = VECTOR_ELT(b, i);
	int *dims = INTEGER(getAttrib(Bi, R_DimSymbol));
	int nlev = dims[0], nci = dims[1];
	int ncisqr = nci * nci, ntot = nlev * nci;
	double *bcopy = Calloc(ntot, double),
	    *bi = REAL(Bi),
	    *bhati = REAL(VECTOR_ELT(bhat, i)),
	    *bpropi = REAL(VECTOR_ELT(bprop, i)),
	    *bVari = REAL(VECTOR_ELT(bVar, i)),
	    *chol = Calloc(ncisqr, double),
	    *delta = Calloc(nci, double),
	    *omgfac = Memcpy(Calloc(ncisqr, double),
			     REAL(VECTOR_ELT(Omega, i)),
			     ncisqr);

				/* subtract quadratic form in
				 * Omega[[i]] at b  */
	F77_CALL(dpotrf)("U", &nci, omgfac, &nci, &j);
	if (j)
	    error(_("Leading %d minor of Omega[[%d]] not positive definite"),
                      j, i + 1);
	Memcpy(bcopy, bi, ntot);
	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
			omgfac, &nci, bcopy, &nlev);
	for (k = 0; k < ntot; k++) devr -= bcopy[k] * bcopy[k];
				/* form bprop and proposal density */
	for (k = 0; k < nlev; k++) {
				/* proposal density at b */
	    for (j = 0; j < nci; j++)
		delta[j] = bi[k + j * nlev] - bhati[k + j * nlev];
	    Memcpy(chol, &(bVari[k * ncisqr]), ncisqr);
	    F77_CALL(dpotrf)("U", &nci, chol, &nci, &j);
	    if (j)
		error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
	    F77_CALL(dtrsv)("U", "T", "N", &nci, chol, &nci,
			    delta, &ione);
	    for (j = 0; j < nci; j++) {
		double nrm = norm_rand(); /* proposal deviate */
		devr += delta[j] * delta[j] - nrm * nrm;
		delta[j] = nrm;
	    }
				/* scale by Cholesky inverse */
	    F77_CALL(dtrmv)("U", "T", "N", &nci, chol, &nci,
			    delta, &ione);
				/* add mean */
	    for (j = 0; j < nci; j++)
		bpropi[k + j * nlev] = bhati[k + j * nlev] + delta[j];
	}
				/* add quadratic form in
				 * Omega[[i]] at bprop  */
	Memcpy(bcopy, bpropi, ntot);
	F77_CALL(dtrmm)("R", "U", "T", "N", &nlev, &nci, &one,
			omgfac, &nci, bcopy, &nlev);
	for (k = 0; k < ntot; k++) devr += bcopy[k] * bcopy[k];

	Free(bcopy); Free(chol); Free(delta); Free(omgfac);
    }
				/* add deviance at bprop */
    devr += random_effects_deviance(GS, bprop);

    if (unif_rand() < exp(-0.5 * devr))
	for (i = 0; i < GS->nf; i++) { /* copy each face of b */
	    SEXP Bi = VECTOR_ELT(b, i);
	    int *dims = INTEGER(getAttrib(Bi, R_DimSymbol));

	    Memcpy(REAL(Bi), REAL(VECTOR_ELT(bprop, i)),
		   dims[0] * dims[1]);
	}
    
    if (asLogical(Matrix_getElement(GS->cv, "msVerbose"))) {
	double *b0 = REAL(VECTOR_ELT(bprop, 0));
	Rprintf("%5.3f:", exp(-0.5 * devr));
	for (k = 0; k < 5; k++) Rprintf("%#10g ", b0[k]);
	Rprintf("\n");
    }
	
    UNPROTECT(2);
}

    
/** 
 * Determine the conditional modes and the conditional variance of the
 * fixed effects given the data and the current random effects.
 * Create a Metropolis-Hasting proposal step from the multivariate
 * normal density, determine the acceptance probability and, if the
 * step is to be accepted, overwrite the contents of fixed with the
 * new contents.
 * 
 * @param GS a GlmerStruct
 * @param b list of random effects
 * @param fixed current value of the fixed effects
 * 
 * @return updated value of the fixed effects
 */
static double *
internal_glmer_fixef_update(GlmerStruct GS, SEXP b,
			    double fixed[])
{
    SEXP dmu_deta, var;
    int i, ione = 1, it, j, lwork = -1;
    double *ans = Calloc(GS->p, double), /* proposal point */
	*md = Calloc(GS->p, double), /* conditional modes */
	*w = Calloc(GS->n, double), *work,
	*wtd = Calloc(GS->n * GS->p, double),
	*z = Calloc(GS->n, double),
	crit, devr, one = 1, tmp, zero = 0;
    
    if (!isNewList(b) || LENGTH(b) != GS->nf)
	error(_("%s must be a %s of length %d"), "b", "list", GS->nf);
    for (i = 0; i < GS->nf; i++) {
	SEXP bi = VECTOR_ELT(b, i);
	if (!isReal(bi) || !isMatrix(bi))
	    error(_("b[[%d]] must be a numeric matrix"), i);
    }
    AZERO(z, GS->n);		/* -Wall */
    Memcpy(md, fixed, GS->p);
				/* calculate optimal size of work array */
    F77_CALL(dgels)("N", &(GS->n), &(GS->p), &ione, wtd, &(GS->n),
		    z,  &(GS->n), &tmp, &lwork, &j);
    if (j)			/* shouldn't happen */
	error(_("%s returned error code %d"), "dgels", j);
    lwork = (int) tmp;
    work = Calloc(lwork, double);
				
    AZERO(GS->off, GS->n); /* fitted values from random effects */
/*     fitted_ranef(GET_SLOT(GS->mer, Matrix_flistSym), GS->unwtd, b, */
/* 		 INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)), GS->off); */
    for (i = 0; i < GS->n; i++)
	(GS->etaold)[i] = ((GS->off)[i] += (GS->offset)[i]);
    
    for (it = 0, crit = GS->tol + 1;
	 it < GS->maxiter && crit > GS->tol; it++) {
				/* fitted values from current beta */
	F77_CALL(dgemv)("N", &(GS->n), &(GS->p), &one,
			GS->X, &(GS->n), md,
			&ione, &zero, REAL(GS->eta), &ione);
				/* add in random effects and offset */
	vecIncrement(REAL(GS->eta), (GS->off), GS->n);
				/* check for convergence */
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
				/* obtain mu, dmu_deta, var */
	eval_check_store(GS->linkinv, GS->rho, GS->mu);
	dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				      REALSXP, GS->n));
	var = PROTECT(eval_check(GS->var, GS->rho, REALSXP, GS->n));
				/* calculate weights and working residual */
	for (i = 0; i < GS->n; i++) {
	    w[i] = GS->wts[i] * REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	    z[i] = w[i] * (REAL(GS->eta)[i] - (GS->off)[i] +
			   ((GS->y)[i] - REAL(GS->mu)[i]) /
			   REAL(dmu_deta)[i]);
	}
	UNPROTECT(2);
				/* weighted copy of the model matrix */
	for (j = 0; j < GS->p; j++)
	    for (i = 0; i < GS->n; i++)
		wtd[i + j * GS->n] = GS->X[i + j * GS->n] * w[i];
				/* weighted least squares solution */
	F77_CALL(dgels)("N", &(GS->n), &(GS->p), &ione, wtd, &(GS->n),
			z, &(GS->n), work, &lwork, &j);
	if (j) error(_("%s returned error code %d"), "dgels", j);
	Memcpy(md, z, GS->p);
    }
				/* wtd contains the Cholesky factor of
				 * the precision matrix */
    devr = normal_kernel(GS->p, md, wtd, GS->n, fixed);
    devr -= fixed_effects_deviance(GS, fixed);
    for (i = 0; i < GS->p; i++) {
	double var = norm_rand();
	ans[i] = var;
	devr -= var * var;
    }
    F77_CALL(dtrsv)("U", "N", "N", &(GS->p), wtd, &(GS->n), ans, &ione);
    for (i = 0; i < GS->p; i++) ans[i] += md[i];
    devr += fixed_effects_deviance(GS, ans);
    crit = exp(-0.5 * devr);	/* acceptance probability */
    tmp = unif_rand();
    if (asLogical(Matrix_getElement(GS->cv, "msVerbose"))) {
	Rprintf("%5.3f: ", crit);
	for (j = 0; j < GS->p; j++) Rprintf("%#10g ", ans[j]);
	Rprintf("\n");
    }
    if (tmp < crit) Memcpy(fixed, ans, GS->p);
    Free(ans); Free(md); Free(w);
    Free(work); Free(wtd); Free(z);
    return fixed;
}

/* Gauss-Hermite Quadrature x positions and weights */
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
    GHQ_w5[3] = {0.0112574109895360, 0.222075915334214,
		 0.533333317311434},
    GHQ_x6[3] = {3.32425737665988, 1.88917584542184,
		 0.61670657963811},
    GHQ_w6[3] = {0.00255578432527774, 0.0886157433798025,
		 0.408828457274383},
    GHQ_x7[4] = {3.7504396535397, 2.36675937022918,
		 1.15440537498316, 0},
    GHQ_w7[4] = {0.000548268839501628, 0.0307571230436095,
		 0.240123171391455, 0.457142843409801},
    GHQ_x8[4] = {4.14454711519499, 2.80248581332504,
		 1.63651901442728, 0.539079802125417},
    GHQ_w8[4] = {0.000112614534992306, 0.0096352198313359,
		 0.117239904139746, 0.373012246473389},
    GHQ_x9[5] = {4.51274578616743, 3.20542894799789,
		 2.07684794313409, 1.02325564627686, 0},
    GHQ_w9[5] = {2.23458433364535e-05, 0.00278914123744297,
		 0.0499164052656755, 0.244097495561989,
		 0.406349194142045},
    GHQ_x10[5] = {4.85946274516615, 3.58182342225163,
		  2.48432579912153, 1.46598906930182,
		  0.484935699216176},
    GHQ_w10[5] = {4.31065250122166e-06, 0.000758070911538954,
		  0.0191115799266379, 0.135483698910192,
		  0.344642324578594},
    GHQ_x11[6] = {5.18800113558601, 3.93616653976536,
		  2.86512311160915, 1.87603498804787,
		  0.928868981484148, 0},
    GHQ_w11[6] = {8.12184954622583e-07, 0.000195671924393029,
		  0.0067202850336527, 0.066138744084179,
		  0.242240292596812, 0.36940835831095};

static const double
    *GHQ_x[12] = {(double *) NULL, GHQ_x1, GHQ_x2, GHQ_x3, GHQ_x4,
		  GHQ_x5, GHQ_x6, GHQ_x7, GHQ_x8, GHQ_x9, GHQ_x10,
		  GHQ_x11},
    *GHQ_w[12] = {(double *) NULL, GHQ_w1, GHQ_w2, GHQ_w3, GHQ_w4,
		  GHQ_w5, GHQ_w6, GHQ_w7, GHQ_w8, GHQ_w9, GHQ_w10,
		  GHQ_w11};

/* Externally accessible functions */

/** 
 * Simulate a sample of random matrices from a Wishart distribution
 * 
 * @param ns Number of samples to generate
 * @param dfp Degrees of freedom
 * @param scal Positive-definite scale matrix
 * 
 * @return 
 */
SEXP
Matrix_rWishart(SEXP ns, SEXP dfp, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), j,
	n = asInteger(ns), psqr;
    double *scCp, *tmp, df = asReal(dfp), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Calloc(psqr, double);
    AZERO(tmp, psqr);
    scCp = Memcpy(Calloc(psqr, double), REAL(scal), psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &j);
    if (j)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3Darray(REALSXP, dims[0], dims[0], n));
  
    GetRNGstate();
    for (j = 0; j < n; j++) {
	double *ansj = REAL(ans) + j * psqr;
	std_rWishart_factor(df, dims[0], tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));
	internal_symmetrize(ansj, dims[0]);
	}

    PutRNGstate();
    Free(scCp); Free(tmp);
    UNPROTECT(1);
    return ans;
}

/** 
 * Perform the PQL optimization
 * 
 * @param GSp pointer to a GlmerStruct object
 * 
 * @return R_NilValue
 */
SEXP glmer_PQL(SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int i; double crit;

    Memcpy(GS->etaold, REAL(GS->eta), GS->n);
    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	internal_weights(GS);
	internal_weight_ZXy(GS->mer, GS->X, GS->Zt, GS->w, GS->z);
	mer_update_ZXy(GS->mer);
	if (!i) mer_initial(GS->mer); /* initialize first fit */
	mer_factor(GS->mer);
	mer_secondary(GS->mer);
	internal_ECMEsteps(GS->mer, i ? 2 : GS->niterEM,
			   GS->EMverbose);
	eval(GS->LMEopt, GS->rho);
	Memcpy(REAL(GS->eta), GS->offset, GS->n);
	internal_mer_fitted(GS->mer, GS->X, GS->Zt, REAL(GS->eta));
	crit = conv_crit(GS->etaold, REAL(GS->eta), GS->n);
    }
    if (crit > GS->tol)
	warning(_("IRLS iterations for PQL did not converge"));

    return R_NilValue;
}

/** 
 * Compute the approximation to the deviance using adaptive
 * Gauss-Hermite quadrature (AGQ).  When nAGQ == 1 this is the Laplace
 * approximation.
 * 
 * @param pars pointer to a numeric vector of parameters
 * @param GSp pointer to a GlmerStruct object
 * @param nAGQp pointer to a scalar integer representing the number of
 * points in AGQ to use
 * 
 * @return the approximation to the deviance as computed using AGQ
 */
SEXP glmer_devAGQ(SEXP pars, SEXP GSp, SEXP nAGQp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    SEXP Omega, bVar;
    int i, j, k, nAGQ = asInteger(nAGQp);
    int n2 = (nAGQ + 1)/2;
    double *f0, LaplaceDev = 0, AGQadjst = 0,
	*bhat = REAL(GET_SLOT(GS->mer, Matrix_ranefSym));
	
    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (GS->nf > 1 && nAGQ > 1) {
	warning(_("AGQ not available for multiple grouping factors - using Laplace"));
	nAGQ = 1;
    }
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);
    return R_NilValue;
    bVar = GET_SLOT(GS->mer, Matrix_bVarSym);
    Omega = GET_SLOT(GS->mer, Matrix_OmegaSym);

    for (i = 0; i < GS->nf; i++) {
	int *dims = INTEGER(getAttrib(VECTOR_ELT(bVar, i),
				      R_DimSymbol));
	int nci = dims[0], nlev = dims[2];
	int ncip1 = nci + 1, ncisqr = nci * nci;
	double *omg = Memcpy(Calloc(ncisqr, double),
			     REAL(VECTOR_ELT(Omega, i)), ncisqr),
	    *bvar = Memcpy(Calloc(ncisqr * nlev, double),
			   REAL(VECTOR_ELT(bVar, i)), ncisqr * nlev);

	LaplaceDev = 0;
				/* Calculate difference of determinants */
        F77_CALL(dpotrf)("U", &nci, omg, &nci, &j);
        if (j)
            error(_("Leading %d minor of Omega[[%d]] not positive definite"),
                  j, i + 1);
        for (j = 0; j < nci; j++) { /* nlev * logDet(Omega_i) */
            /* Note: we subtract because Omega has been inverted */
            LaplaceDev -= 2 * nlev * log(omg[j * ncip1]);
        }
        for (k = 0; k < nlev; k++) {
	    double *bVik = bvar + k * ncisqr;
            F77_CALL(dpotrf)("U", &nci, bVik, &nci, &j);
            if (j)
                error(_("Leading %d minor of bVar[[%d]][,,%d] not positive definite"),
                      j, i + 1, k + 1);
            for (j = 0; j < nci; j++) LaplaceDev -= 2 * log(bVik[j * ncip1]);
        }

	f0 = Calloc(nlev, double);
	rel_dev_1(GS, bhat, nlev, nci, i, (double *) NULL,
		  omg, bvar, f0);
	for (k = 0; k < nlev; k++) LaplaceDev += f0[k];
	if (nAGQ > 1) {
	    double *fx = Calloc(nlev, double),
		*rellik = Calloc(nlev, double),
		*delb = Calloc(nci, double);

	    if (nci > 1) error(_("code not yet written"));
	    AZERO(rellik, nlev);	/* zero accumulator */
	    for (k = 0; k < n2; k++) {	
		delb[0] = GHQ_x[nAGQ][k];
		if (delb[0]) {
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omg, bvar, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		    delb[0] *= -1;
		    rel_dev_1(GS, bhat, nlev, nci, i, delb,
			      omg, bvar, fx);
		    for (j = 0; j < nlev; j++) {
			rellik[j] += GHQ_w[nAGQ][k] *
			    exp(-(fx[j] - f0[j])/2);
		    }
		} else {
		    for (j = 0; j < nlev; j++)
			rellik[j] += GHQ_w[nAGQ][k];
		}
	    }
	    for (j = 0; j < nlev; j++)
		AGQadjst -= 2 * log(rellik[j]);
	    Free(fx); Free(rellik);
	}
	Free(f0); Free(omg); Free(bvar);
    }
    UNPROTECT(1);
    return ScalarReal(LaplaceDev + AGQadjst);
}

/** 
 * Release the storage for a GlmerStruct
 * 
 * @param GSp External pointer to a  GlmerStruct
 * 
 * @return R_NilValue
 */
SEXP glmer_finalize(SEXP GSp) {
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    
    Free(GS->w); Free(GS->X); Free(GS->y); Free(GS->z);
    Free(GS->Zt); Free(GS->off); Free(GS->offset); Free(GS->wts);
    Free(GS->etaold);
    Free(GS);
    return R_NilValue;
}

/** 
 * Determine the conditional modes and the conditional variance of the
 * fixed effects given the data and the current random effects.
 * Create a Metropolis-Hasting proposal step from the multivariate
 * normal density, determine the acceptance probability and return the
 * current value or the proposed value.
 * 
 * @param GSp pointer to a GlmerStruct
 * @param b list of random effects
 * @param fixed current value of the fixed effects
 * 
 * @return updated value of the fixed effects
 */
SEXP glmer_fixed_update(SEXP GSp, SEXP b, SEXP fixed)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
 
    if (!isReal(fixed) || LENGTH(fixed) != GS->p)
	error(_("%s must be a %s of length %d"), "fixed",
		"numeric vector", GS->p);
    GetRNGstate();
    internal_glmer_fixef_update(GS, b, REAL(fixed));
    PutRNGstate();
    return fixed;
}

/** 
 * Return an external pointer object to a GlmerStruct created in
 * environment rho
 * 
 * @param rho An environment
 * 
 * @return An external pointer to a GlmerStruct
 */
SEXP glmer_init(SEXP rho) {
    GlmerStruct GS;
    SEXP tmp, y, Ztx;
    
    
    GS = (GlmerStruct) Calloc(1, glmer_struct);
    if (!isEnvironment(rho))
	error(_("`rho' must be an environment"));
    GS->rho = rho;
    GS->mer = find_and_check(rho, install("mer"), VECSXP, 0);
    y = GET_SLOT(GS->mer, Matrix_ySym);
    GS->n = LENGTH(y);
    GS->p = LENGTH(GET_SLOT(GS->mer, Matrix_rXySym));
    GS->X = Memcpy(Calloc(GS->n * GS->p, double),
		   REAL(GET_SLOT(GS->mer, Matrix_XSym)), GS->n * GS->p);
    GS->y = Memcpy(Calloc(GS->n, double), REAL(y), GS->n);
    GS->w = Calloc(GS->n, double);
    GS->z = Calloc(GS->n, double);
    Ztx = GET_SLOT(GET_SLOT(GS->mer, Matrix_ZtSym), Matrix_xSym);
    GS->Zt = Memcpy(Calloc(LENGTH(Ztx), double), REAL(Ztx), LENGTH(Ztx));
    GS->mu = find_and_check(rho, install("mu"), REALSXP, GS->n);
    tmp = find_and_check(rho, install("offset"), REALSXP, GS->n);
    GS->offset = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    tmp = find_and_check(rho, install("wts"), REALSXP, GS->n);
    GS->wts = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    GS->off = Calloc(GS->n, double);
    GS->etaold = Calloc(GS->n, double);
    GS->cv = find_and_check(rho, install("cv"), VECSXP, 0);
    GS->niterEM = asInteger(Matrix_getElement(GS->cv, "niterEM"));
    GS->EMverbose = asLogical(Matrix_getElement(GS->cv, "EMverbose"));
    GS->tol = asReal(Matrix_getElement(GS->cv, "tolerance"));
    GS->maxiter = asInteger(Matrix_getElement(GS->cv, "maxIter"));
    GS->nf = LENGTH(GET_SLOT(GS->mer, Matrix_flistSym));
    GS->npar = GS->p +
	coef_length(GS->nf, INTEGER(GET_SLOT(GS->mer, Matrix_ncSym)));
    GS->eta = find_and_check(rho, install("eta"), REALSXP, GS->n);

    GS->linkinv = find_and_check(rho, install("linkinv"),
				 LANGSXP, 0);
    GS->mu_eta = find_and_check(rho, install("mu.eta"),
				LANGSXP, 0);
    GS->var = find_and_check(rho, install("variance"),
			     LANGSXP, 0);
    GS->LMEopt = find_and_check(rho, install("doLMEopt"),
				LANGSXP, 0);
    GS->dev_resids = find_and_check(rho, install("dev.resids"),
				    LANGSXP, 0);
    return R_MakeExternalPtr(GS, R_NilValue, GS->mer);
}

/** 
 * Determine the conditional modes and the conditional variance of the
 * random effects given the data and the current fixed effects and
 * variance components.
 *
 * Create a Metropolis-Hasting proposal step from a multivariate
 * normal density centered at bhat with variance-covariance matrix
 * from bVar, determine the acceptance probability and return the
 * current value or the proposed value.
 * 
 * @param GSp pointer to a GlmerStruct
 * @param fixed pointer to a numeric vector of the fixed effects
 * @param varc pointer to a numeric vector of the variance components
 * @param varc pointer to current values of b
 * 
 * @return updated b from the Metropolis-Hastings step
 */
SEXP glmer_ranef_update(SEXP GSp, SEXP fixed, SEXP varc, SEXP b)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    int nvarc = GS->npar - GS->p;

    if (!isReal(fixed) || LENGTH(fixed) != GS->p)
	error(_("`%s' must be a numeric vector of length %d"),
	      "fixed", GS->p);
    if (INTEGER(GET_SLOT(GS->mer, Matrix_ncSym))[GS->nf] > 0)
	error(_("the mer object must be set to skip fixed effects"));
    if (!isReal(varc) || LENGTH(varc) != nvarc)
	error(_("`%s' must be a numeric vector of length %d"),
	      "varc", nvarc);

    GetRNGstate();
    /* Don't check for convergence failure after internal_bhat.
     * It is determining the mean of the proposal density and
     * does not need exact convergence. */
    internal_bhat(GS, REAL(fixed), REAL(varc));
    internal_glmer_ranef_update(GS, b);
    PutRNGstate();

    UNPROTECT(1);
    return b;
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
    internal_ECMEsteps(x, asInteger(nsteps), asLogical(Verbp));
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
	D = GET_SLOT(x, Matrix_DSym),
	Omega = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	val;
    int *dRZX = INTEGER(getAttrib(RZXP, R_DimSymbol)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	Q, Qsqr, RZXpos, facepos,
	i, ione = 1, j, nf = length(Omega), p = dRZX[1] - 1, pos;
    SEXP gradComp = mer_gradComp(x);
    double
	*RZX = REAL(RZXP),
	*b = REAL(RZXP) + dRZX[0] * p,
	*bbface,		/* vec of second faces of gradComp elts */
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
	int nci = nc[i];
	int ncisqr = nci * nci;
	double *fDi = REAL(VECTOR_ELT(gradComp, i)),
	    mult = 1./((double)(Gp[i + 1] - Gp[i])/nci);

	Memcpy(bbface + pos, fDi + ncisqr, ncisqr);
	/* outer product of the third face of gradComp on the diagonal
	 * of the third face of val */
	F77_CALL(dsyr)("U", &ncisqr, &mult, fDi + 2 * ncisqr, &ione,
		       REAL(val) + 2 * Qsqr + pos * Q, &Q);
	pos += ncisqr;
    }
				/* fifth face is outer product of bbface */
    F77_CALL(dsyr)("U", &Q, &one, bbface, &ione, REAL(val) + 4 * Qsqr, &Q);
				/* fourth face from \bb\trans\der\vb\der\bb */
    AZERO(REAL(val) + 3 * Qsqr, Qsqr); /* zero accumulator */
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
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 * 
 * @param mer pointer to an mer object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 * 
 * @return a matrix
 */
SEXP
mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp)
{
    SEXP ans, Omega = GET_SLOT(x, Matrix_OmegaSym),
	Omegacp = PROTECT(duplicate(Omega));
    cholmod_factor *L =
      (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))), "REML"),
	i, ione = 1, j, n = LENGTH(GET_SLOT(x, Matrix_ySym)),
	nf = LENGTH(Omega), nsamp = asInteger(nsampp),
	p = LENGTH(GET_SLOT(x, Matrix_rXySym)),
	q = LENGTH(GET_SLOT(x, Matrix_rZySym)),
	saveb = asLogical(savebp),
	trans = asLogical(transp);
    double 
	*RXX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	*bhat = REAL(GET_SLOT(x, Matrix_ranefSym)),
	*betahat = REAL(GET_SLOT(x, Matrix_fixefSym)), 
	*bnew = Calloc(q, double), *betanew = Calloc(p, double),
	*dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	df = n - (REML ? p : 0), m1[] = {-1,0}, one[] = {1,0};
    int nrbase = p + 1 + coef_length(nf, nc); /* rows always included */
    int nrtot = nrbase + (saveb ? q : 0);
    cholmod_dense *chb, *chbnew = numeric_as_chm_dense(bnew, q);
    
    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    GetRNGstate();
    for (i = 0; i < nsamp; i++) {
	double *col = REAL(ans) + i * nrtot, sigma;
				/* factor x and get secondary values */
	mer_factor(x);
	mer_secondary(x);
				/* simulate and store new value of sigma */
	sigma = exp(dcmp[3]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
				/* simulate scaled, independent normals */
	for (j = 0; j < p; j++) betanew[j] = sigma * norm_rand();
	for (j = 0; j < q; j++) bnew[j] = sigma * norm_rand();
				/* betanew := RXX^{-1} %*% betanew */
	F77_CALL(dtrsv)("U", "N", "N", &p, RXX, &p, betanew, &ione);
				/* bnew := bnew - RZX %*% betanew */
	F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, betanew, &ione,
			one, bnew, &ione);
				/* chb := L^{-T} %*% bnew */
	chb = cholmod_solve(CHOLMOD_Lt, L, chbnew, &c);
				/* Copy chb to bnew and free chb */
	Memcpy(bnew, (double *)(chb->x), q);
 	cholmod_free_dense(&chb, &c);
				/* Add conditional modes and store beta */
	for (j = 0; j < p; j++) {
	    col[j] = (betanew[j] += betahat[j]);
	}
				/* Add conditional modes and
				 * optionally store b */
	for (j = 0; j < q; j++) {
	    bnew[j] += bhat[j];
	    if (saveb) col[nrbase + j] = bnew[j];
	}
	/* Update and store variance-covariance of the random effects */
	internal_Omega_update(Omega, bnew, sigma, nf, nc, Gp, col + p + 1, trans);
    }
    PutRNGstate();
    Free(betanew); Free(bnew); free(chbnew);
				/* Restore original Omega */
    SET_SLOT(x, Matrix_OmegaSym, Omegacp);
    mer_factor(x);

    UNPROTECT(2);
    return ans;
}

/** 
 * Create a Markov Chain Monte Carlo sample from a fitted generalized
 * linear mixed model
 * 
 * @param GSpt External pointer to a GlmerStruct
 * @param b Conditional modes of the random effects at the parameter
 * estimates
 * @param fixed Estimates of the fixed effects
 * @param varc Estimates of the variance components
 * @param savebp Logical indicator of whether or not to save the
 * random effects in the MCMC sample
 * @param nsampp Integer value of the number of samples to generate
 * 
 * @return 
 */
SEXP 
glmer_MCMCsamp(SEXP GSpt, SEXP b, SEXP fixed, SEXP varc,
	       SEXP savebp, SEXP nsampp) 
{ 
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSpt);
    int i, j, nf = LENGTH(b), nsamp = asInteger(nsampp),
	p = LENGTH(fixed), q = LENGTH(varc),
	saveb = asLogical(savebp);
    int *mcpar, nc = p + q;
    SEXP ans, mcparSym = install("mcpar");
    
    if (nsamp <= 0) nsamp = 1;
    nc = p + q;
    if (saveb)
	for (i = 0; i < nf; i++) {
	    int *dims = INTEGER(getAttrib(VECTOR_ELT(b, i),
					  R_DimSymbol));
	    nc += dims[0] * dims[1];
	}
    ans = PROTECT(allocMatrix(REALSXP, nsamp, nc));
    GetRNGstate();
    for (i = 0; i < nsamp; i++) {
	internal_glmer_fixef_update(GS, b, REAL(fixed));
	internal_bhat(GS, REAL(fixed), REAL(varc));
	internal_glmer_ranef_update(GS, b);
/* 	internal_Omega_update(GS->mer, b); */
	internal_mer_coef(GS->mer, 2, REAL(varc));
	for (j = 0; j < p; j++)
	    REAL(ans)[i + j * nsamp] = REAL(fixed)[j];
	for (j = 0; j < q; j++)
	    REAL(ans)[i + (j + p) * nsamp] = REAL(varc)[j];
	if (saveb) {
	    int base = p + q, k;
	    for (k = 0; k < nf; k++) {
		SEXP bk = VECTOR_ELT(b, k);
		int *dims = INTEGER(getAttrib(bk, R_DimSymbol));
		int klen = dims[0] * dims[1];

		for (j = 0; j < klen; j++)
		    REAL(ans)[i + (j + base) * nsamp] = REAL(bk)[j];
		base += klen;
	    }
	}
    }
    PutRNGstate();
    UNPROTECT(1);
				/* set (S3) class and mcpar attribute */
    setAttrib(ans, R_ClassSymbol, mkString("mcmc"));
    setAttrib(ans, mcparSym, allocVector(INTSXP, 3));
    mcpar = INTEGER(getAttrib(ans, mcparSym));
    mcpar[0] = mcpar[2] = 1;
    mcpar[1] = nsamp;

    return ans;
} 

/**
 * Create an mer object from a list of grouping factors and a list of model
 * matrices. 
 *
 * @param fl named list of grouping factors
 * @param Ztl list of transposes of model matrices
 * @param Xp model matrix for the fixed effects
 * @param yp response vector
 * @param method character vector describing the estimation method
 *
 * @return pointer to an mer object
 */
SEXP mer_create(SEXP fl, SEXP ZZt, SEXP Xp, SEXP yp, SEXP method,
		 SEXP ncp, SEXP cnames, SEXP useS, SEXP call,
		 SEXP family)
{
    SEXP Omega, bVar, gradComp, fnms = getAttrib(fl, R_NamesSymbol),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer"))), xnms;
    cholmod_sparse *ts1, *ts2, *Zt;
    cholmod_factor *F;
    int *nc = INTEGER(ncp), *Gp, *xdims, i,
	nf = LENGTH(fl), nobs = LENGTH(yp), p, q;
    char *devCmpnms[] = {"n", "p", "yty", "logryy2", "logDetL2",
			 "logDetOmega", "logDetRXX", ""};
    char *devnms[] = {"ML", "REML", ""};
    double *dcmp;
				/* Check arguments to be duplicated */
    if (!isReal(yp)) error(_("yp must be a real vector"));
    SET_SLOT(val, Matrix_ySym, duplicate(yp));
    if (!isMatrix(Xp) || !isReal(Xp))
	error(_("Xp must be a real matrix"));
    xdims = INTEGER(getAttrib(Xp, R_DimSymbol));
    if (xdims[0] != nobs) error(_("Xp must have %d rows"), nobs);
    p = xdims[1];
    xnms = VECTOR_ELT(getAttrib(Xp, R_DimNamesSymbol), 1);
    SET_SLOT(val, Matrix_XSym, duplicate(Xp));
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);
    }
    SET_SLOT(val, Matrix_flistSym, duplicate(fl));
    if (!isString(method) || LENGTH(method) != 1)
	error(_("method must be a character vector of length 1"));
    SET_SLOT(val, Matrix_methodSym, duplicate(method));
    if (!isLogical(useS) || LENGTH(useS) != 1)
	error(_("useS must be a logical vector of length 1"));
    SET_SLOT(val, Matrix_useScaleSym, duplicate(useS));
    if (!isNewList(cnames) || LENGTH(cnames) != nf + 1)
	error(_("cnames must be a list of length %d"), nf + 1);
    SET_SLOT(val, Matrix_cnamesSym, duplicate(cnames));
    if (!isInteger(ncp) || LENGTH(ncp) != nf)
	error(_("ncp must be an integer vector of length %d"), nf);
    SET_SLOT(val, Matrix_callSym, duplicate(call));
    SET_SLOT(val, Matrix_familySym, duplicate(family));
    SET_SLOT(val, Matrix_ncSym, duplicate(ncp));
    Gp = INTEGER(ALLOC_SLOT(val, Matrix_GpSym, INTSXP, nf + 1));
    Gp[0] = 0;
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);

    }
    SET_SLOT(val, Matrix_ZtSym, duplicate(ZZt));
    Zt = as_cholmod_sparse(GET_SLOT(val, Matrix_ZtSym));
				/* analyze Zt */
    q = Zt->nrow;
    i = c.supernodal;
    c.supernodal = CHOLMOD_SUPERNODAL;
    F = cholmod_analyze(Zt, &c);
    c.supernodal = i;
/* FIXME: Need to set up a finalizer for F. */
/* Right now that storage is not being released. */
    SET_SLOT(val, Matrix_LSym, R_MakeExternalPtr(F, R_NilValue, val));
    ts1 = cholmod_aat(Zt, (int*)NULL/* fset */,(size_t)0,
		      CHOLMOD_PATTERN, &c);
    ts2 = cholmod_copy(ts1, 1/*upper triangle*/, CHOLMOD_PATTERN, &c);
    SET_SLOT(val, Matrix_ZtZSym,
	     alloc_dsCMatrix(q, cholmod_nnz(ts2, &c), "U", R_NilValue,
			     R_NilValue));
    cholmod_free_sparse(&ts1, &c); cholmod_free_sparse(&ts2, &c);
				/* allocate other slots */
    SET_SLOT(val, Matrix_devianceSym, Matrix_make_named(REALSXP, devnms));
    SET_SLOT(val, Matrix_devCompSym, Matrix_make_named(REALSXP, devCmpnms));
    dcmp = REAL(GET_SLOT(val, Matrix_devCompSym));
    AZERO(dcmp, 7);		/* cosmetic */
    dcmp[0] = (double) nobs;
    dcmp[1] = (double) p;
				/* allocate and populate list slots */
    Omega = ALLOC_SLOT(val, Matrix_OmegaSym, VECSXP, nf);
    bVar = ALLOC_SLOT(val, Matrix_bVarSym, VECSXP, nf);
    gradComp = ALLOC_SLOT(val, Matrix_gradCompSym, VECSXP, nf);
    setAttrib(Omega, R_NamesSymbol, duplicate(fnms));
    setAttrib(bVar, R_NamesSymbol, duplicate(fnms));
    setAttrib(gradComp, R_NamesSymbol, duplicate(fnms));
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
				/* create ZtX, RZX, XtX, RXX */
    SET_SLOT(val, Matrix_ZtXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, Matrix_RZXSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
    SET_SLOT(val, Matrix_XtXSym, alloc_dpoMatrix(p, "U", xnms, xnms));
    SET_SLOT(val, Matrix_RXXSym, alloc_dtrMatrix(p, "U", "N", xnms, xnms));
    SET_SLOT(val, Matrix_ZtySym, allocVector(REALSXP, q));
    SET_SLOT(val, Matrix_rZySym, allocVector(REALSXP, q));
    SET_SLOT(val, Matrix_XtySym, allocVector(REALSXP, p));
    SET_SLOT(val, Matrix_rXySym, allocVector(REALSXP, p));
    mer_update_ZXy(val);
    Free(Zt);
				/* secondary slots */
    SET_SLOT(val, Matrix_ranefSym, allocVector(REALSXP, q));
    SET_SLOT(val, Matrix_RZXinvSym, alloc_dgeMatrix(q, p, R_NilValue, xnms));
				/* initialize */
    mer_initial(val);
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
    int	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym));
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
    int clen = coef_length(LENGTH(GET_SLOT(x, Matrix_flistSym)),
			   INTEGER(GET_SLOT(x, Matrix_ncSym)));   
    if (LENGTH(coef) != clen || !isReal(coef))
	error(_("coef must be a numeric vector of length %d"), clen);
    internal_mer_coefGets(x, REAL(coef), asInteger(pType));
    return x;
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
    cholmod_factor *L;
    cholmod_sparse *Lm;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	nz, q;

    L = cholmod_copy_factor((cholmod_factor *)
			    R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym)),
			    &c);
    dims[0] = dims[1] = q = (int)(L->n);
    Lm = cholmod_factor_to_sparse(L, &c); cholmod_free_factor(&L, &c);
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, q + 1)),
	   (int *) Lm->p, q + 1);
    nz = ((int *)(Lm->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)),
	   (int *) Lm->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)),
	   (double *) Lm->x, nz);
    cholmod_free_sparse(&Lm, &c);
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
    cholmod_factor
	*L = (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    cholmod_sparse
	*b = cholmod_allocate_sparse(L->n, L->n, L->n, 1, 1,
				     0, CHOLMOD_REAL, &c),
	*Linv;
    double *bx = (double *)(b->x);
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtCMatrix")));
    int *bi = (int *) (b->i), *bp = (int *) (b->p),
	*dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	j, nz, q;

    dims[0] = dims[1] = q = (int)(L->n);
    for (j = 0; j < q; j++) {
	bp[j] = bi[j] = j;
	bx[j] = 1;
    }
    bp[q] = q;
    Linv = cholmod_spsolve(CHOLMOD_L, L, b, &c);
    cholmod_free_sparse(&b, &c);
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    SET_SLOT(ans, Matrix_diagSym, mkString("N"));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, q + 1)),
	   (int *) Linv->p, q + 1);
    nz = ((int *)(Linv->p))[q];
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nz)),
	   (int *) Linv->i, nz);
    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nz)),
	   (double *) Linv->x, nz);
    cholmod_free_sparse(&Linv, &c);
    UNPROTECT(1);
    return ans;
}

/**
 * Create and factor Z'Z+Omega.  Also create RZX and RXX, the deviance
 * components, and the value of the deviance for both ML and REML.
 *
 * @param x pointer to an lmer object
 *
 * @return NULL
 */
SEXP mer_factor(SEXP x)
{
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    cholmod_sparse *A,
	*zz = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtZSym));
    cholmod_factor *L =
	(cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    cholmod_dense *ZtX = as_cholmod_dense(GET_SLOT(x, Matrix_ZtXSym)),
	*Zty = numeric_as_chm_dense(REAL(GET_SLOT(x, Matrix_ZtySym)), L->n),
	*rZy, *RZX;
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)), i, info, ione = 1,
	nf = LENGTH(Omega), p = ZtX->ncol, q = L->n;
    double *RXX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
	*rXy = REAL(GET_SLOT(x, Matrix_rXySym)),
	*dcmp = REAL(GET_SLOT(x, Matrix_devCompSym)),
	*dev = REAL(GET_SLOT(x, Matrix_devianceSym)),
	one[2] = {1, 0}, m1[2] = {-1, 0};
    double nml = dcmp[0], nreml = dcmp[0] - dcmp[1];
	    
    dcmp[5] = Omega_log_det(Omega, nf, nc, Gp); /* logDet(Omega) */
    A = ZZ_inflate(zz, nf, Omega, nc, Gp); free(zz);
    if (!cholmod_factorize(A, L, &c))
	error(_("rank_deficient Z'Z+Omega"));
    cholmod_free_sparse(&A, &c);
    dcmp[4] = 2 * chm_log_abs_det(L); /* 2 * logDet(L) */
				/* calculate and store RZX and rZy */
    RZX = cholmod_solve(CHOLMOD_L, L, ZtX, &c); free(ZtX);
    rZy = cholmod_solve(CHOLMOD_L, L, Zty, &c); free(Zty);
    Memcpy(REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXSym), Matrix_xSym)),
	   (double *) RZX->x, q * p);
    Memcpy(REAL(GET_SLOT(x, Matrix_rZySym)), (double *) rZy->x, q);
				/* downdate XtX and factor */
    Memcpy(RXX, REAL(GET_SLOT(GET_SLOT(x, Matrix_XtXSym), Matrix_xSym)), p * p);
    F77_CALL(dsyrk)("U", "T", &p, &q, m1, (double*)RZX->x, &q, one, RXX, &p);
    F77_CALL(dpotrf)("U", &p, RXX, &p, &info);
    if (info) {
	error(_("Leading minor of order %d in downdated X'X is not positive definite"),
	      info);
	dcmp[3] = dcmp[6] = dev[0] = dev[1] = NA_REAL;
    } else {
	for (dcmp[6] = 0, i = 0; i < p; i++) /* 2 * logDet(RXX) */
	    dcmp[6] += 2. * log(RXX[i * (p + 1)]);
				/* solve for rXy  and ryy^2 */
	Memcpy(rXy, REAL(GET_SLOT(x, Matrix_XtySym)), p);
	F77_CALL(dgemv)("T", &q, &p, m1, (double*) RZX->x, &q,
			(double*) rZy->x, &ione, one, rXy, &ione);
	F77_CALL(dtrsv)("U", "T", "N", &p, RXX, &p, rXy, &ione);
	dcmp[3] = log(dcmp[2] /* dcmp[3] = log(ryy^2); dcmp[2] = y'y; */
		      - F77_CALL(ddot)(&p, rXy, &ione, rXy, &ione)
		      - F77_CALL(ddot)(&q, (double*)rZy->x, &ione,
				       (double*)rZy->x, &ione));
				/* evaluate ML and REML deviance */
	dev[0] = dcmp[4] - dcmp[5] +
	    nml*(1.+dcmp[3]+log(2.*PI/nml));
	dev[1] = dcmp[4] - dcmp[5] + dcmp[6] +
	    nreml*(1.+dcmp[3]+log(2.*PI/nreml));
    }
	    
    cholmod_free_dense(&RZX, &c); cholmod_free_dense(&rZy, &c);
				/* signal that secondary slots are not valid */
    SET_SLOT(x, Matrix_fixefSym, allocVector(REALSXP, 0));
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

SEXP mer_fitted(SEXP x, SEXP useFe, SEXP useRe)
{
    int n = LENGTH(GET_SLOT(x, Matrix_ySym));
    SEXP ans = PROTECT(allocVector(REALSXP, n));

    AZERO(REAL(ans), n);
    internal_mer_fitted(x,
			(asLogical(useFe)
			 ? REAL(GET_SLOT(x, Matrix_XSym))
			 : (double *) NULL),
			(asLogical(useRe)
			 ? REAL(GET_SLOT(GET_SLOT(x, Matrix_ZtSym), Matrix_xSym))
			 : (double *) NULL),
			REAL(ans));
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
    int nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym));
    SEXP ans;
    
    mer_secondary(x);
    ans = PROTECT(duplicate(GET_SLOT(x, Matrix_fixefSym)));
    setAttrib(ans, R_NamesSymbol,
	      duplicate(VECTOR_ELT(GET_SLOT(x, Matrix_cnamesSym), nf)));
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
    SEXP bVarP = GET_SLOT(x, Matrix_bVarSym),
	OmegaP = GET_SLOT(x, Matrix_OmegaSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	gradComp = GET_SLOT(x, Matrix_gradCompSym),
	ranefP = GET_SLOT(x, Matrix_ranefSym);
    cholmod_factor
	*L = (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    cholmod_dense *RZX = as_cholmod_dense(RZXP), *tmp1;
    int q = LENGTH(ranefP), p = LENGTH(GET_SLOT(x, Matrix_rXySym));
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*Perm = (int *)(L->Perm),
	*iperm = Calloc(q, int),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, j, k, nf = length(OmegaP);
    double *b = REAL(GET_SLOT(x, Matrix_ranefSym)), m1[] = {-1, 0};

    for (j = 0; j < q; j++) iperm[Perm[j]] = j;
    mer_secondary(x);
    /* FIXME: Store this array (and perhaps calculate it as part of mer_secondary) */
    tmp1 = cholmod_solve(CHOLMOD_Lt, L, RZX, &c); free(RZX);
    F77_CALL(dtrsm)("R", "U", "N", "N", &q, &p, m1,
		    REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym),
				  Matrix_xSym)), &p,
		    (double *)(tmp1->x), &q);
    for (i = 0; i < nf; i++) {
	SEXP bVPi = VECTOR_ELT(bVarP, i);
	int nci = nc[i], RZXrows = Gp[i + 1] - Gp[i];
	cholmod_sparse
	    *sm = cholmod_allocate_sparse(q, nci, nci, TRUE, TRUE, 0,
					  CHOLMOD_REAL, &c),
	    *tsm1, *tsm2, *tsm3;
	int *si = (int*)(sm->i), *sp = (int*)(sm->p),
	    ncisq = nci * nci, nlev = RZXrows/nci;
	double *bVi = REAL(bVPi),
	    *bi = b + Gp[i], *mm = REAL(VECTOR_ELT(gradComp, i)),
	    *sx = (double*)(sm->x),
	    *tmp = Memcpy(Calloc(ncisq, double),
			  REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(OmegaP, i)),
			       Matrix_xSym)), ncisq),
	    *RZXi = (double *)(tmp1->x) + Gp[i],
	    dlev = (double) nlev,
	    one[] = {1,0}, zero[] = {0,0};

	for (j = 0; j < nci; j++) {
	    sp[j] = j;
	    sx[j] = 1;
	}
	sp[nci] = nci;
	
	for (k = 0; k < nlev; k++) {
	    int rk0 = Gp[i] + k * nci, nnz;
	    for (j = 0; j < nci; j++) si[j] = iperm[rk0 + j];
	    tsm1 = cholmod_spsolve(CHOLMOD_Lt, L, sm, &c);
	    nnz = cholmod_nnz(tsm1, &c);
	    if (nci == 1) {
		double *xp = (double*)(tsm1->x);
		int ione = 1;
		bVi[k] = F77_CALL(ddot)(&nnz, xp, &ione, xp, &ione);
	    } else {
		/* FIXME: Create a cholmod_ata that uses dsyrk and an array to hold the i indices. */
		tsm2 = cholmod_transpose(tsm1, 1, &c);
		tsm3 = cholmod_aat(tsm2, (int*) NULL, 0, 1, &c);
		/* Now copy the contents to bVi */
	    }
	}
 	if (nci == 1) {
	    int ione = 1;
 	    mm[0] = ((double) nlev)/(tmp[0] * tmp[0]);
 	    mm[1] = F77_CALL(ddot)(&nlev, bi, &ione, bi, &ione);
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
	    F77_CALL(dsyrk)("U", "N", &nci, &nlev, one, bi, &nci,
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
    SEXP Omega = GET_SLOT(x, Matrix_OmegaSym);
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	dind, i, ifour = 4, info, ione = 1, nf = length(Omega),
	odind, ptyp = asInteger(pType);
    SEXP
	gradComp = mer_gradComp(x),
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
			REAL(VECTOR_ELT(gradComp, i)), &ncisqr,
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
 * Create and insert initial values for Omega.
 *
 * @param x pointer to an mer object
 *
 * @return NULL
 */
SEXP mer_initial(SEXP x)
{
    SEXP Omg = GET_SLOT(x, Matrix_OmegaSym),
	ZtZ = GET_SLOT(x, Matrix_ZtZSym);
    int	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*p = INTEGER(GET_SLOT(ZtZ, Matrix_pSym)),
	i, nf = length(Omg);
    double *xx = REAL(GET_SLOT(ZtZ, Matrix_xSym));

    for (i = 0; i < nf; i++) {
	double *omgi = REAL(GET_SLOT(VECTOR_ELT(Omg, i),
				     Matrix_xSym));
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
	dpoMatrix_chol(VECTOR_ELT(Omg, i));
    }
    mer_factor(x);
    return R_NilValue;
}

/** 
 * Return the permutation of the columns of Z as a pMatrix object
 * 
 * @param x pointer to an mer object
 * 
 * @return the permutation as an pMatrix object
 */
SEXP mer_pMatrix(SEXP x)
{
    cholmod_factor *L =
      (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("pMatrix")));
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2)),
	*perm = INTEGER(ALLOC_SLOT(ans, Matrix_permSym, INTSXP, L->n)), i;

    dims[0] = dims[1] = (int) L->n;
    for (i = 0; i < (int) L->n; i++) perm[i] = ((int *)(L->Perm))[i] + 1;
    UNPROTECT(1);
    return ans;
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
    SEXP cnames = GET_SLOT(x, Matrix_cnamesSym),
	flist = GET_SLOT(x, Matrix_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	*nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	i, ii, jj,
	nf = length(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, Matrix_ranefSym));

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
 * Update the secondary slots (fixef, ranef, RZXinv, bVar, gradComp)
 * 
 * @param x pointer to a mer object
 * 
 */
SEXP mer_secondary(SEXP x)
{
    SEXP fixef = GET_SLOT(x, Matrix_fixefSym);

    if (!LENGTH(fixef)) {
	SEXP rXy = GET_SLOT(x, Matrix_rXySym),
	    RZXP = GET_SLOT(x, Matrix_RZXSym),
	    ranef = GET_SLOT(x, Matrix_ranefSym);
	int ione = 1, p = LENGTH(rXy), q = LENGTH(ranef);
	cholmod_factor *L =
	    (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
	cholmod_dense *td1, *td2,
	    *chRZX = as_cholmod_dense(RZXP),
	    *chranef = numeric_as_chm_dense(REAL(ranef), q);
	double *RXX = REAL(GET_SLOT(GET_SLOT(x, Matrix_RXXSym), Matrix_xSym)),
	    *RZX = REAL(GET_SLOT(RZXP, Matrix_xSym)),
	    *RZXinv = REAL(GET_SLOT(GET_SLOT(x, Matrix_RZXinvSym), Matrix_xSym)),
	    m1[] = {-1,0}, one[] = {1,0};
				/* allocate and copy */
	fixef = ALLOC_SLOT(x, Matrix_fixefSym, REALSXP, p);
	Memcpy(REAL(fixef), REAL(rXy), p);
	Memcpy(REAL(ranef), REAL(GET_SLOT(x, Matrix_rZySym)), q);
				/* fixef, ranef, RZXinv */
	F77_CALL(dtrsv)("U", "N", "N", &p, RXX,	&p, REAL(fixef), &ione);
	F77_CALL(dgemv)("N", &q, &p, m1, RZX, &q, REAL(fixef), &ione,
			one, REAL(ranef), &ione);
	td1 = cholmod_solve(CHOLMOD_Lt, L, chranef, &c);
	td2 = cholmod_solve(CHOLMOD_Pt, L, td1, &c);
	Memcpy(REAL(ranef), (double *)(td2->x), q);
	free(chranef); cholmod_free_dense(&td1, &c); cholmod_free_dense(&td2, &c);
	td1 = cholmod_solve(CHOLMOD_Lt, L, chRZX, &c); free(chRZX);
	F77_CALL(dtrsm)("R", "U", "N", "N", &q, &p, m1, RXX, &p,
			(double *)(td1->x), &q);
	Memcpy(RZXinv, (double *)(td1->x), p * q);
	cholmod_free_dense(&td1, &c);
	/* FIXME: Should bVar and gradComp be updated unconditionally? */
    }
    return R_NilValue;
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
 * @param useScale Logical indicator of whether to use the scale
 * 
 * @return a matrix of simulated linear predictors
 */
SEXP mer_simulate(SEXP x, SEXP nsimP)
{
    int *nc = INTEGER(GET_SLOT(x, Matrix_ncSym)),
	*Gp = INTEGER(GET_SLOT(x, Matrix_GpSym)),
	REML = !strcmp(CHAR(asChar(GET_SLOT(x, Matrix_methodSym))),"REML"),
	i, ii, j, nsim = asInteger(nsimP),
	nf = LENGTH(GET_SLOT(x, Matrix_OmegaSym)),
	n = LENGTH(GET_SLOT(x, Matrix_ySym)),
	q = LENGTH(GET_SLOT(x, Matrix_ZtySym));
    SEXP ans = PROTECT(allocMatrix(REALSXP, n, nsim)),
	Omega = GET_SLOT(x, Matrix_OmegaSym);
    cholmod_dense *cha = as_cholmod_dense(ans),
	*chb = cholmod_allocate_dense(q, nsim, q, CHOLMOD_REAL, &c);
    double one[] = {1,0}, zero[] = {0,0},
	scale = (asLogical(GET_SLOT(x, Matrix_useScaleSym)) ?
		 internal_mer_sigma(x, REML) : 1);
    cholmod_sparse *Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
	
    GetRNGstate();
    for (ii = 0; ii < nsim; ii++) {
	for (i = 0; i < nf; i++) {
	    int nci = nc[i], relen = Gp[i + 1] - Gp[i];
	    int nlev = relen/nci;
	    double *bi = (double *)(chb->x) + ii * q + Gp[i],
		*Rmat = REAL(GET_SLOT(dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				      Matrix_xSym));

	    for (j = 0; j < relen; j++) bi[j] = norm_rand();
	    F77_CALL(dtrsm)("L", "U", "N", "N", &nci, &nlev, &scale,
			    Rmat, &nci, bi, &nci);
	}
    }
    PutRNGstate();

    if (!cholmod_sdmult(Zt, 1, one, zero, chb, cha, &c))
	error(_("cholmod_sdmult failed"));
    cholmod_free_dense(&chb, &c);
    free(Zt); free(cha);
    UNPROTECT(1);
    return ans;
}

/** 
 * Update the derived quantities (RZX, rZy, RXX, rxy, and dcmp[2] = y'y
 * when Z, X or y has been changed.
 * 
 * @param x pointer to an mer object
 * 
 * @return NULL
 */
SEXP mer_update_ZXy(SEXP x)
{
    SEXP Xp = GET_SLOT(x, Matrix_XSym), ZtZ = GET_SLOT(x, Matrix_ZtZSym);
    SEXP ZtZx = GET_SLOT(ZtZ, Matrix_xSym),
	ZtZp = GET_SLOT(ZtZ, Matrix_pSym), ZtZi = GET_SLOT(ZtZ, Matrix_iSym);
    int *dims = INTEGER(getAttrib(Xp, R_DimSymbol)), i, ione = 1, j;
    cholmod_factor
	*L = (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    int *perm = (int *)L->Perm,
	n = dims[0], nnz, p = dims[1], q = (int)L->n;
    cholmod_sparse *ts1, *ts2,
	*Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
    double *X = REAL(Xp),
	*XtX = REAL(GET_SLOT(GET_SLOT(x, Matrix_XtXSym), Matrix_xSym)),
	*Xty = REAL(GET_SLOT(x, Matrix_XtySym)),
	*ZtX = REAL(GET_SLOT(GET_SLOT(x, Matrix_ZtXSym), Matrix_xSym)),
	*Zty = REAL(GET_SLOT(x, Matrix_ZtySym)),
	*y = REAL(GET_SLOT(x, Matrix_ySym)),
	one[] = {1, 0}, zero[] = {0,0};
    cholmod_dense *td1, *Xd = as_cholmod_dense(Xp),
	*yd = numeric_as_chm_dense(y, n);
				/* y'y */
    REAL(GET_SLOT(x, Matrix_devCompSym))[2] =
	F77_CALL(ddot)(&n, y, &ione, y, &ione); 
				/* ZtZ */
    ts1 = cholmod_aat(Zt, (int *) NULL, (size_t) 0, 1/* mode */, &c);
    /* cholmod_aat returns stype == 0; copy to set stype == 1 */ 
    ts2 = cholmod_copy(ts1, 1/* stype */, 1/* mode */, &c);
    nnz = cholmod_nnz(ts2, &c);
    if (((int)(ts2->ncol) + 1) != LENGTH(ZtZp))
	error(_("Order Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZp) - 1, (int)(ts2->ncol));
    Memcpy(INTEGER(ZtZp), (int*)(ts2->p), LENGTH(ZtZp));
    if (nnz != LENGTH(ZtZx))
	error(_("Number of nonzeros in Z'Z has changed - was %d, now %d"),
	      LENGTH(ZtZx), nnz);
    Memcpy(INTEGER(ZtZi), (int*)(ts2->i), nnz);
    Memcpy(REAL(ZtZx), (double *)(ts2->x), nnz);
    cholmod_free_sparse(&ts1, &c); cholmod_free_sparse(&ts2, &c);
				/* PZ'X into ZtX */
    td1 = cholmod_allocate_dense(q, p, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Zt, 0, one, zero, Xd, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (j = 0; j < p; j++) { 	/* permute the columns */
	double *dcol = ZtX + j * q,
	    *scol = (double *)(td1->x) + j * q;
	for (i = 0; i < q; i++) dcol[i] = scol[perm[i]];
    }
    cholmod_free_dense(&td1, &c); Free(Xd);
				/* PZ'y into Zty */
    td1 = cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Zt, 0, one, zero, yd, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (i = 0; i < q; i++) Zty[i] = ((double *)(td1->x))[perm[i]];
    cholmod_free_dense(&td1, &c); Free(yd); Free(Zt);
				/* XtX and Xty */
    AZERO(XtX, p * p);		
    F77_CALL(dsyrk)("U", "T", &p, &n, one, X, &n, zero, XtX, &p);
    F77_CALL(dgemv)("T", &n, &p, one, X, &n, y, &ione, zero, Xty, &ione);
    SET_SLOT(x, Matrix_fixefSym, allocVector(REALSXP, 0));
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
    SEXP y = GET_SLOT(x, Matrix_ySym),
	Xty = GET_SLOT(x, Matrix_XtySym),
	Zty = GET_SLOT(x, Matrix_ZtySym);
    cholmod_factor *L =
      (cholmod_factor *) R_ExternalPtrAddr(GET_SLOT(x, Matrix_LSym));
    int *perm = (int*)(L->Perm), i, ione = 1,
	n = LENGTH(y), p = LENGTH(Xty), q = LENGTH(Zty);
    cholmod_sparse *Zt = as_cholmod_sparse(GET_SLOT(x, Matrix_ZtSym));
    cholmod_dense *td1, *yd = numeric_as_chm_dense(REAL(y), n);
    double one[] = {1,0}, zero[] = {0,0};

    if (!isReal(ynew) || LENGTH(ynew) != n)
	error(_("ynew must be a numeric vector of length %d"), n);
    Memcpy(REAL(y), REAL(ynew), n);
    				/* y'y */
    REAL(GET_SLOT(x, Matrix_devCompSym))[2] =
	F77_CALL(ddot)(&n, REAL(y), &ione, REAL(y), &ione); 
				/* PZ'y into Zty */
    td1 = cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    if (!cholmod_sdmult(Zt, 0, one, zero, yd, td1, &c))
	error(_("cholmod_sdmult failed"));
    for (i = 0; i < q; i++) REAL(Zty)[i] = ((double *)(td1->x))[perm[i]];
    cholmod_free_dense(&td1, &c); free(yd); free(Zt);
    				/* Xty */
    F77_CALL(dgemv)("T", &n, &p, one, REAL(GET_SLOT(x, Matrix_XSym)),
		    &n, REAL(y), &ione, zero, REAL(Xty), &ione);
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
	/* ZZxP = GET_SLOT(x, Matrix_ZZxSym), */
	ZtXP = GET_SLOT(x, Matrix_ZtXSym),
	XtXP = GET_SLOT(x, Matrix_XtXSym),
	RZXP = GET_SLOT(x, Matrix_RZXSym),
	RXXP = GET_SLOT(x, Matrix_RXXSym)
	/* , cnames = GET_SLOT(x, Matrix_cnamesSym) */
	;
    int *ZtXd = INTEGER(getAttrib(ZtXP, Matrix_DimSym)),
	*XtXd = INTEGER(getAttrib(XtXP, Matrix_DimSym));

    if (!match_mat_dims(ZtXd, INTEGER(getAttrib(RZXP, Matrix_DimSym))))
	return mkString(_("Dimensions of slots ZtX and RZX must match"));
    if (!match_mat_dims(XtXd, INTEGER(getAttrib(RXXP, Matrix_DimSym))))
	return mkString(_("Dimensions of slots XtX and RXX must match"));
    if (ZtXd[1] != XtXd[0] || XtXd[0] != XtXd[1])
	return mkString(_("Slot XtX must be a square matrix with same ncol as ZtX"));
    return ScalarLogical(1);
}

/** 
 * Create the sparse Zt matrix from a factor list and list of model matrices
 * 
 * @param fl list of factors
 * @param Ztl list of transposes of model matrices
 * 
 * @return a freshly created sparse Zt object
 */
SEXP Zt_create(SEXP fl, SEXP Ztl)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), fi, tmmat;
    int *dims, *p, *ii, i, nrtot = 0, nf = LENGTH(fl), nobs;
    int *Gp = Calloc(nf + 1, int), *nr = Calloc(nf, int),
	*offset = Calloc(nf + 1, int);
    double *x;
    
    if (!isNewList(fl) || nf < 1)
	error(_("fl must be a non-null list"));
    if (!isNewList(Ztl) || LENGTH(Ztl) != nf)
	error(_("Ztl must be a list of length %d"), nf);
    fi = VECTOR_ELT(fl, 0);
    nobs = LENGTH(fi);
    if (!isFactor(fi) || nobs < 1)
	error(_("fl[[1]] must be a non-empty factor"));
    offset[0] = Gp[0] = 0;
    for (i = 0; i < nf; i++) {	/* check consistency and establish dimensions */
	fi = VECTOR_ELT(fl, i);	/* grouping factor */
	if (!isFactor(fi) || LENGTH(fi) != nobs)
	    error(_("fl[[%d]] must be a factor of length %d"), i + 1, nobs);
	tmmat = VECTOR_ELT(Ztl, i); /* transpose of model matrix */
	if (!isMatrix(tmmat) || !isReal(tmmat))
	    error(_("Ztl[[%d]] must be real matrix"), i + 1);
	dims = INTEGER(getAttrib(tmmat, R_DimSymbol));
	if (dims[1] != nobs)
	    error(_("Ztl[[%d]] must have %d columns"), i + 1, nobs);
	nrtot += (nr[i] = dims[0]);
	offset[i + 1] = offset[i] + nr[i];
	Gp[i + 1] = Gp[i] + dims[0] * LENGTH(getAttrib(fi, R_LevelsSymbol));
    }
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = Gp[nf]; dims[1] = nobs;
    p = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, nobs + 1));
    ii = INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nrtot * nobs));
    x = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nrtot * nobs));
    p[0] = 0; for(i = 0; i < nobs; i++) p[i + 1] = p[i] + nrtot;

    for (i = 0; i < nf; i++) {	/* fill ans */
	int *vals = INTEGER(VECTOR_ELT(fl, i)), j;
	double *xvals = REAL(VECTOR_ELT(Ztl, i));

	for (j = 0; j < nobs; j++) {
	    int jbase = Gp[i] + nr[i] * (vals[j] - 1), k;
	    for (k = 0; k < nr[i]; k++) {
		int ind = j * nrtot + offset[i] + k;
		ii[ind] = jbase + k;
		x[ind] = xvals[j * nr[i] + k];
	    }
	}
    }

    Free(offset); Free(Gp); Free(nr);
    UNPROTECT(1);
    return ans;
}
