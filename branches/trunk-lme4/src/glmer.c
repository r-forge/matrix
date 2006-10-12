#include "glmer.h"

/**
 * Evaluate the quadratic form in b defined by Omega
 *
 * @param b vector of random effects
 * @param Omega - list of dpoMatrix objects defining the pattern for Omega
 * @param nf - number of grouping factors
 * @param Gp - group pointers
 * @param nc - number of columns per factor
 *
 * @return value of the quadratic form in b
 */
static double
b_quadratic(const double b[], SEXP Omega, const int Gp[], const int nc[])
{
    int i, ione = 1, nf = LENGTH(Omega);
    double ans = 0., one[] = {1.,0.};

    for (i = 0; i < nf; i++) {
	int nci = nc[i], ntot = Gp[i + 1] - Gp[i];
	int nlev = ntot/nci;
	double *bcp = Memcpy(Calloc(ntot, double), b + Gp[i], ntot),
	    *omgf = REAL(GET_SLOT(M_dpoMatrix_chol(VECTOR_ELT(Omega, i)),
				  lme4_xSym));

	F77_CALL(dtrmm)("L", "U", "N", "N", &nci, &nlev, one, omgf,
			&nci, bcp, &nci);
	ans += F77_CALL(ddot)(&ntot, bcp, &ione, bcp, &ione);
	Free(bcp);
    }
    return ans;
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
 * Find a variable of a given name in a given environment and check
 * that its length and mode are correct.
 *
 * @param rho Environment in which to find the variable
 * @param nm Name of the variable to find
 * @param mode Desired mode
 * @param len Desired length (ignored if <= 0)
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
    if (len > 0 && LENGTH(ans) != len)
	error(_("object `%s' must be of length %d"),
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

static double
internal_Gaussian_deviance(int p, int q, cholmod_factor *L,
			   double RZX[], double RXX[], double betahat[],
			   double bhat[], double beta[], double b[])
{
    int i, ione = 1;
    double *bb = Calloc(q, double), *betab = Calloc(p, double), ans = 0;
    cholmod_sparse *Lm;
    cholmod_dense *chb = M_numeric_as_chm_dense(bb, q),
	*Ltb = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    cholmod_factor *Lcp;
    double one[] = {1,0}, zero[] = {0,0};

    for (i = 0; i < p; i++) betab[i] = beta[i] - betahat[i];
    for (i = 0; i < q; i++) bb[i] = b[i] - bhat[i];
    Lcp = M_cholmod_copy_factor(L, &c); /* next call changes Lcp */
    Lm = M_cholmod_factor_to_sparse(Lcp, &c); M_cholmod_free_factor(&Lcp, &c);
    if (!M_cholmod_sdmult(Lm, 1 /* transpose */, one, zero, chb, Ltb, &c))
	error(_("Error return from cholmod_sdmult"));
    Memcpy(bb, (double *)(Ltb->x), q);
    M_cholmod_free_sparse(&Lm, &c); M_cholmod_free_dense(&Ltb, &c);
    F77_CALL(dgemv)("N", &q, &p, one, RZX, &q, betab, &ione, one, bb, &ione);
    for (i = 0; i < q; i++) ans += bb[i] * bb[i];

    F77_CALL(dtrmv)("U", "N", "N", &p, RXX, &p, betab, &ione);
    for (i = 0; i < p; i++) ans += betab[i] * betab[i];
    Free(chb); Free(bb); Free(betab);
    return ans;
}

/**
 * Evaluate new weights and working residuals.
 *
 * @param GS a GlmerStruct object
 */
static void
internal_glmer_reweight(GlmerStruct GS) {
    SEXP dmu_deta, var;
    int i;
    double *w = REAL(GET_SLOT(GS->mer, lme4_wtsSym)),
	*y = REAL(GET_SLOT(GS->mer, lme4_ySym)),
	*z = REAL(GET_SLOT(GS->mer, lme4_wrkresSym));

				/* reweight mer */
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    dmu_deta = PROTECT(eval_check(GS->mu_eta, GS->rho,
				  REALSXP, GS->n));
    var = PROTECT(eval_check(GS->var, GS->rho,
			     REALSXP, GS->n));
    /* FIXME: speedup doing REAL(.) outside of loop :*/
    for (i = 0; i < GS->n; i++) {
	w[i] = GS->wts[i] * REAL(dmu_deta)[i]/sqrt(REAL(var)[i]);
	z[i] = REAL(GS->eta)[i] - GS->offset[i] +
	    (y[i] - REAL(GS->mu)[i])/REAL(dmu_deta)[i];
    }
    UNPROTECT(2);
    mer_update_ZXy(GS->mer);
}

/**
 * Iterate to determine the conditional modes of the random effects.
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
    SEXP fixef = GET_SLOT(GS->mer, lme4_fixefSym);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(GS->mer, lme4_LSym));
    int i;
    double crit = GS->tol + 1, *etap = REAL(GS->eta);

    if (varc)	  /* skip this step if varc == (double*) NULL */
	internal_mer_coefGets(GS->mer, varc, 2);
    Memcpy(REAL(fixef), fixed, LENGTH(fixef));
    internal_mer_Zfactor(GS->mer, L);
    internal_mer_ranef(GS->mer);
    internal_mer_fitted(GS->mer, GS->offset, etap);
    Memcpy(GS->etaold, etap, GS->n);

    for (i = 0; i < GS->maxiter && crit > GS->tol; i++) {
	internal_glmer_reweight(GS);
	internal_mer_Zfactor(GS->mer, L);
	internal_mer_ranef(GS->mer);
	internal_mer_fitted(GS->mer, GS->offset, etap);
	crit = conv_crit(GS->etaold, etap, GS->n);
    }
    internal_mer_Xfactor(GS->mer);
    Free(L);
    return (crit > GS->tol) ? 0 : i;
}

/**
 * Evaluate the conditional deviance for the stored random effects.
 *
 * @param GS Pointer to a GlmerStruct
 *
 * @return conditional deviance
 */
static double
random_effects_deviance(GlmerStruct GS)
{
    int i;
    double ans, *devs;

    internal_mer_fitted(GS->mer, GS->offset, REAL(GS->eta));
    eval_check_store(GS->linkinv, GS->rho, GS->mu);
    devs = REAL(PROTECT(eval_check(GS->dev_resids, GS->rho, REALSXP, GS->n)));
    for (i = 0, ans = 0; i < GS->n; i++) ans += devs[i];
    UNPROTECT(1);
    return ans;
}

/**
 * Calculate the deviance for a generalizedlinear mixed model at
 * arbitrary parameter values.  This version restores the original
 * values of the fixef and ranef slots after evaluating at arbitrary
 * beta and b.
 *
 * @param GS a generalized mixed-effects model pointer
 * @param beta fixed-effects parameter vector
 * @param b random-effects vector
 *
 * @return deviance
 */
static
double glmm_deviance(GlmerStruct GS, const double beta[], const double b[])
{
    SEXP x = GS->mer;
    SEXP fixefp = GET_SLOT(x, lme4_fixefSym),
	ranefp = GET_SLOT(x, lme4_ranefSym);
    int p = LENGTH(fixefp), q = LENGTH(ranefp);
    double *fixcp = Memcpy(Calloc(p, double), REAL(fixefp), p),
	*rancp = Memcpy(Calloc(q, double), REAL(ranefp), q), ans;

    mer_factor(x);
    Memcpy(REAL(fixefp), beta, p);
    Memcpy(REAL(ranefp), b, q);
    ans = random_effects_deviance(GS) +
	b_quadratic(b, GET_SLOT(x, lme4_OmegaSym),
		    INTEGER(GET_SLOT(x, lme4_GpSym)),
		    INTEGER(GET_SLOT(x, lme4_ncSym)));
    Memcpy(REAL(fixefp), fixcp, p);
    Memcpy(REAL(ranefp), rancp, q);
    Free(fixcp); Free(rancp);
    return ans;
}


/* Externally accessible functions */
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
    int i; double crit, *etap = REAL(GS->eta);

    Memcpy(GS->etaold, etap, GS->n);
    for (i = 0, crit = GS->tol + 1;
	 i < GS->maxiter && crit > GS->tol; i++) {
	internal_glmer_reweight(GS);
	if (!i) mer_initial(GS->mer); /* initialize first fit */
	internal_ECMEsteps(GS->mer, i ? 2 : GS->niterEM,
			   GS->EMverbose);
	eval(GS->LMEopt, GS->rho);
	internal_mer_fitted(GS->mer, GS->offset, etap);
	crit = conv_crit(GS->etaold, etap, GS->n);
    }
    if (crit > GS->tol)
	warning(_("IRLS iterations for PQL did not converge"));

    return R_NilValue;
}

/**
 * Compute the Laplace approximation to the deviance.
 *
 * @param pars pointer to a numeric vector of parameters
 * @param GSp pointer to a GlmerStruct object
 *
 * @return the Laplace approximation to the deviance
 */
SEXP glmer_devLaplace(SEXP pars, SEXP GSp)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSp);
    SEXP Omega = GET_SLOT(GS->mer, lme4_OmegaSym);
    int *Gp = INTEGER(GET_SLOT(GS->mer, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(GS->mer, lme4_ncSym));
    double *bhat = REAL(GET_SLOT(GS->mer, lme4_ranefSym)),
	*dcmp = REAL(GET_SLOT(GS->mer, lme4_devCompSym)),
	*dev = REAL(GET_SLOT(GS->mer, lme4_devianceSym));

    if (!isReal(pars) || LENGTH(pars) != GS->npar)
	error(_("`%s' must be a numeric vector of length %d"),
	      "pars", GS->npar);
    if (!internal_bhat(GS, REAL(pars), REAL(pars) + (GS->p)))
	return ScalarReal(DBL_MAX);
    dev[0] = dcmp[4] - dcmp[5] + random_effects_deviance(GS) +
	b_quadratic(bhat, Omega, Gp, nc);
    dev[1] = NA_REAL;
    return ScalarReal(dev[0]);
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

    Free(GS->offset); Free(GS->wts); Free(GS->etaold);
    Free(GS);
    return R_NilValue;
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
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 4, 0)
    GS->mer = find_and_check(rho, install("mer"), S4SXP, 0);
#else
    GS->mer = find_and_check(rho, install("mer"), VECSXP, 0);
#endif /* S4SXP */
    y = GET_SLOT(GS->mer, lme4_ySym);
    GS->n = LENGTH(y);
    GS->p = LENGTH(GET_SLOT(GS->mer, lme4_rXySym));
    GS->y = Memcpy(Calloc(GS->n, double), REAL(y), GS->n);
    Ztx = GET_SLOT(GET_SLOT(GS->mer, lme4_ZtSym), lme4_xSym);
    GS->eta = find_and_check(rho, install("eta"), REALSXP, GS->n);
    GS->mu = find_and_check(rho, install("mu"), REALSXP, GS->n);
    tmp = find_and_check(rho, install("offset"), REALSXP, GS->n);
    GS->offset = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    tmp = find_and_check(rho, install("weights"), REALSXP, GS->n);
    GS->wts = Memcpy(Calloc(GS->n, double), REAL(tmp), GS->n);
    GS->etaold = Calloc(GS->n, double);
    GS->cv = find_and_check(rho, install("cv"), VECSXP, 0);
    GS->niterEM = asInteger(internal_getElement(GS->cv, "niterEM"));
    GS->EMverbose = asLogical(internal_getElement(GS->cv, "EMverbose"));
    GS->tol = asReal(internal_getElement(GS->cv, "tolerance"));
    GS->maxiter = asInteger(internal_getElement(GS->cv, "maxIter"));
    GS->nf = LENGTH(GET_SLOT(GS->mer, lme4_flistSym));
    GS->npar = GS->p +
	coef_length(GS->nf, INTEGER(GET_SLOT(GS->mer, lme4_ncSym)));

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

/* MCMCsamp for generalized linear mixed models.
 *
 * 1) Detect if there  is a scale factor or not (not done yet).
 * 2) Sample beta and b according to the normal approximation
 * 3) Evaluate the Metropolis-Hastings ratio and accept or reject
 * 4) If step is accepted then reweight/update
 * 5) Sample from the Wishart for the variance matrix.
 * 6) If necessary, sample from the scale factor distribution (not done yet).
*/

/**
 * Create a Markov Chain Monte Carlo sample from a fitted generalized
 * linear mixed model
 *
 * @param GSpt External pointer to a GlmerStruct
 * @param savebp Logical indicator of whether or not to save the
 *   random effects in the MCMC sample
 * @param nsampp number of samples to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP
glmer_MCMCsamp(SEXP GSpt, SEXP savebp, SEXP nsampp, SEXP transp, SEXP verbosep)
{
    GlmerStruct GS = (GlmerStruct) R_ExternalPtrAddr(GSpt);
    SEXP ans, x = GS->mer;
    SEXP Omega = GET_SLOT(x, lme4_OmegaSym),
	Omegacp = PROTECT(duplicate(Omega));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*nc = INTEGER(GET_SLOT(x, lme4_ncSym)),
	i, j, nf = LENGTH(Omega), nsamp = asInteger(nsampp),
	p = LENGTH(GET_SLOT(x, lme4_rXySym)),
	q = LENGTH(GET_SLOT(x, lme4_rZySym)),
	saveb = asLogical(savebp),
	trans = asLogical(transp),
	verbose = asLogical(verbosep);
    double
	*RXX = REAL(GET_SLOT(GET_SLOT(x, lme4_RXXSym), lme4_xSym)),
	*RZX = REAL(GET_SLOT(GET_SLOT(x, lme4_RZXSym), lme4_xSym)),
	*bhat = REAL(GET_SLOT(x, lme4_ranefSym)),
	*betahat = REAL(GET_SLOT(x, lme4_fixefSym)),
	*bcur = Calloc(q, double), *betacur = Calloc(p, double), /* current */
	*bnew = Calloc(q, double), *betanew = Calloc(p, double), /* proposed */
	*ansp, MHratio;
    int nrbase = p + 1 + coef_length(nf, nc); /* rows always included */
    int nrtot = nrbase + (saveb ? q : 0);

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    GetRNGstate();
    /* initialize current values of b and beta to estimates */
    Memcpy(betacur, REAL(GET_SLOT(x, lme4_fixefSym)), p);
    Memcpy(bcur, REAL(GET_SLOT(x, lme4_ranefSym)), q);
    if(verbose) Rprintf("%12s\n", "MHratio");
    for (i = 0; i < nsamp; i++) {
	double odev, ndev, ogaus, ngaus, *col = ansp + i * nrtot;
				/* update bhat */
/* 	internal_bhat(GS, betahat, (double *) NULL); */
				/* generate proposal for b and beta */
				/* save Gaussian deviance at proposal */
	ngaus = internal_betab_update(p, q, 1, L, RZX, RXX,
				      betahat, bhat, betanew, bnew);
				/* Gaussian deviance at current */
	ogaus = internal_Gaussian_deviance(p, q, L, RZX, RXX, betahat, bhat,
					 betacur, bcur);
				/*  glmm deviance at proposal */
	ndev = glmm_deviance(GS, betanew, bnew);
				/*  glmm deviance at current */
	odev = glmm_deviance(GS, betacur, bcur);
				/* evaluate Metropolis-Hastings ratio */
	MHratio = exp((ngaus - ogaus - ndev + odev)/2);
	if(verbose) Rprintf("%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n",
			    ogaus, ngaus, odev, ndev, ngaus - ogaus, odev - ndev,
			    MHratio);
	/* Accept proposal with probability min(MHratio, 1) */
	if (unif_rand() < MHratio) {
	    Memcpy(betacur, betanew, p);
	    Memcpy(bcur, bnew, q);
	}
				/* Store beta */
	for (j = 0; j < p; j++) col[j] = betacur[j];
				/* Optionally store b */
	if (saveb) for (j = 0; j < q; j++) col[nrbase + j] = bcur[j];
				/* Update and store Omega */
	internal_Omega_update(Omega, bcur, 1, nf, nc, Gp, col + p, trans);
	internal_mer_refactor(x);
	mer_secondary(x);

	col[nrbase - 1] = glmm_deviance(GS, betacur, bcur); /* store deviance */
    }
    PutRNGstate();
    Free(bcur); Free(betacur); Free(betanew); Free(bnew);
				/* Restore original Omega */
    SET_SLOT(x, lme4_OmegaSym, Omegacp);
    internal_bhat(GS, betahat, (double *) NULL);
    internal_mer_refactor(x);

    Free(L);
    UNPROTECT(2);
    return ans;
}
