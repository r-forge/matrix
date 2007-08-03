#include "lmer.h"

#if R_VERSION < R_Version(2,6,0)

/**
 * Allocate a 3-dimensional array
 *
 * @param mode The R mode (e.g. INTSXP)
 * @param nrow number of rows
 * @param ncol number of columns
 * @param nface number of faces
 *
 * @return A 3-dimensional array of the indicated dimensions and mode
 */
static SEXP alloc3DArray(SEXPTYPE mode, int nrow, int ncol, int nface)
{
    SEXP s, t;
    int n;

    if (nrow < 0 || ncol < 0 || nface < 0)
	error(_("negative extents to 3D array"));
    if ((double)nrow * (double)ncol * (double)nface > INT_MAX)
	error(_("alloc3Darray: too many elements specified"));
    n = nrow * ncol * nface;
    PROTECT(s = allocVector(mode, n));
    PROTECT(t = allocVector(INTSXP, 3));
    INTEGER(t)[0] = nrow;
    INTEGER(t)[1] = ncol;
    INTEGER(t)[2] = nface;
    setAttrib(s, R_DimSymbol, t);
    UNPROTECT(2);
    return s;
}

#endif

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SET_SLOT(obj, nm, allocVector(type, length));
    return GET_SLOT(obj, nm);
}

/* Functions for the lmer representation */

				/* positions in the deviance vector */
enum devP {ML_POS=0, REML_POS, ldL2_POS, ldRX2_POS, lpdisc_POS, disc_POS, bqd_POS};
/* ("ML", "REML", "ldL2", "ldRX2", "lpdisc", "disc", "bqd") */
				/* positions in the dims vector */
enum dimP {nf_POS=0, n_POS, p_POS, q_POS, s_POS, np_POS, isREML_POS, fTyp_POS, nest_POS};
/* ("nf", "n", "p", "q", "s", "np", "REML", "fTyp", "nest") */

#define isREML(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isREML_POS]
#define L_SLOT(x) AS_CHM_FR(GET_SLOT(x, lme4_LSym))

/**
 * Return the sum of squares of the first n elements of x
 *
 * @param n
 * @param x
 *
 * @return sum of squares
 */
static double R_INLINE lme4_sumsq(const double *x, int n)
{
    double ans = 0;
    int i;

    for (i = 0; i < n; i++) ans += x[i] * x[i];
    return ans;
}

/**
 * Evaluate the logarithm of the square of the determinant of L
 * (i.e. the logarithm of the determinant of LL')
 *
 * @param L
 *
 * @return logarithm of the square of the determinant of L
 */
static double
chm_log_det2(CHM_FR L)
{
    double ans = 0;
    int i;

    if (L->is_super) {
	for (i = 0; i < L->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(L->pi))[i + 1] - ((int *)(L->pi))[i],
		nc = ((int *)(L->super))[i + 1] - ((int *)(L->super))[i];
	    double *x = (double *)(L->x) + ((int *)(L->px))[i];

	    for (j = 0; j < nc; j++) {
		ans += 2 * log(fabs(x[j * nrp1]));
	    }
	}
    } else {
	int *li = (int*)(L->i), *lp = (int*)(L->p), j, p;
	double *lx = (double *)(L->x);
	
	for (j = 0; j < L->n; j++) {
	    for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++) {};
	    if (li[p] != j) break; /* what happened to the diagonal element? */
	    ans += log(lx[p] * ((L->is_ll) ? lx[p] : 1.));
	}
    }
    return ans;
}

/**
 * Populate the st, nc and nlev arrays.  Return the maximum element of nc.
 *
 * @param ST pointer to a list (length nf) of matrices
 * @param Gp pointer to integer vector (length nf + 1) of group
 * pointers
 * @param st length nf array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST.  Not used if NULL.
 * @param nc length nf array of ints to be filled with the number of
 * columns
 * @param nlev length nf array of ints to be filled with the number of
 * levels of the grouping factor for each term
 * 
 * @return maximum element of nc
 */
static int			/* populate the st, nc and nlev arrays */
ST_nc_nlev(const SEXP ST, const int *Gp, double **st, int *nc, int *nlev)
{
    int ans = 0, i, nf = LENGTH(ST);

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int nci = *INTEGER(getAttrib(STi, R_DimSymbol));

	if (nci > ans) ans = nci;
	if (st) st[i] = REAL(STi);
	nc[i] = nci;
	nlev[i] = (Gp[i + 1] - Gp[i])/nci;
    }
    return ans;
}

/**
 * Return the group in the (nf, Gp) combination to which ind belongs
 *
 * @param ind a row number
 * @param nf number of groups of rows
 * @param Gp group pointers
 *
 */
static int R_INLINE Gp_grp(int ind, int nf, const int *Gp)
{
    int i;
    for (i = 0; i < nf; i++) if (ind < Gp[i + 1]) return i;
    error(_("invalid row index %d (max is %d)"), ind, Gp[nf]);
    return -1;			/* -Wall */
}

	
/**
 * Determine the nonzero positions in the jth column of Vt
 *
 * @param nz array to hold the answer
 * @param j column index
 * @param nf number of groups of rows
 * @param Gp group pointers
 * @param zi row indices in Zt
 * @param zp column pointers for Zt
 *
 * @return count of nonzeros in the jth column of Vt
 */
static int
Vt_nz_col(int *nz, int j, int nf, const int *Gp, const int *nc,
       const int *nlev, const int *zi, const int *zp)
{
    int ans, i, p;

    AZERO(nz, Gp[nf]);
    for (p = zp[j]; p < zp[j + 1]; p++) {
	int zrow = zi[p];
	int k = Gp_grp(zrow, nf, Gp);

	nz[zrow] = 1;		/* T contains the identity */
	if (nc[k] > 1) {
	    int nextra = (zrow - Gp[k]) / nlev[k];
	    for (i = 1; i <= nextra; i++)
		nz[zrow - i * nlev[k]] = 1;
	}
    }
    for (i = 0, ans = 0; i < Gp[nf]; i++) if (nz[i]) ans++;
    return ans;
}

/**
 * Create the Vt matrix pattern from Zt, ST and Gp.  Partition the
 * columns into s groups and overlay them.
 * 
 *
 * @param Zt 
 * @param ST
 * @param GpP pointer to the Gp array
 *
 * @return Vt
 */
SEXP mer_create_Vt(SEXP Zt, SEXP ST, SEXP GpP)
{
    SEXP ans;
    int *Gp = INTEGER(GpP), *nnz, *nz, *vi, *vp, *zdims, *zi, *zp,
	Vnnz, ZtOK, j, nf = LENGTH(ST);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    R_CheckStack();

				/* Trivial case, all nc == 1 */
    if (ST_nc_nlev(ST, Gp, (double**)NULL, nc, nlev) <= 1)
	return duplicate(Zt);
				/* Check the nonzero pattern in Zt */
    zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym));
    nnz = Alloca(zdims[1], int);
    nz = Alloca(zdims[0], int);
    R_CheckStack();

    zi = INTEGER(GET_SLOT(Zt, lme4_iSym));
    zp = INTEGER(GET_SLOT(Zt, lme4_pSym));
    for (j = 0, ZtOK = 1; j < zdims[1]; j++) {
	nnz[j] = Vt_nz_col(nz, j, nf, Gp, nc, nlev, zi, zp);
	if (nnz[j] != (zp[j + 1] - zp[j])) ZtOK = 0;
    }
    if (ZtOK) return duplicate(Zt);
				/* Must create a new dgCMatrix object */ 
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    SET_SLOT(ans, lme4_DimSym, duplicate(GET_SLOT(Zt, lme4_DimSym)));
    SET_SLOT(ans, lme4_DimNamesSym, allocVector(VECSXP, 2));
				/* create and evaluate the p slot */
    vp = INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, zdims[1] + 1));
    vp[0] = 0;
    for (j = 0; j < zdims[1]; j++) vp[j + 1] = vp[j] + nnz[j];
    Vnnz = vp[zdims[1]];
    vi = INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, Vnnz));
    AZERO(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, Vnnz)), Vnnz);
    for (j = 0; j < zdims[1]; j++) { /* fill in the i slot */
	int i, pos = vp[j];
	Vt_nz_col(nz, j, nf, Gp, nc, nlev, zi, zp);
	for (i = 0; i < zdims[0]; i++) if (nz[i]) vi[pos++] = i;
    }

    UNPROTECT(1); 
    return ans;
}

/**
 * Change the x slot in Vt according to ST'Zt
 *
 * @param x An object that inherits from class mer
 */
SEXP mer_update_Vt(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym));
    CHM_SP cVt = AS_CHM_SP(GET_SLOT(x, lme4_VtSym)),
	cZt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym));
    int *vi = (int*)(cVt->i), *vp = (int*)(cVt->p),
	*zi = (int*)(cZt->i), *zp = (int*)(cZt->p),
	j, nf = LENGTH(ST), p;
    int vnnz = vp[cVt->ncol], znnz = zp[cZt->ncol];
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*),
	*vx = (double*)(cVt->x), *zx = (double*)(cZt->x),
	one[] = {1,0};
    R_CheckStack();

    if (vnnz == znnz) Memcpy(vx, zx, znnz);
    else {
	AZERO(vx, vnnz);
	for (j = 0; j < cVt->ncol; j++) {
	    int pv = vp[j];
	    for (p = zp[j]; p < zp[j + 1]; p++) {
		while (vi[pv] < zi[p]) pv++;
		if (vi[pv] != zi[p])
		    error(_("nonconforming Zt and Vt structures, j = %d"), j);
		vx[pv] = zx[p];
	    }
	}
    }
	
    if (ST_nc_nlev(ST, Gp, st, nc, nlev) > 1) /* T' == I when ncmax == 1 */
	for (j = 0; j < cVt->ncol; j++) /* multiply column j by T' */
	    for (p = vp[j]; p < vp[j + 1];) {
		int i = Gp_grp(vi[p], nf, Gp);

		if (nc[i] <= 1) p++;
		else {
		    int nr = p;	/* number of rows in `B' in dtrmm call */
		    while ((vi[nr] - Gp[i]) < nlev[i]) nr++;
		    nr -= p;	/* nr == 1 except in models with carry-over */
		    F77_CALL(dtrmm)("R", "L", "N", "U", &nr, nc + i,
				    one, st[i], nc + i, vx + p, &nr);
		    p += (nr * nc[i]);
		}
	    }
    for (p = 0; p < vnnz; p++) {
	int i = Gp_grp(vi[p], nf, Gp);
	vx[p] *= st[i][((vi[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
    return R_NilValue;
}

/**
 * dest = P  %*% S %*% t(T) %*% src
 *
 * @param dest matrix whose contents are overwritten
 * @param ST ST slot
 * @param Gp group pointers
 * @param src originating numeric matrix
 * @param Perm permutation to be applied
 *
 */
static void
TSPp_dense_mult(double *dx, int nf, int m, int n, double **st,
		const int *nc, const int *nlev, const int *Gp,
		const double *sx, const int *Perm)
{
    int i, j, k;
    double *tmp = Alloca(m, double), one = 1;
    R_CheckStack(); 

    for (j = 0; j < n; j++) {
				/* apply permutation if given */
	if (Perm) for (i = 0; i < m; i++) tmp[Perm[i]] = sx[j * m + i];
	else Memcpy(tmp, sx + j * m, m);
	for (i = 0; i < nf; i++) {
	    for (k = 0; k < nc[i]; k++) { /* multiply by \tilde{S}_i */
		double dd = st[i][k * (nc[i] + 1)];
		int base = Gp[i] + k * nlev[i], kk;
		for (kk = 0; kk < nlev[i]; kk++) tmp[base + kk] *= dd;
	    }
	    if (nc[i] > 1) {	/* multiply by \tilde{T}_i */
		F77_CALL(dtrmm)("R", "L", "T", "U", nlev + i, nc + i, &one,
				st[i], nc + i, tmp + Gp[i], nlev + i);
	    }
	}
	Memcpy(dx + j * m, tmp, m);
    }
}

static void
TSPp_SEXP_mult(SEXP dest, SEXP ST, const int *Gp, SEXP src, const int *Perm)
{
    int isM = isMatrix(src), nf = LENGTH(ST);
    int m = isM ? INTEGER(getAttrib(src, R_DimSymbol))[0] : LENGTH(src),
	n = isM ? INTEGER(getAttrib(src, R_DimSymbol))[1] : 1;
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack(); 

    ST_nc_nlev(ST, Gp, st, nc, nlev);
    TSPp_dense_mult(REAL(dest), nf, m, n, st, nc, nlev, Gp, REAL(src), Perm);
}

/**
 * Update the ranef slot from the uvec slot.
 *
 * @param x An object that inherits from class lmer
 */
SEXP mer_update_b(SEXP x)
{
    TSPp_SEXP_mult(GET_SLOT(x, lme4_ranefSym),
		   GET_SLOT(x, lme4_STSym),
		   INTEGER(GET_SLOT(x, lme4_GpSym)),
		   GET_SLOT(x, lme4_uvecSym),
		   INTEGER(GET_SLOT(GET_SLOT(x, lme4_LSym),
				    lme4_permSym)));
    return R_NilValue;
}

/**
 * Evaluate starting estimates for the elements of ST
 *
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param Gpp length nf+1 vector of group pointers for the rows of Zt
 * @param Zt transpose of Z matrix
 *
 */
SEXP ST_initialize(SEXP ST, SEXP Gpp, SEXP Zt)
{
    int *Gp = INTEGER(Gpp),
	*Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	i, j, k, nf = LENGTH(ST);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int),
	nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = Alloca(Zdims[0], double),
	**st = Alloca(nf, double*),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    
    ST_nc_nlev(ST, Gp, st, nc, nlev);
    AZERO(rowsqr, Zdims[0]);
    for (i = 0; i < nnz; i++) rowsqr[zi[i]] += zx[i] * zx[i];
    for (i = 0; i < nf; i++) {
	AZERO(st[i], nc[i] * nc[i]);
	for (j = 0; j < nc[i]; j++) {
	    double *stij = st[i] + j * (nc[i] + 1);
	    for (k = 0; k < nlev[i]; k++)
		*stij += rowsqr[Gp[i] + j * nlev[i] + k];
	    *stij = sqrt(nlev[i]/(0.375 * *stij));
	}
    }
    return R_NilValue;
}

/**
 * dest = P  %*% S %*% t(T) %*% src
 *
 * @param dest matrix whose contents are overwritten
 * @param ST ST slot
 * @param Gp group pointers
 * @param src originating numeric matrix
 * @param Perm permutation to be applied
 *
 */
static void
PSTp_dense_mult(SEXP dest, SEXP ST, const int *Gp, SEXP src, int *Perm)
{
    double *dx = REAL(dest), *sx = REAL(src);
    int *dims = INTEGER(getAttrib(src, R_DimSymbol)), i, j, k;
    int m = dims[0], n = dims[1], ncmax, nf = LENGTH(ST);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*), *tmp = Alloca(m, double), one = 1;
    R_CheckStack(); 

    ncmax = ST_nc_nlev(ST, Gp, st, nc, nlev);
    for (j = 0; j < n; j++) {
	Memcpy(tmp, sx + j * m, m);
	for (i = 0; i < nf; i++) {
	    if (nc[i] > 1) {	/* multiply by \tilde{T}_i' */
		F77_CALL(dtrmm)("R", "L", "N", "U", nlev + i, nc + i, &one,
				st[i], nc + i, tmp + Gp[i], nlev + i);
	    }
	    for (k = 0; k < nc[i]; k++) { /* multiply by \tilde{S}_i */
		double dd = st[i][k * (nc[i] + 1)];
		int base = Gp[i] + k * nlev[i], kk;
		for (kk = 0; kk < nlev[i]; kk++) tmp[base + kk] *= dd;
	    }
	}
				/* apply permutation if given */
	if (Perm) for (i = 0; i < m; i++) dx[j * m + i] = tmp[Perm[i]];
	else Memcpy(dx + j * m, tmp, m);
    }
}

/**
 * Evaluate the elements of the deviance slot in an lmer object
 *
 * @param x pointer to an lmer object
 *
 */
SEXP
lmer_update_dev(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), j, info;
    double *d = REAL(GET_SLOT(x, lme4_devianceSym)),
	*rxy = REAL(GET_SLOT(x, lme4_RXySym)), dn = (double)(dims[n_POS]),
	dnmp = (double)(dims[n_POS] - dims[p_POS]), mone = -1, one = 1;
    SEXP RVXy = GET_SLOT(x, lme4_RVXySym);
    CHM_SP cVt = AS_CHM_SP(GET_SLOT(x, lme4_VtSym));
    CHM_DN cRVXy = AS_CHM_DN(RVXy), ans;
    CHM_FR L = L_SLOT(x);
    int pp1 = cRVXy->ncol, q = cRVXy->nrow;
    R_CheckStack();

    mer_update_Vt(x);
    j = c.final_ll;
    c.final_ll = L->is_ll;
    if (!M_cholmod_factorize_p(cVt, &one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    c.final_ll = j;
    d[ldL2_POS] = chm_log_det2(L);
				/* Evaluate PST'Z'[X:y] in RVXy */
    PSTp_dense_mult(RVXy, GET_SLOT(x, lme4_STSym),
		    INTEGER(GET_SLOT(x, lme4_GpSym)),
		    GET_SLOT(x, lme4_ZtXySym), (int*)(L->Perm));
				/* solve for RVXy */
    ans = M_cholmod_solve(CHOLMOD_L, L, cRVXy, &c);
    Memcpy(REAL(RVXy), (double*)(ans->x), q * pp1);
    M_cholmod_free_dense(&ans, &c);
    				/* downdate XytXy  */
    F77_CALL(dlacpy)("U", &pp1, &pp1, REAL(GET_SLOT(x, lme4_XytXySym)),
		     &pp1, rxy, &pp1);
    F77_CALL(dsyrk)("U", "T", &pp1, &q, &mone, REAL(RVXy),
		    &q, &one, rxy, &pp1);
    F77_CALL(dpotrf)("U", &pp1, rxy, &pp1, &info);
    if (info)
	error(_("the leading minor of order %d of the downdated [X:y]'[X:y] is not positive definite"),
	      info);
    d[lpdisc_POS] = 2 * log(rxy[pp1 * pp1 - 1]);
    for (j = 0, d[ldRX2_POS] = 0; j < (pp1 - 1); j++)
	d[ldRX2_POS] += 2 * log(rxy[j * (pp1 + 1)]);
    d[ML_POS] = d[ldL2_POS] + dn * (1. + d[lpdisc_POS] + log(2. * PI / dn));
    d[REML_POS] = d[ldL2_POS] + d[ldRX2_POS] + dnmp *
	(1. + d[lpdisc_POS] + log(2. * PI / dnmp));
    d[bqd_POS] = d[disc_POS] = NA_REAL;
    return R_NilValue;
}

/**
 * Extract the parameters from ST list
 *
 * @param ST ST slot from an mer object
 * @param pars vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double
*internal_ST_getPars(SEXP ST, double *pars)
{
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++)
	    pars[pos++] = st[j * ncp1];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
}

/**
 * Extract the parameters from the ST slot of an lmer object
 *
 * @param x an lmer object
 *
 * @return pointer to a REAL vector
 */
SEXP ST_getPars(SEXP x)
{
    int np = INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS];
    SEXP ST = GET_SLOT(x, lme4_STSym);
    SEXP ans = PROTECT(allocVector(REALSXP, np));

    internal_ST_getPars(ST, REAL(ans));
    UNPROTECT(1); 
    return ans;
}

/**
 * Update the ST slot of an lmer object
 *
 * @param pars double vector of the appropriate length
 * @param ST ST slot from an lmer object
 *
 * @return pointer to the updated ST object
 */
static SEXP
ST_SetPars(const double *pars, SEXP ST)
{
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++)
	    st[j * ncp1] = pars[pos++];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		st[k + j * nci] = pars[pos++];
    }
    return ST;
}

/**
 * Update the ST slot of an lmer object from a REAL vector of
 * parameters and update the Cholesky factorization
 *
 * @param x an lmer object
 * @param pars a REAL vector of the appropriate length
 *
 */
SEXP ST_setPars(SEXP x, SEXP pars)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int npar = INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS];

    if (!isReal(pars) || LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    ST_SetPars(REAL(pars), ST);
    return R_NilValue;
}

/**
 * Update the contents of the fixef, ranef and uvec slots in an lmer
 * object.
 *
 * @param x an lmer object
 *
 * @return R_NilValue
 */
SEXP lmer_update_effects(SEXP x)
{
    SEXP uvec = GET_SLOT(x, lme4_uvecSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), ione = 1;
    int p = dims[p_POS],  pp1 = dims[p_POS] + 1, q = dims[q_POS];
    double *RXy = REAL(GET_SLOT(x, lme4_RXySym)),
	*RVXy = REAL(GET_SLOT(x, lme4_RVXySym)),
	*fixef = REAL(GET_SLOT(x, lme4_fixefSym)),
	*u = REAL(uvec), mone = -1, one = 1;
    CHM_FR L = L_SLOT(x);
    CHM_DN cu = AS_CHM_DN(uvec), sol;
    R_CheckStack();

    lmer_update_dev(x);
    Memcpy(fixef, RXy + p * pp1, p);
    F77_CALL(dtrsv)("U", "N", "N", &p, RXy, &pp1, fixef, &ione);
    Memcpy(u, RVXy + p * q, q);
    F77_CALL(dgemv)("N", &q, &p, &mone, RVXy, &q, fixef, &ione, &one, u, &ione);
    if (!(sol = M_cholmod_solve(CHOLMOD_Lt, L, cu, &c)))
	error(_("cholmod_solve (CHOLMOD_Lt) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    Memcpy(u, (double*)(sol->x), q);
    M_cholmod_free_dense(&sol, &c);
    mer_update_b(x);
    return R_NilValue;
}

/**
 * Return the REML or ML conditional estimate of sigma, the standard
 * deviation of the per-observation noise term.
 *
 * @param REML non-zero for REML estimate, 0 for ML estimate
 * @param dims vector of dimensions
 * @param deviance vector of deviance components
 */
static R_INLINE double
Mer_sigma(int REML, const int* dims, const double* deviance)
{
    return sqrt(exp(deviance[lpdisc_POS])/
		((double)(dims[n_POS] - (REML ? dims[p_POS] : 0))));
}


/**
 * Extract the estimate of the scale factor from an lmer object
 *
 * @param x an lmer object
 * @param which scalar integer (< 0 => REML, 0 => native, > 0 => ML)
 *
 * @return scalar REAL value
 */
SEXP mer_sigma(SEXP x, SEXP which)
{
    int w = asInteger(which);
		
    return ScalarReal(Mer_sigma(w < 0 || (!w && isREML(x)),
				INTEGER(GET_SLOT(x, lme4_dimsSym)),
				REAL(GET_SLOT(x, lme4_devianceSym))));
}

/**
 * Extract the posterior variances of the random effects in an lmer object
 *
 * @param x pointer to an lmer object
 *
 * @return pointer to a list of arrays
 */
SEXP mer_postVar(SEXP x, SEXP uS)
{
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int nf = dims[nf_POS], q = dims[q_POS];
    SEXP ans = PROTECT(allocVector(VECSXP, nf));
    int i, j, k;
    double *vv, one[] = {1,0},
	sc = asLogical(uS) ?
	     Mer_sigma(isREML(x), dims, REAL(GET_SLOT(x, lme4_devianceSym)))
	    : 1;
    CHM_SP sm1, sm2;
    CHM_DN dm1;
    CHM_FR L = L_SLOT(x);
    int *Perm = (int*)(L->Perm), *iperm = Alloca(q, int),
	*nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    double **st = Alloca(nf, double*);
    R_CheckStack();

    ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
    for (j = 0; j < q; j++) iperm[Perm[j]] = j; /* inverse permutation */
    for (i = 0; i < nf; i++) {
	int ncisqr = nc[i] * nc[i];
	CHM_SP rhs = M_cholmod_allocate_sparse(q, nc[i], nc[i],
					       1/*sorted*/, 1/*packed*/,
					       0/*stype*/, CHOLMOD_REAL, &c);

	SET_VECTOR_ELT(ans, i, alloc3DArray(REALSXP, nc[i], nc[i], nlev[i]));
	vv = REAL(VECTOR_ELT(ans, i));
	for (j = 0; j <= nc[i]; j++) ((int *)(rhs->p))[j] = j;
	for (j = 0; j < nc[i]; j++)
	    ((double *)(rhs->x))[j] = st[i][j * (nc[i] + 1)] * sc;
	for (k = 0; k < nlev[i]; k++) {
	    double *vvk = vv + k * ncisqr;
	    for (j = 0; j < nc[i]; j++)
		((int*)(rhs->i))[j] = iperm[Gp[i] + j * nlev[i]];
	    sm1 = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
	    sm2 = M_cholmod_transpose(sm1, 1 /*values*/, &c);
	    M_cholmod_free_sparse(&sm1, &c);
	    sm1 = M_cholmod_aat(sm2, (int*)NULL, (size_t)0, 1 /*mode*/, &c);
	    dm1 = M_cholmod_sparse_to_dense(sm1, &c);
	    M_cholmod_free_sparse(&sm1, &c); M_cholmod_free_sparse(&sm2, &c);
	    Memcpy(vvk, (double*)(dm1->x), ncisqr);
	    M_cholmod_free_dense(&dm1, &c);
	    if (nc[i] > 1) {
		F77_CALL(dtrmm)("L", "L", "N", "U", nc + i, nc + i,
				one, st[i], nc + i, vvk, nc + i);
		F77_CALL(dtrmm)("R", "L", "T", "U", nc + i, nc + i,
				one, st[i], nc + i, vvk, nc + i);
	    }
	}
	M_cholmod_free_sparse(&rhs, &c);
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Create and initialize L
 *
 * @param Zt pointer to the random-effects model matrix (transposed)
 *
 * @return L
 */
SEXP mer_create_L(SEXP Vt)
{
    double one[] = {1, 0};
    int fll = c.final_ll;
    CHM_SP cVt = AS_CHM_SP(Vt);
    CHM_FR L;
    R_CheckStack();

    c.final_ll = 1;
    L = M_cholmod_analyze(cVt, &c);
    if (!M_cholmod_factorize_p(cVt, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    c.final_ll = fll;

    return M_chm_factor_to_SEXP(L, 1);
}

/**
 * Check validity of an object that inherits from the mer class.
 *
 * @param x Pointer to an lmer or glmer or nlmer object
 *
 * @return TRUE if the object is a valid mer object, otherwise a string
 *         that describes the violation.
 */
SEXP mer_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
	XP = GET_SLOT(x, lme4_XSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	fixefP = GET_SLOT(x, lme4_fixefSym),
	flistP = GET_SLOT(x, lme4_flistSym),
	offsetP = GET_SLOT(x, lme4_offsetSym),
	ranefP = GET_SLOT(x, lme4_ranefSym),
	uvecP = GET_SLOT(x, lme4_uvecSym),
	weightsP = GET_SLOT(x, lme4_weightsSym),
	y = GET_SLOT(x, lme4_ySym);
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int i, n = dd[n_POS], nf = dd[nf_POS], nq,
	p = dd[p_POS], q = dd[q_POS], s = dd[s_POS];
    CHM_SP Zt = AS_CHM_SP(GET_SLOT(x, lme4_ZtSym)),
	Vt =  AS_CHM_SP(GET_SLOT(x, lme4_VtSym));
    CHM_FR L = L_SLOT(x);
    R_CheckStack();
				/* check lengths */
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length nf"));
    if (LENGTH(GpP) != nf + 1)
	return mkString(_("Slot Gp must have length dims['nf'] + 1"));
    if (LENGTH(y) != n)
	return mkString(_("Slot y must have length dims['n']"));
    if (Gp[0] != 0 || Gp[nf] != q)
	return mkString(_("Gp[1] != 0 or Gp[dims['nf'] + 1] != dims['q']"));
    if (LENGTH(fixefP) != p)
	return mkString(_("Slot fixef must have length ['p']"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length dims['q']"));
    if (LENGTH(uvecP) != q)
	return mkString(_("Slot uvec must have length dims['q']"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or dims['n']"));
    if (LENGTH(offsetP) && LENGTH(offsetP) != n)
	return mkString(_("Slot offset must have length 0 or dims['n']"));
    if (LENGTH(devianceP) != (bqd_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (bqd_POS + 1))
	return mkString(_("deviance slot not named or incorrect length"));
    if (LENGTH(dimsP) != (nest_POS + 1) ||
	LENGTH(getAttrib(dimsP, R_NamesSymbol)) != (nest_POS + 1))
	return mkString(_("dims slot not named or incorrect length"));
    if (L->n != q || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size dims['q']"));
    if (Zt->nrow != q || Zt->ncol != n * s)
	return mkString(_("Slot Zt must by dims['q']  by dims['n']*dims['s']"));
    if (Vt->nrow != q || Vt->ncol != n * s)
	return mkString(_("Slot Zt must be dims['q']  by dims['n']*dims['s']"));
    dd = INTEGER(getAttrib(XP, R_DimSymbol)); 
    if (!isReal(XP) || dd[1] != p)
	return mkString(_("Slot X must be a numeric matrix with dims['p'] columns"));
    if (dd[0] != 0 && dd[0] != (n * s)) /* special case */
	return mkString(_("Slot X must be have 0 or dims['n']*dims['s'] rows"));

    nq = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i), fli = VECTOR_ELT(flistP, i);
	int *dm = INTEGER(getAttrib(STi, R_DimSymbol));
	if (!isMatrix(STi) || !isReal(STi) || dm[0] != dm[1])
	    return
		mkString(_("Slot ST must be a list of square numeric matrices"));
	if (Gp[i] > Gp[i + 1])
	    return mkString(_("Gp must be non-decreasing"));
	if (!isFactor(fli))
	    return mkString(_("flist must be a list of factors"));
	nq += dm[0] * LENGTH(getAttrib(fli, R_LevelsSymbol));
    }
    if (q != nq)
	return mkString(_("q is not sum of columns by levels"));

    return ScalarLogical(1);
}

#define BUF_SIZE 127

/**
 * Check the dimensions of the matrix pointer MP.
 *
 * @param buffer character buffer of length BUF_SIZE + 1
 * @param nm name of matrix - used in the error message
 * @param MP pointer to a matrix
 *
 * @return TRUE if the object is a valid lmer object, else a string
 * describing the nature of the violation.
 */
static char R_INLINE
*chkDims(char *buffer, char *nm, SEXP MP, int nr, int nc)
{
    int *dd = isMatrix(MP) ?
	INTEGER(getAttrib(MP, R_DimSymbol)) :
	(int *)NULL;

    if (!dd) error(_("Argument MP to chkDims is not a matrix"));
    buffer[0] = '\0';
    if (!isReal(MP) || dd[0] != nr || dd[1] != nc)
	snprintf(buffer, BUF_SIZE,
		 "Matrix %s must be a %d by %d numeric matrix",
		nm, nr, nc);
    return buffer;
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
    int *dd = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int p = dd[p_POS], q = dd[q_POS];
    int pp1 = p + 1;
    char *buf = Alloca(BUF_SIZE + 1, char);
    R_CheckStack();

    /* Note: We don't need to check the matrix slots to see if they are matrices. */
    if (strlen(chkDims(buf, "ZtXy", GET_SLOT(x, lme4_ZtXySym), q, pp1)))
	return mkString(buf);
    if (strlen(chkDims(buf, "XytXy", GET_SLOT(x, lme4_XytXySym), pp1, pp1)))
	return mkString(buf);
    if (strlen(chkDims(buf, "RXy", GET_SLOT(x, lme4_RXySym), pp1, pp1)))
	return mkString(buf);
    if (strlen(chkDims(buf, "RVXy", GET_SLOT(x, lme4_RVXySym), q, pp1)))
	return mkString(buf);
    return ScalarLogical(1);
}

/**
 * Check validity of an glmer object.
 *
 * @param x Pointer to an glmer object
 *
 * @return TRUE if the object is a valid glmer object, else a string
 *         describing the nature of the violation.
 */
SEXP glmer_validate(SEXP x)
{
    return ScalarLogical(1);
}

/**
 * Check validity of an nlmer object.
 *
 * @param x Pointer to an nlmer object
 *
 * @return TRUE if the object is a valid nlmer object, else a string
 *         describing the nature of the violation.
 */
SEXP nlmer_validate(SEXP x)
{
    if (LENGTH(GET_SLOT(x, lme4_muSym)) != LENGTH(GET_SLOT(x, lme4_ySym)))
	return mkString(_("Slot mu must have length dims['n']"));
    return ScalarLogical(1);
}

/* Functions common to nonlinear and generalized linear mixed models */

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


#define IRLS_MAXITER  40
#define IRLS_TOL      1e-8
#define IRLS_SMIN     1e-4

/**
 * Evaluate the convergence criterion and copy eta to
 * etaold
 *
 * @param eta current value of the linear predictors
 * @param etaold previous values of the linear predictors
 *
 * @return convergence criterion
 */
static double conv_crit(double etaold[], double eta[], int n)
{
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
 * Evaluate the linear predictor as model offset + X \beta + V u
 *
 * @param x pointer to a glmer or a nlmer object
 * @param ans pointer to memory to be overwritten with the values
 *
 * @return ans
 */
static double *nglmer_eta(SEXP x, double *ans)
{
    SEXP X = GET_SLOT(x, lme4_XSym),
	moff = GET_SLOT(x, lme4_offsetSym);
    int *xd = INTEGER(getAttrib(X, R_DimSymbol)), ione = 1;
    double one[] = {1,0};
    CHM_SP cVt = AS_CHM_SP(GET_SLOT(x, lme4_VtSym));
    CHM_DN cans = N_AS_CHM_DN(ans, xd[0], 1),
	cu = AS_CHM_DN(GET_SLOT(x, lme4_uvecSym));
    R_CheckStack();

    AZERO(ans, xd[0]);
    if (LENGTH(moff)) Memcpy(ans, REAL(moff), xd[0]);
				/* eta := eta + X \beta */
    F77_CALL(dgemv)("N", xd, xd + 1, one, REAL(X), xd,
		    REAL(GET_SLOT(x, lme4_fixefSym)), &ione, one, ans, &ione);
				/* eta := eta + V u */
    if (!M_cholmod_sdmult(cVt, 1 /* trans */, one, one, cu, cans, &c))
	error(_("cholmod_sdmult error returned"));
    return ans;
}    

/* Nonlinear mixed models */

/**
 * Create a dgCMatrix object in the form of the transpose of M
 *
 * @param Vt V' as a dgCMatrix object
 * @param sP pointer to the INTEGER object s
 *
 */
SEXP nlmer_create_Mt(SEXP Vt, SEXP sP)
{
    int n, p, s = asInteger(sP);
    CHM_SP ans;
    CHM_TR tMt = M_cholmod_sparse_to_triplet(AS_CHM_SP(Vt), &c);
    R_CheckStack();

    if (s <= 0) error(_("s must be > 0"));
    if (tMt->ncol % s)
	error(_("Number of columns in Vt, %d, is not a multiple of s = %d"),
	      tMt->ncol, s);
    n = tMt->ncol /= s;
    for (p = 0; p < tMt->nnz; p++) ((int*)(tMt->j))[p] %= tMt->ncol;
    ans = M_cholmod_triplet_to_sparse(tMt, tMt->nnz, &c);
    M_cholmod_free_triplet(&tMt, &c);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "", R_NilValue);
}

/**
 * Evaluate the nonlinear model, its gradient matrix and the residual
 * sum of squares
 *
 * @param x pointer to an nlmer object
 *
 */
SEXP nlmer_eval_model(SEXP x)
{
    SEXP gg, pnames = GET_SLOT(x, lme4_pnamesSym),
	rho = GET_SLOT(x, lme4_envSym), vv;
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *gdims;
    int i, n = dims[n_POS], s = dims[s_POS];
    double *y = REAL(GET_SLOT(x, lme4_ySym)), *mu, ans;
    double *Phi = Alloca(n * s, double);
    R_CheckStack();

    nglmer_eta(x, Phi);
    /* distribute the parameters in the environment */
    for (i = 0; i < dims[s_POS]; i++) {
	vv = findVarInFrame(rho, install(CHAR(STRING_ELT(pnames, i))));
	if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("Parameter %s in the environment must be a length %d numeric vector"),
		  CHAR(STRING_ELT(pnames, i)), n);
	Memcpy(REAL(vv), Phi + i * n, n);
    }
    vv = PROTECT(eval(GET_SLOT(x, lme4_modelSym), rho));
    if (!isReal(vv) || LENGTH(vv) != n)
	error(_("evaluated model is not a numeric vector of length %d"), n);
    gg = getAttrib(vv, lme4_gradientSym);
    if (!isReal(gg) || !isMatrix(gg))
	error(_("gradient attribute of evaluated model must be a numeric matrix"));
    gdims = INTEGER(getAttrib(gg, R_DimSymbol));
    if (gdims[0] != n ||gdims[1] != s)
	error(_("gradient matrix must be of size %d by %d"), n, s);
    SET_SLOT(x, lme4_muSym, vv);
    mu = REAL(vv);
    for (i = 0, ans = 0; i < n; i++) {
	double res = y[i] - mu[i];
	ans += res * res;
    }
    REAL(GET_SLOT(x, lme4_devianceSym))[disc_POS] = ans;
    UNPROTECT(1);
    return R_NilValue;
}

/**
 * Update the transpose of M from the current gradient and Vt
 *
 * @param x pointer to an nlmer object
 *
 * @return R_NilValue
 */
SEXP nlmer_update_Mt(SEXP x)
{
    SEXP Mt = GET_SLOT(x, lme4_MtSym), Vt = GET_SLOT(x, lme4_VtSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	*mi = INTEGER(GET_SLOT(Mt, lme4_iSym)),
	*mp = INTEGER(GET_SLOT(Mt, lme4_pSym)),
	*vi = INTEGER(GET_SLOT(Vt, lme4_iSym)),
	*vp = INTEGER(GET_SLOT(Vt, lme4_pSym)),
	cfll = c.final_ll, jv;
    int n = dims[n_POS], s = dims[s_POS];
    double *grad = REAL(getAttrib(GET_SLOT(x, lme4_muSym), lme4_gradientSym)),
	*mx = REAL(GET_SLOT(Mt, lme4_xSym)),
	*vx = REAL(GET_SLOT(Vt, lme4_xSym)), one[] = {1,0};
    CHM_SP cMt = AS_CHM_SP(Mt);
    CHM_FR L = L_SLOT(x);
    R_CheckStack();

    AZERO(mx, mp[n]);
    for (jv = 0; jv < n * s; jv++) {
	int iv, jm = jv % n, im;
	for (iv = vp[jv]; iv < vp[jv + 1]; iv++) {
	    for (im = mp[jm]; im < mp[jm + 1]; im++)
		if (mi[im] == vi[iv]) break;
	    if (im == mp[jm + 1])
		error(_("Structure of Mt incompatible with Vt, jv = %d, iv = %d"),
		      jv, iv);
	    mx[im] += grad[jv] * vx[iv];
	}
    }
    c.final_ll = L->is_ll;
    if (!M_cholmod_factorize_p(cMt, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    c.final_ll = cfll;
    REAL(GET_SLOT(x, lme4_devianceSym))[ldL2_POS] = chm_log_det2(L);
    return R_NilValue;
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an nlmer object
 *
 * @return An indicator of whether the iterations converged
 */
static int Nlmer_condMode(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *perm;
    int i, j, n = dims[n_POS], q = dims[q_POS];
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym)), *del,
	*mu, *u = REAL(GET_SLOT(x, lme4_uvecSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym)), crit = DOUBLE_XMAX,
	lpd_old = DOUBLE_XMAX, one[] = {1,0}, step, zero[] = {0,0};
    double *Mtres = Alloca(q, double),
	*res = Alloca(n, double),
	*rhs = Alloca(q, double),
	*uold = Alloca(q, double);
    CHM_FR L = L_SLOT(x);
    CHM_DN cMtres = N_AS_CHM_DN(Mtres, q, 1), cdel,
	cres = N_AS_CHM_DN(res, n, 1),
	crhs = N_AS_CHM_DN(rhs, q, 1);
    CHM_SP Mt = AS_CHM_SP(GET_SLOT(x, lme4_MtSym));
    R_CheckStack();

    mer_update_Vt(x);
    if (!(L->is_ll)) error(_("L must be LL', not LDL'"));
    if (q > n) error(_("q = %d > n = %d"), q, n);
    perm = (int*)(L->Perm);
    Memcpy(uold, u, q);
    nlmer_eval_model(x);
    for (i = 0; i < IRLS_MAXITER; i++) {
	mu = REAL(GET_SLOT(x, lme4_muSym)); /* nlmer_eval_model changes this */
	nlmer_update_Mt(x);	/* also updates L and dev[ldL2_POS] */
	for (j = 0; j < n; j++) res[j] = y[j] - mu[j]; /* residual then M'(y-mu) */
	if (!(M_cholmod_sdmult(Mt, 0 /* trans */, one, zero, cres, cMtres, &c)))
	    error(_("cholmod_sdmult returned error code"));
	crit = ((double)(n - q))/((lme4_sumsq(res, n) + lme4_sumsq(uold, q)) * q);
				/* apply permutation and subtract uold */
	for (j = 0; j < q; j++) rhs[j] = Mtres[perm[j]] - uold[j];
	if (!(cdel = M_cholmod_solve(CHOLMOD_L, L, crhs, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
	Memcpy(rhs, (double*)(cdel->x), q);
	M_cholmod_free_dense(&cdel, &c);
	crit *= lme4_sumsq(rhs, q);
	if (crit < IRLS_TOL) break;
	if (!(cdel = M_cholmod_solve(CHOLMOD_Lt, L, crhs, &c)))
	    error(_("cholmod_solve (CHOLMOD_L) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
	del = (double*)(cdel->x);
	for (step = 1; step > IRLS_SMIN; step /= 2) { /* step halving */
	    for (j = 0; j < q; j++) u[j] = uold[j] + step * del[j];
	    dev[bqd_POS] = lme4_sumsq(u, q);
	    nlmer_eval_model(x); /* also updates dev[disc_POS] */
	    dev[lpdisc_POS] = log(dev[disc_POS] + dev[bqd_POS]);
#ifdef DEBUG_NLMER
	    Rprintf("S %2d:%8.5f %8.5f %20.15g %20.15g %20.15g\n", i, crit, step,
		    dev[disc_POS], dev[bqd_POS], dev[lpdisc_POS]);
#endif
	    if (lpd_old > dev[lpdisc_POS]) break;
	}
	M_cholmod_free_dense(&cdel, &c);
	lpd_old = dev[lpdisc_POS];
	Memcpy(uold, u, q);
    }
    dev[ML_POS] = dev[REML_POS] =
	dev[ldL2_POS] + n * (1. + dev[lpdisc_POS] + log(2. * PI /((double)n)));
#ifdef DEBUG_NLMER
    Rprintf("C %2d:%20.15g %20.15g %20.15g %20.15f %20.15f\n", i, crit,
	    dev[disc_POS], dev[bqd_POS], dev[lpdisc_POS], dev[ML_POS]);
#endif
    M_cholmod_free_dense(&cdel, &c);
    return (crit > IRLS_TOL) ? 0 : i;
}

SEXP nlmer_condMode(SEXP x) {
    return ScalarInteger(Nlmer_condMode(x));
}

/* Generalized linear mixed models */
				/* utilities */
static const double LTHRESH = 30.;
static const double MLTHRESH = -30.;
static double MPTHRESH = 0;
static double PTHRESH = 0;
static const double INVEPS = 1/DOUBLE_EPS;

/** 
 * Evaluate x/(1 - x). An inline function is used so that x is
 * only evaluated once. 
 * 
 * @param x input in the range (0, 1)
 * 
 * @return x/(1 - x) 
 */
static R_INLINE double x_d_omx(double x) {
    if (x < 0 || x > 1)
	error(_("Value %d out of range (0, 1)"), x);
    return x/(1 - x);
}

/** 
 * Evaluate x/(1 + x). An inline function is used so that x is
 * evaluated once only.
 * 
 * @param x input
 * 
 * @return x/(1 + x) 
 */
static R_INLINE double x_d_opx(double x) {return x/(1 + x);}

static R_INLINE double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

/**
 * Evaluate the inverse link function at eta storing the result in mu
 *
 * @param x pointer to a glmer object
 */
SEXP glmer_linkinv(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[fTyp_POS];
    double *eta = REAL(findVarInFrame(rho, lme4_etaSym)),
	*mu = REAL(findVarInFrame(rho, lme4_muSym));

    switch(fltype) {
    case 1: 			/* binomial with logit link */
	for (i = 0; i < n; i++) {
	    double etai = eta[i], tmp;
	    tmp = (etai < MLTHRESH) ? DOUBLE_EPS :
		((etai > LTHRESH) ? INVEPS : exp(etai));
	    mu[i] = x_d_opx(tmp);
	}
	break;
    case 2:			/* binomial with probit link */
	if (!MPTHRESH) {
	    MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
	    PTHRESH = -MPTHRESH;
	}
	for (i = 0; i < n; i++) {
	    double etai = eta[i];
	    mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
		((etai > PTHRESH) ? 1 - DOUBLE_EPS :
		 pnorm5(etai, 0, 1, 1, 0));
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    double tmp = exp(eta[i]);
	    mu[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	}
 	break;
    default:
	error(_("General form of glmer_linkinv not yet written"));
     } 
    return R_NilValue;
}

/**
 * Evaluate the variance function for the link
 *
 * @param x pointer to a glmer object
 * @param var pointer to positions to hold computed values
 *
 * @return var
 */
static double *glmer_var(SEXP x, double *var)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[fTyp_POS];
    double *mu = REAL(findVarInFrame(rho, lme4_muSym));

    switch(fltype) {
    case 1: 			/* binomial family with logit or probit link */
    case 2:
	for (i = 0; i < n; i++) {
	    double mui = mu[i];
	    var[i] = mui * (1 - mui);
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    var[i] = mu[i];
	}
	break;
    default:
	error(_("General form of glmer_var not yet written"));
    }
    return var;
}

/**
 * Evaluate the derivative of mu wrt eta for the link
 *
 * @param x pointer to a glmer object
 * @param dmu_deta pointer to positions to hold computed values
 *
 * @return dmu_deta
 */
static double *glmer_dmu_deta(SEXP x, double *dmu_deta)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[fTyp_POS];
    double *eta = REAL(findVarInFrame(rho, lme4_etaSym));

    switch(fltype) {
    case 1: 			/* binomial with logit link */
	for (i = 0; i < n; i++) {
	    double etai = eta[i];
	    double opexp = 1 + exp(etai);
	    
	    dmu_deta[i] = (etai > LTHRESH || etai < MLTHRESH) ?
		DOUBLE_EPS : exp(etai)/(opexp * opexp);
	}
	break;
    case 2:			/* binomial with probit link */
	for (i = 0; i < n; i++) {
	    double tmp = dnorm4(eta[i], 0, 1, 0);
	    dmu_deta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    double tmp = exp(eta[i]);
	    dmu_deta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
	}
	break;
/*     default: { */
/* 	SEXP ans = PROTECT(eval_check(mu_eta, rho, REALSXP, n)); */
/* 	Memcpy(dmu_deta, REAL(ans), n); */
/* 	UNPROTECT(1); */
/*     } */
    }
    return dmu_deta;
}

/**
 * Evaluate the deviance residuals
 *
 * @param x pointer to a glmer object
 * @param dev_res pointer to an area to hold the result
 *
 * @return sum of the deviance residuals
 */
SEXP glmer_dev_resids(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, fltype = dims[fTyp_POS], n = dims[n_POS];
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	*dev_res = REAL(findVarInFrame(rho, lme4_devResidSym)),
	*mu = REAL(findVarInFrame(rho, lme4_muSym)),
	*wts = REAL(GET_SLOT(x, lme4_weightsSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym));

    switch(fltype) {
    case 1: 			/* binomial with logit or probit link */
    case 2:
	for (i = 0; i < n; i++) {
	    double mui = mu[i], yi = y[i];
	    
	    dev_res[i] = 2 * wts[i] *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
	}
	break;
    case 3:			/* Poisson with log link */
	for (i = 0; i < n; i++) {
	    double mui = mu[i], yi = y[i];
	    dev_res[i] = 2 * wts[i] * (y_log_y(yi, mui) - (yi - mui));
	}
	break;
    default:
	error(_("General form of glmer_dev_resids not yet written"));
    }
    for (i = 0, dev[disc_POS] = 0; i < n; i++)
	dev[disc_POS] += dev_res[i];
    return R_NilValue;
}

/**
 * Evaluate the linear predictor as model offset + X \beta + V u
 *
 * @param x pointer to a glmer object
 *
 * @return R_NilValue
 */
SEXP glmer_eta(SEXP x)
{
    nglmer_eta(x, REAL(findVarInFrame(GET_SLOT(x, lme4_envSym),
				      lme4_etaSym)));
    return R_NilValue;
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an glmer object
 *
 * @return An indicator of whether the iterations converged
 */
static int Glmer_condMode(SEXP x)
{
    SEXP env = GET_SLOT(x, lme4_envSym),
	moff = GET_SLOT(x, lme4_offsetSym);
    CHM_SP cVt = AS_CHM_SP(GET_SLOT(x, lme4_VtSym));
    int *vp = (int*)(cVt->p), cfll = c.final_ll, i, j, p;
    int n = cVt->ncol, nnz = vp[cVt->ncol], q = cVt->nrow;
    CHM_DN td1, td2;
    double crit = IRLS_TOL + 1, dn = (double)n,
	m_one[] = {-1, 0}, one[] = {1,0};
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	*eta = REAL(findVarInFrame(env, lme4_etaSym)),
	*mo = LENGTH(moff) ? REAL(moff) : (double*) NULL,
	*u = REAL(GET_SLOT(x, lme4_uvecSym)),
	*vx = (double*)(cVt->x),
	*w = REAL(GET_SLOT(x, lme4_weightsSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym));
    CHM_FR L = L_SLOT(x);
    double *dmu_deta = Alloca(n, double),
	*uold = Alloca(q, double), *var = Alloca(n, double),
	*wx = Alloca(nnz, double); /* weighted x slot of cVt */
    CHM_DN res = N_AS_CHM_DN(Alloca(n, double), n, 1),
	rhs = N_AS_CHM_DN(Alloca(q, double), q, 1);
    R_CheckStack();

    mer_update_Vt(x);
    Memcpy(uold, u, q); c.final_ll = L->is_ll;
    nglmer_eta(x, eta); glmer_linkinv(x);
    glmer_dmu_deta(x, dmu_deta); glmer_var(x, var);
    for (i = 0; i < IRLS_MAXITER; i++) {
	for (j = 0; j < n; j++) { /* calculate weighted Vt and residual */
	    double ww = sqrt(w[j] * dmu_deta[j] * dmu_deta[j]/var[j]);
	    ((double*)(res->x))[j] = y[j] - dmu_deta[j] - (mo ? mo[j] : 0);
	    for (p = vp[j]; p < vp[j + 1]; p++) wx[p] = vx[p] * ww;
	}
	cVt->x = (void*)wx;	/* Factor V'WV + I */
	if (!M_cholmod_factorize_p(cVt, one, (int*) NULL, (size_t) 0, L, &c)) { 
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	}
				/* evaluate V'[y-dmu_deta]-u */
	cVt->x = (void*)vx; Memcpy((double*)(rhs->x), u, q);
	if (!(j = M_cholmod_sdmult(cVt, 0/*trans*/, one, m_one, res, rhs, &c)))
	    error(_("cholmod_sdmult returned error code %d"), j);
				/* solve for the increment */
	td1 = M_cholmod_solve(CHOLMOD_L, L, rhs, &c);
	td2 = M_cholmod_solve(CHOLMOD_Lt, L, td1, &c);    
	M_cholmod_free_dense(&td1, &c);
	for (j = 0; j < q; j++) u[j] = uold[j] + ((double*)(td2->x))[j];
	M_cholmod_free_dense(&td2, &c);
	dev[bqd_POS] = lme4_sumsq(u, q);
	nglmer_eta(x, eta); glmer_linkinv(x); glmer_dev_resids(x);
	glmer_dmu_deta(x, dmu_deta); glmer_var(x, var);
 	crit = conv_crit(uold, u, q);
#ifdef DEBUG_NLMER
	Rprintf("%3d: %20.15g %20.15g", i, dev[disc_POS], dev[bqd_POS]);
	Rprintf(" %20.15g\n", dev[disc_POS] + dev[bqd_POS]);
#endif
	if (crit < IRLS_TOL) break;
    }
    dev[lpdisc_POS] = log(dev[disc_POS] + dev[bqd_POS]);
    dev[ldL2_POS] = chm_log_det2(L);
    dev[ML_POS] = dev[REML_POS] =
	dev[ldL2_POS] + dn * (1. + dev[lpdisc_POS] + log(2. * PI / dn));

    c.final_ll = cfll;
    return (crit > IRLS_TOL) ? 0 : i;
}

SEXP glmer_condMode(SEXP x) {
    return ScalarInteger(Glmer_condMode(x));
}


/* Functions common to all three forms of mixed models */

/**
 * Set the parameters in an lmer object and evaluate the deviance of
 * an lmm or the Laplace approximation to the deviance of a glmm.
 *
 * @param x pointer to a glmer object
 * @param xv vector of parameter values
 * @param mtype model type: 0 -> lmm, 1 -> nlmm, 2 -> glmm
 *
 * @return deviance
 */
static double
update_deviance(SEXP x, const double *xv, int mtype)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym));

    ST_SetPars(xv, ST);		/* common parameters */
    switch(mtype) {
    case 0: {			/* linear mixed model */
	lmer_update_dev(x);
	break;
    }
    case 1: {			/* nonlinear mixed model */
	Memcpy(REAL(GET_SLOT(x, lme4_fixefSym)),
	       xv + dims[np_POS], dims[p_POS]);
	Nlmer_condMode(x);
	break;
    }
    case 2: {			/* generalized linear mixed model */
	Memcpy(REAL(GET_SLOT(x, lme4_fixefSym)),
	       xv + dims[np_POS], dims[p_POS]);
	Glmer_condMode(x);
	break;
    }
    default:
	error(_("Unknown form of model for update_deviance"));
    }
    return dev[dims[isREML_POS] ? REML_POS : ML_POS];
}

/**
 * Optimize the deviance of an lmer object or the Laplace approximation
 * to the deviance of a nlmer or glmer object.
 *
 * @param x pointer to an lmer, glmer or nlmer object
 * @param verbp indicator of verbose output
 * @param mtypep model type 0 -> lmer, 1 -> nlmer, 2 -> glmer
 *
 * @return indicator of convergence from nlminb_iterate
 */
SEXP 
mer_optimize(SEXP x, SEXP verbp, SEXP mtypep)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, j,
	mtype = asInteger(mtypep), pos, verb = asLogical(verbp);
    int nf = dims[nf_POS], nv = dims[np_POS] + (mtype?dims[p_POS]:0);
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    double *g = (double*)NULL, *h = (double*)NULL, fx = R_PosInf;
    int *iv = Alloca(liv, int);
    double *b = Alloca(2 * nv, double), *d = Alloca(nv, double),
	*v = Alloca(lv, double), *xv = Alloca(nv, double);
    R_CheckStack();

    internal_ST_getPars(ST, xv);
    if (mtype) {
	SEXP fixef = GET_SLOT(x, lme4_fixefSym);
	Memcpy(xv + dims[np_POS], REAL(fixef), LENGTH(fixef));
    }
				/* initialize the state vectors v and iv */
    S_Rf_divset(OPT, iv, liv, lv, v);
    if (verb) iv[OUTLEV] = 1;
				/* set the bounds to plus/minus Infty  */
    for (i = 0; i < nv; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
				/* reset lower bounds on elements of S */
    for (i = 0, pos = 0; i < nf; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }
    S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    while (iv[0] == 1 || iv[0] == 2) {
	fx = update_deviance(x, xv, mtype); 
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, nv, v, xv);
    }
    i = iv[0];
    return ScalarInteger(i);
}

/* Functions for sampling from a Wishart distribution */
/* FIXME: Move these to the R sources */

/**
 * Simulate the Cholesky factor of a standardized Wishart variate with
 * dimension p and df degrees of freedom.
 *
 * @param df degrees of freedom
 * @param p dimension of the Wishart distribution
 * @param upper if 0 the result is lower triangular, otherwise upper
                triangular
 * @param ans array of size p * p to hold the result
 *
 * @return ans
 */
static double
*std_rWishart_factor(double df, int p, int upper, double ans[])
{
    int i, j, pp1 = p + 1;

    if (df < (double) p || p <= 0)
	error("inconsistent degrees of freedom and dimension");
    for (j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(df - (double) j));
	for (i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = norm_rand();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

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
lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal)
{
    SEXP ans;
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), i, j, k,
	n = asInteger(ns), psqr;
    double *scCp, *ansp, *tmp, df = asReal(dfp), one = 1, zero = 0;

    if (!isMatrix(scal) || !isReal(scal) || dims[0] != dims[1])
	error("scal must be a square, real matrix");
    if (n <= 0) n = 1;
    psqr = dims[0] * dims[0];
    tmp = Alloca(psqr, double);
    scCp = Alloca(psqr, double);
    R_CheckStack();

    Memcpy(scCp, REAL(scal), psqr);
    AZERO(tmp, psqr);
    F77_CALL(dpotrf)("U", &(dims[0]), scCp, &(dims[0]), &j);
    if (j)
	error("scal matrix is not positive-definite");
    PROTECT(ans = alloc3DArray(REALSXP, dims[0], dims[0], n));
    ansp = REAL(ans);
    GetRNGstate();
    for (j = 0; j < n; j++) {
	double *ansj = ansp + j * psqr;
	std_rWishart_factor(df, dims[0], 1, tmp);
	F77_CALL(dtrmm)("R", "U", "N", "N", dims, dims,
			&one, scCp, dims, tmp, dims);
	F77_CALL(dsyrk)("U", "T", &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));

	for (i = 1; i < dims[0]; i++)
	    for (k = 0; k < i; k++)
		ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
    }

    PutRNGstate();
    UNPROTECT(1);
    return ans;
}

/* MCMC sampling functions */


/**
 * Update the fixed effects and the random effects in an MCMC sample
 * from an lmer model.
 *
 * @param x an lmer model with mer_update_b
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
static void lmer_MCMC_betab(SEXP x, double sigma)
{
    SEXP fixefP = GET_SLOT(x, lme4_fixefSym),
	uvecP = GET_SLOT(x, lme4_uvecSym);
    int ione = 1, j, p = LENGTH(fixefP), q = LENGTH(uvecP);
    int pp1 = p + 1;
    double *fix = REAL(fixefP), *uv = REAL(uvecP), mone = -1, one = 1;
    CHM_FR L = L_SLOT(x);
    double *del1 = Alloca(q, double), *del2 = Alloca(p, double);
    CHM_DN sol, rhs = N_AS_CHM_DN(del1, q, 1);
    R_CheckStack();

    lmer_update_effects(x);
				/* Update the fixed effects */
    for (j = 0; j < p; j++) del2[j] = sigma * norm_rand();
    F77_CALL(dtrsv)("U", "N", "N", &p, REAL(GET_SLOT(x, lme4_RXySym)),
		    &pp1, del2, &ione);
    for (j = 0; j < p; j++) fix[j] += del2[j];
				/* Update the orthogonal random effects */
    for (j = 0; j < q; j++) del1[j] = sigma * norm_rand();
    F77_CALL(dgemv)("N", &q, &p, &mone, REAL(GET_SLOT(x, lme4_RVXySym)),
		    &q, del2, &ione, &one, del1, &ione);
    sol = M_cholmod_solve(CHOLMOD_Lt, L, rhs, &c);
    for (j = 0; j < q; j++) uv[j] += ((double*)(sol->x))[j];
    M_cholmod_free_dense(&sol, &c);
				/* Update the random effects */
    TSPp_SEXP_mult(GET_SLOT(x, lme4_ranefSym), GET_SLOT(x, lme4_STSym),
		   INTEGER(GET_SLOT(x, lme4_GpSym)), uvecP, (int*)(L->Perm));
}

static void
safe_pd_matrix(double x[], const char uplo[], int n, double thresh)
{
    int info, lwork = 3 * n, nm1 = n - 1;
    double *work = Alloca(3 * n, double),
	*w = Alloca(n, double),
	*xcp = Memcpy(Alloca(n * n, double), x, n * n);

    F77_CALL(dsyev)("N", uplo, &n, xcp, &n, w, work, &lwork, &info);
    if (info) error(_("dsyev returned %d"), info);
    if (w[nm1] <= 0) error(_("no positive eigenvalues!"));
    if ((w[0]/w[nm1]) < thresh) {
	int i, np1 = n + 1;
	double incr = w[nm1] * thresh;
	for (i = 0; i < n; i++) x[i * np1] += incr;
    }
}

static void mer_MCMC_ST(SEXP x, double sigma, double *vals, int trans)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int i, info, j, k, ncmax, nf = LENGTH(ST),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym));
    double *rr = REAL(GET_SLOT(x, lme4_ranefSym)),
	*scal, *var, *wfac, one = 1, zero = 0;
    double **st = Alloca(nf, double*);
    int *nc = Alloca(nf, int), *nlev = Alloca(nf, int);
    R_CheckStack();

    ncmax = ST_nc_nlev(ST, Gp, st, nc, nlev);
    scal = Alloca(ncmax * ncmax, double);
    wfac = Alloca(ncmax * ncmax, double);
    var = Alloca(ncmax * ncmax, double);
    R_CheckStack();
    for (i = 0; i < nf; i++) {
	int nci = nc[i], nl = nlev[i];
	int ncisq = nci * nci, ncip1 = nci + 1;
	
	AZERO(scal, ncisq);
	F77_CALL(dsyrk)("U", "T", &nci, &nl, &one, rr + Gp[i], &nl,
			&zero, scal, &nci); 
	if (nci == 1) *scal = sqrt(*scal);
	else {
	    safe_pd_matrix(scal, "L", nci, 1e-7);
	    F77_CALL(dpotrf)("U", &nci, scal, &nci, &info);
	    if (info)
		error(_("scal matrix is not pd after safe_pd_matrix!"));
	}
				/* generate random factor from std Wishart */
	AZERO(wfac, ncisq);
	std_rWishart_factor((double) (nl - nci + 1), nci, 1, wfac);
				/* form and store variance-covariance matrix */
	F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &nci,
			&one, wfac, &nci, scal, &nci);
	F77_CALL(dsyrk)("U", "T", &nci, &nci, &one, scal, &nci,
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
	if (nci == 1) {
	    *st[i] = sqrt(*var)/sigma;
	} else {
	    double *sti = st[i];

	    for (j = 0; j < nci; j++) /* copy var to st[i] transposing */
		for (k = j; k < nci; k++)
		    sti[k + j * nci] = var[j + k * nci];
				/* factor the variance-covariance matrix */
	    F77_CALL(dpotrf)("L", &nci, sti, &nci, &info);
	    if (info)
		error(_("new var-cov matrix is not positive definite!"));
	    for (j = 0; j < nci; j++) { /* convert to LDL' then TSST' */
		for (k = j + 1; k < nci; k++)
		    sti[k + j * nci] /= sti[j * ncip1];
		sti[j * ncip1] /= sigma;
	    }
	}
    }
}

/**
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 *
 * @param x pointer to an lmer object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP lmer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbosep, SEXP deviancep)
{
    SEXP Pars = PROTECT(ST_getPars(x)), ans;
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int dPOS = dims[REML_POS] ? REML_POS : ML_POS, i, j,
	n = dims[n_POS], nsamp = asInteger(nsampp),
	p = dims[p_POS], q = dims[q_POS],
	saveb = asLogical(savebp), trans = asLogical(transp),
	verbose = asLogical(verbosep), dev = asLogical(deviancep);
    double *ansp, *deviance = REAL(GET_SLOT(x, lme4_devianceSym)),
	*fixef = REAL(GET_SLOT(x, lme4_fixefSym)),
	*ranef = (double*)NULL, df = n - (dims[REML_POS] ? p : 0);
    int nrbase = p + 1 + dims[np_POS]; /* rows always included */
    int nrtot = nrbase + dev + (saveb ? q : 0);

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    GetRNGstate();
    if (verbose) Rprintf("%12s %14s\n", "sigma", "deviance");
    if (saveb) ranef = REAL(GET_SLOT(x, lme4_ranefSym));

    for (i = 0; i < nsamp; i++) {
	double *col = ansp + i * nrtot, sigma;
				/* simulate and store new value of sigma */
	sigma = exp(deviance[lpdisc_POS]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
	lmer_update_effects(x);	/* Conditional estimates of beta, u and b */
	lmer_MCMC_betab(x, sigma); /* Updated beta, u and b */
	for (j = 0; j < p; j++) col[j] = fixef[j]; /* Store beta */
	if (saveb)				   /* store b */
	    for (j = 0; j < q; j++) col[nrbase + dev + j] = ranef[j];
	mer_MCMC_ST(x, sigma, col + p + 1, trans);
	lmer_update_dev(x);	/* Refactor and evaluate deviance */
	/* FIXME: This deviance should be for sigma, beta, b, ST */
	if (dev) col[nrbase] = deviance[dPOS]; 
	if (verbose) Rprintf("%12.6g %14.8g\n", sigma, deviance[dPOS]);
    }
    PutRNGstate();
				/* Restore pars, refactor, etc. */
    ST_setPars(x, Pars);
    lmer_update_dev(x);
    lmer_update_effects(x);
    UNPROTECT(2);
    return ans;
}

/* Functions for the pedigree class */

/**
 * Create the left Cholesky factor of the numerator relationship
 * matrix from a pedigree.
 *
 * @param x a pedigree object
 * @param ans T stored as a non-unit, lower dtCMatrix
 *
 * @return ans with elements modified to incorporate D
 */
SEXP pedigree_chol(SEXP x, SEXP ans)
{
    SEXP Sire = GET_SLOT(x, install("sire"));
    int *ai = INTEGER(GET_SLOT(ans, lme4_iSym)),
	*ap = INTEGER(GET_SLOT(ans, lme4_pSym)),
	*dam = INTEGER(GET_SLOT(x, install("dam"))),
	*sire = INTEGER(Sire), 
	i, j, n = LENGTH(Sire);
    double *ax = REAL(GET_SLOT(ans, lme4_xSym)), *F, Di, tmp;

    setAttrib(ans, install("F"), allocVector(REALSXP, n));
    F = REAL(getAttrib(ans, install("F")));
    for (i = 0; i < n; i++) {
	int k, p = sire[i] - 1, q = dam[i] - 1;
	if (sire[i] == NA_INTEGER) {
	    F[i] = 0;
	    Di = (dam[i] == NA_INTEGER) ? 1 : sqrt(0.75 - 0.25 * F[q]);
	} else {
	    if (dam[i] == NA_INTEGER) { /* sire only */
		F[i] = 0;
		Di = sqrt(0.75 - 0.25 * F[p]);
	    } else {		/* both parents in pedigree */
		Di = sqrt(0.5 - 0.25 * (F[p] + F[q]));
		F[i] = NA_REAL;
		if ((ap[i + 1] - ap[i]) > 1) {	  /* skip if no progeny */
		    if (p > q) {j = p; p = q; q = j;} /* ensure p <= q */
		    for (j = 0, F[i] = 0; j <= p; j++) {
			for (k = ap[j], tmp = 0;
			     k < ap[j + 1] && ai[k] <= q; k++) {
			    int ii = ai[k];
			    if (ii == p) tmp = ax[k];
			    if (ii == q) F[i] += tmp * ax[k]/2;
			}
		    }
		}
	    }
	}
	for (j = ap[i]; j < ap[i + 1]; j++) ax[j] *= Di;
    }
    return ans;
}

/* NOTE: This function requires that missing parents be coded as zero */
/**
 * Create the inbreeding coefficients according to the algorithm given
 * in "Comparison of four direct algorithms for computing inbreeding
 * coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
 * Science Journal (2005) 76, 401--406.  This function is a modified
 * version of the code published in an appendix to that paper.
 *
 * @param x a pedigree object
 *
 * @return a list of the inbreeding coefficients
 */

SEXP pedigree_inbreeding(SEXP x)
{
    SEXP ans, sp = GET_SLOT(x, install("sire"));
    int i, j, t, n = LENGTH(sp), S, D;
    int *SI, *MI,		/* start and minor */
	*sire = INTEGER(sp),
	*dam = INTEGER(GET_SLOT(x, install("dam")));
    double *F = Alloca(n + 1, double), /* inbreeding coefficients */
      *L = Alloca(n + 1, double), *B = Alloca(n + 1, double);
    int *Anc = Alloca(n + 1, int),	/* ancestor */
	*LAP = Alloca(n + 1, int); 	/* longest ancestoral path */
    R_CheckStack();
    
    F[0] =-1; LAP[0] =-1; /* set F and lap for unknown parents */
    for(i = 1, t = -1; i <= n; i++) { 	/* evaluate LAP and its maximum */
	S = sire[i]; D = dam[i]; /* parents of animal i */
	LAP[i] = ((LAP[S] < LAP[D]) ? LAP[D] : LAP[S]) + 1;
	if (LAP[i] > t) t = LAP[i];
    }
    SI = Alloca(t + 1, int);
    MI = Alloca(t + 1, int);
    for(i = 0; i <= t ; ++i) SI[i] = MI[i] = 0; /* initialize start and minor */
    for(i = 1; i <= n; i++) { 	/* evaluate F */
	S = sire[i]; D = dam[i]; /* parents of animal i */
	B[i] = 0.5 - 0.25 * (F[S] + F[D]); 
				/* adjust start and minor */
	for (j = 0; j < LAP[i]; j++) {++SI[j]; ++MI[j];} 
	if (S == 0 || D == 0) { /* both parents unknown */
	    F[i] = L[i] = 0; continue;
	}
	if(S == sire[i-1] && D == dam[i-1]) { /* full-sib with last animal */
	    F[i] = F[i-1]; L[i] = L[i-1]; continue;
	}
    
	F[i] = -1; L[i] = 1; 
	t = LAP[i]; /* largest lap group number in the animal's pedigree */
	Anc[MI[t]++] = i; /* initialize Anc and increment MI[t] */
	while(t > -1) { /* from the largest lap group to zero */
	    j = Anc[--MI[t]]; /* next ancestor */
	    S = sire[j]; D = dam[j]; /* parents of the ancestor */
	    if (S) {
		if (!L[S]) Anc[MI[LAP[S]]++] = S; 
				/* add sire in its lap group in Anc
				 * array if it is not added yet and
				 * increment the minor index for the group */ 
		L[S] += 0.5 * L[j]; /* contribution to sire */
	    }
	    if (D) {
		if (!L[D]) Anc[MI[LAP[D]]++] = D;
		L[D] += 0.5 * L[j]; /* contribution to dam */
	    }
	    F[i] += L[j] * L[j] * B[j];
	    L[j] = 0; /*clear L[j] for the evaluation of the next animal */
	    if (MI[t] == SI[t]) --t; /* move to the next lap group when
				      * all ancestors in group t have been
				      * evaluated */
	} 
    }
    ans = PROTECT(allocVector(REALSXP, n));
    Memcpy(REAL(ans), F + 1, n);
    UNPROTECT(1);
    return ans;
}


#if 0

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


/**
 * Update the working residual as y - mu + M u
 *
 * @param x pointer to an nlmer object
 * @param w array to hold the updated working residual
 *
 * @return w pointer
 */
static double *Nlmer_update_wrkres(SEXP x, double *w)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i;
    CHM_SP Mt = AS_CHM_SP(GET_SLOT(x, lme4_MtSym));
    CHM_DN cwrk = N_AS_CHM_DN(w, dims[n_POS], 1),
	cu = AS_CHM_DN(GET_SLOT(x, lme4_uvecSym));
    double *mu = REAL(GET_SLOT(x, lme4_muSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym)), one[] = {1,0};

    R_CheckStack();
    Memcpy(w, y, dims[n_POS]);
    if (!(i = M_cholmod_sdmult(Mt, 1 /* trans */, one, one, cu, cwrk, &c)))
	error(_("cholmod_sdmult returned error code %d"), i);
    for (i = 0; i < dims[n_POS]; i++) w[i] -= mu[i];
    return w;
}

/**
 * Externally callable function to return the working residual
 *
 * @param x pointer to an nlmer object
 *
 * @return working residual
 */
SEXP nlmer_update_wrkres(SEXP x)
{
    SEXP ans = PROTECT(allocVector(REALSXP, LENGTH(GET_SLOT(x, lme4_ySym))));

    Nlmer_update_wrkres(x, REAL(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Update the relative variance-covariance matrices by sampling from a
 * Wishart distribution with scale factor determined by the current
 * sample of random effects.
 *
 * @param sigma current value of sigma
 * @param Gp vector of group pointers
 * @param unew current sample from the orthogonal random effects
 * @param ST list of factorizations to update
 */
static void
internal_ST_update(double sigma, int trans, const int *Gp,
		   const double *unew, SEXP ST, double *dest)
{
    int i, j, k, info, nf = LENGTH(ST), pos = 0;
    double one = 1, zero = 0, sigsq = sigma * sigma;
    double sigsqinv = 1/sigsq;

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi), sd;
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int nlev = (Gp[i + 1] - Gp[i])/nci;

	if (nci == 1) {		/* fast update for scalar */
	    sd = sqrt(lme4_sumsq(unew + Gp[i], nlev)/rchisq((double)nlev)) * st[0];
	    st[0] = sd/sigma;
	    dest[pos++] = (trans ? 2 * log(sd) : sd * sd);
	} else {
	    int ncip1 = nci + 1, ncisqr = nci * nci;
	    double *scal = Calloc(ncisqr, double), /* factor of scale matrix */
		*wfac = Calloc(ncisqr, double); /* factor of Wishart variate */
				/* create the scale matrix (close to identity) */
	    AZERO(scal, ncisqr);
	    F77_CALL(dsyrk)("L", "N", &nci, &nlev, &sigsqinv, unew + Gp[i], &nci,
			    &zero, scal, &nci);
	    F77_CALL(dpotrf)("L", &nci, scal, &nci, &info);
	    if (info) error(_("Crossprod of b*[[%d]] not positive definite"), i + 1);
	    for (i = 0; i < nci; i++) /* premultiply by S */
		for (j = 0; j <= i; j++) scal[i + j * nci] *= st[i * ncip1];
	    F77_CALL(dtrmm)("L", "L", "N", "U", /* premultiply by T */
			    &nci, &nci, &one, st, &nci, scal, &nci);
				/* generate (lower) random factor from std Wishart */
	    std_rWishart_factor((double)(nlev - nci + 1), nci, 0, wfac);
	    F77_CALL(dtrsm)("L", "L", "N", "N",/* form a scaled variance factor */
			    &nci, &nci, &one, wfac, &nci, scal, &nci);
	    Memcpy(st, scal, ncisqr);
	    for (j = 0; j < nci; j++) { /* form the ST representation */
		for (i = j + 1; i < nci ; i++)
		    st[i + j * nci] /= st[j * ncip1];
	    }
				/* Overwrite wfac with variance-covariance */
	    F77_CALL(dsyrk)("L", "N", &nci, &nci, &sigsq, scal, &nci,
			    &zero, wfac, &nci);
	    if (trans) {
		for (j = 0; j < nci; j++) {
		    double vj = wfac[j * ncip1]; /* jth variance */
		    double sdj = sqrt(vj);

		    for (k = 0; k < j; k++) /* jth row in lower tri */
			wfac[k * nci + j] = atanh(wfac[k * nci + j]/sdj);
		    for (k = j + 1; k < nci; k++)
			wfac[j * nci + k] /= sdj; /* jth col in lower tri */
		    wfac[j * ncip1] = log(vj);
		}
	    }
	    for (j = 0; j < nci; j++) dest[pos++] = wfac[j * ncip1];
	    for (j = 1; j < nci; j++) {
		for (k = 0; k < j; k++)
		    dest[pos++] = wfac[k * nci + j];
	    }
	    Free(scal); Free(wfac);
	}
    }
}

#endif
