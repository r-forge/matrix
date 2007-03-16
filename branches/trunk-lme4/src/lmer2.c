#include "lmer2.h"

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

/* Functions for the mer2 representation */

				/* positions in the deviance vector */
enum devP {ML_POS=0, REML_POS, ldZ_POS, ldX_POS, lr2_POS};
			/* {"ML", "REML", "ldZ", "ldX", "lr2", ""} */
				/* positions in the dims vector */
enum dimP {nf_POS=0, n_POS, p_POS, q_POS, isREML_POS, isGLMM_POS, isNest_POS};
	      /* {"nf", "n", "p", "q", "isREML", "isGLMM", "isNested"} */

#define isREML(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isREML_POS]
#define isGLMM(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isGLMM_POS]
#define isNested(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isNest_POS]

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
static int check_nesting(int nf, const int nc[], const int Gp[], const int p[])
{
    int **cnz = Calloc(nf, int*), ans = 1, i, j, k, nct;

    for (i = 0, nct = 0; i < nf; i++) { /* total number of columns */
	nct += nc[i];
	cnz[i] = Calloc(nc[i], int);
    }
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
    for (i = 0, nct = 0; i < nf; i++) Free(cnz[i]);
    Free(cnz);
    return ans;
}

/**
 * Multiply a vector by the virtual T and S matrices represented by ST
 *
 * @param Gp vector of group pointers
 * @param ST compounded matrices S* and T*
 * @param b vector of random effects to be transformed from b* to b
 *
 * @return b after transformation
 */
static double *TS_mult(const int *Gp, SEXP ST, double *b)
{
    int i, ione = 1, nf = LENGTH(ST);

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncip1 = nci + 1,
	    nlev = (Gp[i + 1] - Gp[i])/nci;

	for (j = 0; j < nlev; j++) {
	    int base = Gp[i] + j * nci;
	    for (k = 0; k < nci; k++) /* multiply by S_i */
		b[base + k] *= st[k * ncip1];
	    if (nci > 1) 	/* multiply by T_i */
		F77_CALL(dtrmv)("L", "N", "U", &nci, st, &nci,
				b + base, &ione);
	}
    }
    return b;
}

/**
 * Evaluate the effects in an mer2 representation
 *
 * @param L factorization
 *
 * @return cholmod_dense object with \hat{b*} in the first q positions
 * and \hat{\beta} in the next p positions.
 */
static cholmod_dense
*internal_mer2_effects(cholmod_factor *L)
{
    int i, nt = (L->n);
    cholmod_dense *X,
	*B = M_cholmod_allocate_dense((size_t)nt, 1, (size_t)nt, CHOLMOD_REAL, &c);
    double *bx = (double*)(B->x);
    
    for (i = 0; i < nt; i++) bx[i] = 0;
    if (L->is_super) {
	int ns = (L->nsuper);
	int nr = ((int *)(L->pi))[ns] - ((int *)(L->pi))[ns - 1],
	    nc = ((int *)(L->super))[ns] - ((int *)(L->super))[ns - 1];
	double *x = (double *)(L->x) + ((int *)(L->px))[ns - 1];

	bx[nt - 1] = x[(nc - 1) * (nr + 1)];
    } else {
	bx[nt - 1] = (L->is_ll) ? ((double*)(L->x))[((int*)(L->p))[nt - 1]] : 1;
    }
    if (!(X = M_cholmod_solve(CHOLMOD_Lt, L, B, &c)))
	error(_("cholmod_solve (CHOLMOD_Lt) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    M_cholmod_free_dense(&B, &c);
    if (!(B = M_cholmod_solve(CHOLMOD_Pt, L, X, &c)))
	error(_("cholmod_solve (CHOLMOD_Pt) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    M_cholmod_free_dense(&X, &c);
    return B;
}

/**
 * Evaluate the logarithm of the square of the determinant of selected
 * sections of a sparse Cholesky factor.
 *
 * @param ans vector of doubles of sufficient length to hold the result
 * @param nans number of values to calculate
 * @param c vector of length nans+1 containing the cut points
 * @param F factorization
 *
 * @return ans
 */
static double*
chm_log_abs_det2(double *ans, int nans, const int *c, const cholmod_factor *F)
{
    int i, ii  = 0, jj = 0;
    for (i = 0; i < nans; i++) ans[i] = 0;
    if (F->is_super) {
	for (i = 0; i < F->nsuper; i++) {
	    int j, nrp1 = 1 + ((int *)(F->pi))[i + 1] - ((int *)(F->pi))[i],
		nc = ((int *)(F->super))[i + 1] - ((int *)(F->super))[i];
	    double *x = (double *)(F->x) + ((int *)(F->px))[i];

	    for (j = 0; j < nc; j++) {
		int col = j + jj;
		if (col < c[ii]) continue;
		while (col >= c[ii + 1] && ++ii < nans) {};
		if (ii >= nans) break;
		ans[ii] += 2 * log(fabs(x[j * nrp1]));
	    }
	    jj += nc;
	}
    } else {
	int *fi = (int*)(F->i), *fp = (int*)(F->p), j, k;
	double *fx = (double *)(F->x);
	
	for (j = 0; ii < nans && j < F->n; j++) {
	    if (j < c[ii]) continue;
	    for (k = fp[j]; fi[k] != j && k < fp[j + 1]; k++) {};
	    if (fi[k] != j) break; /* what happened to the diagonal element? */
	    while (j >= c[ii + 1] && ++ii < nans) {};
	    if (ii >= nans) break;
	    ans[ii] += log(fx[k] * ((F->is_ll) ? fx[k] : 1.));
	}
    }
    return ans;
}

/**
 * Evaluate the elements of the deviance slot given a factorization of
 * A* and the dimensions vector.
 *
 * @param d pointer to the contents of the slot
 * @param dims dimensions
 * @param L factor of the current A*
 *
 * @return d
 */
static double*
internal_deviance(double *d, const int *dims, const cholmod_factor *L)
{
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
    int c[] = {0,  q, p + q, p + q + 1};
    double dn = (double) n, dnmp = (double)(n - p);
    
    chm_log_abs_det2(d + ldZ_POS, lr2_POS + 1 - ldZ_POS, c, L);
    d[0] = d[2] + dn * (1. + d[4] + log(2. * PI / dn));
    d[1] = d[2] + d[3] + dnmp * (1. + d[4] + log(2. * PI / dnmp));
    return d;
}

/**
 * Update A to A* and evaluate its numeric factorization in L.
 *
 * @param deviance Hold the result
 * @param dims dimensions {
 * @param nc length nf vector of number of random effects per factor
 * @param Gp length nf+3 vector of group pointers for the rows of A
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param A symmetric matrix of size Gp[nf+2]
 * @param F factorization to be modified
 *
 */
static void
internal_update_L(double *deviance, int *dims, const int *Gp,
		  SEXP ST, cholmod_sparse *A, cholmod_factor *L)
{
    cholmod_sparse *Ac = M_cholmod_copy_sparse(A, &c);
    int *ai = (int *)(Ac->i), *ap = (int *)(Ac->p), nf = *dims,
	i, j, ione = 1;
    double *ax = (double*)(Ac->x) , one[] = {1, 0};
    

    if ((!Ac->sorted) || Ac->stype <= 0) {
	M_cholmod_free_sparse(&Ac, &c);
	error(_("A must be a sorted cholmod_sparse object with stype > 0"));
    }

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int base = Gp[i], k, kk;
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - Gp[i])/nci;

	/* if nci == 1 then T == I and the next section is skipped */
	if (nci > 1) {	/* evaluate ith section of SAST */
	    int maxrows = -1;
	    double *db = Calloc(nci * nci, double), /* diagonal blcok */
		*wrk = (double *) NULL;
	    
/* FIXME: calculate and store maxrows in mer2_create */
	    if (nci > 2) {	/* calculate scratch storage size */
		for (j = 0; j < nlev; j++) {
		    int cj = base + j * nci; /* first column in this group */
		    int nnzm1 = ap[cj + 1] - ap[cj] - 1;
		    if (nnzm1 > maxrows) maxrows = nnzm1;
		}
		wrk = Calloc(maxrows * nci, double);
	    }
	    
	    for (j = 0; j < nlev; j++) {
		int cj = base + j * nci; /* first column in this block */
		int nnz = ap[cj + 1] - ap[cj]; /* nonzeros in column cj */
		int nnzm1 = nnz - 1;
		
		if (nnzm1) { /* elements above diagonal block to update */
		    if (nci == 2)
			F77_CALL(dtrmm)("R", "L", "N", "U", &nnzm1,
					&nci, one, st, &nci,
					ax + ap[cj], &nnz);
		    else {	
			for (k = 0; k < nci; k++) /* copy columns to wrk */
			    Memcpy(wrk + k * maxrows, ax + ap[cj + k],
				   nnzm1);
			F77_CALL(dtrmm)("R", "L", "N", "U", &nnzm1,
					&nci, one, st, &nci, wrk,
					&maxrows);
			for (k = 0; k < nci; k++) /* copy results back */
			    Memcpy(ax + ap[cj + k], wrk + k * maxrows,
				   nnzm1);
		    }
		    /* evaluate T'A for rows and columns in this block */
		    for (k = 0; k < nci; k++) {
			for (kk = 0; kk < nnzm1;) {
			    int ind = ap[cj + k] + kk;
			    if (Gp[i] <= ai[ind]) {
				F77_CALL(dtrmv)("L", "T", "U", &nci,
						st, &nci, ax + ind,
						&ione);
				kk += nci;
			    } else kk++;
			}
		    }
		}
				/* update the diagonal block */
		for (k = 0; k < nci; k++) /* copy upper triangle */
		    Memcpy(db + k * nci, ax + ap[cj + k] + nnzm1, k + 1);
		for (k = 1; k < nci; k++) /* symmetrize */
		    for (kk = 0; kk < k; kk++)
			db[k + kk * nci] = db[kk + k * nci];
		F77_CALL(dtrmm)("L", "L", "T", "U", &nci, &nci, one,
				st, &nci, db, &nci);
		F77_CALL(dtrmm)("R", "L", "N", "U", &nci, &nci, one,
				st, &nci, db, &nci);
		for (k = 0; k < nci; k++) /* restore updated upper triangle */
		    Memcpy(ax + ap[cj + k] + nnzm1, db + k * nci, k + 1);
	    }
	    if (nci > 2) Free(wrk);
				/* evaluate T'AT for all blocks to the right */
	    for (j = Gp[i+1]; j < Gp[nf + 2]; j++) {
		for (k = ap[j]; k < ap[j + 1]; ) {
		    if (ai[k] >= Gp[i + 1]) break;
		    if (ai[k] < Gp[i]) {
			k++;
			continue;
		    }
		    F77_CALL(dtrmv)("L", "T", "U", &nci, st, &nci,
				    ax + k, &ione);
		    k += nci;
		}
	    }
	    Free(db);
	}
				/* Multiply by S from left. */
	for (j = Gp[i]; j < A->ncol; j++)
	    for (k = ap[j]; k < ap[j+1]; k++) {
		int row = ai[k];
		if (row < Gp[i]) continue;
		if (Gp[i + 1] <= row) break;
		ax[k] *= st[((row - Gp[i]) % nci) * ncip1];
	    }
				/* Multiply by S from right */
	for (j = Gp[i]; j < Gp[i + 1]; j += nci) {
	    for (k = 0; k < nci; k++)
		for (kk = ap[j + k]; kk < ap[j + k + 1]; kk++)
		    ax[kk] *= st[k * ncip1];
	}
				/* Increment diagonal */
	for (j = Gp[i]; j < Gp[i + 1]; j++) {
	    k = ap[j + 1] - 1;
	    if (ai[k] != j) error(_("Logic error"));
	    ax[k]++;
	}
    }
				/* force LL' decomp and sorted row indices */
    i = c.final_ll; c.final_ll = TRUE;
    j = c.final_monotonic; c.final_monotonic = TRUE;
    if (!M_cholmod_factorize(Ac, L, &c)) { 
	error(_("cholmod_factorize failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    }	
				/* restore previous settings */
    c.final_ll = i; c.final_monotonic = j;
    internal_deviance(deviance, dims, L);
    M_cholmod_free_sparse(&Ac, &c);
}

static void
internal_mer2_initial(SEXP ST, int *Gp, cholmod_sparse *A)
{
    int *ai = (int*)(A->i), *ap = (int*)(A->p), i, nf = LENGTH(ST);
    double *ax = (double*)(A->x);
    
    if (!(A->sorted) || (A->stype <= 0))
	error(_("A should be upper triangular and sorted"));
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int bb = Gp[i], j, k;
	int ncip1 = nci + 1, nlev = (Gp[i + 1] - bb)/nci;
	
	AZERO(st, nci * nci);
	for (j = 0; j < nlev; j++) {
	    int base = bb + j * nci;
	    for (k = 0; k < nci; k++) {
		int cl = ap[base + k + 1] - 1; /* last element in the column */
		if (ai[cl] != (base + k))
		    error(_("Weird structure in A.  Expecting row %d, got %d"),
			  base + k, ai[cl]);
		st[k * ncip1] += ax[cl];
	    }
	}
	for (k = 0; k < nci; k++)
	    st[k * ncip1] = sqrt(((double) nlev)/(0.375*st[k * ncip1]));
    }
}

static R_INLINE void
make_cholmod_sparse_sorted(cholmod_sparse *A)
{
    if(!A->sorted) {
	int i = M_cholmod_sort(A, &c); 
	if(!i)
	    error(_("cholmod_sort returned error code %d"),i);
    }
}

static void
internal_update_A(cholmod_sparse *ZXyt, SEXP wtP, SEXP offP,
		  cholmod_sparse *A)
{
    cholmod_sparse *ts1 = M_cholmod_copy_sparse(ZXyt, &c), *ts2;
    int *ai, *ap, *zi = (int*)(ts1->i), *zp = (int*)(ts1->p),
	i, j, m = ts1->nrow, n = ts1->ncol,
	ol = LENGTH(offP), wl = LENGTH(wtP);
    double *ax, *zx = (double*)(ts1->x), *off = REAL(offP), *wts = REAL(wtP);

    make_cholmod_sparse_sorted(ts1);
    if (ol) {		/* skip if length 0 */
	if (ol != n) {
	    M_cholmod_free_sparse(&ts1, &c);
	    error(_("Length of offset is %d, should be %d"), ol, n);
	}
	for (j = 0; j < n; j++) { /* iterate over columns */
	    int ind = zp[j + 1] - 1;
	    if (zi[ind] != (m - 1)) {
		M_cholmod_free_sparse(&ts1, &c);
		error(_("missing y position in ZXyt at column %d"), j+1);
	    }
	    zx[ind] -= off[j];
	}
    }
    if (wl) {		/* skip if length 0 */
	if (wl != n) {
	    M_cholmod_free_sparse(&ts1, &c);
	    error(_("Length of offset is %d, should be %d"), wl, n);
	}
	for (j = 0; j < n; j++) { /* iterate over columns */
	    double wt = sqrt(wts[j]);
	    for (i = zp[j]; i < zp[j + 1]; i++) zx[i] *= wt;
	}
    }
    
    ts2 = M_cholmod_aat(ts1, (int*)NULL, (size_t)0, 1, &c);
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = M_cholmod_copy(ts2, +1/*upper triangle*/, +1/*values*/, &c);
    M_cholmod_free_sparse(&ts2, &c);
    make_cholmod_sparse_sorted(ts1);
    make_cholmod_sparse_sorted(A);
    ai = (int*)(A->i); zi = (int*)(ts1->i);
    ap = (int*)(A->p); zp = (int*)(ts1->p);
    ax = (double*)(A->x); zx = (double*)(ts1->x);
    for (j = 0; j < m; j++) {
	if (ap[j + 1] != zp[j + 1]) {
	    M_cholmod_free_sparse(&ts1, &c);
	    error(_("pattern mismatch for A and update at column %d"),
		  j + 1);
	}
	for (i = ap[j]; i < ap[j + 1]; i++) {
	    if (ai[i] != zi[i]) {
		M_cholmod_free_sparse(&ts1, &c);
		error(_("pattern mismatch for A and update at column %d"),
		      j + 1);
	    }
	    ax[i] = zx[i];
	}
    }
    M_cholmod_free_sparse(&ts1, &c);
}

/**
 * Create an mer2 object from a list of grouping factors and a list of model
 * matrices.
 *
 * @param fl named list of grouping factors
 * @param ZZt transpose of Z as a sparse matrix
 * @param Xtp transpose of model matrix for the fixed effects
 * @param yp response vector
 * @param REMLp logical scalar indicating if REML is to be used
 * @param ncp integer vector of the number of random effects per level
 *        of each grouping factors
 * @param cnames list of column names of model matrices
 * @param offset numeric vector of offsets
 * @param weights numeric vector of prior weights
 *
 * @return pointer to an mer2 object
 */
SEXP mer2_create(SEXP fl, SEXP ZZt, SEXP Xtp, SEXP yp, SEXP REMLp,
		 SEXP ncp, SEXP cnames, SEXP offset, SEXP wts)
{
    SEXP ST, fnms = getAttrib(fl, R_NamesSymbol),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("mer2"))), xnms;
    cholmod_sparse *A, *Zt, *ZXyt, *ts1, *ts2;
    cholmod_dense *Xy;
    cholmod_factor *L;
    int *Perm, *Gp, *nc = INTEGER(ncp), *dims, *xdims, *zdims,
	i, j, nf = LENGTH(fl), nobs = LENGTH(yp), p, q;
    double *Xt, *fixef, *ranef, *y;
    char *DEVIANCE_NAMES[]={"ML","REML","ldZ","ldX","lr2",""};
    char *DIMS_NAMES[]={"nf","n","p","q","isREML","isGLMM","isNested",""};
				/* record dimensions */
    SET_SLOT(val, lme4_dimsSym, internal_make_named(INTSXP, DIMS_NAMES));
    dims = INTEGER(GET_SLOT(val, lme4_dimsSym));
    dims[nf_POS] = nf;
    dims[n_POS] = nobs;
    dims[isREML_POS] = asLogical(REMLp);
    dims[isGLMM_POS] = FALSE;
    dims[isNest_POS] = TRUE;
				/* Check arguments to be duplicated */
    if (!isReal(yp)) error(_("yp must be a real vector"));
    y = REAL(yp);
    if (!isReal(offset) || (LENGTH(offset) && LENGTH(offset) != nobs))
	error(_("offset must be a real vector of length %d"), nobs);
    SET_SLOT(val, lme4_offsetSym, duplicate(offset));
    if (!isReal(wts) || (LENGTH(wts) && LENGTH(wts) != nobs))
	error(_("wts must be a real vector of length %d"), nobs);
    SET_SLOT(val, lme4_weightsSym, duplicate(wts));
				/*  check flist and create Gp*/
    if (!isNewList(cnames) || LENGTH(cnames) != nf + 2)
	error(_("cnames must be a list of length %d"), nf + 2);
    SET_SLOT(val, lme4_cnamesSym, duplicate(cnames));
    if (!isInteger(ncp) || LENGTH(ncp) != nf)
	error(_("ncp must be an integer vector of length %d"), nf);
    if (!isNewList(fl) || nf < 1) error(_("fl must be a nonempty list"));
    Gp = INTEGER(ALLOC_SLOT(val, lme4_GpSym, INTSXP, nf + 3));
    Gp[0] = 0;
    for (i = 0; i < nf; i++) {
	SEXP fli = VECTOR_ELT(fl, i);
	if (!isFactor(fli) || LENGTH(fli) != nobs)
	    error(_("fl[[%d] must be a factor of length %d"), i+1, nobs);
	Gp[i + 1] = Gp[i] + LENGTH(getAttrib(fli, R_LevelsSymbol)) * nc[i];
    }
    SET_SLOT(val, lme4_flistSym, duplicate(fl));
				/* construct ZXyt matrix */
    if (!isMatrix(Xtp) || !isReal(Xtp))
	error(_("Xtp must be a real matrix"));
    xdims = INTEGER(getAttrib(Xtp, R_DimSymbol));
    if (xdims[1] != nobs) error(_("Xtp must have %d rows"), nobs);
    dims[p_POS] = p = xdims[0]; Xt = REAL(Xtp);
    fixef = REAL(ALLOC_SLOT(val, lme4_fixefSym, REALSXP, p));
    AZERO(fixef, p);
    Gp[nf + 1] = Gp[nf] + p;	/* fixed effects */
    Gp[nf + 2] = Gp[nf + 1] + 1; /* response */
    xnms = VECTOR_ELT(getAttrib(Xtp, R_DimNamesSymbol), 0);
    zdims = INTEGER(GET_SLOT(ZZt, lme4_DimSym));
    if (zdims[1] != nobs) error(_("Zt must have %d columns"), nobs);
    dims[q_POS] = q = zdims[0];
    ranef = REAL(ALLOC_SLOT(val, lme4_ranefSym, REALSXP, q));
    AZERO(ranef, q);
				/* determine permutation */
    Perm = Calloc(q + p + 1, int);
    for (j = 0; j <= (p + q); j++) Perm[j] = j; /* initialize to identity */
    Zt = M_as_cholmod_sparse(ZZt);
				/* check for non-trivial perm */
    if (nf > 1) {
	ts1 = M_cholmod_aat(Zt, (int*)NULL/* fset */,(size_t)0,
			    CHOLMOD_PATTERN, &c);
	ts2 = M_cholmod_copy(ts1, -1/*lower triangle*/, CHOLMOD_PATTERN, &c);
	M_cholmod_free_sparse(&ts1, &c);
	if (!check_nesting(nf, nc, Gp, (int*)(ts2->p))) {
	    dims[isNest_POS] = FALSE;
	    L = M_cholmod_analyze(Zt, &c);
	    if (!L)
		error(_("cholmod_analyze returned with status %d"), c.status);
	    Memcpy(Perm, (int*)(L->Perm), q);
	    M_cholmod_free_factor(&L, &c);
	}
	M_cholmod_free_sparse(&ts2, &c);
    }
				/* create [X;-y]' */
    Xy = M_cholmod_allocate_dense(p + 1, nobs, p + 1, CHOLMOD_REAL, &c);
    for (j = 0; j < nobs; j++) {
	Memcpy(((double*)(Xy->x)) + j * (p + 1), Xt + j * p, p);
	((double*)(Xy->x))[j * (p + 1) + p] = -y[j];
    }
    ts1 = M_cholmod_dense_to_sparse(Xy, TRUE, &c);
    M_cholmod_free_dense(&Xy, &c);
    ts2 = M_cholmod_vertcat(Zt, ts1, TRUE, &c);
    M_cholmod_free_sparse(&ts1, &c);
    SET_SLOT(val, lme4_ZXytSym,
	     M_chm_sparse_to_SEXP(ts2, 0, 0, 0, "", R_NilValue));
    M_cholmod_free_sparse(&ts2, &c);
    ZXyt = M_as_cholmod_sparse(GET_SLOT(val, lme4_ZXytSym));
				/*  Create and store A's pattern */
    ts1 = M_cholmod_aat(ZXyt, (int*)NULL, (size_t)0, 1, &c);
    ts2 = M_cholmod_copy(ts1, +1/*upper triangle*/, +1/*values*/, &c);
    M_cholmod_free_sparse(&ts1, &c);
    make_cholmod_sparse_sorted(ts2);
    SET_SLOT(val, lme4_ASym,
	     M_chm_sparse_to_SEXP(ts2, 0, 0, 0, "", R_NilValue));
    M_cholmod_free_sparse(&ts2, &c);
				/* update A using weights and offset */
    A = M_as_cholmod_sparse(GET_SLOT(val, lme4_ASym));
    internal_update_A(ZXyt, wts, offset, A);
				/* allocate, populate and initialize ST */
    ST = ALLOC_SLOT(val, lme4_STSym, VECSXP, nf);
    setAttrib(ST, R_NamesSymbol, duplicate(fnms));
    for (i = 0; i < nf; i++)
	SET_VECTOR_ELT(ST, i, allocMatrix(REALSXP, nc[i], nc[i]));
    internal_mer2_initial(ST, Gp, A);
    j = c.supernodal;		/* force supernodal for non-nested */
    if (!dims[isNest_POS]) c.supernodal = CHOLMOD_SUPERNODAL;
    i = c.nmethods;
    c.nmethods = 1;		/* force user-specified permutation */
    				/* Create L  with user-specified perm */
    L = M_cholmod_analyze_p(A, Perm, (int*)NULL, (size_t)0, &c);
    if (!L)
	error(_("cholmod_analyze_p returned with status %d"), c.status);
				/* restore previous settings */
    c.nmethods = i;
    c.supernodal = j;
    				/* initialize and store L */
    SET_SLOT(val, lme4_devianceSym,
	     internal_make_named(REALSXP, DEVIANCE_NAMES));
    internal_update_L(REAL(GET_SLOT(val, lme4_devianceSym)),
			  dims, Gp, ST, A, L);
    SET_SLOT(val, lme4_LSym, M_chm_factor_to_SEXP(L, 1));

    Free(A); Free(ZXyt); Free(Perm); Free(Zt);
    UNPROTECT(1);
    return val;
}

/**
 * Calculate the number of parameters in ST
 *
 * @param nf number of factors
 * @param nc number of random effects per level of the factor
 *
 * @return the total number of parameters in ST
 */
static R_INLINE int mer2_npar(SEXP ST)
{
    int ans = 0, i, nf = LENGTH(ST);

    for (i = 0; i < nf; i++) {
	int nci = INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol))[0];
	ans += (nci * (nci + 1))/2;
    }
    return ans;
}

/**
 * Extract the parameters from ST
 *
 * @param ST ST slot from an mer2 object
 * @param pars double vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double
*internal_mer2_getPars(SEXP ST, double *pars)
{
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++) pars[pos++] = st[j * ncp1];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		pars[pos++] = st[k + j * nci];
    }
    return pars;
}

/**
 * Extract the parameters from the ST slot of an mer2 object
 *
 * @param x an mer2 object
 *
 * @return pointer to a REAL vector
 */
SEXP mer2_getPars(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    SEXP ans = PROTECT(allocVector(REALSXP, mer2_npar(ST)));

    internal_mer2_getPars(ST, REAL(ans));
    UNPROTECT(1); 
    return ans;
}

/**
 * Update the ST slot of an mer2 object
 *
 * @param pars double vector of the appropriate length
 * @param ST ST slot from an mer2 object
 *
 * @return pointer to the updated ST object
 */
static SEXP
internal_mer2_setPars(const double *pars, SEXP ST)
{
    int i, nf = LENGTH(ST), pos = 0;
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *st = REAL(STi);
	int nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int j, k, ncp1 = nci + 1;

	for (j = 0; j < nci; j++) st[j * ncp1] = pars[pos++];
	for (j = 0; j < (nci - 1); j++)
	    for (k = j + 1; k < nci; k++)
		st[k + j * nci] = pars[pos++];
    }
    return ST;
}

/**
 * Update the ST slot of an mer2 object from a REAL vector of
 * parameters and update the Cholesky factorization
 *
 * @param x an mer2 object
 * @param pars a REAL vector of the appropriate length
 *
 * @return x
 */
SEXP mer2_setPars(SEXP x, SEXP pars)
{
    cholmod_sparse *A = M_as_cholmod_sparse(GET_SLOT(x, lme4_ASym));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int npar = mer2_npar(ST);

    if (!isReal(pars) || LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    internal_mer2_setPars(REAL(pars), ST);
    internal_update_L(REAL(GET_SLOT(x, lme4_devianceSym)),
		      INTEGER(GET_SLOT(x, lme4_dimsSym)), 
		      INTEGER(GET_SLOT(x, lme4_GpSym)), ST, A, L);
    Free(A); Free(L);
    return x;
}

/* FIXME: Should I start using the name "discrepancy" instead of
   "deviance"? */
/**
 * Extract the deviance from an mer2 object
 *
 * @param x an mer2 object
 * @param which scalar integer < 0 => REML, 0 => native, > 0 => ML
 *
 * @return scalar REAL value
 */
SEXP mer2_deviance(SEXP x, SEXP which)
{
    int w = asInteger(which);
    int POS = (w < 0 || (!w && isREML(x))) ? REML_POS : ML_POS; 

    return ScalarReal(REAL(GET_SLOT(x, lme4_devianceSym))[POS]);
}

    
/**
 * Update the contents of the fixef and ranef slots
 *
 * @param x an mer2 object
 *
 * @return R_NilValue
 */
SEXP mer2_update_effects(SEXP x)
{
    int *dd = INTEGER(GET_SLOT(x, lme4_dimsSym));
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)), *bstbx;
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    cholmod_dense *bstarb;

    bstarb = internal_mer2_effects(L);
    bstbx = (double*)(bstarb->x);
    Memcpy(b, bstbx, dd[q_POS]);
    Memcpy(REAL(GET_SLOT(x, lme4_fixefSym)), bstbx + dd[q_POS], dd[p_POS]);
    M_cholmod_free_dense(&bstarb, &c);
    TS_mult(INTEGER(GET_SLOT(x, lme4_GpSym)),
			 GET_SLOT(x, lme4_STSym), b);
    Free(L);
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
internal_mer2_sigma(int REML, const int* dims, const double* deviance)
{
    return sqrt(exp(deviance[lr2_POS])/
		((double)(dims[n_POS] - (REML ? dims[p_POS] : 0))));
}

static int internal_mer2_optimize(SEXP x, int verb)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    cholmod_sparse *A = M_as_cholmod_sparse(GET_SLOT(x, lme4_ASym));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, j,
	n = mer2_npar(ST), nf = length(ST), pos = 0;
    int REML = dims[isREML_POS],
	liv = S_iv_length(OPT, n), lv = S_v_length(OPT, n);
    int *iv = Calloc(liv, int);
    double *b = Calloc(2 * n, double), *d = Calloc(n, double),
	*deviance = REAL(GET_SLOT(x, lme4_devianceSym)),
	*g = (double*)NULL, *h = (double*)NULL,
	*v = Calloc(lv, double),
	*xv = internal_mer2_getPars(ST, Calloc(n, double)),
	fx = R_PosInf;

    S_Rf_divset(OPT, iv, liv, lv, v);
    if (verb) iv[OUTLEV] = 1;
    for (i = 0; i < n; i++) {
	b[2*i] = R_NegInf; b[2*i+1] = R_PosInf; d[i] = 1;
    }
    for (i = 0; i < nf; i++) {
	int nc = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
	for (j = 0; j < nc; j++) b[pos + 2*j] = 0;
	pos += nc * (nc + 1);
    }
    do {
	S_nlminb_iterate(b, d, fx, g, h, iv, liv, lv, n, v, xv);
	internal_mer2_setPars(xv, ST);
	internal_update_L(deviance, dims, Gp, ST, A, L);
	fx = deviance[REML ? REML_POS : ML_POS];
    } while (iv[0] < 3);
    i = iv[0];
    Free(iv); Free(v); Free(xv); Free(b); Free(d);
    return i;
}

SEXP mer2_optimize(SEXP x, SEXP verb)
{
    return ScalarInteger(internal_mer2_optimize(x, asInteger(verb)));
}

/**
 * Extract the estimate of the scale factor from an mer2 object
 *
 * @param x an mer2 object
 * @param which scalar integer (< 0 => REML, 0 => native, > 0 => ML)
 *
 * @return scalar REAL value
 */
SEXP mer2_sigma(SEXP x, SEXP which)
{
    int w = asInteger(which);
		
    return ScalarReal(internal_mer2_sigma(w < 0 || (!w && isREML(x)),
					  INTEGER(GET_SLOT(x, lme4_dimsSym)),
					  REAL(GET_SLOT(x, lme4_devianceSym))));
}

/**
 * Extract the unscaled lower Cholesky factor of the relative
 * variance-covariance matrix for the fixed-effects in an mer2 object.
 *
 * @param x an mer2 object
 *
 * @return a REAL p by p lower triangular matrix (it's a matrix, not a Matrix)
 */

SEXP mer2_vcov(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *select;
    int i, p = dims[p_POS], q = dims[q_POS];
    SEXP ans = PROTECT(allocMatrix(REALSXP, p, p));

    if (p) {
	cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
	 /* need a copy because factor_to_sparse changes 1st arg */
	cholmod_factor *Lcp = M_cholmod_copy_factor(L, &c);
	cholmod_sparse *Lsp, *Lred; /* sparse and reduced-size sparse */
	cholmod_dense *Ld;

	if (!(Lcp->is_ll))
	    if (!M_cholmod_change_factor(Lcp->xtype, 1, 0, 1, 1, Lcp, &c))
		error(_("cholmod_change_factor failed with status %d"), c.status);
	Lsp = M_cholmod_factor_to_sparse(Lcp, &c);
	M_cholmod_free_factor(&Lcp, &c);
	select = Calloc(p, int);
	for (i = 0; i < p; i++) select[i] = q + i;
	Lred = M_cholmod_submatrix(Lsp, select, p, select, p,
				 1 /* values */, 1 /* sorted */, &c);
	M_cholmod_free_sparse(&Lsp, &c);
	Ld = M_cholmod_sparse_to_dense(Lred, &c);
	M_cholmod_free_sparse(&Lred, &c);
	Memcpy(REAL(ans), (double*)(Ld->x), p * p);
	M_cholmod_free_dense(&Ld, &c);
/* FIXME: This does not allow for a possible P_X permutation  */
	F77_CALL(dtrtri)("L", "N", &p, REAL(ans), &p, &i);
	if (i)
	    error(_("Lapack routine dtrtri returned error code %d"), i);
	Free(L); Free(select);
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Extract the conditional modes of the random effects as a list of matrices
 *
 * @param x Pointer to an mer object
 *
 * @return a list of matrices containing the conditional modes of the
 * random effects
 */
SEXP mer2_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
        cnames = GET_SLOT(x, lme4_cnamesSym),
	flist = GET_SLOT(x, lme4_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ii, jj, nf = LENGTH(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, lme4_ranefSym));

    mer2_update_effects(x);
    setAttrib(val, R_NamesSymbol,
	      duplicate(getAttrib(flist, R_NamesSymbol)));
    for (i = 0; i < nf; i++) {
	SEXP nms, rnms = getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol);
	int nci = INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol))[0];
	int mi = length(rnms);
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
 * Extract the posterior variances of the random effects
 *
 * @param x pointer to a mer object
 *
 * @return pointer to a list of arrays
 */
SEXP mer2_postVar(SEXP x)
{
    double *deviance = REAL(GET_SLOT(x, lme4_devianceSym)), one = 1;
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, j, nf = dims[nf_POS], p = dims[p_POS], q = dims[q_POS];
    int ppq = p + q;
    double sc = internal_mer2_sigma(isREML(x), dims, deviance);
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym)),
	*Lcp = (cholmod_factor*)NULL;
    cholmod_sparse *rhs, *B, *Bt, *BtB;
    cholmod_dense *BtBd;
    int *Perm = (int*)(L->Perm), *iperm = Calloc(ppq, int),
	*fset = Calloc(ppq, int);
    SEXP ST = GET_SLOT(x, lme4_STSym),
	ans = PROTECT(allocVector(VECSXP, nf));
    
    for (j = 0; j < ppq; j++) {
	iperm[Perm[j]] = j;
	fset[j] = j;
    }
    if (!L->is_ll) {
	Lcp = M_cholmod_copy_factor(L, &c);
	Free(L);
	L = Lcp;
	j = M_cholmod_change_factor(CHOLMOD_REAL, TRUE/*ll*/,
				    FALSE/*super*/, TRUE/*packed*/,
				    TRUE/*sorted*/, L, &c);
	if (!j) error(_("cholmod_change_factor failed"));
    }
    sc = sc * sc;		/* variance scale factor */
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int j, k, kk, nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	int nlev = (Gp[i + 1] - Gp[i])/nci;
	SEXP ansi = PROTECT(alloc3Darray(REALSXP, nci, nci, nlev));
	int ncip1 = nci + 1, ncisqr = nci * nci;
	double *vv = REAL(ansi),
	    *st = Memcpy(Calloc(ncisqr, double), REAL(STi), ncisqr);

	SET_VECTOR_ELT(ans, i, ansi); UNPROTECT(1);
	AZERO(vv, ncisqr * nlev);
	rhs = M_cholmod_allocate_sparse((size_t)(ppq + 1),
					(size_t) nci, (size_t) nci,
					1/*sorted*/, 1/*packed*/,
					0/*stype*/, CHOLMOD_REAL, &c);
	((int*)(rhs->p))[0] = 0;
	for (k = 0; k < nci; k++) {
	    ((double*)(rhs->x))[k] = 1.;
	    ((int*)(rhs->p))[k + 1] = k + 1;
	}
	for (k = 0; k < nci; k++) {
	    double mult = st[k * ncip1];
	    for (kk = k + 1; kk < nci; kk++)
		st[kk + k * nci] *= mult;
	}
	for (j = 0; j < nlev; j++) {
	    int *ip, *pp, base = Gp[i] + j * nci;
	    double *xp;
	    
	    for (k = 0; k < nci; k++)
		((int*)(rhs->i))[k] = iperm[base + k];
	    B = M_cholmod_spsolve(CHOLMOD_L, L, rhs, &c);
	    ip = (int*)(B->i);
	    pp = (int*)(B->p);
	    xp = (double*)(B->x);
	    if (nci == 1) {
		for (k = 0; k < pp[1]; k++)
		    if (ip[k] < ppq) vv[j] += xp[k] * xp[k];
		vv[j] *= sc * st[0] * st[0];
	    } else {
		double *vvj = vv + j * ncisqr;
		Bt = M_cholmod_transpose(B, TRUE/*values*/, &c);
		BtB = M_cholmod_aat(Bt, fset, (size_t)ppq, 1/*mode*/,&c);
		M_cholmod_free_sparse(&Bt, &c);
		BtBd = M_cholmod_sparse_to_dense(BtB, &c);
		M_cholmod_free_sparse(&BtB, &c);
		Memcpy(vvj, (double*)(BtBd->x), ncisqr);
		M_cholmod_free_dense(&BtBd, &c);
		F77_CALL(dtrmm)("L", "L", "N", "N", &nci, &nci,
				&one, st, &nci, vvj, &nci);
		F77_CALL(dtrmm)("R", "L", "T", "N", &nci, &nci,
				&sc, st, &nci, vvj, &nci);
	    }
	    M_cholmod_free_sparse(&B, &c);
	}
	M_cholmod_free_sparse(&rhs, &c);
	Free(st);
    }
    if (L == Lcp) M_cholmod_free_factor(&L, &c); else Free(L);
    Free(iperm); Free(fset);
    UNPROTECT(1);
    return ans;
}

/**
 * Update the coefficients beta and bstar given the current L
 *
 * @param sigma current standard deviation of the noise term
 * @param L current factorization of A*
 * @param bbsthat current conditional modes
 * @param bbstnew array to receive updated coefficients
 *
 * @return squared length of spherical Gaussian deviate
 */
static double
internal_betabst_update(double sigma, cholmod_factor *L,
		       cholmod_dense *bbsthat,
		       cholmod_dense *bbstnew)
{
    int i, ppq1 = L->n;
    int ppq = ppq1 - 1;
    double *bx, *hx = (double*)(bbsthat->x),
	*nx = (double*)(bbstnew->x), ans = 0;
    cholmod_dense *chb;

    nx[ppq] = 0.;
    for (i = 0; i < ppq; i++) {
	double nr = norm_rand();
	ans += nr * nr;
	nx[i] = sigma * nr;
    }
    				/* chb := L^{-T} %*% bnew */
    chb = M_cholmod_solve(CHOLMOD_Lt, L, bbstnew, &c);
    bx = (double*)(chb->x);
    for (i = 0; i < ppq; i++) nx[i] = bx[i] + hx[i];
    M_cholmod_free_dense(&chb, &c);

    return ans;
}

/**
 * Update the relative variance-covariance matrices by sampling from a
 * Wishart distribution with scale factor determined by the current
 * sample of random effects.
 *
 * @param sigma current value of sigma
 * @param Gp vector of group pointers
 * @param bnew current sample from the random effects
 * @param ST list of factorizations to update
 */
static void
internal_ST_update(double sigma, int trans, const int *Gp,
		   const double *bnew, SEXP ST, double *dest)
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
	    double ssq = 0;
	    for (j = 0; j < nlev; j++) ssq += bnew[Gp[i] + j] * bnew[Gp[i] + j];
	    sd = sqrt(ssq/rchisq((double)nlev)) * st[0];
	    st[0] = sd/sigma;
	    dest[pos++] = (trans ? 2 * log(sd) : sd * sd);
	} else {
	    int ncip1 = nci + 1, ncisqr = nci * nci;
	    double *scal = Calloc(ncisqr, double), /* factor of scale matrix */
		*wfac = Calloc(ncisqr, double); /* factor of Wishart variate */
				/* create the scale matrix (close to identity) */
	    AZERO(scal, ncisqr);
	    F77_CALL(dsyrk)("L", "N", &nci, &nlev, &sigsqinv, bnew + Gp[i], &nci,
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

/**
 * Generate a Markov-Chain Monte Carlo sample from a fitted
 * linear mixed model.
 *
 * @param x pointer to an mer2 object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP mer2_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbosep, SEXP deviancep)
{
    SEXP ST = GET_SLOT(x, lme4_STSym), Pars = PROTECT(mer2_getPars(x)), ans;
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    cholmod_sparse *A = M_as_cholmod_sparse(GET_SLOT(x, lme4_ASym));
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	*Gp = INTEGER(GET_SLOT(x, lme4_GpSym));
    int dPOS = dims[REML_POS] ? REML_POS : ML_POS, i, j,
	n = dims[n_POS], nsamp = asInteger(nsampp),
	p = dims[p_POS], q = dims[q_POS],
	saveb = asLogical(savebp), trans = asLogical(transp),
	verbose = asLogical(verbosep), dev = asLogical(deviancep);
    int qpp1 = q + p + 1;
    double
	*deviance = REAL(GET_SLOT(x, lme4_devianceSym)),
	*ansp, df = n - (dims[REML_POS] ? p : 0);
    int nrbase = p + 1 + mer2_npar(ST); /* rows always included */
    int nrtot = nrbase + dev + (saveb ? q : 0);
    cholmod_dense *chhat,
	*chnew = M_cholmod_allocate_dense(qpp1, 1, qpp1, CHOLMOD_REAL, &c);
    double *bstar = (double*) NULL, *nx = (double*)(chnew->x);

    if (nsamp <= 0) nsamp = 1;
    ans = PROTECT(allocMatrix(REALSXP, nrtot, nsamp));
    ansp = REAL(ans);
    for (i = 0; i < nrtot * nsamp; i++) ansp[i] = NA_REAL;
    if (saveb) bstar = Calloc(q, double);
    GetRNGstate();
    if (verbose) Rprintf("%12s %14s\n", "sigma", "deviance");

    for (i = 0; i < nsamp; i++) {
	double *col = ansp + i * nrtot, sigma;
				/* simulate and store new value of sigma */
	sigma = exp(deviance[lr2_POS]/2)/sqrt(rchisq(df));
	col[p] = (trans ? 2 * log(sigma) : sigma * sigma);
	/* simulate new fixed and random effects */
				/* Evaluate conditional estimates */
	chhat = internal_mer2_effects(L);
	internal_betabst_update(sigma, L, chhat, chnew);
	M_cholmod_free_dense(&chhat, &c);
				/* Store beta */
	for (j = 0; j < p; j++) col[j] = nx[q + j];
	if (saveb) {		/* Optionally store b */
	    TS_mult(Gp, ST, Memcpy(bstar, nx, q));
	    for (j = 0; j < q; j++) col[nrbase + dev + j] = bstar[j];
	}
	internal_ST_update(sigma, trans, Gp, nx, ST, col + p + 1);
				/* Refactor and evaluate deviance */
	internal_update_L(deviance, dims, Gp, ST, A, L);
				/* store new variance-covariance parameters */
				/* optionally store deviance */
	/* FIXME: This deviance should be for sigma, beta, b, ST */
	if (dev) col[nrbase] = deviance[dPOS]; 
	if (verbose) Rprintf("%12.6g %14.8g\n", sigma, deviance[dPOS]);
    }
    PutRNGstate();
    if (saveb) Free(bstar);
    Free(L);
    M_cholmod_free_dense(&chnew, &c);
				/* Restore pars, refactor, etc. */
    mer2_setPars(x, Pars);
    UNPROTECT(2);
    return ans;
}

/**
 * Check validity of an mer2 object.
 *
 * @param x Pointer to an mer2 object
 *
 * @return TRUE if the object is a valid mer2 object, else a string
 * describing the nature of the violation.
 */
SEXP mer2_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	fixefP = GET_SLOT(x, lme4_fixefSym),
	flistP = GET_SLOT(x, lme4_flistSym),
	offsetP = GET_SLOT(x, lme4_offsetSym),
	ranefP = GET_SLOT(x, lme4_ranefSym),
	weightsP = GET_SLOT(x, lme4_weightsSym) ;
    cholmod_sparse *A = M_as_cholmod_sparse(GET_SLOT(x, lme4_ASym)),
	*ZXyt = M_as_cholmod_sparse(GET_SLOT(x, lme4_ZXytSym));
    cholmod_factor *L = M_as_cholmod_factor(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nf = dd[nf_POS], n = dd[n_POS], p = dd[p_POS], q = dd[q_POS];
    int i, nq, ppq1 = p + q + 1;

				/* check lengths */
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length nf"));
    if (LENGTH(GpP) != (nf + 3))
	return mkString(_("Slot Gp must have length nf + 3"));
    if (Gp[0] != 0 || Gp[nf + 2] != ppq1)
	return mkString(_("Gp[1] != 0 or Gp[nf+3] != p+q+1"));
    if (LENGTH(fixefP) != p)
	return mkString(_("Slot fixef must have length p"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length q"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or n"));
    if (LENGTH(offsetP) && LENGTH(offsetP) != n)
	return mkString(_("Slot offset must have length 0 or n"));
    if (LENGTH(devianceP) != (lr2_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (lr2_POS + 1))
	return mkString(_("deviance slot not named or incorrect length"));
    if (ZXyt->nrow != ppq1 || ZXyt->ncol != n)
	return mkString(_("Slot ZXyt must have dimensions (p+q+1) by n"));
    if (A->nrow != ppq1 || A->ncol != ppq1 || A->stype <= 0)
	return mkString(_("Slot A must be symmetric (upper) of size (p+q+1)"));
    if (L->n != ppq1 || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size (p+q+1)"));

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

    Free(A); Free(L); Free(ZXyt);
    return ScalarLogical(1);
}

SEXP lmer2_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp,
		  SEXP nc, SEXP cnames, SEXP offset, SEXP wts, SEXP fr,
		  SEXP terms, SEXP call)
{
    SEXP mer = PROTECT(mer2_create(fl, Zt, Xp, yp, REMLp, nc, cnames, offset, wts)),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("lmer2")));
    SEXP nms = getAttrib(GET_SLOT(MAKE_CLASS("mer2"), install("slots")), R_NamesSymbol);
    int i, n = LENGTH(nms);

    for (i = 0; i < n; i++) {
	SEXP sym = install(CHAR(STRING_ELT(nms, i)));
	SET_SLOT(val, sym, GET_SLOT(mer, sym));
    }
    SET_SLOT(val, install("frame"), fr);
    SET_SLOT(val, install("terms"), terms);
    SET_SLOT(val, install("call"), call);

    UNPROTECT(2);
    return val;
}
