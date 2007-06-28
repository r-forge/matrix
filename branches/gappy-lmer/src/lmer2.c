#include "lmer2.h"

/* Functions for the lmer2 representation */

				/* positions in the deviance vector */
enum devP {ML_POS=0, REML_POS, ldZ_POS, ldX_POS, lr2_POS, bqd_POS, Sdr_POS};
			/* {"ML", "REML", "ldZ", "ldX", "lr2", "bQuad" "sumDevR" ""} */
				/* positions in the dims vector */
enum dimP {nf_POS=0, n_POS, p_POS, q_POS, s_POS, np_POS, isREML_POS, famType_POS, isNest_POS};
	      /* {"nf", "n", "p", "q", "s", "np", "isREML", "famType", "isNested"} */

#define isREML(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isREML_POS]
#define isGLMM(x) (INTEGER(GET_SLOT(x, lme4_dimsSym))[famType_POS] >= 0)
#define isNested(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isNest_POS]

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
static SEXP
alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface)
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

/**
 * Return the element of a given name from a named list
 *
 * @param list
 * @param nm name of desired element
 *
 * @return element of list with name nm
 */
static SEXP R_INLINE getListElement(SEXP list, char *nm) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    if (!isNull(names))
	for (i = 0; i < LENGTH(names); i++)
	    if (!strcmp(CHAR(STRING_ELT(names, i)), nm))
		return(VECTOR_ELT(list, i));
    return R_NilValue;
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
static int check_nesting(int nf, SEXP ST, const int Gp[], const int p[])
{
    int **cnz = Calloc(nf, int*), *nc = Calloc(nf, int), ans = 1, i, j, k, nct;

    for (i = 0, nct = 0; i < nf; i++) { /* total number of columns */
	nc[i] = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
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
    Free(cnz); Free(nc);
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
 * Populate the st, nc and nlev arrays.  Return an indicator of
 * non-trivial T.
 *
 * @param nf number of random-effects terms
 * @param ST pointer to a list (length nf) of matrices
 * @param Gp pointer to integer vector (length nf + 1) of group
 * pointers
 * @param st length nf array of (double*) pointers to be filled with
 * pointers to the contents of the matrices in ST
 * @param nc length nf array of ints to be filled with the number of
 * columns
 * @param nlev length nf array of ints to be filled with the number of
 * levels of the grouping factor for each term
 * 
 * @return 0 when T == the identity matrix (all nc[i] are 1) else 1
 */
static int			/* populate the st, nc and nlev arrays */
st_nc_nlev(const SEXP ST, const int *Gp, double **st, int *nc, int *nlev)
{
    int ans = 0, i, nf = LENGTH(ST);

    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	int nci = *INTEGER(getAttrib(STi, R_DimSymbol));

	if (nci > 1) ans = 1;
	st[i] = REAL(STi); nc[i] = nci; nlev[i] = (Gp[i + 1] - Gp[i])/nci;
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
 * Change the x slot in Vt according to ST'Zt
 *
 * @param Zt transpose of Z as a dgCMatrix object
 * @param ST compounded matrices S* and T*
 * @param Gp nf + 1 vector of group pointers
 * @param Vt transpose of V as a dgCMatrix object
 *
 * @return square of determinant of L
 */
static double
internal_VtL_update(SEXP Zt, SEXP ST, const int *Gp,
		    SEXP Vt, CHM_FR L)
{
    CHM_SP cVt = AS_CHM_SP(Vt);
    int *vi = (int*)(cVt->i), *vp = (int*)(cVt->p),
	fll = c.final_ll, nf = LENGTH(ST), p;
    int *nc = (int*)alloca(nf * sizeof(int)),
	*nlev = (int*)alloca(nf * sizeof(int)),
	nnz = vp[cVt->ncol];
    double **st = (double**)alloca(nf * sizeof(double*)),
	*vx = (double*)(cVt->x), one[] = {1,0};
    
    R_CheckStack();
    if (st_nc_nlev(ST, Gp, st, nc, nlev)) { /* multiply by T' */
	error(_("Code not yet written"));
    } else Memcpy(vx, REAL(GET_SLOT(Zt, lme4_xSym)), nnz);
    for (p = 0; p < nnz; p++) {
	int i = Gp_grp(vi[p], nf, Gp);
	vx[p] *= st[i][((vi[p] - Gp[i]) / nlev[i]) * (nc[i] + 1)];
    }
    c.final_ll = L->is_ll;
    if (!M_cholmod_factorize_p(cVt, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    c.final_ll = fll;

    return chm_log_det2(L);
}

SEXP VtL_update(SEXP x)
{
    return ScalarReal(
	internal_VtL_update(GET_SLOT(x, install("Zt")),
			    GET_SLOT(x, lme4_STSym),
			    INTEGER(GET_SLOT(x, lme4_GpSym)),
			    GET_SLOT(x, install("Vt")),
			    AS_CHM_FR(GET_SLOT(x, lme4_LSym))));
}

/**
 * Multiply a dense matrix by the virtual T and S matrices represented by ST
 *
 * @param Gp vector of group pointers
 * @param ST compounded matrices S* and T*
 * @param b vector of random effects to be transformed from u to b
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
    int *nc = (int*)alloca(nf * sizeof(int)),
	*nlev = (int*)alloca(nf * sizeof(int)),
	nnz = INTEGER(GET_SLOT(Zt, lme4_pSym))[Zdims[1]];
    double *rowsqr = (double*)alloca(Zdims[0] * sizeof(double)),
	**st = (double**)alloca(nf * sizeof(double*)),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));
    
    st_nc_nlev(ST, Gp, st, nc, nlev);
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
 * Evaluate the effects in an lmer2 representation
 *
 * @param L factorization
 *
 * @return cholmod_dense object with \hat{b*} in the first q positions
 * and \hat{\beta} in the next p positions.
 */
static CHM_DN
internal_lmer2_effects(CHM_FR L)
{
    int i, nt = (L->n);
    CHM_DN X, B = M_cholmod_allocate_dense(L->n, 1, L->n, CHOLMOD_REAL, &c);
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
chm_log_abs_det2(double *ans, int nans, const int *c, const CHM_FR F)
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
internal_deviance(double *d, const int *dims, const CHM_FR L)
{
    int n = dims[n_POS], p = dims[p_POS], q = dims[q_POS];
    int c[] = {0,  q, p + q, p + q + 1};
    double dn = (double) n, dnmp = (double)(n - p);
    
    chm_log_abs_det2(d + ldZ_POS, lr2_POS + 1 - ldZ_POS, c, L);
    d[ML_POS] = d[ldZ_POS] + dn * (1. + d[lr2_POS] + log(2. * PI / dn));
    d[REML_POS] = d[ldZ_POS] + d[ldX_POS] + dnmp *
	(1. + d[lr2_POS] + log(2. * PI / dnmp));
    d[bqd_POS] = d[Sdr_POS] = 0.;
    return d;
}

/**
 * Update A to A* and evaluate its numeric factorization in L.
 *
 * @param deviance Hold the result
 * @param dims dimensions
 * @param Gp length nf+3 vector of group pointers for the rows of A
 * @param ST pointers to the nf ST factorizations of the diagonal
 *     elements of Sigma 
 * @param A symmetric matrix of size Gp[nf+2]
 * @param F factorization to be modified
 *
 */
static void
internal_update_L(double *deviance, int *dims, const int *Gp,
		  SEXP ST, CHM_SP A, CHM_FR L)
{
    CHM_SP Ac = M_cholmod_copy_sparse(A, &c);
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
	    
/* FIXME: calculate and store maxrows in lmer2_create */
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

static R_INLINE void
make_cholmod_sparse_sorted(CHM_SP A)
{
    if(!A->sorted) {
	int i = M_cholmod_sort(A, &c); 
	if(!i)
	    error(_("cholmod_sort returned error code %d"),i);
    }
}

/**
 * Extract the parameters from ST
 *
 * @param ST ST slot from an lmer2 object
 * @param pars double vector of the appropriate length
 *
 * @return pointer to the parameter vector
 */
static double
*internal_lmer2_getPars(SEXP ST, double *pars)
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
 * Extract the parameters from the ST slot of an lmer2 object
 *
 * @param x an lmer2 object
 *
 * @return pointer to a REAL vector
 */
SEXP lmer2_getPars(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    SEXP ans = PROTECT(allocVector(REALSXP,
				   INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS]));

    internal_lmer2_getPars(ST, REAL(ans));
    UNPROTECT(1); 
    return ans;
}

/**
 * Update the ST slot of an lmer2 object
 *
 * @param pars double vector of the appropriate length
 * @param ST ST slot from an lmer2 object
 *
 * @return pointer to the updated ST object
 */
static SEXP
internal_lmer2_setPars(const double *pars, SEXP ST)
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
 * Update the ST slot of an lmer2 object from a REAL vector of
 * parameters and update the Cholesky factorization
 *
 * @param x an lmer2 object
 * @param pars a REAL vector of the appropriate length
 *
 * @return x
 */
SEXP lmer2_setPars(SEXP x, SEXP pars)
{
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int npar = INTEGER(GET_SLOT(x, lme4_dimsSym))[np_POS];

    R_CheckStack();
    if (!isReal(pars) || LENGTH(pars) != npar)
	error(_("pars must be a real vector of length %d"), npar);
    internal_lmer2_setPars(REAL(pars), ST);
    internal_update_L(REAL(GET_SLOT(x, lme4_devianceSym)),
		      INTEGER(GET_SLOT(x, lme4_dimsSym)), 
		      INTEGER(GET_SLOT(x, lme4_GpSym)), ST, A, L);
    return x;
}

/**
 * Extract the deviance from an lmer2 object
 *
 * @param x an lmer2 object
 * @param which scalar integer < 0 => REML, 0 => native, > 0 => ML
 *
 * @return scalar REAL value
 */
SEXP lmer2_deviance(SEXP x, SEXP which)
{
    int w = asInteger(which);
    int POS = (w < 0 || (!w && isREML(x))) ? REML_POS : ML_POS; 

    return ScalarReal(REAL(GET_SLOT(x, lme4_devianceSym))[POS]);
}

    
/**
 * Update the contents of the fixef and ranef slots
 *
 * @param x an lmer2 object
 *
 * @return R_NilValue
 */
SEXP lmer2_update_effects(SEXP x)
{
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i;
    int q = dims[q_POS];
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)), *bstbx,
	*dev = REAL(GET_SLOT(x, lme4_devianceSym));
    CHM_DN bstarb;

    R_CheckStack();
    bstarb = internal_lmer2_effects(L);
    bstbx = (double*)(bstarb->x);
    Memcpy(b, bstbx, q);
    for (i = 0, dev[bqd_POS] = 0; i < q; i++) /* accumulate ssqs of bstar */
	dev[bqd_POS] += bstbx[i] * bstbx[i];
    /* FIXME: apply the permutation when copying */
    Memcpy(REAL(GET_SLOT(x, lme4_fixefSym)), bstbx + q, dims[p_POS]);
    M_cholmod_free_dense(&bstarb, &c);
    TS_mult(INTEGER(GET_SLOT(x, lme4_GpSym)),
	    GET_SLOT(x, lme4_STSym), b);
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
internal_lmer2_sigma(int REML, const int* dims, const double* deviance)
{
    return sqrt(exp(deviance[lr2_POS])/
		((double)(dims[n_POS] - (REML ? dims[p_POS] : 0))));
}


/**
 * Extract the estimate of the scale factor from an lmer2 object
 *
 * @param x an lmer2 object
 * @param which scalar integer (< 0 => REML, 0 => native, > 0 => ML)
 *
 * @return scalar REAL value
 */
SEXP lmer2_sigma(SEXP x, SEXP which)
{
    int w = asInteger(which);
		
    return ScalarReal(internal_lmer2_sigma(w < 0 || (!w && isREML(x)),
					  INTEGER(GET_SLOT(x, lme4_dimsSym)),
					  REAL(GET_SLOT(x, lme4_devianceSym))));
}

/**
 * Extract the unscaled lower Cholesky factor of the relative
 * variance-covariance matrix for the fixed-effects in an lmer2 object.
 *
 * @param x an lmer2 object
 *
 * @return a REAL p by p lower triangular matrix (it's a matrix, not a Matrix)
 */

SEXP lmer2_vcov(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *select;
    int i, p = dims[p_POS], q = dims[q_POS];
    SEXP ans = PROTECT(allocMatrix(REALSXP, p, p));

    if (p) {
	CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
	 /* need a copy because factor_to_sparse changes 1st arg */
	CHM_FR Lcp = M_cholmod_copy_factor(L, &c);
	CHM_SP Lsp, Lred; /* sparse and reduced-size sparse */
	CHM_DN Ld;

	R_CheckStack();
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
	Free(select);
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
SEXP lmer2_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
        cnames = GET_SLOT(x, lme4_cnamesSym),
	flist = GET_SLOT(x, lme4_flistSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ii, jj, nf = LENGTH(flist);
    SEXP val = PROTECT(allocVector(VECSXP, nf));
    double *b = REAL(GET_SLOT(x, lme4_ranefSym));

    lmer2_update_effects(x);
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
SEXP lmer2_postVar(SEXP x)
{
    double *deviance = REAL(GET_SLOT(x, lme4_devianceSym)), one = 1;
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, j, nf = dims[nf_POS], p = dims[p_POS], q = dims[q_POS];
    int ppq = p + q;
    double sc = internal_lmer2_sigma(isREML(x), dims, deviance);
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym)), Lcp = (CHM_FR) NULL;
    CHM_SP rhs, B, Bt, BtB;
    CHM_DN BtBd;
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
internal_betabst_update(double sigma, CHM_FR L, CHM_DN bbsthat,
		       CHM_DN bbstnew)
{
    int i, ppq1 = L->n;
    int ppq = ppq1 - 1;
    double *bx, *hx = (double*)(bbsthat->x),
	*nx = (double*)(bbstnew->x), ans = 0;
    CHM_DN chb;

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
    int *dims = INTEGER(getAttrib(scal, R_DimSymbol)), j,
	n = asInteger(ns), psqr;
    double *scCp, *ansp, *tmp, df = asReal(dfp), one = 1, zero = 0;

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
	internal_symmetrize(ansj, dims[0]);
    }

    PutRNGstate();
    Free(scCp); Free(tmp);
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
 * @param x pointer to an lmer2 object
 * @param savebp pointer to a logical scalar indicating if the
 * random-effects should be saved
 * @param nsampp pointer to an integer scalar of the number of samples
 * to generate
 * @param transp pointer to an logical scalar indicating if the
 * variance components should be transformed.
 *
 * @return a matrix
 */
SEXP lmer2_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbosep, SEXP deviancep)
{
    SEXP ST = GET_SLOT(x, lme4_STSym), Pars = PROTECT(lmer2_getPars(x)), ans;
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
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
    int nrbase = p + 1 + dims[np_POS]; /* rows always included */
    int nrtot = nrbase + dev + (saveb ? q : 0);
    CHM_DN chhat, chnew = M_cholmod_allocate_dense(qpp1, 1, qpp1, CHOLMOD_REAL, &c);
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
	chhat = internal_lmer2_effects(L);
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
    M_cholmod_free_dense(&chnew, &c);
				/* Restore pars, refactor, etc. */
    lmer2_setPars(x, Pars);
    UNPROTECT(2);
    return ans;
}

/**
 * Create and initialize L
 *
 * @param Zt pointer to the random-effects model matrix (transposed)
 *
 * @return L
 */
SEXP lmer2_createL(SEXP Zt)
{
    CHM_SP cZt = AS_CHM_SP(Zt);
    CHM_FR L;
    double one[] = {1, 0};
    int fll = c.final_ll;

    L = M_cholmod_analyze(cZt, &c);
    if (!M_cholmod_factorize_p(cZt, one, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    c.final_ll = fll;

    return M_chm_factor_to_SEXP(L, 1);
}

/**
 * Check validity of an lmer2 object.
 *
 * @param x Pointer to an lmer2 object
 *
 * @return TRUE if the object is a valid lmer2 object, else a string
 * describing the nature of the violation.
 */
SEXP lmer2_validate(SEXP x)
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
    CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym)),
	ZXyt = AS_CHM_SP(GET_SLOT(x, lme4_ZXytSym));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nf = dd[nf_POS], n = dd[n_POS], p = dd[p_POS], q = dd[q_POS];
    int i, nq, ppq1 = p + q + 1;

    R_CheckStack();
				/* check lengths */
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length nf"));
    if (LENGTH(GpP) != (nf + 3))
	return mkString(_("Slot Gp must have length nf + 3"));
    if (Gp[0] != 0 || Gp[nf + 2] != ppq1)
	return mkString(_("Gp[1] != 0 or Gp[nf+3] != p+q+1"));
    if (p && LENGTH(fixefP) != p) /* p == 0 is a special case */
	return mkString(_("Slot fixef must have length p"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length q"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or n"));
    if (LENGTH(offsetP) && LENGTH(offsetP) != n)
	return mkString(_("Slot offset must have length 0 or n"));
    if (LENGTH(devianceP) != (Sdr_POS + 1) ||
	LENGTH(getAttrib(devianceP, R_NamesSymbol)) != (Sdr_POS + 1))
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

    return ScalarLogical(1);
}

/**
 * Remove the names attribute from an object
 *
 * @param x 
 *
 * @return x without the names attribute
 */
static SEXP R_INLINE unname(SEXP x) {
    setAttrib(x, R_NamesSymbol, R_NilValue);
    return x;
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


#define IRLS_MAXITER  60
#define IRLS_TOL      1e-9

/**
 * Evaluate the convergence criterion and copy eta to
 * etaold
 *
 * @param eta current value of the linear predictors
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

static int*
nz_col(int *nz, int j, int nf, const int *Gp,
       const int *nc, const int *zi, const int *zp)
{
    int i, ii, k, nextra, zrow;
    for (ii = zp[j]; ii < zp[j + 1]; ii++) {
	k = Gp_grp(zrow = zi[ii], nf, Gp);
	nextra = (zrow - Gp[k]) % nc[k];
	for (i = 0; i <= nextra; i++) nz[zrow - i] = 1;
    }
    return nz;
}

SEXP nlmer_create_Vt(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	Zt = GET_SLOT(x, lme4_ZtSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    SEXP Zdims = GET_SLOT(Zt, lme4_DimSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)), *adims = INTEGER(Zdims),
	*vi, *vp, *zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	*zp = INTEGER(GET_SLOT(Zt, lme4_pSym)),
	i, j, nf = LENGTH(ST), nnz;
    int *nc = Calloc(nf, int), *nz = Calloc(adims[0], int);
    
    SET_SLOT(ans, lme4_DimSym, duplicate(Zdims));
    SET_SLOT(ans, lme4_DimNamesSym, allocVector(VECSXP, 2));
    for (i = 0; i < nf; i++) 	/* populate nc */
	nc[i] = *INTEGER(getAttrib(VECTOR_ELT(ST, i), R_DimSymbol));
				/* create and evaluate the p slot */
    vp = INTEGER(ALLOC_SLOT(ans, lme4_pSym, INTSXP, adims[1] + 1));
    vp[0] = 0;
    for (j = 0; j < adims[1]; j++) {
	AZERO(nz, adims[0]);
	nz_col(nz, j, nf, Gp, nc, zi, zp);
	for (i = 0, nnz = 0; i < adims[0]; i++)
	    if (nz[i]) nnz++;
	vp[j + 1] = vp[j] + nnz;
    }
    vi = INTEGER(ALLOC_SLOT(ans, lme4_iSym, INTSXP, vp[adims[1]]));
    AZERO(REAL(ALLOC_SLOT(ans, lme4_xSym, REALSXP, vp[adims[1]])), vp[adims[1]]);
				/* fill in the i slot */
    for (j = 0; j < adims[1]; j++) {
	int pos = vp[j];
	AZERO(nz, adims[0]);
	nz_col(nz, j, nf, Gp, nc, zi, zp);
	for (i = 0; i < adims[0]; i++) if (nz[i]) vi[pos++] = i;
    }

    UNPROTECT(1); Free(nc); Free(nz);
    return ans;
}

/**
 * Update the Vt slot in a glmer or nlmer object.
 *
 * @param x pointer to a glmer or nlmer object
 *
 */
SEXP nlmer_update_Vt(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	Vt = GET_SLOT(x, install("Vt")),
	Zt =  GET_SLOT(x, lme4_ZtSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	*vi = INTEGER(GET_SLOT(Vt, lme4_iSym)),
	*vp = INTEGER(GET_SLOT(Vt, lme4_pSym)),
	*zi = INTEGER(GET_SLOT(Zt, lme4_iSym)),
	*zp = INTEGER(GET_SLOT(Zt, lme4_pSym)),
	i, ione = 1, iv, iz, j, mnc, nf = LENGTH(ST),
	ncol = INTEGER(GET_SLOT(Zt, lme4_DimSym))[1];
    int *nc = Calloc(nf, int);
    double **st = Calloc(nf, double*), *tmp,
	*vx = REAL(GET_SLOT(Vt, lme4_xSym)),
	*zx = REAL(GET_SLOT(Zt, lme4_xSym));

    for (i = 0, mnc = 0; i < nf; i++) {	/* populate nc and st */
	SEXP STi = VECTOR_ELT(ST, i);
	nc[i] = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	if (nc[i] > mnc) mnc = nc[i]; /* max num of cols */
	st[i] = REAL(STi);
    }
    tmp = Calloc(mnc, double);
    for (j = 0; j < ncol; j++) {
	int iz2 = zp[j + 1];
				/* premultiply by T' */
	for (iz = zp[j], iv = vp[j]; iz < iz2; iz++) {
	    int k = Gp_grp(zi[iz], nf, Gp);
	    if (nc[k] > 1) {
		int itmp = (zi[iz] - Gp[k]) % nc[k];
		AZERO(tmp, mnc);
		tmp[itmp] = zx[iz];
		for (i = 1; i < nc[k] && (iz + 1) < iz2; i++) {
		    if (zi[iz + 1] != zi[iz] + 1) break;
		    tmp[itmp++] = zx[++iz];
		}
		F77_CALL(dtrmv)("L", "T", "U", &(nc[k]), st[k],
				&(nc[k]), tmp, &ione);
		for (i = 0; i < nc[k] && iv < vp[j + 1]; i++, iv++) {
		    vx[iv] = tmp[i];
		    if (vi[iv + 1] != vi[iv] + 1) break;
		}
	    } else vx[iv++] = zx[iz++];
	}
	for (iv = vp[j]; iv < vp[j + 1]; iv++) {
	    int k = Gp_grp(vi[iv], nf, Gp);
	    vx[iv] *= st[k][((vi[iv] - Gp[k]) % nc[k]) * (nc[k] + 1)];
	}
    }    
    Free(st); Free(nc); Free(tmp);
    return R_NilValue;
}
/* FIXME: Use TS_mult instead */
/**
 * Update the ranef slot, b=TSu, in a glmer or nlmer object.
 *
 * @param x pointer to a glmer or nlmer object
 *
 */
SEXP nlmer_update_ranef(SEXP x)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *Gp = INTEGER(GET_SLOT(x, lme4_GpSym)),
	i, ione = 1, nf = LENGTH(ST);
    double *b = REAL(GET_SLOT(x, lme4_ranefSym)),
	*u = REAL(GET_SLOT(x, install("uvec")));
    
    for (i = 0; i < nf; i++) {
	SEXP STi = VECTOR_ELT(ST, i);
	double *sti = REAL(STi);
	int base = Gp[i], j, k,
	    nci = INTEGER(getAttrib(STi, R_DimSymbol))[0];
	
	for (j = base; j < Gp[i+1]; j += nci) {
	    for (k = 0; k < nci; k++) { /* premultiply  by S */
		int jj = j + k;
		b[jj] = u[jj] * sti[k];
	    }
	    if (nci > 1)	/* premultiply  by T */
		F77_CALL(dtrmv)("L", "N", "U", &nci, sti,
				&nci, &(u[j]), &ione);
	}
    }
    return R_NilValue;
}

/* Nonlinear mixed models */

SEXP nlmer_validate(SEXP x)
{
    SEXP GpP = GET_SLOT(x, lme4_GpSym),
	ST = GET_SLOT(x, lme4_STSym),
 	devianceP = GET_SLOT(x, lme4_devianceSym),
	dimsP = GET_SLOT(x, lme4_dimsSym),
	fixefP = GET_SLOT(x, lme4_fixefSym),
	flistP = GET_SLOT(x, lme4_flistSym),
	ranefP = GET_SLOT(x, lme4_ranefSym),
	weightsP = GET_SLOT(x, lme4_weightsSym) ;
    CHM_SP Xt = AS_CHM_SP(GET_SLOT(x, install("Xt"))),
	Zt = AS_CHM_SP(GET_SLOT(x, install("Zt")));
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    int *Gp = INTEGER(GpP), *dd = INTEGER(dimsP);
    int nf = dd[nf_POS], n = dd[n_POS], p = dd[p_POS], q = dd[q_POS], s = dd[s_POS];

    R_CheckStack();
    if (!LENGTH(devianceP))
	return mkString(_("deviance slot must have positive length"));
    if (nf < 1 || LENGTH(flistP) != nf || LENGTH(ST) != nf)
	return mkString(_("Slots ST, and flist must have length nf"));
    if (LENGTH(GpP) != (nf + 1))
	return mkString(_("Slot Gp must have length nf + 1"));
    if (Gp[0] != 0 || Gp[nf] != q)
	return mkString(_("Gp[1] != 0 or Gp[nf+1] != q"));
    if (LENGTH(ranefP) != q)
	return mkString(_("Slot ranef must have length q"));
    if (LENGTH(fixefP) != p)
	return mkString(_("Slot fixef must have length p"));
    if (LENGTH(weightsP) && LENGTH(weightsP) != n)
	return mkString(_("Slot weights must have length 0 or n"));
    if (Zt->nrow != q || Zt->ncol != n * s)
	return mkString(_("Slot Zt must have dimensions q by n*s"));
    if (Xt->nrow != p || Xt->ncol != n * s)
	return mkString(_("Slot Xt must have dimensions p by n*s"));
    if (L->n != q || !L->is_ll || !L->is_monotonic)
	return mkString(_("Slot L must be a monotonic LL' factorization of size q"));
    
    return ScalarLogical(1);
}

/**
 * Evaluate the nonlinear model, its gradient matrix and the residual
 * sum of squares
 *
 * @param x pointer to an nlmer object
 * @param uform logical indicator of whether to use V'u (as opposed to Z'b)
 *
 * @return the residual sum of squares
 */
static double
internal_nlmer_eval_model(SEXP x, int uform)
{
    SEXP gg, pnames = GET_SLOT(x, install("pnames")),
	rho = GET_SLOT(x, install("env")), vv;
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), *gdims;
    int i, n = dims[n_POS], s = dims[s_POS];
    double *y = REAL(GET_SLOT(x, lme4_ySym)), *mu,
	one[] = {1, 0}, sumsq, zero[] = {0, 0};
    CHM_DN Phi = M_cholmod_allocate_dense(n * s, 1, n * s, CHOLMOD_REAL, &c);
    
    R_CheckStack();
				/* Evaluate Phi */
    if (!(i = M_cholmod_sdmult(AS_CHM_SP(GET_SLOT(x, lme4_XtSym)), 1 /*trans*/, one,
			       zero, AS_CHM_DN(GET_SLOT(x, lme4_fixefSym)), Phi, &c)))
	error(_("cholmod_sdmult returned error code %d"), i);
    i = (uform) ? M_cholmod_sdmult(AS_CHM_SP(GET_SLOT(x, install("Vt"))), 1 /*trans*/,
				   one, one, AS_CHM_DN(GET_SLOT(x, install("uvec"))),
				   Phi, &c)
	: M_cholmod_sdmult(AS_CHM_SP(GET_SLOT(x, lme4_ZtSym)), 1 /*trans*/, one, one,
			   AS_CHM_DN(GET_SLOT(x, lme4_ranefSym)), Phi, &c);
    R_CheckStack();
    if (i) error(_("cholmod_sdmult returned error code %d"), i);

    /* distribute the parameters in the environment */
    for (i = 0; i < dims[s_POS]; i++) {
	vv = findVarInFrame(rho, install(CHAR(STRING_ELT(pnames, i))));
	if (!isReal(vv) || LENGTH(vv) != n)
	    error(_("Parameter %s in the environment must be a length %d numeric vector"),
		  CHAR(STRING_ELT(pnames, i)), n);
	Memcpy(REAL(vv), ((double*)(Phi->x)) + i * n, n);
    }
    vv = PROTECT(eval(GET_SLOT(x, install("model")), rho));
    if (!isReal(vv) || LENGTH(vv) != n)
	error(_("evaluated model is not a numeric vector of length %d"), n);
    gg = getAttrib(vv, install("gradient"));
    if (!isReal(gg) || !isMatrix(gg))
	error(_("gradient attribute of evaluated model must be a numeric matrix"));
    gdims = INTEGER(getAttrib(gg, R_DimSymbol));
    if (gdims[0] != n ||gdims[1] != s)
	error(_("gradient matrix must be of size %d by %d"), n, s);
    SET_SLOT(x, lme4_muSym, vv);
    mu = REAL(vv);
    for (i = 0, sumsq = 0; i < n; i++) {
	double res = y[i] - mu[i];
	sumsq += res * res;
    }
    UNPROTECT(1);
    return sumsq;
}

SEXP nlmer_eval_model(SEXP x, SEXP uform)
{
    return ScalarReal(internal_nlmer_eval_model(x, asLogical(uform)));
}


/**
 * Create an nlmer object
 *
 * @param env environment in which to evaluate the model
 * @param model nonlinear model expression as a call
 * @param frame model frame
 * @param pnames character vector of parameter names
 * @param call matched call to the R nlmer function
 * @param flist factor list
 * @param Xt transpose of fixed-effects model matrix
 * @param Zt transpose of random-effects model matrix
 * @param y response vector
 * @param weights weights vector (may have length 0)
 * @param cnames list of column names
 * @param Gp integer vector of group pointers
 * @param fixef numeric vector of fixed effects
 */
SEXP nlmer_create(SEXP env, SEXP model, SEXP frame, SEXP pnames,
		  SEXP call, SEXP flist, SEXP Xt, SEXP Zt, SEXP y,
		  SEXP weights, SEXP cnames, SEXP Gp, SEXP fixef)
{
    SEXP ST, ans = PROTECT(NEW_OBJECT(MAKE_CLASS("nlmer")));
    char *DEVIANCE_NAMES[]={"ML","REML","ldZ","ldX","lr2", "bqd", "Sdr", ""};
    char *DIMS_NAMES[]={"nf","n","p","q", "s", "np","isREML","famType","isNested",""};
    int *Gpp = INTEGER(Gp), *Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)), *dims, i, iT;
    CHM_SP cVt, ts1, ts2;
    CHM_FR L;
    double one[] = {1,0};

    SET_SLOT(ans, install("env"), env);
    SET_SLOT(ans, install("model"), model);
    SET_SLOT(ans, install("frame"), frame);
    SET_SLOT(ans, install("pnames"), pnames);
    SET_SLOT(ans, install("call"), call);
    SET_SLOT(ans, lme4_flistSym, flist);
    SET_SLOT(ans, install("Xt"), Xt);
    SET_SLOT(ans, lme4_ZtSym, Zt);
    SET_SLOT(ans, lme4_ySym, y);
    SET_SLOT(ans, lme4_weightsSym, weights);
    SET_SLOT(ans, lme4_cnamesSym, cnames);
    SET_SLOT(ans, lme4_GpSym, Gp);
    SET_SLOT(ans, lme4_fixefSym, fixef);
    SET_SLOT(ans, lme4_dimsSym,
	     internal_make_named(INTSXP, DIMS_NAMES));
    dims = INTEGER(GET_SLOT(ans, lme4_dimsSym));
    dims[nf_POS] = LENGTH(flist);
    dims[n_POS] = LENGTH(y);
    dims[np_POS] = 0;
    dims[p_POS] = LENGTH(fixef);
    dims[q_POS] = Zdims[0];
    dims[s_POS] = Zdims[1]/dims[n_POS];
    dims[isREML_POS] = FALSE;
    dims[famType_POS] = -1;
    dims[isNest_POS] = TRUE;
    SET_SLOT(ans, lme4_ranefSym, allocVector(REALSXP, Zdims[0]));
    AZERO(REAL(GET_SLOT(ans, lme4_ranefSym)), Zdims[0]);
    SET_SLOT(ans, install("uvec"), allocVector(REALSXP, Zdims[0]));
    AZERO(REAL(GET_SLOT(ans, install("uvec"))), Zdims[0]);
    internal_nlmer_eval_model(ans, 0); /* check the model evaluation */

    SET_SLOT(ans, lme4_devianceSym,
	     internal_make_named(REALSXP, DEVIANCE_NAMES));
    AZERO(REAL(GET_SLOT(ans, lme4_devianceSym)), 7);

    SET_SLOT(ans, lme4_STSym, allocVector(VECSXP, dims[nf_POS]));
    ST = GET_SLOT(ans, lme4_STSym);
    iT = TRUE;			/* is T the identity? */
    for (i = 0; i < dims[nf_POS]; i++) {
	int nci = (Gpp[i + 1] - Gpp[i])/
	    LENGTH(getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol));
	SET_VECTOR_ELT(ST, i, allocMatrix(REALSXP, nci, nci));
	dims[np_POS] += (nci*(nci + 1))/2;
	if (nci > 1) iT = FALSE;
    }
    ST_initialize(ST, Gp, Zt);	/* initialize ST */
    SET_SLOT(ans, install("Vt"), iT ? duplicate(Zt) : nlmer_create_Vt(ans));
    nlmer_update_Vt(ans);
    cVt = AS_CHM_SP(GET_SLOT(ans, install("Vt")));
    
    /* Create Mt beginning with s identity matrices concatenated horizontally */
    ts1 = M_cholmod_allocate_sparse((size_t) dims[n_POS],
				    (size_t) Zdims[1], (size_t) Zdims[1],
				    1/*sorted*/, 1/*packed*/,
				    0/*stype*/, CHOLMOD_REAL, &c);
    for (i = 0; i < Zdims[1]; i++) {
	((int*)(ts1->p))[i] = i;
	((int*)(ts1->i))[i] = i % dims[n_POS];
	((double*)(ts1->x))[i] = 1;
    }
    ((int*)(ts1->p))[Zdims[1]] = Zdims[1];
    ts2 = M_cholmod_transpose(ts1, TRUE/*values*/, &c);
    M_cholmod_free_sparse(&ts1, &c);
    ts1 = M_cholmod_ssmult(cVt, ts2, /* Create pattern for Mt */
			   0 /*stype*/, 1 /*values*/, 1 /*sorted*/, &c);
    M_cholmod_free_sparse(&ts2, &c);
    SET_SLOT(ans, install("Mt"),
	     M_chm_sparse_to_SEXP(ts1, 0, 0, 0, "", R_NilValue));
/* FIXME: Check for nesting in Mt here? */
    i = c.final_ll;
    c.final_ll = 1;
    L = M_cholmod_analyze(ts1, &c); /* Create pattern for L */
    if (!M_cholmod_factorize_p(ts1, one, (int*) NULL, 0, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    if (!M_cholmod_change_factor(CHOLMOD_REAL, 1 /* to_ll */,
				 L->is_super, 1 /* packed */,
				 1 /* monotonic */, L, &c))
	error(_("cholmod_change_factor failed"));
    c.final_ll = i;
    SET_SLOT(ans, lme4_LSym, M_chm_factor_to_SEXP(L, 0));
    M_cholmod_free_factor(&L, &c);

    M_cholmod_free_sparse(&ts1, &c);
    nlmer_update_Vt(ans);

    UNPROTECT(1);
    return ans;
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
    SEXP Mt = GET_SLOT(x, install("Mt")),
	Vt = GET_SLOT(x, install("Vt"));
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)),
	*mi = INTEGER(GET_SLOT(Mt, lme4_iSym)),
	*mp = INTEGER(GET_SLOT(Mt, lme4_pSym)),
	*vi = INTEGER(GET_SLOT(Vt, lme4_iSym)),
	*vp = INTEGER(GET_SLOT(Vt, lme4_pSym)), jv;
    double *grad = REAL(getAttrib(GET_SLOT(x, lme4_muSym),
				  install("gradient"))),
	*mx = REAL(GET_SLOT(Mt, lme4_xSym)),
	*vx = REAL(GET_SLOT(Vt, lme4_xSym));
    int n = dims[n_POS], s = dims[s_POS];

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
    return R_NilValue;
}

/**
 * Update the working residual as y - mu + M u
 *
 * @param x pointer to an nlmer object
 * @param wrkres array to hold the updated working residual
 *
 * @return wrkres pointer
 */
static
double *internal_nlmer_update_wrkres(SEXP x, double *wrkres)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i;
    CHM_SP Mt = AS_CHM_SP(GET_SLOT(x, install("Mt")));
    CHM_DN cwrk = M_numeric_as_chm_dense(Alloca(1, cholmod_dense),
					 wrkres, dims[n_POS]),
	cu = AS_CHM_DN(GET_SLOT(x, install("uvec")));
    double *mu = REAL(GET_SLOT(x, lme4_muSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym)), one[] = {1,0};

    R_CheckStack();
    Memcpy(wrkres, y, dims[n_POS]);
    if (!(i = M_cholmod_sdmult(Mt, 1 /* trans */, one, one, cu, cwrk, &c)))
	error(_("cholmod_sdmult returned error code %d"), i);
    for (i = 0; i < dims[n_POS]; i++) wrkres[i] -= mu[i];
    return wrkres;
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

    internal_nlmer_update_wrkres(x, REAL(ans));
    UNPROTECT(1);
    return ans;
}

/**
 * Iterate to determine the conditional modes of the random effects.
 *
 * @param x pointer to an nlmer object
 *
 * @return An indicator of whether the iterations converged
 */
static int internal_nbhat(SEXP x)
{
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, j;
    int n = dims[n_POS], q = dims[q_POS], zq[] = {0, dims[q_POS]};
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    CHM_DN cz = M_cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, &c),
	cu, cMtz = M_cholmod_allocate_dense(q, 1, q, CHOLMOD_REAL, &c);
    CHM_SP Mt = AS_CHM_SP(GET_SLOT(x, install("Mt")));
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym)),
	*u = REAL(GET_SLOT(x, install("uvec"))),
	*uold = Calloc(q, double), *z = ((double*)(cz->x)),
	crit = IRLS_TOL + 1, dn = (double)n, one[] = {1,0}, zero[] = {0,0};

    nlmer_update_Vt(x);
    Memcpy(uold, u, q);
    for (i = 0; i < IRLS_MAXITER && crit > IRLS_TOL; i++) {
	dev[Sdr_POS] = internal_nlmer_eval_model(x, 1);
	for (j = 0, dev[bqd_POS] = 0; j < q; j++) dev[bqd_POS] += u[j] * u[j];
#ifdef DEBUG_NLMER
	Rprintf("%3d: %20.15g %20.15g %20.15g\n", i, dev[Sdr_POS], dev[bqd_POS],
		dev[Sdr_POS] + dev[bqd_POS]);
#endif
	nlmer_update_Mt(x);
	j = c.final_ll;
	c.final_ll = L->is_ll;
	if (!M_cholmod_factorize_p(Mt, one, (int*) NULL, 0, L, &c))
	    error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
		  c.status, L->minor, L->n);
	c.final_ll = j;
	internal_nlmer_update_wrkres(x, z);
	if (!(j = M_cholmod_sdmult(Mt, 0 /* trans */, one, zero, cz, cMtz, &c)))
	    error(_("cholmod_sdmult returned error code %d"), j);
	if (!(cu = M_cholmod_solve(CHOLMOD_A, L, cMtz, &c)))
	    error(_("cholmod_solve (CHOLMOD_A) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
	Memcpy(u, (double*)(cu->x), q);
	M_cholmod_free_dense(&cu, &c); 
/* FIXME: replace this by an orthogonality convergence criterion */
 	crit = conv_crit(uold, u, q);
    }
    dev[lr2_POS] = log(dev[Sdr_POS] + dev[bqd_POS]);
    chm_log_abs_det2(&(dev[ldZ_POS]), 1, zq, L);
    dev[ML_POS] = dev[REML_POS] =
	dev[ldZ_POS] + dn * (1. + dev[lr2_POS] + log(2. * PI / dn));

    Free(uold);
    M_cholmod_free_dense(&cu, &c); M_cholmod_free_dense(&cMtz, &c);
    return (crit > IRLS_TOL) ? 0 : i;
}

SEXP nlmer_bhat(SEXP x) {
    return ScalarInteger(internal_nbhat(x));
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
 * evaluated once only. 
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

static int flType(SEXP famName)
{
    const char *fam = CHAR(STRING_ELT(famName, 0)),
	*lnk = CHAR(STRING_ELT(famName, 1));

    if ((!strcmp(fam, "gaussian")) && (!strcmp(lnk, "identity")))
	return -1;
    if ((!strcmp(fam, "binomial")) && (!strcmp(lnk, "logit")))
	return 1;
    if ((!strcmp(fam, "binomial")) && (!strcmp(lnk, "probit")))
	return 2;
    if ((!strcmp(fam, "poisson")) && (!strcmp(lnk, "log")))
	return 3;
    return 0;
}

/**
 * Evaluate the inverse link function at eta storing the result in mu
 *
 * @param x pointer to a glmer2 object
 */
SEXP glmer_linkinv(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[famType_POS];
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
 * @param x pointer to a glmer2 object
 * @param var pointer to positions to hold computed values
 *
 * @return var
 */
static double *glmer_var(SEXP x, double *var)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[famType_POS];
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
 * @param x pointer to a glmer2 object
 * @param dmu_deta pointer to positions to hold computed values
 *
 * @return dmu_deta
 */
static double *glmer_dmu_deta(SEXP x, double *dmu_deta)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS], fltype = dims[famType_POS];
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
 * @param x pointer to a glmer2 object
 * @param dev_res pointer to an area to hold the result
 *
 * @return sum of the deviance residuals
 */
SEXP glmer_dev_resids(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, fltype = dims[famType_POS], n = dims[n_POS];
    double *dev_res = REAL(findVarInFrame(rho, install("devResid"))),
	*mu = REAL(GET_SLOT(x, lme4_muSym)),
	*wts = REAL(findVarInFrame(rho, install("pwts"))),
	*y = REAL(findVarInFrame(rho, lme4_ySym)), sum;

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
    for (i = 0, sum = 0; i < n; i++) sum += dev_res[i];
    return R_NilValue;
}

/**
 * Evaluate the linear predictor as model offset + X \beta + V u
 *
 * @param x pointer to a glmer2 object
 *
 * @return R_NilValue
 */
SEXP glmer_eta(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    SEXP etap = findVarInFrame(rho, lme4_etaSym),
	moff = findVarInFrame(rho, install("offset"));
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int ione = 1, n = dims[n_POS], p = dims[p_POS];
    double *eta = REAL(etap), one[] = {1,0};
    CHM_SP cVt = AS_CHM_SP(GET_SLOT(x, install("Vt")));
    CHM_DN ceta = AS_CHM_DN(etap),
	cu = AS_CHM_DN(GET_SLOT(x, install("uvec")));

    if (moff != R_NilValue)	/* initialize eta to offset */
	Memcpy(eta, REAL(moff), n);
    else			/* initialize eta to zero */
	AZERO(eta, n);
				/* eta := eta + X \beta */
    F77_CALL(dgemv)("N", &n, &p, one, REAL(GET_SLOT(x, lme4_XSym)), &n,
		    REAL(GET_SLOT(x, lme4_fixefSym)), &ione, one, eta, &ione);
				/* eta := eta + V u */
    if (!M_cholmod_sdmult(cVt, 1 /* trans */, one, one, cu, ceta, &c))
	error(_("cholmod_sdmult error returned"));
    return R_NilValue;
}

/**
 * Evaluate new weights and working residuals.
 *
 * @param x pointer to a glmer2 object
 */
SEXP glmer_reweight(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    SEXP moff = findVarInFrame(rho, lme4_offsetSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    int i, n = dims[n_POS];
    double
	*dmu_deta = Calloc(n, double),
	*eta = REAL(findVarInFrame(rho, lme4_etaSym)),
	*mo = (moff == R_NilValue) ? (double*) NULL : REAL(moff),
	*mu = REAL(findVarInFrame(rho, lme4_muSym)),
	*var = Calloc(n, double),
	*w = REAL(findVarInFrame(rho, lme4_weightsSym)),
	*y = REAL(GET_SLOT(x, lme4_ySym)),
	*z = REAL(findVarInFrame(rho, install("z")));

				/* initialize weights to prior wts */
    Memcpy(w, REAL(findVarInFrame(rho, install("pwts"))), n); 
    glmer_linkinv(x);		/* evaluate mu */
    glmer_dmu_deta(x, dmu_deta);
    glmer_var(x, var);
    for (i = 0; i < n; i++) {
	w[i] *= dmu_deta[i] * dmu_deta[i]/var[i];
	/* adjusted working variate */
 	z[i] = (mo ? mo[i] : 0) - eta[i] - (y[i] - mu[i])/dmu_deta[i];
    }
    Free(dmu_deta); Free(var);
    return R_NilValue;
}

/**
 * Create a glmer object
 *
 * @param env environment in which to evaluate the linkinv, dmu_deta, etc. functions
 * @param frame model frame
 * @param famName name of the generalized linear model family and link
 * @param call matched call to the R nlmer function
 * @param flist factor list
 * @param X fixed-effects model matrix
 * @param Zt transpose of random-effects model matrix
 * @param y response vector
 * @param cnames list of column names
 * @param Gp integer vector of group pointers
 * @param fixef numeric vector of fixed effects
 */
SEXP glmer_create(SEXP env, SEXP frame, SEXP famName, SEXP call,
		  SEXP flist, SEXP X, SEXP Zt, SEXP y, SEXP cnames,
		  SEXP Gp, SEXP fixef)
{
    SEXP ST, ans = PROTECT(NEW_OBJECT(MAKE_CLASS("glmer2")));
    char *DEVIANCE_NAMES[]={"ML","REML","ldZ","ldX","lr2", "bqd", "Sdr", ""};
    char *DIMS_NAMES[]={"nf","n","p","q", "s", "np","isREML","famType","isNested",""};
    int *Gpp = INTEGER(Gp), *Zdims = INTEGER(GET_SLOT(Zt, lme4_DimSym)), *dims, i, iT;
    CHM_SP cVt;
    CHM_FR L;
    double one[] = {1,0};

    SET_SLOT(ans, install("env"), env);
    SET_SLOT(ans, install("famName"), famName);
    SET_SLOT(ans, install("frame"), frame);
    SET_SLOT(ans, install("call"), call);
    SET_SLOT(ans, lme4_flistSym, flist);
    SET_SLOT(ans, lme4_XSym, X);
    SET_SLOT(ans, lme4_ZtSym, Zt);
    SET_SLOT(ans, lme4_ySym, y);
    SET_SLOT(ans, lme4_cnamesSym, cnames);
    SET_SLOT(ans, lme4_GpSym, Gp);
    SET_SLOT(ans, lme4_fixefSym, fixef);
    SET_SLOT(ans, lme4_dimsSym,
	     internal_make_named(INTSXP, DIMS_NAMES));
    dims = INTEGER(GET_SLOT(ans, lme4_dimsSym));
    dims[nf_POS] = LENGTH(flist);
    dims[n_POS] = LENGTH(y);
    dims[np_POS] = 0;
    dims[p_POS] = LENGTH(fixef);
    dims[q_POS] = Zdims[0];
    dims[s_POS] = 0;
    dims[isREML_POS] = FALSE;
    dims[famType_POS] = flType(famName);
    dims[isNest_POS] = TRUE;
    SET_SLOT(ans, lme4_ranefSym, allocVector(REALSXP, Zdims[0]));
    AZERO(REAL(GET_SLOT(ans, lme4_ranefSym)), Zdims[0]);
    SET_SLOT(ans, install("uvec"), allocVector(REALSXP, Zdims[0]));
    AZERO(REAL(GET_SLOT(ans, install("uvec"))), Zdims[0]);
    SET_SLOT(ans, lme4_devianceSym,
	     internal_make_named(REALSXP, DEVIANCE_NAMES));
    AZERO(REAL(GET_SLOT(ans, lme4_devianceSym)), 7);

    SET_SLOT(ans, lme4_STSym, allocVector(VECSXP, dims[nf_POS]));
    ST = GET_SLOT(ans, lme4_STSym);
    iT = TRUE;			/* is T the identity? */
    for (i = 0; i < dims[nf_POS]; i++) {
	int nci = (Gpp[i + 1] - Gpp[i])/
	    LENGTH(getAttrib(VECTOR_ELT(flist, i), R_LevelsSymbol));
	SET_VECTOR_ELT(ST, i, allocMatrix(REALSXP, nci, nci));
	dims[np_POS] += (nci*(nci + 1))/2;
	if (nci > 1) iT = FALSE;
    }
    ST_initialize(ST, Gp, Zt); /* initialize ST */
    SET_SLOT(ans, install("Vt"), iT ? duplicate(Zt) : nlmer_create_Vt(ans));
    nlmer_update_Vt(ans);
    cVt = AS_CHM_SP(GET_SLOT(ans, install("Vt")));

    i = c.final_ll;
    c.final_ll = 1;
    L = M_cholmod_analyze(cVt, &c); /* Create pattern for L */
    if (!M_cholmod_factorize_p(cVt, one, (int*) NULL, 0, L, &c))
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    if (!M_cholmod_change_factor(CHOLMOD_REAL, 1 /* to_ll */,
				 L->is_super, 1 /* packed */,
				 1 /* monotonic */, L, &c))
	error(_("cholmod_change_factor failed"));
    c.final_ll = i;
    SET_SLOT(ans, lme4_LSym, M_chm_factor_to_SEXP(L, 0));
    M_cholmod_free_factor(&L, &c);

    nlmer_update_Vt(ans);
    UNPROTECT(1);
    return ans;
}

/**
 * Factor (Vt WW (Vt)' + I) as L, solve
 * (LL')u = (Vt W z) where W is diagonal
 *
 * @param Vt - sparse representation of V'
 * @param w - diagonal of matrix W
 * @param z - working residual vector
 * @param L - factor to update
 * @param u - holds the solution
 *
 * @return squared length of u
 */
static double
internal_glmer_update_L(SEXP x)
{
    SEXP rho = GET_SLOT(x, lme4_envSym);
    CHM_SP Vt = AS_CHM_SP(GET_SLOT(x, install("Vt")));
    CHM_SP wV = M_cholmod_copy_sparse(Vt, &c);
    CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
    double *u = REAL(GET_SLOT(x, install("uvec"))), *w, *z;
    int *wi = (int*)(wV->i), *wp = (int*)(wV->p),
	i, j, m = wV->nrow, n = wV->ncol;
    CHM_DN td,
	rhs = M_cholmod_allocate_dense(m, 1, m, CHOLMOD_REAL, &c);
    double *wx = (double*)(wV->x), *rh = (double*)(rhs->x),
	one[] = {1, 0}, val;

    /* reweight the model then extract the pointers (which may have changed) */
    glmer_reweight(x);
    w = REAL(findVarInFrame(rho, lme4_weightsSym));
    z = REAL(findVarInFrame(rho, install("z")));

    AZERO(rh, m);
    for (j = 0; j < n; j++) {
	for (i = wp[j]; i < wp[j + 1]; i++) {
	    wx[i] *= w[j];	       /* weight jth column */
	    rh[wi[i]] += wx[i] * z[j]; /* inner product */
	}
    }
    if (!M_cholmod_factorize_p(wV, one, (int*) NULL, (size_t) 0, L, &c)) { 
	error(_("cholmod_factorize_p failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    }
    M_cholmod_free_sparse(&wV, &c);
				/* CHOLMOD_A also applies any permutation */
    td = M_cholmod_solve(CHOLMOD_A, L, rhs, &c);
    M_cholmod_free_dense(&rhs, &c);
    for (i = 0, val = 0; i < m; i++) {
	double tt = ((double*)(td->x))[i];
	u[i] = tt;
	val += tt * tt;
    }
    M_cholmod_free_dense(&td, &c);
    return val;
}

SEXP glmer_update_L(SEXP x)
{
    return ScalarReal(internal_glmer_update_L(x));
}

static void internal_gbhat(SEXP x)
{
    internal_glmer_update_L(x);
}

/* Functions common to all three forms of mixed models */

/**
 * Set the parameters in an lmer2 object and evaluate the deviance of
 * an lmm or the Laplace approximation to the deviance of a glmm.
 *
 * @param x pointer to a glmer2 object
 * @param xv vector of parameter values
 * @param mtype model type: 0 -> lmm, 1 -> nlmm, 2 -> glmm
 *
 * @return deviance
 */
static double
update_deviance(SEXP x, const double *xv, int mtype)
{
    SEXP ST = GET_SLOT(x, lme4_STSym),
	fixefp = GET_SLOT(x, lme4_fixefSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym));
    double *dev = REAL(GET_SLOT(x, lme4_devianceSym));

    internal_lmer2_setPars(xv, ST); /* common parameters */
    switch(mtype) {
    case 0: {			  /* linear mixed model */
	CHM_SP A = AS_CHM_SP(GET_SLOT(x, lme4_ASym));
	CHM_FR L = AS_CHM_FR(GET_SLOT(x, lme4_LSym));
	internal_update_L(dev, dims, INTEGER(GET_SLOT(x, lme4_GpSym)),
			  ST, A, L);
	break;
    }
    case 1: {			  /* nonlinear mixed model */
	Memcpy(REAL(fixefp), xv + dims[np_POS], dims[p_POS]);
	internal_nbhat(x);
	break;
    }
    case 2: {
	Memcpy(REAL(fixefp), xv + dims[np_POS], dims[p_POS]);
	internal_gbhat(x);
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
 * @param x pointer to an lmer2, glmer2 or nlmer object
 * @param verbp indicator of verbose output
 * @param mtypep model type 0 -> lmer, 1 -> nlmer, 2 -> glmer
 *
 * @return indicator of convergence from nlminb_iterate
 */
SEXP 
lmer2_optimize(SEXP x, SEXP verbp, SEXP mtypep)
{
    SEXP ST = GET_SLOT(x, lme4_STSym);
    int *dims = INTEGER(GET_SLOT(x, lme4_dimsSym)), i, j,
	mtype = asInteger(mtypep), pos, verb = asLogical(verbp);
    int nf = dims[nf_POS], nv = dims[np_POS] + (mtype?dims[p_POS]:0);
    int liv = S_iv_length(OPT, nv), lv = S_v_length(OPT, nv);
    int *iv = Calloc(liv, int);
    double *b = Calloc(2 * nv, double), *d = Calloc(nv, double),
	*g = (double*)NULL, *h = (double*)NULL,
	*v = Calloc(lv, double),
	*xv = internal_lmer2_getPars(ST, Calloc(nv, double)),
	fx = R_PosInf;

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
    Free(iv); Free(v); Free(xv); Free(b); Free(d);
    return ScalarInteger(i);
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
#endif
