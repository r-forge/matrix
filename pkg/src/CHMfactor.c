				/* CHOLMOD factors */
#include "CHMfactor.h"

SEXP CHMfactor_to_sparse(SEXP x)
{
    CHM_FR L = AS_CHM_FR(x), Lcp;
    CHM_SP Lm;
    R_CheckStack();

    /* cholmod_factor_to_sparse changes its first argument. Make a copy */
    Lcp = cholmod_copy_factor(L, &c);
    if (!(Lcp->is_ll))
	if (!cholmod_change_factor(Lcp->xtype, 1, 0, 1, 1, Lcp, &c))
	    error(_("cholmod_change_factor failed with status %d"), c.status);
    Lm = cholmod_factor_to_sparse(Lcp, &c); cholmod_free_factor(&Lcp, &c);
    return chm_sparse_to_SEXP(Lm, 1/*do_free*/, -1/*uploT*/, 0/*Rkind*/,
			      "N"/*non_unit*/, R_NilValue/*dimNames*/);
}

SEXP CHMfactor_solve(SEXP a, SEXP b, SEXP system)
{
    CHM_FR L = AS_CHM_FR(a);
    SEXP bb = PROTECT(dup_mMatrix_as_dgeMatrix(b));
    CHM_DN B = AS_CHM_DN(bb), X;
    int sys = asInteger(system);
    R_CheckStack();

    if (!(sys--))		/* -- align with CHOLMOD defs */
	error(_("system argument is not valid"));

    X = cholmod_solve(sys, L, B, &c);
    UNPROTECT(1);
    return chm_dense_to_SEXP(X, 1/*do_free*/, 0/*Rkind*/,
			     GET_SLOT(bb, Matrix_DimNamesSym));
}

SEXP CHMfactor_spsolve(SEXP a, SEXP b, SEXP system)
{
    CHM_FR L = AS_CHM_FR(a);
    CHM_SP B = AS_CHM_SP(b);
    int sys = asInteger(system);
    R_CheckStack();

    if (!(sys--))		/* -- align with CHOLMOD defs */
	error(_("system argument is not valid"));
    return chm_sparse_to_SEXP(cholmod_spsolve(sys, L, B, &c),
			      1/*do_free*/, 0/*uploT*/, 0/*Rkind*/,
			      "", GET_SLOT(b, Matrix_DimNamesSym));
}

/**
 * Evaluate the logarithm of the square of the determinant of L
 *
 * @param x pointer to a CHMfactor object
 *
 * @return log(det(L)^2)
 *
 */
double chm_factor_ldetL2(CHM_FR f)
{
    int i, j, p;
    double ans = 0;

    if (f->is_super) {
	int *lpi = (int*)(f->pi), *lsup = (int*)(f->super);
	for (i = 0; i < f->nsuper; i++) { /* supernodal block i */
	    int nrp1 = 1 + lpi[i + 1] - lpi[i],
		nc = lsup[i + 1] - lsup[i];
	    double *x = (double*)(f->x) + ((int*)(f->px))[i];

	    for (j = 0; j < nc; j++) {
		ans += 2 * log(fabs(x[j * nrp1]));
	    }
	}
    } else {
	int *li = (int*)(f->i), *lp = (int*)(f->p);
	double *lx = (double *)(f->x);
	
	for (j = 0; j < f->n; j++) {
	    for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++) {};
	    if (li[p] != j) {
		error(_("%d diagonal element of Cholesky factor is missing"), j);
		break;		/* -Wall */
	    }
	    ans += log(lx[p] * ((f->is_ll) ? lx[p] : 1.));
	}
    }
    return ans;
}

SEXP CHMfactor_ldetL2(SEXP x)
{
    CHM_FR L = AS_CHM_FR(x); R_CheckStack();

    return ScalarReal(chm_factor_ldetL2(L));
}
