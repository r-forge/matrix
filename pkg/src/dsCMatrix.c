#include "dsCMatrix.h"

SEXP dsCMatrix_chol(SEXP x, SEXP pivot)
{
    CHM_FR N = AS_CHM_FR(dsCMatrix_Cholesky(x, pivot, ScalarLogical(FALSE),
					    ScalarLogical(FALSE),
					    ScalarReal(0.)));
    /* Must use a copy; cholmod_factor_to_sparse modifies first arg. */
    CHM_FR Ncp = cholmod_copy_factor(N, &c);
    CHM_SP L = cholmod_factor_to_sparse(Ncp, &c);
    CHM_SP R = cholmod_transpose(L, /*values*/ 1, &c);
    SEXP ans = PROTECT(chm_sparse_to_SEXP(R, 1/*do_free*/, 1/*uploT*/, 0/*Rkind*/,
					  "N"/*diag*/, GET_SLOT(x, Matrix_DimNamesSym)));
    R_CheckStack();

    cholmod_free_factor(&Ncp, &c);
    cholmod_free_sparse(&L, &c);
    if (asLogical(pivot)) {
	SEXP piv = PROTECT(allocVector(INTSXP, N->n));
	int *dest = INTEGER(piv), *src = (int*)N->Perm, i;

	for (i = 0; i < N->n; i++) dest[i] = src[i] + 1;
	setAttrib(ans, install("pivot"), piv);
	/* FIXME: Because of the cholmod_factor -> S4 obj ->
	 * cholmod_factor conversions, the value of N->minor will
	 * always be N->n.  Change as_cholmod_factor and
	 * chm_factor_to_SEXP to keep track of Minor.
	 */
	setAttrib(ans, install("rank"), ScalarInteger((size_t) N->minor));
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}

SEXP dsCMatrix_Cholesky(SEXP Ap, SEXP permP, SEXP LDLp, SEXP superP, SEXP ImultP)
{
    char fname[12] = "spdCholesky"; /* template for factorization name */
    /* S|s : super or not
     * P|p : permuted or not
     * D|d :  LDL' or not (= LL')
     */
    int perm = asLogical(permP), LDL = asLogical(LDLp), super = asLogical(superP);
    double mult = asReal(ImultP);
    SEXP Chol;
    CHM_SP A;
    CHM_FR L;
    int sup, ll;

    if (super) fname[0] = 'S';
    if (perm) fname[1] = 'P';
    if (LDL) fname[2] = 'D';
    Chol = get_factors(Ap, fname);
    if (Chol != R_NilValue) return Chol; /* return a cached factor */

    A = AS_CHM_SP(Ap);
    R_CheckStack();
    if (!A->stype)
	error("Non-symmetric matrix passed to dsCMatrix_chol");

    sup = c.supernodal;
    ll = c.final_ll;

    c.final_ll = !LDL;	/* leave as LL' or form LDL' */
    c.supernodal = super ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL;

    if (perm) {
	L = cholmod_analyze(A, &c); /* get fill-reducing permutation */
    } else {			    /* require identity permutation */
	int nmethods = c.nmethods, ord0 = c.method[0].ordering,
	    postorder = c.postorder;
	c.nmethods = 1;
	c.method[0].ordering = CHOLMOD_NATURAL;
	c.postorder = FALSE;
	L = cholmod_analyze(A, &c);
	c.nmethods = nmethods; c.method[0].ordering = ord0;
	c.postorder = postorder;
    }
    if (!cholmod_factorize_p(A, &mult, (int*)NULL, 0 /*fsize*/, L, &c))
	error(_("Cholesky factorization failed"));
    c.supernodal = sup;	/* restore previous setting */
    c.final_ll = ll;
    Chol = set_factors(Ap, chm_factor_to_SEXP(L, 1), fname);
    return Chol;
}

/**
 * Fast version of getting at the diagonal matrix D of the
 * (generalized) simplicial Cholesky LDL' decomposition of a
 * (sparse symmetric) dsCMatrix.
 *
 * @param Ap  symmetric CsparseMatrix
 * @param permp  logical indicating if permutation is allowed
 *
 * @return SEXP containing either the vector diagonal entries of D,
 *         or just  sum_i D[i], prod_i D[i] or  sum_i log(D[i]).
 */
/* Started as copy + modification of dsCMatrix_Cholesky();
 * builds strongly on diag_tC(...., SEXP resultKind) in Csparse.c
*/
SEXP dsCMatrix_LDL_D(SEXP Ap, SEXP permP, SEXP resultKind)
{
    char fname[12] = "spDCholesky"; /* template for factorization name */
    /* S|s : super or not
     * P|p : permuted or not
     * D|d :  LDL' or not (= LL')
     */
    const int perm = asLogical(permP), LDL = 1, super = 0;
    SEXP Chol;
    CHM_SP A;
    CHM_FR L;
    int sup, ll;

    /* if (super) fname[0] = 'S'; */
    if (perm) fname[1] = 'P';
    /* if (LDL) fname[2] = 'D'; */
    Chol = get_factors(Ap, fname);
    if (Chol != R_NilValue) { /* use the *cached* Cholesky() factor */
	return diag_tC(GET_SLOT(Chol, Matrix_pSym),
		       GET_SLOT(Chol, Matrix_xSym),
		       GET_SLOT(Chol, Matrix_permSym),
		       resultKind);
    }
    else {
	A = AS_CHM_SP(Ap);
	R_CheckStack();
	if (!A->stype)
	    error("Non-symmetric matrix passed to dsCMatrix_LDL_D");

	sup = c.supernodal;
	ll = c.final_ll;

	c.final_ll = !LDL;	/* leave as LL' or form LDL' */
	c.supernodal = super ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL;

	if (perm) {
	    L = cholmod_analyze(A, &c); /* get fill-reducing permutation */
	} else {			/* require identity permutation */
	    int nmethods = c.nmethods,
		ord0 = c.method[0].ordering, postorder = c.postorder;
	    c.nmethods = 1;
	    c.method[0].ordering = CHOLMOD_NATURAL; c.postorder = FALSE;
	    L = cholmod_analyze(A, &c);
	    c.nmethods = nmethods;
	    c.method[0].ordering = ord0;            c.postorder = postorder;
	}
	if (!cholmod_factorize(A, L, &c))
	    error(_("Cholesky factorization failed"));
	/* restore previous setting */
	c.supernodal = sup;
	c.final_ll = ll;

	/* Now, instead of
	 *   Chol = set_factors(Ap, chm_factor_to_SEXP(L, 1), fname);
	 *   return Chol;
	 *
	 * get the correct entries from the CHOLMOD 'L' struct directly : */

	return diag_tC_ptr(L->n,
			   L->p,
			   L->x,
			   L->Perm,
			   resultKind);
    }
    return R_NilValue;/*just for now --- FIXME */
}

static
SEXP get_factor_pattern(SEXP obj, char *pat, int offset)
{
    SEXP facs = GET_SLOT(obj, Matrix_factorSym), nms;
    int i;

    if (!LENGTH(facs)) /* Fast return for empty factor slot: */
	return R_NilValue;
    nms = getAttrib(facs, R_NamesSymbol);
    for (i = 0; i < LENGTH(nms); i++) {
	const char *nm = CHAR(STRING_ELT(nms, i));
	if (strlen(nm) > offset && !strcmp(pat + offset, nm + offset))
	    return VECTOR_ELT(facs, i);
    }
    return R_NilValue;
}

SEXP dsCMatrix_Csparse_solve(SEXP a, SEXP b)
{
    SEXP Chol = get_factor_pattern(a, "...Cholesky", 3);
    CHM_FR L;
    CHM_SP cx, cb = AS_CHM_SP(b);
    R_CheckStack();

    if (Chol == R_NilValue) /* compute (and cache) "sPDCholesky" */
	Chol = dsCMatrix_Cholesky(a,
				  ScalarLogical(1),  /* permuted  : "P" */
				  ScalarLogical(1),  /* LDL'	  : "D" */
				  ScalarLogical(0),  /* simplicial: "s" */
				  ScalarReal(0.));   /* Imult */
    L = AS_CHM_FR(Chol);
    R_CheckStack();
    cx = cholmod_spsolve(CHOLMOD_A, L, cb, &c);
    return chm_sparse_to_SEXP(cx, /*do_free*/ 1, /*uploT*/ 0,
			      /*Rkind*/ 0, /*diag*/ "N",
			      /*dimnames = */ R_NilValue);
}

SEXP dsCMatrix_matrix_solve(SEXP a, SEXP b)
{
    SEXP Chol = get_factor_pattern(a, "...Cholesky", 3);
    CHM_FR L;
    CHM_DN cx, cb = AS_CHM_DN(PROTECT(mMatrix_as_dgeMatrix(b)));
    R_CheckStack();

    if (Chol == R_NilValue) /* compute (and cache) "sPDCholesky" */
	Chol = dsCMatrix_Cholesky(a,
				  ScalarLogical(1),  /* permuted  : "P" */
				  ScalarLogical(1),  /* LDL'      : "D" */
				  ScalarLogical(0),  /* simplicial: "s" */
				  ScalarReal(0.));   /* Imult */
    L = AS_CHM_FR(Chol);
    R_CheckStack();
    cx = cholmod_solve(CHOLMOD_A, L, cb, &c);
    UNPROTECT(1);
    return chm_dense_to_SEXP(cx, 1, 0, /*dimnames = */ R_NilValue);
}

/* Needed for printing dsCMatrix objects */
/* FIXME: Create a more general version of this operation: also for lsC, (dsR?),..
*         e.g. make  compressed_to_dgTMatrix() in ./dgCMatrix.c work for dsC */
SEXP dsCMatrix_to_dgTMatrix(SEXP x)
{
    CHM_SP A = AS_CHM_SP(x);
    CHM_SP Afull = cholmod_copy(A, /*stype*/ 0, /*mode*/ 1, &c);
    CHM_TR At = cholmod_sparse_to_triplet(Afull, &c);
    R_CheckStack();

    if (!A->stype)
	error("Non-symmetric matrix passed to dsCMatrix_to_dgTMatrix");
    cholmod_free_sparse(&Afull, &c);
    return chm_triplet_to_SEXP(At, 1, /*uploT*/ 0, /*Rkind*/ 0, "",
			       GET_SLOT(x, Matrix_DimNamesSym));
}
