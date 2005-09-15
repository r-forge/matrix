#include "Pattern.h"

cholmod_sparse *factor_to_pattern(SEXP fact)
{
    int *vals = INTEGER(fact), i, n = LENGTH(fact),
	nlev = LENGTH(getAttrib(fact, R_LevelsSymbol));
    cholmod_sparse *A;

    if (!isFactor(fact))
	error("`fact' must be a factor");
    A = cholmod_allocate_sparse((size_t) nlev, (size_t) n, (size_t) n,
				TRUE, TRUE, 0, CHOLMOD_PATTERN, &c);
    for (i = 0; i <= n; i++) ((int *) (A->p))[i] = i;
    for (i = 0; i < n; i++) ((int *) (A->i))[i] = vals[i] - 1;
    if (cholmod_check_sparse(A, &c) < 0)
	error("sparse Pattern matrix has invalid structure");
    return A;
}

SEXP factor_prod(SEXP f1, SEXP f2)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    int *dims, n = LENGTH(f1), i, nl1, nl2, nmax, nnz, super;
    cholmod_triplet *A;
    cholmod_sparse *B;
    cholmod_factor *F;

    if (!isFactor(f1) || !isFactor(f2) || LENGTH(f2) != n)
	error("f1 and f2 must be factors of the same length");
    nl1 = LENGTH(getAttrib(f1, R_LevelsSymbol));
    nl2 = LENGTH(getAttrib(f2, R_LevelsSymbol));
    
    A = cholmod_allocate_triplet((size_t) nl1, (size_t) nl2, (size_t) n,
				 0, CHOLMOD_PATTERN, &c);
    for (i = 0; i < n; i++) {
	((int *)(A->i))[i] = INTEGER(f1)[i] - 1;
	((int *)(A->j))[i] = INTEGER(f2)[i] - 1;
    }
    A->nnz = n;
    nmax = nl1 * nl2;
    if (n < nmax) nmax = n;
    B = cholmod_triplet_to_sparse(A, nmax, &c);
    cholmod_free_triplet(&A, &c);
				/*  force a simplicial factorization */
    super = c.supernodal;
    c.supernodal = CHOLMOD_SIMPLICIAL;
    F = cholmod_analyze(B, &c); 
    c.supernodal = super;

    SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, 2));
    dims = INTEGER(VECTOR_ELT(ans, 0));
    dims[0] = nl1; dims[1] = nl2;
    SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, nl1));
    Memcpy(INTEGER(VECTOR_ELT(ans, 1)), (int *) F->Perm, nl1);
    SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, nl1));
    Memcpy(INTEGER(VECTOR_ELT(ans, 2)), (int *) F->ColCount, nl1);
    Rprintf("Ordering used: %d\n", F->ordering);
    Rprintf("is_super: %d\n", F->is_super);
    
    cholmod_free_sparse(&B, &c);
    cholmod_free_factor(&F, &c);

    UNPROTECT(1);
    return ans;
}

