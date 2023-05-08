#include <Rmath.h> /* math.h, logspace_add, logspace_sub */
#include "factorizations.h"

static cs *dgC2cs(SEXP obj)
{
    cs *A = (cs *) R_alloc(1, sizeof(cs)); 
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
	i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    A->nzmax = LENGTH(i);
    A->m = INTEGER(dim)[0];
    A->n = INTEGER(dim)[1];
    A->p = INTEGER(p);
    A->i = INTEGER(i);
    A->x = REAL(x);
    A->nz = -1;
    UNPROTECT(4);
    return A;
}

static SEXP cs2dgC(cs *A, const char *cl)
{
    int nnz = A->p[A->n];
    R_xlen_t np1 = (R_xlen_t) A->n + 1;
    SEXP obj = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dim = PROTECT(allocVector(INTSXP, 2)),
	p = PROTECT(allocVector(INTSXP, np1)),
	i = PROTECT(allocVector(INTSXP, nnz)),
	x = PROTECT(allocVector(REALSXP, nnz));
    INTEGER(dim)[0] = A->m;
    INTEGER(dim)[1] = A->n;
    Matrix_memcpy(INTEGER(p), A->p, np1, sizeof(int));
    Matrix_memcpy(INTEGER(i), A->i, nnz, sizeof(int));
    Matrix_memcpy(REAL(x), A->x, nnz, sizeof(double));
    SET_SLOT(obj, Matrix_DimSym, dim);
    SET_SLOT(obj, Matrix_pSym, p);
    SET_SLOT(obj, Matrix_iSym, i);
    SET_SLOT(obj, Matrix_xSym, x);
    UNPROTECT(5);
    return obj;
}

SEXP dgeMatrix_trf_(SEXP obj, int warn)
{
    SEXP val = get_factor(obj, "LU");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("denseLU"));
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
    int *pdim = INTEGER(dim), r = (pdim[0] < pdim[1]) ? pdim[0] : pdim[1];
    SET_SLOT(val, Matrix_DimSym, dim);
    SET_SLOT(val, Matrix_DimNamesSym, dimnames);
    if (r > 0) {
	PROTECT_INDEX pid;
	SEXP perm = PROTECT(allocVector(INTSXP, r)), x;
	PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pid);
	REPROTECT(x = duplicate(x), pid);
	int *pperm = INTEGER(perm), info;
	double *px = REAL(x);
	
	F77_CALL(dgetrf)(pdim, pdim + 1, px, pdim, pperm, &info);
	
	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dgetrf", info);
	else if (info > 0 && warn > 0) {
	    /* MJ: 'dgetrf' does not distinguish between singular, */
	    /*     finite matrices and matrices containing NaN ... */
	    /*     hence this message can mislead                  */
	    if (warn > 1)
		error  (_("LAPACK '%s': matrix is exactly singular, U[i,i]=0, i=%d"),
			"dgetrf", info);
	    else 
		warning(_("LAPACK '%s': matrix is exactly singular, U[i,i]=0, i=%d"),
			"dgetrf", info);
	}
	
	SET_SLOT(val, Matrix_permSym, perm);
	SET_SLOT(val, Matrix_xSym, x);
	UNPROTECT(2); /* x, perm */
    }
    set_factor(obj, "LU", val);
    UNPROTECT(3); /* dimnames, dim, val */
    return val;
}

SEXP dsyMatrix_trf_(SEXP obj, int warn)
{
    SEXP val = get_factor(obj, "BunchKaufman");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("BunchKaufman"));
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_DimSym, dim);
    set_symmetrized_DimNames(val, dimnames, -1);
    SET_SLOT(val, Matrix_uploSym, uplo);
    if (n > 0) {
	R_xlen_t nn;
	SEXP perm = PROTECT(allocVector(INTSXP, n)),
	    x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	    y = PROTECT(allocVector(REALSXP, nn = XLENGTH(x)));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int *pperm = INTEGER(perm), lwork = -1, info;
	double *px = REAL(x), *py = REAL(y), tmp, *work;

	Matrix_memset(py, 0, nn, sizeof(double));
	F77_CALL(dlacpy)(&ul, pdim, pdim, px, pdim, py, pdim FCONE);
	F77_CALL(dsytrf)(&ul, pdim, py, pdim, pperm, &tmp, &lwork, &info FCONE);
	lwork = (int) tmp;
	Matrix_Calloc(work, lwork, double);
	F77_CALL(dsytrf)(&ul, pdim, py, pdim, pperm, work, &lwork, &info FCONE);
	Matrix_Free(work, lwork);
	
	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dsytrf", info);
	else if (info > 0 && warn > 0) {
	    /* MJ: 'dsytrf' does not distinguish between singular, */
	    /*     finite matrices and matrices containing NaN ... */
	    /*     hence this message can mislead                  */
	    if (warn > 1)
		error  (_("LAPACK '%s': matrix is exactly singular, D[i,i]=0, i=%d"),
			"dsytrf", info);
	    else
		warning(_("LAPACK '%s': matrix is exactly singular, D[i,i]=0, i=%d"),
			"dsytrf", info);
	}
	
	SET_SLOT(val, Matrix_permSym, perm);
	SET_SLOT(val, Matrix_xSym, y);
	UNPROTECT(3); /* y, x, perm */
    }
    set_factor(obj, "BunchKaufman", val);
    UNPROTECT(4); /* uplo, dimnames, dim, val */
    return val;
}

SEXP dspMatrix_trf_(SEXP obj, int warn)
{
    SEXP val = get_factor(obj, "pBunchKaufman");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("pBunchKaufman"));
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_DimSym, dim);
    set_symmetrized_DimNames(val, dimnames, -1);
    SET_SLOT(val, Matrix_uploSym, uplo);
    if (n > 0) {
	PROTECT_INDEX pid;
	SEXP perm = PROTECT(allocVector(INTSXP, n)), x;
	PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pid);
	REPROTECT(x = duplicate(x), pid);
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int *pperm = INTEGER(perm), info;
	double *px = REAL(x);
    
	F77_CALL(dsptrf)(&ul, pdim, px, pperm, &info FCONE);
    
	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dsptrf", info);
	else if (info > 0 && warn > 0) {
	    /* MJ: 'dsptrf' does not distinguish between singular, */
	    /*     finite matrices and matrices containing NaN ... */
	    /*     hence this message can mislead                  */
	    if (warn > 1)
		error  (_("LAPACK '%s': matrix is exactly singular, D[i,i]=0, i=%d"),
			"dsptrf", info);
	    else
		warning(_("LAPACK '%s': matrix is exactly singular, D[i,i]=0, i=%d"),
			"dsptrf", info);
	}

	SET_SLOT(val, Matrix_permSym, perm);
	SET_SLOT(val, Matrix_xSym, x);
	UNPROTECT(2); /* x, perm */
    }
    set_factor(obj, "pBunchKaufman", val);
    UNPROTECT(4); /* uplo, dimnames, dim, val */
    return val;
}

SEXP dpoMatrix_trf_(SEXP obj, int warn)
{
    SEXP val = get_factor(obj, "Cholesky");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("Cholesky"));
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_DimSym, dim);
    set_symmetrized_DimNames(val, dimnames, -1);
    SET_SLOT(val, Matrix_uploSym, uplo);
    if (n > 0) {
	R_xlen_t nn;
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	    y = PROTECT(allocVector(REALSXP, nn = XLENGTH(x)));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int info;
	double *px = REAL(x), *py = REAL(y);

	Matrix_memset(py, 0, nn, sizeof(double));
	F77_CALL(dlacpy)(&ul, pdim, pdim, px, pdim, py, pdim FCONE);
	F77_CALL(dpotrf)(&ul, pdim, py, pdim, &info FCONE);

	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dpotrf", info);
	else if (info > 0) {
	    if (warn > 1)
		error  (_("LAPACK '%s': leading minor of order %d is not positive definite"),
			"dpotrf", info);
	    else if (warn > 0)
		warning(_("LAPACK '%s': leading minor of order %d is not positive definite"),
			"dpotrf", info);
	    UNPROTECT(6); /* y, x, uplo, dimnames, dim, val */
	    return ScalarInteger(info);
	}
	
	SET_SLOT(val, Matrix_xSym, y);
	UNPROTECT(2); /* y, x */
    }
    set_factor(obj, "Cholesky", val);
    UNPROTECT(4); /* uplo, dimnames, dim, val */
    return val;
}

SEXP dppMatrix_trf_(SEXP obj, int warn)
{
    SEXP val = get_factor(obj, "pCholesky");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("pCholesky"));
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    SET_SLOT(val, Matrix_DimSym, dim);
    set_symmetrized_DimNames(val, dimnames, -1);
    SET_SLOT(val, Matrix_uploSym, uplo);
    if (n > 0) {
	PROTECT_INDEX pid;
	SEXP x;
	PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pid);
	REPROTECT(x = duplicate(x), pid);
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int info;
	double *px = REAL(x);
	
	F77_CALL(dpptrf)(&ul, pdim, px, &info FCONE);

	if (info < 0)
	    error(_("LAPACK '%s' gave error code %d"),
		  "dpptrf", info);
	else if (info > 0) {
	    if (warn > 1)
		error  (_("LAPACK '%s': leading minor of order %d is not positive definite"),
			"dpptrf", info);
	    else if (warn > 0)
		warning(_("LAPACK '%s': leading minor of order %d is not positive definite"),
			"dpptrf", info);
	    UNPROTECT(5); /* x, uplo, dimnames, dim, val */
	    return ScalarInteger(info);
	}
	
	SET_SLOT(val, Matrix_xSym, x);
	UNPROTECT(1); /* x */
    }
    set_factor(obj, "pCholesky", val);
    UNPROTECT(4); /* uplo, dimnames, dim, val */
    return val;
}

SEXP dgeMatrix_trf(SEXP obj, SEXP warn)
{
    return dgeMatrix_trf_(obj, asInteger(warn));
}

SEXP dsyMatrix_trf(SEXP obj, SEXP warn)
{
    return dsyMatrix_trf_(obj, asInteger(warn));
}

SEXP dspMatrix_trf(SEXP obj, SEXP warn)
{
    return dspMatrix_trf_(obj, asInteger(warn));
}

SEXP dpoMatrix_trf(SEXP obj, SEXP warn)
{
    return dpoMatrix_trf_(obj, asInteger(warn));
}

SEXP dppMatrix_trf(SEXP obj, SEXP warn)
{
    return dppMatrix_trf_(obj, asInteger(warn));
}

int dgCMatrix_trf_(cs *A, css **S, csn **N, int doError,
		   int order, double tol)
{
    if (A->m != A->n)
	error(_("LU factorization of m-by-n dgCMatrix requires m == n"));
    /* Symbolic analysis : */
    *S = cs_sqr(order, A, 0);
    /* Numeric factorization : */
    *N = cs_lu(A, *S, tol);
    if (*N) {
	cs *T;
	/* Drop zeros from L and sort it : */
	cs_dropzeros((*N)->L);
	T = cs_transpose((*N)->L, 1);
	cs_spfree((*N)->L);
	(*N)->L = cs_transpose(T, 1);
	cs_spfree(T);
	/* Drop zeros from U and sort it : */
	cs_dropzeros((*N)->U);
	T = cs_transpose((*N)->U, 1);
	cs_spfree((*N)->U);
	(*N)->U = cs_transpose(T, 1);
	cs_spfree(T);
	return 0;
    } else {
	if (*S) *S = cs_sfree(*S);
	if(doError)
	    error(_("LU factorization of dgCMatrix failed: out of memory or near-singular"));
	return 1;
    }
}

SEXP dgCMatrix_trf(SEXP obj, SEXP doError, SEXP keepDimNames,
		   SEXP order, SEXP tol)
{
    SEXP val = get_factor(obj, "LU");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("sparseLU"));

    double tol_ = asReal(tol);
    if (ISNAN(tol_))
	error(_("'tol' is not a number"));
    
    int order_ = asInteger(order);
    if (order_ == NA_INTEGER)
	order_ = (tol_ == 1.0) ? 2 : 1;
    else if (order_ < 0)
	order_ = 0;
    else if (order_ > 3)
	order_ = 3;

    cs *A = dgC2cs(obj);
    css *S;
    csn *N;
    if (dgCMatrix_trf_(A, &S, &N, asLogical(doError), order_, tol_)) {
	if (N) cs_nfree(N);
	if (S) cs_sfree(S);
	UNPROTECT(1); /* val */
	/* Defensive code will check with isS4 : */
	return ScalarLogical(NA_LOGICAL);
    }
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    SET_SLOT(val, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */
    if (asLogical(keepDimNames)) {
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(val, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */
    }
    
    SEXP L = PROTECT(cs2dgC(N->L, "dtCMatrix")),
	U = PROTECT(cs2dgC(N->U, "dtCMatrix")),
	uplo = PROTECT(mkString("L"));
    SET_SLOT(L, Matrix_uploSym, uplo);
    SET_SLOT(val, Matrix_LSym, L);
    SET_SLOT(val, Matrix_USym, U);
    UNPROTECT(3); /* uplo, U, L */

    SEXP p = PROTECT(allocVector(INTSXP, A->m));
    int *pp = cs_pinv(N->pinv, A->m);
    Matrix_memcpy(INTEGER(p), pp, A->m, sizeof(int));
    SET_SLOT(val, Matrix_pSym, p);
    UNPROTECT(1); /* p */
    if (order_ > 0) {
	SEXP q = PROTECT(allocVector(INTSXP, A->n));
	int *pq = S->q;
	Matrix_memcpy(INTEGER(q), pq, A->n, sizeof(int));
	SET_SLOT(val, Matrix_qSym, q);
	UNPROTECT(1); /* q */
    }

    cs_free(pp);
    cs_nfree(N);
    cs_sfree(S);

    set_factor(obj, "LU", val);
    UNPROTECT(1); /* val */
    return val;
}

int dgCMatrix_orf_(cs *A, css **S, csn **N, int doError,
		   int order)
{
    if (A->m < A->n)
	error(_("QR factorization of m-by-n dgCMatrix requires m >= n"));
    /* Symbolic analysis : */
    *S = cs_sqr(order, A, 1);
    /* Numeric factorization : */
    *N = cs_qr(A, *S);
    if (*N) {
	cs *T;
	/* Drop zeros from V and sort it : */
	cs_dropzeros((*N)->L);
	T = cs_transpose((*N)->L, 1);
	cs_spfree((*N)->L);
	(*N)->L = cs_transpose(T, 1);
	cs_spfree(T);
	/* Drop zeros from R and sort it : */
	cs_dropzeros((*N)->U);
	T = cs_transpose((*N)->U, 1);
	cs_spfree((*N)->U);
	(*N)->U = cs_transpose(T, 1);
	cs_spfree(T);
	return 0;
    } else {
	if (*S) *S = cs_sfree(*S);
	if(doError)
	    error(_("QR factorization of dgCMatrix failed: out of memory"));
	return 1;
    }
}

SEXP dgCMatrix_orf(SEXP obj, SEXP doError, SEXP keepDimNames,
		   SEXP order)
{
    SEXP val = get_factor(obj, "QR");
    if (!isNull(val))
	return val;
    PROTECT(val = NEW_OBJECT_OF_CLASS("sparseQR"));

    int order_ = asInteger(order);
    if (order_ < 0)
	order_ = 0;
    else if (order_ > 3)
	order_ = 3;

    cs *A = dgC2cs(obj);
    css *S;
    csn *N;
    if (dgCMatrix_orf_(A, &S, &N, asLogical(doError), order_)) {
	if (N) cs_nfree(N);
	if (S) cs_sfree(S);
	UNPROTECT(1); /* val */
	/* Defensive code will check with isS4 : */
	return ScalarLogical(NA_LOGICAL);
    }
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    SET_SLOT(val, Matrix_DimSym, dim);
    UNPROTECT(1); /* dim */
    if (asLogical(keepDimNames)) {
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(val, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */
    }
    
    SEXP V = PROTECT(cs2dgC(N->L, "dgCMatrix")),
	R = PROTECT(cs2dgC(N->U, "dgCMatrix"));
    SET_SLOT(val, Matrix_VSym, V);
    SET_SLOT(val, Matrix_RSym, R);
    UNPROTECT(2); /* R, V */

    SEXP beta = PROTECT(allocVector(REALSXP, A->n));
    double *pbeta = N->B;
    Matrix_memcpy(REAL(beta), pbeta, A->n, sizeof(double));
    SET_SLOT(val, Matrix_betaSym, beta);
    UNPROTECT(1); /* beta */
    
    SEXP p = PROTECT(allocVector(INTSXP, S->m2));
    int *pp = cs_pinv(S->pinv, S->m2);
    Matrix_memcpy(INTEGER(p), pp, S->m2, sizeof(int));
    SET_SLOT(val, Matrix_pSym, p);
    UNPROTECT(1); /* p */
    if (order_ > 0) {
	SEXP q = PROTECT(allocVector(INTSXP, A->n));
	int *pq = S->q;
	Matrix_memcpy(INTEGER(q), pq, A->n, sizeof(int));
	SET_SLOT(val, Matrix_qSym, q);
	UNPROTECT(1); /* q */
    }

    cs_free(pp);
    cs_nfree(N);
    cs_sfree(S);

    set_factor(obj, "QR", val);
    UNPROTECT(1); /* val */
    return val;
}

int dpCMatrix_trf_(CHM_SP A, CHM_FR *L,
		   int perm, int ldl, int super, double mult)
{
    CHM_store_common();
    
    double beta[2];
    beta[0] = mult;
    beta[1] = 0.0;
    
    if (!perm) {
	/* Require identity permutation : */
	c.nmethods = 1;
	c.method[0].ordering = CHOLMOD_NATURAL;
	c.postorder = FALSE;
    }
    c.final_ll = ldl == 0;
    c.supernodal = (super == NA_LOGICAL || super < 0) ? CHOLMOD_AUTO :
	((super > 0) ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL);
    
    *L = cholmod_analyze(A, &c);
    int res = cholmod_factorize_p(A, beta, (int *) NULL, 0, *L, &c);
    
    CHM_restore_common();

    return res;
}

SEXP dpCMatrix_trf(SEXP obj,
		   SEXP perm, SEXP ldl, SEXP super, SEXP mult)
{
    int perm_ = asLogical(perm), ldl_ = asLogical(ldl),
	super_ = asInteger(super);
    double mult_ = asReal(mult);
    if (!R_FINITE(mult_))
	error(_("'mult' is not a number or not finite"));
    char nm[] = "spdCholesky";
    if (perm_)
	nm[1] = 'P';
    if (ldl_)
	nm[2] = 'D';
    SEXP ch = R_NilValue;
    if (super_ == NA_LOGICAL || !super_)
	ch = get_factor(obj, nm);
    if (isNull(ch) && (super_ == NA_LOGICAL || super_)) {
	nm[0] = 'S';
	ch = get_factor(obj, nm);
    }
    int cached = !isNull(ch);
    if (cached && mult_ == 0.0)
	return ch;

    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(ch, &pid);
    CHM_SP A = AS_CHM_SP__(obj);
    CHM_FR L;
    R_CheckStack();

    if (cached) {
	double beta[2];
	beta[0] = mult_;
	beta[1] = 0.0;
	L = AS_CHM_FR(ch);
	R_CheckStack();
	L = cholmod_copy_factor(L, &c);
	cholmod_factorize_p(A, beta, (int *) NULL, 0, L, &c);
    } else {
	dpCMatrix_trf_(A, &L, perm_, ldl_, super_, mult_);
	if (super_ == NA_LOGICAL)
	    nm[0] = (L->is_super) ? 'S' : 's';
    }
    REPROTECT(ch = chm_factor_to_SEXP(L, 1), pid);
    
    SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
    set_symmetrized_DimNames(ch, dimnames, -1);
    UNPROTECT(1); /* dimnames */
    
    if (!cached && mult_ == 0.0)
	set_factor(obj, nm, ch);
    UNPROTECT(1); /* ch */
    return ch;
}

SEXP denseLU_expand(SEXP obj)
{
    /* A = P L U   <=>   P' A = L U   where ... 
       
       A -> [m,n]
       P -> [m,m], permutation
       L -> if m >= n then [m,n] else [m,m], lower trapezoidal, unit diagonal
       U -> if m <= n then [m,n] else [n,n], upper trapezoidal
       
       square L,U given as dtrMatrix with appropriate 'uplo', 'diag' slots,
       non-square L,U given as dgeMatrix
    */
    
    const char *nms[] = {"P", "L", "U", ""};
    PROTECT_INDEX pidA, pidB;
    SEXP res = PROTECT(Rf_mkNamed(VECSXP, nms)),
	P = PROTECT(NEW_OBJECT_OF_CLASS("pMatrix")),
	dim, x;
    PROTECT_WITH_INDEX(dim = GET_SLOT(obj, Matrix_DimSym), &pidA);
    PROTECT_WITH_INDEX(x = GET_SLOT(obj, Matrix_xSym), &pidB);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n, j;
    
    if (m == n) {
	SEXP L = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    U = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    uplo = PROTECT(mkString("L")),
	    diag = PROTECT(mkString("U"));
	SET_SLOT(L, Matrix_DimSym, dim);
	SET_SLOT(U, Matrix_DimSym, dim);
	SET_SLOT(P, Matrix_DimSym, dim);
	SET_SLOT(L, Matrix_uploSym, uplo);
	SET_SLOT(L, Matrix_diagSym, diag);
	SET_SLOT(L, Matrix_xSym, x);
	SET_SLOT(U, Matrix_xSym, x);
	SET_VECTOR_ELT(res, 1, L);
	SET_VECTOR_ELT(res, 2, U);
	UNPROTECT(4); /* diag, uplo, U, L */
    } else {
	SEXP G = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
	    T = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    y = PROTECT(allocVector(REALSXP, (R_xlen_t) r * r));
	REPROTECT(x = duplicate(x), pidB);
	double *px = REAL(x), *py = REAL(y);
	
	SET_SLOT(G, Matrix_DimSym, dim);
	REPROTECT(dim = allocVector(INTSXP, 2), pidA);
	pdim = INTEGER(dim);
	pdim[0] = pdim[1] = r;
	SET_SLOT(T, Matrix_DimSym, dim);
	REPROTECT(dim = allocVector(INTSXP, 2), pidA);
	pdim = INTEGER(dim);
	pdim[0] = pdim[1] = m;
	SET_SLOT(P, Matrix_DimSym, dim);
	
	if (m < n) {
            /* G is upper trapezoidal, T is unit lower triangular */
	    SEXP uplo = PROTECT(mkString("L")),
		diag = PROTECT(mkString("U"));
	    SET_SLOT(T, Matrix_uploSym, uplo);
	    SET_SLOT(T, Matrix_diagSym, diag);
	    UNPROTECT(2); /* diag, uplo */

	    Matrix_memcpy(py, px, (R_xlen_t) m * m, sizeof(double));
	    ddense_unpacked_make_triangular(px, m, n, 'U', 'N');
	} else {
            /* G is unit lower trapezoidal, T is upper triangular */
	    double *tmp = px;
	    for (j = 0; j < n; ++j, px += m, py += r)
		Matrix_memcpy(py, px, j+1, sizeof(double));
	    ddense_unpacked_make_triangular(tmp, m, n, 'L', 'U');
	}
	SET_SLOT(G, Matrix_xSym, x);
	SET_SLOT(T, Matrix_xSym, y);
	
	SET_VECTOR_ELT(res, (m < n) ? 2 : 1, G);
	SET_VECTOR_ELT(res, (m < n) ? 1 : 2, T);
	UNPROTECT(3); /* y, T, G */
    }

    SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	perm = PROTECT(allocVector(INTSXP, m));
    int *ppivot = INTEGER(pivot), *pperm = INTEGER(perm), *pinvperm, pos, tmp;
    Matrix_Calloc(pinvperm, m, int);

    for (j = 0; j < m; ++j) /* initialize inverse permutation */
	pinvperm[j] = j;
    for (j = 0; j < r; ++j) { /* generate inverse permutation */
	pos = ppivot[j] - 1;
	if (pos != j) {
	    tmp = pinvperm[j];
	    pinvperm[j] = pinvperm[pos];
	    pinvperm[pos] = tmp;
	}
    }
    for (j = 0; j < m; ++j) /* invert inverse permutation (0->1-based) */
	pperm[pinvperm[j]] = j + 1;
    Matrix_Free(pinvperm, m);

    SET_SLOT(P, Matrix_permSym, perm);
    SET_VECTOR_ELT(res, 0, P);
    
#define DO_DIMNAMES							\
    do {								\
	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),		\
	    ndn = PROTECT(getAttrib(dn, R_NamesSymbol));		\
	for (j = 0; j < 2; ++j) {					\
	    SEXP rn = VECTOR_ELT(dn, j);				\
	    if (!isNull(rn) || !isNull(ndn)) {				\
		SEXP X = VECTOR_ELT(res, (j == 0) ? 0 : LENGTH(res) - 1), \
		    dnX = PROTECT(allocVector(VECSXP, 2));		\
		SET_VECTOR_ELT(dnX, j, rn);				\
		if (!isNull(ndn)) {					\
		    SEXP ndnX = PROTECT(allocVector(STRSXP, 2));	\
		    SET_STRING_ELT(ndnX, j, STRING_ELT(ndn, j));	\
		    setAttrib(dnX, R_NamesSymbol, ndnX);		\
		    UNPROTECT(1);					\
		}							\
		SET_SLOT(X, Matrix_DimNamesSym, dnX);			\
		UNPROTECT(1);						\
	    }								\
	}								\
	UNPROTECT(2);							\
    } while (0)
    
    DO_DIMNAMES;
    UNPROTECT(6); /* perm, pivot, x, dim, P, res */
    return res;
}

SEXP BunchKaufman_expand(SEXP obj)
{
    SEXP P_ = PROTECT(NEW_OBJECT_OF_CLASS("pMatrix")),
	T_ = PROTECT(NEW_OBJECT_OF_CLASS("dtCMatrix")),
	D_ = PROTECT(NEW_OBJECT_OF_CLASS("dsCMatrix")),
	dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int i, j, s, n = INTEGER(dim)[0];
    R_xlen_t n1a = (R_xlen_t) n + 1;
    if (n > 0) {
	SET_SLOT(P_, Matrix_DimSym, dim);
	SET_SLOT(T_, Matrix_DimSym, dim);
	SET_SLOT(D_, Matrix_DimSym, dim);
    }
    UNPROTECT(1); /* dim */

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    if (!upper) {
	SET_SLOT(T_, Matrix_uploSym, uplo);
	SET_SLOT(D_, Matrix_uploSym, uplo);
    }
    UNPROTECT(1); /* uplo */
    
    SEXP diag = PROTECT(mkString("U"));
    SET_SLOT(T_, Matrix_diagSym, diag);
    UNPROTECT(1); /* diag */
    
    SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	D_p = PROTECT(allocVector(INTSXP, n1a));
    int *ppivot = INTEGER(pivot), *D_pp = INTEGER(D_p),
	b = n, dp = (upper) ? 1 : 2;
    D_pp[0] = 0;
    j = 0;
    while (j < n) {
	if (ppivot[j] > 0) {
	    D_pp[j+1] = D_pp[j] + 1;
	    j += 1;
	} else {
	    D_pp[j+1] = D_pp[j] + dp;
	    D_pp[j+2] = D_pp[j] + 3;
	    j += 2;
	    --b;
	}
    }
    SET_SLOT(D_, Matrix_pSym, D_p);
    UNPROTECT(1); /* D_p */

    SEXP P, P_perm, T, T_p, T_i, T_x,
	D_i = PROTECT(allocVector(INTSXP, D_pp[n])),
	D_x = PROTECT(allocVector(REALSXP, D_pp[n])),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    int *P_pperm, *T_pp, *T_pi, *D_pi = INTEGER(D_i);
    double *T_px, *D_px = REAL(D_x), *px = REAL(x);

    int unpacked = (double) n * n <= R_XLEN_T_MAX &&
	(R_xlen_t) n * n == XLENGTH(x);

    R_xlen_t len = (R_xlen_t) 2 * b + 1, k = (upper) ? len - 1 : 0;
    SEXP res = PROTECT(allocVector(VECSXP, len));

    j = 0;
    while (b--) {
	s = (ppivot[j] > 0) ? 1 : 2;
	dp = (upper) ? j : n - j - s;
	
	PROTECT(P = duplicate(P_));
	PROTECT(P_perm = allocVector(INTSXP, n));
	PROTECT(T = duplicate(T_));
	PROTECT(T_p = allocVector(INTSXP, n1a));
	PROTECT(T_i = allocVector(INTSXP, (R_xlen_t) s * dp));
	PROTECT(T_x = allocVector(REALSXP, (R_xlen_t) s * dp));
	
	P_pperm = INTEGER(P_perm);
	T_pp = INTEGER(T_p);
	T_pi = INTEGER(T_i);
	T_px = REAL(T_x);
	T_pp[0] = 0;
	
	for (i = 0; i < j; ++i) {
	    T_pp[i+1] = 0;
	    P_pperm[i] = i + 1;
	}
	for (i = j; i < j+s; ++i) {
	    T_pp[i+1] = T_pp[i] + dp;
	    P_pperm[i] = i + 1;
	}
	for (i = j+s; i < n; ++i) {
	    T_pp[i+1] = T_pp[i];
	    P_pperm[i] = i + 1;
	}
	
	if (s == 1) {
	    P_pperm[j] = ppivot[j];
	    P_pperm[ppivot[j]-1] = j + 1;
	} else if (upper) {
	    P_pperm[j] = -ppivot[j];
	    P_pperm[-ppivot[j]-1] = j + 1;
	} else {
	    P_pperm[j+1] = -ppivot[j];
	    P_pperm[-ppivot[j]-1] = j + 2;
	}

	if (upper) {
	    for (i = 0; i < j; ++i) {
		*(T_pi++) = i;
		*(T_px++) = *(px++);
	    }
	    *(D_pi++) = j;
	    *(D_px++) = *(px++);
	    ++j;
	    if (unpacked)
		px += n - j;
	    if (s == 2) {
		for (i = 0; i < j-1; ++i) {
		    *(T_pi++) = i;
		    *(T_px++) = *(px++);
		}
		*(D_pi++) = j - 1;
		*(D_pi++) = j;
		*(D_px++) = *(px++);
		*(D_px++) = *(px++);
		++j;
		if (unpacked)
		    px += n - j;
	    }
	} else {
	    if (s == 2) {
		*(D_pi++) = j;
		*(D_pi++) = j + 1;
		*(D_px++) = *(px++);
		*(D_px++) = *(px++);
		for (i = j+2; i < n; ++i) {
		    *(T_pi++) = i;
		    *(T_px++) = *(px++);
		}
		++j;
		if (unpacked)
		    px += j;
	    }
	    *(D_pi++) = j;
	    *(D_px++) = *(px++);
	    for (i = j+1; i < n; ++i) {
		*(T_pi++) = i;
		*(T_px++) = *(px++);
	    }
	    ++j;
	    if (unpacked)
		px += j;
	}

	SET_SLOT(P, Matrix_permSym, P_perm);
	SET_SLOT(T, Matrix_pSym, T_p);
	SET_SLOT(T, Matrix_iSym, T_i);
	SET_SLOT(T, Matrix_xSym, T_x);

	if (upper) {
	    SET_VECTOR_ELT(res, k-1, P);
	    SET_VECTOR_ELT(res, k  , T);
	    k -= 2;
	} else {
	    SET_VECTOR_ELT(res, k  , P);
	    SET_VECTOR_ELT(res, k+1, T);
	    k += 2;
	}
	UNPROTECT(6); /* T_x, T_i, T_p, T, P_perm, P */
    }
    
    SET_SLOT(D_, Matrix_iSym, D_i);
    SET_SLOT(D_, Matrix_xSym, D_x);
    SET_VECTOR_ELT(res, k, D_);

    DO_DIMNAMES;
    
#undef DO_DIMNAMES
    
    UNPROTECT(8); /* res, x, D_x, D_i, pivot, D_, T_, P_ */ 
    return res;
}

SEXP denseLU_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("determinant of non-square matrix is undefined"));
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = 0.0; /* result for n == 0 */
    if (n > 0) {
	SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	    x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int j, *ppivot = INTEGER(pivot);
	R_xlen_t n1a = (R_xlen_t) n + 1;
	double *px = REAL(x);
	
	for (j = 0; j < n; ++j, px += n1a, ++ppivot) {
	    if (*px < 0.0) {
		modulus += log(-(*px));
		if (*ppivot == j + 1)
		    sign = -sign;
	    } else {
		/* incl. 0, NaN cases */
		modulus += log(*px);
		if (*ppivot != j + 1)
		    sign = -sign;
	    }
	}
	UNPROTECT(2); /* x, pivot */
    }
    if (!givelog)
	modulus = exp(modulus);
    return as_det_obj(modulus, givelog, sign);
}

SEXP sparseLU_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = 0.0; /* result for n == 0 */
    if (n > 0) {
	SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym)),
	    p = PROTECT(GET_SLOT(U, Matrix_pSym)),
	    i = PROTECT(GET_SLOT(U, Matrix_iSym)),
	    x = PROTECT(GET_SLOT(U, Matrix_xSym));
	int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
	double *px = REAL(x);

	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    if (kend > k && pi[kend - 1] == j) {
		if (px[kend - 1] < 0.0) {
		    modulus += log(-px[kend - 1]);
		    sign = -sign;
		} else {
		    /* incl. 0, NaN cases */
		    modulus += log(px[kend - 1]);
		}
	    } else {
		UNPROTECT(4); /* x, i, p, U */
		return as_det_obj((givelog) ? 0.0 : 1.0, givelog, 1);
	    }
	    k = kend;
	}
	UNPROTECT(4); /* x, i, p, U */

	PROTECT(p = GET_SLOT(obj, Matrix_pSym));
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
	    sign = -sign;
	UNPROTECT(1); /* p */
	PROTECT(p = GET_SLOT(obj, Matrix_qSym));
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
	    sign = -sign;
	UNPROTECT(1); /* p */
    }
    if (!givelog)
	modulus = exp(modulus);
    return as_det_obj(modulus, givelog, sign);
}

SEXP sparseQR_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("determinant of non-square matrix is undefined"));
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = 0.0; /* result for n == 0 */
    if (n > 0) {
	SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
	PROTECT(dim = GET_SLOT(R, Matrix_DimSym));
	if (INTEGER(dim)[0] > n)
	    error(_("determinant(<sparseQR>) does not support structurally rank deficient case"));
	UNPROTECT(1); /* dim */
	
	SEXP p = PROTECT(GET_SLOT(R, Matrix_pSym)),
	    i = PROTECT(GET_SLOT(R, Matrix_iSym)),
	    x = PROTECT(GET_SLOT(R, Matrix_xSym));
	int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
	double *px = REAL(x);
	
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    if (kend > k && pi[kend - 1] == j) {
		if (px[kend - 1] < 0.0) {
		    modulus += log(-px[kend - 1]);
		    sign = -sign;
		} else {
		    /* incl. 0, NaN cases */
		    modulus += log(px[kend - 1]);
		}
	    } else {
		UNPROTECT(4); /* x, i, p, R */
		return as_det_obj((givelog) ? 0.0 : 1.0, givelog, 1);
	    }
	    k = kend;
	}
	UNPROTECT(4); /* x, i, p, U */

	PROTECT(p = GET_SLOT(obj, Matrix_pSym));
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
	    sign = -sign;
	UNPROTECT(1); /* p */
	PROTECT(p = GET_SLOT(obj, Matrix_qSym));
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
	    sign = -sign;
	UNPROTECT(1); /* p */
    }
    if (!givelog)
	modulus = exp(modulus);
    return as_det_obj(modulus, givelog, sign);
}

SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = 0.0; /* result for n == 0 */
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    if (n > 0) {
	int unpacked = (double) n * n <= R_XLEN_T_MAX &&
	    (R_xlen_t) n * n == XLENGTH(x);	

	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1); /* uplo */
	
	SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym));
	int j = 0, *ppivot = INTEGER(pivot);
	R_xlen_t n1a = (R_xlen_t) n + 1;
	double *px = REAL(x), a, b, c, logab, logcc;	
	while (j < n) {
	    if (ppivot[j] > 0) {
		if (*px < 0.0) {
		    modulus += log(-(*px));
		    sign = -sign;
		} else {
		    /* incl. 0, NaN cases */
		    modulus += log(*px);
		}
		px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
		j += 1;
	    } else {
		a = *px;
		if (upper) {
		    px += (unpacked) ? n1a : j + 2;
		    b = *px;
		    c = *(px - 1);
		    px += (unpacked) ? n1a : j + 3;
		} else {
		    c = *(px + 1);
		    px += (unpacked) ? n1a : n - j;
		    b = *px;
		    px += (unpacked) ? n1a : n - j - 1;
		}
		logab = log(fabs(a)) + log(fabs(b));
		logcc = 2.0 * log(fabs(c));
		if ((a < 0.0) != (b < 0.0)) {
		    /* det = ab - cc = -(abs(ab) + cc) < 0 */
		    modulus += logspace_add(logab, logcc);
		    sign = -sign;
		} else if (logab < logcc) {
		    /* det = ab - cc = -(cc - ab) < 0 */
		    modulus += logspace_sub(logcc, logab);
		    sign = -sign;
		} else {
		    /* det = ab - cc > 0 */
		    modulus += logspace_sub(logab, logcc);
		}
		j += 2;
	    }
	}
	UNPROTECT(2); /* x, pivot */
    }
    if (!givelog)
	modulus = exp(modulus);
    return as_det_obj(modulus, givelog, sign);
}

SEXP Cholesky_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = 0.0; /* result for n == 0 */
    if (n > 0) {
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	int unpacked = (double) n * n <= R_XLEN_T_MAX &&
	    (R_xlen_t) n * n == XLENGTH(x);	
	SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
	int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1); /* uplo */
	
	R_xlen_t n1a = (R_xlen_t) n + 1;
	double *px = REAL(x);
	for (int j = 0; j < n; ++j) {
	    if (*px < 0.0) {
		modulus += log(-(*px));
		sign = -sign;
	    } else {
		/* incl. 0, NaN cases */
		modulus += log(*px);
	    }
	    px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
	}
	UNPROTECT(1); /* x */
    }
    if (!givelog)
	modulus = exp(modulus);
    return as_det_obj(modulus, givelog, sign);
}

SEXP CHMfactor_determinant(SEXP obj, SEXP logarithm)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int n = INTEGER(dim)[0];
    UNPROTECT(1); /* dim */
    int givelog = asLogical(logarithm) != 0, sign = 1;
    double modulus = 0.0; /* result for n == 0 */
    if (n > 0) {
	CHM_FR L = AS_CHM_FR(obj);
	R_CheckStack();
	cholmod_change_factor(CHOLMOD_REAL, L->is_ll, 0, 0, 0, L, &c);
	int j, *pp = L->p;
	double *px = L->x;
	for (j = 0; j < n; ++j) {
	    /* incl. 0, NaN cases */
	    modulus += log(px[pp[j]]);
	}
	if (L->is_ll)
	    modulus *= 2.0;
    }
    if (!givelog)
	modulus = exp(modulus);
    return as_det_obj(modulus, givelog, sign);
}

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0

SEXP LU_validate(SEXP obj)
{
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (!isReal(x))
	return mkString(_("'x' slot is not of type \"double\""));
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (XLENGTH(x) != pdim[0] * (double) pdim[1])
	return mkString(_("length of 'x' slot is not prod(Dim)"));
    return DimNames_validate(obj, pdim);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer denseLU_expand() */
#if 0

SEXP LU_expand(SEXP x)
{
    const char *nms[] = {"L", "U", "P", ""};
    // x[,] is  m x n    (using LAPACK dgetrf notation)
    SEXP L, U, P, val = PROTECT(Rf_mkNamed(VECSXP, nms)),
	lux = GET_SLOT(x, Matrix_xSym),
	dd = GET_SLOT(x, Matrix_DimSym);
    int *iperm, *perm, *pivot = INTEGER(GET_SLOT(x, Matrix_permSym)),
	*dim = INTEGER(dd), m = dim[0], n = dim[1], nn = m, i;
    size_t m_ = (size_t) m; // to prevent integer (multiplication) overflow
    Rboolean is_sq = (n == m), L_is_tri = TRUE, U_is_tri = TRUE;

    // nn :=  min(n,m) ==  length(pivot[])
    if(!is_sq) {
	if(n < m) { // "long"
	    nn = n;
	    L_is_tri = FALSE;
	} else { // m < n : "wide"
	    U_is_tri = FALSE;
	}
    }

    SET_VECTOR_ELT(val, 0, NEW_OBJECT_OF_CLASS(L_is_tri ? "dtrMatrix":"dgeMatrix"));
    SET_VECTOR_ELT(val, 1, NEW_OBJECT_OF_CLASS(U_is_tri ? "dtrMatrix":"dgeMatrix"));
    SET_VECTOR_ELT(val, 2, NEW_OBJECT_OF_CLASS("pMatrix"));
    L = VECTOR_ELT(val, 0);
    U = VECTOR_ELT(val, 1);
    P = VECTOR_ELT(val, 2);
    if(is_sq || !L_is_tri) {
	SET_SLOT(L, Matrix_xSym, duplicate(lux));
	SET_SLOT(L, Matrix_DimSym, duplicate(dd));
    } else { // !is_sq && L_is_tri -- m < n -- "wide" -- L is  m x m
	size_t m2 = m_ * m;
	double *Lx = REAL(ALLOC_SLOT(L, Matrix_xSym, REALSXP, m2));
	int *dL = INTEGER(ALLOC_SLOT(L, Matrix_DimSym, INTSXP, 2));
	dL[0] = dL[1] = m;
	// fill lower-diagonal (non-{0,1}) part -- remainder by ddense_unpacked_make_*() below:
	Memcpy(Lx, REAL(lux), m2);
    }
    if(is_sq || !U_is_tri) {
	SET_SLOT(U, Matrix_xSym, duplicate(lux));
	SET_SLOT(U, Matrix_DimSym, duplicate(dd));
    } else { // !is_sq && U_is_tri -- m > n -- "long" -- U is  n x n
	double *Ux = REAL(ALLOC_SLOT(U, Matrix_xSym, REALSXP, ((size_t) n) * n)),
	       *xx = REAL(lux);
	int *dU = INTEGER(ALLOC_SLOT(U, Matrix_DimSym, INTSXP, 2));
	dU[0] = dU[1] = n;
	/* fill upper-diagonal (non-0) part -- remainder by ddense_unpacked_make_*() below:
	 * this is more complicated than in the L case, as the x / lux part we need
	 * is  *not*  continguous:  Memcpy(Ux, REAL(lux), n * n); -- is  WRONG */
	for (size_t j = 0; j < n; j++) {
	    Memcpy(Ux+j*n, xx+j*m, j+1);
	    // for (i = 0; i <= j; i++)
	    //   Ux[i + j*n] = xx[i + j*m];
	}
    }
    if(L_is_tri) {
	SET_SLOT(L, Matrix_uploSym, mkString("L"));
	SET_SLOT(L, Matrix_diagSym, mkString("U"));
    }
    // fill the upper right part with 0  *and* the diagonal with 1
    ddense_unpacked_make_triangular(REAL(GET_SLOT(L, Matrix_xSym)),
				    m, (is_sq || !L_is_tri) ? n : m, 'L', 'U');

    if(U_is_tri) {
	SET_SLOT(U, Matrix_uploSym, mkString("U"));
	SET_SLOT(U, Matrix_diagSym, mkString("N"));
	
    }
    // fill the lower left part with 0
    ddense_unpacked_make_triangular(REAL(GET_SLOT(U, Matrix_xSym)),
				    (is_sq || !U_is_tri) ? m : n, n, 'U', 'N');
    
    SET_SLOT(P, Matrix_DimSym, duplicate(dd));
    if(!is_sq) // m != n -- P is  m x m
	INTEGER(GET_SLOT(P, Matrix_DimSym))[1] = m;
    perm = INTEGER(ALLOC_SLOT(P, Matrix_permSym, INTSXP, m));
    Matrix_Calloc(iperm, m, int);

    for (i = 0; i < m; i++) iperm[i] = i + 1; /* initialize permutation*/
    for (i = 0; i < nn; i++) {	/* generate inverse permutation */
	int newp = pivot[i] - 1;
	if (newp != i) { // swap
	    int tmp = iperm[i]; iperm[i] = iperm[newp]; iperm[newp] = tmp;
	}
    }
    // invert the inverse
    for (i = 0; i < m; i++) perm[iperm[i] - 1] = i + 1;

    Matrix_Free(iperm, m);
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
