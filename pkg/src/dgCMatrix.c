#include "dgCMatrix.h"

/* for Csparse_transpose() : */
#include "Csparse.h"
#include "chm_common.h"

/* FIXME -- we "forget" about dimnames almost everywhere : */

/* for dgCMatrix  _and_ lgCMatrix and others  (but *not*  ngC...) : */
SEXP xCMatrix_validate(SEXP x)
{
    /* Almost everything now in Csparse_validate ( ./Csparse.c )
     * *but* the checking of the 'x' slot : */
    if (length(GET_SLOT(x, Matrix_iSym)) !=
	length(GET_SLOT(x, Matrix_xSym)))
	return mkString(_("lengths of slots 'i' and 'x' must match"));

    return ScalarLogical(1);
}

/* for dgRMatrix  _and_ lgRMatrix and others  (but *not*  ngC...) : */
SEXP xRMatrix_validate(SEXP x)
{
    /* Almost everything now in Rsparse_validate ( ./Csparse.c )
     * *but* the checking of the 'x' slot : */
    if (length(GET_SLOT(x, Matrix_jSym)) !=
	length(GET_SLOT(x, Matrix_xSym)))
	return mkString(_("lengths of slots 'j' and 'x' must match"));

    return ScalarLogical(1);
}

/* This and the following R_to_CMatrix() lead to memory-not-mapped seg.faults
 * only with {32bit + R-devel + enable-R-shlib} -- no idea why */
SEXP compressed_to_TMatrix(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    /* however, for Csparse, we now effectively use the cholmod-based
     * Csparse_to_Tsparse() in ./Csparse.c ; maybe should simply write
     * an  as_cholmod_Rsparse() function and then do "as there" ...*/
    SEXP indSym = col ? Matrix_iSym : Matrix_jSym,
	ans,	indP = GET_SLOT(x, indSym),
	pP = GET_SLOT(x, Matrix_pSym);
    int npt = length(pP) - 1;
    char *ncl = strdup(class_P(x));
    char *valid[] = {
	"dgCMatrix", "dsCMatrix", "dtCMatrix", /* 0: 0:2 */
	"lgCMatrix", "lsCMatrix", "ltCMatrix", /* 1: 3:5 */
	"ngCMatrix", "nsCMatrix", "ntCMatrix", /* 2: 6:8 */
	"zgCMatrix", "zsCMatrix", "ztCMatrix", /* 3: 9:11 */

	"dgRMatrix", "dsRMatrix", "dtRMatrix", /* 4: 12:14 */
	"lgRMatrix", "lsRMatrix", "ltRMatrix", /* 5: 15:17 */
	"ngRMatrix", "nsRMatrix", "ntRMatrix", /* 6: 18:20 */
	"zgRMatrix", "zsRMatrix", "ztRMatrix", /* 7: 21:23 */
	""};
    int ctype = Matrix_check_class(ncl, valid);

    if (ctype < 0)
	error(_("invalid class(x) '%s' in compressed_to_TMatrix(x)"), ncl);

    /* replace 'C' or 'R' with 'T' :*/
    ncl[2] = 'T';
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(ncl)));

    slot_dup(ans, x, Matrix_DimSym);
    if((ctype / 3) % 4 != 2) /* not n..Matrix */
	slot_dup(ans, x, Matrix_xSym);
    if(ctype % 3) { /* s(ymmetric) or t(riangular) : */
	slot_dup(ans, x, Matrix_uploSym);
	if(ctype % 3 == 2) /* t(riangular) : */
	    slot_dup(ans, x, Matrix_diagSym);
    }
    SET_DimNames(ans, x);
    SET_SLOT(ans, indSym, duplicate(indP));
    expand_cmprPt(npt, INTEGER(pP),
		  INTEGER(ALLOC_SLOT(ans, col ? Matrix_jSym : Matrix_iSym,
				     INTSXP, length(indP))));
    free(ncl);
    UNPROTECT(1);
    return ans;
}

SEXP R_to_CMatrix(SEXP x)
{
    SEXP ans, tri = PROTECT(allocVector(LGLSXP, 1));
    char *ncl = strdup(class_P(x));
    char *valid[] = {
	"dgRMatrix", "dsRMatrix", "dtRMatrix",
	"lgRMatrix", "lsRMatrix", "ltRMatrix",
	"ngRMatrix", "nsRMatrix", "ntRMatrix",
	"zgRMatrix", "zsRMatrix", "ztRMatrix",
	""};
    int ctype = Matrix_check_class(ncl, valid);
    int *x_dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), *a_dims;

    if (ctype < 0)
	error(_("invalid class(x) '%s' in R_to_CMatrix(x)"), ncl);

    /* replace 'R' with 'C' : */
    ncl[2] = 'C';
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(ncl)));

    a_dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    /* reversed dim() since we will transpose: */
    a_dims[0] = x_dims[1];
    a_dims[1] = x_dims[0];

    /* triangular: */ LOGICAL(tri)[0] = 0;
    if((ctype / 3) != 2) /* not n..Matrix */
	slot_dup(ans, x, Matrix_xSym);
    if(ctype % 3) { /* s(ymmetric) or t(riangular) : */
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((*uplo_P(x) == 'U') ? "L" : "U"));
	if(ctype % 3 == 2) { /* t(riangular) : */
	    LOGICAL(tri)[0] = 1;
	    slot_dup(ans, x, Matrix_diagSym);
	}
    }
    SET_SLOT(ans, Matrix_iSym, duplicate(GET_SLOT(x, Matrix_jSym)));
    slot_dup(ans, x, Matrix_pSym);
    ans = Csparse_transpose(ans, tri);
    SET_DimNames(ans, x);
    free(ncl);
    UNPROTECT(2);
    return ans;
}

SEXP compressed_non_0_ij(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    SEXP ans, indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP indP = GET_SLOT(x, indSym),
	pP = GET_SLOT(x, Matrix_pSym);
    int n_el = length(indP), i, *ij;

    ij = INTEGER(ans = PROTECT(allocMatrix(INTSXP, n_el, 2)));
    /* expand the compressed margin to 'i' or 'j' : */
    expand_cmprPt(length(pP) - 1, INTEGER(pP), &ij[col ? n_el : 0]);
    /* and copy the other one: */
    if (col)
	for(i = 0; i < n_el; i++)
	    ij[i] = INTEGER(indP)[i];
    else /* row compressed */
	for(i = 0; i < n_el; i++)
	    ij[i + n_el] = INTEGER(indP)[i];

    UNPROTECT(1);
    return ans;
}

SEXP dgCMatrix_lusol(SEXP x, SEXP y)
{
    SEXP ycp = PROTECT(duplicate(y));
    CSP xc = AS_CSP(x);
    R_CheckStack();

    if (xc->m != xc->n || xc->m <= 0)
	error(_("dgCMatrix_lusol requires a square, non-empty matrix"));
    if (!isReal(ycp) || LENGTH(ycp) != xc->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    if (!cs_lusol(/*order*/ 1, xc, REAL(ycp), /*tol*/ 1e-7))
	error(_("cs_lusol failed"));

    UNPROTECT(1);
    return ycp;
}

SEXP dgCMatrix_qrsol(SEXP x, SEXP y)
{
    SEXP ycp = PROTECT(duplicate(y));
    CSP xc = AS_CSP(x);
    R_CheckStack();

    if (xc->m < xc->n || xc->n <= 0)
	error(_("dgCMatrix_qrsol requires a 'tall' rectangular matrix"));
    if (!isReal(ycp) || LENGTH(ycp) != xc->m)
	error(_("Dimensions of system to be solved are inconsistent"));
    if (!cs_qrsol(/*order*/ 1, xc, REAL(ycp)))
	error(_("cs_qrsol failed"));

    UNPROTECT(1);
    return ycp;
}

/* Modified version of Tim Davis's cs_qr_mex.c file for MATLAB */
SEXP dgCMatrix_QR(SEXP Ap, SEXP order)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sparseQR")));
    CSP A = AS_CSP(Ap), D;
    css *S;
    csn *N;
    int m = A->m, n = A->n, ord = asLogical(order) ? 3 : 0, *p;
    int *dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    R_CheckStack();

    if (m < n) error("A must have # rows >= # columns") ;
    dims[0] = m; dims[1] = n;
    S = cs_sqr(ord, A, 1);	/* symbolic QR ordering & analysis*/
    if (!S) error("cs_sqr failed");
    N = cs_qr(A, S);		/* numeric QR factorization */
    if (!N) error("cs_qr failed") ;
    cs_dropzeros(N->L);		/* drop zeros from V and sort */
    D = cs_transpose(N->L, 1);
    cs_spfree(N->L);
    N->L = cs_transpose(D, 1);
    cs_spfree(D);
    cs_dropzeros(N->U);		/* drop zeros from R and sort */
    D = cs_transpose(N->U, 1);
    cs_spfree(N->U) ;
    N->U = cs_transpose(D, 1);
    cs_spfree(D);
    m = N->L->m;		/* m may be larger now */
    p = cs_pinv(S->pinv, m);	/* p = pinv' */
    SET_SLOT(ans, install("V"),
	     Matrix_cs_to_SEXP(N->L, "dgCMatrix", 0));
    Memcpy(REAL(ALLOC_SLOT(ans, install("beta"),
			   REALSXP, n)), N->B, n);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym,
			      INTSXP, m)), p, m);
    SET_SLOT(ans, install("R"),
	     Matrix_cs_to_SEXP(N->U, "dgCMatrix", 0));
    if (ord)
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("q"),
				  INTSXP, n)), S->q, n);
    else
	ALLOC_SLOT(ans, install("q"), INTSXP, 0);
    cs_nfree(N);
    cs_sfree(S);
    cs_free(p);
    UNPROTECT(1);
    return ans;
}

/* Modified version of Tim Davis's cs_lu_mex.c file for MATLAB */
SEXP dgCMatrix_LU(SEXP Ap, SEXP orderp, SEXP tolp)
{
    /* Is currently only called as  .Call(dgCMatrix_LU, x, TRUE, 1)) */
    SEXP ans = get_factors(Ap, "LU");
    CSP A = AS_CSP(Ap), D;
    css *S;
    csn *N;
    int n, order = asInteger(orderp), *p;
    double tol = asReal(tolp);
    R_CheckStack();

    /* FIXME: dgCMatrix_LU should check ans for consistency in
     * permutation type with the requested value - Should have two
     * classes or two different names in the factors list for LU with
     * permuted columns or not. */

    if (ans != R_NilValue) return ans;
    n = A->n;
    if (A->m != n)
	error("LU decomposition applies only to square matrices");
    if (order) {		/* not using natural order */
	order = (tol == 1) ? 2	/* amd(S'*S) w/dense rows or I */
	    : 1;		/* amd (A+A'), or natural */
    }
    S = cs_sqr (order, A, 0) ;	/* symbolic ordering, no QR bound */
    N = cs_lu (A, S, tol) ;	/* numeric factorization */
    if (!N) {
	/*WAS: error ("cs_lu failed (singular, or out of memory)") ; */
	return R_NilValue;
	/* and the caller can warn() or stop() ... */
    }
    cs_dropzeros (N->L) ;	/* drop zeros from L and sort it */
    D = cs_transpose (N->L, 1) ;
    cs_spfree (N->L) ;
    N->L = cs_transpose (D, 1) ;
    cs_spfree (D) ;
    cs_dropzeros (N->U) ;	/* drop zeros from U and sort it */
    D = cs_transpose (N->U, 1) ;
    cs_spfree (N->U) ;
    N->U = cs_transpose (D, 1) ;
    cs_spfree (D) ;
    p = cs_pinv (N->pinv, n) ;	/* p=pinv' */
    ans = PROTECT(NEW_OBJECT(MAKE_CLASS("sparseLU")));
    SET_SLOT(ans, install("L"),
	     Matrix_cs_to_SEXP(N->L, "dtCMatrix", 0));
    SET_SLOT(ans, install("U"),
	     Matrix_cs_to_SEXP(N->U, "dtCMatrix", 0));
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym,
			      INTSXP, n)), p, n);
    if (order)
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("q"),
				  INTSXP, n)), S->q, n);
    cs_nfree(N);
    cs_sfree(S);
    cs_free(p);
    UNPROTECT(1);
    return set_factors(Ap, ans, "LU");
}

SEXP dgCMatrix_matrix_solve(SEXP Ap, SEXP b)
{
    /* b is dense or NULL [ <--> solve(A) */
    SEXP lu = dgCMatrix_LU(Ap, ScalarLogical(1), ScalarReal(1));
    SEXP qslot = GET_SLOT(lu, install("q"));
    CSP L = AS_CSP(GET_SLOT(lu, install("L"))),
	U = AS_CSP(GET_SLOT(lu, install("U")));
    SEXP ans = PROTECT( !isNull(b) ? dup_mMatrix_as_dgeMatrix(b)
			: new_dgeMatrix(U->n, U->n));
    int *bdims = INTEGER(GET_SLOT(ans, Matrix_DimSym));
    int j, n = bdims[0], nrhs = bdims[1];
    int *p = INTEGER(GET_SLOT(lu, Matrix_pSym)),
	*q = LENGTH(qslot) ? INTEGER(qslot) : (int *) NULL;
    double *ax = REAL(GET_SLOT(ans, Matrix_xSym)),
	*x = Alloca(n, double);
    R_CheckStack();

    if (U->n != n || nrhs < 1 || n < 1)
	error(_("Dimensions of system to be solved are inconsistent"));
    for (j = 0; j < nrhs; j++) {
	if(!isNull(b))
	    cs_pvec(p, ax + j * n, x, n);  /* x = b(p) */
	else { /* solve(A): (RHS) B = I_n,  hence  b = e_j (j-th unit vector) */
	    int i;
	    for(i=0; i < n; i++) x[i] = (p[i] == j) ? 1. : 0.;
	}
	cs_lsolve(L, x);	       /* x = L\x */
	cs_usolve(U, x);	       /* x = U\x */
	if (q)			       /* b(q) = x */
	    cs_ipvec(q, x, ax + j * n, n);
	else
	    Memcpy(ax + j * n, x, n);
    }
    UNPROTECT(1);
    return ans;
}

SEXP dgCMatrix_cholsol(SEXP x, SEXP y)
{
    CHM_SP cx = AS_CHM_SP(x);
    CHM_FR L;
    CHM_DN cy = AS_CHM_DN(y), rhs, cAns;
    double one[] = {1,0}, zero[] = {0,0};
    SEXP ans = PROTECT(allocVector(VECSXP, 3));
    R_CheckStack();

    if (cx->ncol < cx->nrow || cx->ncol <= 0)
	error(_("dgCMatrix_cholsol requires a 'short, wide' rectangular matrix"));
    if (cy->nrow != cx->ncol)
	error(_("Dimensions of system to be solved are inconsistent"));
    rhs = cholmod_allocate_dense(cx->nrow, 1, cx->nrow, CHOLMOD_REAL, &c);
    if (!(cholmod_sdmult(cx, 0 /* trans */, one, zero, cy, rhs, &c)))
	error(_("cholmod_sdmult error"));
    L = cholmod_analyze(cx, &c);
    if (!cholmod_factorize(cx, L, &c))
	error(_("cholmod_factorize failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
/* FIXME: Do this in stages so an "effects" vector can be calculated */
    if (!(cAns = cholmod_solve(CHOLMOD_A, L, rhs, &c)))
	error(_("cholmod_solve (CHOLMOD_A) failed: status %d, minor %d from ncol %d"),
	      c.status, L->minor, L->n);
    SET_VECTOR_ELT(ans, 0, chm_factor_to_SEXP(L, 0));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, cx->nrow));
    Memcpy(REAL(VECTOR_ELT(ans, 1)), (double*)(cAns->x), cx->nrow);
/* FIXME: Change this when the "effects" vector is available */
    SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, cx->nrow));
    Memcpy(REAL(VECTOR_ELT(ans, 1)), (double*)(rhs->x), cx->nrow);

    cholmod_free_factor(&L, &c);
    cholmod_free_dense(&rhs, &c);
    cholmod_free_dense(&cAns, &c);
    UNPROTECT(1);
    return ans;
}


/* Define all of
 *  dgCMatrix_colSums(....)
 *  igCMatrix_colSums(....)
 *  lgCMatrix_colSums_d(....)
 *  lgCMatrix_colSums_i(....)
 *  ngCMatrix_colSums_d(....)
 *  ngCMatrix_colSums_i(....)
 */
#define _dgC_
#include "t_gCMatrix_colSums.c"

#define _igC_
#include "t_gCMatrix_colSums.c"

#define _lgC_
#include "t_gCMatrix_colSums.c"

#define _ngC_
#include "t_gCMatrix_colSums.c"

#define _lgC_mn
#include "t_gCMatrix_colSums.c"

#define _ngC_mn
#include "t_gCMatrix_colSums.c"


SEXP lgCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means)
{
    if(asLogical(means)) /* ==> result will be "double" / "dsparseVector" */
	return lgCMatrix_colSums_d(x, NArm, spRes, trans, means);
    else
	return lgCMatrix_colSums_i(x, NArm, spRes, trans, means);
}

SEXP ngCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means)
{
    if(asLogical(means)) /* ==> result will be "double" / "dsparseVector" */
	return ngCMatrix_colSums_d(x, NArm, spRes, trans, means);
    else
	return ngCMatrix_colSums_i(x, NArm, spRes, trans, means);
}
