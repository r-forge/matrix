#include "factorizations.h"

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0

SEXP LU_validate(SEXP obj)
{
    /* NB: 'Dim' already checked by MatrixFactorization_validate() */
    
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (!isReal(x))
	return mkString(_("'x' slot is not of type \"double\""));
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (XLENGTH(x) != pdim[0] * (double) pdim[1])
	return mkString(_("length of 'x' slot is not prod(Dim)"));
    return DimNames_validate(obj, pdim);
}

#endif /* MJ */

SEXP denseLU_expand(SEXP obj)
{
    /* A = P L U, where ... 
       
       A -> [m,n]
       P -> [m,m], permutation
       L -> if m <= n then [m,n] else [m,m], lower trapezoidal, unit diagonal
       U -> if m >= n then [m,n] else [n,n], upper trapezoidal
       
       square L,U given as dtrMatrix with appropriate 'uplo', 'diag' slots,
       non-square L,U given as dgeMatrix
    */
    
    const char *nms[] = {"L", "U", "P", ""};
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
	SET_VECTOR_ELT(res, 0, L);
	SET_VECTOR_ELT(res, 1, U);
	UNPROTECT(4); /* diag, uplo, U, L */
    } else {
	SEXP G = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
	    T = PROTECT(NEW_OBJECT_OF_CLASS("dtrMatrix")),
	    y = PROTECT(allocVector(REALSXP, (R_xlen_t) r * r));
	REPROTECT(x = duplicate(x), pidB);
	double *px = REAL(x), *py = REAL(y);
	int whichT = (m < n) ? 0 : 1;
	
	SET_SLOT(G, Matrix_DimSym, dim);
	REPROTECT(dim = allocVector(INTSXP, 2), pidA);
	pdim = INTEGER(dim);
	pdim[0] = pdim[1] = r;
	SET_SLOT(T, Matrix_DimSym, dim);
	REPROTECT(dim = allocVector(INTSXP, 2), pidA);
	pdim = INTEGER(dim);
	pdim[0] = pdim[1] = m;
	SET_SLOT(P, Matrix_DimSym, dim);
	
	if (whichT == 0) {
            /* G is upper trapezoidal, T is unit lower triangular */
	    SEXP uplo = PROTECT(mkString("L")),
		diag = PROTECT(mkString("U"));
	    SET_SLOT(T, Matrix_uploSym, uplo);
	    SET_SLOT(T, Matrix_diagSym, diag);
	    UNPROTECT(2); /* diag, uplo */

	    Memcpy(py, px, (size_t) m * m);
	    ddense_unpacked_make_triangular(px, m, n, 'U', 'N');
	} else {
            /* G is unit lower trapezoidal, T is upper triangular */
	    double *tmp = px;
	    for (j = 0; j < n; ++j, px += m, py += r)
		Memcpy(py, px, j+1);
	    ddense_unpacked_make_triangular(tmp, m, n, 'L', 'U');
	}
	SET_SLOT(G, Matrix_xSym, x);
	SET_SLOT(T, Matrix_xSym, y);
	
	SET_VECTOR_ELT(res, !whichT, G);
	SET_VECTOR_ELT(res,  whichT, T);
	UNPROTECT(3); /* y, T, G */
    }

    SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
	perm = PROTECT(allocVector(INTSXP, m));
    int *ppivot = INTEGER(pivot), *pperm = INTEGER(perm), *pinvperm, pos, tmp;
    Calloc_or_Alloca_TO(pinvperm, m, int);

    /* MJ: inversion step below can be skipped once class indMatrix
           is generalized to include _column_ index matrices ...
	   but may need to be kept anyway for backwards compatibility
    */
    
    for (j = 0; j < m; ++j) /* initialize column permutation */
	pinvperm[j] = j;
    for (j = 0; j < r; ++j) { /* generate column permutation */
	pos = ppivot[j] - 1;
	if (pos != j) {
	    tmp = pinvperm[j];
	    pinvperm[j] = pinvperm[pos];
	    pinvperm[pos] = tmp;
	}
    }
    for (j = 0; j < m; ++j) /* invert column permutation (0->1 based) */
	pperm[pinvperm[j]] = j + 1;
    Free_FROM(pinvperm, m);

    SET_SLOT(P, Matrix_permSym, perm);
    SET_VECTOR_ELT(res, 2, P);
    UNPROTECT(6); /* perm, pivot, x, dim, P, res */
    return res;
}

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
    Calloc_or_Alloca_TO(iperm, m, int);

    for (i = 0; i < m; i++) iperm[i] = i + 1; /* initialize permutation*/
    for (i = 0; i < nn; i++) {	/* generate inverse permutation */
	int newp = pivot[i] - 1;
	if (newp != i) { // swap
	    int tmp = iperm[i]; iperm[i] = iperm[newp]; iperm[newp] = tmp;
	}
    }
    // invert the inverse
    for (i = 0; i < m; i++) perm[iperm[i] - 1] = i + 1;

    Free_FROM(iperm, m);
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
