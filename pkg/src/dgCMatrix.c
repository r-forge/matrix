#include "dgCMatrix.h"

#include "chm_common.h"

/* FIXME -- we "forget" about dimnames almost everywhere : */

SEXP dgCMatrix_validate(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym),
	xslot = GET_SLOT(x, Matrix_xSym);
    int j,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow = dims[0],
	ncol = dims[1],
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    if (length(islot) != length(xslot))
	return mkString(_("lengths of slots i and x must match"));
    if (length(pslot) != ncol + 1)
	return mkString(_("slot p must have length ncol + 1"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(islot) != xp[ncol])
	return mkString(_("last element of slot p must match length of slots i and x"));
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j+1])
	    return mkString(_("slot p must be non-decreasing"));
    }
    for (j = 0; j < length(islot); j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return mkString(_("all row indices must be between 0 and nrow-1"));
    }
    if (csc_unsorted_columns(ncol, xp, xi))
	csc_sort_columns(ncol, xp, xi, REAL(xslot));

    return ScalarLogical(1);
}

#if 0 		/* more general methods now defined in ./Csparse.c */
SEXP csc_crossprod(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))), tmp;
    int *xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

    int j, *iVal, ncol = length(pslot) - 1, maxnz, nnz = 0, *pVal;
    double *xVal;

    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(ans, Matrix_DimSym, allocVector(INTSXP, 2));
    SET_SLOT(ans, Matrix_uploSym, mkString("L"));
    maxnz = (ncol * (ncol + 1))/2;
    iVal = Calloc(maxnz, int); xVal = Calloc(maxnz, double);
    SET_SLOT(ans, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    tmp = GET_SLOT(ans, Matrix_pSym);
    pVal = INTEGER(tmp);
    for (j = 0; j < ncol; j++) {
	pVal[j] = nnz;
	if (xp[j] < xp[j+1]) {	/* column j contains some non-zeros */
	    int ind, jj;
	    double accum = 0.;
				/* diagonal elements */
	    for (ind = xp[j]; ind < xp[j+1]; ind++)
		accum += xx[ind] * xx[ind];
	    iVal[nnz] = j;
	    xVal[nnz] = accum;
	    nnz++;
				/* off-diagonals (lower triangle only) */
	    for (jj = j+1; jj < ncol; jj++) {
		int ind2;

		ind = xp[j];
		ind2 = xp[jj];
		accum = 0.;
		while (ind < xp[j+1] && ind2 < xp[jj+1]) {
		    if (xi[ind] < xi[ind2]) ind++;
		    else {
			if (xi[ind] > xi[ind2]) ind2++;
			else {
			    accum += xx[ind] * xx[ind2];
			    ind++; ind2++;
			}
		    }
		}
		if (accum != 0.) {
		    iVal[nnz] = jj;
		    xVal[nnz] = accum;
		    nnz++;
		}
	    }
	}
    }
    pVal[ncol] = nnz;

    SET_SLOT(ans, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(ans, Matrix_iSym)), iVal, nnz);
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(ans, Matrix_xSym)), xVal, nnz);
    Free(iVal); Free(xVal); UNPROTECT(1);
    return dgCMatrix_set_Dim(ans, ncol);
}

SEXP csc_tcrossprod(SEXP x)
{
    cholmod_sparse *chx = as_cholmod_sparse(x);
    cholmod_sparse *cha = cholmod_aat(chx, (int *) NULL, 0, 1, &c);
    cholmod_sparse *chas = cholmod_copy(cha, 1, 1, &c);

    Free(chx);
    cholmod_free_sparse(&cha, &c);
    return chm_sparse_to_SEXP(chas, 1, 0, "", R_NilValue);
}

SEXP csc_matrix_crossprod(SEXP x, SEXP y, SEXP classed)
{
    int cl = asLogical(classed);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *xdims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*ydims = INTEGER(cl ? GET_SLOT(y, Matrix_DimSym) :
			 getAttrib(y, R_DimSymbol)),
	*vdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
    int *xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	*xp = INTEGER(GET_SLOT(x, Matrix_pSym));
    int j, k = xdims[0], m = xdims[1], n = ydims[1];
    double *vx, *xx = REAL(GET_SLOT(x, Matrix_xSym)),
	*yx = REAL(cl ? GET_SLOT(y, Matrix_xSym) : y);

    if (!cl && !(isMatrix(y) && isReal(y)))
	error(_("y must be a numeric matrix"));
    if (ydims[0] != k)
	error(_("x and y must have the same number of rows"));
    if (m < 1 || n < 1 || k < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    vdims[0] = m; vdims[1] = n;
    vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n));
    for (j = 0; j < n; j++) {
	int i; double *ypt = yx + j * k;
	for(i = 0; i < m; i++) {
	    int ii; double accum = 0.;
	    for (ii = xp[i]; ii < xp[i+1]; ii++) {
		accum += xx[ii] * ypt[xi[ii]];
	    }
	    vx[i + j * m] = accum;
	}
    }
    UNPROTECT(1);
    return val;
}
#endif	/* more general methods */

SEXP compressed_to_dgTMatrix(SEXP x, SEXP colP)
{
    int col = asLogical(colP); /* 1 if "C"olumn compressed;  0 if "R"ow */
    SEXP indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgTMatrix"))),
	indP = GET_SLOT(x, indSym),
	pP = GET_SLOT(x, Matrix_pSym);
    int npt = length(pP) - 1;

    SET_SLOT(ans, Matrix_DimSym, duplicate(GET_SLOT(x, Matrix_DimSym)));
    SET_SLOT(ans, Matrix_xSym,  duplicate(GET_SLOT(x, Matrix_xSym)));
    SET_SLOT(ans, indSym, duplicate(indP));
    expand_cmprPt(npt, INTEGER(pP),
		  INTEGER(ALLOC_SLOT(ans, col ? Matrix_jSym : Matrix_iSym,
				     INTSXP, length(indP))));
    UNPROTECT(1);
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

#if 0 				/* more general methods in ./Csparse.c */
SEXP csc_to_matrix(SEXP x)
{
    SEXP ans, pslot = GET_SLOT(x, Matrix_pSym);
    int j, ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	*xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *ax;

    ax = REAL(ans = PROTECT(allocMatrix(REALSXP, nrow, ncol)));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0.;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++) {
	    ax[j * nrow + xi[ind]] = xx[ind];
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP csc_to_dgeMatrix(SEXP x)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	Dimslot = GET_SLOT(x, Matrix_DimSym);
    int *dims = INTEGER(Dimslot),
	*xp = INTEGER(GET_SLOT(x, Matrix_pSym)),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *ax;
    int j, nrow = dims[0], ncol = dims[1];

    SET_SLOT(ans, Matrix_DimSym, duplicate(Dimslot));
    SET_SLOT(ans, Matrix_xSym, allocVector(REALSXP, nrow*ncol));
    SET_SLOT(ans, Matrix_factorSym, allocVector(VECSXP, 0));
    ax = REAL(GET_SLOT(ans, Matrix_xSym));
    for (j = 0; j < (nrow * ncol); j++) ax[j] = 0.;
    for (j = 0; j < ncol; j++) {
	int ind;
	for (ind = xp[j]; ind < xp[j+1]; ind++) {
	    ax[j * nrow + xi[ind]] = xx[ind];
	}
    }
    UNPROTECT(1);
    return ans;
}

SEXP double_to_csc(double *a, int *dim_a)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    int j, maxnz, nrow, ncol, nnz, *vp, *vi;
    double *vx;

    nrow = dim_a[0]; ncol = dim_a[1];
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    SET_SLOT(val, Matrix_pSym, allocVector(INTSXP, ncol + 1));
    vp = INTEGER(GET_SLOT(val, Matrix_pSym));
    maxnz = nrow * ncol;
    vi = Calloc(maxnz, int); vx = Calloc(maxnz, double);
    nnz = 0;
    for (j = 0; j < ncol; j++) {
	int i;
	vp[j] = nnz;
	for (i = 0; i < nrow; i++) {
	    double val = a[i + j * nrow];
	    if (val != 0.) {
		vi[nnz] = i;
		vx[nnz] = val;
		nnz++;
	    }
	}
    }
    vp[ncol] = nnz;
    SET_SLOT(val, Matrix_iSym, allocVector(INTSXP, nnz));
    Memcpy(INTEGER(GET_SLOT(val, Matrix_iSym)), vi, nnz);
    SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, nnz));
    Memcpy(REAL(GET_SLOT(val, Matrix_xSym)), vx, nnz);
    Free(vi); Free(vx);
    UNPROTECT(1);
    return dgCMatrix_set_Dim(val, nrow);
}
#endif	/* more general methods */

#if 0 				/* use dense_to_Csparse in ./dense.c */
SEXP matrix_to_csc(SEXP A)
{
    if (!(isMatrix(A) && isReal(A)))
	error(_("A must be a numeric matrix"));
    return double_to_csc(REAL(A),
			 INTEGER(getAttrib(A, R_DimSymbol)));
}

SEXP dgeMatrix_to_csc(SEXP x)
{
    return double_to_csc(   REAL(GET_SLOT(x, Matrix_xSym)),
			 INTEGER(GET_SLOT(x, Matrix_DimSym)));
}
#endif	/* use dense_to_Csparse */

#if 0 				/* more general method in Tsparse */
SEXP dgTMatrix_to_csc(SEXP dgTMatrix)
{
    SEXP Tisl = GET_SLOT(dgTMatrix, Matrix_iSym);
    int *Ti = INTEGER(Tisl),
	*Tj = INTEGER(GET_SLOT(dgTMatrix, Matrix_jSym)),
	i, nrow, ncol,
	nz = length(Tisl);

    nrow = ncol = -1;
    for(i = 0; i < nz; i++) {
	if (Ti[i] > nrow) nrow = Ti[i];
	if (Tj[i] > ncol) ncol = Tj[i];
    }
    return triple_as_SEXP(nrow + 1, ncol + 1, nz, Ti, Tj,
			  REAL(GET_SLOT(dgTMatrix, Matrix_xSym)),
			  "dgCMatrix");
}
#endif 				/* more general method */

#if 0 				/* use Csparse_band and R code */
SEXP csc_getDiag(SEXP x)
{
    SEXP pslot = GET_SLOT(x, Matrix_pSym), ans;
    int *xp = INTEGER(pslot),
	*xi = INTEGER(GET_SLOT(x, Matrix_iSym)),
	j,
	ncol = length(pslot) - 1,
	nrow = INTEGER(GET_SLOT(x, Matrix_DimSym))[0],
	ndiag;
    double *xx = REAL(GET_SLOT(x, Matrix_xSym)), *diag;

    ndiag = (nrow < ncol) ? nrow : ncol;
    ans = PROTECT(allocVector(REALSXP, ndiag));
    diag = REAL(ans);
    for (j = 0; j < ndiag; j++) {
	int ind;
	diag[j] = 0.;
	for (ind = xp[j]; ind < xp[j+1]; ind++) {
	    if (xi[ind] == j) diag[j] = xx[ind];
	}
    }
    UNPROTECT(1);
    return ans;
}
#endif	/* use Csparse_band */

#if 0 				/* more general method in Csparse.c */
SEXP csc_transpose(SEXP x)
{
    cholmod_sparse *chx = as_cholmod_sparse(x);
    SEXP ans =
	chm_sparse_to_SEXP(cholmod_transpose(chx, 1, &c), 1, 0, "", R_NilValue);
    Free(chx);
    return ans;
}
#endif	/* more general method */

SEXP csc_matrix_mm(SEXP a, SEXP b, SEXP classed, SEXP right)
{
    int cl = asLogical(classed), rt = asLogical(right);
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*ai = INTEGER(GET_SLOT(a, Matrix_iSym)),
	*ap = INTEGER(GET_SLOT(a, Matrix_pSym)),
	*bdims = INTEGER(cl ? GET_SLOT(b, Matrix_DimSym) :
			 getAttrib(b, R_DimSymbol)),
	*cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2)),
	chk, ione = 1, j, jj, k, m, n;
    double *ax = REAL(GET_SLOT(a, Matrix_xSym)),
	*bx = REAL(cl ? GET_SLOT(b, Matrix_xSym) : b), *cx;

    if (rt) {
	m = bdims[0]; n = adims[1]; k = bdims[1]; chk = adims[0];
    } else {
	m = adims[0]; n = bdims[1]; k = adims[1]; chk = bdims[0];
    }
    if (chk != k)
	error(_("Matrices are not conformable for multiplication"));
    if (m < 1 || n < 1 || k < 1)
	error(_("Matrices with zero extents cannot be multiplied"));
    cx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n));
    AZERO(cx, m * n); /* zero the accumulators */
    for (j = 0; j < n; j++) { /* across columns of c */
	if (rt) {
	    int kk, k2 = ap[j + 1];
	    for (kk = ap[j]; kk < k2; kk++) {
		F77_CALL(daxpy)(&m, &ax[kk], &bx[ai[kk]*m],
				&ione, &cx[j*m], &ione);
	    }
	} else {
	    double *ccol = cx + j * m,
		*bcol = bx + j * k;

	    for (jj = 0; jj < k; jj++) { /* across columns of a */
		int kk, k2 = ap[jj + 1];
		for (kk = ap[jj]; kk < k2; kk++) {
		    ccol[ai[kk]] += ax[kk] * bcol[jj];
		}
	    }
	}
    }
    cdims[0] = m; cdims[1] = n;
    UNPROTECT(1);
    return val;
}

#if 0				/* no longer used */
SEXP csc_col_permute(SEXP x, SEXP perm)
{
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), tmp;
    int *iperm, *prm, *vi, *vp, *xi, *xp, j, k, ncol, pos;
    double *vx, *xx;

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    tmp = GET_SLOT(x, Matrix_DimSym);
    SET_SLOT(val, Matrix_DimSym, duplicate(tmp));
    ncol = INTEGER(tmp)[1];
    if (!(isInteger(perm) && length(perm) == ncol))
	error(_("perm must be an integer vector of length %d"),
	      ncol);
    prm = INTEGER(perm);
    if (!R_ldl_valid_perm(ncol, prm))
	error(_("perm is not a valid 0-based permutation"));
    iperm = Calloc(ncol, int);
    for (j = 0; j < ncol; j++) iperm[prm[j]] = j;
    tmp = GET_SLOT(x, Matrix_pSym);
    xp = INTEGER(tmp);
    SET_SLOT(val, Matrix_pSym, duplicate(tmp));
    vp = INTEGER(GET_SLOT(val, Matrix_pSym));
    tmp = GET_SLOT(x, Matrix_iSym);
    xi = INTEGER(tmp);
    SET_SLOT(val, Matrix_iSym, duplicate(tmp));
    vi = INTEGER(GET_SLOT(val, Matrix_iSym));
    tmp = GET_SLOT(x, Matrix_xSym);
    xx = REAL(tmp);
    SET_SLOT(val, Matrix_xSym, duplicate(tmp));
    vx = REAL(GET_SLOT(val, Matrix_xSym));

    pos = vp[0] = 0;
    for (j = 0; j < ncol; j++) {
	int jj = iperm[j];
	int j1 = xp[jj], j2 = xp[jj+1];
	vp[j + 1] = vp[j] + (j2 - j1);
	for (k = j1; k < j2; k++) {
	    vi[pos] = xi[k];
	    vx[pos] = xx[k];
	    pos++;
	}
    }
    Free(iperm);
    UNPROTECT(1);
    return val;
}
#endif	/* no longer used */
