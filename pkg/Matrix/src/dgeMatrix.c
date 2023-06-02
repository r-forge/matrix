#include <float.h> /* DBL_EPSILON */
#include "dgeMatrix.h"

double get_norm_dge(SEXP obj, const char *typstr)
{
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    R_xlen_t i, nx = XLENGTH(x);
    double *px = REAL(x);
    for (i = 0; i < nx; ++i) {
	if (ISNAN(px[i])) {
	    UNPROTECT(1);
	    return NA_REAL;
	}
    }
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim);
    double norm, *work = NULL;

    if (typstr[0] == 'I')
	work = (double *) R_alloc((size_t) pdim[0], sizeof(double));
    norm = F77_CALL(dlange)(typstr, pdim, pdim + 1, px, pdim,
			    work FCONE);
    
    UNPROTECT(2);
    return norm;
}

SEXP dgeMatrix_norm(SEXP obj, SEXP type)
{
    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_norm_type(CHAR(type));
    double norm = get_norm_dge(obj, typstr);
    UNPROTECT(1);
    return ScalarReal(norm);
}

SEXP dgeMatrix_rcond(SEXP obj, SEXP type)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
    int *pdim = INTEGER(dim);
    if (pdim[0] != pdim[1] || pdim[0] < 1)
	error(_("'rcond' requires a square, nonempty matrix"));

    char typstr[] = {'\0', '\0'};
    PROTECT(type = asChar(type));
    typstr[0] = La_rcond_type(CHAR(type));

    SEXP trf = PROTECT(dgeMatrix_trf_(obj, 0)),
	x = PROTECT(GET_SLOT(trf, Matrix_xSym));
    int info;
    double *px = REAL(x), norm = get_norm_dge(obj, typstr), rcond;
    
    F77_CALL(dgecon)(typstr, pdim, px, pdim, &norm, &rcond,
		     (double *) R_alloc((size_t) 4 * pdim[0], sizeof(double)),
		     (int *) R_alloc((size_t) pdim[0], sizeof(int)),
		     &info FCONE);

    UNPROTECT(4);
    return ScalarReal(rcond);
}

/* MJ: no longer needed ... prefer more general unpackedMatrix_diag_[gs]et() */
#if 0

SEXP dgeMatrix_getDiag(SEXP x)
{
#define geMatrix_getDiag_1					\
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));		\
    int i, m = dims[0], nret = (m < dims[1]) ? m : dims[1];	\
    R_xlen_t m1 = m + 1;					\
    SEXP x_x = GET_SLOT(x, Matrix_xSym)

    geMatrix_getDiag_1;
    SEXP ret = PROTECT(allocVector(REALSXP, nret));
    double *rv = REAL(ret),
	   *xv = REAL(x_x);

#define geMatrix_getDiag_2			\
    for (i = 0; i < nret; i++) {		\
	rv[i] = xv[i * m1];			\
    }						\
    UNPROTECT(1);				\
    return ret

    geMatrix_getDiag_2;
}

SEXP lgeMatrix_getDiag(SEXP x)
{
    geMatrix_getDiag_1;

    SEXP ret = PROTECT(allocVector(LGLSXP, nret));
    int *rv = LOGICAL(ret),
	*xv = LOGICAL(x_x);

    geMatrix_getDiag_2;
}

#undef geMatrix_getDiag_1
#undef geMatrix_getDiag_2

SEXP dgeMatrix_setDiag(SEXP x, SEXP d)
{
#define geMatrix_setDiag_1					\
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));		\
    int m = dims[0], nret = (m < dims[1]) ? m : dims[1];	\
    SEXP ret = PROTECT(duplicate(x));				\
    SEXP r_x = GET_SLOT(ret, Matrix_xSym);			\
    int l_d = LENGTH(d); Rboolean d_full = (l_d == nret);	\
    if (!d_full && l_d != 1)					\
	error(_("replacement diagonal has wrong length"))

    geMatrix_setDiag_1;
    double *dv = REAL(d), *rv = REAL(r_x);

#define geMatrix_setDiag_2			\
    R_xlen_t m1 = m + 1;                        \
    if(d_full) for (int i = 0; i < nret; i++)	\
	rv[i * m1] = dv[i];			\
    else for (int i = 0; i < nret; i++)		\
	rv[i * m1] = *dv;			\
    UNPROTECT(1);				\
    return ret

    geMatrix_setDiag_2;
}

SEXP lgeMatrix_setDiag(SEXP x, SEXP d)
{
    geMatrix_setDiag_1;
    int *dv = INTEGER(d), *rv = INTEGER(r_x);

    geMatrix_setDiag_2;
}

#undef geMatrix_setDiag_1
#undef geMatrix_setDiag_2

/* was unused, not replaced: */
SEXP dgeMatrix_addDiag(SEXP x, SEXP d)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	m = dims[0], nret = (m < dims[1]) ? m : dims[1];
    R_xlen_t m1 = m + 1;
    SEXP ret = PROTECT(duplicate(x)),
	r_x = GET_SLOT(ret, Matrix_xSym);
    double *dv = REAL(d), *rv = REAL(r_x);
    int l_d = LENGTH(d); Rboolean d_full = (l_d == nret);
    if (!d_full && l_d != 1)
	error(_("diagonal to be added has wrong length"));

    if(d_full) for (int i = 0; i < nret; i++) rv[i * m1] += dv[i];
    else for (int i = 0; i < nret; i++)	      rv[i * m1] += *dv;
    UNPROTECT(1);
    return ret;
}

#endif /* MJ */

SEXP dgeMatrix_svd(SEXP x, SEXP nnu, SEXP nnv)
{
    int /* nu = asInteger(nnu),
	   nv = asInteger(nnv), */
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));
    SEXP val = PROTECT(allocVector(VECSXP, 3));

    if (dims[0] && dims[1]) {
	int m = dims[0], n = dims[1], mm = (m < n)?m:n,
	    lwork = -1, info;
	double tmp, *work;
	int *iwork, n_iw = 8 * mm;
	if(8 * (double)mm != n_iw) // integer overflow
	    error(_("dgeMatrix_svd(x,*): dim(x)[j] = %d is too large"), mm);
	Matrix_Calloc(iwork, n_iw, int);

	SET_VECTOR_ELT(val, 0, allocVector(REALSXP, mm));
	SET_VECTOR_ELT(val, 1, allocMatrix(REALSXP, m, mm));
	SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, mm, n));
	F77_CALL(dgesdd)("S", &m, &n, xx, &m,
			 REAL(VECTOR_ELT(val, 0)),
			 REAL(VECTOR_ELT(val, 1)), &m,
			 REAL(VECTOR_ELT(val, 2)), &mm,
			 &tmp, &lwork, iwork, &info FCONE);
	lwork = (int) tmp;
	Matrix_Calloc(work, lwork, double);

	F77_CALL(dgesdd)("S", &m, &n, xx, &m,
			 REAL(VECTOR_ELT(val, 0)),
			 REAL(VECTOR_ELT(val, 1)), &m,
			 REAL(VECTOR_ELT(val, 2)), &mm,
			 work, &lwork, iwork, &info FCONE);

	Matrix_Free(iwork, n_iw);
	Matrix_Free(work, lwork);
    }
    UNPROTECT(1);
    return val;
}

const static double padec [] = /* for matrix exponential calculation. */
{
  5.0000000000000000e-1,
  1.1666666666666667e-1,
  1.6666666666666667e-2,
  1.6025641025641026e-3,
  1.0683760683760684e-4,
  4.8562548562548563e-6,
  1.3875013875013875e-7,
  1.9270852604185938e-9,
};

/**
 * Matrix exponential - based on the _corrected_ code for Octave's expm function.
 *
 * @param x real square matrix to exponentiate
 *
 * @return matrix exponential of x
 */
SEXP dgeMatrix_exp(SEXP x)
{
    const double one = 1.0, zero = 0.0;
    const int i1 = 1;
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    const int n = Dims[1];
    const R_xlen_t n_ = n, np1 = n + 1, nsqr = n_ * n_; // nsqr = n^2

    SEXP val = PROTECT(duplicate(x));
    int i, ilo, ilos, ihi, ihis, j, sqpow;
    int *pivot = R_Calloc(n, int);
    double *dpp = R_Calloc(nsqr, double), /* denominator power Pade' */
	*npp = R_Calloc(nsqr, double), /* numerator power Pade' */
	*perm = R_Calloc(n, double),
	*scale = R_Calloc(n, double),
	*v = REAL(GET_SLOT(val, Matrix_xSym)),
	*work = R_Calloc(nsqr, double), inf_norm, m1_j/*= (-1)^j */, trshift;
    R_CheckStack();

    if (n < 1 || Dims[0] != n)
	error(_("Matrix exponential requires square, non-null matrix"));
    if(n == 1) {
	v[0] = exp(v[0]);
	UNPROTECT(1);
	return val;
    }

    /* Preconditioning 1.  Shift diagonal by average diagonal if positive. */
    trshift = 0;		/* determine average diagonal element */
    for (i = 0; i < n; i++) trshift += v[i * np1];
    trshift /= n;
    if (trshift > 0.) {		/* shift diagonal by -trshift */
	for (i = 0; i < n; i++) v[i * np1] -= trshift;
    }

    /* Preconditioning 2. Balancing with dgebal. */
    F77_CALL(dgebal)("P", &n, v, &n, &ilo, &ihi, perm, &j FCONE);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);
    F77_CALL(dgebal)("S", &n, v, &n, &ilos, &ihis, scale, &j FCONE);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);

    /* Preconditioning 3. Scaling according to infinity norm */
    inf_norm = F77_CALL(dlange)("I", &n, &n, v, &n, work FCONE);
    sqpow = (inf_norm > 0) ? (int) (1 + log(inf_norm)/log(2.)) : 0;
    if (sqpow < 0) sqpow = 0;
    if (sqpow > 0) {
	double scale_factor = 1.0;
	for (i = 0; i < sqpow; i++) scale_factor *= 2.;
	for (R_xlen_t i = 0; i < nsqr; i++) v[i] /= scale_factor;
    }

    /* Pade' approximation. Powers v^8, v^7, ..., v^1 */
    Matrix_memset(npp, 0, nsqr, sizeof(double));
    Matrix_memset(dpp, 0, nsqr, sizeof(double));
    m1_j = -1;
    for (j = 7; j >=0; j--) {
	double mult = padec[j];
	/* npp = m * npp + padec[j] *m */
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, npp, &n,
			&zero, work, &n FCONE FCONE);
	for (R_xlen_t i = 0; i < nsqr; i++) npp[i] = work[i] + mult * v[i];
	/* dpp = m * dpp + (m1_j * padec[j]) * m */
	mult *= m1_j;
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, dpp, &n,
			&zero, work, &n FCONE FCONE);
	for (R_xlen_t i = 0; i < nsqr; i++) dpp[i] = work[i] + mult * v[i];
	m1_j *= -1;
    }
    /* Zero power */
    for (R_xlen_t i = 0; i < nsqr; i++) dpp[i] *= -1.;
    for (j = 0; j < n; j++) {
	npp[j * np1] += 1.;
	dpp[j * np1] += 1.;
    }

    /* Pade' approximation is solve(dpp, npp) */
    F77_CALL(dgetrf)(&n, &n, dpp, &n, pivot, &j);
    if (j) error(_("dgeMatrix_exp: dgetrf returned error code %d"), j);
    F77_CALL(dgetrs)("N", &n, &n, dpp, &n, pivot, npp, &n, &j FCONE);
    if (j) error(_("dgeMatrix_exp: dgetrs returned error code %d"), j);
    Memcpy(v, npp, nsqr);

    /* Now undo all of the preconditioning */
    /* Preconditioning 3: square the result for every power of 2 */
    while (sqpow--) {
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, v, &n,
			&zero, work, &n FCONE FCONE);
	Memcpy(v, work, nsqr);
    }
    /* Preconditioning 2: apply inverse scaling */
    for (j = 0; j < n; j++) {
	R_xlen_t jn = j * n_;
	for (i = 0; i < n; i++)
	    v[i + jn] *= scale[i]/scale[j];
    }


    /* 2 b) Inverse permutation  (if not the identity permutation) */
    if (ilo != 1 || ihi != n) {
	/* Martin Maechler's code */

#define SWAP_ROW(I,J) F77_CALL(dswap)(&n, &v[(I)], &n, &v[(J)], &n)

#define SWAP_COL(I,J) F77_CALL(dswap)(&n, &v[(I)*n_], &i1, &v[(J)*n_], &i1)

#define RE_PERMUTE(I)				\
	int p_I = (int) (perm[I]) - 1;		\
	SWAP_COL(I, p_I);			\
	SWAP_ROW(I, p_I)

	/* reversion of "leading permutations" : in reverse order */
	for (i = (ilo - 1) - 1; i >= 0; i--) {
	    RE_PERMUTE(i);
	}

	/* reversion of "trailing permutations" : applied in forward order */
	for (i = (ihi + 1) - 1; i < n; i++) {
	    RE_PERMUTE(i);
	}
    }

    /* Preconditioning 1: Trace normalization */
    if (trshift > 0.) {
	double mult = exp(trshift);
	for (R_xlen_t i = 0; i < nsqr; i++) v[i] *= mult;
    }

    /* Clean up */
    R_Free(work); R_Free(scale); R_Free(perm); R_Free(npp); R_Free(dpp); R_Free(pivot);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_Schur(SEXP x, SEXP vectors, SEXP isDGE)
{
// 'x' is either a traditional matrix or a  dgeMatrix, as indicated by isDGE.
    int *dims, n, vecs = asLogical(vectors), is_dge = asLogical(isDGE),
	info, izero = 0, lwork = -1, nprot = 1;

    if(is_dge) {
	dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    } else { // traditional matrix
	dims = INTEGER(getAttrib(x, R_DimSymbol));
	if(!isReal(x)) { // may not be "numeric" ..
	    x = PROTECT(coerceVector(x, REALSXP)); // -> maybe error
	    nprot++;
	}
    }
    double *work, tmp;
    const char *nms[] = {"WR", "WI", "T", "Z", ""};
    SEXP val = PROTECT(Rf_mkNamed(VECSXP, nms));

    n = dims[0];
    if (n != dims[1] || n < 1)
	error(_("dgeMatrix_Schur: argument x must be a non-null square matrix"));
    const R_xlen_t n2 = ((R_xlen_t)n) * n; // = n^2

    SET_VECTOR_ELT(val, 0, allocVector(REALSXP, n));
    SET_VECTOR_ELT(val, 1, allocVector(REALSXP, n));
    SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, n, n));
    Memcpy(REAL(VECTOR_ELT(val, 2)),
	   REAL(is_dge ? GET_SLOT(x, Matrix_xSym) : x),
	   n2);
    SET_VECTOR_ELT(val, 3, allocMatrix(REALSXP, vecs ? n : 0, vecs ? n : 0));
    F77_CALL(dgees)(vecs ? "V" : "N", "N", NULL, dims, (double *) NULL, dims, &izero,
		    (double *) NULL, (double *) NULL, (double *) NULL, dims,
		    &tmp, &lwork, (int *) NULL, &info FCONE FCONE);
    if (info) error(_("dgeMatrix_Schur: first call to dgees failed"));
    lwork = (int) tmp;
    Matrix_Calloc(work, lwork, double);

    F77_CALL(dgees)(vecs ? "V" : "N", "N", NULL, dims, REAL(VECTOR_ELT(val, 2)), dims,
		    &izero, REAL(VECTOR_ELT(val, 0)), REAL(VECTOR_ELT(val, 1)),
		    REAL(VECTOR_ELT(val, 3)), dims, work, &lwork,
		    (int *) NULL, &info FCONE FCONE);
    Matrix_Free(work, lwork);
    if (info) error(_("dgeMatrix_Schur: dgees returned code %d"), info);
    UNPROTECT(nprot);
    return val;
} // dgeMatrix_Schur

/* MJ: no longer needed ... prefer more general R_dense_(col|row)Sums() */
#if 0

// colSums(), colMeans(), rowSums() and rowMeans() -- called from ../R/colSums.R
SEXP dgeMatrix_colsums(SEXP x, SEXP naRmP, SEXP cols, SEXP mean)
{
    int keepNA = !asLogical(naRmP); // <==>  na.rm = FALSE,  the default
    int doMean = asLogical(mean);
    int useCols = asLogical(cols);
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int i, j, m = dims[0], n = dims[1];
    R_xlen_t m_ = (R_xlen_t) m; // ==> use "long-integer" arithmetic for indices
    SEXP ans = PROTECT(allocVector(REALSXP, (useCols) ? n : m));
    double *aa = REAL(ans), *xx = REAL(GET_SLOT(x, Matrix_xSym));

    if (useCols) {  /* col(Sums|Means) : */
	R_xlen_t cnt = m_; // := number of 'valid' entries in current column
	for (j = 0; j < n; j++) { // column j :
	    double *x_j = xx + m_ * j, s = 0.;

	    if (keepNA)
		for (i = 0; i < m; i++) s += x_j[i];
	    else {
		cnt = 0;
		for (i = 0; i < m; i++)
		    if (!ISNAN(x_j[i])) {cnt++; s += x_j[i];}
	    }
	    if (doMean) {
		if (cnt > 0) s /= cnt; else s = NA_REAL;
	    }
	    aa[j] = s;
	}
    } else { /* row(Sums|Means) : */
	Rboolean do_count = (!keepNA) && doMean;
	int *cnt = (int*) NULL;
	if (do_count) Matrix_Calloc(cnt, m, int);

	// (taking care to access x contiguously: vary i inside j)
	for (i = 0; i < m; i++) {
	    aa[i] = 0.;
	    if(do_count) cnt[i] = 0;
	}
	for (j = 0; j < n; j++) {
	    if (keepNA)
		for (i = 0; i < m; i++) aa[i] += xx[i + j * m_];
	    else
		for (i = 0; i < m; i++) {
		    double el = xx[i + j * m_];
		    if (!ISNAN(el)) {
			aa[i] += el;
			if (doMean) cnt[i]++;
		    }
		}
	}
	if (doMean) {
	    if (keepNA)
		for (i = 0; i < m; i++) aa[i] /= n;
	    else
		for (i = 0; i < m; i++)
		    aa[i] = (cnt[i] > 0) ? aa[i]/cnt[i] : NA_REAL;
	}
	if (do_count) Matrix_Free(cnt, m);
    }

    SEXP nms = VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), useCols ? 1 : 0);
    if(!isNull(nms))
	setAttrib(ans, R_NamesSymbol, duplicate(nms));
    UNPROTECT(1);
    return ans;
}

#endif /* MJ */
