#include "cscBlocked.h"
/* TODO
 *  - code for trans = 'T' in cscb_syrk
 *  - code for non-trivial cscb_trmm and cscb_ldl
 */
SEXP cscBlocked_validate(SEXP x)
{
    SEXP pp = GET_SLOT(x, Matrix_pSym),
	ip = GET_SLOT(x, Matrix_iSym),
	xp = GET_SLOT(x, Matrix_xSym),
	dp = getAttrib(xp, R_DimSymbol);
    int *pv = INTEGER(pp),
	*dim = INTEGER(dp),
	ncol = length(pp) - 1;
    int nnz = pv[ncol];

    if (!(isReal(xp) && isArray(xp)))
	return mkString("slot x should be a real array");
    if (length(dp) != 3)
	return mkString("slot x should be a 3-dimensional array");
    if (length(ip) != nnz)
	return mkString("length of slot i does not matck last element of slot p");
    if (dim[2] != nnz)
	return
	    mkString("third dimension of slot x does not match number of nonzeros");
    return ScalarLogical(1);
}

static int*
expand_column_pointers(int ncol, const int mp[], int mj[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int j2 = mp[j+1], jj;
	for (jj = mp[j]; jj < j2; jj++) mj[jj] = j;
    }
    return mj;
}

static R_INLINE
int Tind(const int rowind[], const int colptr[], int i, int j)
{
    int k, k2 = colptr[j + 1];
    for (k = colptr[j]; k < k2; k++)
	if (rowind[k] == i) return k;
    error("row %d and column %d not defined in rowind and colptr",
	  i, j);
    return -1;			/* -Wall */
}

/** 
 * Perform one of the matrix operations 
 *  C := alpha*op(A)*B + beta*C
 * or
 *  C := alpha*B*op(A) + beta*C
 * where A is a compressed, sparse, blocked matrix and
 * B and C are dense matrices.
 * 
 * @param side 'L' or 'l' for left, otherwise right
 * @param transa 'T' or 't' for transpose.
 * @param m number of rows in C
 * @param n number of columns in C
 * @param k number of rows in B if side = 'L', otherwise
 *        number of columns in B.
 * @param alpha
 * @param nr number of rows per block of A
 * @param nc number of columns per block of A
 * @param ap vector of column pointers in A
 * @param ai vector of row indices in A
 * @param ax contents of non-zero blocks of A
 * @param b matrix to be multiplied
 * @param ldb leading dimension of b as declared in the calling
 *        routine
 * @param beta scalar multiplier of c
 * @param c product matrix to be modified
 * @param ldc leading dimension of c as declared in the calling
 *        routine
 */
void
cscBlocked_mm(char side, char transa, int m, int n, int k,
	      double alpha, int nr, int nc,
	      const int ap[], const int ai[],
	      const double ax[],
	      const double b[], int ldb,
	      double beta, double c[], int ldc)
{
    int j, kk, lside = (side == 'L' || side == 'l');
    int ncb, nrb, sz = nr * nc, tra = (transa == 'T' || transa == 't');

    if (nr < 1 || nc < 1 || m < 0 || n < 0 || k < 0)
	error("improper dims m=%d, n=%d, k=%d, nr=%d, nc=%d",
		  m, n, k, nr, nc);
    if (ldc < n) error("incompatible dims n=%d, ldc=%d", n, ldc);
    if (lside) {
	if (ldb < k)
	    error("incompatible L dims k=%d, ldb=%d, n=%d, nr=%d, nc=%d",
		  k, ldb, n, nr, nc);
	if (tra) {
	    if (m % nc || k % nr)
		error("incompatible LT dims m=%d, k = %d, nr=%d, nc=%d",
		      m, k, nr, nc);
	    nrb = k/nr; ncb = m/nc;
	} else {
	    if (m % nr || k % nc)
		error("incompatible LN dims m=%d, k = %d, nr=%d, nc=%d",
		      m, k, nr, nc);
	    nrb = m/nr; ncb = k/nc;
	}
	for (j = 0; j < ncb; j++) {
	    int j2 = ap[j + 1];
	    for (kk = ap[j]; kk < j2; kk++) {
		int ii = ai[kk];
		if (ii < 0 || ii >= nrb)
		    error("improper row index ii=%d, nrb=%d", ii, nrb);
		if (tra) {
		    F77_CALL(dgemm)("T", "N", &nc, &n, &nr,
				    &alpha, ax + kk * sz, &nr,
				    b + ii * nr, &ldb,
				    &beta, c + j * nc, &ldc);
		} else {
		    F77_CALL(dgemm)("N", "N", &nr, &n, &nc,
				    &alpha, ax + kk * sz, &nr,
				    b + j*nc, &ldb,
				    &beta, c + ii * nr, &ldc);
		}
	    }
	}
    } else {
	error("Call to cscBlocked_mm must have side == 'L'");
    }
}

/** 
 * Perform one of the matrix operations 
 *  C := alpha*op(A)*B + beta*C
 * or
 *  C := alpha*B*op(A) + beta*C
 * where A is a compressed, sparse, blocked matrix and
 * B and C are dense matrices.
 * 
 * @param side 'L' or 'l' for left, otherwise right
 * @param transa 'T' or 't' for transpose.
 * @param m number of rows in C
 * @param n number of columns in C
 * @param k number of rows in B if side = 'L', otherwise
 *        number of columns in B.
 * @param alpha
 * @param A pointer to a cscBlocked object
 * @param B matrix to be multiplied
 * @param ldb leading dimension of b as declared in the calling
 *        routine
 * @param beta scalar multiplier of c
 * @param C product matrix to be modified
 * @param ldc leading dimension of c as declared in the calling
 *        routine
 */
void
cscb_mm(char side, char transa, int m, int n, int k,
	double alpha, SEXP A,
	const double B[], int ldb,
	double beta, double C[], int ldc)
{
    int lside = (side == 'L' || side == 'l'),
	tra = (transa == 'T' || transa == 't');
    SEXP AxP = GET_SLOT(A, Matrix_xSym),
	ApP = GET_SLOT(A, Matrix_pSym);
    int *adims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	*Ap = INTEGER(ApP),
	*Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	ancb = length(ApP) - 1, /* number of column blocks */
	anrb;			/* number of row blocks */
    int absz = adims[0] * adims[1]; /* block size */
    int j;
    double *Ax = REAL(AxP);

    if (ldc < m) error("incompatible dims m=%d, ldc=%d", m, ldc);
    if (lside) {
	/* B is of size k by n */
	if (ldb < k)
	    error("incompatible L dims k=%d, ldb=%d, n=%d, nr=%d, nc=%d",
		  k, ldb, n, adims[0], adims[1]);
	if (tra) {
	    if (m % adims[1] || k % adims[0])
		error("incompatible LT dims m=%d, k = %d, nr=%d, nc=%d",
		      m, k, adims[0], adims[1]);
	    if (ancb != m/adims[1])
		error("incompatible LT dims m=%d, ancb=%d, adims=[%d,%d,%d]",
		      m, ancb, adims[0], adims[1], adims[3]);
	    anrb = k/adims[0];
	} else {
	    if (m % adims[0] || k % adims[1])
		error("incompatible LN dims m=%d, k = %d, nr=%d, nc=%d",
		      m, k, adims[0], adims[1]);
	    if (ancb != k/adims[1])
		error("incompatible LN dims k=%d, ancb=%d, adims=[%d,%d,%d]",
		      k, ancb, adims[0], adims[1], adims[3]);
	    anrb = m/adims[0];
	}
	for (j = 0; j < ancb; j++) {
	    int kk, j2 = Ap[j + 1];
	    for (kk = Ap[j]; kk < j2; kk++) {
		int ii = Ai[kk];
		if (ii < 0 || ii >= anrb)
		    error("improper row index ii=%d, anrb=%d", ii, anrb);
		if (tra) {
		    F77_CALL(dgemm)("T", "N", adims+1, &n, adims,
				    &alpha, Ax + kk * absz, adims,
				    B + ii * adims[0], &ldb,
				    &beta, C + j * adims[1], &ldc);
		} else {
		    F77_CALL(dgemm)("N", "N", adims, &n, adims+1,
				    &alpha, Ax + kk * absz, adims,
				    B + j*adims[1], &ldb,
				    &beta, C + ii * adims[0], &ldc);
		}
	    }
	}
    } else {
	/* B is of size m by k */
	error("Call to cscb_mm must have side == 'L'");
    }
}

/** 
 * Invert a triangular sparse blocked matrix.  This is not done in
 * place because the number of non-zero blocks in A-inverse can be
 * different than the number of non-zero blocks in A.
 * 
 * @param upper 'U' indicates upper triangular, 'L' lower
 * @param unit 'U' indicates unit diagonal, 'N' non-unit
 * @param A Pointer to a triangular cscBlocked object
 * @param Ai Pointer to a triangular cscBlocked object
 */
void
cscb_tri(char upper, char unit, SEXP A, const int Parent[], SEXP AI)
{
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym),
	AIpP = GET_SLOT(AI, Matrix_pSym),
	AIxP = GET_SLOT(AI, Matrix_xSym);
    int *Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*AIi = INTEGER(GET_SLOT(AI, Matrix_iSym)),
	*AIp = INTEGER(AIpP),
	*adim = INTEGER(getAttrib(AxP, R_DimSymbol)),
	*aidim = INTEGER(getAttrib(AxP, R_DimSymbol)),
	ancb = length(ApP) - 1,
	aincb = length(AIpP) - 1,
	iup = (upper == 'U' || upper == 'u'),
	iunit = (unit == 'U' || unit == 'u');

    if (aidim[0] != adim[0] || aidim[1] != adim[1] || aidim[2] != adim[2] ||
	adim[0] != adim[1] || ancb != aincb)
	error("Incompatible dimensions A[%d,%d,%d]x%d; AI[%d,%d,%d]x%d",
	      adim[0], adim[1], adim[2], ancb, 
	      aidim[0], aidim[1], aidim[2], aincb);
    if (!iunit)
	error("Code for non-unit triangular matrices not yet written");
    if (adim[2] > 0)
	error("Code for non-trivial unit inverse not yet written");
}

/** 
 * Search for the element in a compressed sparse matrix at a given row and column
 * 
 * @param row row index
 * @param col column index
 * @param cp column pointers
 * @param ci row indices
 * 
 * @return index of element in ci, if it exists, else -1
 */
static R_INLINE
int Ind(int row, int col, const int cp[], const int ci[])
{
    int i, i2 = cp[col + 1];
    for (i = cp[col]; i < i2; i++)
	if (ci[i] == row) return i;
    return -1;
}

/** 
 * Perform one of the matrix operations 
 *  C := alpha*A*A' + beta*C,
 * or
 *  C := alpha*A'*A + beta*C,
 * where A is a compressed, sparse, blocked matrix and
 * C is a compressed, sparse, symmetric blocked matrix.
 * 
 * @param uplo 'U' or 'u' for upper triangular storage, else lower.
 * @param trans 'T' or 't' for transpose.
 * @param alpha scalar multiplier of outer product
 * @param A compressed sparse blocked matrix
 * @param beta scalar multiplier of c
 * @param C compressed sparse blocked symmetric matrix to be updated
 */
void cscb_syrk(char uplo, char trans, double alpha, SEXP A, double beta, SEXP C)
{
    SEXP AxP = GET_SLOT(A, Matrix_xSym),
	ApP = GET_SLOT(A, Matrix_pSym),
	CxP = GET_SLOT(C, Matrix_xSym),
	CpP = GET_SLOT(C, Matrix_pSym);
    int *adims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	*Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*cdims = INTEGER(getAttrib(CxP, R_DimSymbol)),
	*Ci = INTEGER(GET_SLOT(C, Matrix_iSym)),
	*Cp = INTEGER(CpP),
	iup = (uplo == 'U' || uplo == 'u'),
	itr = (trans == 'T' || trans == 't'),
	j, k;
    double *Ax = REAL(AxP), *Cx = REAL(CxP), one = 1.;
    int scalar = (adims[0] == 1 && adims[1] == 1),
	anc = length(ApP) - 1,
	asz = adims[0] * adims[1],
	csz = cdims[0] * cdims[1];


    if (cdims[0] != cdims[1]) error("blocks in C must be square");
    if (itr) {
	error("Code for trans == 'T' not yet written");
				/* FIXME: Write this part */
    } else {
	if (adims[0] != cdims[0])
	    error("Inconsistent dimensions: A[%d,%d,%d], C[%d,%d,%d]",
		  adims[0], adims[1], adims[2],
		  cdims[0], cdims[1], cdims[2]);
				/* check the row indices */
	for (k = 0; k < adims[2]; k++) {
	    int aik = Ai[k];
	    if (aik < 0 || aik >= cdims[2])
		error("Row index %d = %d is out of range [0, %d]",
		      k, aik, cdims[2] - 1);
	}
				/* multiply C by beta */
	if (beta != 1.)
	    for (j = 0; j < csz * cdims[2]; j++) Cx[j] *= beta;
				/* individual products */
	for (j = 0; j < anc; j++) {
	    int k, kk, k2 = Ap[j+1];
	    for (k = Ap[j]; k < k2; k++) {
		int ii = Ai[k], K = Ind(ii, ii, Cp, Ci);

		if (K < 0) error("cscb_syrk: C[%d,%d] not defined", ii, ii);
		if (scalar) Cx[K] += alpha * Ax[k] * Ax[k];
		else F77_CALL(dsyrk)(&uplo, "N", cdims, adims + 1,
				     &alpha, Ax + k * asz, adims,
				     &one, Cx + K * csz, cdims);

		for (kk = k+1; kk < k2; kk++) {
		    int jj = Ai[kk];
		    K = (iup) ? Ind(ii, jj, Cp, Ci) : Ind(jj, ii, Cp, Ci);
		    
		    if (K < 0) error("cscb_syrk: C[%d,%d] not defined", ii, jj);
		    if (scalar) Cx[K] += alpha * Ax[k] * Ax[kk];
		    else F77_CALL(dgemm)("N", "T", cdims, cdims + 1, adims + 1,
					 &alpha, Ax + ((iup)?kk:k) * asz, adims,
					 Ax + ((iup)?k:kk) * asz, adims,
					 &one, Cx + K * asz, cdims);
		}
	    }
	}
    }
}

/** 
 * Create the LDL' decomposition of the positive definite symmetric
 * cscBlocked matrix A (upper triangle stored) in L and D.
 * 
 * @param A pointer to a cscBlocked object containing the upper
 * triangle of a positive definite symmetric matrix.
 * @param Parent the parent array for A
 * @param L pointer to a cscBlocked object to hold L
 * @param D pointer to a 3D array to hold D
 * 
 * @return n the number of column blocks in A for success.  A value
 * less than n indicates the first column block whose diagonal was not
 * positive definite.
 */
int
cscb_ldl(SEXP A, const int Parent[], SEXP L, SEXP D)
{
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym);
    int *adims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	diag, info, j, k, n = length(ApP) - 1;
    int *Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*Li = INTEGER(GET_SLOT(L, Matrix_iSym)),
	*Lp = INTEGER(GET_SLOT(L, Matrix_pSym)), nci = adims[0];
    double *Lx = REAL(GET_SLOT(L, Matrix_xSym)),
	*Ax = REAL(AxP), *Dx = REAL(D), minus1 = -1., one = 1.;
    
    if (adims[1] != nci || nci < 1)
	error("cscb_ldl: dim(A) is [%d, %d, %d]", adims[0], adims[1], adims[2]);
    for (j = 0, diag = 1; j < n; j++) { /* check for trivial structure */
	if (Parent[j] >= 0) {diag = 0; break;}
    }
    if (diag) {
	int ncisqr = nci * nci;
	Memcpy(Dx, Ax, ncisqr * n);
	for (j = 0; j < n; j++) { /* form D_i^{1/2} */
	    F77_CALL(dpotrf)("U", &nci, Dx + j * ncisqr, &nci, &k);
	    if (k) error("D[ , , %d], level %d, is not positive definite",
			 j + 1);
	}
	return n;
    }
    if (nci == 1) {
	k = R_ldl_numeric(n, Ap, Ai, Ax, Lp, Parent, Li, Lx, Dx,
			  (int *) NULL, (int *) NULL);
	if (k < n) error("cscb_ldl: error code %k from R_ldl_numeric", k);
	for (j = 0; j < n; j++) Dx[j] = sqrt(Dx[j]);
	return n;
    } else {		   /* Copy of ldl_numeric from the LDL package
			    * modified for blocked sparse matrices */ 
	int i, k, p, p2, len, nci = adims[0], ncisqr = adims[0]*adims[0], top;
	int *Lnz = Calloc(n, int),
	    *Pattern = Calloc(n, int),
	    *Flag = Calloc(n, int);
	double *Y = Calloc(n * ncisqr, double), *Yi = Calloc(ncisqr, double);

	for (k = 0; k < n; k++) {
	    /* compute nonzero Pattern of kth row of L, in topological order */
	    AZERO(Y + k*ncisqr, ncisqr); /* Y[,,0:k] is now all zero */
	    top = n;		/* stack for pattern is empty */
	    Flag[k] = k;	/* mark node k as visited */
	    Lnz[k] = 0;		/* count of nonzeros in column k of L */
	    p2 = Ap[k+1];
	    for (p = Ap[k]; p < p2; p++) {
		i = Ai[p];	/* get A[i,k] */
		if (i <= k) {	/* [i,k] in upper triangle? Should always be true */
				/* copy A(i,k) into Y */ 
		    Memcpy(Y + i * ncisqr, Ax + p * ncisqr, ncisqr); 
		    /* follow path from i to root of etree, stop at flagged node */
		    for (len = 0; Flag[i] != k; i = Parent[i]) {
			Pattern[len++] = i; /* L[k,i] is nonzero */
			Flag[i] = k; /* mark i as visited */
		    }
		    while (len > 0) { /* push path on top of stack */
			Pattern[--top] = Pattern[--len];
		    }
		}
	    }
	    /* Pattern [top ... n-1] now contains nonzero pattern of L[,k] */
	    /* compute numerical values in kth row of L (a sparse triangular solve) */
	    Memcpy(Dx + k * ncisqr, Y + k * ncisqr, ncisqr); /* get D[,,k] */
	    AZERO(Y + k*ncisqr, ncisqr); /* clear Y[,,k] */
	    for (; top < n; top++) {
		i = Pattern[top];
		Memcpy(Yi, Y + i*ncisqr, ncisqr); /* copy Y[,,i] */
		AZERO(Y + i*ncisqr, ncisqr); /* clear Y[,,i] */
		p2 = Lp[i] + Lnz[i];
		for (p = Lp[i]; p < p2; p++) {
		    F77_CALL(dgemm)("N", "N", &nci, &nci, &nci, &minus1,
				    Lx + p*ncisqr, &nci, Yi, &nci,
				    &one, Y + Li[p]*ncisqr, &nci);
		}
		/* FIXME: Check that this is the correct order and transposition */
		F77_CALL(dtrsm)("L", "U", "N", "N", &nci, &nci,
				&one, Dx + i*ncisqr, &nci, Yi, &nci);
		F77_CALL(dsyrk)("U", "N", &nci, &nci, &minus1, Yi, &nci,
				&one, Dx + k*ncisqr, &nci);
		F77_CALL(dtrsm)("L", "U", "T", "N", &nci, &nci,
				&one, Dx + i*ncisqr, &nci, Yi, &nci);
		Li[p] = k;	/* store L[k,i] in column form of L */
		Memcpy(Lx + p * ncisqr, Yi, ncisqr);
		Lnz[i]++;	/* increment count of nonzeros in col i */
	    }
	    F77_CALL(dpotrf)("U", &nci, Dx + k*ncisqr, &nci, &info);
	    if (info) {
		Free(Y); Free(Yi); Free(Pattern); Free(Flag); Free(Lnz); 
		return k;	    /* failure, D[,,k] is not positive definite */
	    }
	}
	Free(Y); Free(Yi); Free(Pattern); Free(Flag); Free(Lnz);
	return n;	/* success, diagonal of D is all nonzero */
    }
    return -1;			/* keep -Wall happy */
}

/** 
 * Perform one of the cscBlocked-matrix operations B := alpha*op(A)*B
 * or B := alpha*B*op(A)
 * 
 * @param side 'L' or 'R' for left or right
 * @param uplo 'U' or 'L' for upper or lower
 * @param transa 'T' or 'N' for transpose or no transpose
 * @param diag 'U' or 'N' for unit diagonal or non-unit
 * @param A pointer to a triangular cscBlocked object
 * @param B pointer to the contents of the matrix B
 * @param m number of rows in B
 * @param n number of columns in B
 * @param ldb leading dimension of B as declared in the calling function
 */
void cscb_trmm(char side, char uplo, char transa, char diag,
	       double alpha, SEXP A, double B[], int m, int n, int ldb)
{
    int ileft = (side == 'L' || side == 'l'),
	iup = (uplo == 'U' || uplo == 'u'),
	itr = (transa == 'T' || transa == 't'),
	iunit = (diag == 'U' || diag == 'u');
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym);
    int *Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*xdims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	i, j, nb = length(ApP) - 1;
    
    if (xdims[0] != xdims[1])
	error("Argument A to cscb_trmm is not triangular");
    if (alpha != 1.0) {
	for (j = 0; j < n; j++) { /* scale by alpha */
	    for (i = 0; i < m; i++)
		B[i + j * ldb] *= alpha;
	}
    }
    if (iunit && xdims[2] < 1) return; /* A is the identity */
    error("Code for non-identity cases of cscb_trmm not yet written");
    iup += 0;			/* keep -Wall happy */
    ileft += 0;
    itr += 0;
}

/** 
 * Solve a triangular system of the form op(A)*X = alpha*B where A
 * is a cscBlocked triangular matrix and B is a dense matrix.
 * 
 * @param uplo 'U' or 'L' for upper or lower
 * @param trans 'T' or 'N' for transpose or no transpose
 * @param diag 'U' or 'N' for unit diagonal or non-unit
 * @param A pointer to a triangular cscBlocked object
 * @param B pointer to the contents of the matrix B
 * @param m number of rows in B
 * @param n number of columns in B
 * @param ldb leading dimension of B as declared in the calling function
 */
void cscb_trsm(char uplo, char transa, char diag,
	       double alpha, SEXP A, double B[], int m, int n, int ldb)
{
    int iup = (uplo == 'U' || uplo == 'u'),
	itr = (transa == 'T' || transa == 't'),
	iunit = (diag == 'U' || diag == 'u');
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym);
    int *Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*xdims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	i, j, nb = length(ApP) - 1;
    double *Ax = REAL(GET_SLOT(A, Matrix_xSym));
    
    if (xdims[0] != xdims[1])
	error("Argument A to cscb_trsm is not triangular");
    if (ldb < m || ldb <= 0 || n <= 0)
	error("cscb_trsm: inconsistent dims m = %d, n = %d, ldb = %d",
	      m, n, ldb);
    if (m != (nb * xdims[0]))
	error("cscb_trsm: inconsistent dims m = %d, A[%d,%d,]x%d",
	      m, xdims[0], xdims[1], xdims[2]);
    if (alpha != 1.0) {
	for (j = 0; j < n; j++) { /* scale by alpha */
	    for (i = 0; i < m; i++)
		B[i + j * ldb] *= alpha;
	}
    }
    if (iunit) {
	if (xdims[2] < 1) return; /* A is the identity */
	if (xdims[0] == 1) {
	    if (iup) error("Code for upper triangle not yet written");
	    if (itr) {
		for (j = 0; j < n; j++)
		    R_ldl_ltsolve(m, B + j * ldb, Ap, Ai, Ax);
	    } else {
		for (j = 0; j < n; j++)
		    R_ldl_lsolve(m, B + j * ldb, Ap, Ai, Ax);
	    }
	    return;
	}
	error("Code for non-scalar cscBlocked objects not yet written");
    }
    error("Code for non-unit cases of cscb_trsm not yet written");
}

/** 
 * Perform one of the operations B := alpha*op(A)*B or
 * B := alpha*B*op(A) where A and B are both cscBlocked.
 * 
 * @param side 'L' or 'R' for left or right
 * @param uplo 'U' or 'L' for upper or lower
 * @param transa 'T' or 'N' for transpose or no transpose
 * @param diag 'U' or 'N' for unit diagonal or non-unit
 * @param alpha scalar multiplier
 * @param A pointer to a triangular cscBlocked object
 * @param B pointer to a general cscBlocked matrix
 */
void cscb_trcbm(char side, char uplo, char transa, char diag,
		double alpha, SEXP A, SEXP B)
{
    int ileft = (side == 'L' || side == 'l'),
	iup = (uplo == 'U' || uplo == 'u'),
	itr = (transa == 'T' || transa == 't'),
	iunit = (diag == 'U' || diag == 'u');
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym),
	BpP = GET_SLOT(B, Matrix_pSym),
	BxP = GET_SLOT(B, Matrix_xSym);
    int *Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*Bi = INTEGER(GET_SLOT(B, Matrix_iSym)),
	*Bp = INTEGER(BpP),
	*axdims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	*bxdims = INTEGER(getAttrib(BxP, R_DimSymbol)),
	ncbA = length(ApP) - 1,
	ncbB = length(BpP) - 1;
    int i, nbx = bxdims[0] * bxdims[1] * bxdims[2];

    if (axdims[0] != axdims[1])
	error("Argument A to cscb_trcbm is not triangular");
    if (alpha != 1.0) {
	for (i = 0; i < nbx; i++) { /* scale by alpha */
	    REAL(BxP)[i] *= alpha;
	}
    }
    if (iunit && axdims[2] < 1) return; /* A is the identity */
    error("Code for non-trivial cscb_trcbm not yet written");
}

/** 
 * Expand a column of a compressed, sparse, column-oriented matrix.
 * 
 * @param dest array to hold the result
 * @param m number of rows in the matrix
 * @param j index (0-based) of column to expand
 * @param Ap array of column pointers
 * @param Ai array of row indices
 * @param Ax array of non-zero values
 * 
 * @return dest
 */
static
double *expand_column(double *dest, int m, int j,
		      const int Ap[], const int Ai[], const double Ax[])
{
    int k, k2 = Ap[j + 1];

    for (k = 0; k < m; k++) dest[k] = 0.;
    for (k = Ap[j]; k < k2; k++) dest[Ai[k]] = Ax[k];
    return dest;
}

/** 
 * Solve one of the systems op(A)*X = alpha*B or
 * X*op(A) = alpha*B where A cscBlocked triangular and B is cscBlocked.
 * 
 * @param side 'L' or 'R' for left or right
 * @param uplo 'U' or 'L' for upper or lower
 * @param transa 'T' or 'N' for transpose or no transpose
 * @param diag 'U' or 'N' for unit diagonal or non-unit
 * @param alpha scalar multiplier
 * @param A pointer to a triangular cscBlocked object
 * @param B pointer to a general cscBlocked matrix
 */
void cscb_trcbsm(char side, char uplo, char transa, char diag,
		 double alpha, SEXP A, const int Parent[], SEXP B)
{
    int ileft = (side == 'L' || side == 'l'),
	iup = (uplo == 'U' || uplo == 'u'),
	itr = (transa == 'T' || transa == 't'),
	iunit = (diag == 'U' || diag == 'u');
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym),
	BpP = GET_SLOT(B, Matrix_pSym),
	BxP = GET_SLOT(B, Matrix_xSym);
    int *Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Ap = INTEGER(ApP),
	*Bi = INTEGER(GET_SLOT(B, Matrix_iSym)),
	*Bp = INTEGER(BpP),
	*axdims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	*bxdims = INTEGER(getAttrib(BxP, R_DimSymbol)),
	ncbA = length(ApP) - 1,
	ncbB = length(BpP) - 1;
    int i, j, nbx = bxdims[0] * bxdims[1] * bxdims[2];
    double *Ax = REAL(AxP), *Bx = REAL(BxP);

    if (axdims[0] != axdims[1])
	error("Argument A to cscb_trcbm is not triangular");
    if (alpha != 1.0) {
	for (i = 0; i < nbx; i++) { /* scale by alpha */
	    REAL(BxP)[i] *= alpha;
	}
    }
    if (iunit && axdims[2] < 1) return;	/* A is the identity */
    if (iunit && axdims[0] == 1) { /* can use R_ldl code */
	if (!ileft && itr) {	/* case required for lmer */
	    int *BTp, nnz = bxdims[2], nrbB;
	    int *tmp = expand_column_pointers(ncbB, Bp, Calloc(nnz, int));
	    int *BTi = Calloc(nnz, int);
	    double *BTx = Calloc(nnz, double), *rhs;

				/* transpose B */
	    for (i = 0, nrbB = -1; i < nnz; i++) if (Bi[i] > nrbB) nrbB = Bi[i];
	    BTp = Calloc(nrbB, int);
	    triplet_to_col(ncbB, nrbB, nnz, tmp, Bi, Bx, BTp, BTi, BTx);
				/* sanity check */
	    if (BTp[nrbB] != nnz) error("cscb_trcbsm: transpose operation failed");
	    Free(tmp);
				/* Solve one column at a time */
	    rhs = Calloc(ncbB, double);
	    AZERO(Bx, nnz);	/* zero the result */
	    for (i = 0; i < nrbB; i++) {
		R_ldl_lsolve(ncbB, expand_column(rhs, ncbB, i, BTp, BTi, BTx),
			     Ap, Ai, Ax);
		for (j = 0; j < ncbB; j++) { /* write non-zeros in sol'n into B */
		    if (BTx[j]) Bx[Tind(Bi, Bp, j, i)] = BTx[j];
		}
		Free(rhs); Free(BTp); Free(BTx); Free(BTi);
	    }
	}
	error("cscb_trcbsm: method not yet written");
    }
    error("cscb_trcbsm: method not yet written");
}

/** 
 * Perform one of the matrix-matrix operations 
 *      C := alpha*op(A)*op(B) + beta*C
 * on compressed, sparse, blocked matrices.
 * 
 * @param transa 'T' for transpose of A, else 'N'
 * @param transb 'T' for transpose of B, else 'N'
 * @param alpha scalar multiplier
 * @param A pointer to a cscBlocked object
 * @param B pointer to a cscBlocked object
 * @param beta scalar multiplier
 * @param C pointer to a cscBlocked object
 */
void cscb_cscbm(char transa, char transb, double alpha, SEXP A, SEXP B,
	     double beta, SEXP C)
{
    int ta = (transa == 'T' || transa == 't') ? 1 : 0,
	tb = (transb == 'T' || transb == 't') ? 1 : 0;
    SEXP ApP = GET_SLOT(A, Matrix_pSym),
	AxP = GET_SLOT(A, Matrix_xSym),
	BpP = GET_SLOT(B, Matrix_pSym),
	BxP = GET_SLOT(B, Matrix_xSym),
	CxP = GET_SLOT(C, Matrix_xSym);
    int *Ap = INTEGER(ApP),
	*Ai = INTEGER(GET_SLOT(A, Matrix_iSym)),
	*Bp = INTEGER(BpP),
	*Bi = INTEGER(GET_SLOT(B, Matrix_iSym)),
	*Cp = INTEGER(GET_SLOT(C, Matrix_pSym)),
	*Ci = INTEGER(GET_SLOT(C, Matrix_iSym)),
	*adims = INTEGER(getAttrib(AxP, R_DimSymbol)),
	*bdims = INTEGER(getAttrib(BxP, R_DimSymbol)),
	*cdims = INTEGER(getAttrib(CxP, R_DimSymbol)),
	nca = length(ApP) - 1,
	ncb = length(BpP) - 1;
    int ablk = adims[0] * adims[1],
	bblk = bdims[0] * bdims[1],
	cblk = cdims[0] * cdims[1];
    double *Ax = REAL(AxP),
	*Bx = REAL(BxP),
	*Cx = REAL(CxP),
	one = 1.0;

    if ((!ta) && tb) {		/* transposed crossproduct */
	int jj;

	if (adims[1] != bdims[1] ||
	    adims[0] != cdims[0] ||
	    bdims[0] != cdims[1])
	    error("C[%d,%d,%d] := A[%d,%d,%d] %*% t(B[%d,%d,%d])",
		  cdims[0], cdims[1], cdims[2],
		  adims[0], adims[1], adims[2],
		  bdims[0], bdims[1], bdims[2]);
	if (nca != ncb)
	    error("C := A(ncblocks = %d) %*% t(B(ncblocks = %d)", nca, ncb);
	if (beta != 1.) {	/* scale C by beta */
	    int ctot = cdims[0] * cdims[1] * cdims[2];
	    for (jj = 0; jj < ctot; jj++) Cx[jj] *= beta;
	}
	for (jj = 0; jj < nca; jj++) {
	    int ia, ib, a2 = Ap[jj + 1], b2 = Bp[jj + 1];
	    for (ia = Ap[jj]; ia < a2; ia++) {
		for (ib = Bp[jj]; ib < b2; ib++) {
		    int cind = Ind(Ai[ia], Bi[ib], Cp, Ci);
		    if (cind < 0)
			error("Invalid index [%d,%d]", Ai[ia], Bi[ib]);
		    F77_CALL(dgemm)("N", "T", cdims, cdims + 1, adims + 1,
				    &alpha, Ax + ia * ablk, adims,
				    Bx + ib * bblk, bdims, &one,
				    Cx + cind * cblk, cdims);
		}
	    }
	}
	return;
    }
    error("Code not yet written");
}
