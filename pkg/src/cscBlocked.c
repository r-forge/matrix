#include "cscBlocked.h"

SEXP cscBlocked_validate(SEXP x)
{
    SEXP pp = GET_SLOT(x, Matrix_pSym),
	ip = GET_SLOT(x, Matrix_iSym),
	xp = GET_SLOT(x, Matrix_xSym),
	dp = getAttrib(xp, R_DimSymbol);
    int *pv = INTEGER(pp),
	*iv = INTEGER(ip),
	*dim = INTEGER(dp),
	ncol = length(pp) - 1;
    int nnz = pv[ncol];

    if (!isReal(xp))
	return ScalarString(mkChar("slot x should be a real array"));
    if (length(dp) != 3)
	return ScalarString(mkChar("slot x should be a 3-dimensional array"));
    if (length(ip) != nnz)
	return ScalarString(mkChar("length of slot i does not matck last element of slot p"));
    if (dim[2] != nnz)
	return ScalarString(mkChar("third dimension of slot x does not match number of nonzeros"));
    return ScalarLogical(1);
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
 * @param k number of columns in B if side = 'L', otherwise
 *        number of rows in B.
 * @param alpha
 * @param nr number of rows per block of A
 * @param nc number of columns per block of A
 * @param ap vector of column pointers in A
 * @param ai vector of row indices in A
 * @param ax contents of non-zero blocks of A
 * @param b matrix to be multiplied
 * @param ldb leading dimension of b as declared in the calling
 *        routine
 * @param c product matrix to be modified
 * @param ldc leading dimension of b as declared in the calling
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
    int lside = (side == 'L' || side == 'l') ? 1 : 0;
    int tra = (transa == 'T' || transa == 't') ? 1 : 0;

    if (lside && !tra) {	/* only case written so far */
	int j, ncb = k/nc, nrb = m/nr, sz = nc * nr;

	if (nr < 1 || nc < 1 || m < 0 || n < 0 || k < 0 ||
	    m % nr || k % nc)
	    error("incompatible dims m=%d, n=%d, k=%d, nr=%d, nc=%d",
		  m, n, k, nr, nc);
	if (ldb < k) error("incompatible dims k=%d, ldb=%d", k, ldb);
	if (ldc < n) error("incompatible dims n=%d, ldc=%d", n, ldc);
	for (j = 0; j < ncb; j++) {
	    int i, i2 = ap[j+1];
	    for (i = ap[j]; i < i2; i++) {
		int rb = ai[i];
		if (rb < 0 || rb >= nrb)
		    error("incompatible dims m=%d, nr=%d, rb=%d",
			  m, nr, rb);
		F77_CALL(dgemm)("N", "N", &nr, &n, &nc,
				&alpha, ax + i * sz, &nr,
				b + j * nc, &ldb,
				&beta, c + rb * nr, &ldc);
	    }
	}
	    
    } else {
	error("Call to cscBlocked_mm must have side == 'L' and transa == 'N'");
    }
}
