#include "bCrosstab.h"
#include "ldl.h"

/** 
 * Calculate the zero-based index in a row-wise packed lower
 * triangular matrix.  This is used for the arrays of blocked sparse matrices.
 * 
 * @param i row number (0-based)
 * @param k column number (0-based)
 * 
 * @return The 0-based index of the (i,k) element of a row-wise packed lower
 * triangular matrix.
 */    
static R_INLINE
int Lind(int i, int k)
{
    return (i * (i + 1))/2 + k;
}

/** 
 * Permute an index vector
 * 
 * @param i vector of 0-based indices
 * @param nnz length of vector i
 * @param perm 0-based permutation vector of length max(i)
 */
static R_INLINE
void ind_permute(int i[], int nnz, int perm[])
{
    int j;
    for (j = 0; j < nnz; j++) i[j] = perm[i[j]];
}

/** 
 * Force indices to be in the lower triangle of a matrix
 * 
 * @param i vector of 0-based row indices
 * @param j vector of 0-based column indices
 * @param nnz length of index vectors
 */
static R_INLINE
void make_lower_triangular(int i[], int j[], int nnz)
{
    int k;
    for (k = 0; k < nnz; k++) {
	if (i[k] < j[k]) {
	    int tmp = i[k];
	    i[k] = j[k];
	    j[k] = tmp;
	}
    }
}

/** 
 * Create a named list of length n
 * 
 * @param n length of list to return
 * @param names names of list elements
 * 
 * @return pointer to a named list of length n
 */
static
SEXP make_named_list(int n, char **names)
{
    SEXP ans = PROTECT(allocVector(VECSXP, n)),
	nms = PROTECT(allocVector(CHARSXP, n));
    int i;

    for (i = 0; i < n; i++) nms[i] = mkChar(names[i]);
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/** 
 * Update a diagonal block
 * 
 * @param ctab pointer to a blocked crosstabulation object
 * @param jm1 index of updating column
 * @param i index of diagonal block to be updated
 */
static void diag_update(ctab, jm1, i)
{
    SEXP db = VECTOR_ELT(ctab, Lind(i, i)),
	jb = VECTOR_ELT(ctab, Lind(i, jm1));
    SEXP dpp = GET_SLOT(db, Matrix_pSym),
	jpp = GET_SLOT(jb, Matrix_pSym);
    int nci = length(dpp) - 1,
	ncj = length(jpp) - 1,
	*dp = INTEGER(dpp),
	*di = INTEGER(GET_SLOT(db), Matrix_iSym),
	*jp = INTEGER(jpp),
	*ji = INTEGER(GET_SLOT(jb), Matrix_iSym);
	
}

/** 
 * Project a column of a pairwise crosstabulation onto the remaining columns
 * 
 * @param ctab pointer to a blocked crosstabulation object
 * @param i index of column to project
 * 
 * @return a list of the projected ctab, the permutation of factor
 * j+1, and the j+1 diagonal element of L-inverse
 */
SEXP bCrosstab_project(SEXP ctab, SEXP jp)
{
    char anms[][] = {"ctab", "perm", "Linv"};
    SEXP ans = PROTECT(make_named_list(3, anms));
    int ctbl = length(ctab), i, j = asInteger(jp), jm1, k;
    int nf = ((int) (-1 + sqrt((double)(1 + 8*ctbl))))/2;

    if (ctbl != (nf*(nf + 1))/2)
	error("length of ctab = %d is not permisable", ctbl);
    if (j < 1 || j > (nf - 1))
	error("j must be in the range [1,%d]", nf - 1);
    jm1 = j - 1;		/* 0-based index */
    for (i = j; i < nf; i++) {	
	diag_update(ctab, jm1, i);
	for (k = i; k < nf; k++) 
	    offdiag_update(ctab, jm1, i, k);
    }
    
}

/** 
 * Permute the levels of one of the grouping factors in a bCrosstab object
 * 
 * @param ctab Pointer to a bCrosstab object
 * @param i index (1-based) of the factor levels to permute
 * @param perm permutation (0-based) to apply
 * 
 * @return the ctab object with the levels permuted
 */
SEXP bCrosstab_permute(SEXP ctab, SEXP ip, SEXP permp)
{
    SEXP cs, trip, ipt, jpt;
    int *perm, ind = asInteger(ip) - 1,
	ctbl = length(ctab), plen = length(permp), j;
    int nf = ((int) (-1 + sqrt((double)(1 + 8*ctbl))))/2;
    int *iperm = (int *) Calloc(plen, int);
    
    if (ctbl != (nf*(nf + 1))/2)
	error("length of ctab = %d is not permisable", ctbl);
    if (0 < ind || nf <= ind)
	error("i must be in the range [1,%d]", nf);
    if (!isInteger(permp) || plen !=
	(length(GET_SLOT(VECTOR_ELT(ctab, Lind(ind, ind)), Matrix_pSym)) - 1))
	error("perm must be an integer vector of length %d",
	      length(GET_SLOT(VECTOR_ELT(ctab, Lind(ind, ind)), Matrix_pSym)) - 1);
    if (!ldl_valid_perm(plen, INTEGER(permp), iperm))
	error("perm is not a valid 0-based permutation vector");
    for (j = 0; j < ind; j++) {
	trip = csc_to_triplet(VECTOR_ELT(ctab, Lind(ind, j)));
	ipt = GET_SLOT(trip, Matrix_iSym);
	ind_permute(INTEGER(ipt), length(ipt), perm);
	SET_VECTOR_ELT(ctab, Lind(ind, j), triplet_to_csc(trip));
    }
    trip = csc_to_triplet(VECTOR_ELT(ctab, Lind(ind, ind)));
    ipt = GET_SLOT(trip, Matrix_iSym);
    ind_permute(INTEGER(ipt), length(ipt), perm);
    jpt = GET_SLOT(trip, Matrix_jSym);
    ind_permute(INTEGER(jpt), length(jpt), perm);
    make_lower_triangular(INTEGER(ipt), INTEGER(jpt), length(jpt));
    SET_VECTOR_ELT(ctab, Lind(ind, ind), triplet_to_csc(trip));
    for (j = ind + 1; j < nf; j++) {
	trip = csc_to_triplet(VECTOR_ELT(ctab, Lind(j, ind)));
	jpt = GET_SLOT(trip, Matrix_jSym);
	ind_permute(INTEGER(jpt), length(jpt), perm);
	SET_VECTOR_ELT(ctab, Lind(j, ind), triplet_to_csc(trip));
    }
    return ctab;
}

	
