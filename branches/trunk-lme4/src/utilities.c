#include "utilities.h"

/** 
 * Replace the value of a slot or subslot of an object in place.  This
 * routine purposely does not copy the value of obj.  Use with caution.
 * 
 * @param obj object with slot to be replaced
 * @param names vector of names.  The last element is the name of the slot to replace.  The leading elements are the names of slots and subslots of obj.
 * @param value the replacement value for the slot
 * 
 * @return obj, with the named slot modified in place.
 */
SEXP
nlme_replaceSlot(SEXP obj, SEXP names, SEXP value)
{
    int lnm1 = length(names) - 1;

    if (lnm1 >= 0) {
	SEXP comp = obj;
	int i;

	for (i = 0; i < lnm1; i++) {
	    comp = GET_SLOT(comp, install(CHAR(STRING_ELT(names, i))));
	}
	SET_SLOT(comp, install(CHAR(STRING_ELT(names, lnm1))), value);
    }
    return obj;
}

/** 
 * Produce a weighted copy of the matrices in MLin in the storage
 * allocated to MLout
 * 
 * @param MLin input matrix list
 * @param wts real vector of weights
 * @param adjst adjusted response
 * @param MLout On input a list of matrices of the same dimensions as MLin.  
 * 
 * @return MLout with its contents overwritten by a weighted copy of
 * MLin according to wts with adjst overwriting the response.
 */
SEXP nlme_weight_matrix_list(SEXP MLin, SEXP wts, SEXP adjst, SEXP MLout)
{
    int i, j, n, nf;
    SEXP lastM;
    
    if (!(isNewList(MLin) && isReal(wts) && isReal(adjst) && isNewList(MLout)))
	error("Incorrect argument type");
    nf = length(MLin);
    if (length(MLout) != nf)
	error("Lengths of MLin (%d) and MLout (%d) must match", nf,
	      length(MLout));
    n = length(wts);
    if (length(adjst) != n)
	error("Expected adjst to have length %d, got %d", n, length(adjst));
    for (i = 0; i < nf; i++) {
	SEXP Min = VECTOR_ELT(MLin, i),
	    Mout = VECTOR_ELT(MLout, i);
	int *din, *dout, k, nc;

	if (!(isMatrix(Min) && isReal(Min)))
	    error("component %d of MLin is not a numeric matrix", i + 1);
	din = INTEGER(getAttrib(Min, R_DimSymbol));
	nc = din[1];
	if (din[0] != n)
	    error("component %d of MLin has %d rows, expected %d", i + 1,
		  din[0], n);
	if (!(isMatrix(Mout) && isReal(Mout)))
	    error("component %d of MLout is not a numeric matrix", i + 1);
	dout = INTEGER(getAttrib(Mout, R_DimSymbol));
	if (dout[0] != n)
	    error("component %d of MLout has %d rows, expected %d", i + 1,
		  dout[0], n);
	if (dout[1] != nc)
	    error("component %d of MLout has %d columns, expected %d", i + 1,
		  dout[1], nc);
	for (k = 0; k < nc; k++) {
	    for (j = 0; j < n; j++) {
		REAL(Mout)[j + k * n] = REAL(Min)[j + k * n] * REAL(wts)[j];
	    }
	}
    }
    lastM = VECTOR_ELT(MLout, nf - 1);
    j = INTEGER(getAttrib(lastM, R_DimSymbol))[1] - 1;
    for (i = 0; i < n; i++)
	REAL(lastM)[j*n + i] = REAL(adjst)[i] * REAL(wts)[i];
    return MLout;
}
