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

SEXP nlme_weight_matrix_list(SEXP MLin, SEXP wts, SEXP MLout)
{
    int i, n, nf;
    
    if (!(isNewList(MLin) && isReal(wts) && isNewList(MLout)))
	error("Incorrect argument type");
    nf = length(MLin);
    if (length(MLout) != nf)
	error("Lengths of MLin (%d) and MLout (%d) must match", nf,
	      length(MLout));
    n = length(wts);
    for (i = 0; i < nf; i++) {
	SEXP Min = VECTOR_ELT(MLin, i),
	    Mout = VECTOR_ELT(MLout, i);
	int *din, *dout, j, k, nc;

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
    return MLout;
}
