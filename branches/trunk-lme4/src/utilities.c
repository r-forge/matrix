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
nlme_replaceSlot(SEXP obj, const SEXPREC* names, const SEXPREC* value)
{
    int lnm1 = length((SEXP) names) - 1;

    if (lnm1 >= 0) {
	SEXP comp = obj;
	int i;

	for (i = 0; i < lnm1; i++) {
	    comp = GET_SLOT(comp, install(CHAR(STRING_ELT((SEXP) names, i))));
	}
	SET_SLOT(comp, install(CHAR(STRING_ELT((SEXP) names, lnm1))),
		 (SEXP) value);
    }
    return obj;
}

