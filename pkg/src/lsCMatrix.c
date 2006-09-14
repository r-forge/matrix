#include "lsCMatrix.h"

/**
 * Check the validity of the slots of an lsCMatrix object
 *
 * @param x Pointer to an lsCMatrix object
 *
 * @return an SEXP that is either TRUE or a character string
 * describing the way in which the object failed the validity check
 */
SEXP lsCMatrix_validate(SEXP x)
{
    SEXP val = symmetricMatrix_validate(x);
    if(isString(val))
	return(val);
    else {
	/* FIXME needed? ltC* inherits from lgC* which does this in validate*/
	SEXP pslot = GET_SLOT(x, Matrix_pSym),
	    islot = GET_SLOT(x, Matrix_iSym);
	/* column sorting now done in Csparse_validate */
/* 	int */
/* 	    ncol = length(pslot) - 1, */
/* 	    *xp = INTEGER(pslot), */
/* 	    *xi = INTEGER(islot); */
/* 	if (csc_unsorted_columns(ncol, xp, xi)) */
/* 	    csc_sort_columns(ncol, xp, xi, (double *) NULL); */

	return ScalarLogical(1);
    }
}
