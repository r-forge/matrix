				/* Sparse triangular logical matrices */
#include "ltCMatrix.h"

/**
 * Check the validity of the slots of an ltCMatrix object
 *
 * @param x Pointer to an ltCMatrix object
 *
 * @return an SEXP that is either TRUE or a character string
 * describing the way in which the object failed the validity check
 */
SEXP ltCMatrix_validate(SEXP x)
{
    /* Because ltCMatrix inherits from triangularMatrix this is not necessary */
/*     SEXP val = triangularMatrix_validate(x); */
/*     if(isString(val)) */
/* 	return(val); */
/*     else { */
	/* FIXME needed? ltC* inherits from lgC* which does this in validate*/
/* 	SEXP pslot = GET_SLOT(x, Matrix_pSym), */
/* 	    islot = GET_SLOT(x, Matrix_iSym); */
/* 	int */
/* 	    ncol = length(pslot) - 1, */
/* 	    *xp = INTEGER(pslot), */
/* 	    *xi = INTEGER(islot); */
	/* column sorting now done in Csparse_validate */
/* 	if (csc_unsorted_columns(ncol, xp, xi)) */
/* 	    csc_sort_columns(ncol, xp, xi, (double *) NULL); */
	return ScalarLogical(1);
/*     } */
}
