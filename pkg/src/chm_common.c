#include "chm_common.h"
#include "Mutils.h"

static R_INLINE int
check_class(char *class, char **valid)
{
    int ans;
    for (ans = 0; ; ans++) {
	if (!strlen(valid[ans])) return -1;
	if (!strcmp(class, valid[ans])) return ans;
    }
}
	
cholmod_sparse *as_cholmod_sparse(SEXP x)
{
    cholmod_sparse *ans = (cholmod_sparse*) malloc(sizeof(cholmod_sparse));
    char *valid[] = {"dgCMatrix", "dsCMatrix", "dtCMatrix",
		     "lgCMatrix", "lsCMatrix", "ltCMatrix",
		     "zgCMatrix", "zsCMatrix", "ztCMatrix",
		     ""};
    int *dims, ctype = check_class(CHAR(asChar(getAttrib(x, R_ClassSymbol))),
				   valid);
    SEXP islot;

    if (ctype < 0) error("invalid class of object to as_cholmod_sparse");
				/* characteristics of the system */
    ans->itype = CHOLMOD_INT;
    ans->dtype = CHOLMOD_DOUBLE;
    ans->packed = TRUE;
    ans->sorted = TRUE;
    ans->x = ans->z = ans->nz = (void *) NULL;
				/* dimensions and nzmax */
    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    ans->nrow = dims[0];
    ans->ncol = dims[1];
    islot = GET_SLOT(x, Matrix_iSym);
    ans->nzmax = LENGTH(islot);
				/* slots always present */
    ans->i = (void *) INTEGER(islot);
    ans->p = (void *) INTEGER(GET_SLOT(x, Matrix_pSym));
				/* set the xtype and any elements */
    switch(ctype) {
    case 0:
    case 1:
    case 2:
	ans->xtype = CHOLMOD_REAL;
	break;
    case 3:
    case 4:
    case 5:
	ans->xtype = CHOLMOD_PATTERN;
	ans->x = (void *) REAL(GET_SLOT(x, Matrix_xSym));
	break;
    case 6:
    case 7:
    case 8:
	ans->xtype = CHOLMOD_COMPLEX;
	ans->x = (void *) COMPLEX(GET_SLOT(x, Matrix_xSym));
	break;
    }
				/* set the stype */
    switch(ctype % 3) {
    case 0: ans->stype = 0; break;
    case 1:
	ans->stype =
	    (!strcmp(CHAR(asChar(getAttrib(x, Matrix_uploSym))), "U")) ?
	    1 : -1;
	break;
    case 2: error("triangular matrices not yet mapped");
    }
    
    return ans;
}

	    
	
	    
