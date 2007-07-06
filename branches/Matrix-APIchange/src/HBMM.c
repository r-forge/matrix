#include "HBMM.h"

#include "mmio.h"

SEXP Matrix_writeMatrixMarket(SEXP obj, SEXP file, SEXP typep)
{
    const char *type = CHAR(asChar(typep));
    char *ff = strdup(CHAR(asChar(file)));
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	*ii = (int *) NULL, *jj = (int *) NULL;
    int M = dims[0], N = dims[1], i, nz = -1, *src;
    MM_typecode matcode;
    double *xx = (double *) NULL;

    mm_set_matrix(&matcode);
    if (type[2] == 'C' || type[2] == 'T') {
	SEXP islot = GET_SLOT(obj, Matrix_iSym);
	nz = LENGTH(islot);
	src = INTEGER(islot);
	ii = Alloca(nz, int);
	R_CheckStack();
	for (i = 0; i < nz; i++) ii[i] = src[i] + 1;
	mm_set_coordinate(&matcode);
    } else error("Only types 'C' and 'T' allowed");

    if (type[0] == 'D') {
	xx = REAL(GET_SLOT(obj, Matrix_xSym));
	mm_set_real(&matcode);
    } else error("Only real matrices allowed");

    if (type[1] == 'S') {
	if (*uplo_P(obj) != 'L')
	    error("Symmetric matrices must be stored in lower triangle");
	mm_set_symmetric(&matcode);
    }
    if (type[1] == 'G') mm_set_general(&matcode);

    if (type[2] == 'C') {
	jj = Alloca(nz, int);
	R_CheckStack();
	expand_cmprPt(N, INTEGER(GET_SLOT(obj, Matrix_pSym)), jj);
	for (i = 0; i < nz; i++) jj[i] += 1; /* 1-based indices */
    }
    if (type[2] == 'T') {
	src = INTEGER(GET_SLOT(obj, Matrix_jSym));
	jj = Alloca(nz, int);
	R_CheckStack();
	for (i = 0; i < nz; i++) jj[i] = src[i] + 1;
    }
    if (!jj) error("storage mode must be T or C");

    mm_write_mtx_crd(ff, M, N, nz, ii, jj, xx, matcode);

    free(ff);
    return R_NilValue;

}
