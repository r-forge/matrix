#include "Mdefines.h"
#include "cs-etc.h"

cs *dgC2cs(SEXP obj, int values)
{
	cs *A = (cs *) R_alloc(1, sizeof(cs));
	memset(A, 0, sizeof(cs));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	A->m = INTEGER(dim)[0];
	A->n = INTEGER(dim)[1];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->nzmax = LENGTH(i);
	A->nz = -1;
	if (values && HAS_SLOT(obj, Matrix_xSym)) {
		SEXP x = GET_SLOT(obj, Matrix_xSym);
		switch (TYPEOF(x)) {
		case CPLXSXP:
			A->x = (double *) COMPLEX(x);
			break;
		case REALSXP:
			A->x = REAL(x);
			break;
		default:
			ERROR_INVALID_TYPE(x, __func__);
			break;
		}
	}
	UNPROTECT(3);
	return A;
}

SEXP cs2dgC(const cs *A, const char *class, int values)
{
	int nnz = ((int *) A->p)[A->n];
	R_xlen_t np1 = (R_xlen_t) A->n + 1;
	SEXP obj = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(allocVector(INTSXP, np1)),
		i = PROTECT(allocVector(INTSXP, nnz));
	INTEGER(dim)[0] = A->m;
	INTEGER(dim)[1] = A->n;
	Matrix_memcpy(INTEGER(p), A->p, np1, sizeof(int));
	Matrix_memcpy(INTEGER(i), A->i, nnz, sizeof(int));
	SET_SLOT(obj, Matrix_pSym, p);
	SET_SLOT(obj, Matrix_iSym, i);
	if (values && A->x) {
		SEXP x;
		if (A->nz == -2) {
			PROTECT(x = allocVector(CPLXSXP, nnz));
			Matrix_memcpy(COMPLEX(x), A->x, nnz, sizeof(Rcomplex));
		} else {
			PROTECT(x = allocVector(REALSXP, nnz));
			Matrix_memcpy(REAL(x), A->x, nnz, sizeof(double));
		}
		SET_SLOT(obj, Matrix_xSym, x);
		UNPROTECT(1);
	}
	UNPROTECT(4);
	return obj;
}
