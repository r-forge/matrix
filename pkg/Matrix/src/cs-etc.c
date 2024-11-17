#include "Mdefines.h"
#include "cs-etc.h"

int Matrix_cs_xtype; /* flag indicating use of cs_di_* or cs_ci_* */

Matrix_cs *M2CXS(SEXP obj, int values)
{
	const char *class = Matrix_class(obj, valid_sparse_compressed, 6, __func__);
	values = values && (class[0] == 'd' || class[0] == 'z');
	char trans = (class[2] == 'C') ? 'N' : 'C';
	Matrix_cs *A = (Matrix_cs *) R_alloc(1, sizeof(Matrix_cs));
	memset(A, 0, sizeof(Matrix_cs));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, (trans == 'N') ? Matrix_iSym : Matrix_jSym));
	A->nzmax = LENGTH(i);
	A->m = INTEGER(dim)[(trans == 'N') ? 0 : 1];
	A->n = INTEGER(dim)[(trans == 'N') ? 1 : 1];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->nz = -1;
	A->xtype = CXSPARSE_PATTERN;
	if (values) {
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	switch (class[0]) {
	case 'd':
		A->x = REAL(x);
		A->xtype = CXSPARSE_REAL;
		break;
	case 'z':
		A->x = COMPLEX(x);
		A->xtype = CXSPARSE_COMPLEX;
		break;
	default:
		break;
	}
	}
	UNPROTECT(3);
	return A;
}

SEXP CXS2M(Matrix_cs *A, int values, char shape)
{
	values = values && (A->xtype == CXSPARSE_REAL || A->xtype == CXSPARSE_COMPLEX);
	char class[] = "..CMatrix";
	class[0] = (!values) ? 'n' : ((A->xtype == CXSPARSE_REAL) ? 'd' : 'z');
	class[1] = shape;
	int nnz = A->p[A->n];
	SEXP obj = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) A->n + 1)),
		i = PROTECT(Rf_allocVector(INTSXP, nnz));
	INTEGER(dim)[0] = A->m;
	INTEGER(dim)[1] = A->n;
	memcpy(INTEGER(p), A->p, sizeof(int) * ((size_t) A->n + 1));
	memcpy(INTEGER(i), A->i, sizeof(int) * (size_t) nnz);
	SET_SLOT(obj, Matrix_pSym, p);
	SET_SLOT(obj, Matrix_iSym, i);
	if (values) {
	SEXP x;
	if (A->xtype == CXSPARSE_REAL) {
		PROTECT(x = Rf_allocVector(REALSXP, nnz));
		memcpy(REAL(x), A->x, sizeof(double) * (size_t) nnz);
	} else {
		PROTECT(x = Rf_allocVector(CPLXSXP, nnz));
		memcpy(COMPLEX(x), A->x, sizeof(Rcomplex) * (size_t) nnz);
	}
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(1);
	}
	UNPROTECT(4);
	return obj;
}

/* Wrappers for the functions that we use at least once : */

Matrix_csd *Matrix_cs_dfree(Matrix_csd *D)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_csd *) cs_ci_dfree((cs_cid *) D);
	else
	return (Matrix_csd *) cs_di_dfree((cs_did *) D);
#else
	return (Matrix_csd *)    cs_dfree((csd    *) D);
#endif
}

Matrix_csd *Matrix_cs_dmperm(const Matrix_cs *A, int seed)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_csd *) cs_ci_dmperm((cs_ci *) A, seed);
	else
	return (Matrix_csd *) cs_di_dmperm((cs_di *) A, seed);
#else
	return (Matrix_csd *)    cs_dmperm((cs    *) A, seed);
#endif
}

int Matrix_cs_dropzeros(Matrix_cs *A)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_dropzeros((cs_ci *) A);
	else
	return cs_di_dropzeros((cs_di *) A);
#else
	return    cs_dropzeros((cs    *) A);
#endif
}

void *Matrix_cs_free(void *p)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_free(p);
	else
	return cs_di_free(p);
#else
	return    cs_free(p);
#endif
}

int Matrix_cs_happly(const Matrix_cs *V, int i, double beta, void *x)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_happly((cs_ci *) V, i, beta, (double _Complex *) x);
	else
	return cs_di_happly((cs_di *) V, i, beta, (double          *) x);
#else
	return    cs_happly((cs    *) V, i, beta, (double          *) x);
#endif
}

int Matrix_cs_ipvec(const int *p, const void *b, void *x, int n)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_ipvec(p, (double _Complex *) b, (double _Complex *) x, n);
	else
	return cs_di_ipvec(p, (double          *) b, (double          *) x, n);
#else
	return    cs_ipvec(p, (double          *) b, (double          *) x, n);
#endif
}

int Matrix_cs_lsolve(const Matrix_cs *L, void *x)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_lsolve((cs_ci *) L, (double _Complex *) x);
	else
	return cs_di_lsolve((cs_di *) L, (double          *) x);
#else
	return    cs_lsolve((cs    *) L, (double          *) x);
#endif
}

Matrix_csn *Matrix_cs_lu(const Matrix_cs *A, const Matrix_css *S, double tol)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_csn *) cs_ci_lu((cs_ci *) A, (cs_cis *) S, tol);
	else
	return (Matrix_csn *) cs_di_lu((cs_di *) A, (cs_dis *) S, tol);
#else
	return (Matrix_csn *)    cs_lu((cs    *) A, (css    *) S, tol);
#endif
}

int Matrix_cs_lusol(int order, const Matrix_cs *A, void *b, double tol)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_lusol(order, (cs_ci *) A, (double _Complex *) b, tol);
	else
	return cs_di_lusol(order, (cs_di *) A, (double          *) b, tol);
#else
	return    cs_lusol(order, (cs    *) A, (double          *) b, tol);
#endif
}

Matrix_csn *Matrix_cs_nfree(Matrix_csn *N)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_csn *) cs_ci_nfree((cs_cin *) N);
	else
	return (Matrix_csn *) cs_di_nfree((cs_din *) N);
#else
	return (Matrix_csn *)    cs_nfree((csn    *) N);
#endif
}

Matrix_cs *Matrix_cs_permute(const Matrix_cs *A, const int *pinv, const int *q, int values)
{
	Matrix_cs *B = NULL;
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX) {
	cs_ci *tmp = cs_ci_permute((cs_ci *) A, pinv, q, values);
	B = (Matrix_cs *) cs_ci_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs_ci));
	tmp = cs_ci_free(tmp);
	} else {
	cs_di *tmp = cs_di_permute((cs_di *) A, pinv, q, values);
	B = (Matrix_cs *) cs_di_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs_di));
	tmp = cs_di_free(tmp);
	}
#else
	cs    *tmp =    cs_permute((cs    *) A, pinv, q, values);
	B = (Matrix_cs *)    cs_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs   ));
	tmp =    cs_free(tmp);
#endif
	B->xtype = CXSPARSE_XTYPE_GET();
	return B;
}

int *Matrix_cs_pinv(const int *p, int n)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_pinv(p, n);
	else
	return cs_di_pinv(p, n);
#else
	return    cs_pinv(p, n);
#endif
}

int Matrix_cs_pvec(const int *p, const void *b, void *x, int n)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_pvec(p, (double _Complex *) b, (double _Complex *) x, n);
	else
	return cs_di_pvec(p, (double          *) b, (double          *) x, n);
#else
	return    cs_pvec(p, (double          *) b, (double          *) x, n);
#endif
}

Matrix_csn *Matrix_cs_qr(const Matrix_cs *A, const Matrix_css *S)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_csn *) cs_ci_qr((cs_ci *) A, (cs_cis *) S);
	else
	return (Matrix_csn *) cs_di_qr((cs_di *) A, (cs_dis *) S);
#else
	return (Matrix_csn *)    cs_qr((cs    *) A, (css    *) S);
#endif
}

int Matrix_cs_qrsol(int order, const Matrix_cs *A, void *b)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_qrsol(order, (cs_ci *) A, (double _Complex *) b);
	else
	return cs_di_qrsol(order, (cs_di *) A, (double          *) b);
#else
	return    cs_qrsol(order, (cs    *) A, (double          *) b);
#endif
}

Matrix_css *Matrix_cs_sfree(Matrix_css *S)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_css *) cs_ci_sfree((cs_cis *) S);
	else
	return (Matrix_css *) cs_di_sfree((cs_dis *) S);
#else
	return (Matrix_css *)    cs_sfree((css    *) S);
#endif
}

Matrix_cs *Matrix_cs_spalloc(int m, int n, int nzmax, int values, int triplet)
{
	Matrix_cs *B = NULL;
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX) {
	cs_ci *tmp = cs_ci_spalloc(m, n, nzmax, values, triplet);
	B = (Matrix_cs *) cs_ci_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs_ci));
	tmp = cs_ci_free(tmp);
	} else {
	cs_di *tmp = cs_di_spalloc(m, n, nzmax, values, triplet);
	B = (Matrix_cs *) cs_di_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs_di));
	tmp = cs_di_free(tmp);
	}
#else
	cs    *tmp =    cs_spalloc(m, n, nzmax, values, triplet);
	B = (Matrix_cs *)    cs_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs   ));
	tmp =    cs_free(tmp);
#endif
	B->xtype = CXSPARSE_XTYPE_GET();
	return B;
}

Matrix_cs *Matrix_cs_speye(int m, int n, int values, int triplet)
{
	int k, d = (m < n) ? m : n;
	Matrix_cs *B = Matrix_cs_spalloc(m, n, d, values, triplet);
	if (!B)
		return NULL;
	int *B__p = B->p, *B__i = B->i;
	for (k = 0; k < d; ++k)
		B__p[k] = B__i[k] = k;
	if (!triplet)
		for (k = d; k <= n; ++k)
			B__p[k] = d;
	if (values) {
#ifdef CXSPARSE
		if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX) {
		double _Complex *B__x = (double _Complex *) B->x;
		for (k = 0; k < d; ++k)
			B__x[k] = 1.0;
		} else {
#endif
		double          *B__x = (double          *) B->x;
		for (k = 0; k < d; ++k)
			B__x[k] = 1.0;
#ifdef CXSPARSE
		}
#endif
	}
	return B;
}

Matrix_cs *Matrix_cs_spfree(Matrix_cs *A)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_cs *) cs_ci_spfree((cs_ci *) A);
	else
	return (Matrix_cs *) cs_di_spfree((cs_di *) A);
#else
	return (Matrix_cs *)    cs_spfree((cs    *) A);
#endif
}

int Matrix_cs_sprealloc(Matrix_cs *A, int nzmax)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_sprealloc((cs_ci *) A, nzmax);
	else
	return cs_di_sprealloc((cs_di *) A, nzmax);
#else
	return    cs_sprealloc((cs    *) A, nzmax);
#endif
}

int Matrix_cs_spsolve(Matrix_cs *L, const Matrix_cs *B, int k, int *xi, void *x, const int *pinv, int lo)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_spsolve((cs_ci *) L, (cs_ci *) B, k, xi, (double _Complex *) x, pinv, lo);
	else
	return cs_di_spsolve((cs_di *) L, (cs_di *) B, k, xi, (double          *) x, pinv, lo);
#else
	return    cs_spsolve((cs    *) L, (cs    *) B, k, xi, (double          *) x, pinv, lo);
#endif
}

Matrix_css *Matrix_cs_sqr(int order, const Matrix_cs *A, int qr)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return (Matrix_css *) cs_ci_sqr(order, (cs_ci *) A, qr);
	else
	return (Matrix_css *) cs_di_sqr(order, (cs_di *) A, qr);
#else
	return (Matrix_css *)    cs_sqr(order, (cs    *) A, qr);
#endif
}

Matrix_cs *Matrix_cs_transpose(const Matrix_cs *A, int values)
{
	Matrix_cs *B = NULL;
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX) {
	cs_ci *tmp = cs_ci_transpose((cs_ci *) A, values);
	B = (Matrix_cs *) cs_ci_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs_ci));
	tmp = cs_ci_free(tmp);
	} else {
	cs_di *tmp = cs_di_transpose((cs_di *) A, values);
	B = (Matrix_cs *) cs_di_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs_di));
	tmp = cs_di_free(tmp);
	}
#else
	cs    *tmp =    cs_transpose((cs    *) A, values);
	B = (Matrix_cs *)    cs_calloc(1, sizeof(Matrix_cs));
	memcpy(B, tmp, sizeof(cs   ));
	tmp =    cs_free(tmp);
#endif
	B->xtype = CXSPARSE_XTYPE_GET();
	return B;
}

int Matrix_cs_usolve(const Matrix_cs *U, void *x)
{
#ifdef CXSPARSE
	if (CXSPARSE_XTYPE_GET() == CXSPARSE_COMPLEX)
	return cs_ci_usolve((cs_ci *) U, (double _Complex *) x);
	else
	return cs_di_usolve((cs_di *) U, (double          *) x);
#else
	return    cs_usolve((cs    *) U, (double          *) x);
#endif
}
