#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

SEXP sparse_aggregate(SEXP from, const char *class)
{
	if (class[2] != 'T')
		return from;

	SEXP i0 = GET_SLOT(from, Matrix_iSym);
	if (XLENGTH(i0) < 2)
		return from;
	if (XLENGTH(i0) > INT_MAX)
		Rf_error(_("number of triplets to be aggregated exceeds %s"),
		         "2^31-1");
	PROTECT(i0);
	SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym)), i1, j1, to;
	int *pi0 = INTEGER(i0), *pi1 = NULL,
		*pj0 = INTEGER(j0), *pj1 = NULL,
		*pdim = DIM(from), m = pdim[0], n = pdim[1],
		nnz = (int) XLENGTH(i0), nnz_ = nnz, *iwork = NULL;
	size_t liwork = (size_t) ((int_fast64_t) n + 1 + n + m + nnz),
		lwork = (size_t) nnz;
	Matrix_Calloc(iwork, liwork, int);

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = NULL, *px1 = NULL, *work = NULL; \
		c##IF_NPATTERN( \
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
		px0 = c##PTR(x0); \
		Matrix_Calloc(work, lwork, c##TYPE); \
		); \
		c##tspaggr(pj1, pi1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
		if (nnz != nnz_) { \
		PROTECT(to = newObject(class)); \
		PROTECT(i1 = Rf_allocVector(INTSXP, nnz)), \
		PROTECT(j1 = Rf_allocVector(INTSXP, nnz)); \
		pi1 = INTEGER(i1); \
		pj1 = INTEGER(j1); \
		SET_SLOT(to, Matrix_iSym, i1); \
		SET_SLOT(to, Matrix_jSym, j1); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
		px1 = c##PTR(x1); \
		SET_SLOT(to, Matrix_xSym, x1); \
		); \
		c##tspaggr(pj1, pi1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
		c##IF_NPATTERN( \
		UNPROTECT(1); /* x1 */ \
		); \
		UNPROTECT(3); /* j1, i1, to */ \
		} \
		c##IF_NPATTERN( \
		Matrix_Free(work, lwork); \
		UNPROTECT(1); /* x0 */ \
		); \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

	Matrix_Free(iwork, liwork);
	UNPROTECT(2); /* j0, i0 */

	if (nnz == nnz_)
		return from;

	PROTECT(to);
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);
	UNPROTECT(1); /* to */
	return to;
}

SEXP R_sparse_aggregate(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_aggregate(s_from, class);
}
