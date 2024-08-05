#include "Mdefines.h"
#include "cholmod-api.h"

/**
 * Coerce from sparseCholesky to (cholmod_factor *)
 *
 * Sets the members of a pointed-to cholmod_factor struct, using "data"
 * obtained from slots of a sparseCholesky.  The result should _not_ be
 * freed using cholmod_free_factor, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param L a pointer to a cholmod_factor struct, to be modified in-place.
 * @param from an S4 object inheriting from virtual class sparseCholesky.
 *
 * @return L.
 */
/* NB: mostly parallel to M2CHF in ./cholmod-etc.c */
cholmod_factor *sexp_as_cholmod_factor(cholmod_factor *L, SEXP from)
{
	static const char *valid[] = {
		"nsimplicialCholesky", "nsupernodalCholesky",
		"dsimplicialCholesky", "dsupernodalCholesky",
		"zsimplicialCholesky", "zsupernodalCholesky", "" };
	const char *class = Matrix_class(from, valid, -1, __func__);
	memset(L, 0, sizeof(cholmod_factor));
	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
		perm = PROTECT(GET_SLOT(from, Matrix_permSym)),
		colcount = PROTECT(GET_SLOT(from, Matrix_colcountSym)),
		ordering = PROTECT(GET_SLOT(from, Matrix_orderingSym));
	L->ordering = INTEGER(ordering)[0];
	L->is_super = class[2] == 'u';
	L->n = (size_t) INTEGER(dim)[0];
	L->minor = L->n;
	if (L->ordering != CHOLMOD_NATURAL)
		L->Perm = INTEGER(perm);
	else {
		/* cholmod_check_factor allows L->Perm == NULL,
		   but cholmod_copy_factor does not test, so it segfaults ...
		*/
		int n = (int) L->n, *Perm = (int *) R_alloc(L->n, sizeof(int));
		for (int j = 0; j < n; ++j)
			Perm[j] = j;
		L->Perm = Perm;
	}
	L->ColCount = INTEGER(colcount);
	if (L->is_super) {
		SEXP maxcsize = PROTECT(GET_SLOT(from, Matrix_maxcsizeSym)),
			maxesize = PROTECT(GET_SLOT(from, Matrix_maxesizeSym)),
			super = PROTECT(GET_SLOT(from, Matrix_superSym)),
			pi = PROTECT(GET_SLOT(from, Matrix_piSym)),
			px = PROTECT(GET_SLOT(from, Matrix_pxSym)),
			s = PROTECT(GET_SLOT(from, Matrix_sSym));
		L->nsuper = (size_t) (LENGTH(super) - 1);
		L->ssize = (size_t) INTEGER(pi)[L->nsuper];
		L->xsize = (size_t) INTEGER(px)[L->nsuper];
		L->maxcsize = (size_t) INTEGER(maxcsize)[0];
		L->maxesize = (size_t) INTEGER(maxesize)[0];
		L->super = INTEGER(super);
		L->pi = INTEGER(pi);
		L->px = INTEGER(px);
		L->s = INTEGER(s);
		L->is_ll = 1;
		L->is_monotonic = 1;
		UNPROTECT(6);
	} else {
		if (class[0] != 'n') {
		SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i = PROTECT(GET_SLOT(from, Matrix_iSym)),
			nz = PROTECT(GET_SLOT(from, Matrix_nzSym)),
			next = PROTECT(GET_SLOT(from, Matrix_nextSym)),
			prev = PROTECT(GET_SLOT(from, Matrix_prevSym)),
			is_ll = PROTECT(GET_SLOT(from, Matrix_isllSym)),
			is_monotonic = PROTECT(GET_SLOT(from, Matrix_ismtSym));
		L->nzmax = (size_t) INTEGER(p)[L->n];
		L->p = INTEGER(p);
		L->i = INTEGER(i);
		L->nz = INTEGER(nz);
		L->next = INTEGER(next);
		L->prev = INTEGER(prev);
		L->is_ll = LOGICAL(is_ll)[0] != 0;
		L->is_monotonic = LOGICAL(is_monotonic)[0] != 0;
		UNPROTECT(7);
		}
	}
	L->itype = CHOLMOD_INT;
	L->xtype = CHOLMOD_PATTERN;
	L->dtype = CHOLMOD_DOUBLE;
	if (class[0] != 'n') {
	SEXP minor = GET_SLOT(from, Matrix_minorSym);
	L->minor = (size_t) INTEGER(minor)[0];
	SEXP x = GET_SLOT(from, Matrix_xSym);
	switch (class[0]) {
	case 'd':
		L->x = REAL(x);
		L->xtype = CHOLMOD_REAL;
		break;
	case 'z':
		L->x = COMPLEX(x);
		L->xtype = CHOLMOD_COMPLEX;
		break;
	default:
		break;
	}
	}
	UNPROTECT(4);
	return L;
}

/**
 * Coerce from [CR]sparseMatrix to (cholmod_sparse *)
 *
 * Sets the members of a pointed-to cholmod_sparse struct, using "data"
 * obtained from slots of a [CR]sparseMatrix.  The result should _not_ be
 * freed using cholmod_free_sparse, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param A a pointer to a cholmod_sparse struct, to be modified in-place.
 * @param from an S4 object inheriting from virtual class [CR]sparseMatrix.
 * @param allocUnit a boolean indicating if the unit diagonal of formally
 *     unit triangular [CR]sparseMatrix should be allocated.
 * @param sortInPlace a boolean indicating if unsorted [CR]sparseMatrix
 *     should be sorted in place to avoid an allocation.
 *
 * @return A.
 */
/* NB: mostly parallel to M2CHS in ./cholmod-etc.c */
cholmod_sparse *sexp_as_cholmod_sparse(cholmod_sparse *A, SEXP from,
                                       Rboolean allocUnit, Rboolean sortInPlace)
{
	/* defined in ./Csparse.c : */
	SEXP checkpi(SEXP dim, SEXP p, SEXP i);

	const char *class = Matrix_class(from, valid_sparse_compressed, 6, __func__);
	memset(A, 0, sizeof(cholmod_sparse));
	int mg = (class[2] == 'C') ? 1 : 0;
	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(from, Matrix_pSym)),
		i = PROTECT(GET_SLOT(from, (mg == 1) ? Matrix_iSym : Matrix_jSym));
	A->nrow = (size_t) INTEGER(dim)[(mg == 1) ? 0 : 1];
	A->ncol = (size_t) INTEGER(dim)[(mg == 1) ? 1 : 0];
	A->nzmax = (size_t) INTEGER(p)[A->ncol];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;
	A->sorted = LOGICAL(checkpi(dim, p, i))[0] != 0;
	A->packed = 1;
	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		A->stype = ((CHAR(STRING_ELT(uplo, 0))[0] == 'U') == (mg == 1)) ? 1 : -1;
	}
	if (!A->sorted && !sortInPlace) {
		void *tmp;
		tmp = A->p;
		A->p = (void *) R_alloc(A->ncol + 1, sizeof(int));
		memcpy(A->p, tmp, sizeof(int) * (A->ncol + 1));
		tmp = A->i;
		A->i = (void *) R_alloc(A->nzmax, sizeof(int));
		memcpy(A->i, tmp, sizeof(int) * A->nzmax);
	}
	if (class[0] != 'n') {
	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	switch (class[0]) {
	case 'l':
	case 'i':
	{
		int *px = (TYPEOF(x) == LGLSXP) ? LOGICAL(x) : INTEGER(x);
		double *Ax = (double *) R_alloc(A->nzmax, sizeof(double));
		for (size_t k = 0; k < A->nzmax; ++k)
			Ax[k] = (px[k] == NA_INTEGER) ? NA_REAL : (double) px[k];
		A->x = Ax;
		A->xtype = CHOLMOD_REAL;
		break;
	}
	case 'd':
		A->x = REAL(x);
		A->xtype = CHOLMOD_REAL;
		if (!A->sorted && !sortInPlace) {
			void *tmp = A->x;
			A->x = (void *) R_alloc(A->nzmax, sizeof(double));
			memcpy(A->x, tmp, sizeof(double) * A->nzmax);
		}
		break;
	case 'z':
		A->x = COMPLEX(x);
		A->xtype = CHOLMOD_COMPLEX;
		if (!A->sorted && !sortInPlace) {
			void *tmp = A->x;
			A->x = (void *) R_alloc(A->nzmax, sizeof(Rcomplex));
			memcpy(A->x, tmp, sizeof(Rcomplex) * A->nzmax);
		}
		break;
	default:
		break;
	}
	UNPROTECT(1); /* x */
	}
	if (!A->sorted && !cholmod_sort(A, &c))
		error(_("'%s' failed in '%s'"), "cholmod_sort", __func__);
	if (class[1] == 't' && A->ncol > 0 && allocUnit) {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		char nu = CHAR(STRING_ELT(diag, 0))[0];
		if (nu != 'N') {
			double one[] = { 1.0, 0.0 };
			cholmod_sparse *I1 = cholmod_speye(A->nrow, A->ncol, A->xtype, &c);
			I1->stype = A->stype;
			cholmod_sparse *A1 = cholmod_add(A, I1, one, one, 1, 1, &c);
			A->nzmax = (size_t) ((int *) A1->p)[A->ncol];
			A->p = (void *) R_alloc(A->ncol + 1, sizeof(int));
			memcpy(A->p, A1->p, sizeof(int) * (A->ncol + 1));
			A->i = (void *) R_alloc(A->nzmax, sizeof(int));
			memcpy(A->i, A1->i, sizeof(int) * A->nzmax);
			if (A->xtype != CHOLMOD_PATTERN) {
			if (A->xtype == CHOLMOD_REAL) {
			A->x = (void *) R_alloc(A->nzmax, sizeof(double));
			memcpy(A->x, A1->x, sizeof(double) * A->nzmax);
			} else {
			A->x = (void *) R_alloc(A->nzmax, sizeof(Rcomplex));
			memcpy(A->x, A1->x, sizeof(double) * A->nzmax);
			}
			}
			cholmod_free_sparse(&I1, &c);
			cholmod_free_sparse(&A1, &c);
		}
	}
	UNPROTECT(3); /* i, p, dim */
	return A;
}

/**
 * Coerce from TsparseMatrix to (cholmod_triplet *)
 *
 * Sets the members of a pointed-to cholmod_triplet struct, using "data"
 * obtained from slots of a TsparseMatrix.  The result should _not_ be
 * freed using cholmod_free_sparse, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param A a pointer to a cholmod_triplet struct, to be modified in-place.
 * @param from an S4 object inheriting from virtual class TsparseMatrix.
 * @param allocUnit a boolean indicating if the unit diagonal of formally
 *     unit triangular TsparseMatrix should be allocated.
 *
 * @return A.
 */
cholmod_triplet *sexp_as_cholmod_triplet(cholmod_triplet *A, SEXP from,
                                         Rboolean allocUnit)
{
	const char *class = Matrix_class(from, valid_sparse_triplet, 6, __func__);
	memset(A, 0, sizeof(cholmod_triplet));
	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym)),
		i = PROTECT(GET_SLOT(from, Matrix_pSym)),
		j = PROTECT(GET_SLOT(from, Matrix_iSym));;
	A->nrow = (size_t) INTEGER(dim)[0];
	A->ncol = (size_t) INTEGER(dim)[1];
	A->nzmax = A->nnz = (size_t) XLENGTH(i);
	A->i = INTEGER(i);
	A->j = INTEGER(j);
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;
	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		A->stype = (CHAR(STRING_ELT(uplo, 0))[0] == 'U') ? 1 : -1;
	}
	if (class[0] != 'n') {
	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	switch (class[0]) {
	case 'l':
	case 'i':
	{
		int *px = (TYPEOF(x) == LGLSXP) ? LOGICAL(x) : INTEGER(x);
		double *Ax = (double *) R_alloc(A->nzmax, sizeof(double));
		for (size_t k = 0; k < A->nzmax; ++k)
			Ax[k] = (px[k] == NA_INTEGER) ? NA_REAL : (double) px[k];
		A->x = Ax;
		A->xtype = CHOLMOD_REAL;
		break;
	}
	case 'd':
		A->x = REAL(x);
		A->xtype = CHOLMOD_REAL;
		break;
	case 'z':
		A->x = COMPLEX(x);
		A->xtype = CHOLMOD_COMPLEX;
		break;
	default:
		break;
	}
	UNPROTECT(1); /* x */
	}
	if (class[1] == 't' && A->ncol > 0 && allocUnit) {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		char nu = CHAR(STRING_ELT(diag, 0))[0];
		if (nu != 'N') {
			A->nzmax += A->ncol;
			void *tmp;
			tmp = A->i;
			A->i = (void *) R_alloc(A->nzmax, sizeof(int));
			memcpy(A->i, tmp, sizeof(int) * A->nnz);
			tmp = A->j;
			A->j = (void *) R_alloc(A->nzmax, sizeof(int));
			memcpy(A->j, tmp, sizeof(int) * A->nnz);
			if (A->xtype != CHOLMOD_PATTERN) {
			tmp = A->x;
			if (A->xtype == CHOLMOD_REAL) {
			A->x = (void *) R_alloc(A->nzmax, sizeof(double));
			memcpy(A->x, tmp, sizeof(double) * A->nnz);
			} else {
			A->x = (void *) R_alloc(A->nzmax, sizeof(Rcomplex));
			memcpy(A->x, tmp, sizeof(Rcomplex) * A->nnz);
			}
			}
			int n = (int) A->ncol,
				*Ai = ((int *) A->i) + A->nnz,
				*Aj = ((int *) A->j) + A->nnz;
			for (int j = 0; j < n; ++j)
				Ai[j] = Aj[j] = j;
			if (A->xtype != CHOLMOD_PATTERN) {
			if (A->xtype == CHOLMOD_REAL) {
			double *Ax = ((double *) A->x) + A->nnz;
			for (int j = 0; j < n; ++j)
				Ax[j] = 1.0;
			} else {
			Rcomplex *Ax = ((Rcomplex *) A->x) + A->nnz;
			for (int j = 0; j < n; ++j)
				Ax[j] = Matrix_zunit;
			}
			}
			A->nnz += A->ncol;
		}
	}
	UNPROTECT(3); /* j, i, dim */
	return A;
}

/**
 * Coerce from .geMatrix or vector to (cholmod_dense *)
 *
 * Sets the members of a pointed-to cholmod_dense struct, using "data"
 * obtained from slots of a .geMatrix.  The result should _not_ be
 * freed using cholmod_free_dense, as the resulting members point to
 * memory controlled by R, not by CHOLMOD.
 *
 * @param A a pointer to a cholmod_dense struct, to be modified in-place.
 * @param from an S4 object inheriting from class .geMatrix _or_
 *     a traditional vector of type "logical", "integer", "double", or
 *     "complex" (to be handled as a 1-column matrix if not a matrix).
 *
 * @return A.
 */
/* NB: mostly parallel to M2CHD in ./cholmod-etc.c */
cholmod_dense *sexp_as_cholmod_dense(cholmod_dense *A, SEXP from)
{
	static const char *valid[] = {
		"ngeMatrix", "lgeMatrix", "igeMatrix", "dgeMatrix", "zgeMatrix", "" };
	const char *class = Matrix_class(from, valid, 0, NULL);
	int m, n;
	if (class) {
		SEXP dim = GET_SLOT(from, Matrix_DimSym);
		m = INTEGER(dim)[0];
		n = INTEGER(dim)[1];
		from = GET_SLOT(from, Matrix_xSym);
	} else {
		switch (TYPEOF(from)) {
		case LGLSXP:
		case INTSXP:
		case REALSXP:
		case CPLXSXP:
			break;
		default:
			ERROR_INVALID_TYPE(from, __func__);
			break;
		}
		SEXP dim = getAttrib(from, R_DimSymbol);
		if (TYPEOF(dim) == INTSXP && LENGTH(dim) == 2) {
			m = INTEGER(dim)[0];
			n = INTEGER(dim)[1];
		} else {
			m = LENGTH(from);
			n = 1;
		}
	}
	PROTECT(from);
	memset(A, 0, sizeof(cholmod_dense));
	A->nrow = (size_t) m;
	A->ncol = (size_t) n;
	A->nzmax = A->nrow * A->ncol;
	A->d = A->nrow;
	A->dtype = CHOLMOD_DOUBLE;
	switch (TYPEOF(from)) {
	case LGLSXP:
	case INTSXP:
	{
		int pattern = class[0] == 'n';
		int *px = (TYPEOF(from) == LGLSXP) ? LOGICAL(from) : INTEGER(from);
		double *Ax = (double *) R_alloc(A->nzmax, sizeof(double));
		for (size_t k = 0; k < A->nzmax; ++k)
			Ax[k] = (px[k] == NA_INTEGER)
				? ((pattern) ? 1.0 : NA_REAL) : (double) px[k];
		A->x = Ax;
		A->xtype = CHOLMOD_REAL;
		break;
	}
	case REALSXP:
		A->x = REAL(from);
		A->xtype = CHOLMOD_REAL;
		break;
	case CPLXSXP:
		A->x = COMPLEX(from);
		A->xtype = CHOLMOD_COMPLEX;
		break;
	default:
		break;
	}
	UNPROTECT(1); /* from */
	return A;
}

/**
 * Coerce from (double *) to (cholmod_dense *) with given dimensions
 *
 * An analogue of base::matrix(data, nrow, ncol),
 * where typeof(data)=="double" and length(data)==nrow*ncol.
 *
 * @param A a pointer to a cholmod_dense struct, to be modified in-place.
 * @param data a pointer to an nrow*ncol*sizeof(double) block of memory.
 * @param nrow the desired number of rows.
 * @param ncol the desired number of columns.
 *
 * @return A.
 */
cholmod_dense *numeric_as_cholmod_dense(cholmod_dense *A,
                                        double *data, int nrow, int ncol)
{
	memset(A, 0, sizeof(cholmod_dense));
	A->nrow = (size_t) nrow;
	A->ncol = (size_t) ncol;
	A->nzmax = A->nrow * A->ncol;
	A->d = A->nrow;
	A->x = data;
	A->xtype = CHOLMOD_REAL;
	A->dtype = CHOLMOD_DOUBLE;
	return A;
}

/**
 * Coerce from (cholmod_factor *) to sparseCholesky
 *
 * Allocates an S4 object inheriting from virtual class sparseCholesky
 * and copies into the slots from members of a pointed-to cholmod_factor
 * struct.  The specific class of the result is determined by struct
 * members xtype and is_super.
 *
 * @param L a pointer to a cholmod_factor struct.
 * @param doFree a flag indicating if and how to free L before returning.
 *     (0) don't free, (>0) free with cholmod_free_factor, (<0) free with
 *     R_Free.
 *
 * @return A sparseCholesky.
 */
/* NB: mostly parallel to CHF2M in ./cholmod-etc.c */
SEXP cholmod_factor_as_sexp(cholmod_factor *L, int doFree)
{

#define errorFree(...) \
	do { \
		MAYBE_FREE; \
		error(__VA_ARGS__); \
	} while (0)

#define MAYBE_FREE \
	do { \
		if (doFree != 0) { \
			if (doFree < 0) \
				R_Free(L); \
			else if (L->itype == CHOLMOD_INT) \
				cholmod_free_factor(&L, &c); \
			else \
				cholmod_l_free_factor(&L, &cl); \
		} \
	} while (0)

	if (L->itype != CHOLMOD_INT)
		errorFree(_("wrong '%s'"), "itype");
	if (L->xtype != CHOLMOD_PATTERN &&
	    L->xtype != CHOLMOD_REAL && L->xtype != CHOLMOD_COMPLEX)
		errorFree(_("wrong '%s'"), "xtype");
	if (L->dtype != CHOLMOD_DOUBLE)
		errorFree(_("wrong '%s'"), "dtype");
	if (L->n > INT_MAX)
		errorFree(_("dimensions cannot exceed %s"), "2^31-1");
	if (L->super) {
		if (L->maxcsize > INT_MAX)
		errorFree(_("'%s' would overflow type \"%s\""),
			      "maxcsize", "integer");
	} else {
		if (L->n == INT_MAX)
		errorFree(_("n+1 would overflow type \"%s\""),
			      "integer");
	}
	if (L->minor < L->n) {
		if (L->is_ll)
		errorFree(_("leading principal minor of order %d is not positive"),
			      (int) L->minor + 1);
		else
		errorFree(_("leading principal minor of order %d is zero"),
			      (int) L->minor + 1);
	}
	char class[] = "...........Cholesky";
	class[0] = (L->xtype == CHOLMOD_PATTERN)
		? 'n' : ((L->xtype == CHOLMOD_REAL) ? 'd' : 'z');
	memcpy(class + 1, (L->is_super) ? "supernodal" : "simplicial", 10);
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym)),
		ordering = PROTECT(GET_SLOT(to, Matrix_orderingSym));
	INTEGER(ordering)[0] = L->ordering;
	INTEGER(dim)[0] = INTEGER(dim)[1] = (int) L->n;
	if (L->ordering != CHOLMOD_NATURAL) {
	SEXP perm = PROTECT(allocVector(INTSXP, (R_xlen_t) L->n));
	memcpy(INTEGER(perm), L->Perm, sizeof(int) * L->n);
	SET_SLOT(to, Matrix_permSym, perm);
	UNPROTECT(1);
	}
	SEXP colcount = PROTECT(allocVector(INTSXP, (R_xlen_t) L->n));
	memcpy(INTEGER(colcount), L->ColCount, sizeof(int) * L->n);
	SET_SLOT(to, Matrix_colcountSym, colcount);
	UNPROTECT(1);
	if (L->is_super) {
		SEXP maxcsize = PROTECT(GET_SLOT(to, Matrix_maxcsizeSym)),
			maxesize = PROTECT(GET_SLOT(to, Matrix_maxesizeSym)),
			super = PROTECT(allocVector(INTSXP, (R_xlen_t) (L->nsuper + 1))),
			pi = PROTECT(allocVector(INTSXP, (R_xlen_t) (L->nsuper + 1))),
			px = PROTECT(allocVector(INTSXP, (R_xlen_t) (L->nsuper + 1))),
			s = PROTECT(allocVector(INTSXP, (R_xlen_t) L->ssize));
		INTEGER(maxcsize)[0] = (int) L->maxcsize;
		INTEGER(maxesize)[0] = (int) L->maxesize;
		memcpy(INTEGER(super), L->super, sizeof(int) * (L->nsuper + 1));
		memcpy(INTEGER(pi), L->pi, sizeof(int) * (L->nsuper + 1));
		memcpy(INTEGER(px), L->px, sizeof(int) * (L->nsuper + 1));
		memcpy(INTEGER(s), L->s, sizeof(int) * L->ssize);
		SET_SLOT(to, Matrix_superSym, super);
		SET_SLOT(to, Matrix_piSym, pi);
		SET_SLOT(to, Matrix_pxSym, px);
		SET_SLOT(to, Matrix_sSym, s);
		UNPROTECT(6);
	} else {
		if (L->xtype != CHOLMOD_PATTERN) {
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) (L->n + 1))),
			i = PROTECT(allocVector(INTSXP, (R_xlen_t) L->nzmax)),
			nz = PROTECT(allocVector(INTSXP, (R_xlen_t) L->n)),
			next = PROTECT(allocVector(INTSXP, (R_xlen_t) (L->n + 2))),
			prev = PROTECT(allocVector(INTSXP, (R_xlen_t) (L->n + 2))),
			is_ll = PROTECT(GET_SLOT(to, Matrix_isllSym)),
			is_monotonic = PROTECT(GET_SLOT(to, Matrix_ismtSym));
		memcpy(INTEGER(p), L->p, sizeof(int) * (L->n + 1));
		memcpy(INTEGER(i), L->i, sizeof(int) * L->nzmax);
		memcpy(INTEGER(nz), L->nz, sizeof(int) * L->n);
		memcpy(INTEGER(next), L->next, sizeof(int) * (L->n + 2));
		memcpy(INTEGER(prev), L->prev, sizeof(int) * (L->n + 2));
		LOGICAL(is_ll)[0] = L->is_ll != 0;
		LOGICAL(is_monotonic)[0] = L->is_monotonic != 0;
		SET_SLOT(to, Matrix_pSym, p);
		SET_SLOT(to, Matrix_iSym, i);
		SET_SLOT(to, Matrix_nzSym, nz);
		SET_SLOT(to, Matrix_nextSym, next);
		SET_SLOT(to, Matrix_prevSym, prev);
		UNPROTECT(7);
		}
	}
	if (L->xtype != CHOLMOD_PATTERN) {
	SEXP minor = GET_SLOT(to, Matrix_minorSym);
	INTEGER(minor)[0] = (int) L->minor;
	SEXP x;
	size_t nnz = (L->is_super) ? L->xsize : L->nzmax;
	if (L->xtype == CHOLMOD_REAL) {
		PROTECT(x = allocVector(REALSXP, (R_xlen_t) nnz));
		memcpy(REAL(x), L->x, sizeof(double) * nnz);
	} else {
		PROTECT(x = allocVector(CPLXSXP, (R_xlen_t) nnz));
		memcpy(COMPLEX(x), L->x, sizeof(Rcomplex) * nnz);
	}
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
	}
	MAYBE_FREE;

#undef MAYBE_FREE

	UNPROTECT(3);
	return to;
}

/**
 * Coerce from (cholmod_sparse *) to CsparseMatrix
 *
 * Allocates an S4 object inheriting from virtual class CsparseMatrix
 * and copies into the slots from members of a pointed-to cholmod_sparse
 * struct.  The specific class of the result is determined by struct
 * members xtype and stype and by arguments ttype and doLogic.
 *
 * @param A a pointer to a cholmod_sparse struct.
 * @param doFree a flag indicating if and how to free A before returning.
 *     (0) don't free, (>0) free with cholmod_free_sparse, (<0) free with
 *     R_Free.
 * @param ttype a flag indicating if the result should be a .tCMatrix.
 *     (0) not .tCMatrix, (>0) .tCMatrix with uplo="U", (<0) .tCMatrix
 *     with uplo="L".  If ttype=0, then the result is a .gCMatrix or
 *     .sCMatrix depending on stype.  (0) .gCMatrix, (>0) .sCMatrix with
 *     uplo="U", (<0) .sCMatrix with uplo="L".
 * @param doLogic a flag indicating if the result should be a l.CMatrix
 *     if xtype=CHOLMOD_REAL.
 * @param diagString a null-terminated string or NULL.  The diag slot
 *     of a .tCMatrix result is "N" if and only if diagString is NULL
 *     or diagString[0] is 'N'.
 * @param dimnames an R object specifying the Dimnames slot of the result,
 *     unused if not a list of length 2.
 *
 * @return A CsparseMatrix.
 */
/* NB: mostly parallel to CHS2M in ./cholmod-etc.c */
SEXP cholmod_sparse_as_sexp(cholmod_sparse *A, int doFree,
                            int ttype, int doLogic, const char *diagString,
                            SEXP dimnames)
{

#define MAYBE_FREE \
	do { \
		if (doFree != 0) { \
		if (doFree < 0) \
			R_Free(A_); \
		else if (A_->itype == CHOLMOD_INT) \
			cholmod_free_sparse(&A_, &c); \
		else \
			cholmod_l_free_sparse(&A_, &cl); \
		} \
	} while (0)

	cholmod_sparse *A_ = A;
	if (A->itype != CHOLMOD_INT)
		errorFree(_("wrong '%s'"), "itype");
	if (A->xtype != CHOLMOD_PATTERN &&
	    A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		errorFree(_("wrong '%s'"), "xtype");
	if (A->dtype != CHOLMOD_DOUBLE)
		errorFree(_("wrong '%s'"), "dtype");
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		errorFree(_("dimensions cannot exceed %s"), "2^31-1");
	if (!A->sorted)
		cholmod_sort(A, &c);
	if (!A->packed || A->stype != 0)
		A = cholmod_copy(A, A->stype, 2, &c);
	char class[] = "..CMatrix";
	class[0] = (A->xtype == CHOLMOD_PATTERN)
		? 'n' : ((A->xtype == CHOLMOD_REAL) ? ((doLogic) ? 'l' : 'd') : 'z');
	class[1] = (ttype != 0) ? 't' : ((A->stype != 0) ? 's' : 'g');
	int nnz = ((int *) A->p)[A->ncol];
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym)),
		p = PROTECT(allocVector(INTSXP, (R_xlen_t) (A->ncol + 1))),
		i = PROTECT(allocVector(INTSXP, (R_xlen_t) nnz));
	INTEGER(dim)[0] = (int) A->nrow;
	INTEGER(dim)[1] = (int) A->ncol;
	memcpy(INTEGER(p), A->p, sizeof(int) * (A->ncol + 1));
	memcpy(INTEGER(i), A->i, sizeof(int) * (size_t) nnz);
	SET_SLOT(to, Matrix_pSym, p);
	SET_SLOT(to, Matrix_iSym, i);
	if (A->xtype != CHOLMOD_PATTERN) {
	SEXP x;
	if (A->xtype == CHOLMOD_REAL) {
	if (doLogic) {
		PROTECT(x = allocVector(LGLSXP, nnz));
		int *px = LOGICAL(x);
		double *Ax = (double *) A->x;
		for (int k = 0; k < nnz; ++k)
			px[k] = (ISNAN(Ax[k])) ? NA_LOGICAL : (Ax[k] != 0.0);
	} else {
		PROTECT(x = allocVector(REALSXP, nnz));
		memcpy(REAL(x), A->x, sizeof(double) * (size_t) nnz);
	}
	} else {
		PROTECT(x = allocVector(CPLXSXP, nnz));
		memcpy(COMPLEX(x), A->x, sizeof(Rcomplex) * (size_t) nnz);
	}
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
	}
	if (ttype < 0 || A->stype < 0) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1);
	}
	if (ttype != 0 && diagString && diagString[0] != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1);
	}
	if (TYPEOF(dimnames) == VECSXP && XLENGTH(dimnames) == 2)
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	if (A != A_)
		cholmod_free_sparse(&A, &c);
	MAYBE_FREE;

#undef MAYBE_FREE

	UNPROTECT(4);
	return to;
}

/**
 * Coerce from (cholmod_triplet *) to TsparseMatrix
 *
 * Allocates an S4 object inheriting from virtual class TsparseMatrix
 * and copies into the slots from members of a pointed-to cholmod_triplet
 * struct.  The specific class of the result is determined by struct
 * members xtype and stype and by arguments ttype and doLogic.
 *
 * @param A a pointer to a cholmod_triplet struct.
 * @param doFree a flag indicating if and how to free A before returning.
 *     (0) don't free, (>0) free with cholmod_free_triplet, (<0) free with
 *     R_Free.
 * @param ttype a flag indicating if the result should be a .tTMatrix.
 *     (0) not .tTMatrix, (>0) .tTMatrix with uplo="U", (<0) .tTMatrix
 *     with uplo="L".  If ttype=0, then the result is a .gTMatrix or
 *     .sTMatrix depending on stype.  (0) .gTMatrix, (>0) .sTMatrix with
 *     uplo="U", (<0) .sTMatrix with uplo="L".
 * @param doLogic a flag indicating if the result should be an l.TMatrix
 *     if xtype=CHOLMOD_REAL.
 * @param diagString a null-terminated string or NULL.  The diag slot
 *     of a .tTMatrix result is "N" if and only if diagString is NULL
 *     or diagString[0] is 'N'.
 * @param dimnames an R object specifying the Dimnames slot of the result,
 *     unused if not a list of length 2.
 *
 * @return A TsparseMatrix.
 */
SEXP cholmod_triplet_as_sexp(cholmod_triplet *A, int doFree,
                             int ttype, int doLogic, const char *diagString,
                             SEXP dimnames)
{

#define MAYBE_FREE \
	do { \
		if (doFree != 0) { \
		if (doFree < 0) \
			R_Free(A); \
		else if (A->itype == CHOLMOD_INT) \
			cholmod_free_triplet(&A, &c); \
		else \
			cholmod_l_free_triplet(&A, &cl); \
		} \
	} while (0)

	if (A->itype != CHOLMOD_INT)
		errorFree(_("wrong '%s'"), "itype");
	if (A->xtype != CHOLMOD_PATTERN &&
	    A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		errorFree(_("wrong '%s'"), "xtype");
	if (A->dtype != CHOLMOD_DOUBLE)
		errorFree(_("wrong '%s'"), "dtype");
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		errorFree(_("dimensions cannot exceed %s"), "2^31-1");
	char class[] = "..TMatrix";
	class[0] = (A->xtype == CHOLMOD_PATTERN)
		? 'n' : ((A->xtype == CHOLMOD_REAL) ? ((doLogic) ? 'l' : 'd') : 'z');
	class[1] = (ttype != 0) ? 't' : ((A->stype != 0) ? 's' : 'g');
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym)),
		i = PROTECT(allocVector(INTSXP, (R_xlen_t) A->nnz)),
		j = PROTECT(allocVector(INTSXP, (R_xlen_t) A->nnz));
	INTEGER(dim)[0] = (int) A->nrow;
	INTEGER(dim)[1] = (int) A->ncol;
	if (A->stype == 0) {
		memcpy(INTEGER(i), A->i, sizeof(int) * A->nnz);
		memcpy(INTEGER(j), A->j, sizeof(int) * A->nnz);
	} else {
		int *pi = INTEGER(i), *Ai = (int *) A->i,
			*pj = INTEGER(j), *Aj = (int *) A->j;
		for (size_t k = 0; k < A->nnz; ++k) {
			if ((Ai[k] <= Aj[k]) == (A->stype > 0)) {
				pi[k] = Ai[k];
				pj[k] = Aj[k];
			} else {
				pi[k] = Aj[k];
				pj[k] = Ai[k];
			}
		}
	}
	SET_SLOT(to, Matrix_iSym, i);
	SET_SLOT(to, Matrix_jSym, j);
	if (A->xtype != CHOLMOD_PATTERN) {
	SEXP x;
	if (A->xtype == CHOLMOD_REAL) {
	if (doLogic) {
		PROTECT(x = allocVector(LGLSXP, (R_xlen_t) A->nnz));
		int *px = LOGICAL(x);
		double *Ax = (double *) A->x;
		for (size_t k = 0; k < A->nnz; ++k)
			px[k] = (ISNAN(Ax[k])) ? NA_LOGICAL : (Ax[k] != 0.0);
	} else {
		PROTECT(x = allocVector(REALSXP, (R_xlen_t) A->nnz));
		memcpy(REAL(x), A->x, sizeof(double) * A->nnz);
	}
	} else {
		PROTECT(x = allocVector(CPLXSXP, (R_xlen_t) A->nnz));
		memcpy(COMPLEX(x), A->x, sizeof(Rcomplex) * A->nnz);
	}
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
	}
	if (ttype < 0 || A->stype < 0) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1);
	}
	if (ttype != 0 && diagString && diagString[0] != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1);
	}
	if (TYPEOF(dimnames) == VECSXP && XLENGTH(dimnames) == 2)
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);

	MAYBE_FREE;

#undef MAYBE_FREE

	UNPROTECT(4);
	return to;
}

/**
 * Coerce from (cholmod_dense *) to [dz]geMatrix
 *
 * Allocates an S4 object of class [dz]geMatrix
 * and copies into the slots from members of a pointed-to cholmod_dense
 * struct.  The specific class of the result is determined by struct
 * member xtype.
 *
 * @param A a pointer to a cholmod_dense struct.
 * @param doFree a flag indicating if and how to free A before returning.
 *     (0) don't free, (>0) free with cholmod_free_dense, (<0) free with
 *     R_Free.
 *
 * @return A [dz]geMatrix.
 */
/* NB: mostly parallel to CHD2M in ./cholmod-etc.c */
SEXP cholmod_dense_as_sexp(cholmod_dense *A, int doFree)
{

#define MAYBE_FREE \
	do { \
		if (doFree != 0) { \
		if (doFree < 0) \
			R_Free(A); \
		else \
			cholmod_free_dense(&A, &c); \
		} \
	} while (0)

	if (A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		errorFree(_("wrong '%s'"), "xtype");
	if (A->dtype != CHOLMOD_DOUBLE)
		errorFree(_("wrong '%s'"), "dtype");
	if (A->d != A->nrow) /* MJ: currently no need to support this case */
		errorFree(_("leading dimension not equal to number of rows"));
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		errorFree(_("dimensions cannot exceed %s"), "2^31-1");
	if (A->nrow * A->ncol > R_XLEN_T_MAX)
		errorFree(_("attempt to allocate vector of length exceeding %s"),
		          "R_XLEN_T_MAX");
	char class[] = ".geMatrix";
	class[0] = (A->xtype == CHOLMOD_COMPLEX) ? 'z' : 'd';
	SEXP to = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(to, Matrix_DimSym));
	INTEGER(dim)[0] = (int) A->nrow;
	INTEGER(dim)[1] = (int) A->ncol;
	SEXP x;
	if (A->xtype == CHOLMOD_REAL) {
		PROTECT(x = allocVector(REALSXP, (R_xlen_t) (A->nrow * A->ncol)));
		memcpy(REAL(x), A->x, sizeof(double) * (A->nrow * A->ncol));
	} else {
		PROTECT(x = allocVector(CPLXSXP, (R_xlen_t) (A->nrow * A->ncol)));
		memcpy(COMPLEX(x), A->x, sizeof(Rcomplex) * (A->nrow * A->ncol));
	}
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
	MAYBE_FREE;

#undef MAYBE_FREE

	UNPROTECT(2);
	return to;
}

/**
 * Log determinant from Cholesky factorization
 *
 * Computes log(det(A)) given the Cholesky factorization of A as
 * P1 * A * P1' = L1 * D * L1' = L * L', L = L1 * sqrt(D).  The
 * result is computed as sum(log(diag(D))) or 2*sum(log(diag(L))),
 * depending on members is_super and is_ll of the supplied struct.
 * Note that CHOLMOD does not require diag(D) to be positive and
 * that this routine does not check (FIXME).
 *
 * @param L a pointer to a cholmod_factor struct.
 */
double cholmod_factor_ldetA(cholmod_factor *L)
{
	int i, j, p;
	double ans = 0;
	if (L->is_super) {
		int *lpi = (int *) L->pi, *lsup = (int *) L->super;
		for (i = 0; i < L->nsuper; i++) {
			int nrp1 = 1 + lpi[i + 1] - lpi[i],
				nc = lsup[i + 1] - lsup[i];
			double *x = (double *) L->x + ((int *) L->px)[i];
			for (R_xlen_t jn = 0, j = 0; j < nc; j++, jn += nrp1)
				ans += 2.0 * log(fabs(x[jn]));
		}
	} else {
		int *li = (int *) L->i, *lp = (int *) L->p;
		double *lx = (double *) L->x;
		for (j = 0; j < L->n; j++) {
			for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++)
				;
			if (li[p] != j) {
				error(_("invalid simplicial Cholesky factorization: structural zero on main diagonal in column %d"),
				      j);
				break;
			}
			ans += log(lx[p] * ((L->is_ll) ? lx[p] : 1.0));
		}
	}
	return ans;
}

/**
 * Update a Cholesky factorization
 *
 * Updates in-place the Cholesky factorization of a symmetric matrix
 * X+alpha*I with the Cholesky factorization of
 * (1) A+beta*I, where A is a symmetric matrix sharing the nonzero pattern
 * of X, or
 * (2) A*A'+beta*I, where A is a general matrix sharing the nonzero pattern
 * of Y, assuming that X = Y*Y'.
 *
 * @param L a pointer to a cholmod_factor struct, to be modified in-place.
 * @param A a pointer to a cholmod_sparse struct.
 * @param beta a multiplier, typically positive, to guarantee strict
 *     diagonal dominance.
 *
 * @return L.
 */
cholmod_factor *cholmod_factor_update(cholmod_factor *L, cholmod_sparse *A,
                                      double beta)
{
	int ll = L->is_ll;
	double z[2];
	z[0] = beta;
	z[1] = 0.0;
	if (!cholmod_factorize_p(A, z, NULL, 0, L, &c))
		error(_("'%s' failed in '%s'"), "cholmod_factorize_p", __func__);
	if (L->is_ll != ll &&
	    !cholmod_change_factor(L->xtype, ll, L->is_super, 1, 1, L, &c))
		error(_("'%s' failed in '%s'"), "cholmod_change_factor", __func__);
	return L;
}
