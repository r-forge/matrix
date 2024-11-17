#include "Mdefines.h"
#include "idz.h"
#include "cholmod-etc.h"

cholmod_common c ;
cholmod_common cl;

static
void Matrix_cholmod_error_handler(int status, const char *file, int line,
                                  const char *message)
{
	cholmod_defaults(&c);
	if (status < 0)
		Rf_error(_("CHOLMOD error '%s' at file '%s', line %d"),
		         message, file, line);
	else
		Rf_warning(_("CHOLMOD warning '%s' at file '%s', line %d"),
		           message, file, line);
	return;
}

int Matrix_cholmod_start(cholmod_common *Common)
{
	int ans = cholmod_start(Common);
	if (!ans)
		Rf_error(_("'%s' failed in '%s'"), "cholmod_start", __func__);
	Common->error_handler = Matrix_cholmod_error_handler;
	return ans;
}

int Matrix_cholmod_finish(cholmod_common *Common)
{
	int ans = cholmod_finish(Common);
	if (!ans)
		Rf_error(_("'%s' failed in '%s'"), "cholmod_finish", __func__);
	return ans;
}

cholmod_factor *M2CHF(SEXP obj, int values)
{
	static const char *valid[] = {
		"nsimplicialCholesky", "nsupernodalCholesky",
		"dsimplicialCholesky", "dsupernodalCholesky",
		"zsimplicialCholesky", "zsupernodalCholesky", "" };
	const char *class = Matrix_class(obj, valid, -1, __func__);
	cholmod_factor *L = (cholmod_factor *) R_alloc(1, sizeof(cholmod_factor));
	memset(L, 0, sizeof(cholmod_factor));
	values = values && (class[0] == 'd' || class[0] == 'z');
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		perm = PROTECT(GET_SLOT(obj, Matrix_permSym)),
		colcount = PROTECT(GET_SLOT(obj, Matrix_colcountSym)),
		ordering = PROTECT(GET_SLOT(obj, Matrix_orderingSym));
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
		SEXP maxcsize = PROTECT(GET_SLOT(obj, Matrix_maxcsizeSym)),
			maxesize = PROTECT(GET_SLOT(obj, Matrix_maxesizeSym)),
			super = PROTECT(GET_SLOT(obj, Matrix_superSym)),
			pi = PROTECT(GET_SLOT(obj, Matrix_piSym)),
			px = PROTECT(GET_SLOT(obj, Matrix_pxSym)),
			s = PROTECT(GET_SLOT(obj, Matrix_sSym));
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
		if (values) {
		SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
			i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			nz = PROTECT(GET_SLOT(obj, Matrix_nzSym)),
			next = PROTECT(GET_SLOT(obj, Matrix_nextSym)),
			prev = PROTECT(GET_SLOT(obj, Matrix_prevSym)),
			is_ll = PROTECT(GET_SLOT(obj, Matrix_isllSym)),
			is_monotonic = PROTECT(GET_SLOT(obj, Matrix_ismtSym));
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
	if (values) {
	SEXP minor = GET_SLOT(obj, Matrix_minorSym);
	L->minor = (size_t) INTEGER(minor)[0];
	SEXP x = GET_SLOT(obj, Matrix_xSym);
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

cholmod_sparse *M2CHS(SEXP obj, int values)
{
	const char *class = Matrix_class(obj, valid_sparse_compressed, 6, __func__);
	cholmod_sparse *A = (cholmod_sparse *) R_alloc(1, sizeof(cholmod_sparse));
	memset(A, 0, sizeof(cholmod_sparse));
	values = values && (class[0] == 'd' || class[0] == 'z');
	int mg = (class[2] == 'C') ? 1 : 0;
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, (mg == 1) ? Matrix_iSym : Matrix_jSym));
	A->nrow = (size_t) INTEGER(dim)[(mg == 1) ? 0 : 1];
	A->ncol = (size_t) INTEGER(dim)[(mg == 1) ? 1 : 0];
	A->nzmax = (size_t) INTEGER(p)[A->ncol];
	A->p = INTEGER(p);
	A->i = INTEGER(i);
	A->stype = 0;
	A->itype = CHOLMOD_INT;
	A->xtype = CHOLMOD_PATTERN;
	A->dtype = CHOLMOD_DOUBLE;
	A->sorted = 1;
	A->packed = 1;
	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		A->stype = ((*CHAR(STRING_ELT(uplo, 0)) == 'U') == (mg == 1)) ? 1 : -1;
	}
	if (values) {
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	switch (class[0]) {
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
	}
	UNPROTECT(3);
	return A;
}

cholmod_dense *M2CHD(SEXP obj, char trans)
{
	static const char *valid[] = { "dgeMatrix", "zgeMatrix", "" };
	const char *class = Matrix_class(obj, valid, 0, __func__);
	cholmod_dense *A = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(A, 0, sizeof(cholmod_dense));
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	size_t m = (size_t) INTEGER(dim)[0], n = (size_t) INTEGER(dim)[1];
	A->nrow = ((trans == 'N') ? m : n);
	A->ncol = ((trans == 'N') ? n : m);
	A->nzmax = A->nrow * A->ncol;
	A->d = A->nrow;
	A->dtype = CHOLMOD_DOUBLE;
	switch (class[0]) {
	case 'd':
	{
		double *px = REAL(x);
		if (trans != 'N') {
			double *py = (double *) R_alloc(A->nzmax, sizeof(double));
			dtrans2(py, px, m, n, trans);
			px = py;
		}
		A->x = px;
		A->xtype = CHOLMOD_REAL;
		break;
	}
	case 'z':
	{
		Rcomplex *px = COMPLEX(x);
		if (trans != 'N') {
			Rcomplex *py = (Rcomplex *) R_alloc(A->nzmax, sizeof(Rcomplex));
			ztrans2(py, px, m, n, trans);
			px = py;
		}
		A->x = px;
		A->xtype = CHOLMOD_COMPLEX;
		break;
	}
	default:
		break;
	}
	UNPROTECT(2);
	return A;
}

SEXP CHF2M(cholmod_factor *L, int values)
{
	values = values &&
		(L->xtype == CHOLMOD_REAL || L->xtype == CHOLMOD_COMPLEX);
	if (L->itype != CHOLMOD_INT)
		return errorChar(_("wrong '%s'"), "itype");
	if (values && L->dtype != CHOLMOD_DOUBLE)
		return errorChar(_("wrong '%s'"), "dtype");
	if (L->n > INT_MAX)
		return errorChar(_("dimensions cannot exceed %s"), "2^31-1");
	if (L->super) {
		if (L->maxcsize > INT_MAX)
		return errorChar(_("'%s' would overflow type \"%s\""),
		                 "maxcsize", "integer");
	} else {
		if (L->n == INT_MAX)
		return errorChar(_("n+1 would overflow type \"%s\""),
		                 "integer");
	}
	if (L->minor < L->n) {
		if (L->is_ll)
		return errorChar(_("leading principal minor of order %d is not positive"),
		                 (int) L->minor + 1);
		else
		return errorChar(_("leading principal minor of order %d is zero"),
		                 (int) L->minor + 1);
	}
	char class[] = "...........Cholesky";
	class[0] = (!values) ? 'n' : ((L->xtype == CHOLMOD_REAL) ? 'd' : 'z');
	memcpy(class + 1, (L->is_super) ? "supernodal" : "simplicial", 10);
	SEXP obj = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		ordering = PROTECT(GET_SLOT(obj, Matrix_orderingSym));
	INTEGER(ordering)[0] = L->ordering;
	INTEGER(dim)[0] = INTEGER(dim)[1] = (int) L->n;
	if (L->ordering != CHOLMOD_NATURAL) {
	SEXP perm = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) L->n));
	memcpy(INTEGER(perm), L->Perm, sizeof(int) * L->n);
	SET_SLOT(obj, Matrix_permSym, perm);
	UNPROTECT(1);
	}
	SEXP colcount = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) L->n));
	memcpy(INTEGER(colcount), L->ColCount, sizeof(int) * L->n);
	SET_SLOT(obj, Matrix_colcountSym, colcount);
	UNPROTECT(1);
	if (L->is_super) {
		SEXP maxcsize = PROTECT(GET_SLOT(obj, Matrix_maxcsizeSym)),
			maxesize = PROTECT(GET_SLOT(obj, Matrix_maxesizeSym)),
			super = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (L->nsuper + 1))),
			pi = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (L->nsuper + 1))),
			px = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (L->nsuper + 1))),
			s = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) L->ssize));
		INTEGER(maxcsize)[0] = (int) L->maxcsize;
		INTEGER(maxesize)[0] = (int) L->maxesize;
		memcpy(INTEGER(super), L->super, sizeof(int) * (L->nsuper + 1));
		memcpy(INTEGER(pi), L->pi, sizeof(int) * (L->nsuper + 1));
		memcpy(INTEGER(px), L->px, sizeof(int) * (L->nsuper + 1));
		memcpy(INTEGER(s), L->s, sizeof(int) * L->ssize);
		SET_SLOT(obj, Matrix_superSym, super);
		SET_SLOT(obj, Matrix_piSym, pi);
		SET_SLOT(obj, Matrix_pxSym, px);
		SET_SLOT(obj, Matrix_sSym, s);
		UNPROTECT(6);
	} else {
		if (values) {
		SEXP p = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (L->n + 1))),
			i = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) L->nzmax)),
			nz = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) L->n)),
			next = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (L->n + 2))),
			prev = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (L->n + 2))),
			is_ll = PROTECT(GET_SLOT(obj, Matrix_isllSym)),
			is_monotonic = PROTECT(GET_SLOT(obj, Matrix_ismtSym));
		memcpy(INTEGER(p), L->p, sizeof(int) * (L->n + 1));
		memcpy(INTEGER(i), L->i, sizeof(int) * L->nzmax);
		memcpy(INTEGER(nz), L->nz, sizeof(int) * L->n);
		memcpy(INTEGER(next), L->next, sizeof(int) * (L->n + 2));
		memcpy(INTEGER(prev), L->prev, sizeof(int) * (L->n + 2));
		LOGICAL(is_ll)[0] = L->is_ll != 0;
		LOGICAL(is_monotonic)[0] = L->is_monotonic != 0;
		SET_SLOT(obj, Matrix_pSym, p);
		SET_SLOT(obj, Matrix_iSym, i);
		SET_SLOT(obj, Matrix_nzSym, nz);
		SET_SLOT(obj, Matrix_nextSym, next);
		SET_SLOT(obj, Matrix_prevSym, prev);
		UNPROTECT(7);
		}
	}
	if (values) {
	SEXP minor = GET_SLOT(obj, Matrix_minorSym);
	INTEGER(minor)[0] = (int) L->minor;
	SEXP x;
	size_t nnz = (L->is_super) ? L->xsize : L->nzmax;
	if (L->xtype == CHOLMOD_REAL) {
		PROTECT(x = Rf_allocVector(REALSXP, (R_xlen_t) nnz));
		memcpy(REAL(x), L->x, sizeof(double) * nnz);
	} else {
		PROTECT(x = Rf_allocVector(CPLXSXP, (R_xlen_t) nnz));
		memcpy(COMPLEX(x), L->x, sizeof(Rcomplex) * nnz);
	}
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(1);
	}
	UNPROTECT(3);
	return obj;
}

SEXP CHS2M(cholmod_sparse *A, int values, char shape)
{
	cholmod_sparse *A_ = A;
	values = values &&
		(A->xtype == CHOLMOD_REAL || A->xtype == CHOLMOD_COMPLEX);
	if (A->itype != CHOLMOD_INT)
		return errorChar(_("wrong '%s'"), "itype");
	if (values && A->dtype != CHOLMOD_DOUBLE)
		return errorChar(_("wrong '%s'"), "dtype");
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		return errorChar(_("dimensions cannot exceed %s"), "2^31-1");
	if (!A->sorted)
		cholmod_sort(A, &c);
	if (!A->packed || A->stype != 0)
		A = cholmod_copy(A, A->stype, 2, &c);
	char class[] = "..CMatrix";
	class[0] = (!values) ? 'n' : ((A->xtype == CHOLMOD_REAL) ? 'd' : 'z');
	class[1] = shape;
	int nnz = ((int *) A->p)[A->ncol];
	SEXP obj = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
		p = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (A->ncol + 1))),
		i = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nnz));
	INTEGER(dim)[0] = (int) A->nrow;
	INTEGER(dim)[1] = (int) A->ncol;
	memcpy(INTEGER(p), A->p, sizeof(int) * (A->ncol + 1));
	memcpy(INTEGER(i), A->i, sizeof(int) * (size_t) nnz);
	SET_SLOT(obj, Matrix_pSym, p);
	SET_SLOT(obj, Matrix_iSym, i);
	if (values) {
	SEXP x;
	if (A->xtype == CHOLMOD_REAL) {
		PROTECT(x = Rf_allocVector(REALSXP, nnz));
		memcpy(REAL(x), A->x, sizeof(double) * (size_t) nnz);
	} else {
		PROTECT(x = Rf_allocVector(CPLXSXP, nnz));
		memcpy(COMPLEX(x), A->x, sizeof(Rcomplex) * (size_t) nnz);
	}
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(1);
	}
	if (A != A_)
		cholmod_free_sparse(&A, &c);
	UNPROTECT(4);
	return obj;
}

SEXP CHD2M(cholmod_dense *A, char trans, char shape)
{
	if (A->xtype != CHOLMOD_REAL && A->xtype != CHOLMOD_COMPLEX)
		return errorChar(_("wrong '%s'"), "xtype");
	if (A->dtype != CHOLMOD_DOUBLE)
		return errorChar(_("wrong '%s'"), "dtype");
	if (A->d != A->nrow) /* currently no need to support this case */
		return errorChar(_("leading dimension not equal to number of rows"));
	if (A->nrow > INT_MAX || A->ncol > INT_MAX)
		return errorChar(_("dimensions cannot exceed %s"), "2^31-1");
	size_t m = A->nrow, n = A->ncol;
	if (m * n > R_XLEN_T_MAX)
		return errorChar(_("attempt to allocate vector of length exceeding %s"),
		                 "R_XLEN_T_MAX");
	char class[] = "...Matrix";
	class[0] = (A->xtype == CHOLMOD_REAL) ? 'd' : 'z';
	class[1] = shape;
	class[2] = (shape == 'g')
		? 'e' : ((shape == 's') ? 'y' : ((shape == 'p') ? 'o' : 'r'));
	SEXP obj = PROTECT(newObject(class)),
		dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	INTEGER(dim)[0] = (int) ((trans != 'N') ? n : m);
	INTEGER(dim)[1] = (int) ((trans != 'N') ? m : n);
	SEXP x;
	if (A->xtype == CHOLMOD_REAL) {
		PROTECT(x = Rf_allocVector(REALSXP, (R_xlen_t) (m * n)));
		double *px = REAL(x), *py = (double *) A->x;
		dtrans2(px, py, m, n, trans);
	} else {
		PROTECT(x = Rf_allocVector(CPLXSXP, (R_xlen_t) (m * n)));
		Rcomplex *px = COMPLEX(x), *py = (Rcomplex *) A->x;
		ztrans2(px, py, m, n, trans);
	}
	SET_SLOT(obj, Matrix_xSym, x);
	UNPROTECT(3);
	return obj;
}
