#include "Mdefines.h"
#include "Minlines.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Csparse.h"

/* .validateCsparse(x, sort.if.needed = TRUE) */
SEXP CsparseMatrix_validate_maybe_sorting(SEXP x)
{

#define RMKMS(_FORMAT_, ...) \
	return mkString(Matrix_sprintf(_FORMAT_, __VA_ARGS__))

	/* NB: below mostly copied from CsparseMatrix_validate in ./validity.c */

	SEXP dim = GET_SLOT(x, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
		i = PROTECT(GET_SLOT(x, Matrix_iSym));
	UNPROTECT(2); /* i, p */

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[2]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "p");
	int j;
	for (j = 1; j <= n; ++j) {
		if (pp[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[j] < pp[j - 1])
			RMKMS(_("'%s' slot is not nondecreasing"), "p");
		if (pp[j] - pp[j - 1] > m)
			RMKMS(_("first differences of '%s' slot exceed %s"), "p", "Dim[1]");
	}

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (XLENGTH(i) < pp[n])
		RMKMS(_("'%s' slot has length less than %s"), "i", "p[length(p)]");
	int *pi = INTEGER(i), k, kend, ik, i0, sorted = 1;
	for (j = 1, k = 0; j <= n; ++j) {
		kend = pp[j];
		i0 = -1;
		while (k < kend) {
			ik = pi[k];
			if (ik == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (ik < 0 || ik >= m)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (ik < i0)
				sorted = 0;
			else if (ik == i0)
				RMKMS(_("'%s' slot is not increasing within columns after sorting"),
				      "i");
			i0 = ik;
			++k;
		}
	}
	if (!sorted) {
	cholmod_sparse *A = M2CS(x, 1);
	A->sorted = 0;
	if (!cholmod_sort(A, &c)) /* sorting in place !! */
		error(_("'%s' failed"), "cholmod_sort");
	pp = (int *) A->p;
	pi = (int *) A->i;
	for (j = 1, k = 0; j <= n; ++j) {
		kend = pp[j];
		i0 = -1;
		while (k < kend) {
			ik = pi[k];
			if (ik <= i0)
				RMKMS(_("'%s' slot is not increasing within columns after sorting"),
				      "i");
			i0 = ik;
			++k;
		}
	}
	}

	return ScalarLogical(1);
}

static
int strmatch(const char *x, const char **valid)
{
	int i = 0;
	while (valid[i][0] != '\0') {
		if (strcmp(x, valid[i]) == 0)
			return i;
		++i;
	}
	return -1;
}

/* <op>(diag(obj))  where  obj=dCHMsimpl (LDLt)  or  obj=dtCMatrix (nonunit) */
SEXP tCsparse_diag(SEXP obj, SEXP op)
{
	static const char *valid[] = {
		/* 0 */    "trace",
		/* 1 */   "sumLog",
		/* 2 */     "prod",
		/* 3 */      "min",
		/* 4 */      "max",
		/* 5 */    "range",
		/* 6 */     "diag",
		/* 7 */ "diagBack", ""};
	int ivalid = -1;
	if (TYPEOF(op) != STRSXP || LENGTH(op) < 1 ||
	    (op = STRING_ELT(op, 0)) == NA_STRING ||
	    (ivalid = strmatch(CHAR(op), valid)) < 0)
		error(_("invalid '%s' to '%s'"), "op", __func__);

	char ul = 'L';
	if (HAS_SLOT(obj, Matrix_uploSym)) {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = *CHAR(STRING_ELT(uplo, 0));
	}

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p) + 1, j, k = 0, kend, n = (int) (XLENGTH(p) - 1),
		len = (ivalid < 5) ? 1 : ((ivalid < 6) ? 2 : n);

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	double *px = REAL(x), tmp;

	SEXP ans = PROTECT(allocVector(REALSXP, len));
	double *pans = REAL(ans);

	switch (ivalid) {
	case 0: /*    trace */
		pans[0] = 0.0;
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend)
				pans[0] += px[(ul == 'U') ? kend - 1 : k];
			k = kend;
		}
		break;
	case 1: /*   sumLog */
		pans[0] = 0.0;
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend)
				pans[0] += log(px[(ul == 'U') ? kend - 1 : k]);
			k = kend;
		}
		break;
	case 2: /*     prod */
		pans[0] = 1.0;
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend)
				pans[0] *= px[(ul == 'U') ? kend - 1 : k];
			else
				pans[0] *= 0.0;
			k = kend;
		}
		break;
	case 3: /*      min */
		pans[0] = R_PosInf;
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			tmp = (k < kend) ? px[(ul == 'U') ? kend - 1 : k] : 0.0;
			if (ISNAN(tmp)) {
				pans[0] = tmp;
				break;
			}
			else if (tmp < pans[0])
				pans[0] = tmp;
			k = kend;
		}
		break;
	case 4: /*      max */
		pans[0] = R_NegInf;
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			tmp = (k < kend) ? px[(ul == 'U') ? kend - 1 : k] : 0.0;
			if (ISNAN(tmp)) {
				pans[0] = tmp;
				break;
			}
			else if (tmp > pans[0])
				pans[0] = tmp;
			k = kend;
		}
		break;
	case 5: /*    range */
		pans[0] = R_PosInf;
		pans[1] = R_NegInf;
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			tmp = (k < kend) ? px[(ul == 'U') ? kend - 1 : k] : 0.0;
			if (ISNAN(tmp)) {
				pans[0] = pans[1] = tmp;
				break;
			}
			else if (tmp < pans[0])
				pans[0] = tmp;
			else if (tmp > pans[1])
				pans[1] = tmp;
			k = kend;
		}
		break;
	case 6: /*     diag */
	case 7: /* diagBack */
	{
		int *pperm = NULL;
		if (ivalid == 7 && HAS_SLOT(obj, Matrix_permSym)) {
			SEXP perm = GET_SLOT(obj, Matrix_permSym);
			if (LENGTH(perm) == n)
				pperm = INTEGER(perm); /* 0-based by assumption */
		}
		for (j = 0; j < n; ++j) {
			kend = pp[j];
			pans[(pperm) ? pperm[j] : j] =
				(k < kend) ? px[(ul == 'U') ? kend - 1 : k] : 0.0;
			k = kend;
		}
		break;
	}
	default:
		break;
	}

	UNPROTECT(3);
	return ans;
}

enum x_slot_kind {
	x_unknown = -2,  /* NA */
	x_pattern = -1,  /*  n */
	x_double  =  0,  /*  d */
	x_logical =  1,  /*  l */
	x_integer =  2,  /*  i */
	x_complex =  3}; /*  z */

/* x[i, j] <- value  where  x=<.[gt]CMatrix>  and  value=<.sparseVector> */
#define _n_Csp_
#include "t_Csparse_subassign.c"

#define _l_Csp_
#include "t_Csparse_subassign.c"

#define _i_Csp_
#include "t_Csparse_subassign.c"

#define _d_Csp_
#include "t_Csparse_subassign.c"

#define _z_Csp_
#include "t_Csparse_subassign.c"

/* dmperm(x, nAns, seed) */
SEXP Csparse_dmperm(SEXP x, SEXP nans, SEXP seed)
{
	Matrix_cs *A = dgC2cs(x, 0);
	SET_MCS_XTYPE(A->xtype);
	Matrix_csd *D = Matrix_cs_dmperm(A, asInteger(seed));
	if (!D)
		return R_NilValue; /* MJ: why not an error ... ? */
	int len = asInteger(nans);
	if (len < 0)
		len = 0;
	else if (len > 6)
		len = 6;
	SEXP nms = PROTECT(allocVector(STRSXP, len)),
		ans = PROTECT(allocVector(VECSXP, len)), tmp;
	int k = len - 1;
	switch (len) {
	case 6:
		SET_STRING_ELT(nms, k, mkChar("cc5"));
		tmp = allocVector(INTSXP, 5);
		memcpy(INTEGER(tmp), D->cc, 5 * sizeof(int));
		SET_VECTOR_ELT(ans, k, tmp);
		k--;
	case 5:
		SET_STRING_ELT(nms, k, mkChar("rr5"));
		tmp = allocVector(INTSXP, 5);
		memcpy(INTEGER(tmp), D->rr, 5 * sizeof(int));
		SET_VECTOR_ELT(ans, k, tmp);
		k--;
	case 4:
		SET_STRING_ELT(nms, k, mkChar("s"));
		tmp = allocVector(INTSXP, D->nb + 1);
		memcpy(INTEGER(tmp), D->s, (D->nb + 1) * sizeof(int));
		SET_VECTOR_ELT(ans, k, tmp);
		k--;
	case 3:
		SET_STRING_ELT(nms, k, mkChar("r"));
		tmp = allocVector(INTSXP, D->nb + 1);
		memcpy(INTEGER(tmp), D->r, (D->nb + 1) * sizeof(int));
		SET_VECTOR_ELT(ans, k, tmp);
		k--;
	case 2:
		SET_STRING_ELT(nms, k, mkChar("q"));
		tmp = allocVector(INTSXP, A->n);
		for (int j = 0, *q0 = D->q, *q1 = INTEGER(tmp); j < A->n; ++j)
			q1[j] = q0[j] + 1;
		SET_VECTOR_ELT(ans, k, tmp);
		k--;
	case 1:
		SET_STRING_ELT(nms, k, mkChar("p"));
		tmp = allocVector(INTSXP, A->m);
		for (int i = 0, *p0 = D->p, *p1 = INTEGER(tmp); i < A->m; ++i)
			p1[i] = p0[i] + 1;
		SET_VECTOR_ELT(ans, k, tmp);
		k--;
	default:
		break;
	}
	D = Matrix_cs_dfree(D);
	setAttrib(ans, R_NamesSymbol, nms);
	UNPROTECT(2);
	return ans;
}

/* writeMM(obj, file) */
SEXP Csparse_MatrixMarket(SEXP obj, SEXP path)
{
	static const char *valid[] = { VALID_CSPARSE, "" };
	int ivalid = R_check_class_etc(obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(obj, __func__);
	const char *class = valid[ivalid];

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(obj, &pid);
	if (class[0] == 'l' || class[1] == 'i') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
		REPROTECT(obj = sparse_as_kind(obj, class, 'd'), pid);
		class = valid[R_check_class_etc(obj, valid)];
	}
	if (class[1] == 't') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		REPROTECT(obj = sparse_as_general(obj, class), pid);
		class = valid[R_check_class_etc(obj, valid)];
	}

	cholmod_sparse *A = M2CS(obj, 1);
	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
	}

	const char *path_ = CHAR(asChar(path));
	FILE *f = fopen(path_, "w");
	if (!f)
		error(_("failure to open file \"%s\" for writing"), path_);
	if (!cholmod_write_sparse(f, A, (cholmod_sparse *) NULL, (char *) NULL, &c))
		error(_("'%s' failed"), "cholmod_write_sparse");
	fclose(f);
	return R_NilValue;
}
