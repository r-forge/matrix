#include "Mdefines.h"
#include "Minlines.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Csparse.h"

// R: .validateCsparse(x, sort.if.needed = FALSE) :
SEXP Csparse_validate2(SEXP x, SEXP maybe_modify)
{
	/* NB: we do *NOT* check a potential 'x' slot here, at all */
	SEXP pslot = GET_SLOT(x, Matrix_pSym),
		islot = GET_SLOT(x, Matrix_iSym);
	Rboolean sorted, strictly;
	int j, k,
		*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
		nrow = dims[0],
		ncol = dims[1],
		*xp = INTEGER(pslot),
		*xi = INTEGER(islot);

	if (length(pslot) != dims[1] + 1)
		return mkString(_("slot p must have length = ncol(.) + 1"));
	if (xp[0] != 0)
		return mkString(_("first element of slot p must be zero"));
	if (length(islot) < xp[ncol]) /* allow over-allocation ! */
		mkString(_("last element of slot p must match length of slots i and x"));
	for (j = 0; j < xp[ncol]; j++)
		if (xi[j] < 0 || xi[j] >= nrow)
			return mkString(_("all row indices must be between 0 and nrow-1"));
	sorted = TRUE; strictly = TRUE;
	for (j = 0; j < ncol; j++) {
		if (xp[j] > xp[j + 1])
			return mkString(_("slot p must be non-decreasing"));
		if (sorted) /* only act if >= 2 entries in column j : */
			for (k = xp[j] + 1; k < xp[j + 1]; k++) {
				if (xi[k] < xi[k - 1])
					sorted = FALSE;
				else if (xi[k] == xi[k - 1])
					strictly = FALSE;
			}
	}
	if (sorted) {
		if (!strictly)
			return mkString(_("slot i is not *strictly* increasing inside a column"));
	}
	else if (!asLogical(maybe_modify))
		return mkString(_("row indices are not sorted within columns"));
	else {
		cholmod_sparse *A = dgC2cholmod(x, 1);
		A->sorted = 0;
		if (!cholmod_sort(A, &c)) /* sorting in place !! */
			error(_("'%s' returned error code"), "cholmod_sort");
		/* Now re-check that row indices are strictly increasing */
		/* (and not just increasing) within each column :        */
		for (j = 0; j < ncol; j++) {
			for (k = xp[j] + 1; k < xp[j + 1]; k++)
				if (xi[k] == xi[k - 1])
					return mkString(_("slot i is not *strictly* increasing inside a column (even after cholmod_l_sort)"));
		}
	}
	return ScalarLogical(1);
}

enum x_slot_kind {
	x_unknown = -2,  /* NA */
	x_pattern = -1,  /*  n */
	x_double  =  0,  /*  d */
	x_logical =  1,  /*  l */
	x_integer =  2,  /*  i */
	x_complex =  3}; /*  z */

#define _d_Csp_
#include "t_Csparse_subassign.c"

#define _l_Csp_
#include "t_Csparse_subassign.c"

#define _i_Csp_
#include "t_Csparse_subassign.c"

#define _n_Csp_
#include "t_Csparse_subassign.c"

#define _z_Csp_
#include "t_Csparse_subassign.c"

SEXP Csparse_MatrixMarket(SEXP x, SEXP fname)
{
	static const char *valid[] = { VALID_CSPARSE, "" };
	int ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	const char *class = valid[ivalid];

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(x, &pid);
	if (class[0] == 'l' || class[1] == 'i') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_kind(SEXP, const char *, char);
		REPROTECT(x = sparse_as_kind(x, class, 'd'), pid);
		class = valid[R_check_class_etc(x, valid)];
	}
	if (class[1] == 't') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		REPROTECT(x = sparse_as_general(x, class), pid);
		class = valid[R_check_class_etc(x, valid)];
	}

	cholmod_sparse *A = dgC2cholmod(x, 1);
	if (class[1] == 's') {
		SEXP uplo = GET_SLOT(x, Matrix_uploSym);
		char ul = *CHAR(STRING_ELT(uplo, 0));
		A->stype = (ul == 'U') ? 1 : -1;
	}

	const char *fname_ = CHAR(asChar(fname));
	FILE *f = fopen(fname_, "w");
	if (!f)
		error(_("failure to open file \"%s\" for writing"), fname_);
	if (!cholmod_write_sparse(f, A, (cholmod_sparse *) NULL, (char *) NULL, &c))
		error(_("'%s' returned error code"), "cholmod_write_sparse");
	fclose(f);
	return R_NilValue;
}

// seed will *not* be used unless it's -1 (inverse perm.) or  0 ("no" / identity) perm.
static csd* Csparse_dmperm_raw(SEXP mat, SEXP seed)
{
    PROTECT(mat = duplicate(mat));
    cs *matx = dgC2cs(mat, HAS_SLOT(mat, Matrix_xSym));
    int iseed = asInteger(seed);
    UNPROTECT(1);
    return cs_dmperm(matx, iseed); // -> ./cs.c
}

/* MM: should allow to return the full info above
   (Timothy Davis, p.126, explains why it's interesting ..) */

/* Here, return the full *named* list to R */
SEXP Csparse_dmperm(SEXP mat, SEXP seed, SEXP nAns) {
    csd *DMp = Csparse_dmperm_raw(mat, seed);
    if(DMp == NULL) // "failure" in cs_dmperm()
	return(R_NilValue);
    int *dims = INTEGER(GET_SLOT(mat, Matrix_DimSym)),
	m = dims[0],
	n = dims[1],
	n_ans = asInteger(nAns),
	nb = DMp->nb;

    SEXP nms = PROTECT(allocVector(STRSXP, n_ans));
    SEXP ans = PROTECT(allocVector(VECSXP, n_ans));
    R_CheckStack();
    int *ip;
    /* p : */SET_STRING_ELT(nms, 0, mkChar("p"));
             SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, m));
    ip = INTEGER(VECTOR_ELT(ans, 0));
    /* 0-based permutation:
     * Memcpy(ip , (int*)(DMp->p), m); */
    // 1-based permutation:
    for(int i=0; i < m; i++) ip[i] = DMp->p[i] + 1;

    /* q : */SET_STRING_ELT(nms, 1, mkChar("q"));
             SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, n));
    ip = INTEGER(VECTOR_ELT(ans, 1));
    /* 0-based permutation:
     * Memcpy(ip , (int*)(DMp->q), m); */
    // 1-based permutation:
    for(int i=0; i < n; i++) ip[i] = DMp->q[i] + 1;

    if(n_ans > 2) {
      /* r : */  SET_STRING_ELT(nms, 2, mkChar("r"));
		 SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, nb+1));
      Memcpy(INTEGER(VECTOR_ELT(ans, 2)), (int*)(DMp->r),   nb+1);

      /* s : */  SET_STRING_ELT(nms, 3, mkChar("s"));
		 SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, nb+1));
      Memcpy(INTEGER(VECTOR_ELT(ans, 3)), (int*)(DMp->s),   nb+1);
      if(n_ans > 4) {
	/* rr5 :*/ SET_STRING_ELT(nms, 4, mkChar("rr5"));
		   SET_VECTOR_ELT(ans, 4, allocVector(INTSXP, 5));
	Memcpy(INTEGER(VECTOR_ELT(ans, 4)), (int*)(DMp->rr),  5);

	/* cc5 :*/ SET_STRING_ELT(nms, 5, mkChar("cc5"));
		   SET_VECTOR_ELT(ans, 5, allocVector(INTSXP, 5));
	Memcpy(INTEGER(VECTOR_ELT(ans, 5)), (int*)(DMp->cc),  5);
      }
    }
    setAttrib(ans, R_NamesSymbol, nms);
    DMp = cs_dfree(DMp);
    UNPROTECT(2);
    return ans;
}

/**
 * Extract the diagonal entries from *triangular* Csparse matrix  __or__ a
 * cholmod_sparse factor (LDL = TRUE).
 *
 * @param n  dimension of the matrix.
 * @param x_p  'p' (column pointer) slot contents
 * @param x_x  'x' (non-zero entries) slot contents
 * @param perm 'perm' (= permutation vector) slot contents; only used for "diagBack"
 * @param resultKind a (SEXP) string indicating which kind of result is desired.
 *
 * @return  a SEXP, either a (double) number or a length n-vector of diagonal entries
 */
static
SEXP diag_tC_ptr(int n, int *x_p, double *x_x, Rboolean is_U, int *perm,
/*                                ^^^^^^ FIXME[Generalize] to int / ... -- via x_slot_kind ? */
		 SEXP resultKind)
{
    const char* res_ch = CHAR(STRING_ELT(resultKind,0));
    enum diag_kind { diag, diag_backpermuted, trace, prod, sum_log, min, max, range
    } res_kind = ((!strcmp(res_ch, "trace")) ? trace :
		  ((!strcmp(res_ch, "sumLog")) ? sum_log :
		   ((!strcmp(res_ch, "prod")) ? prod :
		    ((!strcmp(res_ch, "min")) ? min :
		     ((!strcmp(res_ch, "max")) ? max :
		      ((!strcmp(res_ch, "range")) ? range :
		       ((!strcmp(res_ch, "diag")) ? diag :
			((!strcmp(res_ch, "diagBack")) ? diag_backpermuted :
			 -1))))))));
    int i, n_x, i_from;
    SEXP ans = PROTECT(allocVector(REALSXP,
/*                                 ^^^^  FIXME[Generalize] */
				   (res_kind == diag ||
				    res_kind == diag_backpermuted) ? n :
				   (res_kind == range ? 2 : 1)));
    double *v = REAL(ans);
/*  ^^^^^^      ^^^^  FIXME[Generalize] */

    i_from = (is_U ? -1 : 0);

#define for_DIAG(v_ASSIGN)					\
    for(i = 0; i < n; i++) {					\
	/* looking at i-th column */				\
	n_x = x_p[i+1] - x_p[i];/* #{entries} in this column */	\
	if( is_U) i_from += n_x;                                \
	v_ASSIGN;						\
	if(!is_U) i_from += n_x;                                \
    }

    /* NOTA BENE: we assume  -- uplo = "L" i.e. lower triangular matrix
     *            for uplo = "U" (makes sense with a "dtCMatrix" !),
     *            should use  x_x[i_from + (n_x - 1)] instead of x_x[i_from],
     *            where n_x = (x_p[i+1] - x_p[i])
     */

    switch(res_kind) {
    case trace: // = sum
	v[0] = 0.;
	for_DIAG(v[0] += x_x[i_from]);
	break;

    case sum_log:
	v[0] = 0.;
	for_DIAG(v[0] += log(x_x[i_from]));
	break;

    case prod:
	v[0] = 1.;
	for_DIAG(v[0] *= x_x[i_from]);
	break;

    case min:
	v[0] = R_PosInf;
	for_DIAG(if(v[0] > x_x[i_from]) v[0] = x_x[i_from]);
	break;

    case max:
	v[0] = R_NegInf;
	for_DIAG(if(v[0] < x_x[i_from]) v[0] = x_x[i_from]);
	break;

    case range:
	v[0] = R_PosInf;
	v[1] = R_NegInf;
	for_DIAG(if(v[0] > x_x[i_from]) v[0] = x_x[i_from];
		 if(v[1] < x_x[i_from]) v[1] = x_x[i_from]);
	break;

    case diag:
	for_DIAG(v[i] = x_x[i_from]);
	break;

    case diag_backpermuted:
	for_DIAG(v[i] = x_x[i_from]);

	warning(_("%s = '%s' (back-permuted) is experimental"),
		"resultKind", "diagBack");
	/* now back_permute : */
	for(i = 0; i < n; i++) {
	    double tmp = v[i]; v[i] = v[perm[i]]; v[perm[i]] = tmp;
	    /*^^^^ FIXME[Generalize] */
	}
	break;

    default: /* -1 from above */
	error(_("diag_tC(): invalid 'resultKind'"));
	/* Wall: */ ans = R_NilValue; v = REAL(ans);
    }

    UNPROTECT(1);
    return ans;
}

/**
 * Extract the diagonal entries from *triangular* Csparse matrix  __or__ a
 * cholmod_sparse factor (LDL = TRUE).
 *
 * @param obj -- now a cholmod_sparse factor or a dtCMatrix
 * @param pslot  'p' (column pointer)   slot of Csparse matrix/factor
 * @param xslot  'x' (non-zero entries) slot of Csparse matrix/factor
 * @param perm_slot  'perm' (= permutation vector) slot of corresponding CHMfactor;
 *		     only used for "diagBack"
 * @param resultKind a (SEXP) string indicating which kind of result is desired.
 *
 * @return  a SEXP, either a (double) number or a length n-vector of diagonal entries
 */
SEXP diag_tC(SEXP obj, SEXP resultKind)
{
    SEXP
	pslot = GET_SLOT(obj, Matrix_pSym),
	xslot = GET_SLOT(obj, Matrix_xSym);
    Rboolean is_U = (R_has_slot(obj, Matrix_uploSym) &&
		     *CHAR(asChar(GET_SLOT(obj, Matrix_uploSym))) == 'U');
    int n = length(pslot) - 1, /* n = ncol(.) = nrow(.) */
	*x_p  = INTEGER(pslot), pp = -1, *perm;
    double *x_x = REAL(xslot);
/*  ^^^^^^        ^^^^ FIXME[Generalize] to INTEGER(.) / LOGICAL(.) / ... xslot !*/

    if(R_has_slot(obj, Matrix_permSym))
	perm = INTEGER(GET_SLOT(obj, Matrix_permSym));
    else perm = &pp;

    return diag_tC_ptr(n, x_p, x_x, is_U, perm, resultKind);
}
