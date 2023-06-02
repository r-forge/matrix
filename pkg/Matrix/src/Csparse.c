/** @file Csparse.c
 * The "CsparseMatrix" class from R package Matrix:
 *
 * Sparse matrices in compressed column-oriented form
 */
#include "Csparse.h"
#include "chm_common.h"
#include "cs_utils.h" /* -> ./cs.h  for cs_dmperm() */


#define _t_Csparse_validate
#include "t_Csparse_validate.c"

#define _t_Csparse_sort
#include "t_Csparse_validate.c"

// R: .validateCsparse(x, sort.if.needed = FALSE) :
SEXP Csparse_validate2(SEXP x, SEXP maybe_modify) {
    return Csparse_validate_(x, asLogical(maybe_modify));
}

// R: Matrix:::.sortCsparse(x) :
SEXP Csparse_sort (SEXP x) {
   int ok = Csparse_sort_2(x, TRUE); // modifying x directly
   if(!ok) warning(_("Csparse_sort(x): x is not a valid (apart from sorting) CsparseMatrix"));
   return x;
}

/** "Cheap" C version of  Csparse_validate() - *not* sorting : */
Rboolean isValid_Csparse(SEXP x)
{
    /* NB: we do *NOT* check a potential 'x' slot here, at all */
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	islot = GET_SLOT(x, Matrix_iSym);
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)), j,
	nrow = dims[0],
	ncol = dims[1],
	*xp = INTEGER(pslot),
	*xi = INTEGER(islot);

    if (length(pslot) != dims[1] + 1)
	return FALSE;
    if (xp[0] != 0)
	return FALSE;
    if (length(islot) < xp[ncol]) /* allow larger slots from over-allocation!*/
	return FALSE;
    for (j = 0; j < xp[ncol]; j++) {
	if (xi[j] < 0 || xi[j] >= nrow)
	    return FALSE;
    }
    for (j = 0; j < ncol; j++) {
	if (xp[j] > xp[j + 1])
	    return FALSE;
    }
    return TRUE;
}

/* MJ: no longer needed ... replacement in ./validity.c */
#if 0

SEXP Csparse_validate(SEXP x)
{
    return Csparse_validate_(x, FALSE);
}

SEXP Rsparse_validate(SEXP x)
{
    /* NB: we do *NOT* check a potential 'x' slot here, at all */
    SEXP pslot = GET_SLOT(x, Matrix_pSym),
	jslot = GET_SLOT(x, Matrix_jSym);
    Rboolean sorted, strictly;
    int i, k,
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	nrow = dims[0],
	ncol = dims[1],
	*xp = INTEGER(pslot),
	*xj = INTEGER(jslot);

    if (length(pslot) != dims[0] + 1)
	return mkString(_("slot p must have length = nrow(.) + 1"));
    if (xp[0] != 0)
	return mkString(_("first element of slot p must be zero"));
    if (length(jslot) < xp[nrow]) /* allow larger slots from over-allocation!*/
	return
	    mkString(_("last element of slot p must match length of slots j and x"));
    for (i = 0; i < length(jslot); i++) {
	if (xj[i] < 0 || xj[i] >= ncol)
	    return mkString(_("all column indices must be between 0 and ncol-1"));
    }
    sorted = TRUE; strictly = TRUE;
    for (i = 0; i < nrow; i++) {
	if (xp[i] > xp[i+1])
	    return mkString(_("slot p must be non-decreasing"));
	if(sorted)
	    for (k = xp[i] + 1; k < xp[i + 1]; k++) {
		if (xj[k] < xj[k - 1])
		    sorted = FALSE;
		else if (xj[k] == xj[k - 1])
		    strictly = FALSE;
	    }
    }
    if (!sorted)
	/* cannot easily use cholmod_sort(.) ... -> "error out" :*/
	return mkString(_("slot j is not increasing inside a column"));
    else if(!strictly) /* sorted, but not strictly */
	return mkString(_("slot j is not *strictly* increasing inside a column"));

    return ScalarLogical(1);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_dense() */
#if 0

/** @brief From a CsparseMatrix, produce a dense one.
 *
 * Directly deals with symmetric, triangular and general.
 * Called from ../R/Csparse.R's  C2dense()
 *
 * @param x a CsparseMatrix: currently all 9 of  "[dln][gst]CMatrix"
 * @param symm_or_tri integer (NA, < 0, > 0, = 0) specifying the knowledge of the caller about x:
 * 	NA  : unknown => will be determined
 *      = 0 : "generalMatrix" (not symm or tri);
 *      < 0 : "triangularMatrix"
 *      > 0 : "symmetricMatrix"
 *
 * @return a "denseMatrix"
 */
SEXP Csparse_to_dense(SEXP x, SEXP symm_or_tri)
{
    Rboolean is_sym, is_tri;
    int is_sym_or_tri = asInteger(symm_or_tri),
	ctype = 0; // <- default = "dgC"
    static const char *valid[] = { MATRIX_VALID_Csparse, ""};
    if(is_sym_or_tri == NA_INTEGER) { // find if  is(x, "symmetricMatrix") :
	ctype = R_check_class_etc(x, valid);
	is_sym = (ctype % 3 == 1);
	is_tri = (ctype % 3 == 2);
    } else {
	is_sym = is_sym_or_tri > 0;
	is_tri = is_sym_or_tri < 0;
	// => both are FALSE  iff  is_.. == 0
	if(is_sym || is_tri)
	    ctype = R_check_class_etc(x, valid);
    }
    CHM_SP chxs = AS_CHM_SP__(x);// -> chxs->stype = +- 1 <==> symmetric // allocated with alloca()
    R_CheckStack();
    Rboolean is_U_tri = is_tri && *diag_P(x) == 'U';
    if(is_U_tri) { // ==>  x := diagU2N(x), directly for chxs;  further: must free chxs
	CHM_SP eye = cholmod_speye(chxs->nrow, chxs->ncol, chxs->xtype, &c);
	double one[] = {1, 0};
	CHM_SP ans = cholmod_add(chxs, eye, one, one,
				 /* values: */ ((ctype / 3) != 2), // TRUE iff not "nMatrix"
				 TRUE, &c);
	cholmod_free_sparse(&eye, &c);
	chxs = cholmod_copy_sparse(ans, &c); // replacing alloca'd chxs with malloc'ed one, which must be freed
	cholmod_free_sparse(&ans, &c);
    }
    /* The following loses the symmetry property, since cholmod_dense has none,
     * BUT, much worse (FIXME!), it also transforms CHOLMOD_PATTERN ("n") matrices
     * to numeric (CHOLMOD_REAL) ones {and we "revert" via chm_dense_to_SEXP()}: */
    CHM_DN chxd = cholmod_sparse_to_dense(chxs, &c);
    /* FIXME: The above FAILS for prod(dim(.)) > INT_MAX
     * ----
     * TODO: use cholmod_l_* but also the 'cl' global ==> many changes in chm_common.[ch]
     * >>>>>>>>>>> TODO <<<<<<<<<<<<
     * CHM_DN chxd = cholmod_l_sparse_to_dense(chxs, &cl); */
    //                   ^^^ important when prod(dim(.)) > INT_MAX
    int chxs_xtype = chxs->xtype;
    int chxs_stype = chxs->stype;
    if(is_U_tri)
       cholmod_free_sparse(&chxs, &c);
    int Rkind = (chxs_xtype == CHOLMOD_PATTERN)? -1 : Real_kind(x);

    SEXP ans = chm_dense_to_SEXP(chxd, 1, Rkind, GET_SLOT(x, Matrix_DimNamesSym),
				 /* transp: */ FALSE);
    // -> a [dln]geMatrix
    if(is_sym) { // ==> want  [dln]syMatrix
	PROTECT(ans);
	const char cl1 = class_P(ans)[0];
	SEXP aa = PROTECT(NEW_OBJECT_OF_CLASS((cl1 == 'd') ? "dsyMatrix" :
					      ((cl1 == 'l') ? "lsyMatrix" : "nsyMatrix")));
	// No need to duplicate() as slots of ans are freshly allocated and ans will not be used
	SET_SLOT(aa, Matrix_xSym,       GET_SLOT(ans, Matrix_xSym));
	SET_SLOT(aa, Matrix_DimSym,     GET_SLOT(ans, Matrix_DimSym));
	SET_SLOT(aa, Matrix_DimNamesSym,GET_SLOT(ans, Matrix_DimNamesSym));
	SET_SLOT(aa, Matrix_uploSym, mkString((chxs_stype > 0) ? "U" : "L"));
	UNPROTECT(2);
	return aa;
    }
    else if(is_tri) { // ==> want  [dln]trMatrix
	PROTECT(ans);
	const char cl1 = class_P(ans)[0];
	SEXP aa = PROTECT(NEW_OBJECT_OF_CLASS((cl1 == 'd') ? "dtrMatrix" :
					      ((cl1 == 'l') ? "ltrMatrix" : "ntrMatrix")));
	// No need to duplicate() as slots of ans are freshly allocated and ans will not be used
	SET_SLOT(aa, Matrix_xSym,       GET_SLOT(ans, Matrix_xSym));
	SET_SLOT(aa, Matrix_DimSym,     GET_SLOT(ans, Matrix_DimSym));
	SET_SLOT(aa, Matrix_DimNamesSym,GET_SLOT(ans, Matrix_DimNamesSym));
	slot_dup(aa, x, Matrix_uploSym);
	/* already by NEW_OBJECT(..) above:
	   SET_SLOT(aa, Matrix_diagSym, mkString("N")); */
	UNPROTECT(2);
	return aa;
    }
    else
	return ans;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_kind() */
#if 0

// FIXME: do not go via CHM (should not be too hard, to just *drop* the x-slot, right?
SEXP Csparse2nz(SEXP x, Rboolean tri)
{
    CHM_SP chxs = AS_CHM_SP__(x);
    CHM_SP chxcp = cholmod_copy(chxs, chxs->stype, CHOLMOD_PATTERN, &c);
    R_CheckStack();

    return chm_sparse_to_SEXP(chxcp, 1/*do_free*/,
			      tri ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			      /* Rkind: pattern */ 0,
			      /* diag = */ tri ? diag_P(x) : "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_nz_pattern(SEXP x, SEXP tri)
{
    int tr_ = asLogical(tri);
    if(tr_ == NA_LOGICAL) {
	warning(_("Csparse_to_nz_pattern(x, tri = NA): 'tri' is taken as TRUE"));
	tr_ = TRUE;
    }
    return Csparse2nz(x, (Rboolean) tr_);
}

// n.CMatrix --> [dli].CMatrix  (not going through CHM!)
SEXP nz_pattern_to_Csparse(SEXP x, SEXP res_kind)
{
    return nz2Csparse(x, asInteger(res_kind));
}

// n.CMatrix --> [dli].CMatrix  (not going through CHM!)
// NOTE: use chm_MOD_xtype(() to change type of  'cholmod_sparse' matrix
SEXP nz2Csparse(SEXP x, enum x_slot_kind r_kind)
{
    const char *cl_x = class_P(x);
    // quick check - if ok, fast
    if(cl_x[0] != 'n' || cl_x[2] != 'C') {
	// e.g. class = "A", from  setClass("A", contains = "ngCMatrix")
	static const char *valid[] = { MATRIX_VALID_nCsparse, ""};
	int ctype = R_check_class_etc(x, valid);
	if(ctype < 0)
	    error(_("not a 'n.CMatrix'"));
	else // fine : get a valid  cl_x  class_P()-like string :
	    cl_x = valid[ctype];
    }
    int nnz = LENGTH(GET_SLOT(x, Matrix_iSym));
    SEXP ans;
    char *ncl = alloca(strlen(cl_x) + 1); /* not much memory required */
    strcpy(ncl, cl_x);
    double *dx_x; int *ix_x;
    ncl[0] = (r_kind == x_double ? 'd' :
	      (r_kind == x_logical ? 'l' :
	       /* else (for now):  r_kind == x_integer : */ 'i'));
    PROTECT(ans = NEW_OBJECT_OF_CLASS(ncl));
    // create a correct 'x' slot:
    switch(r_kind) {
	int i;
    case x_double: // 'd'
	dx_x = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz));
	for (i=0; i < nnz; i++) dx_x[i] = 1.;
	break;
    case x_logical: // 'l'
	ix_x = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, nnz));
	for (i=0; i < nnz; i++) ix_x[i] = TRUE;
	break;
    case x_integer: // 'i'
	ix_x = INTEGER(ALLOC_SLOT(ans, Matrix_xSym, INTSXP, nnz));
	for (i=0; i < nnz; i++) ix_x[i] = 1;
	break;

    default:
	error(_("nz2Csparse(): invalid/non-implemented r_kind = %d"),
	      r_kind);
    }

    // now copy all other slots :
    slot_dup(ans, x, Matrix_iSym);
    slot_dup(ans, x, Matrix_pSym);
    slot_dup(ans, x, Matrix_DimSym);
    slot_dup(ans, x, Matrix_DimNamesSym);
    if(ncl[1] != 'g') { // symmetric or triangular ...
	slot_dup_if_has(ans, x, Matrix_uploSym);
	slot_dup_if_has(ans, x, Matrix_diagSym);
    }
    UNPROTECT(1);
    return ans;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_matrix() */
#if 0

SEXP Csparse_to_matrix(SEXP x, SEXP chk, SEXP symm)
{
    int is_sym = asLogical(symm);
    if(is_sym == NA_LOGICAL) { // find if  is(x, "symmetricMatrix") :
	static const char *valid[] = { MATRIX_VALID_Csparse, ""};
	int ctype = R_check_class_etc(x, valid);
	is_sym = (ctype % 3 == 1);
    }
    return chm_dense_to_matrix(
	cholmod_sparse_to_dense(AS_CHM_SP2(x, asLogical(chk)), &c),
	1 /*do_free*/,
	(is_sym
	 ? get_symmetrized_DimNames(x, -1)
	 : GET_SLOT(x, Matrix_DimNamesSym)));
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_vector() */
#if 0

SEXP Csparse_to_vector(SEXP x)
{
    return chm_dense_to_vector(cholmod_sparse_to_dense(AS_CHM_SP__(x), &c), 1);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer CRsparse_as_Tsparse() */
#if 0

SEXP Csparse_to_Tsparse(SEXP x, SEXP tri)
{
    CHM_SP chxs = AS_CHM_SP__(x);
    CHM_TR chxt = cholmod_sparse_to_triplet(chxs, &c);
    int tr = asLogical(tri);
    int Rkind = (chxs->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    return chm_triplet_to_SEXP(chxt, 1,
			       tr ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			       Rkind, tr ? diag_P(x) : "",
			       GET_SLOT(x, Matrix_DimNamesSym));
}

#endif /* MJ */

/* MJ: unused */
#if 0

SEXP Csparse_to_tCsparse(SEXP x, SEXP uplo, SEXP diag)
{
    CHM_SP chxs = AS_CHM_SP__(x);
    int Rkind = (chxs->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();
    return chm_sparse_to_SEXP(chxs, /* dofree = */ 0,
			      /* uploT = */ (*CHAR(asChar(uplo)) == 'U')? 1: -1,
			       Rkind, /* diag = */ CHAR(STRING_ELT(diag, 0)),
			       GET_SLOT(x, Matrix_DimNamesSym));
}

SEXP Csparse_to_tTsparse(SEXP x, SEXP uplo, SEXP diag)
{
    CHM_SP chxs = AS_CHM_SP__(x);
    CHM_TR chxt = cholmod_sparse_to_triplet(chxs, &c);
    int Rkind = (chxs->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();
    return chm_triplet_to_SEXP(chxt, 1,
			      /* uploT = */ (*CHAR(asChar(uplo)) == 'U')? 1: -1,
			       Rkind, /* diag = */ CHAR(STRING_ELT(diag, 0)),
			       GET_SLOT(x, Matrix_DimNamesSym));
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_as_general() */
#if 0

SEXP Csparse_symmetric_to_general(SEXP x)
{
    CHM_SP chx = AS_CHM_SP__(x), chgx;
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if (!(chx->stype))
	error(_("Nonsymmetric matrix in Csparse_symmetric_to_general"));
    chgx = cholmod_copy(chx, /* stype: */ 0, chx->xtype, &c);
    return chm_sparse_to_SEXP(chgx, 1, 0, Rkind, "",
			      get_symmetrized_DimNames(x, -1));
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_force_symmetric() */
#if 0

// Called from R's forceSymmetricCsparse(), .gC2sC()
SEXP Csparse_general_to_symmetric(SEXP x, SEXP uplo, SEXP sym_dmns)
{
    int *adims = INTEGER(GET_SLOT(x, Matrix_DimSym)), n = adims[0];
    if(n != adims[1]) {
	error(_("Csparse_general_to_symmetric(): matrix is not square!"));
	return R_NilValue; /* -Wall */
    }
    CHM_SP chx = AS_CHM_SP__(x), chgx;
    int uploT = (*CHAR(asChar(uplo)) == 'U');
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();
    chgx = cholmod_copy(chx, /* stype: */ (uploT ? 1 : -1), chx->xtype, &c);

    SEXP ans, dn = GET_SLOT(x, Matrix_DimNamesSym);
    int symDmns = asLogical(sym_dmns);
    /* 3 cases ... 0: do nothing
                   1: symmetrize
          NA_LOGICAL: symmetrize if necessary
    */
    symDmns = (symDmns == 1 ||
	       (symDmns != 0 &&
		((!isNull(VECTOR_ELT(dn, 0)) && !isNull(VECTOR_ELT(dn, 1))) ||
		 !isNull(getAttrib(dn, R_NamesSymbol)))));
    if (symDmns) {
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symmDN(newdn, dn, -1);
	dn = newdn;
    }
    ans = chm_sparse_to_SEXP(chgx, 1, 0, Rkind, "", dn);
    UNPROTECT(symDmns);
    return ans;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_transpose() */
#if 0

SEXP Csparse_transpose(SEXP x, SEXP tri)
{
    /* TODO: lgCMatrix & igC* currently go via double prec. cholmod -
     *       since cholmod (& cs) lacks sparse 'int' matrices */
    CHM_SP chx = AS_CHM_SP__(x);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    CHM_SP chxt = cholmod_transpose(chx, chx->xtype, &c);
    SEXP dn = PROTECT(duplicate(GET_SLOT(x, Matrix_DimNamesSym))), tmp;
    int tr = asLogical(tri);
    R_CheckStack();

    tmp = VECTOR_ELT(dn, 0);	/* swap the dimnames */
    SET_VECTOR_ELT(dn, 0, VECTOR_ELT(dn, 1));
    SET_VECTOR_ELT(dn, 1, tmp);
    tmp = PROTECT(getAttrib(dn, R_NamesSymbol));
    if(!isNull(tmp)) { // swap names(dimnames(.)):
	SEXP nms_dns = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(nms_dns, 1, STRING_ELT(tmp, 0));
	SET_VECTOR_ELT(nms_dns, 0, STRING_ELT(tmp, 1));
	setAttrib(dn, R_NamesSymbol, nms_dns);
	UNPROTECT(1);
    }

    SEXP ans = chm_sparse_to_SEXP(chxt, 1, /* SWAP 'uplo' for triangular */
				  tr ? ((*uplo_P(x) == 'U') ? -1 : 1) : 0,
				  Rkind, tr ? diag_P(x) : "", dn);
    UNPROTECT(2);
    return ans;
}

#endif /* MJ */

/** @brief Csparse_drop(x, tol):  drop entries with absolute value < tol, i.e,
 *  at least all "explicit" zeros. */
SEXP Csparse_drop(SEXP x, SEXP tol)
{
    const char *cl = class_P(x);
    /* dtCMatrix, etc; [1] = the second character =?= 't' for triangular */
    int tr = (cl[1] == 't'); // FIXME - rather  R_check_class_etc(..)
    CHM_SP chx = AS_CHM_SP__(x);
    CHM_SP ans = cholmod_copy(chx, chx->stype, chx->xtype, &c);
    double dtol = asReal(tol);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if(!cholmod_drop(dtol, ans, &c))
	error(_("cholmod_drop() failed"));
   return chm_sparse_to_SEXP(ans, 1,
			      tr ? ((*uplo_P(x) == 'U') ? 1 : -1) : 0,
			      Rkind, tr ? diag_P(x) : "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

/** @brief Horizontal Concatenation -  cbind( <Csparse>,  <Csparse>)
 */
SEXP Csparse_horzcat(SEXP x, SEXP y)
{
#define CSPARSE_CAT(_KIND_)						\
    CHM_SP chx = AS_CHM_SP__(x), chy = AS_CHM_SP__(y);			\
    R_CheckStack();							\
    void* chx_x = chx->x;						\
    void* chx_z = chx->z;						\
    void* chy_x = chy->x;						\
    void* chy_z = chy->z;						\
    int Rk_x = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : x_pattern, \
	Rk_y = (chy->xtype != CHOLMOD_PATTERN) ? Real_kind(y) : x_pattern, Rkind; \
    if(Rk_x == x_pattern || Rk_y == x_pattern) { /* at least one of them is patter"n" */ \
	if(Rk_x == x_pattern && Rk_y == x_pattern) { /* fine */		\
	} else { /* only one is a patter"n"				\
		  * "Bug" in cholmod_horzcat()/vertcat():               \
		  * returns patter"n" matrix if one of them is */	\
	    Rboolean ok;						\
	    if(Rk_x == x_pattern) {					\
		ok = chm_MOD_xtype(CHOLMOD_REAL, chx, &c); Rk_x = 0;	\
	    } else if(Rk_y == x_pattern) {				\
		ok = chm_MOD_xtype(CHOLMOD_REAL, chy, &c); Rk_y = 0;	\
	    } else							\
		error(_("Impossible Rk_x/Rk_y in Csparse_%s(), please report"), _KIND_); \
	    if(!ok)							\
		error(_("chm_MOD_xtype() was not successful in Csparse_%s(), please report"), \
		      _KIND_);						\
	}								\
    }									\
    Rkind = /* logical if both x and y are */ (Rk_x == 1 && Rk_y == 1) ? 1 : 0

    CSPARSE_CAT("horzcat");
    // TODO: currently drops dimnames - and we fix at R level;

    SEXP retval = PROTECT(
	chm_sparse_to_SEXP(cholmod_horzcat(chx, chy, 1, &c),
			   1, 0, Rkind, "", R_NilValue));
/* AS_CHM_SP(x) fills result with points to R-allocated memory but
   chm_MOD_xtype can change ->x and ->z to cholmod_alloc'ed memory.
   The former needs no freeing but the latter does.
   The first 2 arguments to cholmod_free should contain the number
   and size of things being freed, but lying about that is sort of ok. */
#define CSPARSE_CAT_CLEANUP					\
    if (chx_x != chx->x) cholmod_free(0, 0, chx->x, &c);	\
    if (chx_z != chx->z) cholmod_free(0, 0, chx->z, &c);	\
    if (chy_x != chy->x) cholmod_free(0, 0, chy->x, &c);	\
    if (chy_z != chy->z) cholmod_free(0, 0, chy->z, &c);	\
    UNPROTECT(1);

    CSPARSE_CAT_CLEANUP;
    return retval;
}

/** @brief Vertical Concatenation -  rbind( <Csparse>,  <Csparse>)
 */
SEXP Csparse_vertcat(SEXP x, SEXP y)
{
    CSPARSE_CAT("vertcat");
    // TODO: currently drops dimnames - and we fix at R level;

    SEXP retval  = PROTECT(
	chm_sparse_to_SEXP(cholmod_vertcat(chx, chy, 1, &c),
			   1, 0, Rkind, "", R_NilValue));
    CSPARSE_CAT_CLEANUP;
    return retval;
}

/* MJ: no longer needed ... prefer R_sparse_band() */
#if 0

SEXP Csparse_band(SEXP x, SEXP k1, SEXP k2)
{
    CHM_SP chx = AS_CHM_SP__(x);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    CHM_SP ans = cholmod_band(chx, asInteger(k1), asInteger(k2), chx->xtype, &c);
    R_CheckStack();

    return chm_sparse_to_SEXP(ans, 1, /* uploT = */ 0, Rkind, "",
			      GET_SLOT(x, Matrix_DimNamesSym));
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_sparse_diag_(U2N|N2U)() */
#if 0

SEXP Csparse_diagU2N(SEXP x)
{
    const char *cl = class_P(x);
    /* dtCMatrix, etc; [1] = the second character =?= 't' for triangular */
    if (cl[1] != 't' || *diag_P(x) != 'U') {
	/* "trivially fast" when not triangular (<==> no 'diag' slot),
	   or not *unit* triangular */
	return (x);
    }
    else { /* unit triangular (diag='U'): "fill the diagonal" & diag:= "N" */
	CHM_SP chx = AS_CHM_SP__(x);
	CHM_SP eye = cholmod_speye(chx->nrow, chx->ncol, chx->xtype, &c);
	double one[] = {1, 0};
	CHM_SP ans = cholmod_add(chx, eye, one, one, TRUE, TRUE, &c);
	int uploT = (*uplo_P(x) == 'U') ? 1 : -1;
	int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;

	R_CheckStack();
	cholmod_free_sparse(&eye, &c);
	return chm_sparse_to_SEXP(ans, 1, uploT, Rkind, "N",
				  GET_SLOT(x, Matrix_DimNamesSym));
    }
}

SEXP Csparse_diagN2U(SEXP x)
{
    const char *cl = class_P(x);
    /* dtCMatrix, etc; [1] = the second character =?= 't' for triangular */
    if (cl[1] != 't' || *diag_P(x) != 'N') {
	/* "trivially fast" when not triangular (<==> no 'diag' slot),
	   or already *unit* triangular */
	return (x);
    }
    else { /* triangular with diag='N'): now drop the diagonal */
	/* duplicate, since chx will be modified: */
	SEXP xx = PROTECT(duplicate(x));
	CHM_SP chx = AS_CHM_SP__(xx);
	int uploT = (*uplo_P(x) == 'U') ? 1 : -1,
	    Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
	R_CheckStack();

	chm_diagN2U(chx, uploT, /* do_realloc */ FALSE);

	SEXP ans = chm_sparse_to_SEXP(chx, /*dofree*/ 0/* or 1 ?? */,
				      uploT, Rkind, "U",
				      GET_SLOT(x, Matrix_DimNamesSym));
	UNPROTECT(1);// only now !
	return ans;
    }
}

#endif /* MJ */

/* MJ: no longer needed ... replacement in ./subscript.c */
#if 0

/**
 * Indexing aka subsetting : Compute  x[i,j], also for vectors i and j
 * Working via CHOLMOD_submatrix, see ./CHOLMOD/MatrixOps/cholmod_submatrix.c
 * @param x CsparseMatrix
 * @param i row     indices (0-origin), or NULL (R, not C)
 * @param j columns indices (0-origin), or NULL
 *
 * @return x[i,j]  still CsparseMatrix --- currently, this loses dimnames
 */
SEXP Csparse_submatrix(SEXP x, SEXP i, SEXP j)
{
    CHM_SP chx = AS_CHM_SP(x); /* << does diagU2N() when needed */
    int rsize = (isNull(i)) ? -1 : LENGTH(i),
	csize = (isNull(j)) ? -1 : LENGTH(j);
    int Rkind = (chx->xtype != CHOLMOD_PATTERN) ? Real_kind(x) : 0;
    R_CheckStack();

    if (rsize >= 0 && !isInteger(i))
	error(_("Index i must be NULL or integer"));
    if (csize >= 0 && !isInteger(j))
	error(_("Index j must be NULL or integer"));

    /* Must treat 'NA's in i[] and j[] here -- they are *not* treated by Cholmod!
     * haveNA := ...
       if(haveNA) {
	 a. i = removeNA(i); j =removeNA(j), and remember where they were
	 b. ans = CHM_SUB(.., i, j)
	 c. add NA rows and/or columns to 'ans' according to
	    place of NA's in i and/or j.
       } else {
	 ans = CHM_SUB(.....)  // == current code
       }
     */
#define CHM_SUB(_M_, _i_, _j_)					\
    cholmod_submatrix(_M_,					\
		      (rsize < 0) ? NULL : INTEGER(_i_), rsize,	\
		      (csize < 0) ? NULL : INTEGER(_j_), csize,	\
		      TRUE, TRUE, &c)
    CHM_SP ans;
    if (!chx->stype) {/* non-symmetric Matrix */
	ans = CHM_SUB(chx, i, j);
    }
    else { /* symmetric : "dsCMatrix";
	      currently, cholmod_submatrix() only accepts "generalMatrix" */
	CHM_SP tmp = cholmod_copy(chx, /* stype: */ 0, chx->xtype, &c);
	ans = CHM_SUB(tmp, i, j);
	cholmod_free_sparse(&tmp, &c);
    }

    // "FIXME": currently dropping dimnames, and adding them afterwards in R :
    /* // dimnames: */
    /* SEXP x_dns = GET_SLOT(x, Matrix_DimNamesSym), */
    /* 	dn = PROTECT(allocVector(VECSXP, 2)); */
    return chm_sparse_to_SEXP(ans, 1, 0, Rkind, "", /* dimnames: */ R_NilValue);
}
#undef CHM_SUB

#endif /* MJ */

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
    FILE *f = fopen(CHAR(asChar(fname)), "w");

    if (!f)
	error(_("failure to open file \"%s\" for writing"),
	      CHAR(asChar(fname)));
    if (!cholmod_write_sparse(f, AS_CHM_SP(x),
			      (CHM_SP)NULL, (char*) NULL, &c))
	error(_("cholmod_write_sparse returned error code"));
    fclose(f);
    return R_NilValue;
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

/* MJ: unused */
#if 0

/**
 * Create a Csparse matrix object from indices and/or pointers.
 *
 * @param cls name of actual class of object to create
 * @param i optional integer vector of length nnz of row indices
 * @param j optional integer vector of length nnz of column indices
 * @param p optional integer vector of length np of row or column pointers
 * @param np length of integer vector p.  Must be zero if p == (int*)NULL
 * @param x optional vector of values
 * @param nnz length of vectors i, j and/or x, whichever is to be used
 * @param dims optional integer vector of length 2 to be used as
 *     dimensions.  If dims == (int*)NULL then the maximum row and column
 *     index are used as the dimensions.
 * @param dimnames optional list of length 2 to be used as dimnames
 * @param index1 indicator of 1-based indices
 *
 * @return an SEXP of class cls inheriting from CsparseMatrix.
 */
SEXP create_Csparse(char* cls, int* i, int* j, int* p, int np,
		    void* x, int nnz, int* dims, SEXP dimnames,
		    int index1)
{
    SEXP ans;
    int *ij = (int*)NULL, *tri, *trj, nrow = -1, ncol = -1;
    int xtype = -1;		/* -Wall */
    CHM_TR T;
    CHM_SP A;

    if (np < 0 || nnz < 0)
	error(_("negative vector lengths not allowed: np = %d, nnz = %d"),
	      np, nnz);
    int mi = (i == (int*)NULL), // := missing 'i'
	mj = (j == (int*)NULL), // := missing 'j'
	mp = (p == (int*)NULL); // := missing 'p'
    if ((mi + mj + mp) != 1)
	error(_("exactly 1 of 'i', 'j' or 'p' must be NULL"));
    if (mp) {
	if (np) error(_("np = %d, must be zero when p is NULL"), np);
    } else {
	if (np) {		/* Expand p to form i or j */
	    if (!(p[0])) error(_("p[0] = %d, should be zero"), p[0]);
	    for (int ii = 0; ii < np; ii++)
		if (p[ii] > p[ii + 1])
		    error(_("p must be non-decreasing"));
	    if (p[np] != nnz)
		error("p[np] = %d != nnz = %d", p[np], nnz);
	    ij = R_Calloc(nnz, int);
	    if (mi) {
		i = ij;
		nrow = np;
	    } else {
		j = ij;
		ncol = np;
	    }
	    /* Expand p to 0-based indices */
	    for (int ii = 0; ii < np; ii++)
		for (int jj = p[ii]; jj < p[ii + 1]; jj++) ij[jj] = ii;
	} else {
	    if (nnz)
		error(_("Inconsistent dimensions: np = 0 and nnz = %d"),
		      nnz);
	}
    }
    /* calculate nrow and ncol */
    if (nrow < 0) {
	for (int ii = 0; ii < nnz; ii++) {
	    int i1 = i[ii] + (index1 ? 0 : 1); /* 1-based index */
	    if (i1 < 1) error(_("invalid row index at position %d"), ii);
	    if (i1 > nrow) nrow = i1;
	}
    }
    if (ncol < 0) {
	for (int jj = 0; jj < nnz; jj++) {
	    int j1 = j[jj] + (index1 ? 0 : 1);
	    if (j1 < 1) error(_("invalid column index at position %d"), jj);
	    if (j1 > ncol) ncol = j1;
	}
    }
    if (dims != (int*)NULL) {
	if (dims[0] > nrow) nrow = dims[0];
	if (dims[1] > ncol) ncol = dims[1];
    }
    /* check the class name */
    if (strlen(cls) != 8)
	error(_("strlen of cls argument = %d, should be 8"), strlen(cls));
    if (strcmp(cls + 2, "CMatrix"))
	error(_("cls = \"%s\" does not end in \"CMatrix\""), cls);
    switch(cls[0]) {
    case 'd':
    case 'l':
	xtype = CHOLMOD_REAL;
    break;
    case 'n':
	xtype = CHOLMOD_PATTERN;
	break;
    default:
	error(_("cls = \"%s\" must begin with 'd', 'l' or 'n'"), cls);
    }
    if (cls[1] != 'g')
	error(_("Only 'g'eneral sparse matrix types allowed"));
    /* allocate and populate the triplet */
    T = cholmod_allocate_triplet((size_t)nrow, (size_t)ncol, (size_t)nnz, 0,
				 xtype, &c);
    T->x = x;
    tri = (int*)T->i;
    trj = (int*)T->j;
    for (int ii = 0; ii < nnz; ii++) {
	tri[ii] = i[ii] - ((!mi && index1) ? 1 : 0);
	trj[ii] = j[ii] - ((!mj && index1) ? 1 : 0);
    }
    /* create the cholmod_sparse structure */
    A = cholmod_triplet_to_sparse(T, nnz, &c);
    cholmod_free_triplet(&T, &c);
    /* copy the information to the SEXP */
    ans = PROTECT(NEW_OBJECT_OF_CLASS(cls));
// FIXME: This has been copied from chm_sparse_to_SEXP in  chm_common.c
    /* allocate and copy common slots */
    nnz = (int) cholmod_nnz(A, &c);
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = (int) A->nrow;
    dims[1] = (int) A->ncol;
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, A->ncol + 1)), (int*)A->p, A->ncol + 1);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nnz)), (int*)A->i, nnz);
    switch(cls[0]) {
    case 'd':
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)), (double*)A->x, nnz);
	break;
    case 'l':
	error(_("code not yet written for cls = \"lgCMatrix\""));
    }
/* FIXME: dimnames are *NOT* put there yet (if non-NULL) */
    cholmod_free_sparse(&A, &c);
    UNPROTECT(1);
    return ans;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_dense_as_sparse() */
#if 0

/**
 * Create a Csparse matrix object from a traditional R matrix
 *
 * @param x   traditional R matrix (numeric, logical, ...)
 * @param cls class (a string), currently must be  "..CMatrix"
 *
 * @return an SEXP of a class inheriting from CsparseMatrix.
 */
SEXP matrix_to_Csparse(SEXP x, SEXP cls)
{
    if (!isMatrix(x))
	error(_("%s must be (traditional R) matrix"), "'x'");
    SEXP d_x  = getAttrib(x, R_DimSymbol),
	dn_x  = getAttrib(x, R_DimNamesSymbol);
    int nr = INTEGER(d_x)[0],
	nc = INTEGER(d_x)[1];

    if (!(isString(cls) && LENGTH(cls) == 1))
	error(_("%s must be character string"), "'cls'");
    R_xlen_t ii, n = XLENGTH(x);
    int xtype = -1;
    if (n != ((R_xlen_t) nr) * nc)
	error(_("nrow * ncol = %d * %d must equal length(x) = %ld"), nr, nc, n);

    const char *ccls = CHAR(STRING_ELT(cls, 0));
    if (strlen(ccls) != 9)
	error(_("strlen of cls argument = %d, should be 9"), strlen(ccls));
    if (strcmp(ccls + 2, "CMatrix"))
	error(_("cls = \"%s\" does not end in \"CMatrix\""), ccls);
    switch(ccls[0]) {
    case 'd':
    case 'l':
	xtype = CHOLMOD_REAL;
    break;
    case 'n':
	xtype = CHOLMOD_PATTERN;
	break;
    default:
	error(_("cls = \"%s\" must begin with 'd', 'l' or 'n' for now"), ccls);
    }
    /* if (ccls[1] != 'g') */
    /* 	error(_("Only 'g'eneral sparse matrix types allowed")); */

    SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS(ccls));
    SET_SLOT(ans, Matrix_DimSym, d_x);
    SET_SLOT(ans, Matrix_DimNamesSym, (!isNull(dn_x) && LENGTH(dn_x) == 2)
	     ? duplicate(dn_x)
	     : allocVector(VECSXP, 2));

    int nz = 0, // current number of nonzero entries
	nnz = MAXOF(256, MAXOF(nr,nc));/* nnz := final number of nonzero entries, yet unknown;
					   -- must start with guess and then grow */
    int *rp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, nc + 1)),
	*ri = R_Calloc(nnz, int); // to become i slot -- of not-yet-known length nnz

    rp[0] = 0; // always

    switch(TYPEOF(x)) {
    case LGLSXP: {
	if(xtype == CHOLMOD_PATTERN) {
#         define _PATTERN_x
#         include "t_matrix_to_Csp.c"
	} else {
#         define _LOGICAL_x
#         include "t_matrix_to_Csp.c"
	}
	break;
    }
    case REALSXP: {
#       define _DOUBLE_x
#       include "t_matrix_to_Csp.c"
	break;
    }
/* case INTSXP: we would have to use
	x = coerceVector(x, REALSXP));
   and then fall through to REALSXP case, but we must *not* modify 'x' here
   FIXME: use a macro or (inline?) function with argument (y), where
   -----  SEXP y = PROTECT(coerceVector(x, REALSXP))

   ==> give error in INTSXP case, so caller (in R) must set  storage.mode(x) <- "double"
*/
#ifdef _USING_INTEGER_NOT_READY__
    case INTSXP: {
#       define _INTEGER_x
#       include "t_matrix_to_Csp.c"
	break;
    }
#endif
#ifdef _USING_COMPLEX_NOT_READY__
    case CPLXSXP: {
#       define _COMPLEX_x
#       include "t_matrix_to_Csp.c"
	break;
    }
#endif
    default:
	error(_("%s must be a logical or double vector"), "'x'");
	break;
    }

    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym,  INTSXP, nnz)), ri, nnz);
    R_Free(ri);

    UNPROTECT(1);
    return ans;
}

#endif /* MJ */

// seed will *not* be used unless it's -1 (inverse perm.) or  0 ("no" / identity) perm.
static csd* Csparse_dmperm_raw(SEXP mat, SEXP seed)
{
    mat = PROTECT(duplicate(mat));
    CSP matx = AS_CSP__(mat); /* m x n ; compressed column, *double* 'x' or none */
    int iseed = asInteger(seed);
    R_CheckStack();
    UNPROTECT(1);
    return cs_dmperm(matx, iseed); // -> ./cs.c
}

/* NB:  cs.h  defines the 'csd' struct as  (NB: csi :== int in  Matrix, for now)

   typedef struct cs_dmperm_results    // cs_dmperm or cs_scc output
   {
   csi *p ;        // size m, row permutation
   csi *q ;        // size n, column permutation
   csi *r ;        // size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
   csi *s ;        // size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
   csi nb ;        // # of blocks in fine dmperm decomposition
   csi rr [5] ;    // coarse row decomposition
   csi cc [5] ;    // coarse column decomposition
   } csd ;
*/

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
    UNPROTECT(2);
    return ans;
}
