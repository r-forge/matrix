#include "Mdefines.h"
#include "Minlines.h"
#include "chm_common.h"

/* MJ: I'd like to stop using these eventually : */
#define uplo_P(x) \
	CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0))
#define diag_P(x) \
	CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0))

SEXP get_SuiteSparse_version(void)
{
    SEXP ans = allocVector(INTSXP, 3);
    int* version = INTEGER(ans);
    SuiteSparse_version(version);
    return ans;
}

SEXP chm_common_env;
static SEXP dboundSym, grow0Sym, grow1Sym, grow2Sym, maxrankSym,
    supernodal_switchSym, supernodalSym, final_asisSym, final_superSym,
    final_llSym, final_packSym, final_monotonicSym, final_resymbolSym,
    prefer_zomplexSym, prefer_upperSym, quick_return_if_not_posdefSym,
    nmethodsSym, m0_ordSym, postorderSym;

void R_cholmod_common_envset(void) {
    SEXP rho = chm_common_env;
    defineVar(dboundSym, ScalarReal(c.dbound), rho);
    defineVar(grow0Sym, ScalarReal(c.grow0), rho);
    defineVar(grow1Sym, ScalarReal(c.grow1), rho);
    defineVar(grow2Sym, ScalarInteger(c.grow2), rho);
    defineVar(maxrankSym, ScalarInteger(c.maxrank), rho);
    defineVar(supernodal_switchSym,
	      ScalarReal(c.supernodal_switch), rho);
    defineVar(supernodalSym, ScalarInteger(c.supernodal), rho);
    defineVar(final_asisSym, ScalarLogical(c.final_asis), rho);
    defineVar(final_superSym, ScalarLogical(c.final_super), rho);
    defineVar(final_llSym, ScalarLogical(c.final_ll), rho);
    defineVar(final_packSym, ScalarLogical(c.final_pack), rho);
    defineVar(final_monotonicSym, ScalarLogical(c.final_monotonic), rho);
    defineVar(final_resymbolSym, ScalarLogical(c.final_resymbol), rho);
    defineVar(prefer_zomplexSym, ScalarLogical(c.prefer_zomplex), rho);
    defineVar(prefer_upperSym, ScalarLogical(c.prefer_upper), rho);
    defineVar(quick_return_if_not_posdefSym,
	      ScalarLogical(c.quick_return_if_not_posdef), rho);
    defineVar(nmethodsSym, ScalarInteger(c.nmethods), rho);
    defineVar(m0_ordSym, ScalarInteger(c.method[0].ordering), rho);
    defineVar(postorderSym, ScalarLogical(c.postorder), rho);
}

void R_cholmod_common_envget(void) {
    SEXP rho = chm_common_env, var;

#define SET_AS_FROM_FRAME(_V_, _KIND_, _SYM_)	\
    var = PROTECT(findVarInFrame(rho, _SYM_));	\
    _V_ = _KIND_(var);				\
    UNPROTECT(1)

    SET_AS_FROM_FRAME(c.dbound, asReal,    dboundSym);
    SET_AS_FROM_FRAME(c.grow0,  asReal,    grow0Sym);
    SET_AS_FROM_FRAME(c.grow1,  asReal,    grow1Sym);
    SET_AS_FROM_FRAME(c.grow2,  asInteger, grow2Sym);
    SET_AS_FROM_FRAME(c.maxrank,asInteger, maxrankSym);
    SET_AS_FROM_FRAME(c.supernodal_switch, asReal, supernodal_switchSym);
    SET_AS_FROM_FRAME(c.supernodal,     asLogical, supernodalSym);
    SET_AS_FROM_FRAME(c.final_asis,     asLogical, final_asisSym);
    SET_AS_FROM_FRAME(c.final_super,    asLogical, final_superSym);
    SET_AS_FROM_FRAME(c.final_ll,       asLogical, final_llSym);
    SET_AS_FROM_FRAME(c.final_pack,     asLogical, final_packSym);
    SET_AS_FROM_FRAME(c.final_monotonic,asLogical, final_monotonicSym);
    SET_AS_FROM_FRAME(c.final_resymbol, asLogical, final_resymbolSym);
    SET_AS_FROM_FRAME(c.prefer_zomplex, asLogical, prefer_zomplexSym);
    SET_AS_FROM_FRAME(c.prefer_upper,   asLogical, prefer_upperSym);
    SET_AS_FROM_FRAME(c.quick_return_if_not_posdef,
					asLogical, quick_return_if_not_posdefSym);
    SET_AS_FROM_FRAME(c.nmethods,           asInteger, nmethodsSym);
    SET_AS_FROM_FRAME(c.method[0].ordering, asInteger, m0_ordSym);
    SET_AS_FROM_FRAME(c.postorder,          asLogical, postorderSym);
}

SEXP R_cholmod_common_envini(SEXP rho) {
    if (!isEnvironment(rho))
	error(_("Argument rho must be an environment"));
    chm_common_env = rho;
    dboundSym = install("dbound");
    grow0Sym = install("grow0");
    grow1Sym = install("grow1");
    grow2Sym = install("grow2");
    maxrankSym = install("maxrank");
    supernodal_switchSym = install("supernodal_switch");
    supernodalSym = install("supernodal");
    final_asisSym = install("final_asis");
    final_superSym = install("final_super");
    final_llSym = install("final_ll");
    final_packSym = install("final_pack");
    final_monotonicSym = install("final_monotonic");
    final_resymbolSym = install("final_resymbol");
    prefer_zomplexSym = install("final_zomplex");
    prefer_upperSym = install("final_upper");
    quick_return_if_not_posdefSym = install("quick_return_if_not_posdef");
    nmethodsSym = install("nmethods");
    m0_ordSym = install("m0.ord");
    postorderSym = install("postorder");
    R_cholmod_common_envset();
    return R_NilValue;
}

/** @brief stype := "symmetry type".
 *
 *  ./CHOLMOD/Include/cholmod_core.h says about  'int stype' entry of cholmod_sparse_struct:
 *    ------------------------------
 * 0:  matrix is "unsymmetric": use both upper and lower triangular parts
 *     (the matrix may actually be symmetric in pattern and value, but
 *     both parts are explicitly stored and used).  May be square or
 *     rectangular.
 * >0: matrix is square and symmetric, use upper triangular part.
 *     Entries in the lower triangular part are ignored.
 * <0: matrix is square and symmetric, use lower triangular part.
 *     Entries in the upper triangular part are ignored.
 */
static int stype(int ctype, SEXP x)
{
    if ((ctype % 3) == 1) return (*uplo_P(x) == 'U') ? 1 : -1;
    return 0;
}

/** @brief xtype: the _kind_ of numeric (think "x slot") of Cholmod sparse matrices.
  #define CHOLMOD_PATTERN 0	 pattern only, no numerical values
  #define CHOLMOD_REAL    1	 a real matrix
  #define CHOLMOD_COMPLEX 2	 a complex matrix (ANSI C99 compatible)
  #define CHOLMOD_ZOMPLEX 3	 a complex matrix (MATLAB compatible)
*/
static int xtype(int ctype)
{
    switch(ctype / 3) {
    case 0: /* "d" */
    case 1: /* "l" */
	return CHOLMOD_REAL;
    case 2: /* "n" */
	return CHOLMOD_PATTERN;
    case 3: /* "z" */
	return CHOLMOD_COMPLEX;
    }
    return -1;
}

/* coerce a vector to REAL and copy the result to freshly R_alloc'd memory */
static void *RallocedREAL(SEXP x)
{
    SEXP rx = PROTECT(coerceVector(x, REALSXP));
    int lx = LENGTH(rx);
    /* We over-allocate the memory chunk so that it is never NULL. */
    /* The CHOLMOD code checks for a NULL pointer even in the length-0 case. */
    double *ans = Memcpy((double*) R_alloc((size_t) lx + 1, sizeof(double)),
			 REAL(rx), lx);
    UNPROTECT(1);
    return (void*)ans;
}


static void *xpt(int ctype, SEXP x)
{
    switch(ctype / 3) {
    case 0: /* "d" */
	return (void *) REAL(GET_SLOT(x, Matrix_xSym));
    case 1: /* "l" */
	return RallocedREAL(GET_SLOT(x, Matrix_xSym));
    case 2: /* "n" */
	return (void *) NULL;
    case 3: /* "z" */
	return (void *) COMPLEX(GET_SLOT(x, Matrix_xSym));
    }
    return (void *) NULL; 	/* -Wall */
}

Rboolean check_sorted_chm(CHM_SP A)
{
    int *Ai = (int*)(A->i), *Ap = (int*)(A->p);
    int j, p;

    for (j = 0; j < A->ncol; j++) {
	int p1 = Ap[j], p2 = Ap[j + 1] - 1;
	for (p = p1; p < p2; p++)
	    if (Ai[p] >= Ai[p + 1])
		return FALSE;
    }
    return TRUE;
}

/**
   Copy cholmod_sparse, to an R_alloc()ed version of it
 */
static void chm2Ralloc(CHM_SP dest, CHM_SP src)
{
    int np1, nnz;

    /* copy all the characteristics of src to dest */
    memcpy(dest, src, sizeof(cholmod_sparse));

    /* R_alloc the vector storage for dest and copy the contents from src */
    np1 = (size_t) src->ncol + 1;
    nnz = (size_t) cholmod_nnz(src, &c);
    dest->p = (void *) Memcpy((   int *) R_alloc(np1, sizeof(int)),
			     (   int *) (src->p), np1);
    dest->i = (void *) Memcpy((   int *) R_alloc(nnz, sizeof(int)),
			     (   int *) (src->i), nnz);
    if (src->xtype)
    dest->x = (void *) Memcpy((double *) R_alloc(nnz, sizeof(double)),
			     (double *) (src->x), nnz);
}

/**
   Copy cholmod_triplet to an R_alloc()ed version of it
 */
static void chTr2Ralloc(CHM_TR dest, CHM_TR src)
{
    size_t nnz;

    /* copy all the (non-pointer) characteristics of src to dest */
    memcpy(dest, src, sizeof(cholmod_triplet));

    /* R_alloc the vector storage for dest and copy the contents from src */
    nnz = (size_t) src->nnz;
    dest->i = (void *) Memcpy((   int *) R_alloc(nnz, sizeof(int)),
			      (   int *) (src->i), nnz);
    dest->j = (void *) Memcpy((   int *) R_alloc(nnz, sizeof(int)),
			      (   int *) (src->j), nnz);
    if (src->xtype)
    dest->x = (void *) Memcpy((double *) R_alloc(nnz, sizeof(double)),
			      (double *) (src->x), nnz);
}

/** "Cheap" C version of  Csparse_validate() - *not* sorting : */
static Rboolean isValid_Csparse(SEXP x)
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

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macros
 * AS_CHM_SP() or AS_CHM_SP__().  It is unusual to call it directly.
 *
 * @param ans a CHM_SP pointer
 * @param x pointer to an object that inherits from CsparseMatrix
 * @param check_Udiag boolean - should a check for (and consequent
 *  expansion of) a unit diagonal be performed.
 * @param sort_in_place boolean - if the i and x slots are to be sorted
 *  should they be sorted in place?  If the i and x slots are pointers
 *  to an input SEXP they should not be modified.
 *
 * @return ans containing pointers to the slots of x, *unless*
 *	check_Udiag and x is unitriangular.
 */
/* AS_CHM_SP  (x) := as_cholmod_sparse((CHM_SP)alloca(sizeof(cholmod_sparse)), x, TRUE,  FALSE)
 * AS_CHM_SP__(x) := as_cholmod_sparse((CHM_SP)alloca(sizeof(cholmod_sparse)), x, FALSE, FALSE)
 */
CHM_SP as_cholmod_sparse(CHM_SP ans, SEXP x,
			 Rboolean check_Udiag, Rboolean sort_in_place)
{
    static const char *valid[] = { VALID_CSPARSE, ""};
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	ctype = R_check_class_etc(x, valid);
    SEXP islot = GET_SLOT(x, Matrix_iSym);

    if (ctype < 0)
	error(_("invalid class of object to as_cholmod_sparse"));
    if (!isValid_Csparse(x))
	error(_("invalid object passed to as_cholmod_sparse"));

    memset(ans, 0, sizeof(cholmod_sparse)); /* zero the struct */

    ans->itype = CHOLMOD_INT;	/* characteristics of the system */
    ans->dtype = CHOLMOD_DOUBLE;
    ans->packed = TRUE;
				/* slots always present */
    ans->i = INTEGER(islot);
    ans->p = INTEGER(GET_SLOT(x, Matrix_pSym));
				/* dimensions and nzmax */
    ans->nrow = dims[0];
    ans->ncol = dims[1];
    /* Allow for over-allocation of the i and x slots.  Needed for
     * sparse X form in lme4.  Right now it looks too difficult to
     * check for the length of the x slot, because of the xpt
     * utility, but the lengths of x and i should agree. */
    ans->nzmax = LENGTH(islot);
				/* values depending on ctype */
    ans->x     = xpt  (ctype, x);
    ans->stype = stype(ctype, x);
    ans->xtype = xtype(ctype);

    /* are the columns sorted (increasing row numbers) ?*/
    ans->sorted = check_sorted_chm(ans);
    if (!(ans->sorted)) { /* sort columns */
	if(sort_in_place) {
	    if (!cholmod_sort(ans, &c))
		error(_("in_place cholmod_sort returned an error code"));
	    ans->sorted = 1;
	}
	else {
	    CHM_SP tmp = cholmod_copy_sparse(ans, &c);
	    if (!cholmod_sort(tmp, &c))
		error(_("cholmod_sort returned an error code"));

#ifdef DEBUG_Matrix
	    /* This "triggers" exactly for return values of dtCMatrix_sparse_solve():*/
	    /* Don't want to translate this: want it report */
	    Rprintf("Note: as_cholmod_sparse() needed cholmod_sort()ing\n");
#endif
	    chm2Ralloc(ans, tmp);
	    cholmod_free_sparse(&tmp, &c);
	}
    }

    if (check_Udiag && ctype % 3 == 2 /* triangular */ && ans->nrow // fails for Dim = (0,0)
	&& (*diag_P(x) == 'U')) { /* diagU2N(.)  "in place" : */
	double one[] = {1, 0};
	CHM_SP eye = cholmod_speye(ans->nrow, ans->ncol, ans->xtype, &c);
	CHM_SP tmp = cholmod_add(ans, eye, one, one, TRUE, TRUE, &c);

#ifdef DEBUG_Matrix_verbose /* happens quite often, e.g. in ../tests/indexing.R : */
	Rprintf("Note: as_cholmod_sparse(<ctype=%d>) - diagU2N\n", ctype);
#endif
	chm2Ralloc(ans, tmp);
	cholmod_free_sparse(&tmp, &c);
	cholmod_free_sparse(&eye, &c);
    } /* else :
       * NOTE: if(*diag_P(x) == 'U'), the diagonal is lost (!);
       * ---- that may be ok, e.g. if we are just converting from/to Tsparse,
       *      but is *not* at all ok, e.g. when used before matrix products */

    return ans;
}

/**
 * Copy the contents of a to an appropriate CsparseMatrix object and,
 * optionally, free a or free both a and its the pointers to its contents.
 *
 * @param a  (cholmod_sparse) matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 R_Free a
 * @param uploT 0 - not triangular; > 0 upper triangular; < 0 lower
 * @param Rkind - vector type to store for a->xtype == CHOLMOD_REAL,
 *                0 - REAL; 1 - LOGICAL  [unused for other a->xtype]
 * @param diag character string suitable for the diag slot of a
 *          triangular matrix (not accessed if uploT == 0).
 * @param dn either R_NilValue or an SEXP suitable for the Dimnames slot.
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_sparse_to_SEXP(CHM_SP a, int dofree, int uploT, int Rkind,
			const char* diag, SEXP dn)
{
    PROTECT(dn); /* dn is usually UNPROTECTed before the call */

				/* ensure a is sorted and packed */

    Rboolean longi = (a->itype) == CHOLMOD_LONG;
    if (!a->sorted || !a->packed)
        longi ? cholmod_l_sort(a, &cl) : cholmod_sort(a, &c);

    SEXP ans;
    char *cls = "";/* -Wall */
    int *dims, nnz, *ansp, *ansi;
    // if (longi) :
    SuiteSparse_long
	*ail = (SuiteSparse_long*)(a->i),
	*apl = (SuiteSparse_long*)(a->p);
    // else  ((a->itype) == CHOLMOD_INT) :
    int *aii = (int*)(a->i),
	*api = (int*)(a->p);

				/* determine the class of the result */

#define DOFREE_MAYBE							\
    if (dofree > 0)							\
	longi ? cholmod_l_free_sparse(&a, &cl) : cholmod_free_sparse(&a, &c); \
    else if (dofree < 0) R_Free(a)


    switch(a->xtype) {
    case CHOLMOD_PATTERN:
	cls = uploT ? "ntCMatrix": ((a->stype) ? "nsCMatrix" : "ngCMatrix");
	break;
    case CHOLMOD_REAL:
	switch(Rkind) {
	case 0:
	    cls = uploT ? "dtCMatrix": ((a->stype) ? "dsCMatrix" : "dgCMatrix");
	    break;
	case 1:
	    cls = uploT ? "ltCMatrix": ((a->stype) ? "lsCMatrix" : "lgCMatrix");
	    break;
	default:
	    DOFREE_MAYBE;
	    error(_("chm_sparse_to_SEXP(<real>, *): invalid 'Rkind' (real kind code)"));
	}
	break;
    case CHOLMOD_COMPLEX:
	cls = uploT ? "ztCMatrix": ((a->stype) ? "zsCMatrix" : "zgCMatrix");
	break;
    default:
	DOFREE_MAYBE;
	error(_("unknown xtype in cholmod_sparse object"));
    }
    ans = PROTECT(newObject(cls));
				/* allocate and copy common slots */
    nnz = longi ? cholmod_l_nnz(a, &cl) : cholmod_nnz(a, &c);
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->nrow; dims[1] = a->ncol;
    ansp = INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, a->ncol + 1));
    ansi = INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, nnz));
    for (int j = 0; j <= a->ncol; j++) ansp[j] = longi ? (int)(apl[j]) : api[j];
    for (int p = 0; p < nnz; p++)      ansi[p] = longi ? (int)(ail[p]) : aii[p];
				/* copy data slot if present */
    if (a->xtype == CHOLMOD_REAL) {
	int i, *m_x;
	double *a_x = (double *) a->x;
	switch(Rkind) {
	case 0:
	    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nnz)),
		   a_x, nnz);
	    break;
	case 1:
	    m_x = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, nnz));
	    for (i=0; i < nnz; i++)
		m_x[i] = ISNAN(a_x[i]) ? NA_LOGICAL : (a_x[i] != 0);
	    break;
	}
    }
    else if (a->xtype == CHOLMOD_COMPLEX) {
	DOFREE_MAYBE;
	error(_("complex sparse matrix code not yet written"));
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, nnz)), */
/* 	       (complex *) a->x, nnz); */
    }
    if (uploT) {		/* slots for triangularMatrix */
	if (a->stype) error(_("Symmetric and triangular both set"));
	SET_SLOT(ans, Matrix_uploSym, mkString((uploT > 0) ? "U" : "L"));
	SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    }
    if (a->stype)		/* slot for symmetricMatrix */
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((a->stype > 0) ? "U" : "L"));
    DOFREE_MAYBE;
    if (dn != R_NilValue)
	SET_SLOT(ans, Matrix_DimNamesSym, duplicate(dn));

    UNPROTECT(2);
    return ans;
}
#undef DOFREE_MAYBE

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macros
 * AS_CHM_TR() or AS_CHM_TR__().  It is unusual to call it directly.
 *
 * @param ans a CHM_TR pointer
 * @param x pointer to an object that inherits from TsparseMatrix
 * @param check_Udiag boolean - should a check for (and consequent
 *  expansion of) a unit diagonal be performed.
 *
 * @return ans containing pointers to the slots of x, *unless*
 *	check_Udiag and x is unitriangular.
 */
CHM_TR as_cholmod_triplet(CHM_TR ans, SEXP x, Rboolean check_Udiag)
{
    static const char *valid[] = { VALID_TSPARSE, ""};
    int ctype = R_check_class_etc(x, valid),
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    SEXP islot = GET_SLOT(x, Matrix_iSym);
    int m = LENGTH(islot);
    Rboolean do_Udiag = (check_Udiag && ctype % 3 == 2 && (*diag_P(x) == 'U'));
    if (ctype < 0) error(_("invalid class of object to as_cholmod_triplet"));

    memset(ans, 0, sizeof(cholmod_triplet)); /* zero the struct */

    ans->itype = CHOLMOD_INT;	/* characteristics of the system */
    ans->dtype = CHOLMOD_DOUBLE;
				/* nzmax, dimensions, types and slots : */
    ans->nnz = ans->nzmax = m;
    ans->nrow = dims[0];
    ans->ncol = dims[1];
    ans->stype = stype(ctype, x);
    ans->xtype = xtype(ctype);
    ans->i = (void *) INTEGER(islot);
    ans->j = (void *) INTEGER(GET_SLOT(x, Matrix_jSym));
    ans->x = xpt(ctype, x);

    if(do_Udiag) {
	/* diagU2N(.) "in place", similarly to Tsparse_diagU2N [./Tsparse.c]
	   (but without new SEXP): */
	int k = m + dims[0];
	CHM_TR tmp = cholmod_l_copy_triplet(ans, &cl);
	int *a_i, *a_j;

	if(!cholmod_reallocate_triplet((size_t) k, tmp, &cl))
	    error(_("as_cholmod_triplet(): could not reallocate for internal diagU2N()"
		      ));

	/* TODO? instead of copy_triplet() & reallocate_triplet()
	 * ---- allocate to correct length + Memcpy() here, as in
	 * Tsparse_diagU2N() & chTr2Ralloc() below */
	a_i = tmp->i;
	a_j = tmp->j;
	/* add (@i, @j)[k+m] = k, @x[k+m] = 1.   for k = 0,..,(n-1) */
	for(k=0; k < dims[0]; k++) {
	    a_i[k+m] = k;
	    a_j[k+m] = k;

	    switch(ctype / 3) {
	    case 0: { /* "d" */
		double *a_x = tmp->x;
		a_x[k+m] = 1.;
		break;
	    }
	    case 1: { /* "l" */
		int *a_x = tmp->x;
		a_x[k+m] = 1;
		break;
	    }
	    case 2: /* "n" */
		break;
	    case 3: { /* "z" */
		double *a_x = tmp->x;
		a_x[2*(k+m)  ] = 1.;
		a_x[2*(k+m)+1] = 0.;
		break;
	    }
	    }
	} /* for(k) */

	chTr2Ralloc(ans, tmp);
	cholmod_l_free_triplet(&tmp, &c);

    } /* else :
       * NOTE: if(*diag_P(x) == 'U'), the diagonal is lost (!);
       * ---- that may be ok, e.g. if we are just converting from/to Tsparse,
       *      but is *not* at all ok, e.g. when used before matrix products */

    return ans;
}

/**
 * Copy the contents of a to an appropriate TsparseMatrix object and,
 * optionally, free a or free both a and its the pointers to its contents.
 *
 * @param a matrix to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 R_Free a
 * @param uploT 0 - not triangular; > 0 upper triangular; < 0 lower
 * @param Rkind - vector type to store for a->xtype == CHOLMOD_REAL,
 *                0 - REAL; 1 - LOGICAL
 * @param diag character string suitable for the diag slot of a
 *          triangular matrix (not accessed if uploT == 0).
 * @param dn either R_NilValue or an SEXP suitable for the Dimnames slot.
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_triplet_to_SEXP(CHM_TR a, int dofree, int uploT, int Rkind,
			 const char* diag, SEXP dn)
{
    SEXP ans;
    char *cl = "";		/* -Wall */
    int *dims;

    PROTECT(dn);  /* dn is usually UNPROTECTed before the call */
				/* determine the class of the result */

#define DOFREE_MAYBE					\
    if (dofree > 0) cholmod_free_triplet(&a, &c);	\
    else if (dofree < 0) R_Free(a)

    switch(a->xtype) {
    case CHOLMOD_PATTERN:
	cl = uploT ? "ntTMatrix" :
	    ((a->stype) ? "nsTMatrix" : "ngTMatrix");
	break;
    case CHOLMOD_REAL:
	switch(Rkind) {
	case 0:
	    cl = uploT ? "dtTMatrix" :
		((a->stype) ? "dsTMatrix" : "dgTMatrix");
	    break;
	case 1:
	    cl = uploT ? "ltTMatrix" :
		((a->stype) ? "lsTMatrix" : "lgTMatrix");
	    break;
	}
	break;
    case CHOLMOD_COMPLEX:
	cl = uploT ? "ztTMatrix" :
	    ((a->stype) ? "zsTMatrix" : "zgTMatrix");
	break;
    default:
	DOFREE_MAYBE;
	error(_("unknown xtype in cholmod_triplet object"));
    }
    ans = PROTECT(newObject(cl));
				/* allocate and copy common slots */
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = a->nrow; dims[1] = a->ncol;
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, a->nnz)),
	   (int *) a->i, a->nnz);
    Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_jSym, INTSXP, a->nnz)),
	   (int *) a->j, a->nnz);
				/* copy data slot if present */
    if (a->xtype == CHOLMOD_REAL) {
	int i, *m_x;
	double *a_x = (double *) a->x;
	switch(Rkind) {
	case 0:
	    Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, a->nnz)),
		   a_x, a->nnz);
	    break;
	case 1:
	    m_x = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, a->nnz));
	    for (i=0; i < a->nnz; i++)
		m_x[i] = ISNAN(a_x[i]) ? NA_LOGICAL : (a_x[i] != 0);
	    break;
	}
    }
    else if (a->xtype == CHOLMOD_COMPLEX) {
	DOFREE_MAYBE;
	error(_("complex sparse matrix code not yet written"));
/* 	Memcpy(COMPLEX(ALLOC_SLOT(ans, Matrix_xSym, CPLXSXP, a->nnz)), */
/* 	       (complex *) a->x, a->nz); */
    }
    if (uploT) {		/* slots for triangularMatrix */
	if (a->stype) error(_("Symmetric and triangular both set"));
	SET_SLOT(ans, Matrix_uploSym, mkString((uploT > 0) ? "U" : "L"));
	SET_SLOT(ans, Matrix_diagSym, mkString(diag));
    }
				/* set symmetry attributes */
    if (a->stype)
	SET_SLOT(ans, Matrix_uploSym,
		 mkString((a->stype > 0) ? "U" : "L"));
    DOFREE_MAYBE;
    if (dn != R_NilValue)
	SET_SLOT(ans, Matrix_DimNamesSym, duplicate(dn));
    UNPROTECT(2);
    return ans;
}
#undef DOFREE_MAYBE

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macro AS_CHM_DN.
 * It is unusual to call it directly.
 *
 * @param ans a CHM_DN pointer.
 * @param x pointer to an object that inherits from (denseMatrix ^ generalMatrix)
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_DN as_cholmod_dense(CHM_DN ans, SEXP x)
{
#define _AS_cholmod_dense_1						\
    static const char *valid[] = { "dgeMatrix", "lgeMatrix", "ngeMatrix", ""}; \
    int dims[2], ctype = R_check_class_etc(x, valid), nprot = 0;	\
									\
    if (ctype < 0) {		/* not a classed matrix */		\
	if (isMatrix(x)) Memcpy(dims, INTEGER(getAttrib(x, R_DimSymbol)), 2); \
	else {dims[0] = LENGTH(x); dims[1] = 1;}			\
	if (isInteger(x)) {						\
	    x = PROTECT(coerceVector(x, REALSXP));			\
	    nprot++;							\
	}								\
	ctype = (isReal(x) ? 0 :					\
		 (isLogical(x) ? 2 : /* logical -> default to "l", not "n" */ \
		  (isComplex(x) ? 6 : -1)));				\
    } else Memcpy(dims, INTEGER(GET_SLOT(x, Matrix_DimSym)), 2);	\
    if (ctype < 0) error(_("invalid class of object to as_cholmod_dense")); \
    memset(ans, 0, sizeof(cholmod_dense)); /* zero the struct */        \
                                                                        \
    ans->dtype = CHOLMOD_DOUBLE; /* characteristics of the system */	\
    ans->x = ans->z = (void *) NULL;					\
				/* dimensions and nzmax */		\
    ans->d = ans->nrow = dims[0];					\
    ans->ncol = dims[1];						\
    ans->nzmax = ((size_t)dims[0]) * dims[1];				\
				/* set the xtype and any elements */	\
    switch(ctype / 2) {							\
    case 0: /* "d" */							\
	ans->xtype = CHOLMOD_REAL;					\
	ans->x = (void *) REAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x); \
	break

    _AS_cholmod_dense_1;

    case 1: /* "l" */
	ans->xtype = CHOLMOD_REAL;
	ans->x = RallocedREAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x);
	break;
    case 2: /* "n" */
	ans->xtype = CHOLMOD_PATTERN;
	ans->x = (void *) LOGICAL((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x);
	break;

#define _AS_cholmod_dense_2						\
    case 3: /* "z" */							\
	ans->xtype = CHOLMOD_COMPLEX;					\
	ans->x = (void *) COMPLEX((ctype % 2) ? GET_SLOT(x, Matrix_xSym) : x); \
	break;								\
    }									\
    UNPROTECT(nprot);							\
    return ans

    _AS_cholmod_dense_2;
}

#undef _AS_cholmod_dense_1
#undef _AS_cholmod_dense_2

void R_cholmod_error(int status, const char *file, int line, const char *message)
{
    R_cholmod_common_envget(); /* restore any setting that may have been changed */

/* NB: keep in sync with M_R_cholmod_error(), ../inst/include/Matrix_stubs.c */

    /* From CHOLMOD/Include/cholmod_core.h : ...status values.
       zero means success, negative means a fatal error, positive is a warning.
    */
#ifndef R_CHOLMOD_ALWAYS_ERROR
    if(status < 0) {
#endif
	error(_("Cholmod error '%s' at file %s, line %d"), message, file, line);
#ifndef R_CHOLMOD_ALWAYS_ERROR
    }
    else
	warning(_("Cholmod warning '%s' at file %s, line %d"),
		message, file, line);
#endif
}

/* just to get 'int' instead of 'void' as required by CHOLMOD's print_function */
static
int R_cholmod_printf(const char* fmt, ...)
{
    va_list(ap);

    va_start(ap, fmt);
    Rprintf((char *)fmt, ap);
    va_end(ap);
    return 0;
}

/**
 * Initialize the CHOLMOD library and replace the print and error functions
 * by R-specific versions.
 *
 * @param c pointer to a cholmod_common structure to be initialized
 *
 * @return TRUE if successful
 */
int R_cholmod_start(CHM_CM c)
{
    int res;
    if (!(res = cholmod_start(c)))
	error(_("Unable to initialize cholmod: error code %d"), res);
    /*SuiteSparse <= 4.x.y :
     * c->print_function = R_cholmod_printf; /. Rprintf gives warning */
    SuiteSparse_config.printf_func = R_cholmod_printf;/* Rprintf gives warning */
    //                                 ^^^^^^^^^ now is misnomer
    /* Since we provide an error handler, it may not be a good idea to allow CHOLMOD printing,
     * because that's not easily suppressed on the R level :
     * Hence consider, at least temporarily *  c->print_function = NULL; */
    c->error_handler = R_cholmod_error;
    return TRUE;
}

CHM_DN numeric_as_chm_dense(CHM_DN ans, double *v, int nr, int nc)
{
    ans->d = ans->nrow = nr;
    ans->ncol = nc;
    ans->nzmax = nr * nc;
    ans->x = (void *) v;
    ans->xtype = CHOLMOD_REAL;
    ans->dtype = CHOLMOD_DOUBLE;
    return ans;
}

/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * @param ans an CHM_FR object
 * @param x pointer to an object that inherits from CHMfactor
 * @param do_check logical indicating if check for correctness should happen
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_FR as_cholmod_factor3(CHM_FR ans, SEXP x, Rboolean do_check)
{
    static const char *valid[] = {
		"dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl", ""};
    int ctype = R_check_class_etc(x, valid);
    if (ctype < 0)
	error(_("object of invalid class to 'as_cholmod_factor()'"));
    memset(ans, 0, sizeof(cholmod_factor));

    SEXP tmp = GET_SLOT(x, install("type"));
    int *type = INTEGER(tmp);

    ans->ordering = type[0];
    ans->is_super = type[2];

    tmp = GET_SLOT(x, install("colcount"));
    ans->n = LENGTH(tmp);
    ans->minor = ans->n;
    ans->ColCount = INTEGER(tmp);

    if (ans->ordering != CHOLMOD_NATURAL)
	ans->Perm = INTEGER(GET_SLOT(x, Matrix_permSym));
    else {
	int j, n = (int) ans->n, *Perm = (int *) R_alloc(ans->n, sizeof(int));
	for (j = 0; j < n; ++j)
	    Perm[j] = j;
	ans->Perm = Perm;
    }

    ans->itype = CHOLMOD_INT;
    ans->dtype = CHOLMOD_DOUBLE;
    if (ctype >= 2)
	ans->xtype = CHOLMOD_PATTERN;
    else {
	ans->xtype = CHOLMOD_REAL;
	ans->x = REAL(GET_SLOT(x, Matrix_xSym));
    }

    if (ans->is_super) {
	tmp = GET_SLOT(x, install("super"));
	ans->nsuper = LENGTH(tmp) - 1;
	ans->super = INTEGER(tmp);
	ans->pi = INTEGER(GET_SLOT(x, install("pi")));
	ans->px = INTEGER(GET_SLOT(x, install("px")));
	ans->s = INTEGER(GET_SLOT(x, install("s")));
	ans->ssize = ((int *) ans->pi)[ans->nsuper];
	ans->xsize = ((int *) ans->px)[ans->nsuper];
	ans->maxcsize = type[4];
	ans->maxesize = type[5];
	ans->is_ll = 1;
	ans->is_monotonic = 1;
    } else {
	ans->p = INTEGER(GET_SLOT(x, Matrix_pSym));
	ans->i = INTEGER(GET_SLOT(x, Matrix_iSym));
	ans->nz = INTEGER(GET_SLOT(x, install("nz")));
	ans->next = INTEGER(GET_SLOT(x, install("nxt")));
	ans->prev = INTEGER(GET_SLOT(x, install("prv")));
	ans->nzmax = ((int *) ans->p)[ans->n];
	ans->is_ll = type[1];
	ans->is_monotonic = type[3];
    }

    if (do_check && !cholmod_check_factor(ans, &c))
	error(_("failure in as_cholmod_factor"));
    return ans;
}

// This has been in the Matrix API  ( ../inst/include/Matrix.h
/**
 * Populate ans with the pointers from x and modify its scalar
 * elements accordingly. Note that later changes to the contents of
 * ans will change the contents of the SEXP.
 *
 * In most cases this function is called through the macro AS_CHM_FR.
 * It is unusual to call it directly.
 *
 * @param ans an CHM_FR object
 * @param x pointer to an object that inherits from CHMfactor
 *
 * @return ans containing pointers to the slots of x.
 */
CHM_FR as_cholmod_factor(CHM_FR ans, SEXP x) {
    return as_cholmod_factor3(ans, x, /* do_check = */ TRUE);
}


/**
 * Copy the contents of f to an appropriate CHMfactor object and,
 * optionally, free f or free both f and its pointer to its contents.
 *
 * @param f cholmod_factor object to be converted
 * @param dofree 0 - don't free a; > 0 cholmod_free a; < 0 R_Free a
 *
 * @return SEXP containing a copy of a
 */
SEXP chm_factor_to_SEXP(CHM_FR f, int dofree)
{
    SEXP ans;
    int *dims, *type;
    char *class = (char*) NULL;	/* -Wall */

#define DOFREE_MAYBE					\
    if(dofree) {					\
	if (dofree > 0) cholmod_free_factor(&f, &c);	\
	else /* dofree < 0 */ R_Free(f);			\
    }

    if(f->minor < f->n) {
	DOFREE_MAYBE;
	error(_("CHOLMOD factorization was unsuccessful"));
	// error(_("previous CHOLMOD factorization was unsuccessful"));
    }

    switch(f->xtype) {
    case CHOLMOD_REAL:
	class = f->is_super ? "dCHMsuper" : "dCHMsimpl";
	break;
    case CHOLMOD_PATTERN:
	class = f->is_super ? "nCHMsuper" : "nCHMsimpl";
	break;
    default:
	DOFREE_MAYBE;
	error(_("f->xtype of %d not recognized"), f->xtype);
    }

    ans = PROTECT(newObject(class));
    dims = INTEGER(ALLOC_SLOT(ans, Matrix_DimSym, INTSXP, 2));
    dims[0] = dims[1] = f->n;
				/* copy component of known length */
    type = INTEGER(ALLOC_SLOT(ans, install("type"), INTSXP, 6));
    type[0] = f->ordering; type[1] = f->is_ll;
    type[2] = f->is_super; type[3] = f->is_monotonic;
    type[4] = f->maxcsize; type[5] = f->maxesize;
    Memcpy(INTEGER(ALLOC_SLOT(ans, install("colcount"), INTSXP, f->n)),
	   (int*)f->ColCount, f->n);
    if (f->ordering != CHOLMOD_NATURAL) {
	Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_permSym, INTSXP, f->n)),
	       (int*)f->Perm, f->n);
    }
    if (f->is_super) {
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("super"), INTSXP, f->nsuper + 1)),
	       (int*)f->super, f->nsuper+1);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("pi"), INTSXP, f->nsuper + 1)),
	       (int*)f->pi, f->nsuper + 1);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("px"), INTSXP, f->nsuper + 1)),
	       (int*)f->px, f->nsuper + 1);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("s"), INTSXP, f->ssize)),
	       (int*)f->s, f->ssize);
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, f->xsize)),
	       (double*)f->x, f->xsize);
    } else {
	Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP, f->nzmax)),
	       (int*)f->i, f->nzmax);
	Memcpy(INTEGER(ALLOC_SLOT(ans, Matrix_pSym, INTSXP, f->n + 1)),
	       (int*)f->p, f->n + 1);
	Memcpy(REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, f->nzmax)),
	       (double*)f->x, f->nzmax);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("nz"), INTSXP, f->n)),
	       (int*)f->nz, f->n);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("nxt"), INTSXP, f->n + 2)),
	       (int*)f->next, f->n + 2);
	Memcpy(INTEGER(ALLOC_SLOT(ans, install("prv"), INTSXP, f->n + 2)),
	       (int*)f->prev, f->n + 2);

    }
    DOFREE_MAYBE;
    UNPROTECT(1);
    return ans;
}
#undef DOFREE_MAYBE

/**
 * Evaluate the logarithm of the square of the determinant of L
 *
 * @param f pointer to a CHMfactor object
 *
 * @return log(det(L)^2)
 *
 */
double chm_factor_ldetL2(CHM_FR f)
{
    int i, j, p;
    double ans = 0;

    if (f->is_super) {
	int *lpi = (int*)(f->pi), *lsup = (int*)(f->super);
	for (i = 0; i < f->nsuper; i++) { /* supernodal block i */
	    int nrp1 = 1 + lpi[i + 1] - lpi[i],
		nc = lsup[i + 1] - lsup[i];
	    double *x = (double*)(f->x) + ((int*)(f->px))[i];

	    for (R_xlen_t jn = 0, j = 0; j < nc; j++, jn += nrp1) { // jn := j * nrp1
		ans += 2 * log(fabs(x[jn]));
	    }
	}
    } else {
	int *li = (int*)(f->i), *lp = (int*)(f->p);
	double *lx = (double *)(f->x);

	for (j = 0; j < f->n; j++) {
	    for (p = lp[j]; li[p] != j && p < lp[j + 1]; p++) {};
	    if (li[p] != j) {
		error(_("diagonal element %d of Cholesky factor is missing"), j);
		break;		/* -Wall */
	    }
	    ans += log(lx[p] * ((f->is_ll) ? lx[p] : 1.));
	}
    }
    return ans;
}

/**
 * Update the numerical values in the factor f as A + mult * I, if A is
 * symmetric, otherwise AA' + mult * I
 *
 * @param f pointer to a CHM_FR object.  f is updated upon return.
 * @param A pointer to a CHM_SP object, possibly symmetric
 * @param mult multiple of the identity to be added to A or AA' before
 * decomposing.
 *
 * @note: A and f must be compatible.  There is no check on this
 * here.  Incompatibility of A and f will cause the CHOLMOD functions
 * to take an error exit.
 *
 */
CHM_FR chm_factor_update(CHM_FR f, CHM_SP A, double mult)
{
    int ll = f->is_ll;
    double mm[2] = {0, 0};
    mm[0] = mult;
    // NB: Result depends if A is "dsC" or "dgC"; the latter case assumes we mean AA' !!!
    if (!cholmod_factorize_p(A, mm, (int*)NULL, 0 /*fsize*/, f, &c))
	/* -> ./CHOLMOD/Cholesky/cholmod_factorize.c */
	error(_("cholmod_factorize_p failed: status %d, minor %d of ncol %d"),
	      c.status, f->minor, f->n);
    if (f->is_ll != ll)
	if(!cholmod_change_factor(f->xtype, ll, f->is_super, 1 /*to_packed*/,
				  1 /*to_monotonic*/, f, &c))
	    error(_("cholmod_change_factor failed"));
    return f;
}
