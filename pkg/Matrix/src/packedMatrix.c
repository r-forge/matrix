#ifdef __GLIBC__
#define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "packedMatrix.h"

#define PM_AR21_UP(i, j) i + (j * (j + 1)) / 2
#define PM_AR21_LO(i, j, n2) i + (j * (n2 - j - 1)) / 2

/* An alternative to the existing utility 'symmetric_DimNames' 
   enabling users to avoid a copy in cases where it is avoidable ...
   perhaps the existing utility can be rebuilt around this one 
   so that they remain in sync ...
*/
void fast_symmetric_DimNames(SEXP dn, SEXP *vec, SEXP *nm)
{
    *vec = VECTOR_ELT(dn, 0);
    if (isNull(*vec))
    {
	*vec = VECTOR_ELT(dn, 1);
    }
    *nm = getAttrib(dn, R_NamesSymbol);
    if (!isNull(*nm))
    {
	*nm = STRING_ELT(*nm, 0);
	if (LENGTH(*nm) == 0) /* note that LENGTH(NA_STRING) = 2 */
	{
	    *nm = STRING_ELT(*nm, 1);
	}
    }
}

/* FIXME: could avoid some duplication by calling 'symmetricMatrix_validate'
   or 'triangularMatrix_validate' conditional on existence of 'diag' slot ...
   would still need to check length(.@x) though
*/
SEXP packedMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);
    if (LENGTH(val) != 2) {
        return mkString(_("'Dim' slot does not have length 2"));
    }
    int n = INTEGER(val)[0];
    if (INTEGER(val)[1] != n) {
        return mkString(_("Matrix is not square"));
    }
    val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym), "LU", "uplo");
    if (isString(val)) {
        return val;
    }
    if (R_has_slot(obj, Matrix_diagSym)) {
	val = check_scalar_string(GET_SLOT(obj, Matrix_diagSym), "NU", "diag");
	if (isString(val)) {
	    return val;
	}
    }
    val = GET_SLOT(obj, Matrix_xSym);
    if (LENGTH(val) != PACKED_LENGTH(n)) {
        return mkString(_("'x' slot does not have length 'n*(n+1)/2', n=Dim[1]"));
    }
    return ScalarLogical(1);
}

#define PM_T_LOOP(_px0_, _px1_)						\
    /* Loop over [i,j] 0-indices of stored triangle of _result_ */	\
    if (up) {								\
	for (int j = 0, k = 0; j < n; ++j) {				\
	    for (int i = j; i < n; ++i)	{				\
		(_px1_)[k++] = (_px0_)[PM_AR21_UP(j, i)];		\
	    }								\
	}								\
    }									\
    else								\
    {									\
	int n2 = 2 * n;							\
	for (int j = 0, k = 0; j < n; ++j) {				\
	    for (int i = 0; i <= j; ++i) {				\
		(_px1_)[k++] = (_px0_)[PM_AR21_LO(j, i, n2)];		\
	    }								\
	}								\
    }

#define PM_T(_datatype_, _sexptype_, _accessor_)			\
    x1 = PROTECT(allocVector(_sexptype_, LENGTH(x0)));			\
    _datatype_ *px0 = _accessor_(x0);					\
    _datatype_ *px1 = _accessor_(x1);					\
    PM_T_LOOP(px0, px1);						\
    SET_SLOT(res, Matrix_xSym, x1);					\
    UNPROTECT(1)

/* t(x) */
SEXP packedMatrix_t(SEXP obj)
{
    /* Initialize result of same class */
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(class_P(obj)));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? upper -> lower : lower -> upper */
    Rboolean up = uplo_P(obj)[0] == 'U';

    SEXP x0 = GET_SLOT(obj, Matrix_xSym);
    if (n > 1) {
	/* Permute 'x' slot */
	SEXP x1;
	if (isReal(x0)) {
	    PM_T(double, REALSXP, REAL);
	} else {
	    PM_T(int, LGLSXP, LOGICAL);
	}
    } else {
	/* Preserve 'x' slot */
	SET_SLOT(res, Matrix_xSym, x0);
    }

    /* Toggle 'uplo' slot */
    SET_SLOT(res, Matrix_uploSym, mkString(up ? "L" : "U"));
    /* Preserve 'Dim' slot */
    SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
    /* Reverse 'Dimnames' slot and (if not absent) 'names(Dimnames)' */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn1, 0, VECTOR_ELT(dn0, 1));
    SET_VECTOR_ELT(dn1, 1, VECTOR_ELT(dn0, 0));
    SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
    if (!isNull(ndn0)) {
	SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(ndn1, 0, STRING_ELT(ndn0, 1));
	SET_STRING_ELT(ndn1, 1, STRING_ELT(ndn0, 0));
	setAttrib(dn1, R_NamesSymbol, ndn1);
	UNPROTECT(1);
    }
    SET_SLOT(res, Matrix_DimNamesSym, dn1);
    UNPROTECT(2);
    return res;
}

#undef PM_T
#undef PM_T_LOOP

#define PM_D_G_UDIAG(_pres_, _one_)		\
    for (int j = 0; j < n; ++j) {		\
	(_pres_)[j] = _one_;			\
    }

#define PM_D_G_NDIAG(_pres_, _px_)					\
    for (int j = 0, pos = 0; j < n; pos += (up ? 1 + (++j) : n - (j++))) { \
	(_pres_)[j] = (_px_)[pos];					\
    }

#define PM_D_G(_datatype_, _sexptype_, _accessor_, _one_)	\
    res = PROTECT(allocVector(_sexptype_, n));			\
    _datatype_ *pres = _accessor_(res);				\
    if (utr) {							\
	PM_D_G_UDIAG(pres, _one_);				\
    } else {							\
	_datatype_ *px = _accessor_(x);				\
	PM_D_G_NDIAG(pres, px);					\
    }

/* 'packedMatrix_diag_get' handles 'names' precisely as documented
   in '?diag', whereas current '*_getDiag' machinery neglects 'names'.

   It combines the work of 6 existing C utilities:
   '[dl][st]pMatrix_detDiag' and '[dl]_packed_getDiag'.

   It does so at the cost of conditioning on class
   (double/logical, triangular/symmetric), 
   but that cost is minimal
   (1 call to isReal, 1 call to R_has_slot) ...
*/
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL) {
	error(_("'names' must be TRUE or FALSE"));
    }

    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';

    SEXP res;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (isReal(x)) {
	PM_D_G(double, REALSXP, REAL, 1.0);
    } else {
	PM_D_G(int, LGLSXP, LOGICAL, 1);
    }

    if (do_nms) {
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym);
	SEXP rn = VECTOR_ELT(dn, 0);
	SEXP cn = VECTOR_ELT(dn, 1);
	if (isNull(rn)) {
	    if (!tr && !isNull(cn)) {
		setAttrib(res, R_NamesSymbol, cn);
	    }
	} else {
	    if (!tr || R_compute_identical(rn, cn, 16)) {
		setAttrib(res, R_NamesSymbol, rn);
	    }
	}
    }
    UNPROTECT(1);
    return res;
}

#undef PM_D_G
#undef PM_D_G_NDIAG
#undef PM_D_G_UDIAG

#define PM_D_S_ONE(_px_, _d_)						\
    for (int j = 0, pos = 0; j < n; pos += (up ? 1 + (++j) : n - (j++))) { \
	(_px_)[pos] = _d_;						\
    }

#define PM_D_S_FULL(_px_, _pval_)					\
    for (int j = 0, pos = 0; j < n; pos += (up ? 1 + (++j) : n - (j++))) { \
	(_px_)[pos] = (_pval_)[j];					\
    }

#define PM_D_S(_datatype_, _accessor_)		\
    _datatype_ *px = _accessor_(x);		\
    _datatype_ *pval = _accessor_(val);		\
    if (nv1) {					\
	_datatype_ d = pval[0];			\
	PM_D_S_ONE(px, d);			\
    } else {					\
	PM_D_S_FULL(px, pval);			\
    }

/* 'packedMatrix_diag_set' handles coercions similarly to base R's
   'diag<-' and performs explicit tests for type compatibility,
   whereas current '*_setDiag' machinery does neither. It also does
   more to avoid memory allocation by copying slot values only if 
   necessary.

   It combines the work of 8 existing C utilities:
   '[dl][st]pMatrix_setDiag' and '(tr_)?[dl]_packed_setDiag'.

   Like 'packedMatrix_diag_get', it conditions on class, but the cost
   of doing so is minimal ...
*/
SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int nv = LENGTH(val);
    Rboolean nv1 = (nv == 1);
    if (!(nv1 || nv == n)) {
	error(_("replacement diagonal has wrong length"));
    }

    /* Initialize result object of same class */
    SEXP res;
    int nprotect = 0;
    if (MAYBE_REFERENCED(obj)) { /* MAYBE_SHARED seems less safe ... */
	res = PROTECT(NEW_OBJECT_OF_CLASS(class_P(obj))); ++nprotect;
	SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
	slot_dup(res, obj, Matrix_xSym);
    } else {
	res = obj;
    }

    /* Toggle 'diag' slot when assigning to unit diagonal "[dln]tpMatrix",
       _even if_ RHS of assignment is unity
     */
    if (Diag_P(res)[0] == 'U') {
	SET_SLOT(res, Matrix_diagSym, mkString("N"));
    }
    /* Reset 'factors' slot when assigning to "[dln]spMatrix" */
    if (R_has_slot(res, Matrix_factorSym) &&
	LENGTH(GET_SLOT(res, Matrix_factorSym)) > 0) {
	SET_SLOT(res, Matrix_factorSym, allocVector(VECSXP, 0));
    }

    /* ? double result : logical result */
    Rboolean dbl = TRUE;
    /* ? upper : lower */
    Rboolean up = uplo_P(res)[0] == 'U';

    /* Test that LHS and RHS of assignment have compatible types,
       coercing one or the other from logical to double if necessary
     */
    SEXP x = GET_SLOT(res, Matrix_xSym);
    switch (TYPEOF(x)) {
    case LGLSXP:
	switch (TYPEOF(val)) {
	case LGLSXP:
	    dbl = FALSE;
	    break;
	case INTSXP:
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
	case REALSXP:
	{
	    /* [ln][st]pMatrix -> d[st]pMatrix */
	    SEXP strcl = getAttrib(res, R_ClassSymbol);
	    char *cl = strdup(CHAR(STRING_ELT(strcl, 0)));
	    cl[0] = 'd';
	    SET_STRING_ELT(strcl, 0, mkChar(cl));
	    free(cl);
	    x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
	    SET_SLOT(res, Matrix_xSym, x);
	    break;
	}
	default:
	    error(_("replacement diagonal has incompatible type '%s'"),
		  type2char(TYPEOF(val)));
	}
	break;
    case REALSXP:
	switch (TYPEOF(val)) {
	case LGLSXP:
	case INTSXP:
	    /* logical, integer -> double */
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
	case REALSXP:
	    break;
	default:
	    error(_("replacement diagonal has incompatible type '%s'"),
		  type2char(TYPEOF(val)));
	}
	break;
    default:
	error(_("'x' slot is not of type '%s' or '%s', which should never happen; please report"),
	      type2char(LGLSXP), type2char(REALSXP));
    }

    if (dbl) {
	PM_D_S(double, REAL);
    } else {
	PM_D_S(int, LOGICAL);
    }
    UNPROTECT(nprotect);
    return res;
}

#undef PM_D_S
#undef PM_D_S_FULL
#undef PM_D_S_ONE

#define PM_SUB0_START					\
    SEXP x = GET_SLOT(obj, Matrix_xSym);		\
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];	\
    int n2 = 2 * n;					\
    int *pindex = INTEGER(index);			\
    /* ? triangular : symmetric */			\
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);	\
    /* ? unit triangular : other */			\
    Rboolean utr = tr && diag_P(obj)[0] == 'U';		\
    /* ? upper : lower */				\
    Rboolean up = uplo_P(obj)[0] == 'U'

#define PM_SUB0_IJ(i, j, _px_, _zero_, _one_)			\
    (utr && i == j						\
     ? _one_							\
     : (up							\
	? (i <= j						\
	   ? (_px_)[PM_AR21_UP(i, j)]				\
	   : (tr						\
	      ? _zero_						\
	      : (_px_)[PM_AR21_UP(j, i)]))			\
	: (i >= j						\
	   ? (_px_)[PM_AR21_LO(i, j, n2)]			\
	   : (tr						\
	      ? _zero_						\
	      : (_px_)[PM_AR21_LO(j, i, n2)]))))

#define PM_SUB0_LOOP(_pres_, _px_, _na_, _zero_, _one_)			\
    for (int k = 0, i, j, pos; k < nindex; ++k)	{			\
	pos = pindex[k];						\
	if (pos == NA_INTEGER) {					\
	    (_pres_)[k] = _na_;						\
	} else {							\
	    pos -= 1; /* 1-index -> 0-index */				\
	    i = pos % n;						\
	    j = pos / n;						\
	    (_pres_)[k] = PM_SUB0_IJ(i, j, _px_, _zero_, _one_);	\
	}								\
    }

#define PM_SUB0(_datatype_, _sexptype_, _accessor_, _na_, _zero_, _one_) \
    SEXP res = PROTECT(allocVector(_sexptype_, nindex));		\
    _datatype_ *pres = _accessor_(res);					\
    _datatype_ *px = _accessor_(x);					\
    PM_SUB0_LOOP(pres, px, _na_, _zero_, _one_);			\
    UNPROTECT(1);							\
    return res

#define PM_SUB0_END						\
    if (isReal(x)) {						\
	PM_SUB0(double, REALSXP, REAL, NA_REAL, 0.0, 1.0);	\
    } else {							\
	PM_SUB0(int, LGLSXP, LOGICAL, NA_LOGICAL, 0, 1);	\
    }

/* 'x[index]' where 'index' is an integer vector supplying 
   integers in 'c(1:(n*n), NA)' _only_
*/
SEXP packedMatrix_sub0_1ary(SEXP obj, SEXP index)
{
    PM_SUB0_START;
    int nindex = LENGTH(index);
    PM_SUB0_END;
}

#undef PM_SUB0_LOOP

#define PM_SUB0_LOOP(_pres_, _px_, _na_, _zero_, _one_)			\
    for (int k = 0, i, j; k < nindex; ++k) {				\
	i = pindex[k];							\
	j = pindex[k + nindex];						\
	if (i == NA_INTEGER || j == NA_INTEGER)	{			\
	    (_pres_)[k] = _na_;						\
	} else {							\
	    i -= 1;							\
	    j -= 1;							\
	    (_pres_)[k] = PM_SUB0_IJ(i, j, _px_, _zero_, _one_);	\
	}								\
    }

/* 'x[index]' where 'index' is a 2-column integer matrix supplying 
   integers in 'c(1:n, NA)' _only_
*/
SEXP packedMatrix_sub0_2ary(SEXP obj, SEXP index)
{
    PM_SUB0_START;
    int nindex = INTEGER(getAttrib(index, R_DimSymbol))[0];
    PM_SUB0_END;
}

#undef PM_SUB0_END
#undef PM_SUB0
#undef PM_SUB0_LOOP
#undef PM_SUB0_START

#define PM_SUB1_LOOP(_px0_, _px1_, _na_, _zero_, _one_)			\
    int i, incr;							\
    incr = (do_col ? 1 : nindex);					\
    for (int k = 0, pos = 0; k < nindex; ++k) {				\
	i = pindex[k];							\
	if (i == NA_INTEGER) {						\
	    if (do_nms) {						\
		SET_STRING_ELT(destnms, k, NA_STRING);			\
	    }								\
	    for (int j = 0; j < n; ++j, pos += incr) {			\
		px1[pos] = _na_;					\
	    }								\
	} else {							\
	    i -= 1; /* 1-index -> 0-index */				\
	    if (do_nms)	{						\
		SET_STRING_ELT(destnms, k, STRING_ELT(srcnms, i));	\
	    }								\
	    /* FIXME: is there a nice way to reduce repetition */	\
	    /* here without obfuscating too much ... ?? */		\
	    if (tr) {							\
		if (up) {						\
		    if (do_col)	{					\
			for (int j = 0; j <= i; ++j, pos += incr) {	\
			    px1[pos] = px0[PM_AR21_UP(j, i)];		\
			}						\
			for (int j = i + 1; j < n; ++j, pos += incr) {	\
			    px1[pos] = _zero_;				\
			}						\
		    } else {						\
			for (int j = 0; j < i; ++j, pos += incr) {	\
			    px1[pos] = _zero_;				\
			}						\
			for (int j = i; j < n; ++j, pos += incr) {	\
			    px1[pos] = px0[PM_AR21_UP(i, j)];		\
			}						\
		    }							\
		} else {						\
		    if (do_col) {					\
			for (int j = 0; j < i; ++j, pos += incr) {	\
			    px1[pos] = _zero_;				\
			}						\
			for (int j = i; j < n; ++j, pos += incr) {	\
			    px1[pos] = px0[PM_AR21_LO(j, i, n2)];	\
			}						\
		    } else {						\
			for (int j = 0; j <= i; ++j, pos += incr) {	\
			    px1[pos] = px0[PM_AR21_LO(i, j, n2)];	\
			}						\
			for (int j = i + 1; j < n; ++j, pos += incr) {	\
			    px1[pos] = _zero_;				\
			}						\
		    }							\
		}							\
		if (utr) {						\
		    px1[pos - (n - i) * incr] = _one_;			\
		}							\
	    } else {							\
		if (up) {						\
		    for (int j = 0; j < i; ++j, pos += incr) {		\
		        px1[pos] = px0[PM_AR21_UP(j, i)];		\
		    }							\
		    for (int j = i; j < n; ++j, pos += incr) {		\
			px1[pos] = px0[PM_AR21_UP(i, j)];		\
		    }							\
		} else {						\
		    for (int j = 0; j <= i; ++j, pos += incr) {		\
			px1[pos] = px0[PM_AR21_LO(i, j, n2)];		\
		    }							\
		    for (int j = i + 1; j < n; ++j, pos += incr) {	\
			px1[pos] = px0[PM_AR21_LO(j, i, n2)];		\
		    }							\
		}							\
	    }								\
	}								\
	if (!do_col) {							\
	    pos = (pos + 1) % nindex;					\
	}								\
    }

#define PM_SUB1(_datatype_, _sexptype_, _accessor_, _na_, _zero_, _one_) \
    x1 = PROTECT(allocVector(_sexptype_, n * nindex));			\
    _datatype_ *px0 = _accessor_(x0);					\
    _datatype_ *px1 = _accessor_(x1);					\
    PM_SUB1_LOOP(px0, px1, _na_, _zero_, _one_);			\
    SET_SLOT(res, Matrix_xSym, x1);					\
    UNPROTECT(1)

/* 'x[index, ]' and 'x[, index]' where 'i' is an integer vector supplying
   integers in 'c(1:n, NA)' _only_
*/
SEXP packedMatrix_sub1(SEXP obj, SEXP index, SEXP drop, SEXP col)
{
    int do_drop = asLogical(drop);
    if (do_drop == NA_LOGICAL) {
	error(_("'drop' must be TRUE or FALSE"));
    }

    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int n2 = 2 * n;
    int nindex = LENGTH(index);
    int *pindex = INTEGER(index);
    int nprotect = 0;
    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';
    /* ? index columns : index rows */
    Rboolean do_col = LOGICAL(col)[0];
    /* margin of the subset */
    int mar = !do_col;

    /* Initialize result of same type but "general" class */
    char *cl = strdup(class_P(obj));
    cl[1] = 'g';
    cl[2] = 'e';
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl)); ++nprotect;
    free(cl);
    
    /* Set 'Dim' slot */
    SEXP d1 = PROTECT(GET_SLOT(res, Matrix_DimSym));
    INTEGER(d1)[mar] = n;
    INTEGER(d1)[!mar] = nindex;
    UNPROTECT(1);

    /* Set 'Dimnames' slot and (if not absent) 'names(Dimnames)' ... */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = PROTECT(GET_SLOT(res, Matrix_DimNamesSym));
    SEXP srcnms, destnms; /* for the "stretched" vector of names */
    if (tr) {
	SET_VECTOR_ELT(dn1, mar, VECTOR_ELT(dn0, mar));
	srcnms = VECTOR_ELT(dn0, !mar);
	SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
	if (!isNull(ndn0)) {
	    setAttrib(dn1, R_NamesSymbol, ndn0);
	}
    } else {
	SEXP s;
	fast_symmetric_DimNames(dn0, &srcnms, &s);
	SET_VECTOR_ELT(dn1, mar, srcnms);
	if (!isNull(s)) {
	    SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	    SET_STRING_ELT(ndn1, 0, s);
	    SET_STRING_ELT(ndn1, 1, s);
	    setAttrib(dn1, R_NamesSymbol, ndn1);
	    UNPROTECT(1);
	}
    }
    Rboolean do_nms = !isNull(srcnms) && nindex > 0;
    if (do_nms) {
	destnms = PROTECT(allocVector(STRSXP, nindex)); ++nprotect;
	SET_VECTOR_ELT(dn1, !mar, destnms);
    }
    UNPROTECT(1);
    
    SEXP x0, x1;
    x0 = GET_SLOT(obj, Matrix_xSym);
    if (isReal(x0)) {
	PM_SUB1(double, REALSXP, REAL, NA_REAL, 0.0, 1.0);
    } else {
	PM_SUB1(int, LGLSXP, LOGICAL, NA_LOGICAL, 0, 1);
    }

    /* Drop dimensions in this special case */
    if (do_drop && nindex == 1) {
	res = GET_SLOT(res, Matrix_xSym);
	if (do_nms) {
	    SEXP nms = VECTOR_ELT(dn1, mar);
	    if (!isNull(nms)) {
		setAttrib(res, R_NamesSymbol, nms);
	    }
	}
    }
    UNPROTECT(nprotect);
    return res;
}

#undef PM_SUB1
#undef PM_SUB1_LOOP
