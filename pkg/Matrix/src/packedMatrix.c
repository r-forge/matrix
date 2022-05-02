#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "packedMatrix.h"

/* unpack(x) */
SEXP packedMatrix_unpack(SEXP from)
{
    static const char *valid_from[] = {
	"dppMatrix", /* must match before "dspMatrix" */
	"dspMatrix", "lspMatrix", "nspMatrix",
	"pCholesky", "pBunchKaufman", /* must match before "dtpMatrix" */
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    static const char *valid_to[] = {
	"dpoMatrix",
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"Cholesky", "BunchKaufman",
	"dtrMatrix", "ltrMatrix", "ntrMatrix", ""};
    int icl = R_check_class_etc(from, valid_from);
    if (icl < 0) {
	PM_ERROR_INVALID_CLASS(class_P(from), "unpack");
    }
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid_to[icl]));
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int n = INTEGER(dim)[0];
    
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (icl <= 3) {
	/* .syMatrix */
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    } else {
	/* .trMatrix */
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	if (icl == 5) {
	    /* pBunchKaufman */
	    SET_SLOT(to, Matrix_permSym, GET_SLOT(from, Matrix_permSym));
	}
    }

    SEXP x_from = GET_SLOT(from, Matrix_xSym);

#define PM_UNPACK(_PREFIX_, _CTYPE_, _SEXPTYPE_, _PTR_)			\
    do {								\
	R_xlen_t len = n * (R_xlen_t) n;				\
	SEXP x_to = PROTECT(allocVector(_SEXPTYPE_, len));		\
	_PREFIX_ ## dense_unpack(_PTR_(x_to), _PTR_(x_from), n,		\
				 (*uplo_P(to) == 'U') ? UPP : LOW,	\
				 NUN /* for speed */);			\
	SET_SLOT(to, Matrix_xSym, x_to);				\
	UNPROTECT(1);							\
    } while (0)
    
    switch (TYPEOF(x_from)) {
    case REALSXP: /* d..Matrix */
	PM_UNPACK(d, double, REALSXP, REAL);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_UNPACK(i, int, LGLSXP, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	PM_UNPACK(i, int, INTSXP, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
    	PM_UNPACK(z, Rcomplex, CPLXSXP, COMPLEX);
	break;
    default:
	PM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x_from), "unpack");
	break;
    }
    
#undef PM_UNPACK
    
    UNPROTECT(1);
    return to;
}

#define PM_IS_DI(_RES_, _X_, _N_, _UP_, _METHOD_)			\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_packed_is_diagonal(REAL(_X_), _N_, _UP_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_packed_is_diagonal(LOGICAL(_X_), _N_, _UP_); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_packed_is_diagonal(INTEGER(_X_), _N_, _UP_); \
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_packed_is_diagonal(COMPLEX(_X_), _N_, _UP_); \
	    break;							\
	default:							\
	    PM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(_X_), _METHOD_);	\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

/* isSymmetric(x, tol = 0) */ 
SEXP packedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    static const char *valid[] = {
	"dspMatrix", "lspMatrix", "nspMatrix", /* be fast */
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	PM_ERROR_INVALID_CLASS(class_P(obj), "is_symmetric");
	return R_NilValue;
    } else if (icl <= 2) {
	/* .spMatrix: symmetric by definition */
	return ScalarLogical(1);
    } else {
	/* .tpMatrix: symmetric iff diagonal */
	if (asLogical(checkDN) != 0 &&
	    !DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym))) {
	    return ScalarLogical(0);
	}
	Rboolean res = FALSE, up = (*uplo_P(obj) == 'U') ? TRUE : FALSE;
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
	PM_IS_DI(res, x, n, up, "is_symmetric");
	return ScalarLogical(res);
    }
}

#define RETURN_TRUE(_KIND_)						\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

/* isTriangular(x, upper) */
SEXP packedMatrix_is_triangular(SEXP obj, SEXP upper)
{
    static const char *valid[] = {
	"dtpMatrix", "ltpMatrix", "ntpMatrix", /* be fast */
	"dspMatrix", "lspMatrix", "nspMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	PM_ERROR_INVALID_CLASS(class_P(obj), "is_triangular");
    }

    SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
    Rboolean have_upper = (*CHAR(STRING_ELT(uplo, 0)) == 'U') ? TRUE : FALSE;
    int need_upper = asLogical(upper);

#define IF_DIAGONAL							\
    Rboolean res = FALSE;						\
    SEXP x = GET_SLOT(obj, Matrix_xSym);				\
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];			\
    PM_IS_DI(res, x, n, have_upper, "is_triangular");			\
    if (res)
    
    if (icl <= 2) {
	/* .tpMatrix: be fast if 'upper', 'uplo' agree; else need diagonal */
	if (need_upper == NA_LOGICAL) {
	    RETURN_TRUE(uplo);
	} else if ((need_upper != 0 &&  have_upper) ||
		   (need_upper == 0 && !have_upper)) {
	    return ScalarLogical(1);
	} else {
	    IF_DIAGONAL {
		return ScalarLogical(1);
	    }
	}
    } else {
	/* .spMatrix: triangular iff diagonal */
	IF_DIAGONAL {
	    if (need_upper == NA_LOGICAL) {
		RETURN_TRUE(mkString("U"));
	    } else {
		return ScalarLogical(1);
	    }
	}
    }

#undef IF_DIAGONAL

    return ScalarLogical(0);
}

#undef RETURN_TRUE

/* isDiagonal(x) */
SEXP packedMatrix_is_diagonal(SEXP obj)
{
    /* _Not_ checking class of 'obj' */
    Rboolean res = FALSE, up = (*uplo_P(obj) == 'U') ? TRUE : FALSE;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    PM_IS_DI(res, x, n, up, "is_diagonal");
    return ScalarLogical(res);
}

#undef PM_IS_DI

/* t(x) */
SEXP packedMatrix_t(SEXP obj)
{
    /* Initialize result of same class, except for pBunchKaufman */
    const char *cl = class_P(obj);
    SEXP res, x0 = GET_SLOT(obj, Matrix_xSym);
    if (TYPEOF(x0) == REALSXP) {
	static const char *valid[] = { "pBunchKaufman", "" };
	if (!R_check_class_etc(obj, valid)) {
	    cl = "dtpMatrix";
	}
    }
    res = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    const char *uplo = uplo_P(obj), *diag = Diag_P(obj);

#define PM_T(_CTYPE_, _SEXPTYPE_, _PTR_)				\
    do {								\
	SEXP x1 = PROTECT(allocVector(_SEXPTYPE_, XLENGTH(x0)));	\
	_CTYPE_ *px0 = _PTR_(x0);					\
	_CTYPE_ *px1 = _PTR_(x1);					\
	int i, j;							\
	if (*uplo == 'U') {						\
	    for (j = 0; j < n; ++j)					\
		for (i = j; i < n; ++i)					\
		    *(px1++) = *(px0 + PM_AR21_UP(j, i));		\
	} else {							\
	    R_xlen_t n2 = ((R_xlen_t) n) * 2;				\
	    for (j = 0; j < n; ++j)					\
		for (i = 0; i <= j; ++i)				\
		    *(px1++) = *(px0 + PM_AR21_LO(j, i, n2));		\
	}								\
	SET_SLOT(res, Matrix_xSym, x1);					\
	UNPROTECT(1);							\
    } while (0)
    
    if (n > 1) {
	/* Permute 'x' slot */
	switch (TYPEOF(x0)) {
	case REALSXP: /* d..Matrix */
	    PM_T(double, REALSXP, REAL);
	    break;
	case LGLSXP: /* [ln]..Matrix */
	    PM_T(int, LGLSXP, LOGICAL);
	    break;
	case INTSXP: /* i..Matrix */
	    PM_T(int, INTSXP, INTEGER);
	    break;
	case CPLXSXP: /* z..Matrix */
	    PM_T(Rcomplex, CPLXSXP, COMPLEX);
	    break;
	default:
	    PM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x0), "t");
	    break;
	}
    } else {
	/* Preserve 'x' slot */
	SET_SLOT(res, Matrix_xSym, x0);
    }

#undef PM_T
    
    /* Preserve 'Dim' slot */
    SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
    /* Preserve or reverse 'Dimnames' slot and 'names(Dimnames)' 
       (preserving if symmetric) 
    */
    if (*diag == ' ') {
	SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
    } else {
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
	UNPROTECT(1);
    }
    /* Toggle 'uplo' slot */
    SET_SLOT(res, Matrix_uploSym, mkString(*uplo == 'U' ? "L" : "U"));
    /* Preserve 'diag' slot */
    if (*diag == 'U') {
	SET_SLOT(res, Matrix_diagSym, GET_SLOT(obj, Matrix_diagSym));
    }
    /* NB: Nothing to do for 'factors' slot: prototype is already list() ...
       FIXME: However, it would be much better to transpose and reverse
       the factors in each factorization
    */
    
    UNPROTECT(1);
    return res;
}

/* diag(x, names) */
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL) {
	error(_("'names' must be TRUE or FALSE"));
    }

    SEXP res = R_NilValue, x = GET_SLOT(obj, Matrix_xSym);
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    const char *uplo = uplo_P(obj), *diag = Diag_P(obj);
    
#define PM_D_G(_CTYPE_, _SEXPTYPE_, _PTR_, _INIT_ONE_, _ONE_)		\
    do {								\
	res = PROTECT(allocVector(_SEXPTYPE_, n));			\
	_CTYPE_ *pres = _PTR_(res);					\
	if (*diag == 'U') {						\
	    _INIT_ONE_;							\
	    for (int j = 0; j < n; ++j)					\
		*(pres++) = _ONE_;					\
	} else {							\
	    _CTYPE_ *px = _PTR_(x);					\
	    if (*uplo == 'U') {						\
		for (int j = 0; j < n; px += (++j)+1)			\
		    *(pres++) = *px;					\
	    } else {							\
		for (int j = 0; j < n; px += n-(j++))			\
		    *(pres++) = *px;					\
	    }								\
	}								\
    } while (0)

    switch (TYPEOF(x)) {
    case REALSXP: /* d..Matrix */
	PM_D_G(double, REALSXP, REAL, , 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_D_G(int, LGLSXP, LOGICAL, , 1);
	break;
    case INTSXP: /* i..Matrix */
	PM_D_G(int, INTSXP, INTEGER, , 1);
	break;
    case CPLXSXP: /* z..Matrix */
	PM_D_G(Rcomplex, CPLXSXP, COMPLEX,
	       Rcomplex one; one.r = 1.0; one.i = 0.0;, one);
	break;
    default:
	PM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x), "diag_get");
	break;
    }

#undef PM_D_G
    
    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames'
	*/
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	Rboolean symmetric = (*diag == ' ') ? TRUE : FALSE;
	if (isNull(cn)) {
	    if (symmetric && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (symmetric || (!isNull(rn) &&
			      (rn == cn || equal_string_vectors(rn, cn, n)))) {
		setAttrib(res, R_NamesSymbol, cn);
	    }
	}
    }
    UNPROTECT(1);
    return res;
}

/* diag(x) <- value */
/* MJ: Is NO_REFERENCES(x@x) ever true?? (in which case do in-place??) */
SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0], nv = LENGTH(val);
    if (nv != 1 && nv != n) {
	error(_("replacement diagonal has wrong length"));
    }
    
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x), tv = TYPEOF(val);
    if (tx < LGLSXP || tx > CPLXSXP) {
	PM_ERROR_INVALID_SLOT_TYPE("x", tx, "diag_set");
    }
    if (tv < LGLSXP || tv > REALSXP) {
	/* Upper bound can become CPLXSXP once we have proper zdenseMatrix */
	error(_("replacement diagonal has incompatible type \"%s\""),
	      type2char(tv));
    }
    
    static const char *valid[] = {
	"dspMatrix", "dtpMatrix",
	"lspMatrix", "ltpMatrix",
	"nspMatrix", "ntpMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	PM_ERROR_INVALID_CLASS(class_P(obj), "diag_set");
    }
    const char *cl = valid[icl];
    
    SEXP res = R_NilValue;
    int nprotect = 0;

    /* Allocate and coerce as necessary */
    if (tv <= tx) {
	if (tv < tx) {
	    val = PROTECT(coerceVector(val, tx)); ++nprotect;
	    tv = tx;
	}
	res = PROTECT(NEW_OBJECT_OF_CLASS(cl)); ++nprotect;
	SET_SLOT(res, Matrix_xSym, duplicate(x));
    } else { /* tv > tx */
	/* dMatrix result is only possibility until we have proper [iz]Matrix */
	if (tv < REALSXP) {
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
	    tv = REALSXP;
	}
	char *dcl = strdup(cl);
	dcl[0] = 'd';
	res = PROTECT(NEW_OBJECT_OF_CLASS(dcl)); ++nprotect;
	SET_SLOT(res, Matrix_xSym, coerceVector(x, tv));
	free(dcl);
	tx = tv;
    }
    x = PROTECT(GET_SLOT(res, Matrix_xSym)); ++nprotect;
    
    /* Transfer slots other than 'x', 'diag', and 'factors';
       latter two should keep their prototypes "N" and list()
    */
    SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
    SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
    SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
    
#define PM_D_S(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px = _PTR_(x);						\
	_CTYPE_ *pval = _PTR_(val);					\
	if (nv == 1) {							\
	    if (*uplo_P(res) == 'U') {					\
		for (int j = 0; j < n; px += (++j)+1)			\
		    *px = *pval;					\
	    } else {							\
		for (int j = 0; j < n; px += n-(j++))			\
		    *px = *pval;					\
	    }								\
	} else {							\
	    if (*uplo_P(res) == 'U') {					\
		for (int j = 0; j < n; px += (++j)+1)			\
		    *px = *(pval++);					\
	    } else {							\
		for (int j = 0; j < n; px += n-(j++))			\
		    *px = *(pval++);					\
	    }								\
	}								\
    } while (0)

    switch (tx) {
    case REALSXP:
	PM_D_S(double, REAL);
	break;
    case LGLSXP:
	PM_D_S(int, LOGICAL);
	break;
    case INTSXP:
	PM_D_S(int, INTEGER);
	break;
    case CPLXSXP:
	PM_D_S(Rcomplex, COMPLEX);
	break;
    default: /* should never happen */
	PM_ERROR_INVALID_SLOT_TYPE("x", tx, "diag_set");
	break;
    }
    
#undef PM_D_S
    
    UNPROTECT(nprotect);
    return res;
}

#define PM_XIJ_TR_UP_NUN(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ <= _J_							\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : _ZERO_)

#define PM_XIJ_TR_UP_UNT(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ < _J_							\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : (_I_ == _J_ ? _ONE_ : _ZERO_))

#define PM_XIJ_TR_LO_NUN(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : _ZERO_)

#define PM_XIJ_TR_LO_UNT(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ > _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : (_I_ == _J_ ? _ONE_ : _ZERO_))

#define PM_XIJ_SY_UP(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ <= _J_				       			\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))	       			\
     : *(_X_ + PM_AR21_UP(_J_, _I_)))

#define PM_XIJ_SY_LO(_X_, _I_, _J_, _N2_, _ZERO_, _ONE_)	\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _N2_))			\
     : *(_X_ + PM_AR21_LO(_J_, _I_, _N2_)))

#define PM_SUB1_LOOP(_XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	int i, j;							\
	if (TYPEOF(index) == INTSXP) {					\
	    int pos, *pindex = INTEGER(index);				\
	    for (R_xlen_t k = 0; k < nindex; ++k) {			\
		pos = *(pindex++);					\
		if (pos == NA_INTEGER || pos > nn) {			\
		    *(pres++) = _NA_;					\
		} else {						\
		    pos -= 1; /* 1-index -> 0-index */			\
		    i = pos % n;					\
		    j = pos / n;					\
		    *(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);	\
		}							\
	    }								\
	} else {							\
	    double pos, *pindex = REAL(index);				\
	    R_xlen_t truncpos;						\
	    for (R_xlen_t k = 0; k < nindex; ++k) {			\
		pos = *(pindex++);					\
		if (!R_FINITE(pos) || (truncpos = (R_xlen_t) pos) > nn) { \
		    *(pres++) = _NA_;					\
		} else {						\
		    truncpos -= 1; /* 1-index -> 0-index */		\
		    i = truncpos % n;					\
		    j = truncpos / n;					\
		    *(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);	\
		}							\
	    }								\
	}								\
    } while (0)

#define PM_SUB1_END(_CTYPE_, _SEXPTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)	\
    do {								\
	SEXP res = PROTECT(allocVector(_SEXPTYPE_, nindex));		\
	_CTYPE_ *pres = _PTR_(res), *px = _PTR_(x);			\
	if (*diag == ' ') { /* symmetric */				\
	    if (*uplo == 'U') {						\
		PM_SUB1_LOOP(PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB1_LOOP(PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else if (*diag == 'N') { /* non-unit triangular */		\
	    if (*uplo == 'U') {						\
		PM_SUB1_LOOP(PM_XIJ_TR_UP_NUN, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB1_LOOP(PM_XIJ_TR_LO_NUN, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else { /* unit triangular */					\
	    if (*uplo == 'U') {						\
		PM_SUB1_LOOP(PM_XIJ_TR_UP_UNT, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB1_LOOP(PM_XIJ_TR_LO_UNT, _NA_, _ZERO_, _ONE_);	\
	    }								\
	}    								\
	UNPROTECT(1);							\
	return res;							\
    } while (0)

#define PM_SUB1								\
    do {								\
	SEXP x = GET_SLOT(obj, Matrix_xSym);				\
	const char *uplo = uplo_P(obj), *diag = Diag_P(obj);		\
									\
	switch (TYPEOF(x)) {						\
	case REALSXP: /* d..Matrix */					\
	    PM_SUB1_END(double, REALSXP, REAL, NA_REAL, 0.0, 1.0);	\
	    break;							\
	case LGLSXP: /* [ln]..Matrix */					\
	    PM_SUB1_END(int, LGLSXP, LOGICAL, NA_LOGICAL, 0, 1);	\
	    break;							\
	case INTSXP: /* i..Matrix */					\
	    PM_SUB1_END(int, INTSXP, INTEGER, NA_INTEGER, 0, 1);	\
	    break;							\
	case CPLXSXP: /* z..Matrix */					\
	{								\
	    Rcomplex na, zero, one;					\
	    na.r = NA_REAL; zero.r = 0.0; one.r = 1.0;			\
	    na.i = NA_REAL; zero.i = 0.0; one.i = 0.0;			\
	    PM_SUB1_END(Rcomplex, CPLXSXP, COMPLEX, na, zero, one);	\
	    break;							\
	}								\
	default:							\
	    PM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x), "sub1");		\
	    break;							\
	}								\
	return R_NilValue;						\
    } while (0)
	  
/* 'x[i]' where 'i' is an integer or double vector with elements 
   greater than or equal to 1
*/
SEXP packedMatrix_sub1(SEXP obj, SEXP index)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t n2 = ((R_xlen_t) n) * 2, nn = ((R_xlen_t) n) * n,
	nindex = XLENGTH(index);
    PM_SUB1;
}

#undef PM_SUB1_LOOP

#define PM_SUB1_LOOP(_XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	int i, j, *pi = INTEGER(index), *pj = pi + nindex;		\
	for (int k = 0; k < nindex; ++k) {				\
	    i = *(pi++);						\
	    j = *(pj++);						\
	    if (i == NA_INTEGER || j == NA_INTEGER) {			\
		*(pres++) = _NA_;					\
	    } else {							\
		i -= 1; /* 1-index -> 0-index */			\
		j -= 1;							\
		*(pres++) = _XIJ_(px, i, j, n2, _ZERO_, _ONE_);		\
	    }								\
	}								\
    } while (0)

/* 'x[i]' where 'i' is a 2-column integer matrix supplying 
   integers in 'c(1:n, NA)' _only_
*/
SEXP packedMatrix_sub1_mat(SEXP obj, SEXP index)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0],
	nindex = INTEGER(getAttrib(index, R_DimSymbol))[0];
    R_xlen_t n2 = ((R_xlen_t) n) * 2;
    PM_SUB1;
}

#undef PM_SUB1
#undef PM_SUB1_LOOP

#define PM_SUB2_LOOP(_FOR_, _XIJ_, _NA_, _ZERO_, _ONE_)			\
    do {								\
        int i, j;							\
	for (int kj = 0; kj < nj; ++kj) {				\
	    if (mj) {							\
		j = kj;							\
	    } else {							\
		j = pj[kj];						\
		if (j == NA_INTEGER) {					\
		    _FOR_ {						\
			*(px1++) = _NA_;				\
		    }							\
		    if (do_cn) {					\
			SET_STRING_ELT(cn1, kj, NA_STRING);		\
		    }							\
		    continue;						\
		} else {						\
		    j -= 1;						\
		}							\
	    }								\
	    _FOR_ {							\
		if (mi) {						\
		    i = ki;						\
		} else {						\
		    i = pi[ki];						\
		    if (i == NA_INTEGER) {				\
			*(px1++) = _NA_;				\
			continue;					\
		    } else {						\
			i -= 1;						\
		    }							\
		}							\
		*(px1++) = _XIJ_(px0, i, j, n2, _ZERO_, _ONE_);		\
	    }								\
	    if (do_cn) {						\
		SET_STRING_ELT(cn1, kj, STRING_ELT(cn0, j));		\
	    }								\
	}								\
	if (do_rn) {							\
	    for (int ki = 0; ki < ni; ++ki) {				\
		if (mi) {						\
		    i = ki;						\
		} else {						\
		    i = pi[ki];						\
		    if (i == NA_INTEGER) {				\
			SET_STRING_ELT(rn1, ki, NA_STRING);		\
			continue;					\
		    } else {						\
			i -= 1;						\
		    }							\
		}							\
		SET_STRING_ELT(rn1, ki, STRING_ELT(rn0, i));		\
	    }								\
	}								\
    } while (0)

#define PM_SUB2(_CTYPE_, _SEXPTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)	\
    do {								\
	SEXP x1 = PROTECT(allocVector(_SEXPTYPE_, len));		\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	if (*diag == ' ') { /* symmetric */				\
	    if (*uplo == 'U') {						\
		if (do_sp) {						\
		    PM_SUB2_LOOP(for (int ki = 0; ki <= kj; ++ki),	\
				 PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
		} else {						\
		    PM_SUB2_LOOP(for (int ki = 0; ki < ni; ++ki),	\
				 PM_XIJ_SY_UP, _NA_, _ZERO_, _ONE_);	\
		}							\
	    } else {							\
		if (do_sp) {						\
		    PM_SUB2_LOOP(for (int ki = kj; ki < ni; ++ki),	\
				 PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
		} else { 						\
		    PM_SUB2_LOOP(for (int ki = 0; ki < ni; ++ki),	\
				 PM_XIJ_SY_LO, _NA_, _ZERO_, _ONE_);	\
		}							\
	    }								\
	} else if (*diag == 'N') { /* non-unit triangular */		\
	    if (*uplo == 'U') {						\
		PM_SUB2_LOOP(for (int ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_UP_NUN, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB2_LOOP(for (int ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_LO_NUN, _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else { /* unit triangular */					\
	    if (*uplo == 'U') {						\
		PM_SUB2_LOOP(for (int ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_UP_UNT, _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		PM_SUB2_LOOP(for (int ki = 0; ki < ni; ++ki),		\
			     PM_XIJ_TR_LO_UNT, _NA_, _ZERO_, _ONE_);	\
	    }								\
	}								\
	SET_SLOT(res, Matrix_xSym, x1);					\
	UNPROTECT(1);							\
    } while (0)

/* 'x[i, ]', 'x[, j]', and 'x[i, j]' where 'i' and 'j' are integer vectors 
   supplying integers in 'c(1:n, NA)' _only_ ... NULL indicates missingness
*/
SEXP packedMatrix_sub2(SEXP obj, SEXP index1, SEXP index2, SEXP drop)
{
    Rboolean mi = isNull(index1), mj = isNull(index2);
    int *pi, *pj, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t ni, nj, n2 = ((R_xlen_t) n) * 2;
    
    if (mi) {
	ni = n;
    } else {
	ni = XLENGTH(index1);
	if (ni > INT_MAX) {
	    error(_("dimensions of result cannot exceed 2^31-1"));
	}
	pi = INTEGER(index1);
    }
    if (mj) {
	nj = n;
    } else {
	nj = XLENGTH(index2);
	if (nj > INT_MAX) {
	    error(_("dimensions of result cannot exceed 2^31-1"));
	}
	pj = INTEGER(index2);
    }

    static const char *valid[] = {
	"dspMatrix", "lspMatrix", "nspMatrix",
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    int icl = R_check_class_etc(obj, valid), nprotect = 0;
    if (icl < 0) {
	PM_ERROR_INVALID_CLASS(class_P(obj), "sub2");
    }
    const char *uplo = uplo_P(obj), *diag = Diag_P(obj);

    /* Initialize result of same type but "general" class,
       except for symmetric indexing of symmetric matrix, 
       when class is retained also
    */
    SEXP res;
    Rboolean do_sp = (icl <= 2 && !mi && !mj && ni == nj &&
		      !memcmp(pi, pj, ni * sizeof(int)));
    if (do_sp) {
	res = PROTECT(NEW_OBJECT_OF_CLASS(valid[icl])); ++nprotect;
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
    } else {
	char *cl = strdup(valid[icl]);
	cl[1] = 'g';
	cl[2] = 'e';
	res = PROTECT(NEW_OBJECT_OF_CLASS(cl)); ++nprotect;
	free(cl);
    }

    /* Set 'Dim' slot */
    SEXP d1 = PROTECT(GET_SLOT(res, Matrix_DimSym));
    INTEGER(d1)[0] = ni;
    INTEGER(d1)[1] = nj;
    UNPROTECT(1);

    /* Set 'Dimnames' slot and 'names(Dimnames)' */
    SEXP dn0, dn1, rn0, rn1, cn0, cn1;
    dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    dn1 = PROTECT(GET_SLOT(res, Matrix_DimNamesSym)); ++nprotect;
    if (*diag == ' ') { /* symmetric */
    	int j;
	if (isNull(rn0 = cn0 = VECTOR_ELT(dn0, j = 1))) {
	    rn0 = cn0 = VECTOR_ELT(dn0, j = 0);
	}
	SEXP s;
	if (!isNull(s = getAttrib(dn0, R_NamesSymbol))) {
	    SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	    SET_STRING_ELT(ndn1, 0, s = STRING_ELT(s, j));
	    SET_STRING_ELT(ndn1, 1, s);
	    setAttrib(dn1, R_NamesSymbol, ndn1);
	    UNPROTECT(1);
	}
    } else { /* triangular */
	rn0 = VECTOR_ELT(dn0, 0);
	cn0 = VECTOR_ELT(dn0, 1);
	SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
	if (!isNull(ndn0)) {
	    setAttrib(dn1, R_NamesSymbol, ndn0);
	}
    }
    Rboolean has_rn, has_cn, do_rn, do_cn;
    has_rn = !isNull(rn0) && ni > 0;
    has_cn = !isNull(cn0) && nj > 0;
    do_rn = do_cn = FALSE;
    if (has_rn) {
	if (mi) {
	    SET_VECTOR_ELT(dn1, 0, rn0);
	} else {
	    rn1 = PROTECT(allocVector(STRSXP, ni)); ++nprotect;
	    SET_VECTOR_ELT(dn1, 0, rn1);
	    do_rn = TRUE;
	}
    }
    if (has_cn && !do_sp) {
	if (mj) {
	    SET_VECTOR_ELT(dn1, 1, cn0);
	} else {
	    cn1 = PROTECT(allocVector(STRSXP, nj)); ++nprotect;
	    SET_VECTOR_ELT(dn1, 1, cn1);
	    do_cn = TRUE;
	}
    }

    /* Set 'x' slot */
    SEXP x0 = GET_SLOT(obj, Matrix_xSym);
    R_xlen_t len = (do_sp ? (ni * (ni + 1)) / 2 : ni * nj);
    switch (TYPEOF(x0)) {
    case REALSXP: /* d..Matrix */
	PM_SUB2(double, REALSXP, REAL, NA_REAL, 0.0, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PM_SUB2(int, LGLSXP, LOGICAL, NA_LOGICAL, 0, 1);
	break;
    case INTSXP: /* i..Matrix */
	PM_SUB2(int, INTSXP, INTEGER, NA_INTEGER, 0, 1);
	break;
    case CPLXSXP: /* z..Matrix */
    {
	Rcomplex na, zero, one;
	na.r = NA_REAL; zero.r = 0.0; one.r = 1.0;
	na.i = NA_REAL; zero.i = 0.0; one.i = 0.0;
	PM_SUB2(Rcomplex, CPLXSXP, COMPLEX, na, zero, one);
	break;
    }
    default:
	PM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x0), "sub2");
    }

    /* Drop dimensions in this special case */
    if (asLogical(drop) != 0 && (ni == 1 || nj == 1)) {
	res = PROTECT(GET_SLOT(res, Matrix_xSym)); ++nprotect;
	if (has_rn && nj == 1 && ni != 1) {
	    setAttrib(res, R_NamesSymbol, VECTOR_ELT(dn1, 0));
	} else if (has_cn && ni == 1 && nj != 1) {
	    setAttrib(res, R_NamesSymbol, VECTOR_ELT(dn1, 1));
	}
    }
    UNPROTECT(nprotect);
    return res;
}

#undef PM_SUB2
#undef PM_SUB2_LOOP
