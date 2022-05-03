#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "unpackedMatrix.h"

#undef HAVE_PROPER_IMATRIX 
#undef HAVE_PROPER_ZMATRIX

/* pack(x) */
SEXP unpackedMatrix_pack(SEXP from, SEXP tr_if_ge, SEXP up_if_ge)
{
    static const char *valid_from[] = {
	"corMatrix", "dpoMatrix", /* must match before "dsyMatrix" */
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"Cholesky", "BunchKaufman", /* must match before "dtrMatrix" */
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    static const char *valid_to[] = {
	"dppMatrix", "dppMatrix", /* no pcorMatrix _yet_ */
	"dspMatrix", "lspMatrix", "nspMatrix",
	"pCholesky", "pBunchKaufman",
	"dtpMatrix", "ltpMatrix", "ntpMatrix", ""};
    int icl = R_check_class_etc(from, valid_from);
    if (icl < 0) {
	UPM_ERROR_INVALID_CLASS(class_P(from), "pack");
    }

    Rboolean ge = (icl > 9);
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int n = INTEGER(dim)[0];
    if (ge) {
	if (INTEGER(dim)[1] != n) {
	    error(_("attempt to pack non-square matrix"));
	}
	icl -= (asLogical(tr_if_ge)) ? 3 : 3+2+3; /* ? ge->tp : ge->sp */
    }
    
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid_to[icl]));
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (ge) {
	SET_SLOT(to, Matrix_uploSym, mkString(asLogical(up_if_ge) ? "U" : "L"));
    } else {
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
	if (icl <= 4) {
	    /* .syMatrix */
	    SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
	} else {
	    /* .trMatrix */
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	    if (icl == 6) {
		/* BunchKaufman */
		SET_SLOT(to, Matrix_permSym, GET_SLOT(from, Matrix_permSym));
	    }
	}
    }

    SEXP x_from = GET_SLOT(from, Matrix_xSym);

#define UPM_PACK(_PREFIX_, _CTYPE_, _SEXPTYPE_, _PTR_)			\
    do {								\
	R_xlen_t len = (n * ((R_xlen_t) n + 1)) / 2;			\
	SEXP x_to = PROTECT(allocVector(_SEXPTYPE_, len));		\
	_PREFIX_ ## dense_pack(_PTR_(x_to), _PTR_(x_from), n,		\
			       (*uplo_P(to) == 'U') ? UPP : LOW,	\
			       NUN /* for speed */);			\
	SET_SLOT(to, Matrix_xSym, x_to);				\
	UNPROTECT(1);							\
    } while (0)
    
    switch (TYPEOF(x_from)) {
    case REALSXP: /* d..Matrix */
	UPM_PACK(d, double, REALSXP, REAL);
	break;
    case LGLSXP: /* [ln]..Matrix */
	UPM_PACK(i, int, LGLSXP, LOGICAL);
	break;
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP: /* i..Matrix */
	UPM_PACK(i, int, INTSXP, INTEGER);
	break;
#endif
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP: /* z..Matrix */
    	UPM_PACK(z, Rcomplex, CPLXSXP, COMPLEX);
	break;
#endif
    default:
	UPM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x_from), "pack");
	break;
    }

#undef UPM_PACK
    
    UNPROTECT(1);
    return to;
}

SEXP matrix_pack(SEXP from, SEXP tr, SEXP up)
{
    SEXP dim = getAttrib(from, R_DimSymbol);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n) {
	error(_("attempt to pack non-square matrix"));
    }
    char cl[] = "..pMatrix";
    int nprotect = 0;
    
    switch (TYPEOF(from)) {
    case LGLSXP:
	cl[0] = 'l';
	break;
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
	cl[0] = 'i';
	break;
#else
    case INTSXP:
	from = PROTECT(coerceVector(from, REALSXP)); ++nprotect;
#endif
    case REALSXP:
	cl[0] = 'd';
	break;
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	cl[0] = 'z';
	break;
#endif
    default:
	error(_("matrix of invalid type \"%s\" in 'matrix_pack()'"),
	      type2char(TYPEOF(from)));
	break;
    }
    cl[1] = asLogical(tr) ? 't' : 's';
    Rprintf("%s\n", cl);
    
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dimnames = getAttrib(from, R_DimNamesSymbol); ++nprotect;
    SET_SLOT(to, Matrix_DimSym, dim);
    if (!isNull(dimnames)) {
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    }
    SET_SLOT(to, Matrix_uploSym, mkString(asLogical(up) ? "U" : "L"));
    
#define M_PACK(_PREFIX_, _CTYPE_, _SEXPTYPE_, _PTR_)			\
    do {								\
	R_xlen_t len = (n * ((R_xlen_t) n + 1)) / 2;			\
	SEXP x = PROTECT(allocVector(_SEXPTYPE_, len));			\
	_PREFIX_ ## dense_pack(_PTR_(x), _PTR_(from), n,		\
			       (*uplo_P(to) == 'U') ? UPP : LOW,	\
			       NUN);					\
	SET_SLOT(to, Matrix_xSym, x);					\
	UNPROTECT(1);							\
    } while (0)

    switch (TYPEOF(from)) {
    case LGLSXP:
	M_PACK(i, int, LGLSXP, LOGICAL);
	break;
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
	M_PACK(i, int, INTSXP, INTEGER);
	break;
#endif
    case REALSXP:
	M_PACK(d, double, REALSXP, REAL);
	break;
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	M_PACK(z, Rcomplex, CPLXSXP, COMPLEX);
	break;
#endif
    default:
	error(_("matrix of invalid type \"%s\" in 'matrix_pack()'"),
	      type2char(TYPEOF(from)));
	break;
    }

#undef M_PACK
    
    UNPROTECT(nprotect);
    return to;
}

#define UPM_IS_SY(_RES_, _X_, _N_, _METHOD_)				\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_symmetric(REAL(_X_), _N_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = ldense_unpacked_is_symmetric(LOGICAL(_X_), _N_);	\
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_symmetric(INTEGER(_X_), _N_);	\
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_symmetric(COMPLEX(_X_), _N_);	\
	    break;							\
	default:							\
	    UPM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(_X_), _METHOD_);	\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

#define UPM_IS_TR(_RES_, _X_, _N_, _UPPER_, _METHOD_)			\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ =							\
		ddense_unpacked_is_triangular(REAL(_X_), _N_, _UPPER_); \
	    break;							\
	case LGLSXP:							\
	    _RES_ =							\
		idense_unpacked_is_triangular(LOGICAL(_X_), _N_, _UPPER_); \
	    break;							\
	case INTSXP:							\
	    _RES_ =							\
		idense_unpacked_is_triangular(INTEGER(_X_), _N_, _UPPER_); \
	    break;							\
	case CPLXSXP:							\
	    _RES_ =							\
		zdense_unpacked_is_triangular(COMPLEX(_X_), _N_, _UPPER_); \
	    break;							\
	default:							\
	    UPM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(_X_), _METHOD_);	\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

#define UPM_IS_DI(_RES_, _X_, _N_, _METHOD_)				\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_diagonal(REAL(_X_), _N_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_unpacked_is_diagonal(LOGICAL(_X_), _N_);	\
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_diagonal(INTEGER(_X_), _N_);	\
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_diagonal(COMPLEX(_X_), _N_);	\
	    break;							\
	default:							\
	    UPM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(_X_), _METHOD_);	\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

/* isSymmetric(x, tol = 0) */
SEXP unpackedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    static const char *valid[] = {
	"dsyMatrix", "lsyMatrix", "nsyMatrix", /* be fast */
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	UPM_ERROR_INVALID_CLASS(class_P(obj), "is_symmetric");
	return R_NilValue;
    } else if (icl <= 2) {
	/* .syMatrix: symmetric by definition */
	return ScalarLogical(1);
    } else {
	/* .(ge|tr)Matrix */
	int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
	if (icl > 5 && pdim[1] != n) {
	    return ScalarLogical(0);
	}
	if (asLogical(checkDN) != 0 &&
	    !DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym))) {
	    return ScalarLogical(0);
	}
	Rboolean res = FALSE;
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (icl <= 5) {
	    /* .trMatrix: symmetric iff diagonal (upper _and_ lower tri.) */
	    Rboolean upper = (*uplo_P(obj) == 'U') ? FALSE : TRUE;
	    UPM_IS_TR(res, x, n, upper, "is_symmetric");
	} else {
	    /* .geMatrix: need to do a complete symmetry check */
	    UPM_IS_SY(res, x, n, "is_symmetric");
	}
	return ScalarLogical(res);
    }
}

/* isSymmetric(x, tol = 0) */
SEXP matrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n) {
	return ScalarLogical(0);
    }
    if (asLogical(checkDN) != 0) {
	SEXP dn = getAttrib(obj, R_DimNamesSymbol);
	if (!isNull(dn) && !DimNames_is_symmetric(dn)) {
	    return ScalarLogical(0);
	}
    }
    Rboolean res = FALSE;
    UPM_IS_SY(res, obj, n, "is_symmetric");
    /* ^FIXME: wrong function name in error message */
    return ScalarLogical(res);
}

#define RETURN_TRUE(_KIND_)						\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

#define RETURN_GE_IS_TR(_X_, _N_, _UPPER_)				\
    do {								\
	Rboolean res = FALSE;						\
	if (_UPPER_ == NA_LOGICAL) {					\
	    UPM_IS_TR(res, _X_, _N_, TRUE, "is_triangular");		\
	    if (res) {							\
	        RETURN_TRUE(mkString("U"));				\
	    }								\
	    UPM_IS_TR(res, _X_, _N_, FALSE, "is_triangular");		\
	    if (res) {							\
		RETURN_TRUE(mkString("L"));				\
	    }								\
	} else {							\
	    if (_UPPER_ != 0) {						\
		UPM_IS_TR(res, _X_, _N_, TRUE, "is_triangular");	\
	    } else { 							\
		UPM_IS_TR(res, _X_, _N_, FALSE, "is_triangular");	\
	    }								\
	    if (res) {							\
		return ScalarLogical(1);				\
	    }								\
	}								\
	return ScalarLogical(0);					\
    } while (0)

/* isTriangular(x, upper) */
SEXP unpackedMatrix_is_triangular(SEXP obj, SEXP upper)
{
    static const char *valid[] = {
	"dtrMatrix", "ltrMatrix", "ntrMatrix", /* be fast */
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	UPM_ERROR_INVALID_CLASS(class_P(obj), "is_triangular");
    }
    
    SEXP uplo;
    Rboolean have_upper;
    int need_upper = asLogical(upper);
    if (icl <= 5) {
	uplo = GET_SLOT(obj, Matrix_uploSym);
	have_upper = (*CHAR(STRING_ELT(uplo, 0)) == 'U') ? TRUE : FALSE;
    }

#define IF_DIAGONAL							\
    Rboolean res = FALSE;						\
    SEXP x = GET_SLOT(obj, Matrix_xSym);				\
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];			\
    UPM_IS_TR(res, x, n, have_upper ? FALSE : TRUE,			\
	      "is_triangular");						\
    if (res)
    
    if (icl <= 2) {
	/* .trMatrix: be fast if 'upper', 'uplo' agree; else need diagonal */
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
	return ScalarLogical(0);
    } else if (icl <= 5) {
	/* .syMatrix: triangular iff diagonal (upper _and_ lower tri.) */
	IF_DIAGONAL {
	    if (need_upper == NA_LOGICAL) {
		RETURN_TRUE(mkString("U"));
	    } else {
		return ScalarLogical(1);
	    }
	}
	return ScalarLogical(0);

#undef IF_DIAGONAL
	
    } else {
	/* .geMatrix: need to do a complete triangularity check */
	int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
	if (pdim[1] != n) {
	    return ScalarLogical(0);
	}
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	RETURN_GE_IS_TR(x, n, need_upper);
    }
}

/* isTriangular(x, upper) */
SEXP matrix_is_triangular(SEXP obj, SEXP upper)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n) {
	return ScalarLogical(0);
    }
    int need_upper = asLogical(upper);
    RETURN_GE_IS_TR(obj, n, need_upper);
    /* ^FIXME: wrong function name in error message */
}

#undef RETURN_GE_IS_TR
#undef RETURN_TRUE

/* isDiagonal(x) */
SEXP unpackedMatrix_is_diagonal(SEXP obj)
{
    static const char *valid[] = {
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	UPM_ERROR_INVALID_CLASS(class_P(obj), "is_diagonal");
    }
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n) {
	return ScalarLogical(0);
    }
    Rboolean res = FALSE;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (icl <= 5) {
	/* .(sy|tr)Matrix: diagonal iff stored triangle is zero off diagonal */
	Rboolean upper = (*uplo_P(obj) == 'U') ? FALSE : TRUE;
	UPM_IS_TR(res, x, n, upper, "is_diagonal");
    } else {
	/* .geMatrix: need to do a complete diagonality check */
	UPM_IS_DI(res, x, n, "is_diagonal");
    }
    return ScalarLogical(res);
}

/* isDiagonal(x) */
SEXP matrix_is_diagonal(SEXP obj)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n) {
	return ScalarLogical(0);
    }
    Rboolean res = FALSE;
    UPM_IS_DI(res, obj, n, "is_diagonal");
    /* ^FIXME: wrong function name in error message */
    return ScalarLogical(res);
}

#undef UPM_IS_DI
#undef UPM_IS_TR
#undef UPM_IS_SY

/* t(x) */
/* MJ: Is NO_REFERENCES(x@x) ever true?? (in which case do in-place??) */
/* MJ: Technically no need to do full transpose of .(sy|tr)Matrix ...  */
/*     but then identical(x, t(t(x))) can be FALSE ...                 */
SEXP unpackedMatrix_t(SEXP obj)
{
    /* Initialize result of same class, except for BunchKaufman */
    const char *cl = class_P(obj);
    SEXP res, x0 = GET_SLOT(obj, Matrix_xSym);
    if (TYPEOF(x0) == REALSXP) {
	static const char *valid[] = { "BunchKaufman", "" };
	if (!R_check_class_etc(obj, valid)) {
	    cl = "dtrMatrix";
	}
    }
    res = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    int *pdim0 = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim0[0], n = pdim0[1];
    R_xlen_t len = XLENGTH(x0);
    const char *uplo = Uplo_P(obj), *diag = Diag_P(obj);
    
#define UPM_T(_CTYPE_, _SEXPTYPE_, _PTR_)				\
    do {								\
	SEXP x1 = PROTECT(allocVector(_SEXPTYPE_, len));		\
	_CTYPE_ *px0 = _PTR_(x0);					\
	_CTYPE_ *px1 = _PTR_(x1);					\
	int i, j;							\
	/* if (*uplo == ' ') { */					\
	R_xlen_t len1s = len - 1;					\
	for (j = 0; j < m; ++j, px0 -= len1s)				\
	    for (i = 0; i < n; ++i, px0 += m)				\
		*(px1++) = *px0;					\
	/* } else { */							\
	/*     Memzero(px1, n * (size_t) n); */				\
	/*     R_xlen_t upos = 0, lpos = 0; */				\
	/*     if (*uplo == 'U') */					\
	/* 	for (j = 0; j < n; upos = (lpos += (++j))) */		\
	/* 	    for (i = j; i < n; ++i, ++lpos, upos += n) */	\
	/* 		px1[lpos] = px0[upos]; */			\
	/*     else */							\
	/* 	for (j = 0; j < n; upos = (lpos += (++j))) */		\
	/* 	    for (i = j; i < n; ++i, ++lpos, upos += n) */	\
	/* 		px1[upos] = px0[lpos]; */			\
	/* } */								\
	SET_SLOT(res, Matrix_xSym, x1);					\
	UNPROTECT(1);							\
    } while (0)

    if (len > 1) {
	/* Permute 'x' slot */
	switch (TYPEOF(x0)) {
	case REALSXP: /* d..Matrix */
	    UPM_T(double, REALSXP, REAL);
	    break;
	case LGLSXP: /* [ln]..Matrix */
	    UPM_T(int, LGLSXP, LOGICAL);
	    break;
	case INTSXP: /* i..Matrix */
	    UPM_T(int, INTSXP, INTEGER);
	    break;
	case CPLXSXP: /* z..Matrix */
	    UPM_T(Rcomplex, CPLXSXP, COMPLEX);
	    break;
	default:
	    UPM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x0), "t");
	    break;
	}
    } else {
	/* Preserve 'x' slot */
	SET_SLOT(res, Matrix_xSym, x0);
    }

#undef UPM_T

    /* Preserve or reverse 'Dim' slot (preserving if square) */
    if (m == n) {
	SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
    } else {
	SEXP dim1 = PROTECT(allocVector(INTSXP, 2));
	int *pdim1 = INTEGER(dim1);
	pdim1[0] = n;
	pdim1[1] = m;
	SET_SLOT(res, Matrix_DimSym, dim1);
	UNPROTECT(1);
    }
    /* Preserve or reverse 'Dimnames' slot and 'names(Dimnames)' 
       (preserving if symmetric)
    */
    if (*uplo != ' ' && *diag == ' ') {
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
    if (*uplo != ' ') {
	/* Toggle 'uplo' slot */
	SET_SLOT(res, Matrix_uploSym, mkString(*uplo == 'U' ? "L" : "U"));
	if (*diag == 'U') {
	    /* Preserve 'diag' slot */
	    SET_SLOT(res, Matrix_diagSym, GET_SLOT(obj, Matrix_diagSym));
	} else if (*diag == ' ') {
	    static const char *valid[] = { "corMatrix", "" };
	    if (!R_check_class_etc(obj, valid)) {
		/* Preserve 'sd' slot */
		SET_SLOT(res, install("sd"), GET_SLOT(obj, install("sd")));
	    }
	}
    }
    /* NB: Nothing to do for 'factors' slot: prototype is already list() ...
       FIXME: However, it would be much better to transpose and reverse
       the factors in each factorization
    */

    UNPROTECT(1);
    return res;
}

/* diag(x, names) */
SEXP unpackedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL) {
	error(_("'names' must be TRUE or FALSE"));
    }

    SEXP res = R_NilValue, x = GET_SLOT(obj, Matrix_xSym);
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    const char *uplo = Uplo_P(obj), *diag = Diag_P(obj);
    
#define UPM_D_G(_CTYPE_, _SEXPTYPE_, _PTR_, _INIT_ONE_, _ONE_)		\
    do {								\
	res = PROTECT(allocVector(_SEXPTYPE_, r));			\
	_CTYPE_ *pres = _PTR_(res);					\
	if (*diag == 'U') {						\
	    _INIT_ONE_;							\
	    for (int j = 0; j < r; ++j)					\
		*(pres++) = _ONE_;					\
	} else {							\
	    _CTYPE_ *px = _PTR_(x);					\
	    R_xlen_t mp1 = ((R_xlen_t) m) + 1;				\
	    for (int j = 0; j < r; ++j, px += mp1)			\
		*(pres++) = *px;					\
	}								\
    } while (0)
    
    switch (TYPEOF(x)) {
    case REALSXP: /* d..Matrix */
	UPM_D_G(double, REALSXP, REAL, , 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	UPM_D_G(int, LGLSXP, LOGICAL, , 1);
	break;
    case INTSXP: /* i..Matrix */
	UPM_D_G(int, INTSXP, INTEGER, , 1);
	break;
    case CPLXSXP: /* z..Matrix */
	UPM_D_G(Rcomplex, CPLXSXP, COMPLEX,
		Rcomplex one; one.r = 1.0; one.i = 0.0;, one);
	break;
    default:
	UPM_ERROR_INVALID_SLOT_TYPE("x", TYPEOF(x), "diag_get");
	break;
    }

#undef UPM_D_G

    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames'
	*/
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	Rboolean symmetric = (*uplo != ' ' && *diag == ' ') ? TRUE : FALSE;
	if (isNull(cn)) {
	    if (symmetric && !isNull(rn)) {
		setAttrib(res, R_NamesSymbol, rn);
	    }
	} else {
	    if (symmetric) {
		setAttrib(res, R_NamesSymbol, cn);
	    } else if (!isNull(rn) &&
		       (rn == cn || equal_string_vectors(rn, cn, r))) {
		setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
	    }
	}
    }
    UNPROTECT(1);
    return res;
}

/* diag(x) <- value */
/* MJ: Is NO_REFERENCES(x@x) ever true?? (in which case do in-place??) */
SEXP unpackedMatrix_diag_set(SEXP obj, SEXP val)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim[0], n = pdim[1], r = (m < n) ? m : n, nv = LENGTH(val);
    if (nv != 1 && nv != r) {
	error(_("replacement diagonal has wrong length"));
    }

    SEXP x = GET_SLOT(obj, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x), tv = TYPEOF(val);
    if (tx < LGLSXP || tx > CPLXSXP) {
	UPM_ERROR_INVALID_SLOT_TYPE("x", tx, "diag_set");
    }
    if (tv < LGLSXP || tv > REALSXP) {
	/* Upper bound can become CPLXSXP once we have proper zdenseMatrix */
	error(_("replacement diagonal has incompatible type \"%s\""),
	      type2char(tv));
    }
    
    static const char *valid[] = {
	"dgeMatrix", "dsyMatrix", "dtrMatrix",
	"lgeMatrix", "lsyMatrix", "ltrMatrix",
	"ngeMatrix", "nsyMatrix", "ntrMatrix", ""};
    int icl = R_check_class_etc(obj, valid);
    if (icl < 0) {
	UPM_ERROR_INVALID_CLASS(class_P(obj), "diag_set");
    }
    const char *cl = valid[icl];

    SEXP res = obj;
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
    if (R_has_slot(res, Matrix_uploSym)) {
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
    }

#define UPM_D_S(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px = _PTR_(x);						\
	_CTYPE_ *pval = _PTR_(val);					\
	R_xlen_t mp1 = ((R_xlen_t) m) + 1;				\
	if (nv == 1)							\
	    for (int j = 0; j < r; ++j, px += mp1)			\
		*px = *pval;						\
	else								\
	    for (int j = 0; j < r; ++j, px += mp1)			\
		*px = *(pval++);					\
    } while (0)
    
    switch (tx) {
    case REALSXP:
	UPM_D_S(double, REAL);
	break;
    case LGLSXP:
	UPM_D_S(int, LOGICAL);
	break;
    case INTSXP:
	UPM_D_S(int, INTEGER);
	break;
    case CPLXSXP:
	UPM_D_S(Rcomplex, COMPLEX);
	break;
    default: /* should never happen */
	UPM_ERROR_INVALID_SLOT_TYPE("x", tx, "diag_set");
	break;
    }

#undef UPM_D_S
    
    UNPROTECT(nprotect);
    return res;
}
