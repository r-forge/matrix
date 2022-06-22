#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "unpackedMatrix.h"

/* pack(x), returning packedMatrix */
SEXP unpackedMatrix_pack(SEXP from, SEXP strict, SEXP tr_if_ge, SEXP up_if_ge)
{
    static const char *valid_from[] = {
	/*  0 */ "Cholesky", "BunchKaufman", /* must match before dtrMatrix */
	/*  2 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/*  5 */ "corMatrix", "dpoMatrix", /* must match before dsyMatrix */
	/*  7 */ "dsyMatrix", "lsyMatrix", "nsyMatrix",
	/* 10 */ "dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    static const char *valid_to[] = {
	/*  0 */ "pCholesky", "pBunchKaufman",
	/*  2 */ "dtpMatrix", "ltpMatrix", "ntpMatrix",
	/*  5 */ "dppMatrix", "dppMatrix", /* no pcorMatrix _yet_ */
	/*  7 */ "dspMatrix", "lspMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid_from);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_pack");
    if (asLogical(strict) == 0) {
	if (ivalid < 2)
	    ivalid = 2; /* Cholesky,BunchKaufman->dtpMatrix */
	else if (ivalid == 5 || ivalid == 6)
	    ivalid = 7; /* corMatrix,dpoMatrix->dspMatrix */
    }

    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), n = pdim[0], shift = 0;
    if (ivalid >= 10) {
	/* .geMatrix */
	if (pdim[1] != n)
	    error(_("attempt to pack non-square matrix"));
	shift = (asLogical(tr_if_ge) != 0) ?  3+2+3 :      3;
	/*                                 ? ge->tp : ge->sp */
    }

    SEXPTYPE tx;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(valid_to[ivalid - shift])),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x_from = GET_SLOT(from, Matrix_xSym),
	x_to = PROTECT(allocVector(tx = TYPEOF(x_from), PM_LENGTH(n))),
	uplo;

    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    SET_SLOT(to, Matrix_xSym, x_to);
    if (ivalid >= 10) {
	/* .geMatrix */
	uplo = mkString((asLogical(up_if_ge) != 0) ? "U" : "L");
    } else {
	/* .(sy|tr)Matrix */
	uplo = GET_SLOT(from, Matrix_uploSym);
	if (ivalid < 5) {
	    /* .trMatrix */
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	    if (ivalid == 1)
		/* BunchKaufman */
		SET_SLOT(to, Matrix_permSym, GET_SLOT(from, Matrix_permSym));
	} else {
	    /* .syMatrix */
	    SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
	}
    }
    SET_SLOT(to, Matrix_uploSym, uplo);
    char ul = *CHAR(STRING_ELT(uplo, 0));

#define PACK(_PREFIX_, _PTR_)						\
    _PREFIX_ ## dense_pack(_PTR_(x_to), _PTR_(x_from), n, ul, 'N')
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	PACK(d, REAL);
	break;
    case LGLSXP: /* [ln]..Matrix */
	PACK(i, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	PACK(i, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
    	PACK(z, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_pack");
	break;
    }

#undef PACK

    UNPROTECT(2);
    return to;
}

/* pack(x), returning packedMatrix */
SEXP matrix_pack(SEXP from, SEXP tr, SEXP up)
{
    SEXP dim = getAttrib(from, R_DimSymbol);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to pack a non-square matrix"));
    SEXPTYPE tf = TYPEOF(from);
    char cl[] = "..pMatrix";
    int nprotect = 0;
    
    switch (tf) {
    case LGLSXP:
	cl[0] = 'l';
	break;
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
	cl[0] = 'i';
	break;
#else
    case INTSXP:
	PROTECT(from = coerceVector(from, tf = REALSXP));
	++nprotect;
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
	ERROR_INVALID_TYPE("matrix", tf, "matrix_pack");
	break;
    }
    cl[1] = (asLogical(tr) != 0) ? 't' : 's';
    
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dimnames = getAttrib(from, R_DimNamesSymbol),
	x = PROTECT(allocVector(tf, PM_LENGTH(n))),
	uplo = mkString(asLogical(up) ? "U" : "L");
    char ul = *CHAR(STRING_ELT(uplo, 0));
    nprotect += 2;

    SET_SLOT(to, Matrix_DimSym, dim);
    if (!isNull(dimnames))
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    SET_SLOT(to, Matrix_xSym, x);
    SET_SLOT(to, Matrix_uploSym, uplo);

#define PACK(_PREFIX_, _PTR_)					\
    _PREFIX_ ## dense_pack(_PTR_(x), _PTR_(from), n, ul, 'N')
    
    switch (tf) {
    case LGLSXP:
	PACK(i, LOGICAL);
	break;
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
	PACK(i, INTEGER);
	break;
#endif
    case REALSXP:
	PACK(d, REAL);
	break;
#ifdef HAVE_PROPER_ZMATRIX
    case CPLXSXP:
	PACK(z, COMPLEX);
	break;
#endif
    default:
	ERROR_INVALID_TYPE("matrix", tf, "matrix_pack");
	break;
    }

#undef PACK

    UNPROTECT(nprotect);
    return to;
}

/* forceSymmetric(x, uplo), returning .syMatrix */
SEXP unpackedMatrix_force_symmetric(SEXP from, SEXP uplo_to)
{
    static const char *valid[] = {
	"dsyMatrix", "lsyMatrix", "nsyMatrix", /* be fast */
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_force_symmetric");
    const char *clf = valid[ivalid];
    
    SEXP uplo_from;
    char ulf = 'U', ult = *CHAR(asChar(uplo_to));
    if (clf[1] != 'g') {
	/* .(sy|tr)Matrix */
	uplo_from = GET_SLOT(from, Matrix_uploSym);
	ulf = *CHAR(STRING_ELT(uplo_from, 0));
    }
    if (ult == '\0') /* to handle missing(uplo) */
	ult = ulf;
    if (clf[1] == 's') {
	/* .syMatrix */
	if (ulf == ult)
	    return from;
	SEXP to = PROTECT(unpackedMatrix_transpose(from));
	if (clf[0] == 'z') {
	    /* Need _conjugate_ transpose */
	    SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	    conjugate(x);
	    UNPROTECT(1);
	}
	UNPROTECT(1);
	return to;
    }
    
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to symmetrize a non-square matrix"));

    /* Now handling just square .(tr|ge)Matrix ... */
    
    char *clt = strdup(clf);
    clt[1] = 's';
    clt[2] = 'y';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x_from = GET_SLOT(from, Matrix_xSym);
    free(clt);
    
    if (clf[1] == 'g' || ulf == ult) {
	/* .geMatrix or .trMatrix with correct uplo */
	SET_SLOT(to, Matrix_xSym, x_from);
    } else {
	/* .trMatrix with incorrect uplo */
	char di = *diag_P(from);
	SEXPTYPE tx = TYPEOF(x_from);
	R_xlen_t nx = XLENGTH(x_from);
	SEXP x_to = PROTECT(allocVector(tx, nx));

#define COPY_DIAGONAL(_PREFIX_, _PTR_)					\
	do {								\
	    Memzero(_PTR_(x_to), nx);					\
	    _PREFIX_ ## dense_unpacked_copy_diagonal(			\
		_PTR_(x_to), _PTR_(x_from), n, nx, 'U' /* unused */, di); \
	} while (0)
	
	switch (tx) {
	case REALSXP: /* d..Matrix */
	    COPY_DIAGONAL(d, REAL);
	    break;
	case LGLSXP: /* [ln]..Matrix */
	    COPY_DIAGONAL(i, LOGICAL);
	    break;
	case INTSXP: /* i..Matrix */
	    COPY_DIAGONAL(i, INTEGER);
	    break;
	case CPLXSXP: /* z..Matrix */
	    COPY_DIAGONAL(z, COMPLEX);
	    break;
	default:
	   ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_force_symmetric");
	    break;
	}

#undef COPY_DIAGONAL
       
	SET_SLOT(to, Matrix_xSym, x_to);
	UNPROTECT(1);
    }

    SET_SLOT(to, Matrix_DimSym, dim);
    set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_uploSym, mkString((ult == 'U') ? "U" : "L"));    
    UNPROTECT(1);
    return to;
}

/* forceSymmetric(x, uplo), returning .syMatrix */
SEXP matrix_force_symmetric(SEXP from, SEXP uplo_to) {
    SEXP dim = getAttrib(from, R_DimSymbol);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to symmetrize a non-square matrix"));
    SEXPTYPE tf = TYPEOF(from);
    char cl[] = ".syMatrix", ult = *CHAR(asChar(uplo_to));
    int nprotect = 0;

    switch (tf) {
    case LGLSXP:
	cl[0] = 'l';
	break;
#ifdef HAVE_PROPER_IMATRIX
    case INTSXP:
	cl[0] = 'i';
	break;
#else
    case INTSXP:
	PROTECT(from = coerceVector(from, tf = REALSXP));
	++nprotect;
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
	ERROR_INVALID_TYPE("matrix", tf, "matrix_force_symmetric");
	break;
    }

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dimnames = getAttrib(from, R_DimNamesSymbol);
    ++nprotect;
    
    if (MAYBE_REFERENCED(from)) {
	R_xlen_t nn = (R_xlen_t) n * n;
	SEXP x = PROTECT(allocVector(tf, nn));
	++nprotect;
	switch (tf) {
	case LGLSXP:
	    Memcpy(LOGICAL(x), LOGICAL(from), nn);
	    break;
#ifdef HAVE_PROPER_IMATRIX
	case INTSXP:
	    Memcpy(INTEGER(x), INTEGER(from), nn);
	    break;
#endif
	case REALSXP:
	    Memcpy(REAL(x), REAL(from), nn);
	    break;
#ifdef HAVE_PROPER_ZMATRIX
	case CPLXSXP:
	    Memcpy(COMPLEX(x), COMPLEX(from), nn);
	    break;
#endif
	default:
	    ERROR_INVALID_TYPE("matrix", tf, "matrix_force_symmetric");
	    break;
	}
	SET_SLOT(to, Matrix_xSym, x);
    } else {
	SET_ATTRIB(from, R_NilValue);
	if (OBJECT(from))
	    SET_OBJECT(from, 0);
	SET_SLOT(to, Matrix_xSym, from);
    }

    SET_SLOT(to, Matrix_DimSym, dim);
    if (!isNull(dimnames))
	set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_uploSym, mkString((ult == 'U') ? "U" : "L"));
    UNPROTECT(nprotect);
    return to;
}

#define UPM_IS_SY(_RES_, _X_, _N_, _LDENSE_, _WHAT_, _METHOD_)		\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_symmetric(REAL(_X_), _N_);	\
	    break;							\
	case LGLSXP:							\
	    _RES_ = (_LDENSE_						\
		     ? ldense_unpacked_is_symmetric(LOGICAL(_X_), _N_)	\
		     : ndense_unpacked_is_symmetric(LOGICAL(_X_), _N_)); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_symmetric(INTEGER(_X_), _N_);	\
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_symmetric(COMPLEX(_X_), _N_);	\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE(_WHAT_, TYPEOF(_X_), _METHOD_);		\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

#define UPM_IS_TR(_RES_, _X_, _N_, _UPLO_, _WHAT_, _METHOD_)		\
    do {								\
	switch (TYPEOF(_X_)) {						\
	case REALSXP:							\
	    _RES_ = ddense_unpacked_is_triangular(REAL(_X_), _N_, _UPLO_); \
	    break;							\
	case LGLSXP:							\
	    _RES_ = idense_unpacked_is_triangular(LOGICAL(_X_), _N_, _UPLO_); \
	    break;							\
	case INTSXP:							\
	    _RES_ = idense_unpacked_is_triangular(INTEGER(_X_), _N_, _UPLO_); \
	    break;							\
	case CPLXSXP:							\
	    _RES_ = zdense_unpacked_is_triangular(COMPLEX(_X_), _N_, _UPLO_); \
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE(_WHAT_, TYPEOF(_X_), _METHOD_);		\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

#define UPM_IS_DI(_RES_, _X_, _N_, _WHAT_, _METHOD_)			\
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
	    ERROR_INVALID_TYPE(_WHAT_, TYPEOF(_X_), _METHOD_);		\
	    _RES_ = FALSE;						\
	    break;							\
	}								\
    } while (0)

/* isSymmetric(x, tol = 0, checkDN) */
SEXP unpackedMatrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    static const char *valid[] = {
	"dsyMatrix", "lsyMatrix", "nsyMatrix", /* be fast */
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0) {
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_is_symmetric");
	return R_NilValue;
    } else if (ivalid < 3) {
	/* .syMatrix: symmetric by definition */
	return ScalarLogical(1);
    } else {
	/* .(ge|tr)Matrix */
	int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
	if (ivalid >= 6 && pdim[1] != n)
	    return ScalarLogical(0);
	if (asLogical(checkDN) != 0 &&
	    !DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym)))
	    return ScalarLogical(0);
	Rboolean res = FALSE;
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (ivalid < 6) {
	    /* .trMatrix: symmetric iff diagonal (upper _and_ lower tri.) */
	    char ul = (*uplo_P(obj) == 'U') ? 'L' : 'U';
	    UPM_IS_TR(res, x, n, ul,
		      "'x' slot", "unpackedMatrix_is_symmetric");
	} else {
	    /* .geMatrix: need to do a complete symmetry check */
	    UPM_IS_SY(res, x, n, ivalid == 7,
		      "'x' slot", "unpackedMatrix_is_symmetric");
	}
	return ScalarLogical(res);
    }
}

/* isSymmetric(x, tol = 0) */
SEXP matrix_is_symmetric(SEXP obj, SEXP checkDN)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    if (asLogical(checkDN) != 0) {
	SEXP dn = getAttrib(obj, R_DimNamesSymbol);
	if (!isNull(dn) && !DimNames_is_symmetric(dn))
	    return ScalarLogical(0);
    }
    Rboolean res = FALSE;
    UPM_IS_SY(res, obj, n, 1, "matrix", "matrix_is_symmetric");
    return ScalarLogical(res);
}

#define RETURN_TRUE_OF_KIND(_KIND_)					\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

#define RETURN_GE_IS_TR(_X_, _N_, _UPPER_, _WHAT_, _METHOD_)		\
    do {								\
	Rboolean res = FALSE;						\
	if (_UPPER_ == NA_LOGICAL) {					\
	    UPM_IS_TR(res, _X_, _N_, 'U', _WHAT_, _METHOD_);		\
	    if (res)							\
	        RETURN_TRUE_OF_KIND(mkString("U"));			\
	    UPM_IS_TR(res, _X_, _N_, 'L', _WHAT_, _METHOD_);		\
	    if (res)							\
		RETURN_TRUE_OF_KIND(mkString("L"));			\
	} else {							\
	    UPM_IS_TR(res, _X_, _N_, (_UPPER_ != 0) ? 'U' : 'L',	\
		      _WHAT_, _METHOD_);				\
	    if (res)							\
		return ScalarLogical(1);				\
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
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_is_triangular");
    
    SEXP uplo;
    char ul;
    int need_upper = asLogical(upper);
    if (ivalid < 6) {
	uplo = GET_SLOT(obj, Matrix_uploSym);
	ul = *CHAR(STRING_ELT(uplo, 0));
    }

#define IF_DIAGONAL							\
    Rboolean res = FALSE;						\
    SEXP x = GET_SLOT(obj, Matrix_xSym);				\
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];			\
    UPM_IS_TR(res, x, n, (ul == 'U') ? 'L' : 'U',			\
	      "'x' slot", "unpackedMatrix_is_triangular");		\
    if (res)
    
    if (ivalid < 3) {
	/* .trMatrix: be fast if 'upper', 'uplo' agree; else need diagonal */
	if (need_upper == NA_LOGICAL) {
	    RETURN_TRUE_OF_KIND(uplo);
	} else if ((need_upper != 0 && ul == 'U') ||
		   (need_upper == 0 && ul != 'U')) {
	    return ScalarLogical(1);
	} else {
	    IF_DIAGONAL {
		return ScalarLogical(1);
	    }
	}
	return ScalarLogical(0);
    } else if (ivalid < 6) {
	/* .syMatrix: triangular iff diagonal (upper _and_ lower tri.) */
	IF_DIAGONAL {
	    if (need_upper == NA_LOGICAL) {
		RETURN_TRUE_OF_KIND(mkString("U"));
	    } else {
		return ScalarLogical(1);
	    }
	}
	return ScalarLogical(0);

#undef IF_DIAGONAL
	
    } else {
	/* .geMatrix: need to do a complete triangularity check */
	int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
	if (pdim[1] != n)
	    return ScalarLogical(0);
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	RETURN_GE_IS_TR(x, n, need_upper,
			"'x' slot", "unpackedMatrix_is_triangular");
    }
}

/* isTriangular(x, upper) */
SEXP matrix_is_triangular(SEXP obj, SEXP upper)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    int need_upper = asLogical(upper);
    RETURN_GE_IS_TR(obj, n, need_upper,
		    "matrix", "matrix_is_triangular");
}

#undef RETURN_GE_IS_TR
#undef RETURN_TRUE

/* isDiagonal(x) */
SEXP unpackedMatrix_is_diagonal(SEXP obj)
{
    static const char *valid[] = {
	"dtrMatrix", "ltrMatrix", "ntrMatrix",
	"dsyMatrix", "lsyMatrix", "nsyMatrix",
	"dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_is_diagonal");
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    Rboolean res = FALSE;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (ivalid < 6) {
	/* .(sy|tr)Matrix: diagonal iff stored triangle is zero off diagonal */
	char ul = (*uplo_P(obj) == 'U') ? 'L' : 'U';
	UPM_IS_TR(res, x, n, ul, "'x' slot", "unpackedMatrix_is_diagonal");
    } else {
	/* .geMatrix: need to do a complete diagonality check */
	UPM_IS_DI(res, x, n, "'x' slot", "unpackedMatrix_is_diagonal");
    }
    return ScalarLogical(res);
}

/* isDiagonal(x) */
SEXP matrix_is_diagonal(SEXP obj)
{
    int *pdim = INTEGER(getAttrib(obj, R_DimSymbol)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    Rboolean res = FALSE;
    UPM_IS_DI(res, obj, n, "matrix", "matrix_is_diagonal");
    return ScalarLogical(res);
}

#undef UPM_IS_DI
#undef UPM_IS_TR
#undef UPM_IS_SY

/* t(x), typically preserving class */
/* MJ: Technically no need to do full transpose of .(sy|tr)Matrix ...  */
/*     but then identical(.@x, t(t(.))@x) can be FALSE ...             */
SEXP unpackedMatrix_transpose(SEXP from)
{
    static const char *valid[] = {
	/*  0 */ "Cholesky", "BunchKaufman",
	/*  2 */ "dtrMatrix", "ltrMatrix", "ntrMatrix",
	/*  5 */ "corMatrix", "dpoMatrix",
	/*  7 */ "dsyMatrix", "lsyMatrix", "nsyMatrix",
	/* 10 */ "dgeMatrix", "lgeMatrix", "ngeMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "unpackedMatrix_transpose");
    if (ivalid == 1)
	ivalid = 2; /* BunchKaufman->dtrMatrix */
    const char *cl = valid[ivalid];
    
    SEXPTYPE tx;
    R_xlen_t nx;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x0 = GET_SLOT(from, Matrix_xSym),
	x1 = PROTECT(allocVector(tx = TYPEOF(x0), nx = XLENGTH(x0)));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    char ul = (ivalid < 10) ? *uplo_P(from) : '\0';
    
    /* Preserve or reverse 'Dim' slot (preserving if square) */
    if (m == n) {
	SET_SLOT(to, Matrix_DimSym, dim);
    } else {
	pdim = INTEGER(GET_SLOT(to, Matrix_DimSym));
	pdim[0] = n;
	pdim[1] = m;
    }
    /* Reverse or preserve 'Dimnames' slot (preserving if symmetric) */
    if (ivalid < 5 || ivalid >= 10)
	set_reversed_DimNames(to, dimnames);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (ivalid < 10) {
	/* Toggle 'uplo' slot */
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "L" : "U"));
	if (ivalid < 5)
	    /* Preserve 'diag' slot */
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
	else if (ivalid == 5)
	    /* Preserve 'sd' slot */
	    SET_SLOT(to, install("sd"), GET_SLOT(from, install("sd")));
    }
    /* NB: Nothing to do for 'factors' slot: prototype is already list() ...
       FIXME: However, it would be much better to also "transpose" each 
       factorization ... */

#define UPM_T(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	int i, j;							\
	/* if (ivalid < 10) { */					\
	/*     Memzero(px1, (size_t) n * n); */				\
	/*     R_xlen_t upos = 0, lpos = 0; */				\
	/*     if (ul == 'U') */					\
	/* 	for (j = 0; j < n; upos = (lpos += (++j))) */		\
	/* 	    for (i = j; i < n; ++i, ++lpos, upos += n) */	\
	/* 		px1[lpos] = px0[upos]; */			\
	/*     else */							\
	/* 	for (j = 0; j < n; upos = (lpos += (++j))) */		\
	/* 	    for (i = j; i < n; ++i, ++lpos, upos += n) */	\
	/* 		px1[upos] = px0[lpos]; */			\
	/* } else { */							\
	R_xlen_t nx1s = nx - 1;						\
	for (j = 0; j < m; ++j, px0 -= nx1s)				\
	    for (i = 0; i < n; ++i, px0 += m)				\
		*(px1++) = *px0;					\
	/* } */								\
    } while (0)

    /* Permute 'x' slot */
    switch (tx) {
    case REALSXP: /* d..Matrix */
	UPM_T(double, REAL);
	break;
    case LGLSXP: /* [ln]..Matrix */
	UPM_T(int, LOGICAL);
	break;
    case INTSXP: /* i..Matrix */
	UPM_T(int, INTEGER);
	break;
    case CPLXSXP: /* z..Matrix */
	UPM_T(Rcomplex, COMPLEX);
	break;
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_transpose");
	break;
    }
    SET_SLOT(to, Matrix_xSym, x1);

#undef UPM_T

    UNPROTECT(2);
    return to;
}

/* diag(x, names) */
SEXP unpackedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
	error(_("'names' must be TRUE or FALSE"));

    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    SEXPTYPE tx;
    SEXP x = GET_SLOT(obj, Matrix_xSym),
	res = PROTECT(allocVector(tx = TYPEOF(x), r));
    char ul = *Uplo_P(obj), di = *Diag_P(obj);
    
#define UPM_D_G(_CTYPE_, _PTR_, _ONE_)					\
    do {								\
	_CTYPE_ *pres = _PTR_(res);					\
	int j;								\
	if (di == 'U') {						\
	    for (j = 0; j < r; ++j)					\
		*(pres++) = _ONE_;					\
	} else {							\
	    _CTYPE_ *px = _PTR_(x);					\
	    R_xlen_t m1a = ((R_xlen_t) m) + 1;				\
	    for (j = 0; j < r; ++j, px += m1a)				\
		*(pres++) = *px;					\
	}								\
    } while (0)
    
    switch (tx) {
    case REALSXP: /* d..Matrix */
	UPM_D_G(double, REAL, 1.0);
	break;
    case LGLSXP: /* [ln]..Matrix */
	UPM_D_G(int, LOGICAL, 1);
	break;
    case INTSXP: /* i..Matrix */
	UPM_D_G(int, INTEGER, 1);
	break;
    case CPLXSXP: /* z..Matrix */
    {
	Rcomplex one;
	one.r = 1.0;
	one.i = 0.0;
	UPM_D_G(Rcomplex, COMPLEX, one);
	break;
    }
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_diag_get");
	break;
    }

#undef UPM_D_G

    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	if (isNull(cn)) {
	    if (ul != ' ' && di == ' ' && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (ul != ' ' && di == ' ')
		setAttrib(res, R_NamesSymbol, cn);
	    else if (!isNull(rn) &&
		     (rn == cn || equal_string_vectors(rn, cn, r)))
		setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
	}
    }
    UNPROTECT(1);
    return res;
}

/* diag(x) <- value */
SEXP unpackedMatrix_diag_set(SEXP obj, SEXP val)
{
    SEXP dim = GET_SLOT(obj, Matrix_DimSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    R_xlen_t nv = XLENGTH(val);
    if (nv != 1 && nv != r)
	error(_("replacement diagonal has wrong length"));

    SEXP x = GET_SLOT(obj, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x), tv = TYPEOF(val);
    if (tx < LGLSXP || tx > CPLXSXP)
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_diag_set");
    if (tv < LGLSXP || tv > REALSXP)
	/* Upper bound can become CPLXSXP once we have proper zMatrix */
	error(_("replacement diagonal has incompatible type \"%s\""),
	      type2char(tv));
    
    static const char *valid[] = {
	"dgeMatrix", "dsyMatrix", "dtrMatrix",
	"lgeMatrix", "lsyMatrix", "ltrMatrix",
	"ngeMatrix", "nsyMatrix", "ntrMatrix", ""};
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "unpackedMatrix_diag_set");

    SEXP res;
    int nprotect = 0;
    
    /* Allocate and coerce as necessary */
    if (tv <= tx) {
	if (tv < tx) {
	    PROTECT(val = coerceVector(val, tv = tx));
	    ++nprotect;
	}
	PROTECT(res = NEW_OBJECT_OF_CLASS(valid[ivalid]));
	PROTECT(x = duplicate(x));
	nprotect += 2;
    } else { /* tv > tx */
	/* dMatrix result is only possibility until we have proper [iz]Matrix */
	if (tv < REALSXP) {
	    PROTECT(val = coerceVector(val, tv = REALSXP));
	    ++nprotect;
	}
	char *cl = strdup(valid[ivalid]);
	cl[0] = 'd';
	PROTECT(res = NEW_OBJECT_OF_CLASS(cl));
	PROTECT(x = coerceVector(x, tx = tv));
	++nprotect;
    	free(cl);
    }
    SET_SLOT(res, Matrix_xSym, x);
    
    /* Transfer slots other than 'x', 'diag', and 'factors';
       latter two should keep their prototypes "N" and list() */
    SET_SLOT(res, Matrix_DimSym, dim);
    SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
    if (R_has_slot(res, Matrix_uploSym))
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));

#define UPM_D_S(_CTYPE_, _PTR_)						\
    do {								\
	_CTYPE_ *px = _PTR_(x), *pval = _PTR_(val);			\
	R_xlen_t m1a = (R_xlen_t) m + 1;				\
	int j;								\
	if (nv == 1)							\
	    for (j = 0; j < r; ++j, px += m1a)				\
		*px = *pval;						\
	else								\
	    for (j = 0; j < r; ++j, px += m1a)				\
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
    default:
	ERROR_INVALID_TYPE("'x' slot", tx, "unpackedMatrix_diag_set");
	break;
    }

#undef UPM_D_S
    
    UNPROTECT(nprotect);
    return res;
}