#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "sparse.h"

#define SPARSE_ERROR_INVALID_CLASS(_CLASS_, _METHOD_)			\
    error(_("invalid class \"%s\" to 'R_sparse_%s()'"),			\
	  _CLASS_, _METHOD_)

#define DIAGONAL_ERROR_INVALID_CLASS(_CLASS_, _METHOD_)			\
    error(_("invalid class \"%s\" to 'R_diagonal_%s()'"),		\
	  _CLASS_, _METHOD_)

#define NZREAL(_X_) ((_X_) != 0.0)

#define NZINTEGER(_X_) ((_X_) != 0)

#define NZCOMPLEX(_X_) ((_X_).r != 0.0 || (_X_).i != 0.0)

/* as(<sparseMatrix>, "[dlniz](sparse)?Matrix") */
SEXP R_sparse_as_kind(SEXP from, SEXP kind)
{
    char k;
    if ((kind = asChar(kind)) == NA_STRING || (k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_sparse_as_kind()'"));
    static const char *valid[] = {
	"dgCMatrix", "dsCMatrix", "dtCMatrix",
	"dgRMatrix", "dsRMatrix", "dtRMatrix",
	"dgTMatrix", "dsTMatrix", "dtTMatrix",
	"lgCMatrix", "lsCMatrix", "ltCMatrix",
	"lgRMatrix", "lsRMatrix", "ltRMatrix",
	"lgTMatrix", "lsTMatrix", "ltTMatrix",
	"ngCMatrix", "nsCMatrix", "ntCMatrix",
	"ngRMatrix", "nsRMatrix", "ntRMatrix",
	"ngTMatrix", "nsTMatrix", "ntTMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	SPARSE_ERROR_INVALID_CLASS(class_P(from), "as_kind");
    const char *clf = valid[ivalid];
    if (k == '.')
	k = clf[0];
    if (k == clf[0])
	return from;
    SEXPTYPE tx = kind2type(k); /* validating 'k' before doing more */

    if (k == 'n')
	PROTECT(from = (clf[2] == 'C'
			? Csparse_drop0(from)
			: (clf[2] == 'R'
			   ? Rsparse_drop0(from)
			   : Tsparse_drop0(from))));
    
    char *clt = strdup(clf);
    clt[0] = k;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt));
    free(clt);
    
    SEXP i, j;
    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    if (clf[2] != 'T')
	SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
    if (clf[2] != 'R')
	SET_SLOT(to, Matrix_iSym, i = GET_SLOT(from, Matrix_iSym));
    if (clf[2] != 'C')
	SET_SLOT(to, Matrix_jSym, j = GET_SLOT(from, Matrix_jSym));
    
    if (clf[0] == 'n') {
	R_xlen_t ix, nx = (clf[2] == 'C') ? XLENGTH(i) : XLENGTH(j);
	SEXP x = PROTECT(allocVector(tx, nx));

#define SET_ONES(_X_, _I_, _N_, _CTYPE_, _PTR_, _ONE_)	\
	do {						\
	    _CTYPE_ *px = _PTR_(_X_);			\
	    for (_I_ = 0; _I_ < _N_; ++_I_)		\
		*(px++) = _ONE_;			\
	} while (0)

#define SET_ONES_CASES(_SEXPTYPE_, _X_, _I_, _N_)			\
	do {								\
	    switch (_SEXPTYPE_) {					\
	    case REALSXP:						\
		SET_ONES(_X_, _I_, _N_, double, REAL, 1.0);		\
		break;							\
	    case LGLSXP:						\
		SET_ONES(_X_, _I_, _N_, int, LOGICAL, 1);		\
		break;							\
	    case INTSXP:						\
		SET_ONES(_X_, _I_, _N_, int, INTEGER, 1);		\
		break;							\
	    case CPLXSXP:						\
	    {								\
		Rcomplex one; one.r = 1.0; one.i = 0.0;			\
		SET_ONES(_X_, _I_, _N_, Rcomplex, COMPLEX, one);	\
		break;							\
	    }								\
	    default:							\
		break;							\
	    }								\
	} while (0)
	
	SET_ONES_CASES(tx, x, ix, nx);
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
    } else if (k != 'n') {
	SET_SLOT(to, Matrix_xSym,
		 coerceVector(GET_SLOT(from, Matrix_xSym), tx));
    }

    UNPROTECT((k == 'n') ? 2 : 1);
    return to;
}

/* as(<diagonalMatrix>, "(.di|.[gst][CRT]|.sparse|.)Matrix") */
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0)
{
    const char *c;
    char c0, c1, c2, ul;
    if ((code = asChar(code)) == NA_STRING ||
	(c0 = (c = CHAR(code))[0]) == '\0' ||
	(c1 = c[1]) == '\0' ||
	(c2 = c[2]) == '\0')
	error(_("invalid 'code' to 'R_diagonal_as_sparse()'"));
    if ((c1 == 't' || c1 == 's') &&
	(((uplo = asChar(uplo)) == NA_STRING || (ul = *CHAR(uplo)) == '\0')))
	error(_("invalid 'uplo' to 'R_diagonal_as_sparse()'"));
    
    static const char *valid[] = { "ddiMatrix", "ldiMatrix", "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	DIAGONAL_ERROR_INVALID_CLASS(class_P(from), "as_sparse");

    const char *clf = valid[ivalid];
    if (c0 == '.')
	c0 = clf[0];
    if (c0 == clf[0] && c1 == clf[1] && c2 == clf[2])
	return from;
    
    SEXPTYPE tx = kind2type(c0); /* validating 'c0' before doing more */
    char clt[] = "...Matrix";
    clt[0] = c0;
    clt[1] = c1;
    clt[2] = c2;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	diag = GET_SLOT(from, Matrix_diagSym);
    char di = *CHAR(STRING_ELT(diag, 0));

    SET_SLOT(to, Matrix_DimSym, dim);
    if (c1 == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (c1 == 'd' || c1 == 't')
	SET_SLOT(to, Matrix_diagSym, diag);

    /* Be fast when coercing between subclasses of diagonalMatrix ... */

    if (c1 == 'd') {
	SET_SLOT(to, Matrix_xSym,
		 (di != 'N'
		  ? allocVector(tx, 0)
		  : coerceVector(GET_SLOT(from, Matrix_xSym), tx)));
	UNPROTECT(1);
	return to;
    }

    /* Now coercing from diagonalMatrix to .sparseMatrix ... */
    
    if (c1 == 't' || c1 == 's')
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "U" : "L"));

    SEXP p, i, x;
    int k, n = INTEGER(dim)[0], nprotect = 1;
    R_xlen_t n1a = (R_xlen_t) n + 1;
    
    if (di != 'N') { /* unit diagonal */
	
	if (c2 != 'T') {
	    PROTECT(p = allocVector(INTSXP, n1a)); ++nprotect;
	    int *pp = INTEGER(p);
	    if (c1 == 't')
		Memzero(pp, n1a);
	    else
		for (k = 0; k <= n; ++k)
		    *(pp++) = k;
	}
	if (c1 == 't') {
	    PROTECT(i = allocVector(INTSXP, 0)); ++nprotect;
	    if (c0 != 'n') {
		PROTECT(x = allocVector(tx, 0)); ++nprotect;
	    }
	} else {
	    PROTECT(i = allocVector(INTSXP, n)); ++nprotect;
	    int *pi = INTEGER(i);
	    if (c0 == 'n') {
		for (k = 0; k < n; ++k)
		    *(pi++) = k;
	    } else {
		PROTECT(x = allocVector(tx, n)); ++nprotect;

#undef SET_ONES
#define SET_ONES(_X_, _I_, _N_, _CTYPE_, _PTR_, _ONE_)		\
		do {						\
		    _CTYPE_ *px = _PTR_(_X_);			\
		    for (_I_ = 0; _I_ < _N_; ++_I_) {		\
			*(pi++) = _I_;				\
			*(px++) = _ONE_;			\
		    }						\
		} while (0)
		
		SET_ONES_CASES(tx, x, k, n);

#undef SET_ONES_CASES
#undef SET_ONES
		
	    }
	}
	
    } else { /* non-unit diagonal */
	
	if (TYPEOF(x = GET_SLOT(from, Matrix_xSym)) != tx) {
	    PROTECT(x = coerceVector(x, tx)); ++nprotect;
	}
	    
	if (c0 == 'n' || asLogical(drop0) != 0) { /* _do_ drop zeros (if any) */
	    
	    int nnz = 0;
	    
#define DROP0_DIAGONAL(_CTYPE_, _PTR_, _NZ_)				\
	    do {							\
		_CTYPE_ *px = _PTR_(x);					\
		if (c2 == 'T') {					\
		    for (k = 0; k < n; ++k)				\
			if (_NZ_(px[k]))				\
			    ++nnz;					\
		} else {						\
		    PROTECT(p = allocVector(INTSXP, n1a)); ++nprotect;	\
		    int *pp = INTEGER(p);				\
		    *(pp++) = 0;					\
		    for (k = 0; k < n; ++k)				\
			*(pp++) = (_NZ_(px[k])) ? ++nnz : nnz;		\
		}							\
		PROTECT(i = allocVector(INTSXP, nnz)); ++nprotect;	\
		if (nnz != n && c0 != 'n') {				\
		    PROTECT(x = allocVector(tx, nnz)); ++nprotect;	\
		}							\
		if (nnz == 0)						\
		    continue;						\
		int *pi = INTEGER(i);					\
		if (nnz == n) {						\
		    for (k = 0; k < n; ++k)				\
			*(pi++) = k;					\
		} else if (c0 == 'n') {					\
		    for (k = 0; k < n; ++k)				\
			if (_NZ_(px[k]))				\
			    *(pi++) = k;				\
		} else {						\
		    _CTYPE_ *newpx = _PTR_(x);				\
		    for (k = 0; k < n; ++k) {				\
			if (_NZ_(px[k])) {				\
			    *(pi++) = k;				\
			    *(newpx++) = px[k];				\
			}						\
		    }							\
		}							\
	    } while (0)

#define DROP0_CASES(_SEXPTYPE_, _MACRO_)			\
	    do {						\
		switch (_SEXPTYPE_) {				\
		case REALSXP:					\
		    _MACRO_(double, REAL, NZREAL);		\
		    break;					\
		case LGLSXP:					\
		    _MACRO_(int, LOGICAL, NZINTEGER);		\
		    break;					\
		case INTSXP:					\
		    _MACRO_(int, INTEGER, NZINTEGER);		\
		    break;					\
		case CPLXSXP:					\
		    _MACRO_(Rcomplex, COMPLEX, NZCOMPLEX);	\
		    break;					\
		default:					\
		    break;					\
		}						\
	    } while (0)

	    DROP0_CASES(tx, DROP0_DIAGONAL);
	    
	} else { /* _don't_ drop zeros */
	    
	    PROTECT(i = allocVector(INTSXP, n)); ++nprotect;
	    int *pi = INTEGER(i);
	    if (c2 == 'T') {
		PROTECT(p = allocVector(INTSXP, n1a)); ++nprotect;
		int *pp = INTEGER(p);
		*(pp++) = 0;
		for (k = 0; k < n; ++k)
		    *(pi++) = *(pp++) = k;
	    } else {
		for (k = 0; k < n; ++k)
		    *(pi++) = k;
	    }
	    
	}
    }

    if (c2 != 'T')
	SET_SLOT(to, Matrix_pSym, p);
    if (c2 != 'R')
	SET_SLOT(to, Matrix_iSym, i);
    if (c2 != 'C')
	SET_SLOT(to, Matrix_jSym, i);
    if (c0 != 'n')
	SET_SLOT(to, Matrix_xSym, x);
    UNPROTECT(nprotect);
    return to;
}

SEXP Csparse_drop0(SEXP from)
{
    
#define DROP0_START_INIT			\
    SEXP x0 = GET_SLOT(from, Matrix_xSym);	\
    SEXPTYPE tx = TYPEOF(x0);			\
    int nx = LENGTH(x0), i, nnz = 0, nnz_;

    DROP0_START_INIT;
    
#define DROP0_START(_CTYPE_, _PTR_, _NZ_)			\
    do {							\
	_CTYPE_ *px0 = _PTR_(x0);				\
	while (nnz < nx && _NZ_(*px0)) {			\
	    ++nnz; ++px0;					\
	}							\
	if (nnz == nx)						\
	    return from;					\
	nnz_ = nnz;						\
	for (i = nnz_; i < nx; ++i, ++px0)			\
	    if (_NZ_(*px0))					\
		++nnz;						\
    } while (0)
    
    DROP0_CASES(tx, DROP0_START);

#define DROP0_END_INIT(_ISYM_)						\
    R_xlen_t ncol1a, col1a;						\
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(class_P(from))),		\
	p0 = GET_SLOT(from, Matrix_pSym),				\
	i0 = GET_SLOT(from, _ISYM_),					\
	p1 = PROTECT(allocVector(INTSXP, ncol1a = XLENGTH(p0))),	\
	i1 = PROTECT(allocVector(INTSXP, nnz)),				\
	x1 = PROTECT(allocVector(tx, nnz));				\
    int s, t,								\
	*pp0 = INTEGER(p0),						\
	*pi0 = INTEGER(i0),						\
	*pp1 = INTEGER(p1),						\
	*pi1 = INTEGER(i1);
    
    DROP0_END_INIT(Matrix_iSym);
    
#define DROP0_END(_CTYPE_, _PTR_, _NZ_)			\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	Memcpy(pi1, pi0, nnz_);				\
	Memcpy(px1, px0, nnz_);				\
	col1a = 0;					\
	while ((s = pp0[col1a]) <= nnz_)		\
	    pp1[col1a++] = s;				\
	for (i = nnz_; i < s; ++i) {			\
	    if (_NZ_(px0[i])) {				\
		pi1[nnz_] = pi0[i];			\
		px1[nnz_] = px0[i];			\
		++nnz_;					\
	    }						\
	}						\
	pp1[col1a] = nnz_;				\
	while (++col1a < ncol1a) {			\
	    t = pp0[col1a];				\
	    for (i = s; i < t; ++i) {			\
		if (_NZ_(px0[i])) {			\
		    pi1[nnz_] = pi0[i];			\
		    px1[nnz_] = px0[i];			\
		    ++nnz_;				\
		}					\
	    }						\
	    pp1[col1a] = nnz_;				\
	    s = t;					\
	}						\
    } while (0)
    
    DROP0_CASES(tx, DROP0_END);

#define DROP0_SET_COMMON						\
    do {								\
	SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));	\
	SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym)); \
	SET_SLOT(to, Matrix_xSym, x1);					\
	if (R_has_slot(from, Matrix_uploSym))				\
	    SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym)); \
	if (R_has_slot(from, Matrix_diagSym))				\
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym)); \
	else								\
	    SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym)); \
    } while (0)

    DROP0_SET_COMMON;

#define DROP0_SET_EXTRA(_ISYM_)			\
    do {					\
	SET_SLOT(to, Matrix_pSym, p1);		\
	SET_SLOT(to, _ISYM_, i1);		\
    } while (0)
    
    DROP0_SET_EXTRA(Matrix_iSym);
    
    UNPROTECT(4);
    return to;  
}

SEXP Rsparse_drop0(SEXP from)
{
    DROP0_START_INIT;
    DROP0_CASES(tx, DROP0_START);
    DROP0_END_INIT(Matrix_jSym);
    DROP0_CASES(tx, DROP0_END);
    DROP0_SET_COMMON;
    DROP0_SET_EXTRA(Matrix_jSym);
    UNPROTECT(4);
    return to;
}

SEXP Tsparse_drop0(SEXP from)
{
    DROP0_START_INIT;
    DROP0_CASES(tx, DROP0_START);
    
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(class_P(from))),
	i0 = GET_SLOT(from, Matrix_iSym),
	j0 = GET_SLOT(from, Matrix_jSym),
	i1 = PROTECT(allocVector(INTSXP, nnz)),
	j1 = PROTECT(allocVector(INTSXP, nnz)),
	x1 = PROTECT(allocVector(tx, nnz));
    int
	*pi0 = INTEGER(i0),
	*pj0 = INTEGER(j0),
	*pi1 = INTEGER(i1),
	*pj1 = INTEGER(j1);

#undef DROP0_END
#define DROP0_END(_CTYPE_, _PTR_, _NZ_)				\
    do {							\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
	Memcpy(pi1, pi0, nnz_);					\
	Memcpy(pj1, pj0, nnz_);					\
	Memcpy(px1, px0, nnz_);					\
	for (i = nnz_; i < nx; ++i) {				\
	    if (_NZ_(px0[i])) {					\
		pi1[nnz_] = pi0[i];				\
		pj1[nnz_] = pj0[i];				\
		px1[nnz_] = px0[i];				\
		++nnz_;						\
	    }							\
	}							\
    } while (0)

    DROP0_CASES(tx, DROP0_END);
    DROP0_SET_COMMON;
    SET_SLOT(to, Matrix_iSym, i1);
    SET_SLOT(to, Matrix_jSym, j1);
    UNPROTECT(4);
    return to;
}

#undef DROP0_CASES
#undef DROP0_START
#undef DROP0_START_INIT
#undef DROP0_END
#undef DROP0_END_INIT
#undef DROP0_SET_COMMON
#undef DROP0_SET_EXTRA
