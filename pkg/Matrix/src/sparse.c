#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "sparse.h"

/* drop0(<[CRT]sparseMatrix>) */
SEXP R_sparse_drop0(SEXP from)
{
    static const char *valid[] = {
	"dgCMatrix", "dtCMatrix", "dsCMatrix",
	"dgRMatrix", "dtRMatrix", "dsRMatrix",
	"dgTMatrix", "dtTMatrix", "dsTMatrix",
	"lgCMatrix", "ltCMatrix", "lsCMatrix",
	"lgRMatrix", "ltRMatrix", "lsRMatrix",
	"lgTMatrix", "ltTMatrix", "lsTMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_drop0");
    const char *cl = valid[ivalid];
    
    SEXP p0, x0 = GET_SLOT(from, Matrix_xSym);
    SEXPTYPE tx = TYPEOF(x0);
    int *pp0, k, n, nnz_, nnz0, nnz1 = 0;
    R_xlen_t n1a;

    if (cl[2] == 'T') {
	nnz0 = LENGTH(x0);
    } else {
	p0 = GET_SLOT(from, Matrix_pSym);
	pp0 = INTEGER(p0);
	n1a = XLENGTH(p0);
	n = (int) (n1a - 1);
	nnz0 = pp0[n];
    }

#define DROP0_CASES(_SEXPTYPE_, _MACRO_)			\
    do {							\
	switch (_SEXPTYPE_) {					\
	case REALSXP:						\
	    _MACRO_(double, REAL, NZ_REAL);			\
	    break;						\
	case LGLSXP:						\
	    _MACRO_(int, LOGICAL, NZ_INTEGER);			\
	    break;						\
	case INTSXP:						\
	    _MACRO_(int, INTEGER, NZ_INTEGER);			\
	    break;						\
	case CPLXSXP:						\
	    _MACRO_(Rcomplex, COMPLEX, NZ_COMPLEX);		\
	    break;						\
	default:						\
	    break;						\
	}							\
    } while (0)
    
#define DROP0_START(_CTYPE_, _PTR_, _NZ_)			\
    do {							\
	_CTYPE_ *px0 = _PTR_(x0);				\
	while (nnz1 < nnz0 && _NZ_(*px0)) {			\
	    ++nnz1;						\
	    ++px0;						\
	}							\
	if (nnz1 == nnz0)					\
	    return from;					\
	nnz_ = nnz1;						\
	for (k = nnz_; k < nnz0; ++k, ++px0)			\
	    if (_NZ_(*px0))					\
		++nnz1;						\
    } while (0)
    
    DROP0_CASES(tx, DROP0_START);

#undef DROP0_START    

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (cl[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (cl[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));

    /* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

    SEXP isym = (cl[2] == 'R') ? Matrix_jSym : Matrix_iSym,
	i0 = GET_SLOT(from, isym),
	i1 = PROTECT(allocVector(INTSXP, nnz1)),
	x1 = PROTECT(allocVector(tx, nnz1));
    int *pi0 = INTEGER(i0), *pi1 = INTEGER(i1);
    SET_SLOT(to, isym, i1);
    SET_SLOT(to, Matrix_xSym, x1);
    
    if (cl[2] == 'T') {

	SEXP j0 = GET_SLOT(from, Matrix_jSym),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pj0 = INTEGER(j0), *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_jSym, j1);

#define DROP0_END_T(_CTYPE_, _PTR_, _NZ_)			\
	do {							\
	    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
	    Memcpy(pi1, pi0, nnz_);				\
	    Memcpy(pj1, pj0, nnz_);				\
	    Memcpy(px1, px0, nnz_);				\
	    for (k = nnz_; k < nnz0; ++k) {			\
		if (_NZ_(px0[k])) {				\
		    pi1[nnz_] = pi0[k];				\
		    pj1[nnz_] = pj0[k];				\
		    px1[nnz_] = px0[k];				\
		    ++nnz_;					\
		}						\
	    }							\
	} while (0)

	DROP0_CASES(tx, DROP0_END_T);

#undef DROP0_END_T
	
    } else { /* cl[2] == 'C' || cl[2] == 'R' */

	SEXP p1 = PROTECT(allocVector(INTSXP, n1a));
	int *pp1 = INTEGER(p1), j, kend;
	SET_SLOT(to, Matrix_pSym, p1);

#define DROP0_END_CR(_CTYPE_, _PTR_, _NZ_)		\
	do {						\
	    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	    Memcpy(pi1, pi0, nnz_);			\
	    Memcpy(px1, px0, nnz_);			\
	    j = 0;					\
	    while ((kend = pp0[j]) <= nnz_)		\
		pp1[j++] = kend;			\
	    for (k = nnz_; k < kend; ++k) {		\
		if (_NZ_(px0[k])) {			\
		    pi1[nnz_] = pi0[k];			\
		    px1[nnz_] = px0[k];			\
		    ++nnz_;				\
		}					\
	    }						\
	    pp1[j] = nnz_;				\
	    while (++j <= n) {				\
		kend = pp0[j];				\
		while (k < kend) {			\
		    if (_NZ_(px0[k])) {			\
			pi1[nnz_] = pi0[k];		\
			px1[nnz_] = px0[k];		\
			++nnz_;				\
		    }					\
		    ++k;				\
		}					\
		pp1[j] = nnz_;				\
	    }						\
	} while (0)
    
	DROP0_CASES(tx, DROP0_END_CR);

#undef DROP0_END_CR
	
    }

    UNPROTECT(4);
    return to;
}

/* as(<[CRT]sparseMatrix>, "[dlniz]Matrix") */
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0)
{
    char k;
    if ((kind = asChar(kind)) == NA_STRING || (k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_sparse_as_kind()'"));
    static const char *valid[] = {
	"dgCMatrix", "dtCMatrix", "dsCMatrix",
	"dgRMatrix", "dtRMatrix", "dsRMatrix",
	"dgTMatrix", "dtTMatrix", "dsTMatrix",
	"lgCMatrix", "ltCMatrix", "lsCMatrix",
	"lgRMatrix", "ltRMatrix", "lsRMatrix",
	"lgTMatrix", "ltTMatrix", "lsTMatrix",
	"ngCMatrix", "ntCMatrix", "nsCMatrix",
	"ngRMatrix", "ntRMatrix", "nsRMatrix",
	"ngTMatrix", "ntTMatrix", "nsTMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_as_kind");
    const char *clf = valid[ivalid];
    if (k == '.')
	k = clf[0];
    if (k == clf[0])
	return from;
    SEXPTYPE tx = kind2type(k); /* validating 'k' before doing more */

#if 0 /* MJ: what I think makes sense */
    int do_drop0 = (k == 'n' ||
		    asLogical(drop0) != 0);
#else /* MJ: behaviour in Matrix 1.4-1 ?? */
    int do_drop0 = ((clf[0] != 'n' && clf[2] != 'T' && k == 'l') ||
		    asLogical(drop0) != 0);
#endif
    if (do_drop0)
	PROTECT(from = R_sparse_drop0(from));
    
    char *clt = strdup(clf);
    clt[0] = k;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)), p, i;

    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    if (clf[2] != 'T')
	SET_SLOT(to, Matrix_pSym, p = GET_SLOT(from, Matrix_pSym));
    if (clf[2] != 'R')
	SET_SLOT(to, Matrix_iSym, i = GET_SLOT(from, Matrix_iSym));
    if (clf[2] != 'C')
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_jSym));
    if (clf[0] != 'n' && clt[0] != 'n') {
	SET_SLOT(to, Matrix_xSym,
		 coerceVector(GET_SLOT(from, Matrix_xSym), tx));
    } else if (clf[0] == 'n') {
	int ix, nx = (clf[2] == 'T') ? LENGTH(i) : INTEGER(p)[XLENGTH(p) - 1];
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
    }

    free(clt);
    UNPROTECT(do_drop0 ? 2 : 1);
    return to;
}

/* as(<diagonalMatrix>, ".(di|[gst][CRT])Matrix") */
SEXP R_diagonal_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP drop0)
{
    const char *zzz;
    char z0, z1, z2, ul;
    if ((code = asChar(code)) == NA_STRING ||
	(z0 = (zzz = CHAR(code))[0]) == '\0' ||
	(z1 = zzz[1]) == '\0' ||
	(z2 = zzz[2]) == '\0')
	error(_("invalid 'code' to 'R_diagonal_as_sparse()'"));
    if ((z1 == 't' || z1 == 's') &&
	(((uplo = asChar(uplo)) == NA_STRING || (ul = *CHAR(uplo)) == '\0')))
	error(_("invalid 'uplo' to 'R_diagonal_as_sparse()'"));
    
    static const char *valid[] = { "ddiMatrix", "ldiMatrix", "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_diagonal_as_sparse");

    const char *clf = valid[ivalid];
    if (z0 == '.')
	z0 = clf[0];
    if (z0 == clf[0] && z1 == clf[1] && z2 == clf[2])
	return from;
    
    SEXPTYPE tx = kind2type(z0); /* validating 'z0' before doing more */
    char clt[] = "...Matrix";
    clt[0] = z0;
    clt[1] = z1;
    clt[2] = z2;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	diag = GET_SLOT(from, Matrix_diagSym);
    char di = *CHAR(STRING_ELT(diag, 0));

    SET_SLOT(to, Matrix_DimSym, dim);
    if (z1 == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (z1 == 'd' || z1 == 't')
	SET_SLOT(to, Matrix_diagSym, diag);

    /* Be fast when coercing between subclasses of diagonalMatrix ... */

    if (z1 == 'd') {
	SET_SLOT(to, Matrix_xSym,
		 (di != 'N'
		  ? allocVector(tx, 0)
		  : coerceVector(GET_SLOT(from, Matrix_xSym), tx)));
	UNPROTECT(1);
	return to;
    }

    /* Now coercing from diagonalMatrix to [CRT]sparseMatrix ... */
    
    if (z1 == 't' || z1 == 's')
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "U" : "L"));

    SEXP p, i, x;
    int k, n = INTEGER(dim)[0], nprotect = 1;
    R_xlen_t n1a = (R_xlen_t) n + 1;
    
    if (di != 'N') { /* unit diagonal */
	
	if (z2 != 'T') {
	    PROTECT(p = allocVector(INTSXP, n1a));
	    ++nprotect;
	    int *pp = INTEGER(p);
	    if (z1 == 't')
		Memzero(pp, n1a);
	    else
		for (k = 0; k <= n; ++k)
		    *(pp++) = k;
	}
	if (z1 == 't') {
	    PROTECT(i = allocVector(INTSXP, 0));
	    ++nprotect;
	    if (z0 != 'n') {
		PROTECT(x = allocVector(tx, 0));
		++nprotect;
	    }
	} else {
	    PROTECT(i = allocVector(INTSXP, n));
	    ++nprotect;
	    int *pi = INTEGER(i);
	    if (z0 == 'n') {
		for (k = 0; k < n; ++k)
		    *(pi++) = k;
	    } else {
		PROTECT(x = allocVector(tx, n));
		++nprotect;

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
	    PROTECT(x = coerceVector(x, tx));
	    ++nprotect;
	}
	    
	if (z0 == 'n' || asLogical(drop0) != 0) { /* _do_ drop zeros (if any) */
	    
	    int nnz = 0;
	    
#define DROP0_DIAGONAL(_CTYPE_, _PTR_, _NZ_)				\
	    do {							\
		_CTYPE_ *px = _PTR_(x);					\
		if (z2 == 'T') {					\
		    for (k = 0; k < n; ++k)				\
			if (_NZ_(px[k]))				\
			    ++nnz;					\
		} else {						\
		    PROTECT(p = allocVector(INTSXP, n1a));		\
		    ++nprotect;						\
		    int *pp = INTEGER(p);				\
		    *(pp++) = 0;					\
		    for (k = 0; k < n; ++k)				\
			*(pp++) = (_NZ_(px[k])) ? ++nnz : nnz;		\
		}							\
		PROTECT(i = allocVector(INTSXP, nnz));			\
		++nprotect;						\
		if (nnz != n && z0 != 'n') {				\
		    PROTECT(x = allocVector(tx, nnz));			\
		    ++nprotect;						\
		}							\
		if (nnz == 0)						\
		    continue;						\
		int *pi = INTEGER(i);					\
		if (nnz == n) {						\
		    for (k = 0; k < n; ++k)				\
			*(pi++) = k;					\
		} else if (z0 == 'n') {					\
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

	    DROP0_CASES(tx, DROP0_DIAGONAL);
	    
	} else { /* _don't_ drop zeros */
	    
	    PROTECT(i = allocVector(INTSXP, n));
	    ++nprotect;
	    int *pi = INTEGER(i);
	    if (z2 == 'T') {
		PROTECT(p = allocVector(INTSXP, n1a));
		++nprotect;
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

    if (z2 != 'T')
	SET_SLOT(to, Matrix_pSym, p);
    if (z2 != 'R')
	SET_SLOT(to, Matrix_iSym, i);
    if (z2 != 'C')
	SET_SLOT(to, Matrix_jSym, i);
    if (z0 != 'n')
	SET_SLOT(to, Matrix_xSym, x);
    UNPROTECT(nprotect);
    return to;
}

/* band(<sparseMatrix>, k1, k2), tri[ul](<sparseMatrix>, k) */
/* NB: argument validation more or less copied from R_dense_band() */
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2)
{
    static const char *valid[] = {
	"dgCMatrix", "dtCMatrix", "dsCMatrix",
	"dgRMatrix", "dtRMatrix", "dsRMatrix",
	"dgTMatrix", "dtTMatrix", "dsTMatrix",
	"lgCMatrix", "ltCMatrix", "lsCMatrix",
	"lgRMatrix", "ltRMatrix", "lsRMatrix",
	"lgTMatrix", "ltTMatrix", "lsTMatrix",
	"ngCMatrix", "ntCMatrix", "nsCMatrix",
	"ngRMatrix", "ntRMatrix", "nsRMatrix",
	"ngTMatrix", "ntTMatrix", "nsTMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_band");
    const char *clf = valid[ivalid];
    
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], a, b;
    if (isNull(k1))
	a = (m > 0) ? 1 - m : 0;
    else if ((a = asInteger(k1)) == NA_INTEGER || a < -m || a > n)
	error(_("'k1' must be an integer from -Dim[1] to Dim[2]"));
    if (isNull(k2))
	b = (n > 0) ? n - 1 : 0;
    else if ((b = asInteger(k2)) == NA_INTEGER || b < -m || b > n)
	error(_("'k2' must be an integer from -Dim[1] to Dim[2]"));
    else if (b < a)
	error(_("'k1' must be less than or equal to 'k2'"));
    /* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
    if (a <= 1 - m && b >= n - 1 && (clf[1] == 't' || m != n || m > 1 || n > 1))
	return from;

    char ulf = 'U', ult = 'U', di = 'N';
    if (clf[1] != 'g') {
	ulf = *uplo_P(from);
	if (clf[1] == 't') {
	    /* Be fast if band contains entire triangle */
	    if ((ulf == 'U') ? (a <= 0 && b >= n - 1) : (b >= 0 && a <= 1 - m))
		return from;
	    if (a <= 0 && b >= 0)
		di = *diag_P(from);
	}
    }

    char *clt = strdup(clf);
    int ge = 0, tr = 0, sy = 0, nprotect = 0;
    ge = m != n || (!(tr = a >= 0 || b <= 0 || clf[1] == 't') &&
		    !(sy = a == -b && clf[1] == 's'));   
    clt[1] = (ge) ? 'g' : ((tr) ? 't' : 's');

    /* band(<R>, a, b) is equivalent to t(band(t(<R>), -b, -a)) ! */
    if (clf[2] == 'R') {
	int tmp;
	tmp = m; m =  n; n =  tmp;
	tmp = a; a = -b; b = -tmp;
	ulf = (ulf == 'U') ? 'L' : 'U';
	clt[2] = 'C';
	PROTECT(from = tCRsparse_as_RCsparse(from));
	++nprotect;
	if (m != n)
	    dim = GET_SLOT(from, Matrix_DimSym);
    }

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    ++nprotect;
    free(clt);

    SET_SLOT(to, Matrix_DimSym, dim);
    if (!sy && clf[1] == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);

    if (!ge) {
	if (sy || clf[1] == 't') {
	    /* We can tighten the band ... but does it help?? FIXME */
	    if (ulf == 'U') {
		if (b >= 0 && a < 0)
		    a = 0;
	    } else {
		if (a <= 0 && b > 0)
		    b = 0;
	    }
	}
	
	ult = (tr && clf[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ulf;
	if (ult != 'U')
	    SET_SLOT(to, Matrix_uploSym, mkString("L"));
	if (di != 'N')
	    SET_SLOT(to, Matrix_diagSym, mkString("U"));
    }

    /* It remains to set some subset of 'p', 'i', 'j', 'x' ... */
    
    SEXP p0, p1, i0, j0, x0;
    int *pp0, *pp1, *pi0, *pj0, nnz0, nnz1, d, j, k, kend;

    i0 = GET_SLOT(from, Matrix_iSym);
    pi0 = INTEGER(i0);
    if (clf[2] == 'T') {
	j0 = GET_SLOT(from, Matrix_jSym);
	pj0 = INTEGER(j0);
	nnz0 = LENGTH(j0);
    } else {
	p0 = GET_SLOT(from, Matrix_pSym);
	pp0 = INTEGER(p0);
	nnz0 = pp0[n];
	++pp0;
    }
    if (clf[0] != 'n')
	x0 = GET_SLOT(from, Matrix_xSym);    
    
    /* Count number of nonzero elements in band */
    nnz1 = 0;
    if (clf[2] == 'T') {

	if (!sy && clf[1] == 's') {
	    for (k = 0; k < nnz0; ++k) {
		if ((d = pj0[k] - pi0[k]) >= a && d <= b)
		    ++nnz1;
		if (d != 0 && -d >= a && -d <= b)
		    ++nnz1;
	    }
	} else {
	    for (k = 0; k < nnz0; ++k) {
		if ((d = pj0[k] - pi0[k]) >= a && d <= b)
		    ++nnz1;
	    }
	}

    } else {

	PROTECT(p1 = allocVector(INTSXP, (R_xlen_t) n + 1));
	SET_SLOT(to, Matrix_pSym, p1);
	pp1 = INTEGER(p1);
	++nprotect;
	*(pp1++) = 0;
	k = 0;
	
	if (!sy && clf[1] == 's') {
	    Memzero(pp1, n);
	    for (j = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if ((d = j - pi0[k]) >= a && d <= b)
			++pp1[j];
		    if (d != 0 && -d >= a && -d <= b)
			++pp1[pi0[k]];
		    ++k;
		}
	    }
	    for (j = 0; j < n; ++j) {
		nnz1 += pp1[j];
		pp1[j] = nnz1;
	    }
	} else {
	    for (j = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if ((d = j - pi0[k]) >= a && d <= b)
			++nnz1;
		    ++k;
		}
		pp1[j] = nnz1;
	    }
	}
	
    }

    /* Be fast if band contains all nonzero elements */
    if (nnz1 == nnz0 && (sy || clf[1] != 's')) {
	SET_SLOT(to, Matrix_iSym, i0);
	if (clf[2] == 'T')
	    SET_SLOT(to, Matrix_jSym, j0);
	if (clf[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, x0);	
	if (clf[2] == 'R')
	    to = tCRsparse_as_RCsparse(to);
	UNPROTECT(nprotect);
	return to;
    }

    SEXP i1, j1, x1;
    int *pi1, *pj1;
    
    PROTECT(i1 = allocVector(INTSXP, nnz1));
    SET_SLOT(to, Matrix_iSym, i1);
    pi1 = INTEGER(i1);
    ++nprotect;
    if (clf[2] == 'T') {
	PROTECT(j1 = allocVector(INTSXP, nnz1));
	SET_SLOT(to, Matrix_jSym, j1);
	pj1 = INTEGER(j1);
	++nprotect;
    }
    if (clf[0] != 'n') {
	PROTECT(x1 = allocVector(TYPEOF(x0), nnz1));
	SET_SLOT(to, Matrix_xSym, x1);
	++nprotect;
    }

#define SPARSE_MAKE_BANDED(_XSET0_, _XSET1_, _XSET2_)			\
    do {								\
	if (clf[2] == 'T') {						\
	    if (!sy && clf[1] == 's') {					\
		for (k = 0; k < nnz0; ++k) {				\
		    if ((d = pj0[k] - pi0[k]) >= a && d <= b) {		\
			*(pi1++) = pi0[k];				\
			*(pj1++) = pj0[k];				\
			_XSET0_; /* *(px1++) = px0[k]; */		\
		    }							\
		    if (d != 0 && -d >= a && -d <= b) {			\
			*(pi1++) = pj0[k];				\
			*(pj1++) = pi0[k];				\
			_XSET0_; /* *(px1++) = px0[k]; */		\
		    }							\
		}							\
	    } else {							\
		for (k = 0; k < nnz0; ++k) {				\
		    if ((d = pj0[k] - pi0[k]) >= a && d <= b) {		\
			*(pi1++) = pi0[k];				\
			*(pj1++) = pj0[k];				\
			_XSET0_; /* *(px1++) = px0[k]; */		\
		    }							\
		}							\
	    }								\
	} else {							\
	    k = 0;							\
	    if (!sy && clf[1] == 's') {					\
		int *pp1_;						\
		Calloc_or_Alloca_TO(pp1_, n, int);			\
		Memcpy(pp1_, pp1 - 1, n);				\
		for (j = 0; j < n; ++j) {				\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			if ((d = j - pi0[k]) >= a && d <= b) {		\
			    pi1[pp1_[j]] = pi0[k];			\
			    _XSET1_; /* px1[pp1_[j]] = px0[k]; */	\
			    ++pp1_[j];					\
			}						\
			if (d != 0 && -d >= a && -d <= b) {		\
			    pi1[pp1_[pi0[k]]] = j;			\
			    _XSET2_; /* px1[pp1_[pi0[k]]] = px0[k]; */	\
			    ++pp1_[pi0[k]];				\
			}						\
			++k;						\
		    }							\
		}							\
		Free_FROM(pp1_, n);					\
	    } else {							\
		for (j = 0; j < n; ++j) {				\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			if ((d = j - pi0[k]) >= a && d <= b) {		\
			    *(pi1++) = pi0[k];				\
			    _XSET0_; /* *(px1++) = px0[k]; */		\
			}						\
			++k;						\
		    }							\
		}							\
	    }								\
									\
	}								\
    } while (0)

#define SPARSE_MAKE_BANDED_X(_CTYPE_, _PTR_)		\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	SPARSE_MAKE_BANDED(*(px1++) = px0[k],		\
			   px1[pp1_[j]] = px0[k],	\
			   px1[pp1_[pi0[k]]] = px0[k]);	\
    } while (0)

#define SPARSE_MAKE_BANDED_N SPARSE_MAKE_BANDED(, , )

    switch (clf[0]) {
    case 'd':
	SPARSE_MAKE_BANDED_X(double, REAL);
	break;
    case 'l':
	SPARSE_MAKE_BANDED_X(int, LOGICAL);
	break;
    case 'n':
	SPARSE_MAKE_BANDED_N;
	break;
    case 'i':
	SPARSE_MAKE_BANDED_X(int, INTEGER);
	break;
    case 'z':
    	SPARSE_MAKE_BANDED_X(Rcomplex, COMPLEX);
	break;
    default:
	break;
    }

#undef SPARSE_MAKE_BANDED_X
#undef SPARSE_MAKE_BANDED_N
#undef SPARSE_MAKE_BANDED
    
    if (clf[2] == 'R')
	to = tCRsparse_as_RCsparse(to);
    UNPROTECT(nprotect);
    return to;
}

SEXP R_sparse_transpose(SEXP from)
{
    static const char *valid[] = {
	"dgCMatrix", "dtCMatrix", "dsCMatrix",
	"dgRMatrix", "dtRMatrix", "dsRMatrix",
	"dgTMatrix", "dtTMatrix", "dsTMatrix",
	"lgCMatrix", "ltCMatrix", "lsCMatrix",
	"lgRMatrix", "ltRMatrix", "lsRMatrix",
	"lgTMatrix", "ltTMatrix", "lsTMatrix",
	"ngCMatrix", "ntCMatrix", "nsCMatrix",
	"ngRMatrix", "ntRMatrix", "nsRMatrix",
	"ngTMatrix", "ntTMatrix", "nsTMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_transpose");
    const char *cl = valid[ivalid];

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

    if (m == n) {
	SET_SLOT(to, Matrix_DimSym, dim);
    } else {
	pdim = INTEGER(GET_SLOT(to, Matrix_DimSym));
	pdim[0] = n;
	pdim[1] = m;
    }
    if (cl[1] == 's')
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    else
	set_reversed_DimNames(to, dimnames);
    if (cl[1] != 'g') {
	SET_SLOT(to, Matrix_uploSym,
		 mkString((*uplo_P(from) == 'U') ? "L" : "U"));
	if (cl[1] == 't')
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    }

    /* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

    if (cl[2] == 'T') {
	SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_jSym));
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_iSym));
	if (cl[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	UNPROTECT(1);
	return to;
    }

    /* Now dealing only with [CR]sparseMatrix ... */

    int m_, n_;
    SEXP isym;
    if (cl[2] == 'C') {
	m_ = m;
	n_ = n;
	isym = Matrix_iSym;
    } else {
	m_ = n;
	n_ = m;
	isym = Matrix_jSym;
    }
    
    R_xlen_t m1a = (R_xlen_t) m_ + 1;
    SEXP 
	p0 = GET_SLOT(from, Matrix_pSym),
	p1 = PROTECT(allocVector(INTSXP, m1a));
    int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), nnz = pp0[n_]; ++pp0;
    SEXP x0, x1,
	i0 = GET_SLOT(from, isym),
	i1 = PROTECT(allocVector(INTSXP, nnz));
    int *pi0 = INTEGER(i0), *pi1 = INTEGER(i1), i, j, k, kend;

    SET_SLOT(to, Matrix_pSym, p1);
    SET_SLOT(to, isym, i1);
    if (cl[0] != 'n') {
	x0 = GET_SLOT(from, Matrix_xSym);
	PROTECT(x1 = allocVector(TYPEOF(x0), nnz));
	SET_SLOT(to, Matrix_xSym, x1);
    }
    
    /* Count number of nonzero elements, by row */
    Memzero(pp1, m1a); ++pp1;
    for (k = 0; k < nnz; ++k)
	++pp1[pi0[k]];

    /* Compute cumulative sum, in place */
    for (i = 1; i < m_; ++i)
	pp1[i] += pp1[i-1];

    /* Allocate work space */
    int *pp1_;
    Calloc_or_Alloca_TO(pp1_, m_, int);
    Memcpy(pp1_, pp1 - 1, m_);

#define SPARSE_T(_XSET_)				\
    do {						\
	k = 0;						\
	for (j = 0; j < n_; ++j) {			\
	    kend = pp0[j];				\
	    while (k < kend) {				\
		i = pi0[k];				\
		pi1[pp1_[i]] = j;			\
		_XSET_;	/* px1[pp1_[i]] = px0[k] */	\
		++pp1_[i];				\
		++k;					\
	    }						\
	}						\
    } while (0)

#define SPARSE_T_X(_CTYPE_, _PTR_)			\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	SPARSE_T(px1[pp1_[i]] = px0[k]);		\
    } while (0);

#define SPARSE_T_N SPARSE_T()

    /* Do transpose */
    switch (cl[0]) {
    case 'd':
	SPARSE_T_X(double, REAL);
	break;
    case 'l':
	SPARSE_T_X(int, LOGICAL);
	break;
    case 'n':
	SPARSE_T_N;
	break;
    case 'i':
	SPARSE_T_X(int, INTEGER);
	break;
    case 'z':
    	SPARSE_T_X(Rcomplex, COMPLEX);
	break;
    default:
	break;
    }

    Free_FROM(pp1_, m_);
    UNPROTECT((cl[0] == 'n') ? 3 : 4);
    return to;
}

SEXP tCRsparse_as_RCsparse(SEXP from)
{
    static const char *valid[] = {
	"dgCMatrix", "dtCMatrix", "dsCMatrix",
	"dgRMatrix", "dtRMatrix", "dsRMatrix",
	"lgCMatrix", "ltCMatrix", "lsCMatrix",
	"lgRMatrix", "ltRMatrix", "lsRMatrix",
	"ngCMatrix", "ntCMatrix", "nsCMatrix",
	"ngRMatrix", "ntRMatrix", "nsRMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "tCRsparse_as_RCsparse");
    
    char *cl = strdup(valid[ivalid]);
    cl[2] = (cl[2] == 'C') ? 'R' : 'C';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    
    if (m == n)	{
	SET_SLOT(to, Matrix_DimSym, dim);
    } else {
	pdim = INTEGER(GET_SLOT(to, Matrix_DimSym));
	pdim[0] = n;
	pdim[1] = m;
    }
    if (cl[1] == 's')
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    else
	set_reversed_DimNames(to, dimnames);
    SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
    if (cl[2] == 'C')
	SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_jSym));
    else
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_iSym));
    if (cl[0] != 'n')
	SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
    if (cl[1] != 'g') {
	SET_SLOT(to, Matrix_uploSym,
		 mkString((*uplo_P(from) == 'U') ? "L" : "U"));
	if (cl[1] == 't')
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    }

    free(cl);
    UNPROTECT(1);
    return to;
}
