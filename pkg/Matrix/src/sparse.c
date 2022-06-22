#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "sparse.h"

SEXP sparse_as_dense(SEXP from, int packed)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_as_dense");
    const char *clf = valid[ivalid];
    packed = packed && clf[1] != 'g';

    char *clt = strdup(clf);
    clt[2] = (clf[1] == 'g')
	? 'e' : ((packed) ? 'p' : ((clf[1] == 't') ? 'r' : 'y'));
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym),
	x0 = (clf[0] != 'n') ? GET_SLOT(from, Matrix_xSym) : R_NilValue;
    free(clt);

    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    double mn = (double) m * n, nx_ = (packed) ? (mn + n) / 2.0 : mn;
    if (nx_ > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
    if (packed && clf[2] != 'C' && mn > R_XLEN_T_MAX)
	error(_("coercing n-by-n [RT]sparseMatrix to packedMatrix "
		"is not supported for n*n exceeding R_XLEN_T_MAX"));
    
    char ul = 'U', di = 'N';
    SEXPTYPE tx = kind2type(clf[0]);
    R_xlen_t nx = (R_xlen_t) nx_;
    SEXP x1 = PROTECT(allocVector(tx, nx));
    
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    SET_SLOT(to, Matrix_xSym, x1);
    if (clf[1] != 'g') {
	SEXP uplo = GET_SLOT(from, Matrix_uploSym);
	SET_SLOT(to, Matrix_uploSym, uplo);
	ul = *CHAR(STRING_ELT(uplo, 0));
	if (clf[1] == 't') {
	    SEXP diag = GET_SLOT(from, Matrix_diagSym);
	    SET_SLOT(to, Matrix_diagSym, diag);
	    di = *CHAR(STRING_ELT(diag, 0));
	}
    }

#define SAD_CASES(_KIND_, _MACRO_)				\
    switch (_KIND_) {						\
    case 'd':							\
    {								\
	double *px0 = REAL(x0), *px1 = REAL(x1);		\
	_MACRO_(px1, *(px0++));					\
	break;							\
    }								\
    case 'l':							\
    {								\
	int *px0 = LOGICAL(x0), *px1 = LOGICAL(x1);		\
	_MACRO_(px1, *(px0++));					\
	break;							\
    }								\
    case 'n':							\
    {								\
	int *px1 = LOGICAL(x1);					\
	_MACRO_(px1, 1);					\
	break;							\
    }								\
    case 'i':							\
    {								\
	int *px0 = INTEGER(x0), *px1 = INTEGER(x1);		\
	_MACRO_(px1, *(px0++));					\
	break;							\
    }								\
    case 'z':							\
    {								\
	Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);	\
	_MACRO_(px1, *(px0++));					\
	break;							\
    }								\
    default:							\
	ERROR_INVALID_TYPE("'x' slot", tx, "sparse_as_dense");	\
	break;							\
    }

    /* It remains to fill 'x' ... */
    
    if (clf[2] == 'C') {
	
	SEXP p = GET_SLOT(from, Matrix_pSym),
	    i = GET_SLOT(from, Matrix_iSym);
	int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
	++pp;
	
#define SAD_C(_X_, _VAL_)						\
	do {								\
	    Memzero(_X_, nx);						\
	    if (!packed) {						\
		for (j = 0, k = 0; j < n; ++j, _X_ += m) {		\
		    kend = pp[j];					\
		    while (k < kend) {					\
			_X_[*pi] = _VAL_;				\
			++pi; ++k;					\
		    }							\
		}							\
	    } else if (ul == 'U') {					\
		for (j = 0, k = 0; j < n; _X_ += (++j)) {		\
		    kend = pp[j];					\
		    while (k < kend) {					\
			_X_[*pi] = _VAL_;				\
			++pi; ++k;					\
		    }							\
		}							\
	    } else {							\
		for (j = 0, k = 0; j < n; _X_ += n-(j++)) {		\
		    kend = pp[j];					\
		    while (k < kend) {					\
			_X_[*pi - j] = _VAL_;				\
			++pi; ++k;					\
		    }							\
		}							\
	    }								\
	} while (0)

	SAD_CASES(clf[0], SAD_C);

#undef SAD_C
	
    } else if (clf[2] == 'R') {
	
	SEXP p = GET_SLOT(from, Matrix_pSym),
	    j = GET_SLOT(from, Matrix_jSym);
	int *pp = INTEGER(p), *pj = INTEGER(j), i, k, kend;
	++pp;
	
#define SAD_R(_X_, _VAL_)						\
	do {								\
	    Memzero(_X_, nx);						\
	    if (!packed) {						\
		R_xlen_t m1 = (R_xlen_t) m;				\
		for (i = 0, k = 0; i < m; ++i, ++_X_) {			\
		    kend = pp[i];					\
		    while (k < kend) {					\
			_X_[m1 * *pj] = _VAL_;				\
			++pj; ++k;					\
		    }							\
		}							\
	    } else if (ul == 'U') {					\
		for (i = 0, k = 0; i < n; ++i) {			\
		    kend = pp[i];					\
		    while (k < kend) {					\
			_X_[PM_AR21_UP(i, *pj)] = _VAL_;		\
			++pj; ++k;					\
		    }							\
		}							\
	    } else {							\
		R_xlen_t n2 = (R_xlen_t) n * 2;				\
		for (i = 0, k = 0; i < n; ++i) {			\
		    kend = pp[i];					\
		    while (k < kend) {					\
			_X_[PM_AR21_LO(i, *pj, n2)] = _VAL_;		\
			++pj; ++k;					\
		    }							\
		}							\
	    }								\
	} while (0)
	
	SAD_CASES(clf[0], SAD_R);

#undef SAD_R
#undef SAD_CASES
	
    } else {
	
	SEXP i = GET_SLOT(from, Matrix_iSym),
	    j = GET_SLOT(from, Matrix_jSym);
	int *pi = INTEGER(i), *pj = INTEGER(j), nnz = LENGTH(i), k;

#define SAD_SUBCASES(_DO_LOOP_)				\
	do {						\
	    if (!packed) {				\
		R_xlen_t m1 = (R_xlen_t) m;		\
		_DO_LOOP_(m1 * *pj + *pi);		\
	    } else if (ul == 'U') {			\
		_DO_LOOP_(PM_AR21_UP(*pi, *pj));	\
	    } else {					\
		R_xlen_t n2 = (R_xlen_t) n * 2;		\
		_DO_LOOP_(PM_AR21_LO(*pi, *pj, n2));	\
	    }						\
	} while (0)
		
	switch (clf[0]) {
	case 'd':
	{
	    double *px0 = REAL(x0), *px1 = REAL(x1);
	    Memzero(px1, nx);

#define SAD_T_D(_INDEX_)					\
	    do {						\
		for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0)	\
		    px1[_INDEX_] += *px0;			\
	    } while (0)
	    
	    SAD_SUBCASES(SAD_T_D);
	    break;
	}
	case 'l':
	{
	    int *px0 = LOGICAL(x0), *px1 = LOGICAL(x1);
	    R_xlen_t index;
	    Memzero(px1, nx);

#define SAD_T_L(_INDEX_)						\
	    do {							\
		for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) {		\
		    if (*px0 != 0) {					\
			index = _INDEX_;				\
			if (*px0 == NA_LOGICAL) {			\
			    if (px1[index] == 0)			\
				px1[index] = NA_LOGICAL;		\
			} else {					\
			    px1[index] = 1;				\
			}						\
		    }							\
		}							\
	    } while (0)

	    SAD_SUBCASES(SAD_T_L);
	    break;
	}
	case 'n':
	{
	    int *px1 = LOGICAL(x1);
	    Memzero(px1, nx);
	    
#define SAD_T_N(_INDEX_)						\
	    do {							\
		for (k = 0; k < nnz; ++k, ++pi, ++pj)			\
		    px1[_INDEX_] = 1;					\
	    } while (0)

	    SAD_SUBCASES(SAD_T_N);
	    break;
	}
	case 'i':
	{
	    int *px0 = INTEGER(x0), *px1 = INTEGER(x1);
	    Memzero(px1, nx);

/* FIXME: not detecting integer overflow here */
#define SAD_T_I(_INDEX_) SAD_T_D(_INDEX_)

	    SAD_SUBCASES(SAD_T_I);
	    break;
	}
	case 'z':
	{
	    Rcomplex *px0 = COMPLEX(x0), *px1 = COMPLEX(x1);
	    Memzero(px1, nx);
	    R_xlen_t index;

#define SAD_T_Z(_INDEX_)					\
	    do {						\
		for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px0) {	\
		    index = _INDEX_;				\
		    px1[index].r += (*px0).r;			\
		    px1[index].i += (*px0).i;			\
		}						\
	    } while (0)
	    
	    SAD_SUBCASES(SAD_T_Z);
	    break;
	}
	default:
	    ERROR_INVALID_TYPE("'x' slot", tx, "sparse_as_dense");
	    break;
	}

#undef SAD_T_D
#undef SAD_T_L
#undef SAD_T_N
#undef SAD_T_I
#undef SAD_T_Z
#undef SAD_SUBCASES
    
    }
    
    if (di != 'N') {

	int j;
	
#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {						\
	    _CTYPE_ *px1 = _PTR_(x1);			\
	    if (!packed) {				\
		R_xlen_t n1a = (R_xlen_t) n + 1;	\
		for (j = 0; j < n; ++j, px1 += n1a)	\
		    *px1 = _ONE_;			\
	    } else if (ul == 'U') {			\
		for (j = 0; j < n; px1 += (++j)+1)	\
		    *px1 = _ONE_;			\
	    } else {					\
		for (j = 0; j < n; px1 += n-(j++))	\
		    *px1 = _ONE_;			\
	    }						\
	} while (0)
	    
	SPARSE_CASES(tx, SET1);

#undef SET1
	
    }

    UNPROTECT(2);
    return to;
}

/* as(<[CRT]sparseMatrix>,       "denseMatrix") */
/* as(<[CRT]sparseMatrix>, "(un)?packedMatrix") */
SEXP R_sparse_as_dense(SEXP from, SEXP packed)
{
    return sparse_as_dense(from, asLogical(packed));
}

/* as(<[CRT]sparseMatrix>, "matrix") */
SEXP R_sparse_as_matrix(SEXP from)
{
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(from = sparse_as_dense(from, 0), &pid);
    REPROTECT(from = dense_as_general(from, '.', 0, 0), pid);
    SEXP to = PROTECT(GET_SLOT(from, Matrix_xSym));
    setAttrib(to, R_DimSymbol, GET_SLOT(from, Matrix_DimSym));
    setAttrib(to, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(2);
    return to;
}

/* as(<[CRT]sparseMatrix>, "vector") */
SEXP R_sparse_as_vector(SEXP from)
{
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(from = sparse_as_dense(from, 0), &pid);
    REPROTECT(from = dense_as_general(from, '.', 0, 0), pid);
    from = GET_SLOT(from, Matrix_xSym);
    UNPROTECT(1);
    return from;
}

/* as(<[CRT]sparseMatrix>, "[dlniz]Matrix") */
SEXP R_sparse_as_kind(SEXP from, SEXP kind, SEXP drop0)
{
    char k;
    if ((kind = asChar(kind)) == NA_STRING || (k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_sparse_as_kind()'"));
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
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

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {						\
	    _CTYPE_ *px = _PTR_(x);			\
	    for (ix = 0; ix < nx; ++ix)			\
		*(px++) = _ONE_;			\
	} while (0)
	
	SPARSE_CASES(tx, SET1);

#undef SET1
	
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(1);
    }

    free(clt);
    UNPROTECT(do_drop0 ? 2 : 1);
    return to;
}

/* as(<[CRT]sparseMatrix>, "generalMatrix") */
SEXP R_sparse_as_general(SEXP from)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_as_general");
    const char *clf = valid[ivalid];
    if (clf[1] == 'g')
	return from;
    
    char *clt = strdup(clf);
    clt[1] = 'g';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    free(clt);
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (clf[1] == 's') {
	set_symmetrized_DimNames(to, dimnames, -1);
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    } else {
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	if (*diag_P(from) == 'N') {
	    if (clf[2] != 'T')
		SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
	    if (clf[2] != 'R')
		SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_iSym));
	    if (clf[2] != 'C')
		SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_jSym));
	    if (clf[0] != 'n')
		SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	    UNPROTECT(1);
	    return to;
	}
    }

    int n = INTEGER(dim)[0];
    char ul = *uplo_P(from);

    if (clf[2] == 'T') {
	
	SEXP i0 = GET_SLOT(from, Matrix_iSym),
	    j0 = GET_SLOT(from, Matrix_jSym);
	int k, *pi0 = INTEGER(i0), *pj0 = INTEGER(j0),
	    nnz0 = LENGTH(i0), nnz1 = nnz0;
	
	if (clf[1] == 's') {
	    for (k = 0; k < nnz0; ++k)
		if (pi0[k] == pj0[k])
		    --nnz1;
	    nnz1 += nnz0;
	} else {
	    nnz1 += n;
	}

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_iSym, i1);
	SET_SLOT(to, Matrix_jSym, j1);
	Memcpy(pi1, pi0, nnz0);
	Memcpy(pj1, pj0, nnz0);
	pi1 += nnz0;
	pj1 += nnz0;

	SEXP x0, x1;
	SEXPTYPE tx;
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}
	
	if (clf[1] == 's') {

#define SAG_T(_XASSIGN_, _XINCR_)				\
	    do {						\
		for (k = 0; k < nnz0; ++k) {			\
		    if (*pi0 != *pj0) {				\
			*(pi1++) = *pj0;			\
			*(pj1++) = *pi0;			\
			_XASSIGN_; /* *(px1++) = *px0; */	\
		    }						\
		    ++pi0; ++pj0; _XINCR_; /* ++px0; */		\
		}						\
	    } while (0)
	    
#define SAG_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		Memcpy(px1, px0, nnz0);				\
		px1 += nnz0;					\
		SAG_T(*(px1++) = *px0, ++px0);			\
	    } while (0)

	    if (clf[0] == 'n')
		SAG_T(, );
	    else
		SPARSE_CASES(tx, SAG_T_X);

#undef SAG_T_X
#undef SAG_T

	} else {

#define SAG_T(_XASSIGN_ONE_)					\
	    do {						\
		for (k = 0; k < n; ++k) {			\
		    *(pi1++) = *(pj1++) = k;			\
		    _XASSIGN_ONE_; /* *(px1++) = _ONE_; */	\
		}						\
	    } while (0)

#define SAG_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		Memcpy(px1, px0, nnz0);				\
		px1 += nnz0;					\
		SAG_T(*(px1++) = _ONE_);			\
	    } while (0)
	    
	    if (clf[0] == 'n')
		SAG_T();
	    else
		SPARSE_CASES(tx, SAG_T_X);

#undef SAG_T_X
#undef SAG_T
	    
	}
	
    } else {

	SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
	    p0 = GET_SLOT(from, Matrix_pSym),
	    p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
	    i0 = GET_SLOT(from, iSym);
	int j, k, kend, nnz1,
	    *pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0);
	pp0++;
	*(pp1++) = 0;
	
	if (clf[1] == 's') {
	    Memzero(pp1, n);
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if (pi0[k] != j)
			++pp1[pi0[k]];
		    ++k;
		}
	    }
	    for (j = 1; j < n; ++j)
		pp1[j] += pp1[j-1];
	    for (j = 0; j < n; ++j)
		pp1[j] += pp0[j];
	} else {
	    for (j = 0; j < n; ++j)
		pp1[j] = pp0[j] + j + 1;
	}

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz1 = pp1[n-1]));
	int *pi1 = INTEGER(i1);
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to, iSym, i1);

	SEXP x0, x1;
	SEXPTYPE tx;
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}

	if (clf[1] == 's') {

	    int *pp1_;
	    Calloc_or_Alloca_TO(pp1_, n, int);
	    Memcpy(pp1_, pp1 - 1, n);
	    
#define SAG_CR(_XASSIGN_IJ_, _XASSIGN_JI_)				\
	    do {							\
		for (j = 0, k = 0; j < n; ++j) {			\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			pi1[pp1_[j]] = pi0[k];				\
			_XASSIGN_IJ_; /* px1[pp1_[j]] = px0[k]; */	\
			++pp1_[j];					\
			if (pi0[k] != j) {				\
			    pi1[pp1_[pi0[k]]] = j;			\
			    _XASSIGN_JI_; /* px1[pp1_[pi0[k]]] = px0[k]; */ \
			    ++pp1_[pi0[k]];				\
			}						\
			++k;						\
		    }							\
		}							\
	    } while (0)
	    
#define SAG_CR_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		SAG_CR(px1[pp1_[j]] = px0[k],			\
		       px1[pp1_[pi0[k]]] = px0[k]);		\
	    } while (0)

	    if (clf[0] == 'n')
		SAG_CR(, );
	    else
		SPARSE_CASES(tx, SAG_CR_X);

#undef SAG_CR_X
#undef SAG_CR

	    Free_FROM(pp1_, n);
	    
	} else {

#define SAG_CR(_XASSIGN_, _XASSIGN_ONE_)				\
	    do {							\
		if (ul == ((clf[2] == 'C') ? 'U' : 'L')) {		\
		    for (j = 0, k = 0; j < n; ++j) {			\
			kend = pp0[j];					\
			while (k < kend) {				\
			    *(pi1++) = *(pi0++);			\
			    _XASSIGN_; /* *(px1++) = *(px0++); */	\
			    ++k;					\
			}						\
			*(pi1++) = j;					\
			_XASSIGN_ONE_; /* *(px1++) = _ONE_; */		\
		    }							\
		} else {						\
		    for (j = 0, k = 0; j < n; ++j) {			\
			kend = pp0[j];					\
			*(pi1++) = j;					\
			_XASSIGN_ONE_; /* *(px1++) = _ONE_; */		\
			while (k < kend) {				\
			    *(pi1++) = *(pi0++);			\
			    _XASSIGN_; /* *(px1++) = *(px0++); */	\
			    ++k;					\
			}						\
		    }							\
		}							\
	    } while (0)

#define SAG_CR_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		SAG_CR(*(px1++) = *(px0++), *(px1++) = _ONE_);	\
	    } while (0)

	    if (clf[0] == 'n')
		SAG_CR(, );
	    else
		SPARSE_CASES(tx, SAG_CR_X);

#undef SAG_CR_X
#undef SAG_CR

	}
	
    }

    UNPROTECT((clf[0] == 'n') ? 3 : 4);
    return to;
}

/* as(<diagonalMatrix>, ".(di|[gts][CRT])Matrix") */
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
    
    static const char *valid[] = { VALID_DIAGONAL, "" };
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

#define SET1(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		do {						\
		    _CTYPE_ *px = _PTR_(x);			\
		    for (k = 0; k < n; ++k) {			\
			*(pi++) = k;				\
			*(px++) = _ONE_;			\
		    }						\
		} while (0)
		
		SPARSE_CASES(tx, SET1);

#undef SET1

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

#undef DROP0_DIAGONAL
	    
	} else { /* _do not_ drop zeros */
	    
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

/* as(<diagonalMatrix>, ".(ge|tr|sy|tp|sp)Matrix") */
SEXP R_diagonal_as_dense(SEXP from, SEXP code, SEXP uplo)
{
    const char *zzz;
    char z0, z1, z2, ul;
    if ((code = asChar(code)) == NA_STRING ||
	(z0 = (zzz = CHAR(code))[0]) == '\0' ||
	(z1 = zzz[1]) == '\0' ||
	(z2 = zzz[2]) == '\0')
	error(_("invalid 'code' to 'R_diagonal_as_dense()'"));
    if ((z1 == 't' || z1 == 's') &&
	(((uplo = asChar(uplo)) == NA_STRING || (ul = *CHAR(uplo)) == '\0')))
	error(_("invalid 'uplo' to 'R_diagonal_as_dense()'"));
    
    static const char *valid[] = { VALID_DIAGONAL, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_diagonal_as_dense");

    const char *clf = valid[ivalid];
    if (z0 == '.')
	z0 = clf[0];
    
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
    int n = INTEGER(dim)[0];
    double nn = (double) n * n, nx_ = (z2 == 'p') ? (nn + n) / 2.0 : nn;
    if (nx_ > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
    R_xlen_t nx = (R_xlen_t) nx_;
    SEXP x = PROTECT(allocVector(tx, nx));
    
    SET_SLOT(to, Matrix_DimSym, dim);
    if (z1 == 's')
	set_symmetrized_DimNames(to, dimnames, -1);
    else
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (z1 == 't' || z1 == 's')
	SET_SLOT(to, Matrix_uploSym, mkString((ul == 'U') ? "U" : "L"));
    if (z1 == 't')
	SET_SLOT(to, Matrix_diagSym, diag);
    SET_SLOT(to, Matrix_xSym, x);
    
#define DAD_COPY_DIAGONAL(_PREFIX_, _PTR_)				\
    do {								\
	Memzero(_PTR_(x), nx);						\
	if (z1 != 't' || di == 'N') {					\
	    SEXP y = PROTECT(coerceVector(GET_SLOT(from, Matrix_xSym), tx)); \
	    if (z2 == 'p')						\
		_PREFIX_ ## dense_packed_copy_diagonal(			\
		    _PTR_(x), _PTR_(y), n, n, ul, ul /* unused */, di);	\
	    else							\
		_PREFIX_ ## dense_unpacked_copy_diagonal(		\
		    _PTR_(x), _PTR_(y), n, n,     ul /* unused */, di);	\
	    UNPROTECT(1);						\
	}								\
    } while (0)
	
    switch (tx) {
    case REALSXP:
	DAD_COPY_DIAGONAL(d, REAL);
	break;
    case LGLSXP:
	DAD_COPY_DIAGONAL(i, LOGICAL);
	break;
    case INTSXP:
	DAD_COPY_DIAGONAL(i, INTEGER);
	break;
    case CPLXSXP:
	DAD_COPY_DIAGONAL(z, COMPLEX);
	break;
    default:
	break;
    }

    UNPROTECT(2);
    return to;
}

/* drop0(<[CRT]sparseMatrix>) 
   TODO: support 'tol' argument, to be interpreted as modulus for zMatrix
*/
SEXP R_sparse_drop0(SEXP from)
{
    static const char *valid[] = { VALID_DSPARSE, VALID_LSPARSE, "" };
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

    SEXP iSym = (cl[2] == 'R') ? Matrix_jSym : Matrix_iSym,
	i0 = GET_SLOT(from, iSym),
	i1 = PROTECT(allocVector(INTSXP, nnz1)),
	x1 = PROTECT(allocVector(tx, nnz1));
    int *pi0 = INTEGER(i0),
	*pi1 = INTEGER(i1);
    SET_SLOT(to, iSym, i1);
    SET_SLOT(to, Matrix_xSym, x1);
    
    if (cl[2] == 'T') {

	SEXP j0 = GET_SLOT(from, Matrix_jSym),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pj0 = INTEGER(j0),
	    *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_jSym, j1);

#define DROP0_END(_CTYPE_, _PTR_, _NZ_)				\
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

	DROP0_CASES(tx, DROP0_END);

#undef DROP0_END
	
    } else { /* cl[2] == 'C' || cl[2] == 'R' */

	SEXP p1 = PROTECT(allocVector(INTSXP, n1a));
	int *pp1 = INTEGER(p1), j, kend;
	SET_SLOT(to, Matrix_pSym, p1);

#define DROP0_END(_CTYPE_, _PTR_, _NZ_)			\
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
    
	DROP0_CASES(tx, DROP0_END);

#undef DROP0_END
	
    }

    UNPROTECT(4);
    return to;
}

/* band(<[CRT]sparseMatrix>, k1, k2), tri[ul](<[CRT]sparseMatrix>, k) */
/* NB: argument validation more or less copied from R_dense_band() */
SEXP R_sparse_band(SEXP from, SEXP k1, SEXP k2)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
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
    
    SEXP p0, p1, i0, j0;
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
	pp0++;
    }
    
    /* Counting number of nonzero elements in band ... */
    
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
	
	if (!sy && clf[1] == 's') {
	    Memzero(pp1, n);
	    for (j = 0, k = 0; j < n; ++j) {
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
	    for (j = 0, k = 0; j < n; ++j) {
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

    if (nnz1 == nnz0 && (sy || clf[1] != 's')) {
	/* No need to allocate in this case: band has all nonzero elements */
	SET_SLOT(to, Matrix_iSym, i0);
	if (clf[2] == 'T')
	    SET_SLOT(to, Matrix_jSym, j0);
	if (clf[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));	
	if (clf[2] == 'R')
	    to = tCRsparse_as_RCsparse(to);
	UNPROTECT(nprotect);
	return to;
    }
    
    /* Now allocating and filling out slots ... */

    SEXP i1, j1;
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

#define SPARSE_BAND(_XASSIGN_, _XASSIGN_IJ_, _XASSIGN_JI_)		\
    do {								\
	if (clf[2] == 'T') {						\
	    if (!sy && clf[1] == 's') {					\
		for (k = 0; k < nnz0; ++k) {				\
		    if ((d = pj0[k] - pi0[k]) >= a && d <= b) {		\
			*(pi1++) = pi0[k];				\
			*(pj1++) = pj0[k];				\
			_XASSIGN_; /* *(px1++) = px0[k]; */		\
		    }							\
		    if (d != 0 && -d >= a && -d <= b) {			\
			*(pi1++) = pj0[k];				\
			*(pj1++) = pi0[k];				\
			_XASSIGN_; /* *(px1++) = px0[k]; */		\
		    }							\
		}							\
	    } else {							\
		for (k = 0; k < nnz0; ++k) {				\
		    if ((d = pj0[k] - pi0[k]) >= a && d <= b) {		\
			*(pi1++) = pi0[k];				\
			*(pj1++) = pj0[k];				\
			_XASSIGN_; /* *(px1++) = px0[k]; */		\
		    }							\
		}							\
	    }								\
	} else {							\
	    if (!sy && clf[1] == 's') {					\
		int *pp1_;						\
		Calloc_or_Alloca_TO(pp1_, n, int);			\
		Memcpy(pp1_, pp1 - 1, n);				\
		for (j = 0, k = 0; j < n; ++j) {			\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			if ((d = j - pi0[k]) >= a && d <= b) {		\
			    pi1[pp1_[j]] = pi0[k];			\
			    _XASSIGN_IJ_; /* px1[pp1_[j]] = px0[k]; */	\
			    ++pp1_[j];					\
			}						\
			if (d != 0 && -d >= a && -d <= b) {		\
			    pi1[pp1_[pi0[k]]] = j;			\
			    _XASSIGN_JI_; /* px1[pp1_[pi0[k]]] = px0[k]; */ \
			    ++pp1_[pi0[k]];				\
			}						\
			++k;						\
		    }							\
		}							\
		Free_FROM(pp1_, n);					\
	    } else {							\
		for (j = 0, k = 0; j < n; ++j) {			\
		    kend = pp0[j];					\
		    while (k < kend) {					\
			if ((d = j - pi0[k]) >= a && d <= b) {		\
			    *(pi1++) = pi0[k];				\
			    _XASSIGN_; /* *(px1++) = px0[k]; */		\
			}						\
			++k;						\
		    }							\
		}							\
	    }								\
	}								\
    } while (0)

#define SPARSE_BAND_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	SPARSE_BAND(*(px1++) = px0[k],			\
		    px1[pp1_[j]] = px0[k],		\
		    px1[pp1_[pi0[k]]] = px0[k]);	\
    } while (0)

    if (clf[0] == 'n') {
	SPARSE_BAND(, , );
    } else {
	SEXPTYPE tx;
	SEXP x0 = GET_SLOT(from, Matrix_xSym),
	    x1 = PROTECT(allocVector(tx = TYPEOF(x0), nnz1));
	SPARSE_CASES(tx, SPARSE_BAND_X);
	SET_SLOT(to, Matrix_xSym, x1);
	UNPROTECT(1);
    }

#undef SPARSE_BAND_X
#undef SPARSE_BAND
    
    if (clf[2] == 'R')
	to = tCRsparse_as_RCsparse(to);
    UNPROTECT(nprotect);
    return to;
}

/* diag(<[CRT]sparseMatrix>, names) */
SEXP R_sparse_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
	error(_("'names' must be TRUE or FALSE"));

    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(obj, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(obj), "R_sparse_diag_get");
    const char *cl = valid[ivalid];

    char kind = cl[0];
    SEXPTYPE type = kind2type(kind);
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)),
	m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    SEXP res = PROTECT(allocVector(type, r));

    if (cl[1] == 't' && *diag_P(obj) != 'N') {

	/* .t[CRT]Matrix with unit diagonal */

	int i;
	
#define DO_ONES(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
	do {					\
	    _CTYPE_ *pres = _PTR_(res);		\
	    for (i = 0; i < r; ++i)		\
		*(pres++) = _ONE_;		\
	} while (0)

	SPARSE_CASES(type, DO_ONES);

#undef DO_ONES
	
    } else if (cl[2] == 'T') {

	/* ..TMatrix with non-unit diagonal */
	
	SEXP x = (kind != 'n') ? GET_SLOT(obj, Matrix_xSym) : R_NilValue,
	    i = GET_SLOT(obj, Matrix_iSym),
	    j = GET_SLOT(obj, Matrix_jSym);
	int *pi = INTEGER(i), *pj = INTEGER(j);
	R_xlen_t k, nnz = XLENGTH(i);

	switch (kind) {
	case 'd':
	{
	    double *px = REAL(x), *pres = REAL(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px)
		if (*pi == *pj)
		    pres[*pi] += *px;
	    break;
	}
	case 'l':
	{
	    int *px = LOGICAL(x), *pres = LOGICAL(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px) {
		if (*pi == *pj && *px != 0) {
		    if (*px == NA_LOGICAL) {
			if (pres[*pi] == 0)
			    pres[*pi] = NA_LOGICAL;
		    } else {
			pres[*pi] = 1;
		    }
		}
	    }
	    break;
	}
	case 'n':
	{
	    int *pres = LOGICAL(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj)
		if (*pi == *pj)
		    pres[*pi] = 1;
	    break;
	}
	case 'i':
	{
	    /* FIXME: not detecting integer overflow here */
	    int *px = INTEGER(x), *pres = INTEGER(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px)
		if (*pi == *pj)
		    pres[*pi] += *px;
	    break;
	}
	case 'z':
	{
	    Rcomplex *px = COMPLEX(x), *pres = COMPLEX(res);
	    Memzero(pres, r);
	    for (k = 0; k < nnz; ++k, ++pi, ++pj, ++px) {
		if (*pi == *pj) {
		    pres[*pi].r += (*px).r;
		    pres[*pi].i += (*px).i;
		}
	    }
	    break;
	}
	default:
	    break;
	}
	
    } else {

	/* ..[CR]Matrix with non-unit diagonal */
	
	SEXP x = (kind != 'n') ? GET_SLOT(obj, Matrix_xSym) : R_NilValue,
	    iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym;
	int j, k, kend,
	    *pp = INTEGER(GET_SLOT(obj, Matrix_pSym)) + 1,
	    *pi = INTEGER(GET_SLOT(obj, iSym));
	char ul = (cl[1] != 'g') ? *uplo_P(obj) : '\0';
	
#define DO_DIAG(_RES_, _VAL_GENERAL_, _VAL_TRAILING_, _VAL_LEADING_, _ZERO_) \
	do {								\
	    if (ul == '\0') {						\
		/* .g[CR]Matrix */					\
	    	for (j = 0, k = 0; j < r; ++j) {			\
		    _RES_[j] = _ZERO_;					\
		    kend = pp[j];					\
		    while (k < kend) {					\
			if (pi[k] == j) {				\
			    _RES_[j] = _VAL_GENERAL_ /* px[k] */;	\
			    k = kend;					\
			    break;					\
			}						\
			++k;						\
		    }							\
		}							\
	    } else if (ul == ((cl[2] == 'C') ? 'U' : 'L')) {		\
		/* .[ts][CR]Matrix with "trailing" diagonal */		\
	    	for (j = 0, k = 0; j < r; ++j) {			\
		    kend = pp[j];					\
		    _RES_[j] = (kend - k > 0 && pi[kend-1] == j		\
				? _VAL_TRAILING_ /* px[kend-1] */	\
				: _ZERO_);				\
		    k = kend;						\
		}							\
	    } else {							\
		/* .[ts][CR]Matrix with "leading" diagonal */		\
		for (j = 0, k = 0; j < r; ++j) {			\
		    kend = pp[j];					\
		    _RES_[j] = (kend - k > 0 && pi[k] == j		\
				? _VAL_LEADING_	/* px[k] */		\
				: _ZERO_);				\
		    k = kend;						\
		}							\
	    }								\
	} while (0)
	
#define DO_DIAG_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)		\
	do {							\
	    _CTYPE_ *px = _PTR_(x), *pres = _PTR_(res);		\
	    DO_DIAG(pres, px[k], px[kend-1], px[k], _ZERO_);	\
	} while (0)
	    
	if (kind == 'n') {
	    int *pres = LOGICAL(res);
	    DO_DIAG(pres, 1, 1, 1, 0);
	} else {
	    SPARSE_CASES(type, DO_DIAG_X);
	}

#undef DO_DIAG_X
#undef DO_DIAG
	
    }

    if (do_nms) {
	/* NB: The logic here must be adjusted once the validity method 
	       for 'symmetricMatrix' enforces symmetric 'Dimnames' */
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	    rn = VECTOR_ELT(dn, 0),
	    cn = VECTOR_ELT(dn, 1);
	if (isNull(cn)) {
	    if (cl[1] == 's' && !isNull(rn))
		setAttrib(res, R_NamesSymbol, rn);
	} else {
	    if (cl[1] == 's')
		setAttrib(res, R_NamesSymbol, cn);
	    else if (!isNull(rn) &&
		     (rn == cn || equal_string_vectors(rn, cn, r)))
		setAttrib(res, R_NamesSymbol, (r == m) ? rn : cn);
	}
    }

    UNPROTECT(1);
    return res;
}

/* t(<[CRT]sparseMatrix>) */
SEXP R_sparse_transpose(SEXP from)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
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
	/* No need to allocate in this case: need only reverse 'i' and 'j' */
	SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_jSym));
	SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_iSym));
	if (cl[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	UNPROTECT(1);
	return to;
    }

    /* Now dealing only with [CR]sparseMatrix ... */

    int m_, n_;
    SEXP iSym;
    if (cl[2] == 'C') {
	m_ = m;
	n_ = n;
	iSym = Matrix_iSym;
    } else {
	m_ = n;
	n_ = m;
	iSym = Matrix_jSym;
    }
    
    R_xlen_t m1a = (R_xlen_t) m_ + 1;
    SEXP p0 = GET_SLOT(from, Matrix_pSym),
	p1 = PROTECT(allocVector(INTSXP, m1a));
    int i, j, k, kend,
	*pp0 = INTEGER(p0),
	*pp1 = INTEGER(p1),
	nnz = pp0[n_];
    SEXP i0 = GET_SLOT(from, iSym),
	i1 = PROTECT(allocVector(INTSXP, nnz));
    int *pi0 = INTEGER(i0),
	*pi1 = INTEGER(i1);
    pp0++;
    
    SET_SLOT(to, Matrix_pSym, p1);
    SET_SLOT(to, iSym, i1);
    
    /* Counting number of nonzero elements, by "row" */
    Memzero(pp1, m1a);
    ++pp1;
    for (k = 0; k < nnz; ++k)
	++pp1[pi0[k]];

    /* Computing cumulative sum, in place */
    for (i = 1; i < m_; ++i)
	pp1[i] += pp1[i-1];

    /* Allocating work space */
    int *pp1_;
    Calloc_or_Alloca_TO(pp1_, m_, int);
    Memcpy(pp1_, pp1 - 1, m_);

#define SPARSE_T(_XASSIGN_)				\
    do {						\
	for (j = 0, k = 0; j < n_; ++j) {		\
	    kend = pp0[j];				\
	    while (k < kend) {				\
		i = pi0[k];				\
		pi1[pp1_[i]] = j;			\
		_XASSIGN_; /* px1[pp1_[i]] = px0[k] */	\
		++pp1_[i];				\
		++k;					\
	    }						\
	}						\
    } while (0)

#define SPARSE_T_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
    do {						\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
	SPARSE_T(px1[pp1_[i]] = px0[k]);		\
    } while (0);

    if (cl[0] == 'n') {
	SPARSE_T();
    } else {
	SEXPTYPE tx;
	SEXP x0 = GET_SLOT(from, Matrix_xSym),
	    x1 = PROTECT(allocVector(tx = TYPEOF(x0), nnz));
	SET_SLOT(to, Matrix_xSym, x1);
	SPARSE_CASES(tx, SPARSE_T_X);
	UNPROTECT(1);
    }
			 
#undef SPARSE_T_X
#undef SPARSE_T

    Free_FROM(pp1_, m_);
    UNPROTECT(3);
    return to;
}

/* forceSymmetric(<[CRT]sparseMatrix>, uplo) */
SEXP R_sparse_force_symmetric(SEXP from, SEXP uplo_to)
{
    static const char *valid[] = {
	VALID_DSPARSE, VALID_LSPARSE, VALID_NSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_sparse_force_symmetric");
    const char *clf = valid[ivalid];
    
    SEXP uplo_from;
    char ulf = 'U', ult = *CHAR(asChar(uplo_to)), di = 'N';
    if (clf[1] != 'g') {
	uplo_from = GET_SLOT(from, Matrix_uploSym);
	ulf = *CHAR(STRING_ELT(uplo_from, 0));
    }
    if (ult == '\0') /* to handle missing(uplo) */
	ult = ulf;
    if (clf[1] == 's') {
	/* .s[CRT]Matrix */
	if (ulf == ult)
	    return from;
	SEXP to = PROTECT(R_sparse_transpose(from));
	if (clf[0] == 'z') {
	    /* Need _conjugate_ transpose */
	    SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	    conjugate(x);
	    UNPROTECT(1);
	}
	UNPROTECT(1);
	return to;
    }
    if (clf[1] == 't')
	di = *diag_P(from);
    
    SEXP dim = GET_SLOT(from, Matrix_DimSym);
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	error(_("attempt to symmetrize a non-square matrix"));

    /* Now handling just square .[gt][CRT]Matrix ... */

    char *clt = strdup(clf);
    clt[1] = 's';
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    free(clt);

    SET_SLOT(to, Matrix_DimSym, dim);
    set_symmetrized_DimNames(to, dimnames, -1);
    SET_SLOT(to, Matrix_uploSym, mkString((ult == 'U') ? "U" : "L"));

    /* It remains to set some subset of 'p', 'i', 'j', and 'x' ... */

    if (clf[1] == 't' && di == 'N' && ulf == ult) {
	
	/* No need to allocate in this case: we have the triangle we want */
	if (clf[2] != 'T')
	    SET_SLOT(to, Matrix_pSym, GET_SLOT(from, Matrix_pSym));
	if (clf[2] != 'R')
	    SET_SLOT(to, Matrix_iSym, GET_SLOT(from, Matrix_iSym));
	if (clf[2] != 'C')
	    SET_SLOT(to, Matrix_jSym, GET_SLOT(from, Matrix_jSym));
	if (clf[0] != 'n')
	    SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
	UNPROTECT(1);
	return to;
	
    } else if (clf[2] == 'T') {

	/* Symmetrizing square .[gt]TMatrix ... */

	SEXP i0 = GET_SLOT(from, Matrix_iSym),
	    j0 = GET_SLOT(from, Matrix_jSym);
	int k, nnz0 = LENGTH(i0), nnz1 = 0,
	    *pi0 = INTEGER(i0),
	    *pj0 = INTEGER(j0);
	
	/* Counting number of nonzero elements in triangle ... */

	if (clf[1] == 't' && di != 'N') {
	    nnz1 = (ulf == ult) ? n + nnz0 : n;
	} else {
	    if (ult == 'U') {
		for (k = 0; k < nnz0; ++k)
		    if (pi0[k] <= pj0[k])
			++nnz1;
	    } else {
		for (k = 0; k < nnz0; ++k)
		    if (pi0[k] >= pj0[k])
			++nnz1;
	    }
	}

	/* Now allocating and filling out slots ... */

	SEXPTYPE tx;
	SEXP x0, x1,
	    i1 = PROTECT(allocVector(INTSXP, nnz1)),
	    j1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pi1 = INTEGER(i1),
	    *pj1 = INTEGER(j1);
	SET_SLOT(to, Matrix_iSym, i1);
	SET_SLOT(to, Matrix_jSym, j1);
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}

	if (clf[1] == 't' && di != 'N') {
	    if (ulf == ult) {
		Memcpy(pi1, pi0, nnz0);
		Memcpy(pj1, pj0, nnz0);
		pi1 += nnz0;
		pj1 += nnz0;
	    }

#define SPARSE_FS(_XASSIGN_)				\
	    do {					\
		for (k = 0; k < n; ++k) {		\
		    *(pi1++) = *(pj1++) = k;		\
		    _XASSIGN_; /* *(px1++) = _ONE_; */	\
		}					\
	    } while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
	    do {							\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
		if (ulf == ult) {					\
		    Memcpy(px1, px0, nnz0);				\
		    px1 += nnz0;					\
		}							\
		SPARSE_FS(*(px1++) = _ONE_);				\
	    } while (0)

	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		
	} else {

#define SPARSE_FS(_XASSIGN_)					\
	    do {						\
		if (ult == 'U') {				\
		    for (k = 0; k < nnz0; ++k) {		\
			if (pi0[k] <= pj0[k]) {			\
			    *(pi1++) = pi0[k];			\
			    *(pj1++) = pj0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
		    }						\
		} else {					\
		    for (k = 0; k < nnz0; ++k) {		\
			if (pi0[k] <= pj0[k]) {			\
			    *(pi1++) = pi0[k];			\
			    *(pj1++) = pj0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
		    }						\
		}						\
	    } while (0)

#define SPARSE_FS_X_BASIC(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
	    do {						\
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
		SPARSE_FS(*(px1++) = px0[k]);			\
	    } while (0)

	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS
	    
	}
	    
    } else {

	/* Symmetrizing square .[gt][CR]Matrix ... */

	SEXP iSym = (clf[2] == 'C') ? Matrix_iSym : Matrix_jSym,
	    p0 = GET_SLOT(from, Matrix_pSym),
	    p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
	    i0 = GET_SLOT(from, iSym);
	int j, k, kend,
	    *pp0 = INTEGER(p0),
	    *pp1 = INTEGER(p1),
	    *pi0 = INTEGER(i0),
	    nnz0 = pp0[n],
	    nnz1 = 0;
	pp0++;
	*(pp1++) = 0;
	
	/* Counting number of nonzero elements in triangle, by "column" ... */

	if (clf[1] == 't') {
	    if (di != 'N') {
		/* Have triangular matrix with unit diagonal */
		if (ulf != ult) {
		    /* Returning identity matrix */
		    for (j = 0; j < n; ++j)
			pp1[j] = ++nnz1;
		} else {
		    /* Returning symmetric matrix with unit diagonal */
		    for (j = 0; j < n; ++j)
			pp1[j] = ++nnz1 + pp0[j];
		    nnz1 += nnz0;
		}
	    } else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
		/* Have triangular matrix with non-unit "trailing" diagonal
		   and returning diagonal part */
		for (j = 0; j < n; ++j) {
		    if (pp0[j-1] < pp0[j] && pi0[pp0[j]-1] == j)
			++nnz1;
		    pp1[j] = nnz1;
		}
	    } else {
		/* Have triangular matrix with non-unit "leading" diagonal
		   and returning diagonal part */
		for (j = 0; j < n; ++j) {
		    if (pp0[j-1] < pp0[j] && pi0[pp0[j-1]] == j)
			++nnz1;
		    pp1[j] = nnz1;
		}
	    }
	} else if (ult == ((clf[2] == 'C') ? 'U' : 'L')) {
	    /* Have general matrix and returning upper triangle */
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if (pi0[k] <= j)
			++nnz1;
		    ++k;
		}
		pp1[j] = nnz1;
	    }
	} else {
	    /* Have general matrix and returning lower triangle */
	    for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend) {
		    if (pi0[k] >= j)
			++nnz1;
		    ++k;
		}
		pp1[j] = nnz1;
	    }
	}

	/* Now allocating and filling out slots ... */

	SEXPTYPE tx;
	SEXP x0, x1, i1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pi1 = INTEGER(i1);
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to, iSym, i1);
	if (clf[0] != 'n') {
	    x0 = GET_SLOT(from, Matrix_xSym);
	    PROTECT(x1 = allocVector(tx = TYPEOF(x0), nnz1));
	    SET_SLOT(to, Matrix_xSym, x1);
	}
	
	if (clf[1] == 't') {
	    if (di != 'N') {
		/* Have triangular matrix with unit diagonal */
		if (ulf != ult) {
		    /* Returning identity matrix */

#define SPARSE_FS(_XASSIGN_)					\
		    do {					\
			for (j = 0; j < n; ++j) {		\
			    *(pi1++) = j;			\
			    _XASSIGN_; /* *(px1++) = _ONE_; */	\
			}					\
		    } while (0)

#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)	\
		    do {				\
			_CTYPE_ *px1 = _PTR_(x1);	\
			SPARSE_FS(*(px1++) = _ONE_);	\
		    } while (0)

		    if (clf[0] == 'n')
			SPARSE_FS();
		    else
			SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		    
		} else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
		    /* Returning symmetric matrix
		       with unit "trailing" diagonal */

#define SPARSE_FS(_XASSIGN_, _XASSIGN_ONE_)				\
		    do {						\
			for (j = 0, k = 0; j < n; ++j) {		\
			    kend = pp0[j];				\
			    while (k < kend) {				\
				*(pi1++) = pi0[k];			\
				_XASSIGN_; /* *(px1++) = px0[k]; */	\
				++k;					\
			    }						\
			    *(pi1++) = j;				\
			    _XASSIGN_ONE_; /* *(px1++) = _ONE_; */	\
			}						\
		    } while (0)
		    
#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		    do {						\
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);	\
			SPARSE_FS(*(px1++) = px0[k], *(px1++) = _ONE_);	\
		    } while (0)

		    if (clf[0] == 'n')
			SPARSE_FS(, );
		    else
			SPARSE_CASES(tx, SPARSE_FS_X);
		    
#undef SPARSE_FS
		    
		} else {
		    /* Returning symmetric matrix
		       with unit "leading" diagonal */

#define SPARSE_FS(_XASSIGN_, _XASSIGN_ONE_)				\
		    do {						\
			for (j = 0, k = 0; j < n; ++j) {		\
			    *(pi1++) = j;				\
			    _XASSIGN_ONE_; /* *(px1++) = _ONE_; */	\
			    kend = pp0[j];				\
			    while (k < kend) {				\
				*(pi1++) = pi0[k];			\
				_XASSIGN_; /* *(px1++) = px0[k]; */	\
				++k;					\
			    }						\
			}						\
		    } while (0)
		    		    
		    if (clf[0] == 'n')
			SPARSE_FS(, );
		    else
			SPARSE_CASES(tx, SPARSE_FS_X);
		    
#undef SPARSE_FS_X
#undef SPARSE_FS
		    
		}
	    } else if (ulf == ((clf[2] == 'C') ? 'U' : 'L')) {
		/* Have triangular matrix with non-unit "trailing" diagonal
		   and returning diagonal part */

#define SPARSE_FS(_XASSIGN_)						\
		do {							\
		    for (j = 0; j < n; ++j) {				\
			if (pp0[j-1] < pp0[j] && pi0[pp0[j]-1] == j) {	\
			    *(pi1++) = j;				\
			    _XASSIGN_; /* *(px1++) = px0[pp0[j]-1]; */	\
			}						\
		    }							\
		} while (0)
		
#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		do {							\
		    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
		    SPARSE_FS(*(px1++) = px0[pp0[j]-1]);		\
		} while (0)

		if (clf[0] == 'n')
		    SPARSE_FS();
		else
		    SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		
	    } else {
		/* Have triangular matrix with non-unit "leading" diagonal 
		   and returning diagonal part */

#define SPARSE_FS(_XASSIGN_)						\
		do {							\
		    for (j = 0; j < n; ++j) {				\
			if (pp0[j-1] < pp0[j] && pi0[pp0[j-1]] == j) {	\
			    *(pi1++) = j;				\
			    _XASSIGN_; /* *(px1++) = px0[pp0[j-1]]; */	\
			}						\
		    }							\
		} while (0)
		
#define SPARSE_FS_X(_CTYPE_, _PTR_, _ZERO_, _ONE_)			\
		do {							\
		    _CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);		\
		    SPARSE_FS(*(px1++) = px0[pp0[j-1]]);		\
		} while (0)

		if (clf[0] == 'n')
		    SPARSE_FS();
		else
		    SPARSE_CASES(tx, SPARSE_FS_X);

#undef SPARSE_FS_X
#undef SPARSE_FS
		
	    }
	} else if (ult == ((clf[2] == 'C') ? 'U' : 'L')) {
	    /* Have general matrix and returning upper triangle */

#define SPARSE_FS(_XASSIGN_)					\
	    do {						\
		for (j = 0, k = 0; j < n; ++j) {		\
		    kend = pp0[j];				\
		    while (k < kend) {				\
			if (pi0[k] <= j) {			\
			    *(pi1++) = pi0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
			++k;					\
		    }						\
		}						\
	    } while (0)
	    
	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS
	    
	} else {
	    /* Have general matrix and returning lower triangle */

#define SPARSE_FS(_XASSIGN_)					\
	    do {						\
		for (j = 0, k = 0; j < n; ++j) {		\
		    kend = pp0[j];				\
		    while (k < kend) {				\
			if (pi0[k] >= j) {			\
			    *(pi1++) = pi0[k];			\
			    _XASSIGN_; /* *(px1++) = px0[k]; */	\
			}					\
			++k;					\
		    }						\
		}						\
	    } while (0)
	    
	    if (clf[0] == 'n')
		SPARSE_FS();
	    else
		SPARSE_CASES(tx, SPARSE_FS_X_BASIC);

#undef SPARSE_FS_X_BASIC
#undef SPARSE_FS

	}
    }
    
    UNPROTECT((clf[0] == 'n') ? 3 : 4);
    return to;
}

/* as(<[CR]sparseMatrix>, "TsparseMatrix") */
SEXP CRsparse_as_Tsparse(SEXP from)
{
    static const char *valid[] = { VALID_CRSPARSE, "" };
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "CRsparse_as_Tsparse");
    const char *clf = valid[ivalid];
    char *clt = strdup(clf);
    clt[2] = 'T';

    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	dim = GET_SLOT(from, Matrix_DimSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (clf[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
    if (clf[1] == 't')
	SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    else
	SET_SLOT(to, Matrix_factorSym, GET_SLOT(from, Matrix_factorSym));
    if (clf[0] != 'n')
	SET_SLOT(to, Matrix_xSym, GET_SLOT(from, Matrix_xSym));
    
    SEXP iSym, jSym;
    int n, nnz, *pp = INTEGER(GET_SLOT(from, Matrix_pSym));
    if (clf[2] == 'C') {
	iSym = Matrix_iSym;
	jSym = Matrix_jSym;
	n = INTEGER(dim)[1];
    } else {
	iSym = Matrix_jSym;
	jSym = Matrix_iSym;
	n = INTEGER(dim)[0];
    }
    nnz = pp[n];
    
    SEXP j1 = PROTECT(allocVector(INTSXP, nnz));
    int j, k, kend, *pj1 = INTEGER(j1);
    SET_SLOT(to, iSym, GET_SLOT(from, iSym));
    SET_SLOT(to, jSym, j1);
    
    for (j = 0, k = 0; j < n; ++j) {
	kend = *(++pp);
	while (k < kend) {
	    *(pj1++) = j;
	    ++k;
	}
    }
    
    free(clt);
    UNPROTECT(2);
    return to;
}

/* as(t(<[CR]sparseMatrix>), "[RC]sparseMatrix") */
SEXP tCRsparse_as_RCsparse(SEXP from)
{
    static const char *valid[] = { VALID_CRSPARSE, "" };
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

/* isDiagonal(<[CR]sparseMatrix>) */
#define CR_IS_DIAGONAL(_C_, _I_)					\
SEXP _C_ ## sparse_is_diagonal(SEXP obj)				\
{									\
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];	\
    if (pdim[1] != n)							\
	return ScalarLogical(0);					\
    int *pp = INTEGER(GET_SLOT(obj, Matrix_pSym));			\
    if (pp[n] > n)							\
	return ScalarLogical(0);					\
    int d, j,								\
	*pi = INTEGER(GET_SLOT(obj, Matrix_ ## _I_ ## Sym));		\
    for (j = 0; j < n; ++j)						\
	if ((d = pp[j+1] - pp[j]) > 1 || (d == 1 && *(pi++) != j))	\
	    return ScalarLogical(0);					\
    return ScalarLogical(1);						\
}

/* Csparse_is_diagonal() */
CR_IS_DIAGONAL(C, i)
/* Rsparse_is_diagonal() */
CR_IS_DIAGONAL(R, j)

/* isDiagonal(<TsparseMatrix>) */
SEXP Tsparse_is_diagonal(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    int nnz = LENGTH(i);
    if (nnz > n)
	return ScalarLogical(0);
    int k,
	*pi = INTEGER(i),
	*pj = INTEGER(GET_SLOT(obj, Matrix_jSym));
    for (k = 0; k < nnz; ++k)
	if (*(pi++) != *(pj++))
	    return ScalarLogical(0);
    return ScalarLogical(1);
}

#define RETURN_TRUE_OF_KIND(_KIND_)					\
    do {								\
	SEXP ans = PROTECT(allocVector(LGLSXP, 1));			\
	LOGICAL(ans)[0] = 1;						\
	setAttrib(ans, install("kind"), _KIND_);			\
	UNPROTECT(1);							\
	return ans;							\
    } while (0)

/* isTriangular(<.g[CR]Matrix>, upper) */
#define CR_IS_TRIANGULAR(_C_, _I_, _UPPER_, _LOWER_)			\
SEXP _C_ ## sparse_is_triangular(SEXP obj, SEXP upper)			\
{									\
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];	\
    if (pdim[1] != n)							\
	return ScalarLogical(0);					\
    int j, k, kend,							\
	*pp = INTEGER(GET_SLOT(obj, Matrix_pSym)),			\
	*pi = INTEGER(GET_SLOT(obj, Matrix_ ## _I_ ## Sym)),		\
	need_upper = asLogical(upper);					\
    ++pp;								\
    if (need_upper == NA_LOGICAL) {					\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_LOWER_)						\
		    goto opposite;					\
		++k;							\
	    }								\
	}								\
	RETURN_TRUE_OF_KIND(mkString("U"));				\
    opposite:								\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_UPPER_)						\
		    return ScalarLogical(0);				\
		++k;							\
	    }								\
	}								\
	RETURN_TRUE_OF_KIND(mkString("L"));				\
    } else if (need_upper != 0) {					\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_LOWER_)						\
		    return ScalarLogical(0);				\
		++k;							\
	    }								\
	}								\
    } else {								\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if (_UPPER_)						\
		    return ScalarLogical(0);				\
		++k;							\
	    }								\
	}								\
    }									\
    return ScalarLogical(1);						\
}

/* Csparse_is_triangular() */
CR_IS_TRIANGULAR(C, i, pi[k] < j, pi[k] > j)
/* Rsparse_is_triangular() */
CR_IS_TRIANGULAR(R, j, pi[k] > j, pi[k] < j)

/* isTriangular(<.gTMatrix>, upper) */
SEXP Tsparse_is_triangular(SEXP obj, SEXP upper)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return ScalarLogical(0);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    int k, nnz = LENGTH(i),
	*pi = INTEGER(i),
	*pj = INTEGER(GET_SLOT(obj, Matrix_jSym)),
	need_upper = asLogical(upper);
    if (need_upper == NA_LOGICAL) {
	for (k = 0; k < nnz; ++k)
	    if (pi[k] > pj[k])
		goto opposite;
	RETURN_TRUE_OF_KIND(mkString("U"));
    opposite:
	for (k = 0; k < nnz; ++k)
	    if (pi[k] < pj[k])
		return ScalarLogical(0);
	RETURN_TRUE_OF_KIND(mkString("L"));
    } else if (need_upper != 0) {
	for (k = 0; k < nnz; ++k)
	    if (pi[k] > pj[k])
		return ScalarLogical(0);
    } else {
	for (k = 0; k < nnz; ++k)
	    if (pi[k] < pj[k])
		return ScalarLogical(0);
    }
    return ScalarLogical(1);
}

#define CR_IS_SYMMETRIC_LOOP(_XCOND_)					\
    do {								\
	for (j = 0, k = 0; j < n; ++j) {				\
	    kend = pp[j];						\
	    while (k < kend) {						\
		if ((i = pi[k]) >= j) {					\
		    if (i == j)						\
			++pp_[j];					\
		    k = kend;						\
		    break;						\
		}							\
		if (pp_[i] == pp[i] || pi[pp_[i]] != j || (_XCOND_))	\
		    return ScalarLogical(0);				\
		++pp_[i];						\
		++pp_[j];						\
		++k;							\
	    }								\
	}								\
    } while (0)

/* isSymmetric(<.g[CR]Matrix>, tol = 0, checkDN) */
#define CR_IS_SYMMETRIC(_C_, _I_)					\
SEXP _C_ ## sparse_is_symmetric(SEXP obj, SEXP checkDN)			\
{									\
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];	\
    if (pdim[1] != n)							\
	return ScalarLogical(0);					\
    if (asLogical(checkDN) != 0 &&					\
	!DimNames_is_symmetric(GET_SLOT(obj, Matrix_DimNamesSym)))	\
	return ScalarLogical(0);					\
    int i, j, k, kend, *pp_,						\
	*pp = INTEGER(GET_SLOT(obj, Matrix_pSym)),			\
	*pi = INTEGER(GET_SLOT(obj, Matrix_ ## _I_ ## Sym));		\
    Calloc_or_Alloca_TO(pp_, n, int);					\
    Memcpy(pp_, pp, n);							\
    ++pp;								\
									\
    /* For all X[i,j] in "leading" triangle,				\
       need that X[j,i] exists and X[j,i] == X[i,j] */			\
    if (R_has_slot(obj, Matrix_xSym)) {					\
	SEXP x = GET_SLOT(obj, Matrix_xSym);				\
	SEXPTYPE tx = TYPEOF(x);					\
	switch (tx) {							\
	case REALSXP:							\
	{								\
	    double *px = REAL(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		ISNAN(px[pp_[i]])					\
		? !ISNAN(px[k])						\
		: (ISNAN(px[k]) || px[pp_[i]] != px[k]));		\
	    break;							\
	}								\
	case LGLSXP:							\
	{								\
	    int *px = LOGICAL(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		px[pp_[i]] == NA_LOGICAL				\
		? (px[k] != NA_LOGICAL)					\
		: (px[k] == NA_LOGICAL || px[pp_[i]] != px[k]));	\
	    break;							\
	}								\
	case INTSXP:							\
	{								\
	    int *px = INTEGER(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		px[pp_[i]] == NA_INTEGER				\
		? (px[k] != NA_INTEGER)					\
		: (px[k] == NA_INTEGER || px[pp_[i]] != px[k]));	\
	    break;							\
	}								\
	case CPLXSXP:							\
	{								\
	    Rcomplex *px = COMPLEX(x);					\
	    CR_IS_SYMMETRIC_LOOP(					\
		ISNAN(px[pp_[i]].r) || ISNAN(px[pp_[i]].i)		\
		? !(ISNAN(px[k].r) || ISNAN(px[k].i))			\
		: (ISNAN(px[k].r) || ISNAN(px[k].i) ||			\
		   px[pp_[i]].r != px[k].r || px[pp_[i]].i != px[k].i)); \
	    break;							\
	}								\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", tx, "[CR]sparse_is_symmetric"); \
	    break;							\
	}								\
    } else {								\
	CR_IS_SYMMETRIC_LOOP(0);					\
    }									\
									\
    /* Need upper, lower triangles to have same number of nonzero elements */ \
    for (j = 0; j < n; ++j)						\
	if (pp_[j] != pp[j])						\
	    return ScalarLogical(0);					\
									\
    Free_FROM(pp_, n);							\
    return ScalarLogical(1);						\
}

/* Csparse_is_symmetric() */
CR_IS_SYMMETRIC(C, i)
/* Rsparse_is_symmetric() */
CR_IS_SYMMETRIC(R, j)