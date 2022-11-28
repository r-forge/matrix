/* unfinished and not-yet-used, along with ../R/subscript.R */
#include "subscript.h"

#define SAFEMULT(m, n)							\
    (((double) m * n <= R_XLEN_T_MAX) ? (R_xlen_t) m * n : R_XLEN_T_MAX)

#define F_X( _X_)   (_X_)
#define F_ND(_X_) ((_X_) ? 1 : 0)
#define F_NS(_X_)          1

static SEXP unpackedMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    R_xlen_t len = XLENGTH(w);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1);

    int ge = cl[1] == 'g', tr = cl[1] == 't', upper = 1, nonunit = 1;
    if (!ge) {
	SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
	upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1);
	if (tr) {
	    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	    nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	    UNPROTECT(1);
	}
    }

    R_xlen_t w1s, l, i_, j_;

#define SUB1_CASES(_SUB1_N_, _SUB1_X_, _F_N_, _F_X_)			\
    do {								\
	switch (cl[0]) {						\
	case 'n':							\
	    _SUB1_N_(int, LOGICAL, NA_LOGICAL, 0, 1, _F_N_);		\
	    break;							\
	case 'l':							\
	    _SUB1_X_(int, LOGICAL, NA_LOGICAL, 0, 1, _F_X_);		\
	    break;							\
	case 'i':							\
	    _SUB1_X_(int, INTEGER, NA_INTEGER, 0, 1, _F_X_);		\
	    break;							\
	case 'd':							\
	    _SUB1_X_(double, REAL, NA_REAL, 0.0, 1.0, _F_X_);		\
	    break;							\
	case 'z':							\
	    _SUB1_X_(Rcomplex, COMPLEX,					\
		     Matrix_zna, Matrix_zzero, Matrix_zone, _F_X_);	\
	    break;							\
	default:							\
	    break;							\
	}								\
    } while (0)
    
#define SUB1_N(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	_CTYPE_ *pres = _PTR_(res);					\
	R_xlen_t mn = SAFEMULT(m, n);					\
	if (TYPEOF(w) == INTSXP) {					\
	    int *pw = INTEGER(w);					\
	    SUB1_LOOP((pw[l] == NA_INTEGER || (R_xlen_t) pw[l] > mn),	\
		      _NA_, _ZERO_, _ONE_, _F_);			\
	} else {							\
	    double *pw = REAL(w), mn1a = (double) mn + 1.0;		\
	    SUB1_LOOP((ISNAN(pw[l]) || pw[l] >= mn1a),			\
		      _NA_, _ZERO_, _ONE_, _F_);			\
	}								\
    } while (0)

#define SUB1_X(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	PROTECT(x = GET_SLOT(x, Matrix_xSym));			\
	_CTYPE_ *px = _PTR_(x);					\
	SUB1_N(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_);	\
	UNPROTECT(1);						\
    } while (0)
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else {						\
		w1s = (R_xlen_t) pw[l] - 1;			\
		if (ge)						\
		    pres[l] = _F_(px[w1s]);			\
		else {						\
		    i_ = w1s % m;				\
		    j_ = w1s / m;				\
		    if (tr) {					\
			if ((upper) ? i_ > j_ : i_ < j_)	\
			    pres[l] = _ZERO_;			\
			else if (!nonunit && i_ == j_)		\
			    pres[l] = _ONE_;			\
			else					\
			    pres[l] = _F_(px[w1s]);		\
		    } else {					\
			if ((upper) ? i_ > j_ : i_ < j_)	\
			    pres[l] = _F_(px[i_ * m + j_]);	\
			else					\
			    pres[l] = _F_(px[w1s]);		\
		    }						\
		}						\
	    }							\
	}							\
    } while (0)
    
    SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);

#undef SUB1_LOOP
    
    UNPROTECT(1);
    return res;
}

static SEXP packedMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    R_xlen_t len = XLENGTH(w);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int m = INTEGER(dim)[0], n = m;
    UNPROTECT(1);

    int tr = cl[1] == 't', upper = 1, nonunit = 1;
    SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
    upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    UNPROTECT(1);
    if (tr) {
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	UNPROTECT(1);
    }
    
    R_xlen_t w1s, l, i_, j_, m2 = (R_xlen_t) m * 2;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	for (l = 0; l < len; ++l) {					\
	    if (_NA_SUBSCRIPT_)						\
		pres[l] = _NA_;						\
	    else {							\
		w1s = (R_xlen_t) pw[l] - 1;				\
		i_ = w1s % m;						\
		j_ = w1s / m;						\
		if (tr) {						\
		    if (upper) {					\
			if (i_ > j_)					\
			    pres[l] = _ZERO_;				\
			else if (!nonunit && i_ == j_)			\
			    pres[l] = _ONE_;				\
			else						\
			    pres[l] = _F_(px[PM_AR21_UP(i_, j_)]);	\
		    } else {						\
			if (i_ < j_)					\
			    pres[l] = _ZERO_;				\
			else if (!nonunit && i_ == j_)			\
			    pres[l] = _ONE_;				\
			else						\
			    pres[l] = _F_(px[PM_AR21_LO(i_, j_, m2)]);	\
		    }							\
		} else {						\
		    if (upper) {					\
			if (i_ > j_)					\
			    pres[l] = _F_(px[PM_AR21_UP(j_, i_)]);	\
			else						\
			    pres[l] = _F_(px[PM_AR21_UP(i_, j_)]);	\
		    } else {						\
			if (i_ < j_)					\
			    pres[l] = _F_(px[PM_AR21_LO(j_, i_, m2)]);	\
			else						\
			    pres[l] = _F_(px[PM_AR21_LO(i_, j_, m2)]);	\
		    }							\
		}							\
	    }								\
	}								\
    } while (0)

    SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(1);
    return res;
}

static SEXP CsparseMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    R_xlen_t len = XLENGTH(w);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1);
    
    SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
	i = PROTECT(GET_SLOT(x, Matrix_iSym));
    int *pp = INTEGER(p), *pi = INTEGER(i);

    R_xlen_t w1s, l = 0;
    int i_, j_, j, k = 0, kend;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	if (_NA_SUBSCRIPT_)						\
	    j = -1;							\
	else {								\
	    w1s = (R_xlen_t) pw[l] - 1;					\
	    i_ =     (int) (w1s % m);					\
	    j_ = j = (int) (w1s / m);					\
	}								\
	while (j >= 0) {						\
	    k = pp[j];							\
	    kend = pp[j+1];						\
	    while (k < kend && j_ == j) {				\
		if (pi[k] < i_)						\
		    ++k;						\
		else {							\
		    if (pi[k] > i_)					\
			pres[l] = _ZERO_;				\
		    else {						\
			pres[l] = _F_(px[k]);				\
			++k;						\
		    }							\
		    ++l;						\
		    if (l == len || _NA_SUBSCRIPT_)			\
			j_ = -1;					\
		    else {						\
			w1s = (R_xlen_t) pw[l] - 1;			\
			i_ = (int) (w1s % m);				\
			j_ = (int) (w1s / m);				\
		    }							\
		}							\
	    }								\
	    while (j_ == j) {						\
		pres[l] = _ZERO_;					\
		++l;							\
		if (l == len || _NA_SUBSCRIPT_)				\
		    j_ = -1;						\
		else {							\
		    w1s = (R_xlen_t) pw[l] - 1;				\
		    i_ = (int) (w1s % m);				\
		    j_ = (int) (w1s / m);				\
		}							\
	    }								\
	    j = j_;							\
	}								\
	while (l < len) {						\
	    pres[l] = _NA_;						\
	    ++l;							\
	}								\
    } while (0)
    
    SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(3);
    return res;
}

static SEXP RsparseMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    R_xlen_t len = XLENGTH(w);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1);
    
    SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
	j = PROTECT(GET_SLOT(x, Matrix_jSym));
    int *pp = INTEGER(p), *pj = INTEGER(j);

    R_xlen_t w1s, l = 0;
    int i_, j_, i, k = 0, kend;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	if (_NA_SUBSCRIPT_)						\
	    i = -1;							\
	else {								\
	    w1s = (R_xlen_t) pw[l] - 1;					\
	    i_ = i = (int) (w1s % m);					\
	    j_ =     (int) (w1s / m);					\
	}								\
	while (i >= 0) {						\
	    k = pp[i];							\
	    kend = pp[i+1];						\
	    while (k < kend && i_ == i) {				\
		if (pj[k] < j_)						\
		    ++k;						\
		else {							\
		    if (pj[k] > j_)					\
			pres[l] = _ZERO_;				\
		    else {						\
			pres[l] = _F_(px[k]);				\
			++k;						\
		    }							\
		    ++l;						\
		    if (l == len || _NA_SUBSCRIPT_)			\
			i_ = -1;					\
		    else {						\
			w1s = (R_xlen_t) pw[l] - 1;			\
			i_ = (int) (w1s % m);				\
			j_ = (int) (w1s / m);				\
		    }							\
		}							\
	    }								\
	    while (i_ == i) {						\
		pres[l] = _ZERO_;					\
		++l;							\
		if (l == len || _NA_SUBSCRIPT_)				\
		    i_ = -1;						\
		else {							\
		    w1s = (R_xlen_t) pw[l] - 1;				\
		    i_ = (int) (w1s % m);				\
		    j_ = (int) (w1s / m);				\
		}							\
	    }								\
	    i = i_;							\
	}								\
	while (l < len) {						\
	    pres[l] = _NA_;						\
	    ++l;							\
	}								\
    } while (0)

    SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);

#undef SUB1_LOOP
    
    UNPROTECT(3);
    return res;
}

static SEXP diagonalMatrix_subscript_1ary(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    R_xlen_t len = XLENGTH(w);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym)),
	diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
    int m = INTEGER(dim)[0], n = m, nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
    UNPROTECT(2);

    R_xlen_t w1s, l, n1a = (R_xlen_t) n + 1;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else if ((w1s = (R_xlen_t) pw[l] - 1) % n1a != 0)	\
		pres[l] = _ZERO_;				\
	    else if (!nonunit)					\
		pres[l] = _ONE_;				\
	    else						\
		pres[l] = _F_(px[w1s / n1a]);			\
	}							\
    } while (0)
    
    SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(1);
    return res;
}

static SEXP indMatrix_subscript_1ary(SEXP x, SEXP w)
{
    R_xlen_t len = XLENGTH(w);
    SEXP res = PROTECT(allocVector(LGLSXP, len));
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1);
    
    SEXP perm = PROTECT(GET_SLOT(x, Matrix_permSym));
    int *pperm = INTEGER(perm);
    
    R_xlen_t w1s, l;
    int i_, j_;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else {						\
		w1s = (R_xlen_t) pw[l] - 1;			\
		i_ = (int) (w1s % m);				\
		j_ = (int) (w1s / m);				\
		pres[l] = (j_ == pperm[i_] - 1) ? 1 : 0;	\
	    }							\
	}							\
    } while (0)

    SUB1_N(int, LOGICAL, NA_LOGICAL, 0, 1, );

#undef SUB1_LOOP
#undef SUB1_N
    
    UNPROTECT(2);
    return res;
}

/* x[i] with 'i' of type "integer" or "double" {in [1,2^52+1) or NA} */
SEXP R_subscript_1ary(SEXP x, SEXP i)
{
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(x, valid);
    const char *cl = valid[ivalid];
    if (ivalid < 0)
	ERROR_INVALID_CLASS(x, "R_subscript_1ary");
    
    switch (cl[2]) {
    case 'e':
    case 'y':
    case 'r':
	return unpackedMatrix_subscript_1ary(x, i, cl);
    case 'p':
	return   packedMatrix_subscript_1ary(x, i, cl);
	
    /* NB: for [CRT], the caller must preprocess 'x' and/or 'i';
           symmetric and unit triangular 'x' are not handled specially,
	   and it is assumed for speed that 'i' is sorted by row [R] 
	   or column [CT] with NA last
    */
	
    case 'C':
	return  CsparseMatrix_subscript_1ary(x, i, cl);
    case 'R':
	return  RsparseMatrix_subscript_1ary(x, i, cl);
    case 'T':
    {
        PROTECT(x = Tsparse_as_CRsparse(x, ScalarLogical(1)));
	char cl_[] = "..TMatrix";
	cl_[0] = cl[0];
	cl_[1] = cl[1];
	SEXP res = CsparseMatrix_subscript_1ary(x, i, cl_);
	UNPROTECT(1);
	return res;
    }
    case 'i':
	return diagonalMatrix_subscript_1ary(x, i, cl);
    default:
	return      indMatrix_subscript_1ary(x, i);
    }
}

static SEXP unpackedMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    int len = (int) (XLENGTH(w) / 2);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int m = INTEGER(dim)[0];
    UNPROTECT(1);
    
    int ge = cl[1] == 'g', tr = cl[1] == 't', upper = 1, nonunit = 1;
    if (!ge) {
	SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
	upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1);
	if (tr) {
	    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	    nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	    UNPROTECT(1);
	}
    }
    
    int l, i_, j_, *pw0 = INTEGER(w), *pw1 = pw0 + len;
    R_xlen_t m1 = (R_xlen_t) m;
    
#define SUB1_N(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	_CTYPE_ *pres = _PTR_(res);					\
	SUB1_LOOP((pw0[l] == NA_INTEGER || pw1[l] == NA_INTEGER),	\
		  _NA_, _ZERO_, _ONE_, _F_);				\
    } while (0)

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else {						\
		i_ = pw0[l] - 1;				\
		j_ = pw1[l] - 1;				\
		if (ge)						\
		    pres[l] = _F_(px[j_ * m1 + i_]);		\
		else if (tr) {					\
		    if ((upper) ? i_ > j_ : i_ < j_)		\
			pres[l] = _ZERO_;			\
		    else if (!nonunit && i_ == j_)		\
			pres[l] = _ONE_;			\
		    else					\
			pres[l] = _F_(px[j_ * m1 + i_]);	\
		} else {					\
		    if ((upper) ? i_ > j_ : i_ < j_)		\
			pres[l] = _F_(px[i_ * m1 + j_]);	\
		    else					\
			pres[l] = _F_(px[j_ * m1 + i_]);	\
		}						\
	    }							\
	}							\
    } while (0)
    
    SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(1);
    return res;
}

static SEXP packedMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    int len = (int) (XLENGTH(w) / 2);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int m = INTEGER(dim)[0];
    UNPROTECT(1);

    int tr = cl[1] == 't', upper = 1, nonunit = 1;
    SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
    upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    UNPROTECT(1);
    if (tr) {
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	UNPROTECT(1);
    }
    
    int l, i_, j_, *pw0 = INTEGER(w), *pw1 = pw0 + len;
    R_xlen_t m2 = (R_xlen_t) m * 2;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	for (l = 0; l < len; ++l) {					\
	    if (_NA_SUBSCRIPT_)						\
		pres[l] = _NA_;						\
	    else {							\
		i_ = pw0[l] - 1;					\
		j_ = pw1[l] - 1;					\
		if (tr) {						\
		    if (upper) {					\
			if (i_ > j_)					\
			    pres[l] = _ZERO_;				\
			else if (!nonunit && i_ == j_)			\
			    pres[l] = _ONE_;				\
			else						\
			    pres[l] = _F_(px[PM_AR21_UP(i_, j_)]);	\
		    } else {						\
			if (i_ < j_)					\
			    pres[l] = _ZERO_;				\
			else if (!nonunit && i_ == j_)			\
			    pres[l] = _ONE_;				\
			else						\
			    pres[l] = _F_(px[PM_AR21_LO(i_, j_, m2)]);	\
		    }							\
		} else {						\
		    if (upper) {					\
			if (i_ > j_)					\
			    pres[l] = _F_(px[PM_AR21_UP(j_, i_)]);	\
			else						\
			    pres[l] = _F_(px[PM_AR21_UP(i_, j_)]);	\
		    } else {						\
			if (i_ < j_)					\
			    pres[l] = _F_(px[PM_AR21_LO(j_, i_, m2)]);	\
			else						\
			    pres[l] = _F_(px[PM_AR21_LO(i_, j_, m2)]);	\
		    }							\
		}							\
	    }								\
	}								\
    } while (0)

    SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(1);
    return res;
}

static SEXP CsparseMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    int len = (int) (XLENGTH(w) / 2);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
	i = PROTECT(GET_SLOT(x, Matrix_iSym));
    int *pp = INTEGER(p), *pi = INTEGER(i);

    int l = 0, i_, j_, j, k = 0, kend, *pw0 = INTEGER(w), *pw1 = pw0 + len;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	if (_NA_SUBSCRIPT_)						\
	    j = -1;							\
	else {								\
	    i_ =     pw0[l] - 1;					\
	    j_ = j = pw1[l] - 1;					\
	}								\
	while (j >= 0) {						\
	    k = pp[j];							\
	    kend = pp[j+1];						\
	    while (k < kend && j_ == j) {				\
		if (pi[k] < i_)						\
		    ++k;						\
		else {							\
		    if (pi[k] > i_)					\
			pres[l] = _ZERO_;				\
		    else {						\
			pres[l] = _F_(px[k]);				\
			++k;						\
		    }							\
		    ++l;						\
		    if (l == len || _NA_SUBSCRIPT_)			\
			j_ = -1;					\
		    else {						\
			i_ = pw0[l] - 1;				\
			j_ = pw1[l] - 1;				\
		    }							\
		}							\
	    }								\
	    while (j_ == j) {						\
		pres[l] = _ZERO_;					\
		++l;							\
		if (l == len || _NA_SUBSCRIPT_)				\
		    j_ = -1;						\
		else {							\
		    i_ = pw0[l] - 1;					\
		    j_ = pw1[l] - 1;					\
		}							\
	    }								\
	    j = j_;							\
	}								\
	while (l < len) {						\
	    pres[l] = _NA_;						\
	    ++l;							\
	}								\
    } while (0)

    SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(3);
    return res;
}

static SEXP RsparseMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    int len = (int) (XLENGTH(w) / 2);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
	j = PROTECT(GET_SLOT(x, Matrix_jSym));
    int *pp = INTEGER(p), *pj = INTEGER(j);
    
    int l = 0, i_, j_, i, k = 0, kend, *pw0 = INTEGER(w), *pw1 = pw0 + len;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	if (_NA_SUBSCRIPT_)						\
	    i = -1;							\
	else {								\
	    i_ = i = pw0[l] - 1;					\
	    j_ =     pw1[l] - 1;					\
	}								\
	while (i >= 0) {						\
	    k = pp[i];							\
	    kend = pp[i+1];						\
	    while (k < kend && i_ == i) {				\
		if (pj[k] < j_)						\
		    ++k;						\
		else {							\
		    if (pj[k] > j_)					\
			pres[l] = _ZERO_;				\
		    else {						\
			pres[l] = _F_(px[k]);				\
			++k;						\
		    }							\
		    ++l;						\
		    if (l == len || _NA_SUBSCRIPT_)			\
			i_ = -1;					\
		    else {						\
			i_ = pw0[l] - 1;				\
			j_ = pw1[l] - 1;				\
		    }							\
		}							\
	    }								\
	    while (i_ == i) {						\
		pres[l] = _ZERO_;					\
		++l;							\
		if (l == len || _NA_SUBSCRIPT_)				\
		    i_ = -1;						\
		else {							\
		    i_ = pw0[l] - 1;					\
		    j_ = pw1[l] - 1;					\
		}							\
	    }								\
	    i = i_;							\
	}								\
	while (l < len) {						\
	    pres[l] = _NA_;						\
	    ++l;							\
	}								\
    } while (0)

    SUB1_CASES(SUB1_N, SUB1_X, F_NS, F_X);

#undef SUB1_LOOP
    
    UNPROTECT(3);
    return res;
}

static SEXP diagonalMatrix_subscript_1ary_mat(SEXP x, SEXP w, const char *cl)
{
    SEXPTYPE typ = kind2type(cl[0]);
    int len = (int) (XLENGTH(w) / 2);
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
    int nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
    UNPROTECT(1);

    int l, *pw0 = INTEGER(w), *pw1 = pw0 + len;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else if (pw0[l] != pw1[l])				\
		pres[l] = _ZERO_;				\
	    else if (!nonunit)					\
		pres[l] = _ONE_;				\
	    else						\
		pres[l] = _F_(px[pw0[l] - 1]);			\
	}							\
    } while (0)
    
    SUB1_CASES(SUB1_X, SUB1_X, F_ND, F_X);
    
#undef SUB1_LOOP
    
    UNPROTECT(1);
    return res;
}

static SEXP indMatrix_subscript_1ary_mat(SEXP x, SEXP w)
{
    int len = (int) (XLENGTH(w) / 2);
    SEXP res = PROTECT(allocVector(LGLSXP, len));
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP perm = PROTECT(GET_SLOT(x, Matrix_permSym));
    int *pperm = INTEGER(perm);
    
    int l, *pw0 = INTEGER(w), *pw1 = pw0 + len;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	for (l = 0; l < len; ++l) {					\
	    if (_NA_SUBSCRIPT_)						\
		pres[l] = _NA_;						\
	    else							\
		pres[l] = (pw1[l] == pperm[pw0[l] - 1]) ? 1 : 0;	\
	}								\
    } while (0)

    SUB1_N(int, LOGICAL, NA_LOGICAL, 0, 1, );

#undef SUB1_CASES
#undef SUB1_N
#undef SUB1_X
#undef SUB1_LOOP
    
    UNPROTECT(2);
    return res;
}

SEXP R_subscript_1ary_mat(SEXP x, SEXP i)
{
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(x, valid);
    const char *cl = valid[ivalid];
    if (ivalid < 0)
	ERROR_INVALID_CLASS(x, "R_subscript_1ary_mat");

    switch (cl[2]) {
    case 'e':
    case 'y':
    case 'r':
	return unpackedMatrix_subscript_1ary_mat(x, i, cl);
    case 'p':
	return   packedMatrix_subscript_1ary_mat(x, i, cl);
	
    /* NB: for [CRT], the caller must preprocess 'x' and/or 'i';
           symmetric and unit triangular 'x' are not handled specially,
	   and it is assumed for speed that 'i' is sorted by row [R] 
	   or column [CT] with NA last
    */
	
    case 'C':
	return  CsparseMatrix_subscript_1ary_mat(x, i, cl);
    case 'R':
	return  RsparseMatrix_subscript_1ary_mat(x, i, cl);
    case 'T':
    {
        PROTECT(x = Tsparse_as_CRsparse(x, ScalarLogical(1)));
	char cl_[] = "..TMatrix";
	cl_[0] = cl[0];
	cl_[1] = cl[1];
	SEXP res = CsparseMatrix_subscript_1ary_mat(x, i, cl_);
	UNPROTECT(1);
	return res;
    }
    case 'i':
	return diagonalMatrix_subscript_1ary_mat(x, i, cl);
    default:
	return      indMatrix_subscript_1ary_mat(x, i);
    }
}

static int keep_tr(int *pi, int *pj, int n, int upper, int nonunit)
{
    int ki;
    for (ki = 0; ki < n; ++ki) {
	if (pi[ki] == NA_INTEGER)
	    return 0;
    }
    int kj, j_;
    for (kj = 0; kj < n; ++kj) {
	if ((j_ = pj[kj]) == NA_INTEGER)
	    return 0;
	for (ki = 0; ki < kj; ++ki) {
	    if ((upper) ? pi[ki] >= j_ : pi[ki] <= j_)
		return 0;
	}
    }
    return (nonunit || memcmp(pi, pj, n * sizeof(int)) != 0) ? 1 : 2;
}

static int keep_sy(int *pi, int *pj, int n)
{
    return memcmp(pi, pj, n * sizeof(int)) == 0;
}

static int keep_di(int *pi, int *pj, int n, int nonunit)
{
    int ki;
    for (ki = 0; ki < n; ++ki) {
	if (pi[ki] == NA_INTEGER)
	    return 0;
    }
    int kj, j_;
    for (kj = 0; kj < n; ++kj) {
	if ((j_ = pj[kj]) == NA_INTEGER)
	    return 0;
	for (ki = 0; ki < kj; ++ki) {
	    if (pi[ki] == j_)
		return 0;
	}
	for (ki = kj+1; ki < n; ++ki) {
	    if (pi[ki] == j_)
		return 0;
	}
    }
    return (nonunit || memcmp(pi, pj, n * sizeof(int)) != 0) ? 1 : 2;
}

static SEXP unpackedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
					  const char *cl)
{
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1); /* dim */

    int mi = isNull(i), mj = isNull(j),
	*pi = (mi) ? NULL : INTEGER(i),
	*pj = (mj) ? NULL : INTEGER(j);
    R_xlen_t
	ni = (mi) ? (R_xlen_t) m : XLENGTH(i),
	nj = (mj) ? (R_xlen_t) n : XLENGTH(j);

    int upper = 1, nonunit = 1, keep = 0;
    if (cl[1] != 'g') {
	SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
	upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1); /* uplo */
	if (cl[1] == 't') {
	    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	    nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	    UNPROTECT(1); /* diag */
	}
    }
    
    char cl_[] = ".geMatrix";
    cl_[0] = cl[0];
    if (cl[1] != 'g' && !(mi || mj)) {
	if (cl[1] == 't') {
	    keep = keep_tr(pi, pj, m, upper, nonunit);
	    if (keep) {
		cl_[1] = 't';
		cl_[2] = 'r';
	    }
	} else {
	    keep = keep_sy(pi, pj, m);
	    if (keep) {
		cl_[1] = 's';
		cl_[2] = 'y';
	    }
	}
    }
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl_));

    PROTECT(dim = GET_SLOT(res, Matrix_DimSym));
    pdim = INTEGER(dim);
    pdim[0] = (int) ni;
    pdim[1] = (int) nj;
    UNPROTECT(1); /* dim */
    
    if (keep) {
	if (!upper) {
	    SEXP uplo = PROTECT(GET_SLOT(res, Matrix_uploSym)),
		uplo0 = PROTECT(mkChar("L"));
	    SET_STRING_ELT(uplo, 0, uplo0);
	    UNPROTECT(2); /* uplo, uplo0 */
	}
	if (!nonunit && keep > 1) {
	    SEXP diag = PROTECT(GET_SLOT(res, Matrix_diagSym)),
		diag0 = PROTECT(mkChar("U"));
	    SET_STRING_ELT(diag, 0, diag0);
	    UNPROTECT(2); /* diag, diag0 */
	}
    }

    if ((double) ni * nj > R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));

    SEXPTYPE tx1;
    R_xlen_t nx1 = (R_xlen_t) ni * nj, m_ = (R_xlen_t) m;
    SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
	x1 = PROTECT(allocVector(tx1 = TYPEOF(x0), nx1));

#define SUB2_CASES							\
    do {								\
	switch (tx1) {							\
	case LGLSXP:							\
	    SUB2(int, LOGICAL, NA_LOGICAL, 0, 1);			\
	    break;							\
	case INTSXP:							\
	    SUB2(int, INTEGER, NA_INTEGER, 0, 1);			\
	    break;							\
	case REALSXP:							\
	    SUB2(double, REAL, NA_REAL, 0.0, 1.0);			\
	    break;							\
	case CPLXSXP:							\
	    SUB2(Rcomplex, COMPLEX, Matrix_zna, Matrix_zzero, Matrix_zone); \
	    break;							\
	default:							\
	    break;							\
	}								\
    } while (0)

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	if (cl_[1] == 'g') {						\
	    if (cl[1] == 'g')						\
		SUB2_LOOP(for (ki = 0; ki < ni; ++ki),			\
			  XIJ_GE, , , _NA_, _ZERO_, _ONE_);		\
	    else if (cl[1] == 't') {					\
		if (upper) {						\
		    if (nonunit)					\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TR_U_N, , , _NA_, _ZERO_, _ONE_);	\
		    else						\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TR_U_U, , , _NA_, _ZERO_, _ONE_);	\
		} else {						\
		    if (nonunit)					\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TR_L_N, , , _NA_, _ZERO_, _ONE_);	\
		    else						\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TR_L_U, , , _NA_, _ZERO_, _ONE_);	\
		}							\
	    } else {							\
		if (upper)						\
		    SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			      XIJ_SY_U, , , _NA_, _ZERO_, _ONE_);	\
		else							\
		    SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			      XIJ_SY_L, , , _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else if (cl_[1] == 't') {					\
	    Matrix_memset(px1, 0, nx1, sizeof(_CTYPE_));		\
	    if (upper) {						\
		if (nonunit)						\
		    SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),		\
			      XIJ_TR_U_N, , px1 += ni-kj-1, _NA_, _ZERO_, _ONE_); \
		else							\
		    SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),		\
			      XIJ_TR_U_U, , px1 += ni-kj-1, _NA_, _ZERO_, _ONE_); \
	    } else {							\
		if (nonunit)						\
		    SUB2_LOOP(for (ki = kj; ki < ni; ++ki),		\
			      XIJ_TR_U_N, px1 += kj, , _NA_, _ZERO_, _ONE_); \
		else							\
		    SUB2_LOOP(for (ki = kj; ki < ni; ++ki),		\
			      XIJ_TR_U_U, px1 += kj, , _NA_, _ZERO_, _ONE_); \
	    }								\
	} else {							\
	    Matrix_memset(px1, 0, nx1, sizeof(_CTYPE_));		\
	    if (upper)							\
		SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),			\
			  XIJ_SY_U, , px1 += ni-kj-1, _NA_, _ZERO_, _ONE_); \
	    else							\
		SUB2_LOOP(for (ki = kj; ki < ni; ++ki),			\
			  XIJ_SY_L, px1 += kj, , _NA_, _ZERO_, _ONE_);	\
	}								\
    } while (0)

#define SUB2_LOOP(_FOR_, _XIJ_, _JUMP1_, _JUMP2_,			\
		  _NA_, _ZERO_, _ONE_)					\
    do {								\
	int i_, j_, ki, kj;						\
	for (kj = 0; kj < nj; ++kj) {					\
	    if (mj)							\
		j_ = kj;						\
	    else {							\
		j_ = pj[kj];						\
		if (j_ != NA_INTEGER)					\
		    j_ -= 1;						\
		else {							\
		    _JUMP1_;						\
		    _FOR_ {						\
			*(px1++) = _NA_;				\
		    }							\
		    _JUMP2_;						\
		    continue;						\
		}							\
	    }								\
	    _JUMP1_;							\
	    _FOR_ {							\
		if (mi)							\
		    i_ = ki;						\
		else {							\
		    i_ = pi[ki];					\
		    if (i_ != NA_INTEGER)				\
			i_ -= 1;					\
		    else {						\
			*(px1++) = _NA_;				\
			continue;					\
		    }							\
		}							\
		*(px1++) = _XIJ_(px0, i_, j_, m_, _ZERO_, _ONE_);	\
	    }								\
	    _JUMP2_;							\
	}								\
    } while (0)

#define XIJ_GE(    _X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    *(_X_ + _J_ * _M1_ + _I_)

#define XIJ_TR_U_N(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    ((_I_ <= _J_)						\
     ? XIJ_GE(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
     : _ZERO_)

#define XIJ_TR_U_U(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    ((_I_ < _J_)						\
     ? XIJ_GE(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
     : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TR_L_N(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    ((_I_ >= _J_)						\
     ? XIJ_GE(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
     : _ZERO_)

#define XIJ_TR_L_U(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    ((_I_ > _J_)						\
     ? XIJ_GE(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
     : ((_I_ == _J_) ? _ONE_ : _ZERO_))
    
#define XIJ_SY_U(  _X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    ((_I_ <= _J_)						\
     ? XIJ_GE(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
     : XIJ_GE(_X_, _J_, _I_, _M1_, _ZERO_, _ONE_))

#define XIJ_SY_L(  _X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
    (_I_ >= _J_							\
     ? XIJ_GE(_X_, _I_, _J_, _M1_, _ZERO_, _ONE_)		\
     : XIJ_GE(_X_, _J_, _I_, _M1_, _ZERO_, _ONE_))
    
    SUB2_CASES;

#undef SUB2
#undef XIJ_GE
#undef XIJ_TR_U_N
#undef XIJ_TR_U_U
#undef XIJ_TR_L_N
#undef XIJ_TR_L_U
#undef XIJ_SY_U
#undef XIJ_SY_L
    
    SET_SLOT(res, Matrix_xSym, x1);
    
    UNPROTECT(3); /* x1, x0, res */
    return res;
}

static SEXP packedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
					const char *cl)
{
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int *pdim = INTEGER(dim), m = pdim[0];
    UNPROTECT(1); /* dim */

    if ((double) m * m > R_XLEN_T_MAX)
	error(_("x[i,j] is not supported for n-by-n packedMatrix 'x' "
		"with n*n exceeding R_XLEN_T_MAX"));
    
    int mi = isNull(i), mj = isNull(j),
	*pi = (mi) ? NULL : INTEGER(i),
	*pj = (mj) ? NULL : INTEGER(j);
    R_xlen_t
	ni = (mi) ? (R_xlen_t) m : XLENGTH(i),
	nj = (mj) ? (R_xlen_t) m : XLENGTH(j);

    int upper = 1, nonunit = 1, keep = 0;
    SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
    upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    UNPROTECT(1); /* uplo */
    if (cl[1] == 't') {
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	nonunit = *CHAR(STRING_ELT(diag, 0)) == 'N';
	UNPROTECT(1); /* diag */
    }
    
    char cl_[] = ".geMatrix";
    cl_[0] = cl[0];
    if (!(mi || mj)) {
	if (cl[1] == 't') {
	    keep = keep_tr(pi, pj, m, upper, nonunit);
	    if (keep) {
		cl_[1] = 't';
		cl_[2] = 'p';
	    }
	} else {
	    keep = keep_sy(pi, pj, m);
	    if (keep) {
		cl_[1] = 's';
		cl_[2] = 'p';
	    }
	}
    }
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl_));

    PROTECT(dim = GET_SLOT(res, Matrix_DimSym));
    pdim = INTEGER(dim);
    pdim[0] = (int) ni;
    pdim[1] = (int) nj;
    UNPROTECT(1); /* dim */
    
    if (keep) {
	if (!upper) {
	    SEXP uplo = PROTECT(GET_SLOT(res, Matrix_uploSym)),
		uplo0 = PROTECT(mkChar("L"));
	    SET_STRING_ELT(uplo, 0, uplo0);
	    UNPROTECT(2); /* uplo, uplo0 */
	}
	if (!nonunit && keep > 1) {
	    SEXP diag = PROTECT(GET_SLOT(res, Matrix_diagSym)),
		diag0 = PROTECT(mkChar("U"));
	    SET_STRING_ELT(diag, 0, diag0);
	    UNPROTECT(2); /* diag, diag0 */
	}
    }

    if (((keep) ? ((double) ni * ((double) ni + 1.0)) / 2.0 : (double) ni * nj)
	> R_XLEN_T_MAX)
	error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
    
    SEXPTYPE tx1;
    R_xlen_t nx1 = (keep) ? ni + (ni * (ni - 1)) / 2 : ni * nj,
	m_ = (R_xlen_t) m * 2;
    SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
	x1 = PROTECT(allocVector(tx1 = TYPEOF(x0), nx1));

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_)			\
    do {								\
	_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1);			\
	if (cl_[1] == 'g') {						\
	    if (cl[1] == 't') {						\
		if (upper) {						\
		    if (nonunit)					\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_);	\
		    else						\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_);	\
		} else {						\
		    if (nonunit)					\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_);	\
		    else						\
			SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
				  XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_);	\
		}							\
	    } else {							\
		if (upper)						\
		    SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			      XIJ_SP_U, , , _NA_, _ZERO_, _ONE_);	\
		else							\
		    SUB2_LOOP(for (ki = 0; ki < ni; ++ki),		\
			      XIJ_SP_L, , , _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else if (cl_[1] == 't') {					\
	    if (upper) {						\
		if (nonunit)						\
		    SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),		\
			      XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_);	\
		else							\
		    SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),		\
			      XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_);	\
	    } else {							\
		if (nonunit)						\
		    SUB2_LOOP(for (ki = kj; ki < ni; ++ki),		\
			      XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_);	\
		else							\
		    SUB2_LOOP(for (ki = kj; ki < ni; ++ki),		\
			      XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_);	\
	    }								\
	} else {							\
	    if (upper)							\
		SUB2_LOOP(for (ki = 0; ki <= kj; ++ki),			\
			  XIJ_SP_U, , , _NA_, _ZERO_, _ONE_);		\
	    else							\
		SUB2_LOOP(for (ki = kj; ki < ni; ++ki),			\
			  XIJ_SP_L, , , _NA_, _ZERO_, _ONE_);		\
	}								\
    } while (0)
    
#define XIJ_TP_U_N(_X_, _I_, _J_, _M2_, _ZERO_, _ONE_)		\
    ((_I_ <= _J_)						\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : _ZERO_)
    
#define XIJ_TP_U_U(_X_, _I_, _J_, _M2_, _ZERO_, _ONE_)		\
    ((_I_ < _J_)						\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TP_L_N(_X_, _I_, _J_, _M2_, _ZERO_, _ONE_)		\
    ((_I_ >= _J_)						\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _M2_))			\
     : _ZERO_)
    
#define XIJ_TP_L_U(_X_, _I_, _J_, _M2_, _ZERO_, _ONE_)		\
    ((_I_ > _J_)						\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _M2_))			\
     : ((_I_ == _J_) ? _ONE_ : _ZERO_))
    
#define XIJ_SP_U(  _X_, _I_, _J_, _M2_, _ZERO_, _ONE_)		\
    ((_I_ <= _J_)						\
     ? *(_X_ + PM_AR21_UP(_I_, _J_))				\
     : *(_X_ + PM_AR21_UP(_J_, _I_)))

#define XIJ_SP_L(  _X_, _I_, _J_, _M2_, _ZERO_, _ONE_)		\
    (_I_ >= _J_							\
     ? *(_X_ + PM_AR21_LO(_I_, _J_, _M2_))			\
     : *(_X_ + PM_AR21_LO(_J_, _I_, _M2_)))

    SUB2_CASES;

#undef SUB2_CASES
#undef SUB2
#undef SUB2_LOOP
#undef XIJ_TP_U_N
#undef XIJ_TP_U_U
#undef XIJ_TP_L_N
#undef XIJ_TP_L_U
#undef XIJ_SP_U
#undef XIJ_SP_L

    SET_SLOT(res, Matrix_xSym, x1);
    
    UNPROTECT(3); /* x1, x0, res */
    return res;
}

static SEXP  CsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
					  const char *cl)
{
    error(_("not yet implemented"));
    return R_NilValue;
}

static SEXP  RsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
					  const char *cl)
{
    error(_("not yet implemented"));
    return R_NilValue;    
}

static SEXP diagonalMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
					  const char *cl)
{
    error(_("not yet implemented"));
    return R_NilValue;
}

static SEXP      indMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j,
					  const char *cl)
{
    error(_("not yet implemented"));
    return R_NilValue;
}

SEXP R_subscript_2ary(SEXP x, SEXP i, SEXP j)
{
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(x, valid);
    const char *cl = valid[ivalid];
    if (ivalid < 0)
	ERROR_INVALID_CLASS(x, "R_subscript_2ary");

    switch (cl[2]) {
    case 'e':
    case 'y':
    case 'r':
	return unpackedMatrix_subscript_2ary(x, i, j, cl);
    case 'p':
	return   packedMatrix_subscript_2ary(x, i, j, cl);
    case 'C':
	return  CsparseMatrix_subscript_2ary(x, i, j, cl);
    case 'R':
	return  RsparseMatrix_subscript_2ary(x, i, j, cl);
    case 'T':
    {
        PROTECT(x = Tsparse_as_CRsparse(x, ScalarLogical(1)));
	char cl_[] = "..TMatrix";
	cl_[0] = cl[0];
	cl_[1] = cl[1];
	SEXP res = CsparseMatrix_subscript_2ary(x, i, j, cl_);
	UNPROTECT(1);
	return res;
    }
    case 'i':
	return diagonalMatrix_subscript_2ary(x, i, j, cl);
    default:
	return      indMatrix_subscript_2ary(x, i, j, cl);
    }
}
