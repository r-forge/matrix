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

    int ge = cl[1] == 'g', tr = cl[1] == 't', upper = 0, unit = 0;
    if (!ge) {
	SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
	upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1);
	if (tr) {
	    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	    unit = *CHAR(STRING_ELT(diag, 0)) != 'N';
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
			else if (unit && i_ == j_)		\
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

    int tr = cl[1] == 't', upper = 0, unit = 0;
    SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
    upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    UNPROTECT(1);
    if (tr) {
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	unit = *CHAR(STRING_ELT(diag, 0)) != 'N';
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
			else if (unit && i_ == j_)			\
			    pres[l] = _ONE_;				\
			else						\
			    pres[l] = _F_(px[PM_AR21_UP(i_, j_)]);	\
		    } else {						\
			if (i_ < j_)					\
			    pres[l] = _ZERO_;				\
			else if (unit && i_ == j_)			\
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
    int m = INTEGER(dim)[0], n = m, unit = *CHAR(STRING_ELT(diag, 0)) != 'N';
    UNPROTECT(2);

    R_xlen_t w1s, l, n1a = (R_xlen_t) n + 1;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else if ((w1s = (R_xlen_t) pw[l] - 1) % n1a != 0)	\
		pres[l] = _ZERO_;				\
	    else if (unit)					\
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
    R_xlen_t len = XLENGTH(w) / 2;
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int m = INTEGER(dim)[0];
    UNPROTECT(1);
    
    int ge = cl[1] == 'g', tr = cl[1] == 't', upper = 0, unit = 0;
    if (!ge) {
	SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
	upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
	UNPROTECT(1);
	if (tr) {
	    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	    unit = *CHAR(STRING_ELT(diag, 0)) != 'N';
	    UNPROTECT(1);
	}
    }
    
    R_xlen_t l, i_, j_;
    int *pw0 = INTEGER(w), *pw1 = pw0 + len;

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
		i_ = (R_xlen_t) pw0[l] - 1;			\
		j_ = (R_xlen_t) pw1[l] - 1;			\
		if (ge)						\
		    pres[l] = _F_(px[j_ * m + i_]);		\
		else if (tr) {					\
		    if ((upper) ? i_ > j_ : i_ < j_)		\
			pres[l] = _ZERO_;			\
		    else if (unit && i_ == j_)			\
			pres[l] = _ONE_;			\
		    else					\
			pres[l] = _F_(px[j_ * m + i_]);		\
		} else {					\
		    if ((upper) ? i_ > j_ : i_ < j_)		\
			pres[l] = _F_(px[i_ * m + j_]);		\
		    else					\
			pres[l] = _F_(px[j_ * m + i_]);		\
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
    R_xlen_t len = XLENGTH(w) / 2;
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
    int m = INTEGER(dim)[0];
    UNPROTECT(1);

    int tr = cl[1] == 't', upper = 0, unit = 0;
    SEXP uplo = PROTECT(GET_SLOT(x, Matrix_uploSym));
    upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
    UNPROTECT(1);
    if (tr) {
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	unit = *CHAR(STRING_ELT(diag, 0)) != 'N';
	UNPROTECT(1);
    }
    
    R_xlen_t l, i_, j_, m2 = (R_xlen_t) m * 2;
    int *pw0 = INTEGER(w), *pw1 = pw0 + len;

#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)		\
    do {								\
	for (l = 0; l < len; ++l) {					\
	    if (_NA_SUBSCRIPT_)						\
		pres[l] = _NA_;						\
	    else {							\
		i_ = (R_xlen_t) pw0[l] - 1;				\
		j_ = (R_xlen_t) pw1[l] - 1;				\
		if (tr) {						\
		    if (upper) {					\
			if (i_ > j_)					\
			    pres[l] = _ZERO_;				\
			else if (unit && i_ == j_)			\
			    pres[l] = _ONE_;				\
			else						\
			    pres[l] = _F_(px[PM_AR21_UP(i_, j_)]);	\
		    } else {						\
			if (i_ < j_)					\
			    pres[l] = _ZERO_;				\
			else if (unit && i_ == j_)			\
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
    R_xlen_t len = XLENGTH(w) / 2;
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
	i = PROTECT(GET_SLOT(x, Matrix_iSym));
    int *pp = INTEGER(p), *pi = INTEGER(i);

    R_xlen_t l = 0;
    int i_, j_, j, k = 0, kend, *pw0 = INTEGER(w), *pw1 = pw0 + len;

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
    R_xlen_t len = XLENGTH(w) / 2;
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym)),
	j = PROTECT(GET_SLOT(x, Matrix_jSym));
    int *pp = INTEGER(p), *pj = INTEGER(j);
    
    R_xlen_t l = 0;
    int i_, j_, i, k = 0, kend, *pw0 = INTEGER(w), *pw1 = pw0 + len;

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
			i_ = 0;						\
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
    R_xlen_t len = XLENGTH(w) / 2;
    SEXP res = allocVector(typ, len);
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
    int unit = *CHAR(STRING_ELT(diag, 0)) != 'N';
    UNPROTECT(1);

    R_xlen_t l;
    int *pw0 = INTEGER(w), *pw1 = pw0 + len;
    
#define SUB1_LOOP(_NA_SUBSCRIPT_, _NA_, _ZERO_, _ONE_, _F_)	\
    do {							\
	for (l = 0; l < len; ++l) {				\
	    if (_NA_SUBSCRIPT_)					\
		pres[l] = _NA_;					\
	    else if (pw0[l] != pw1[l])				\
		pres[l] = _ZERO_;				\
	    else if (unit)					\
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
    R_xlen_t len = XLENGTH(w) / 2;
    SEXP res = PROTECT(allocVector(LGLSXP, len));
    if (len == 0)
	return res;
    PROTECT(res);
    
    SEXP perm = PROTECT(GET_SLOT(x, Matrix_permSym));
    int *pperm = INTEGER(perm);
    
    R_xlen_t l;
    int *pw0 = INTEGER(w), *pw1 = pw0 + len;
    
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

#undef SUB1_LOOP
#undef SUB1_N
#undef SUB1_X
#undef SUB1_CASES
    
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

#if 0

static SEXP unpackedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{
    
}

static SEXP   packedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP  CsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP  RsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP  TsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP diagonalMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{
    
}

static SEXP      indMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

SEXP R_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop)
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
	return unpackedMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'p':
	return   packedMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'C':
	return  CsparseMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'R':
	return  RsparseMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'T':
	return  TsparseMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'i':
	return diagonalMatrix_subscript_2ary(x, i, j, drop, cl);
    default:
	return      indMatrix_subscript_2ary(x, i, j, drop, cl);
    }
}

#endif
