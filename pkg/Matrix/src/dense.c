#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
#include "dense.h"
/* #include "chm_common.h" */

/* as(<denseMatrix>, "[CRT]sparseMatrix") */
SEXP R_dense_as_sparse(SEXP from, SEXP code, SEXP uplo, SEXP diag)
{
    const char *c;
    char c0, c1, c2;
    if ((code = asChar(code)) == NA_STRING ||
	(c0 = (c = CHAR(code))[0]) == '\0' ||
	((c1 = c[1]) != '.' && c1 != 'g' && c1 != 't' && c1 != 's') ||
	((c2 = c[2]) != 'C' && c2 != 'R' && c2 != 'T'))
	error(_("invalid 'code' to 'R_dense_as_sparse()'"));
    
    SEXPTYPE txf,
	txt = (c0 == '.') ? NILSXP : kind2type(c0); /* before doing more */

    char *cl, ul = 'U', di = 'N';
    static const char *valid[] = {
	"dgeMatrix", "dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix",
	"lgeMatrix", "ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix",
	"ngeMatrix", "ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid), nprotect = 0,
	*pdim = NULL, doDN = 1, packed = 0;
    SEXP dim, dimnames, x_from;
    if (ivalid >= 0) {
	cl = strdup(valid[ivalid]);
	packed = (cl[2] == 'p');
	dim = GET_SLOT(from, Matrix_DimSym);
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
	x_from = GET_SLOT(from, Matrix_xSym);
	txf = TYPEOF(x_from);
	if (cl[1] != 'g') {
	    uplo = GET_SLOT(from, Matrix_uploSym);
	    ul = *CHAR(STRING_ELT(uplo, 0));
	    if (cl[1] == 't') {
		diag = GET_SLOT(from, Matrix_diagSym);
		di = *CHAR(STRING_ELT(diag, 0));
	    }
	}
	pdim = INTEGER(dim);
    } else {
	/* 'from' is a base matrix or base vector, but we behave as though
	   it is the 'x' slot of a .(ge|tr|sy)Matrix (depending on 'code')
	   for efficiency, relying on the user to specify 'uplo' and 'diag'
	   as necessary ...
	*/
	x_from = from;
	txf = TYPEOF(x_from);
	char tmp[] = ".g.Matrix";
	cl = &tmp[0];
	cl[0] = type2kind(txf); /* before doing more */
	if (isMatrix(from)) {
	    dim = getAttrib(from, R_DimSymbol);
	    dimnames = getAttrib(from, R_DimNamesSymbol);
	    pdim = INTEGER(dim);
	    doDN = !isNull(dimnames);
	} else {
	    R_xlen_t len = XLENGTH(from);
	    if (len > INT_MAX)
		error(_("vector of length exceeding 2^31-1 "
			"to 'R_dense_as_sparse()'"));
	    PROTECT(dim = allocVector(INTSXP, 2));
	    ++nprotect;
	    pdim = INTEGER(dim);
	    pdim[0] = (int) len;
	    pdim[1] = 1;
	    SEXP nms = getAttrib(from, R_NamesSymbol);
	    doDN = !isNull(nms);
	    if (doDN) {
		PROTECT(dimnames = allocVector(VECSXP, 2));
		++nprotect;
		SET_VECTOR_ELT(dimnames, 0, nms);
	    }
	}
	if (c1 == 't' || c1 == 's') {
	    if (pdim[0] != pdim[1])
		error(_("attempt to construct triangular or symmetric "
			"%csparseMatrix from non-square matrix"), c2);
	    if ((uplo = asChar(uplo)) == NA_STRING ||
		(ul = *CHAR(uplo)) == '\0')
		error(_("invalid 'uplo' to 'R_dense_as_sparse()'"));
	    uplo = mkString((ul == 'U') ? "U" : "L");
	    if (c1 == 't') {
		if ((diag = asChar(diag)) == NA_STRING ||
		    (di = *CHAR(diag)) == '\0')
		    error(_("invalid 'diag' to 'R_dense_as_sparse()'"));
		diag = mkString((di == 'N') ? "N" : "U");
	    }
	    cl[1] = c1;
	}
#ifndef HAVE_PROPER_IMATRIX
	if (c0 == '.' && txf == INTSXP) {
	    PROTECT(x_from = coerceVector(x_from, txf = REALSXP));
	    ++nprotect;
	}
#endif
    }

    if (c0 != '.')
	cl[0] = c0;
    else
	txt = txf;
    cl[2] = c2;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)), p_to, i_to, j_to, x_to;
    ++nprotect;
    int m = pdim[0], n = pdim[1], i, j, *pp, *pi, *pj;
    R_xlen_t nnz = 0;

    SET_SLOT(to, Matrix_DimSym, dim);
    if (doDN)
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
    if (cl[1] != 'g')
	SET_SLOT(to, Matrix_uploSym, uplo);
    if (cl[1] == 't')
	SET_SLOT(to, Matrix_diagSym, diag);
    if (cl[2] != 'T') {
	PROTECT(p_to = allocVector(INTSXP,
				   (R_xlen_t) ((cl[2] == 'C') ? n : m) + 1));
	++nprotect;
	SET_SLOT(to, Matrix_pSym, p_to);
	pp = INTEGER(p_to);
	*(pp++) = 0;
	if (n > 0 && di != 'N' && ul == ((cl[2] == 'R') ? 'L' : 'U'))
	    *(pp++) = 0; /* first row or column skipped in these loops */
    }

#define DAS_LOOP_GE2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	for (j = 0; j < n; ++j) {				\
	    for (i = 0; i < m; ++i, ++_X_)			\
		if (_NZ_(*_X_)) _DO_INNER_;			\
	    _DO_OUTER_;						\
	}							\
    } while (0)
    
#define DAS_LOOP_GE2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t mn1s = (R_xlen_t) m * n - 1;			\
	for (i = 0; i < m; ++i, _X_ -= mn1s) {			\
	    for (j = 0; j < n; ++j, _X_ += m)			\
		if (_NZ_(*_X_)) _DO_INNER_;			\
	    _DO_OUTER_;						\
	}							\
    } while (0)
    
#define DAS_LOOP_TRN2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    for (j = 0; j < n; _X_ += n-(++j)) {		\
		for (i = 0; i <= j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; _X_ += (++j)) {			\
		for (i = j; i < n; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TRN2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = (R_xlen_t) n * n - 1;				\
	    for (i = 0; i < n; ++i, _X_ -= (d -= n)) {		\
		for (j = i; j < n; ++j, _X_ += n)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    d = -1;						\
	    for (i = 0; i < n; ++i, _X_ -= (d += n)) {		\
		for (j = 0; j <= i; ++j, _X_ += n)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)

#define DAS_LOOP_TRU2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    _X_ += n;						\
	    for (j = 1; j < n; ++j) {				\
		for (i = 0; i < j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
		_X_ += n-j;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; ++j) {				\
		_X_ += j+1;					\
		for (i = j+1; i < n; ++i, ++_X_)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TRU2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = (R_xlen_t) n * (n - 1) - 1;			\
	    for (i = 0; i < n; ++i) {				\
		for (j = i+1; j < n; ++j) {			\
		    _X_ += n;					\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d -= n);				\
	    }							\
	} else {						\
	    ++_X_;						\
	    d = -1;						\
	    for (i = 1; i < n; ++i) {				\
		for (j = 0; j < i; ++j) {			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		    _X_ += n;					\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d += n);				\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TPN2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    for (j = 0; j < n; ++j) {				\
		for (i = 0; i <= j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; ++j) {				\
		for (i = j; i < n; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)
    
#define DAS_LOOP_TPN2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = PM_LENGTH(n) - 1;				\
	    for (i = 0; i < n; _X_ -= (d -= (++i))) {		\
		for (j = i; j < n; _X_ += (++j))		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    d = -1;						\
	    for (i = 0; i < n; _X_ -= (d += n-(++i))) {		\
		for (j = 0; j <= i; _X_ += n-(++j))		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)

#define DAS_LOOP_TPU2C(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	if (ul == 'U') {					\
	    for (j = 1; j < n; ++j) {				\
		++_X_;						\
		for (i = 0; i < j; ++i, ++_X_)			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	} else {						\
	    for (j = 0; j < n; ++j) {				\
		++_X_;						\
		for (i = j+1; i < n; ++i, ++_X_)		\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		_DO_OUTER_;					\
	    }							\
	}							\
    } while (0)

#define DAS_LOOP_TPU2R(_X_, _NZ_, _DO_INNER_, _DO_OUTER_)	\
    do {							\
	R_xlen_t d;						\
	if (ul == 'U') {					\
	    d = PM_LENGTH(n-1) - 1;				\
	    for (i = 0; i < n; ++i) {				\
		for (j = i+1; j < n; ++j) {			\
		    _X_ += j;					\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d -= i+1);				\
	    }							\
	} else {						\
	    ++_X_;						\
	    d = -1;						\
	    for (i = 1; i < n; ++i) {				\
		for (j = 0; j < i; ++j) {			\
		    if (_NZ_(*_X_)) _DO_INNER_;			\
		    _X_ += n-j-1;				\
		}						\
		_DO_OUTER_;					\
		_X_ -= (d += n-i);				\
	    }							\
	}							\
    } while (0)
    
#define DAS_VALID2T						\
    if (nnz > INT_MAX)						\
	error(_("attempt to construct sparse matrix with "	\
		"more than 2^31-1 nonzero elements"))
    
#define DAS_VALID2C						\
    do { DAS_VALID2T; else *(pp++) = (int) nnz; } while (0)
    
#define DAS_SUBSUBCASES(_X_, _NZ_, _LOOP_2C_, _LOOP_2R_)	\
    do {							\
	switch (cl[2]) {					\
	case 'C':						\
	    _LOOP_2C_(_X_, _NZ_, ++nnz, DAS_VALID2C);		\
	    break;						\
	case 'R':						\
	    _LOOP_2R_(_X_, _NZ_, ++nnz, DAS_VALID2C);		\
	    break;						\
	default: /* 'T' */					\
	    _LOOP_2C_(_X_, _NZ_, ++nnz, DAS_VALID2T);		\
	    break;						\
	}							\
    } while (0)
    
#define DAS_SUBCASES(_CTYPE_, _PTR_, _NZ_)				\
    do {								\
	_CTYPE_ *px = _PTR_(x_from);					\
	if (cl[1] == 'g')						\
	    /* .geMatrix */						\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_GE2C,  DAS_LOOP_GE2R);	\
	else if (!packed && di == 'N')					\
	    /* .syMatrix, non-unit diagonal .trMatrix */		\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TRN2C, DAS_LOOP_TRN2R);	\
	else if (!packed)						\
	    /* unit diagonal .trMatrix */				\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TRU2C, DAS_LOOP_TRU2R);	\
	else if (di == 'N')						\
	    /* .spMatrix, non-unit diagonal .tpMatrix */		\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TPN2C, DAS_LOOP_TPN2R);	\
	else								\
	    /* unit diagonal .tpMatrix */				\
	    DAS_SUBSUBCASES(px, _NZ_, DAS_LOOP_TPU2C, DAS_LOOP_TPU2R);	\
    } while (0)
    
#define DAS_CASES(_SEXPTYPE_)						\
    do {								\
	switch (_SEXPTYPE_) {						\
	case REALSXP:							\
	    DAS_SUBCASES(double, REAL, NZ_REAL);			\
	    break;							\
	case LGLSXP:							\
	    DAS_SUBCASES(int, LOGICAL, NZ_INTEGER);			\
	    break;							\
	case INTSXP:							\
	    DAS_SUBCASES(int, INTEGER, NZ_INTEGER);			\
	    break;							\
	case CPLXSXP:							\
	    DAS_SUBCASES(Rcomplex, COMPLEX, NZ_COMPLEX);		\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", _SEXPTYPE_, "R_dense_as_sparse"); \
	    break;							\
	}								\
    } while (0)

    /* First we loop over the _nontrivial part_ of the denseMatrix 'from',
       by row ('R' case) or by column ('C' and 'T' cases), counting the
       nonzero elements and filling the 'p' slot of the result accordingly 
       ('C' and 'R' cases) ... */
    DAS_CASES(txf);

#undef DAS_SUBCASES
#undef DAS_SUBSUBCASES
#undef DAS_VALID2C
#undef DAS_VALID2T

    /* Then we allocate ... */
    if (cl[2] != 'R') {
	PROTECT(i_to = allocVector(INTSXP, nnz));
	++nprotect;
	SET_SLOT(to, Matrix_iSym, i_to);
	pi = INTEGER(i_to);
    }
    if (cl[2] != 'C') {
	PROTECT(j_to = allocVector(INTSXP, nnz));
	++nprotect;
	SET_SLOT(to, Matrix_jSym, j_to);
	pj = INTEGER(j_to);
    }
    if (cl[0] != 'n') {
	PROTECT(x_to = allocVector(txf, nnz));
	++nprotect;
	if (txf == txt)
	    SET_SLOT(to, Matrix_xSym, x_to);
    }

#define DAS_SUBSUBCASES(_X_, _Y_, _NZ_, _LOOP_2C_, _LOOP_2R_)		\
    do {								\
	switch (cl[2]) {						\
	case 'C':							\
	    if (cl[0] == 'n')						\
		_LOOP_2C_(_X_, _NZ_,					\
			  *(pi++) = i, );				\
	    else							\
		_LOOP_2C_(_X_, _NZ_,					\
			  do {						\
			      *(pi++) = i;				\
			      *(_Y_++) = *_X_;				\
			  } while (0), );				\
	    break;							\
	case 'R':							\
	    if (cl[0] == 'n')						\
		_LOOP_2R_(_X_, _NZ_,					\
			  *(pj++) = j, );				\
	    else							\
		_LOOP_2R_(_X_, _NZ_,					\
			  do {						\
			      *(pj++) = j;				\
			      *(_Y_++) = *_X_;				\
			  } while (0), );				\
	    break;							\
	default: /* 'T' */						\
	    if (cl[0] == 'n')						\
		_LOOP_2C_(_X_, _NZ_,					\
			  do {						\
			      *(pi++) = i;				\
			      *(pj++) = j;				\
			  } while (0), );				\
	    else							\
		_LOOP_2C_(_X_, _NZ_,					\
			  do {						\
			      *(pi++) = i;				\
			      *(pj++) = j;				\
			      *(_Y_++) = *_X_;				\
			  } while (0), );				\
	    break;							\
	}								\
    } while (0)
    
#define DAS_SUBCASES(_CTYPE_, _PTR_, _NZ_)				\
    do {								\
	_CTYPE_ *px = _PTR_(x_from), *py;				\
	if (cl[0] != 'n')						\
	    py = _PTR_(x_to);						\
	if (cl[1] == 'g')						\
	    /* .geMatrix */						\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_GE2C,  DAS_LOOP_GE2R); \
	else if (!packed && di == 'N')					\
	    /* .syMatrix, non-unit diagonal .trMatrix */		\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TRN2C, DAS_LOOP_TRN2R); \
	else if (!packed)						\
	    /* unit diagonal .trMatrix */				\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TRU2C, DAS_LOOP_TRU2R); \
	else if (di == 'N')						\
	    /* .spMatrix, non-unit diagonal .tpMatrix */		\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TPN2C, DAS_LOOP_TPN2R); \
	else								\
	    /* unit diagonal .tpMatrix */				\
	    DAS_SUBSUBCASES(px, py, _NZ_, DAS_LOOP_TPU2C, DAS_LOOP_TPU2R); \
    } while (0)

    /* Then we loop back over the same elements in order to fill
       the 'i', 'j', and 'x' slots of the result (whichever exist) ... */
    DAS_CASES(txf);

#undef DAS_CASES
#undef DAS_SUBCASES
#undef DAS_SUBSUBCASES
#undef DAS_LOOP_GE2C
#undef DAS_LOOP_GE2R
#undef DAS_LOOP_TRN2C
#undef DAS_LOOP_TRN2R
#undef DAS_LOOP_TRU2C
#undef DAS_LOOP_TRU2R
#undef DAS_LOOP_TPN2C
#undef DAS_LOOP_TPN2R
#undef DAS_LOOP_TPU2C
#undef DAS_LOOP_TPU2R

    if (cl[0] != 'n' && txf != txt)
	SET_SLOT(to, Matrix_xSym, coerceVector(x_to, txt));
    if (ivalid >= 0)
	free(cl);
    UNPROTECT(nprotect);
    return to;
}

/* as(<denseMatrix>, "[dlniz](dense)?Matrix") */
SEXP R_dense_as_kind(SEXP from, SEXP kind)
{
    char k;
    if ((kind = asChar(kind)) == NA_STRING || (k = *CHAR(kind)) == '\0')
	error(_("invalid 'kind' to 'R_dense_as_kind()'"));
    static const char *valid[] = {
	"dgeMatrix", "dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix",
	"lgeMatrix", "ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix",
	"ngeMatrix", "ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid < 0)
	ERROR_INVALID_CLASS(class_P(from), "R_dense_as_kind");
    const char *clf = valid[ivalid];
    if (k == '.')
	k = clf[0];
    if (k == clf[0])
	return from;
    SEXPTYPE tt = kind2type(k); /* validating 'k' before doing more */
    
    char *clt = strdup(clf);
    clt[0] = k;
    SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(clt)),
	x = GET_SLOT(from, Matrix_xSym);
    SEXPTYPE tf = TYPEOF(x);
    free(clt);
    
    SET_SLOT(to, Matrix_DimSym, GET_SLOT(from, Matrix_DimSym));
    SET_SLOT(to, Matrix_DimNamesSym, GET_SLOT(from, Matrix_DimNamesSym));
    SET_SLOT(to, Matrix_xSym, (tf == tt) ? x : coerceVector(x, tt));
    if (clf[1] != 'g') {
	SET_SLOT(to, Matrix_uploSym, GET_SLOT(from, Matrix_uploSym));
	if (clf[1] == 't')
	    SET_SLOT(to, Matrix_diagSym, GET_SLOT(from, Matrix_diagSym));
    }
    UNPROTECT(1);
    return to;
}

/* as(<denseMatrix>, "matrix") */
SEXP R_dense_as_matrix(SEXP from)
{
    PROTECT(from = dense_as_geMatrix(from, '.', 1, 0));
    SEXP to = PROTECT(GET_SLOT(from, Matrix_xSym));
    setAttrib(to, R_DimSymbol, GET_SLOT(from, Matrix_DimSym));
    setAttrib(to, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(2);
    return to;
}

/* as(<.geMatrix>, "matrix") */
SEXP R_geMatrix_as_matrix(SEXP from)
{
    SEXP to = PROTECT(duplicate(GET_SLOT(from, Matrix_xSym)));
    setAttrib(to, R_DimSymbol, GET_SLOT(from, Matrix_DimSym));
    setAttrib(to, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(1);
    return to;
}

/* band(<denseMatrix>, k1, k2), tri[ul](<denseMatrix>, k)
   band(     <matrix>, k1, k2), tri[ul](     <matrix>, k) */
SEXP R_dense_band(SEXP from, SEXP k1, SEXP k2)
{
    const char *clf;
    static const char *valid[] = {
	"dgeMatrix", "dtrMatrix", "dsyMatrix", "dtpMatrix", "dspMatrix",
	"lgeMatrix", "ltrMatrix", "lsyMatrix", "ltpMatrix", "lspMatrix",
	"ngeMatrix", "ntrMatrix", "nsyMatrix", "ntpMatrix", "nspMatrix", ""};
    int ivalid = R_check_class_etc(from, valid);
    if (ivalid >= 0) {
	clf = valid[ivalid];
    } else {
	/* matrix->.geMatrix with unreferenced 'x' slot ... modify directly */
	PROTECT(from = matrix_as_geMatrix(from, '.', 0, 0));
	clf = class_P(from);
    }
    
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
    
    SEXP uplo_from, diag;
    char ulf = 'U', ult = 'U', di = 'N';
    if (clf[1] != 'g') {
	uplo_from = GET_SLOT(from, Matrix_uploSym);
	ulf = *CHAR(STRING_ELT(uplo_from, 0));
	if (clf[1] == 't') {
	    if ((ulf == 'U') ? (a <= 0 && b >= n - 1) : (b >= 0 && a <= 1 - m))
		return from;
	    diag = GET_SLOT(from, Matrix_diagSym);
	    di = *CHAR(STRING_ELT(diag, 0));
	}
    }
    
    SEXP to, x_to,
	x_from = GET_SLOT(from, Matrix_xSym),
	dimnames = GET_SLOT(from, Matrix_DimNamesSym);
    SEXPTYPE tx = TYPEOF(x_from);
    Rboolean tr, sy;
    
#define UNPACKED_MAKE_BANDED(_PREFIX_, _PTR_)				\
    _PREFIX_ ## dense_unpacked_make_banded(_PTR_(x_to), m, n, a, b, di)

#define PACKED_MAKE_BANDED(_PREFIX_, _PTR_)				\
    _PREFIX_ ## dense_packed_make_banded(_PTR_(x_to), n, a, b, ult, di)
    
#define DENSE_BAND(_MAKE_BANDED_)					\
    do {								\
	switch (tx) {							\
	case REALSXP: /* d..Matrix */					\
	    _MAKE_BANDED_(d, REAL);					\
	    break;							\
	case LGLSXP: /* [ln]..Matrix */					\
	    _MAKE_BANDED_(i, LOGICAL);					\
	    break;							\
	case INTSXP: /* i..Matrix */					\
	    _MAKE_BANDED_(i, INTEGER);					\
	    break;							\
	case CPLXSXP: /* z..Matrix */					\
	    _MAKE_BANDED_(z, COMPLEX);					\
	    break;							\
	default:							\
	    ERROR_INVALID_TYPE("'x' slot", tx, "R_dense_band");		\
	    break;							\
	}								\
    } while (0)

#define DENSE_BAND_CONSTRICT			\
    do {					\
	if (ulf == 'U') {			\
	    if (a < 0) a = 0;			\
	} else {				\
	    if (b > 0) b = 0;			\
	}					\
    } while(0)
    
    if (m != n || (!(tr = a >= 0 || b <= 0 || clf[1] == 't') &&
		   !(sy = a == -b && clf[1] == 's'))) { /* return .geMatrix */
	
	if (clf[1] == 'g') {
	    if (ivalid >= 0) {
		PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
		PROTECT(x_to = duplicate(x_from));
		SET_SLOT(to, Matrix_DimSym, dim);
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
		SET_SLOT(to, Matrix_xSym, x_to);
	    } else {
		to = from;
		PROTECT(x_to = x_from);
	    }
	} else { /* clf[1] == 's' */
	    PROTECT(to = dense_as_geMatrix(from, '.', 1, 0));
	    PROTECT(x_to = GET_SLOT(to, Matrix_xSym));
	}
	DENSE_BAND(UNPACKED_MAKE_BANDED);
	UNPROTECT(2);
	return to;

    } else if (tr) { /* return .t.Matrix */
	
	ult = (a >= 0) ? 'U' : ((b <= 0) ? 'L' : ulf);

	if (clf[1] == 't') {
	    PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
	    SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	    if (di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	    
	    DENSE_BAND_CONSTRICT;
	    if (ulf == ult) {
		PROTECT(x_to = duplicate(x_from));
		if (clf[2] == 'p')
		    DENSE_BAND(PACKED_MAKE_BANDED);
		else
		    DENSE_BAND(UNPACKED_MAKE_BANDED);
	    } else {
		R_xlen_t nx = XLENGTH(x_from);
		PROTECT(x_to = allocVector(tx, nx));
		
#define PACKED_COPY_DIAGONAL(_PREFIX_, _PTR_)				\
		do {							\
		    Memzero(_PTR_(x_to), nx);				\
		    if (a <= b)						\
			_PREFIX_ ## dense_packed_copy_diagonal(		\
			    _PTR_(x_to), _PTR_(x_from),			\
			    n, nx, ult, ulf, di);			\
		} while (0)
		
#define UNPACKED_COPY_DIAGONAL(_PREFIX_, _PTR_)				\
		do {							\
		    Memzero(_PTR_(x_to), nx);				\
		    if (a <= b)						\
			_PREFIX_ ## dense_unpacked_copy_diagonal(	\
			    _PTR_(x_to), _PTR_(x_from),			\
			    n, nx, 'U' /* unused */, di);		\
		} while (0)
		
		if (clf[2] == 'p')
		    DENSE_BAND(PACKED_COPY_DIAGONAL);
		else
		    DENSE_BAND(UNPACKED_COPY_DIAGONAL);
	    }
	    
	} else { /* clf[1] == 'g' || clf[1] == 's' */
	    char *clt = strdup(clf);
	    clt[1] = 't';
	    if (clf[2] != 'p')
		clt[2] = 'r';
	    PROTECT(to = NEW_OBJECT_OF_CLASS(clt));
	    free(clt);
	    if (clf[1] == 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	    else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	    
	    PROTECT(x_to = (clf[1] == 'g' || n <= 1
			    ? (ivalid >= 0
			       ? duplicate(x_from)
			       : x_from)
			    : (ulf == ult
			       ? duplicate(x_from)
			       : (clf[2] == 'p'
				  ? packed_transpose(x_from, n, ulf)
				  : unpacked_force(x_from, n, ulf, '\0')))));
	    
	    if (clf[2] == 'p')
		DENSE_BAND(PACKED_MAKE_BANDED);
	    else
		DENSE_BAND(UNPACKED_MAKE_BANDED);
	}
	
	SET_SLOT(to, Matrix_uploSym, mkString((ult == 'U') ? "U" : "L"));

    } else if (sy) { /* return .s.Matrix */

	PROTECT(to = NEW_OBJECT_OF_CLASS(clf));
	PROTECT(x_to = duplicate(x_from));
	ult = ulf;
	DENSE_BAND_CONSTRICT;
	
	if (clf[2] == 'p')
	    DENSE_BAND(PACKED_MAKE_BANDED);
	else
	    DENSE_BAND(UNPACKED_MAKE_BANDED);
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	SET_SLOT(to, Matrix_uploSym, uplo_from);
	
    } else {
	
	error(_("should never happen; probably a bug in 'R_dense_band()'; "
		"please report!"));
	
    }

#undef PACKED_COPY_DIAGONAL
#undef PACKED_MAKE_BANDED
#undef UNPACKED_COPY_DIAGONAL
#undef UNPACKED_MAKE_BANDED
#undef DENSE_BAND_CONSTRICT
#undef DENSE_BAND

    SET_SLOT(to, Matrix_DimSym, dim);
    SET_SLOT(to, Matrix_xSym, x_to);
    UNPROTECT((ivalid >= 0) ? 2 : 3);
    return to;
}

/**
 * Perform a left cyclic shift of columns j to k in the upper triangular
 * matrix x, then restore it to upper triangular form with Givens rotations.
 * The algorithm is based on the Fortran routine DCHEX from Linpack.
 *
 * The lower triangle of x is not modified.
 *
 * @param x Matrix stored in column-major order
 * @param ldx leading dimension of x
 * @param j column number (0-based) that will be shifted to position k
 * @param k last column number (0-based) to be shifted
 * @param cosines cosines of the Givens rotations
 * @param sines sines of the Givens rotations
 *
 * @return 0 for success
 */
static
int left_cyclic(double x[], int ldx, int j, int k,
		double cosines[], double sines[])
{
    if (j >= k)
	error(_("incorrect left cyclic shift, j (%d) >= k (%d)"), j, k);
    if (j < 0)
	error(_("incorrect left cyclic shift, j (%d) < 0"), j, k);
    if (ldx < k)
	error(_("incorrect left cyclic shift, k (%d) > ldx (%d)"), k, ldx);

    double *lastcol = (double*) R_alloc(k+1, sizeof(double));
    int i;
				/* keep a copy of column j */
    for(i = 0; i <= j; i++) lastcol[i] = x[i + j*ldx];
				/* For safety, zero the rest */
    for(i = j+1; i <= k; i++) lastcol[i] = 0.;
    for(int jj = j+1, ind = 0; jj <= k; jj++, ind++) { /* columns to be shifted */
	int diagind = jj*(ldx+1); //  ind == (jj-j) - 1
	double tmp = x[diagind], cc, ss;
	/* Calculate the Givens rotation. */
				/* This modified the super-diagonal element */
	F77_CALL(drotg)(x + diagind-1, &tmp, cosines + ind, sines + ind);
	cc = cosines[ind]; ss = sines[ind];
				/* Copy column jj+1 to column jj. */
	for(i = 0; i < jj; i++) x[i + (jj-1)*ldx] = x[i+jj*ldx];
				/* Apply rotation to columns up to k */
	for(i = jj; i < k; i++) {
	    tmp = cc*x[(jj-1)+i*ldx] + ss*x[jj+i*ldx];
	    x[jj+i*ldx] = cc*x[jj+i*ldx] - ss*x[(jj-1)+i*ldx];
	    x[(jj-1)+i*ldx] = tmp;
	}
				/* Apply rotation to lastcol */
	lastcol[jj] = -ss*lastcol[jj-1]; lastcol[jj-1] *= cc;
    }
				/* Copy lastcol to column k */
    for(i = 0; i <= k; i++) x[i+k*ldx] = lastcol[i];
    return 0;
}

static
SEXP getGivens(double x[], int ldx, int jmin, int rank)
{
    int shiftlen = (rank - jmin) - 1;
    SEXP ans = PROTECT(allocVector(VECSXP, 4)), nms, cosines, sines;

    SET_VECTOR_ELT(ans, 0, ScalarInteger(jmin));
    SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
    SET_VECTOR_ELT(ans, 2, cosines = allocVector(REALSXP, shiftlen));
    SET_VECTOR_ELT(ans, 3,   sines = allocVector(REALSXP, shiftlen));
    setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 4));
    SET_STRING_ELT(nms, 0, mkChar("jmin"));
    SET_STRING_ELT(nms, 1, mkChar("rank"));
    SET_STRING_ELT(nms, 2, mkChar("cosines"));
    SET_STRING_ELT(nms, 3, mkChar("sines"));
    if (left_cyclic(x, ldx, jmin, rank - 1, REAL(cosines), REAL(sines)))
	error(_("Unknown error in getGivens"));
    UNPROTECT(1);
    return ans;
}

SEXP checkGivens(SEXP X, SEXP jmin, SEXP rank)
{
    SEXP ans = PROTECT(allocVector(VECSXP, 2)),
	Xcp = PROTECT(duplicate(X));
    int  *Xdims;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
    SET_VECTOR_ELT(ans, 1, getGivens(REAL(Xcp), Xdims[0],
				     asInteger(jmin), asInteger(rank)));
    SET_VECTOR_ELT(ans, 0, Xcp);
    UNPROTECT(2);
    return ans;
}

SEXP lsq_dense_Chol(SEXP X, SEXP y)
{
    SEXP ans;
    double d_one = 1., d_zero = 0.;

    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    int *Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP)),
	n = Xdims[0],
	p = Xdims[1];
    if (!(isReal(y) & isMatrix(y)))
	error(_("y must be a numeric (double precision) matrix"));
    int *ydims = INTEGER(coerceVector(getAttrib(y, R_DimSymbol), INTSXP));
    if (ydims[0] != n)
	error(_(
	    "number of rows in y (%d) does not match number of rows in X (%d)"),
	    ydims[0], n);
    int k = ydims[1];
    if (k < 1 || p < 1) return allocMatrix(REALSXP, p, k);
    ans = PROTECT(allocMatrix(REALSXP, p, k));
    F77_CALL(dgemm)("T", "N", &p, &k, &n, &d_one, REAL(X), &n, REAL(y), &n,
		    &d_zero, REAL(ans), &p FCONE FCONE);
    double *xpx = (double *) R_alloc(p * p, sizeof(double));
    F77_CALL(dsyrk)("U", "T", &p, &n, &d_one, REAL(X), &n, &d_zero,
		    xpx, &p FCONE FCONE);
    int info;
    F77_CALL(dposv)("U", &p, &k, xpx, &p, REAL(ans), &p, &info FCONE);
    if (info) error(_("Lapack routine dposv returned error code %d"), info);
    UNPROTECT(1);
    return ans;
}


SEXP lsq_dense_QR(SEXP X, SEXP y)
{
    if (!(isReal(X) & isMatrix(X)))
	error(_("X must be a numeric (double precision) matrix"));
    int *Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP)),
	n = Xdims[0],
	p = Xdims[1];
    if (!(isReal(y) & isMatrix(y)))
	error(_("y must be a numeric (double precision) matrix"));
    int *ydims = INTEGER(coerceVector(getAttrib(y, R_DimSymbol), INTSXP));
    if (ydims[0] != n)
	error(_(
	    "number of rows in y (%d) does not match number of rows in X (%d)"),
	    ydims[0], n);
    int k = ydims[1];
    if (k < 1 || p < 1) return allocMatrix(REALSXP, p, k);
    double tmp, *xvals = (double *) Memcpy(R_alloc(n * p, sizeof(double)), REAL(X), n * p);
    SEXP ans = PROTECT(duplicate(y));
    int lwork = -1, info;
    F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
		    &tmp, &lwork, &info FCONE);
    if (info)
	error(_("First call to Lapack routine dgels returned error code %d"),
	      info);
    lwork = (int) tmp;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgels)("N", &n, &p, &k, xvals, &n, REAL(ans), &n,
		    work, &lwork, &info FCONE);
    if (info)
	error(_("Second call to Lapack routine dgels returned error code %d"),
	      info);
    UNPROTECT(1);
    return ans;
}

/* Rank-Correcting/Adapting LAPACK  QR Decomposition
 * From Doug Bates' initial import; __unused__

 * Provides a qr() with 'rcond' and rank reduction while(rcond < tol),
 * possibly via Givens rotations but WITHOUT PIVOTING

 * .Call(Matrix:::lapack_qr, A, 1e-17) --> ~/R/MM/Pkg-ex/Matrix/qr-rank-deficient.R

 * TODO: export as Matrix::qrNoPiv() or qr1()  or similar
 */
SEXP lapack_qr(SEXP Xin, SEXP tl)
{
    if (!(isReal(Xin) & isMatrix(Xin)))
	error(_("X must be a real (numeric) matrix"));
    double tol = asReal(tl);
    if (tol < 0.) error(_("tol, given as %g, must be non-negative"), tol);
    if (tol > 1.) error(_("tol, given as %g, must be <= 1"), tol);
    SEXP ans = PROTECT(allocVector(VECSXP,5)), X, qraux, pivot;
    SET_VECTOR_ELT(ans, 0, X = duplicate(Xin));
    int *Xdims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP)),
	n = Xdims[0], i,
	p = Xdims[1],
	trsz = (n < p) ? n : p ; /* size of triangular part of decomposition */

    SET_VECTOR_ELT(ans, 2, qraux = allocVector(REALSXP, trsz));
    SET_VECTOR_ELT(ans, 3, pivot = allocVector(INTSXP, p));
    for (i = 0; i < p; i++) INTEGER(pivot)[i] = i + 1;
    SEXP nms,
	Givens = PROTECT(allocVector(VECSXP, trsz - 1));
    setAttrib(ans, R_NamesSymbol, nms = allocVector(STRSXP, 5));
    SET_STRING_ELT(nms, 0, mkChar("qr"));
    SET_STRING_ELT(nms, 1, mkChar("rank"));
    SET_STRING_ELT(nms, 2, mkChar("qraux"));
    SET_STRING_ELT(nms, 3, mkChar("pivot"));
    SET_STRING_ELT(nms, 4, mkChar("Givens"));
    int rank = trsz,
	nGivens = 0;
    double rcond = 0.;
    if (n > 0 && p > 0) {
	int  info, *iwork, lwork;
	double *xpt = REAL(X), *work, tmp;

	lwork = -1;
	F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), &tmp, &lwork, &info);
	if (info)
	    error(_("First call to dgeqrf returned error code %d"), info);
	lwork = (int) tmp;
	work = (double *) R_alloc((lwork < 3*trsz) ? 3*trsz : lwork,
				  sizeof(double));
	F77_CALL(dgeqrf)(&n, &p, xpt, &n, REAL(qraux), work, &lwork, &info);
	if (info)
	    error(_("Second call to dgeqrf returned error code %d"), info);
	iwork = (int *) R_alloc(trsz, sizeof(int));
	F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			 work, iwork, &info FCONE FCONE FCONE);
	if (info)
	    error(_("Lapack routine dtrcon returned error code %d"), info);
	while (rcond < tol) {	/* check diagonal elements */
	    double minabs = (xpt[0] < 0.) ? -xpt[0]: xpt[0];
	    int jmin = 0;
	    for (i = 1; i < rank; i++) {
		double el = xpt[i*n]; // had  i*(n+1)  which looks wrong to MM
		if(el < 0.) el = -el;
		if (el < minabs) {
		    jmin = i;
		    minabs = el;
		}
	    }
	    if (jmin < (rank - 1)) {
		SET_VECTOR_ELT(Givens, nGivens, getGivens(xpt, n, jmin, rank));
		nGivens++;
	    } // otherwise jmin == (rank - 1) , so just "drop that column"
	    rank--;
	    // new  rcond := ... for reduced rank
	    F77_CALL(dtrcon)("1", "U", "N", &rank, xpt, &n, &rcond,
			     work, iwork, &info FCONE FCONE FCONE);
	    if (info)
		error(_("Lapack routine dtrcon returned error code %d"), info);
	}
    }
    SEXP Gcpy, sym;
    SET_VECTOR_ELT(ans, 4, Gcpy = allocVector(VECSXP, nGivens));
    for (i = 0; i < nGivens; i++)
	SET_VECTOR_ELT(Gcpy, i, VECTOR_ELT(Givens, i));
    SET_VECTOR_ELT(ans, 1, ScalarInteger(rank));
    sym = PROTECT(install("useLAPACK")); setAttrib(ans, sym, ScalarLogical(1)); UNPROTECT(1);
    sym = PROTECT(install("rcond"));     setAttrib(ans, sym, ScalarReal(rcond));UNPROTECT(1);
    UNPROTECT(2);
    return ans;
}

/* MJ: no longer needed ... prefer R_dense_as_sparse() above */
#if 0

SEXP dense_to_Csparse(SEXP x)
{
    SEXP ge_x = PROTECT(dense_as_geMatrix(x, '.', 2, 0)),
	Dim = GET_SLOT(ge_x, Matrix_DimSym);
    int *dims = INTEGER(Dim);
    Rboolean longi = (dims[0] * (double)dims[1] > INT_MAX);
    // int itype = longi ? CHOLMOD_LONG : CHOLMOD_INT;
    CHM_DN chxd = AS_CHM_xDN(ge_x); // cholmod_dense (has no itype)
    CHM_SP chxs;
    /* cholmod_dense_to_sparse() in CHOLMOD/Core/ below does only work for
       "REAL" 'xtypes', i.e. *not* for "nMatrix".
       ===> need "_x" in above AS_CHM_xDN() call.

       Also it cannot keep symmetric / triangular, hence the
       as_geMatrix() above.  Note that this is already a *waste* for
       symmetric matrices; However, we could conceivably use an
       enhanced cholmod_dense_to_sparse(), with an extra boolean
       argument for symmetry.
    */
#define DLONG
/* You can try defining DLONG -- then just get a seg.fault :
 *  I think it is because of this in  ./CHOLMOD/Include/cholmod_core.h :
 *
 The itype of all parameters for all CHOLMOD routines must match.
 --- ^^^^^ ------------------------------------------------------
 but then as_cholmod_dense should *not* make a difference: cholmod_dense has *no* itype   (????)
*/
    if(longi) { // calling cholmod_dense_to_sparse() gives wrong matrix
#ifdef DLONG
	chxs = cholmod_l_dense_to_sparse(chxd, 1, &cl);
	// in gdb, I found that 'chxs' seems "basically empty": all
	// p chxs->foo   give ''Cannot access memory at address 0x....''
	// for now rather give error:
	if(cl.status)
	    error(_("dense_to_Csparse(<LARGE>): cholmod_l_dense_to_sparse failure status=%d"),
		  cl.status);
#else
        error(_("Matrix dimension %d x %d (= %g) too large [FIXME calling cholmod_l_dense_to_sparse]"),
	      m,n, m * (double)n);
#endif
    } else { // fits, using integer (instead of long int) 'itype'
	chxs = cholmod_dense_to_sparse(chxd, 1, &c);
    }

    int Rkind = (chxd->xtype == CHOLMOD_REAL) ? Real_KIND2(x) : 0;
    /* Note: when 'x' was integer Matrix, Real_KIND(x) = -1, but *_KIND2(.) = 0 */
    R_CheckStack();

    UNPROTECT(1);
    /* chm_sparse_to_SEXP() *could* deal with symmetric
     * if chxs had such an stype; and we should be able to use uplo below */
    return chm_sparse_to_SEXP(chxs, 1, 0/*TODO: uplo_P(x) if x has an uplo slot*/,
			      Rkind, "",
			      isMatrix(x) ? getAttrib(x, R_DimNamesSymbol)
			      : GET_SLOT(x, Matrix_DimNamesSym));
}

#endif /* MJ */

SEXP ddense_symmpart(SEXP x)
/* Class of the value will be dsyMatrix */
{
    SEXP dx = PROTECT(dense_as_geMatrix(x, 'd', 2, 0));
    int *adims = INTEGER(GET_SLOT(dx, Matrix_DimSym)), n = adims[0];

    if(n != adims[1]) {
	error(_("matrix is not square! (symmetric part)"));
	return R_NilValue; /* -Wall */
    } else {
	SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dsyMatrix"));
	double *xx = REAL(GET_SLOT(dx, Matrix_xSym));

	/* only need to assign the *upper* triangle (uplo = "U");
	 * noting that diagonal remains unchanged */
	R_xlen_t n_ = n;
	for (int j = 0; j < n; j++) {
	    for (int i = 0; i < j; i++) {
		xx[j * n_ + i] = (xx[j * n_ + i] + xx[i * n_ + j]) / 2.;
	    }
	}

/* Copy dx to ans:
   Because slots of dx are freshly allocated and dx will not
   be used, we use the slots themselves and don't duplicate */	
#define MK_SYMMETRIC_DIMNAMES_AND_RETURN(_J_ /* -1|0|1 */)		\
	SET_SLOT(ans, Matrix_DimSym,  GET_SLOT(dx, Matrix_DimSym));	\
	set_symmetrized_DimNames(ans, GET_SLOT(dx, Matrix_DimNamesSym), _J_); \
	SET_SLOT(ans, Matrix_uploSym, mkString((_J_) ? "U" : "L"));	\
	SET_SLOT(ans, Matrix_xSym,    GET_SLOT(dx, Matrix_xSym));	\
	UNPROTECT(2);							\
	return ans

        MK_SYMMETRIC_DIMNAMES_AND_RETURN(-1);
    }
}

SEXP ddense_skewpart(SEXP x)
/* Class of the value will be dgeMatrix */
{
    SEXP dx = PROTECT(dense_as_geMatrix(x, 'd', 2, 0));
    int *adims = INTEGER(GET_SLOT(dx, Matrix_DimSym)), n = adims[0];

    if(n != adims[1]) {
	error(_("matrix is not square! (skew-symmetric part)"));
	return R_NilValue; /* -Wall */
    } else {
	SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix"));
	double *xx = REAL(GET_SLOT(dx, Matrix_xSym));
	R_xlen_t n_ = n;

	for (int j = 0; j < n_; j++) {
	    xx[j * n_ + j] = 0.;
	    for (int i = 0; i < j; i++) {
		double s = (xx[j * n_ + i] - xx[i * n_ + j]) / 2.;
		xx[j * n_ + i] =  s;
		xx[i * n_ + j] = -s;
	    }
	}

        MK_SYMMETRIC_DIMNAMES_AND_RETURN(-1);
    }
}

/* MJ: no longer needed ... prefer (un)?packedMatrix_force_symmetric() */
#if 0

SEXP dense_to_symmetric(SEXP x, SEXP uplo, SEXP symm_test)
/* Class of result will be [dln]syMatrix */
{
/*== FIXME: allow  uplo = NA   and then behave a bit like symmpart():
 *== -----  would use the *dimnames* to determine U or L   (??)
 */
    
    int symm_tst = asLogical(symm_test);
    SEXP dx = PROTECT(dense_as_geMatrix(x, '.', 2, 0));
    SEXP ans;
    const char *cl = class_P(dx);
    /* same as in ..._geMatrix() above:*/
    enum dense_enum M_type = ( (cl[0] == 'd') ? ddense :
			       ((cl[0] == 'l') ? ldense : ndense) );
    int *adims = INTEGER(GET_SLOT(dx, Matrix_DimSym)), n = adims[0];
    if(n != adims[1]) {
	UNPROTECT(1);
	error(_("dense_to_symmetric(): matrix is not square!"));
	return R_NilValue; /* -Wall */
    }

    if(symm_tst) {
	int i, j;
	R_xlen_t n_ = n;
	
#define CHECK_SYMMETRIC							\
	do {								\
	    for (j = 0; j < n; j++) {					\
		for (i = 0; i < j; i++)	{				\
		    if(xx[j * n_ + i] != xx[i * n_ + j]) {		\
			UNPROTECT(1);					\
			error(_("matrix is not symmetric [%d,%d]"), i+1, j+1); \
			return R_NilValue; /* -Wall */			\
		    }							\
		}							\
	    }								\
	} while (0)
	
	if(M_type == ddense) {
	    double *xx = REAL(GET_SLOT(dx, Matrix_xSym));
	    CHECK_SYMMETRIC;
	} else { /* (M_type == ldense || M_type == ndense) */
	    int *xx = LOGICAL(GET_SLOT(dx, Matrix_xSym));
	    CHECK_SYMMETRIC;
	}
    }
#undef CHECK_SYMMETRIC
    
    ans = PROTECT(NEW_OBJECT_OF_CLASS(M_type == ddense
				      ? "dsyMatrix"
				      : (M_type == ldense
					 ? "lsyMatrix"
					 : "nsyMatrix")));
    int uploT = (*CHAR(asChar(uplo)) == 'U');
    MK_SYMMETRIC_DIMNAMES_AND_RETURN(uploT);
}

#endif /* MJ */

/* MJ: no longer needed ... prefer R_dense_band() above */
#if 0

SEXP dense_band(SEXP x, SEXP k1P, SEXP k2P)
/* Always returns a full matrix with entries outside the band zeroed
 * Class of the value can be [dln]trMatrix or [dln]geMatrix
 */
{
    int k1 = asInteger(k1P), k2 = asInteger(k2P);

    if (k1 > k2) {
	error(_("Lower band %d > upper band %d"), k1, k2);
	return R_NilValue; /* -Wall */
    }
    else {
	SEXP ans = PROTECT(dense_as_geMatrix(x, '.', 2, 0));
	int *adims = INTEGER(GET_SLOT(ans, Matrix_DimSym)),
	    j, m = adims[0], n = adims[1],
	    sqr = (adims[0] == adims[1]),
	    tru = (k1 >= 0), trl = (k2 <= 0);
	const char *cl = class_P(ans);
	enum dense_enum M_type = ( (cl[0] == 'd') ? ddense :
			      ((cl[0] == 'l') ? ldense : ndense));


#define SET_ZERO_OUTSIDE				\
	for (j = 0; j < n; j++) {			\
	    int i, i1 = j - k2, i2 = j + 1 - k1;	\
	    R_xlen_t jm = j * (R_xlen_t) m;		\
	    if(i1 > m) i1 = m;				\
	    if(i2 < 0) i2 = 0;				\
	    for (i = 0; i < i1; i++) xx[i + jm] = 0;	\
	    for (i = i2; i < m; i++) xx[i + jm] = 0;	\
	}

	if(M_type == ddense) {
	    double *xx = REAL(GET_SLOT(ans, Matrix_xSym));
	    SET_ZERO_OUTSIDE
	}
	else { /* (M_type == ldense || M_type == ndense) */
	    int *xx = LOGICAL(GET_SLOT(ans, Matrix_xSym));
	    SET_ZERO_OUTSIDE
	}

	if (!sqr || (!tru && !trl)) { /* return the *geMatrix */
	    UNPROTECT(1);
	    return ans;
	}
	else {
	    /* Copy ans to a *trMatrix object (must be square) */
	    SEXP aa = PROTECT(NEW_OBJECT_OF_CLASS(M_type == ddense
						  ? "dtrMatrix"
						  : (M_type == ldense
						     ? "ltrMatrix"
						     : "ntrMatrix")));
	    /* Because slots of ans are freshly allocated and ans will not be
	     * used, we use the slots themselves and don't duplicate */
	    SET_SLOT(aa, Matrix_xSym,        GET_SLOT(ans, Matrix_xSym));
	    SET_SLOT(aa, Matrix_DimSym,      GET_SLOT(ans, Matrix_DimSym));
	    SET_SLOT(aa, Matrix_DimNamesSym, GET_SLOT(ans, Matrix_DimNamesSym));
	    SET_SLOT(aa, Matrix_diagSym,     mkString("N"));
	    SET_SLOT(aa, Matrix_uploSym,     mkString(tru ? "U" : "L"));
	    UNPROTECT(2);
	    return aa;
	}
    }
}

#endif /* MJ */
