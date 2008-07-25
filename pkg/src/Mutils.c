#include <limits.h>

#include "Mutils.h"
#include <R_ext/Lapack.h>

#include <Rversion.h>

#if 0 /* defined(R_VERSION) && R_VERSION >= R_Version(2, 7, 0) *
       * La_norm_type() & La_rcond_type() are now in R_ext/Lapack.h
       * but because of the 'module-mess' that's not sufficient */
#else
char La_norm_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(
	    _("argument type[1]='%s' must be a character string of string length 1"),
	    typstr);
    typup = toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* aliases */
    else if (typup == 'E')
	typup = 'F';
    else if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error(_("argument type[1]='%s' must be one of 'M','1','O','I','F' or 'E'"),
	      typstr);
    return typup;
}

char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error(_("argument type[1]='%s' must be one of '1','O', or 'I'"),
	      typstr);
    return typup;
}
#endif

double get_double_by_name(SEXP obj, char *nm)
{
    SEXP nms = getAttrib(obj, R_NamesSymbol);
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error(_("object must be a named, numeric vector"));
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    return REAL(obj)[i];
	}
    }
    return R_NaReal;
}

SEXP
set_double_by_name(SEXP obj, double val, char *nm)
{
    SEXP nms = getAttrib(obj, R_NamesSymbol);
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error("object must be a named, numeric vector");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    REAL(obj)[i] = val;
	    return obj;
	}
    }
    {
	SEXP nx = PROTECT(allocVector(REALSXP, len + 1)),
	    nnms = allocVector(STRSXP, len + 1);

	setAttrib(nx, R_NamesSymbol, nnms);
	for (i = 0; i < len; i++) {
	    REAL(nx)[i] = REAL(obj)[i];
	    SET_STRING_ELT(nnms, i, duplicate(STRING_ELT(nms, i)));
	}
	REAL(nx)[len] = val;
	SET_STRING_ELT(nnms, len, mkChar(nm));
	UNPROTECT(1);
	return nx;
    }
}

SEXP as_det_obj(double val, int log, int sign)
{
    SEXP det = PROTECT(allocVector(VECSXP, 2)),
	nms = allocVector(STRSXP, 2),
	vv = ScalarReal(val);

    setAttrib(det, R_NamesSymbol, nms);
    SET_STRING_ELT(nms, 0, mkChar("modulus"));
    SET_STRING_ELT(nms, 1, mkChar("sign"));
    setAttrib(vv, install("logarithm"), ScalarLogical(log));
    SET_VECTOR_ELT(det, 0, vv);
    SET_VECTOR_ELT(det, 1, ScalarInteger(sign));
    setAttrib(det, R_ClassSymbol, mkString("det"));
    UNPROTECT(1);
    return det;
}

SEXP get_factors(SEXP obj, char *nm)
{
    SEXP fac = GET_SLOT(obj, Matrix_factorSym),
	nms = getAttrib(fac, R_NamesSymbol);
    int i, len = length(fac);

    if ((!isNewList(fac)) || (length(fac) > 0 && nms == R_NilValue))
	error("factors slot must be a named list");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    return VECTOR_ELT(fac, i);
	}
    }
    return R_NilValue;
}

SEXP set_factors(SEXP obj, SEXP val, char *nm)
{
    SEXP fac = GET_SLOT(obj, Matrix_factorSym),
	nms = getAttrib(fac, R_NamesSymbol), nfac, nnms;
    int i, len = length(fac);

    if ((!isNewList(fac)) || (length(fac) > 0 && nms == R_NilValue))
	error("factors slot must be a named list");
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    SET_VECTOR_ELT(fac, i, duplicate(val));
	    return val;
	}
    }
    nfac = PROTECT(allocVector(VECSXP, len + 1));
    nnms = PROTECT(allocVector(STRSXP, len + 1));
    setAttrib(nfac, R_NamesSymbol, nnms);
    for (i = 0; i < len; i++) {
	SET_VECTOR_ELT(nfac, i, VECTOR_ELT(fac, i));
	SET_STRING_ELT(nnms, i, duplicate(STRING_ELT(nms, i)));
    }
    SET_VECTOR_ELT(nfac, len, duplicate(val));
    SET_STRING_ELT(nnms, len, mkChar(nm));
    SET_SLOT(obj, Matrix_factorSym, nfac);
    UNPROTECT(2);
    return val;
}

#if 0 				/* unused */
/* useful for all the ..CMatrix classes (and ..R by [0] <-> [1]); but unused */
SEXP dgCMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}
#endif	/* unused */

/* Fill in the "trivial remainder" in  n*m  array ;
 *  typically the 'x' slot of a "dtrMatrix" :
 * But should be usable for double/logical/int/complex : */

#define MAKE_TRIANGULAR_BODY(_TO_, _FROM_, _ZERO_, _ONE_)	\
{								\
    int i, j, *dims = INTEGER(GET_SLOT(_FROM_, Matrix_DimSym));	\
    int n = dims[0], m = dims[1];				\
								\
    if (*uplo_P(_FROM_) == 'U') {				\
	for (j = 0; j < n; j++)					\
	    for (i = j+1; i < m; i++)				\
		_TO_[i + j*m] = _ZERO_;				\
    } else {							\
	for (j = 1; j < n; j++)					\
	    for (i = 0; i < j && i < m; i++)			\
		_TO_[i + j*m] = _ZERO_;				\
    }								\
    if (*diag_P(_FROM_) == 'U') {				\
	j = (n < m) ? n : m;					\
	for (i = 0; i < j; i++)					\
	    _TO_[i * (m + 1)] = _ONE_;				\
    }								\
}

void make_d_matrix_triangular(double *to, SEXP from)
    MAKE_TRIANGULAR_BODY(to, from, 0., 1.)
void make_i_matrix_triangular(int *to, SEXP from)
    MAKE_TRIANGULAR_BODY(to, from, 0, 1)


/* Should work for double/logical/int/complex : */
#define MAKE_SYMMETRIC_BODY(_TO_, _FROM_)			\
{								\
    int i, j, n = INTEGER(GET_SLOT(_FROM_, Matrix_DimSym))[0];	\
								\
    if (*uplo_P(_FROM_) == 'U') {				\
	for (j = 0; j < n; j++)					\
	    for (i = j+1; i < n; i++)				\
		_TO_[i + j*n] = _TO_[j + i*n];			\
    } else {							\
	for (j = 1; j < n; j++)					\
	    for (i = 0; i < j && i < n; i++)			\
		_TO_[i + j*n] = _TO_[j + i*n];			\
    }								\
}

void make_d_matrix_symmetric(double *to, SEXP from)
    MAKE_SYMMETRIC_BODY(to, from)

void make_i_matrix_symmetric(int *to, SEXP from)
    MAKE_SYMMETRIC_BODY(to, from)


/**
 * Create a named vector of type TYP
 *
 * @param TYP a vector SEXP type (e.g. REALSXP)
 * @param names names of list elements with null string appended
 *
 * @return pointer to a named vector of type TYP
 */
SEXP
Matrix_make_named(int TYP, char **names)
{
    SEXP ans, nms;
    int i, n;

    for (n = 0; strlen(names[n]) > 0; n++) {}
    ans = PROTECT(allocVector(TYP, n));
    nms = PROTECT(allocVector(STRSXP, n));
    for (i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

#define Matrix_Error_Bufsiz    4096

SEXP check_scalar_string(SEXP sP, char *vals, char *nm)
{
    SEXP val = ScalarLogical(1);
    char *buf;
    /* only allocate when needed: in good case, none is needed */
#define SPRINTF buf = Alloca(Matrix_Error_Bufsiz, char); R_CheckStack(); sprintf

    if (length(sP) != 1) {
	SPRINTF(buf, _("'%s' slot must have length 1"), nm);
    } else {
	const char *str = CHAR(STRING_ELT(sP, 0));
	if (strlen(str) != 1) {
	    SPRINTF(buf, _("'%s' must have string length 1"), nm);
	} else {
	    int i, len, match;
	    for (i = 0, len = strlen(vals), match = 0; i < len; i++) {
		if (str[0] == vals[i])
		    return R_NilValue;
	    }
	    SPRINTF(buf, _("'%s' must be in '%s'"), nm, vals);
	}
    }
    /* 'error' returns : */
    val = mkString(buf);
    return val;
#undef SPRINTF
}

/* FIXME? Something like this should be part of the R API ?
 *        But then, R has the more general  compute_identical()
 * in src/main/identical.c: Rboolean compute_identical(SEXP x, SEXP y);
*/
Rboolean equal_string_vectors(SEXP s1, SEXP s2)
{
    Rboolean n1 = isNull(s1), n2 = isNull(s2);
    if (n1 || n2)
	return (n1 == n2) ? TRUE : FALSE;
    else if (TYPEOF(s1) != STRSXP || TYPEOF(s2) != STRSXP) {
	error(_("'s1' and 's2' must be \"character\" vectors"));
	return FALSE; /* -Wall */
    } else {
	int n = LENGTH(s1), i;
	if (n != LENGTH(s2))
	    return FALSE;
	for(i = 0; i < n; i++) {
	    /* note that compute_identical() code for STRSXP
	       is careful about NA's which we don't need */
	    if(strcmp(CHAR(STRING_ELT(s1, i)),
		      CHAR(STRING_ELT(s2, i))))
		return FALSE;
	}
	return TRUE; /* they *are* equal */
    }
}


SEXP dense_nonpacked_validate(SEXP obj)
{
    int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if ((dims[0] * dims[1]) != length(GET_SLOT(obj, Matrix_xSym)))
	return mkString(_("length of x slot != prod(Dim)"));
    return ScalarLogical(1);
}


#define PACKED_TO_FULL(TYPE)						\
TYPE *packed_to_full_ ## TYPE(TYPE *dest, const TYPE *src,		\
		        int n, enum CBLAS_UPLO uplo)			\
{									\
    int i, j, pos = 0;							\
									\
    AZERO(dest, n*n);							\
    for (j = 0; j < n; j++) {						\
	switch(uplo) {							\
	case UPP:							\
	    for (i = 0; i <= j; i++) dest[i + j * n] = src[pos++];	\
	    break;							\
	case LOW:							\
	    for (i = j; i < n; i++) dest[i + j * n] = src[pos++];	\
	    break;							\
	default:							\
	    error(_("'uplo' must be UPP or LOW"));			\
	}								\
    }									\
    return dest;							\
}

PACKED_TO_FULL(double)
PACKED_TO_FULL(int)

#define FULL_TO_PACKED(TYPE)						\
TYPE *full_to_packed_ ## TYPE(TYPE *dest, const TYPE *src, int n,	\
		      enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag)	\
{									\
    int i, j, pos = 0;							\
									\
    for (j = 0; j < n; j++) {						\
	switch(uplo) {							\
	case UPP:							\
	    for (i = 0; i <= j; i++)					\
		dest[pos++] = (i == j && diag== UNT) ? 1 : src[i + j*n]; \
	    break;							\
	case LOW:							\
	    for (i = j; i < n; i++)					\
		dest[pos++] = (i == j && diag== UNT) ? 1 : src[i + j*n]; \
	    break;							\
	default:							\
	    error(_("'uplo' must be UPP or LOW"));			\
	}								\
    }									\
    return dest;							\
}

FULL_TO_PACKED(double)
FULL_TO_PACKED(int)



/**
 * Copy the diagonal elements of the packed denseMatrix x to dest
 *
 * @param dest vector of length ncol(x)
 * @param x pointer to an object representing a packed array
 * @param n number of columns in the matrix represented by x
 *
 * @return dest
 */
void d_packed_getDiag(double *dest, SEXP x, int n)
{
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));

#define END_packed_getDiag						\
    int j, pos = 0;							\
									\
    if (*uplo_P(x) == 'U') {						\
	for(pos= 0, j=0; j < n; pos += 1+(++j))	 dest[j] = xx[pos];	\
    } else {								\
	for(pos= 0, j=0; j < n; pos += (n - j), j++) dest[j] = xx[pos]; \
    }									\
    return

    END_packed_getDiag;
}

void l_packed_getDiag(int *dest, SEXP x, int n)
{
    int *xx = LOGICAL(GET_SLOT(x, Matrix_xSym));

    END_packed_getDiag;
}

#undef END_packed_getDiag

void tr_d_packed_getDiag(double *dest, SEXP x)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocVector(REALSXP, n));
    double *v = REAL(val);

    if (*diag_P(x) == 'U') {
	int j;
	for (j = 0; j < n; j++) v[j] = 1.;
    } else {
	d_packed_getDiag(v, x, n);
    }
    UNPROTECT(1);
    return;
}

void tr_l_packed_getDiag(   int *dest, SEXP x)
{
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP val = PROTECT(allocVector(LGLSXP, n));
    int *v = LOGICAL(val);

    if (*diag_P(x) == 'U') {
	int j;
	for (j = 0; j < n; j++) v[j] = 1;
    } else {
	l_packed_getDiag(v, x, n);
    }
    UNPROTECT(1);
    return;
}

SEXP Matrix_expand_pointers(SEXP pP)
{
    int n = length(pP) - 1;
    int *p = INTEGER(pP);
    SEXP ans = PROTECT(allocVector(INTSXP, p[n]));

    expand_cmprPt(n, p, INTEGER(ans));
    UNPROTECT(1);
    return ans;
}


/**
 * Return the element of a given name from a named list
 *
 * @param list
 * @param nm name of desired element
 *
 * @return element of list with name nm
 */
SEXP
Matrix_getElement(SEXP list, char *nm) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < LENGTH(names); i++)
	if (!strcmp(CHAR(STRING_ELT(names, i)), nm))
	    return(VECTOR_ELT(list, i));
    return R_NilValue;
}

/**
 * Zero a square matrix of size nc then copy a vector to the diagonal
 *
 * @param dest destination array of length nc * nc
 * @param A pointer to a square Matrix object
 *
 * @return dest
 */
static double *
install_diagonal(double *dest, SEXP A)
{
    int nc = INTEGER(GET_SLOT(A, Matrix_DimSym))[0];
    int i, ncp1 = nc + 1, unit = *diag_P(A) == 'U';
    double *ax = REAL(GET_SLOT(A, Matrix_xSym));

    AZERO(dest, nc * nc);
    for (i = 0; i < nc; i++)
	dest[i * ncp1] = (unit) ? 1. : ax[i];
    return dest;
}

static int *
install_diagonal_int(int *dest, SEXP A)
{
    int nc = INTEGER(GET_SLOT(A, Matrix_DimSym))[0];
    int i, ncp1 = nc + 1, unit = *diag_P(A) == 'U';
    int *ax = INTEGER(GET_SLOT(A, Matrix_xSym));

    AZERO(dest, nc * nc);
    for (i = 0; i < nc; i++)
	dest[i * ncp1] = (unit) ? 1 : ax[i];
    return dest;
}


/**  Duplicate a [dln]denseMatrix or a numeric matrix or vector
 *  as a [dln]geMatrix.
 *  This is for the many *_matrix_{prod,crossprod,tcrossprod,etc.}
 *  functions that work with both classed and unclassed matrices.
 *
 * @param A	  either a denseMatrix object or a matrix object
 */
/* NOTA BENE: If you enlarge this list, do change '14' and '6' below !
 * ---------  ddiMatrix & ldiMatrix  are no longer ddense or ldense on the R level,
 *            ---         ---        but are still dealt with here.
 */

#define ddense_CLASSES							\
		    "dgeMatrix", "dtrMatrix",				\
		    "dsyMatrix", "dpoMatrix", "ddiMatrix",		\
		    "dtpMatrix", "dspMatrix", "dppMatrix",		\
		    /* sub classes of those above:*/			\
		    /* dtr */ "Cholesky", "LDL", "BunchKaufman",	\
		    /* dtp */ "pCholesky", "pBunchKaufman",		\
		    /* dpo */ "corMatrix"

#define ldense_CLASSES					\
		    "lgeMatrix", "ltrMatrix",		\
		    "lsyMatrix", "ldiMatrix",		\
		    "ltpMatrix", "lspMatrix"

#define ndense_CLASSES					\
		    "ngeMatrix", "ntrMatrix",		\
		    "nsyMatrix",			\
		    "ntpMatrix", "nspMatrix"

/* Generalized -- "geMatrix" -- dispatch where needed : */
SEXP dup_mMatrix_as_geMatrix(SEXP A)
{
    SEXP ans, ad = R_NilValue, an = R_NilValue;	/* -Wall */
    const char *cl = class_P(A);
    char *valid[] = {"_NOT_A_CLASS_",/* *_CLASSES defined in ./Mutils.h */
		    ddense_CLASSES, /* 14 */
		    ldense_CLASSES, /* 6  */
		    ndense_CLASSES, /* 5  */
		    ""};
    int sz, ctype = Matrix_check_class(cl, valid), nprot = 1;
    enum dense_enum { ddense, ldense, ndense } M_type = ddense /* -Wall */;

    if (ctype > 0) { /* a [nld]denseMatrix object */
	ad = GET_SLOT(A, Matrix_DimSym);
	an = GET_SLOT(A, Matrix_DimNamesSym);
	M_type = (ctype <= 14) ? ddense :
	    ((ctype <= 14+6) ? ldense : ndense);
    }
    else if (ctype < 0) {	/* not a (recognized) classed matrix */

	if (isReal(A))
	    M_type = ddense;
	else if (isInteger(A)) {
	    A = PROTECT(coerceVector(A, REALSXP));
	    nprot++;
	    M_type = ddense;
	}
	else if (isLogical(A))
	    M_type = ldense;
	else
	    error(_("invalid class '%s' to dup_mMatrix_as_geMatrix"), cl);

#define	DUP_MMATRIX_NON_CLASS						\
	if (isMatrix(A)) {	/* "matrix" */				\
	    ad = getAttrib(A, R_DimSymbol);				\
	    an = getAttrib(A, R_DimNamesSymbol);			\
	} else {/* maybe "numeric" (incl integer,logical) --> (n x 1) */\
	    int* dd = INTEGER(ad = PROTECT(allocVector(INTSXP, 2)));	\
	    nprot++;							\
	    dd[0] = LENGTH(A);						\
	    dd[1] = 1;							\
	    an = R_NilValue;						\
	}								\
	ctype = 0

	DUP_MMATRIX_NON_CLASS;
    }

    ans = PROTECT(NEW_OBJECT(MAKE_CLASS(M_type == ddense ? "dgeMatrix" :
					(M_type == ldense ? "lgeMatrix" :
					 "ngeMatrix"))));
#define DUP_MMATRIX_SET_1					\
    SET_SLOT(ans, Matrix_DimSym, duplicate(ad));		\
    SET_SLOT(ans, Matrix_DimNamesSym, (LENGTH(an) == 2) ? 	\
	     duplicate(an): allocVector(VECSXP, 2));		\
    sz = INTEGER(ad)[0] * INTEGER(ad)[1]

    DUP_MMATRIX_SET_1;

    if(M_type == ddense) { /* ddense -> dge */

	double *ansx;

#define DUP_MMATRIX_ddense_CASES						\
	ansx = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, sz));			\
	switch(ctype) {								\
	case 0:			/* unclassed real matrix */			\
	    Memcpy(ansx, REAL(A), sz);						\
	    break;								\
	case 1:			/* dgeMatrix */					\
	    Memcpy(ansx, REAL(GET_SLOT(A, Matrix_xSym)), sz);			\
	    break;								\
	case 2:			/* dtrMatrix   and subclasses */		\
	case 9: case 10: case 11:   /* ---	Cholesky, LDL, BunchKaufman */	\
	    Memcpy(ansx, REAL(GET_SLOT(A, Matrix_xSym)), sz);			\
	    make_d_matrix_triangular(ansx, A);					\
	    break;								\
	case 3:			/* dsyMatrix */					\
	case 4:			/* dpoMatrix  + subclass */			\
	case 14:	 		/* ---	corMatrix */			\
	    Memcpy(ansx, REAL(GET_SLOT(A, Matrix_xSym)), sz);			\
	    make_d_matrix_symmetric(ansx, A);					\
	    break;								\
	case 5:			/* ddiMatrix */					\
	    install_diagonal(ansx, A);						\
	    break;								\
	case 6:			/* dtpMatrix  + subclasses */			\
	case 12: case 13: 		/* ---	pCholesky, pBunchKaufman */	\
	    packed_to_full_double(ansx, REAL(GET_SLOT(A, Matrix_xSym)),		\
				  INTEGER(ad)[0],				\
				  *uplo_P(A) == 'U' ? UPP : LOW);		\
	    make_d_matrix_triangular(ansx, A);					\
	    break;								\
	case 7:			/* dspMatrix */					\
	case 8:			/* dppMatrix */					\
	    packed_to_full_double(ansx, REAL(GET_SLOT(A, Matrix_xSym)),		\
				  INTEGER(ad)[0],				\
				  *uplo_P(A) == 'U' ? UPP : LOW);		\
	    make_d_matrix_symmetric(ansx, A);					\
	    break;								\
	}  /* switch(ctype) */

	DUP_MMATRIX_ddense_CASES

    }
    else { /* M_type == ldense || M_type = ndense  */
	/* ldense -> lge */
	/* ndense -> nge */
	int *ansx = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, sz));

	switch(ctype) {
	case 0:			/* unclassed logical matrix */
	    Memcpy(ansx, LOGICAL(A), sz);
	    break;

	case 1+14:			/* lgeMatrix */
	case 1+14+6:			/* ngeMatrix */
	    Memcpy(ansx, LOGICAL(GET_SLOT(A, Matrix_xSym)), sz);
	    break;
	case 2+14:			/* ltrMatrix */
	case 2+14+6:			/* ntrMatrix */
	    Memcpy(ansx, LOGICAL(GET_SLOT(A, Matrix_xSym)), sz);
	    make_i_matrix_triangular(ansx, A);
	    break;
	case 3+14:			/* lsyMatrix */
	case 3+14+6:			/* nsyMatrix */
	    Memcpy(ansx, LOGICAL(GET_SLOT(A, Matrix_xSym)), sz);
	    make_i_matrix_symmetric(ansx, A);
	    break;
	case 4+14:			/* ldiMatrix */
	case 4+14+6:			/* ndiMatrix */
	    install_diagonal_int(ansx, A);
	    break;
	case 5+14:			/* ltpMatrix */
	case 5+14+6:			/* ntpMatrix */
	    packed_to_full_int(ansx, LOGICAL(GET_SLOT(A, Matrix_xSym)),
			       INTEGER(ad)[0],
			       *uplo_P(A) == 'U' ? UPP : LOW);
	    make_i_matrix_triangular(ansx, A);
	    break;
	case 6+14:			/* lspMatrix */
	case 6+14+6:			/* nspMatrix */
	    packed_to_full_int(ansx, LOGICAL(GET_SLOT(A, Matrix_xSym)),
			       INTEGER(ad)[0],
			       *uplo_P(A) == 'U' ? UPP : LOW);
	    make_i_matrix_symmetric(ansx, A);
	    break;

	default:
	    error(_("unexpected ctype = %d in dup_mMatrix_as_geMatrix"), ctype);
	}  /* switch(ctype) */

    }  /* if(M_type == .) */

    UNPROTECT(nprot);
    return ans;
}

SEXP dup_mMatrix_as_dgeMatrix(SEXP A)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	ad = R_NilValue , an = R_NilValue;	/* -Wall */
    const char *cl = class_P(A);
    char *valid[] = {"_NOT_A_CLASS_", ddense_CLASSES, ""};
    int ctype = Matrix_check_class(cl, valid), nprot = 1, sz;
    double *ansx;

    if (ctype > 0) {		/* a ddenseMatrix object */
	ad = GET_SLOT(A, Matrix_DimSym);
	an = GET_SLOT(A, Matrix_DimNamesSym);
    }
    else if (ctype < 0) {	/* not a (recognized) classed matrix */

	DUP_MMATRIX_NON_CLASS;

	if (isInteger(A) || isLogical(A)) {
	    A = PROTECT(coerceVector(A, REALSXP));
	    nprot++;
	}
	if (!isReal(A))
	    error(_("invalid class '%s' to dup_mMatrix_as_dgeMatrix"), cl);
    }

    DUP_MMATRIX_SET_1;

    DUP_MMATRIX_ddense_CASES

    UNPROTECT(nprot);
    return ans;
}


SEXP new_dgeMatrix(int nrow, int ncol)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	 ad = PROTECT(allocVector(INTSXP, 2));

    INTEGER(ad)[0] = nrow;
    INTEGER(ad)[1] = ncol;
    SET_SLOT(ans, Matrix_DimSym, ad);
    SET_SLOT(ans, Matrix_DimNamesSym, allocVector(VECSXP, 2));
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, nrow * ncol);

    UNPROTECT(2);
    return ans;
}

/**
 * Return the 0-based index of a string match in a vector of strings
 * terminated by an empty string.  Returns -1 for no match.
 *
 * @param ij: 2-column integer matrix
 * @param di: dim(.), i.e. length 2 integer vector
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd(SEXP ij, SEXP di)
{
    SEXP ans;
    int *ij_di = INTEGER(getAttrib(ij, R_DimSymbol)), n = ij_di[0];
    int *Di = INTEGER(di), *IJ = INTEGER(ij);

    if(!(isMatrix(ij) && isInteger(ij) && ij_di[1] == 2))
	error(_("Argument ij must be 2-column integer matrix"));

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];
	int i, *j_ = IJ+n;
	for(i=0; i < n; i++) ii[i] = IJ[i] + j_[i] * nr;
    } else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];
	int i, *j_ = IJ+n;
	for(i=0; i < n; i++) ii[i] = IJ[i] + j_[i] * nr;
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Return the 0-based index of a string match in a vector of strings
 * terminated by an empty string.  Returns -1 for no match.
 *
 * @param i: integer vector
 * @param j: integer vector of same length as 'i'
 * @param di: dim(.), i.e. length 2 integer vector
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di)
{
    SEXP ans;
    int n = LENGTH(i);
    int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

    if(LENGTH(j) != n || !isInteger(i) || !isInteger(j))
	error(_("i and j must be integer vectors of the same length"));

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];
	int i;
	for(i=0; i < n; i++) ii[i] = i_[i] + j_[i] * nr;
    } else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];
	int i;
	for(i=0; i < n; i++) ii[i] = i_[i] + j_[i] * nr;
    }
    UNPROTECT(1);
    return ans;
}
