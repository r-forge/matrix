#include "Mutils.h"

/* Slot validity methods ===============================================
   Called by various class validity methods (see below).
*/

/**
 * Test that `dim` is a length-2, non-negative integer vector.
 *
 * @param dim A `SEXP`, 
 *     typically the `Dim` slot of a (to be validated) `Matrix`.
 * @param domain A string specifying a domain for message translation.
 *
 * @return Either `TRUE` (indicating success) or a length-1 `STRSXP`
 *     containing an error message.
 */
SEXP Dim_validate(SEXP dim, const char* domain)
{
    /* TODO? coerce from REALSXP to INTSXP?
       // if (TYPEOF(dim) != INTSXP && TYPEOF(dim) != REALSXP)
       //     return mkString(_("'Dim' slot is not numeric"));
       though above is not enough as we must prohibit Dim[i] > INT_MAX

       FIXME? Prohibit is.object(dim) or maybe just inherits(dim, "factor") 
       and return a different error message in that case?
    */
    if (TYPEOF(dim) != INTSXP)
	return mkString(_("'Dim' slot is not of type \"integer\""));
    if (LENGTH(dim) != 2)
	return mkString(_("'Dim' slot does not have length 2"));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (m == NA_INTEGER || n == NA_INTEGER)
	return mkString(_("'Dim' slot contains NA"));
    if (m < 0 || n < 0)
	return mkString(dngettext(domain,
				  "'Dim' slot contains negative value",
				  "'Dim' slot contains negative values",
				  (m < 0 && n < 0) ? 2 : 1));
    return ScalarLogical(1);
}

SEXP R_Dim_validate(SEXP dim)
{
    return Dim_validate(dim, "Matrix");
}

#ifdef Matrix_SupportingCachedMethods

SEXP R_Dim_validate_old(SEXP obj, SEXP domain)
{
    return Dim_validate(GET_SLOT(obj, Matrix_DimSym),
			CHAR(STRING_ELT(domain, 0)));
}

#endif

/**
 * Test that `dimnames` is a valid length-2 list.
 *
 * @param dimnames A `SEXP`,
 *     typically the `Dimnames` slot of a (to be validated) `Matrix`.
 * @param pdim Pointer to a length-2, non-negative `int` array,
 *     typically from the `Dim` slot of a (to be validated) `Matrix`.
 *     Array validity _must_ be checked by the caller.
 *
 * @return Either `TRUE` (indicating success) or a length-1 `STRSXP`
 *     containing an error message.
 */
SEXP DimNames_validate(SEXP dimnames, int *pdim)
{
    char *buf;

    /* Allocate only when needed: in valid case, it is _not_ needed */
#define SPRINTF								\
    buf = Alloca(Matrix_ErrorBufferSize, char); R_CheckStack(); sprintf

    if (TYPEOF(dimnames) != VECSXP) {
	SPRINTF(buf, _("'Dimnames' slot is not a list"));
	return mkString(buf);
    }
    if (LENGTH(dimnames) != 2) {
	SPRINTF(buf, _("'Dimnames' slot does not have length 2"));
	return mkString(buf);
    }
    for (int j = 0; j < 2; ++j) {
	/* Behave as 'do_matrix()' from src/main/array.c:
	   Dimnames[[j]] must be NULL or _coercible to_ character
	   of length Dim[j] or 0 ... see 'R_Dimnames_fixup()' below
	*/
	SEXP s = VECTOR_ELT(dimnames, j);
	if (!isNull(s)) {
	    if (!isVector(s)) {
		SPRINTF(buf, _("Dimnames[[%d]] is not NULL or a vector"), j+1);
		return mkString(buf);
	    }
	    if (LENGTH(s) != pdim[j]) {
		if (LENGTH(s) != 0) {
		    SPRINTF(buf, _("length of Dimnames[[%d]] (%d) "
				   "is not equal to Dim[%d] (%d)"),
			    j+1, LENGTH(s), j+1, pdim[j]);
		    return mkString(buf);
		}
	    }
	}
    }
#undef SPRINTF    
    return ScalarLogical(1);
}

SEXP R_DimNames_validate(SEXP dimnames, SEXP dim)
{
    return DimNames_validate(dimnames, INTEGER(dim));
}

#ifdef Matrix_SupportingCachedMethods

SEXP R_DimNames_validate_old(SEXP obj)
{
    return DimNames_validate(GET_SLOT(obj, Matrix_DimNamesSym),
			     INTEGER(GET_SLOT(obj, Matrix_DimSym)));
}

#endif


/* Class validity methods ==============================================
   NB: These assume that validity methods for superclasses 
   have already been called via validObject() ...
*/

SEXP Matrix_validate(SEXP obj)
{
    SEXP dim = GET_SLOT(obj, Matrix_DimSym), val = Dim_validate(dim, "Matrix");
    return (isString(val))
	? val
	: DimNames_validate(GET_SLOT(obj, Matrix_DimNamesSym), INTEGER(dim));
}

SEXP MatrixFactorization_validate(SEXP obj)
{
    return Dim_validate(GET_SLOT(obj, Matrix_DimSym), "MatrixFactorization");
}

SEXP compMatrix_validate(SEXP obj)
{
    SEXP fac = GET_SLOT(obj, Matrix_factorSym);
    if (TYPEOF(fac) != VECSXP ||
	(XLENGTH(fac) > 0 && isNull(getAttrib(fac, R_NamesSymbol))))
	return mkString(_("'factors' slot is not a named list"));
    else
	return ScalarLogical(1);
}

#define TYPEMATRIX_VALIDATE(_PREFIX_, _SEXPTYPE_, _T2C_SEXPTYPE_)	\
SEXP _PREFIX_ ## Matrix_validate(SEXP obj)				\
{									\
    if (TYPEOF(GET_SLOT(obj, Matrix_xSym)) != _SEXPTYPE_)		\
	return mkString(_("'x' slot is not of type \"" #_T2C_SEXPTYPE_ "\"")); \
    else								\
	return ScalarLogical(1);					\
}
/* dMatrix_validate() */
TYPEMATRIX_VALIDATE(     d, REALSXP,  double)
/* lMatrix_validate() */
TYPEMATRIX_VALIDATE(     l,  LGLSXP, logical)
/* ndenseMatrix_validate() */
/* NB: "nsparseMatrix" has no 'x' slot, only "ndenseMatrix" ... */
TYPEMATRIX_VALIDATE(ndense,  LGLSXP, logical)
/* iMatrix_validate() */
TYPEMATRIX_VALIDATE(     i,  INTSXP, integer)
/* zMatrix_validate() */
TYPEMATRIX_VALIDATE(     z, CPLXSXP, complex)
#undef TYPEMATRIX_VALIDATE

SEXP symmetricMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    
#ifdef ENFORCE_SYMMETRIC_DIMNAMES
    /* This check can be expensive when both rownames and colnames have 
       nonzero length, and even more so when coercions to character are 
       required ... Users can avoid the expense by setting at least one 
       of rownames and colnames to NULL or by ensuring that they are the 
       same object, as testing for pointer equality is fast ...
     */
    
# define ANY_TO_STRING(x)					\
    (isString(x)						\
     ? x							\
     : (inherits(x, "factor")					\
	? asCharacterFactor(x)					\
	: coerceVector(x, STRSXP)))

    SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym),
	ndn = getAttrib(dn, R_NamesSymbol);
    const char *ndn0, *ndn1;
    if (!isNull(ndn) &&
	*(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
	*(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
	strcmp(ndn0, ndn1) != 0)
	return mkString(_("Dimnames[1] differs from Dimnames[2]"));
    if (n > 0) {
	/* NB: It is already known that the length of 'dn[[i]]' is 0 or 'n' */ 
	SEXP rn, cn;
	if (!isNull(rn = VECTOR_ELT(dn, 0)) &&
	    !isNull(cn = VECTOR_ELT(dn, 1)) &&
	    LENGTH(rn) == n &&
	    LENGTH(cn) == n &&
	    rn != cn &&
	    !equal_string_vectors(ANY_TO_STRING(rn), ANY_TO_STRING(cn), n))
	    return mkString(_("Dimnames[1] differs from Dimnames[2]"));
    }
    
# undef ANY_TO_STRING
#endif

    SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
    if (TYPEOF(uplo) != STRSXP)
	return mkString(_("'uplo' slot is not of type \"character\""));
    if (LENGTH(uplo) != 1)
	return mkString(_("'uplo' slot does not have length 1"));
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
	return mkString(_("'uplo' slot is not \"U\" or \"N\""));

    return ScalarLogical(1);
}

SEXP triangularMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (pdim[1] != pdim[0])
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));

    SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
    if (TYPEOF(uplo) != STRSXP)
	return mkString(_("'uplo' slot is not of type \"character\""));
    if (LENGTH(uplo) != 1)
	return mkString(_("'uplo' slot does not have length 1"));
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
	return mkString(_("'uplo' slot is not \"U\" or \"N\""));
    
    SEXP diag = GET_SLOT(obj, Matrix_diagSym);
    if (TYPEOF(diag) != STRSXP)
	return mkString(_("'diag' slot is not of type \"character\""));
    if (LENGTH(diag) != 1)
	return mkString(_("'diag' slot does not have length 1"));
    const char *di = CHAR(STRING_ELT(diag, 0));
    if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
	return mkString(_("'diag' slot is not \"N\" or \"U\""));
    
    return ScalarLogical(1);
}

SEXP diagonalMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    
    SEXP diag = GET_SLOT(obj, Matrix_diagSym);
    if (TYPEOF(diag) != STRSXP)
	return mkString(_("'diag' slot is not of type \"character\""));
    if (LENGTH(diag) != 1)
	return mkString(_("'diag' slot does not have length 1"));
    const char *di = CHAR(STRING_ELT(diag, 0));
    if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
	return mkString(_("'diag' slot is not \"N\" or \"U\""));

    int nx = LENGTH(GET_SLOT(obj, Matrix_xSym));
    if (di[0] == 'N') {
	if (nx != n)
	    return mkString(_("'diag' slot is \"N\" but "
			      "'x' slot does not have length n=Dim[1]"));
    } else {
	if (nx != 0)
	    return mkString(_("'diag' slot is \"U\" (identity matrix) but "
			      "'x' slot does not have length 0"));
    }
    return ScalarLogical(1);
}

SEXP indMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), m = pdim[0], n = pdim[1];
    if (n == 0 && m > 0)
	return mkString(_("m-by-0 indMatrix invalid for positive 'm'"));
    SEXP perm = GET_SLOT(obj, Matrix_permSym);
    if (TYPEOF(perm) != INTSXP)
	return mkString(_("'perm' slot is not of type \"integer\""));
    if (XLENGTH(perm) != m)
	return mkString(_("'perm' slot does not have length Dim[1]"));
    int i, *pperm = INTEGER(perm);
    for (i = 0; i < m; ++i) {
	if (*pperm == NA_INTEGER)
	    return mkString(_("'perm' slot contains NA"));
	if (*pperm < 1 || *pperm > n)
	    return mkString(_("'perm' slot has elements not in {1,...,Dim[2]}"));
	++pperm;
    }
    return ScalarLogical(1);
}

SEXP pMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    if (n > 1) {
	int i, *work, *pperm = INTEGER(GET_SLOT(obj, Matrix_permSym));
	Calloc_or_Alloca_TO(work, n, int);
	Memzero(work, n);
	--work;
	for (i = 0; i < n; ++i) {
	    if (work[*pperm])
		return mkString(_("'perm' slot contains duplicates"));
	    work[*pperm] = 1;
	    ++pperm;
	}
	++work;
	Free_FROM(work, n);
    }
    return ScalarLogical(1);
}

SEXP CsparseMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), m = pdim[0], n = pdim[1];

    SEXP p = GET_SLOT(obj, Matrix_pSym);
    if (TYPEOF(p) != INTSXP)
	return mkString(_("'p' slot is not of type \"integer\""));
    if (XLENGTH(p) - 1 != n)
	return mkString(_("'p' slot does not have length Dim[2]+1"));
    int *pp = INTEGER(p);
    
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    if (TYPEOF(i) != INTSXP)
	return mkString(_("'i' slot is not of type \"integer\""));
    if (pp[n] != NA_INTEGER && XLENGTH(i) < pp[n])
	return mkString(_("'i' slot has length less than p[length(p)]"));
    int *pi = INTEGER(i);

    if (pp[0] != 0)
	return mkString(_("first element of 'p' slot is not 0"));

    int k = 0, kend, ik, imin;
    while ((n--) > 0) {
	kend = *(++pp);
	if (kend == NA_INTEGER)
	    return mkString(_("'p' slot contains NA"));
	if (kend < k)
	    return mkString(_("'p' slot is not nondecreasing"));
	if (kend - k > m)
	    return mkString(_("first differences of 'p' slot exceed Dim[1]"));
	imin = -1;
	while (k < kend) {
	    ik = pi[k];
	    if (ik == NA_INTEGER)
		return mkString(_("'i' slot contains NA"));
	    if (ik < 0 || ik >= m)
		return mkString(_("'i' slot has elements not in {0,...,Dim[1]-1}"));
	    if (ik <= imin)
		return mkString(_("'i' slot is not increasing within columns"));
	    imin = ik;
	    ++k;
	}
    }
    return ScalarLogical(1);
}

SEXP RsparseMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), m = pdim[0], n = pdim[1];

    SEXP p = GET_SLOT(obj, Matrix_pSym);
    if (TYPEOF(p) != INTSXP)
	return mkString(_("'p' slot is not of type \"integer\""));
    if (XLENGTH(p) - 1 != m)
	return mkString(_("'p' slot does not have length Dim[1]+1"));
    int *pp = INTEGER(p);
    
    SEXP j = GET_SLOT(obj, Matrix_jSym);
    if (TYPEOF(j) != INTSXP)
	return mkString(_("'j' slot is not of type \"integer\""));
    if (pp[m] != NA_INTEGER && XLENGTH(j) < pp[m])
	return mkString(_("'j' slot has length less than p[length(p)]"));
    int *pj = INTEGER(j);

    if (pp[0] != 0)
	return mkString(_("p[1] is not 0"));

    int k = 0, kend, jk, jmin;
    while ((m--) > 0) {
	kend = *(++pp);
	if (kend == NA_INTEGER)
	    return mkString(_("'p' slot contains NA"));
	if (kend < k)
	    return mkString(_("'p' slot is not nondecreasing"));
	if (kend - k > n)
	    return mkString(_("first differences of 'p' slot exceed Dim[2]"));
	jmin = -1;
	while (k < kend) {
	    jk = pj[k];
	    if (jk == NA_INTEGER)
		return mkString(_("'j' slot contains NA"));
	    if (jk < 0 || jk >= n)
		return mkString(_("'j' slot has elements not in {0,...,Dim[2]-1}"));
	    if (jk <= jmin)
		return mkString(_("'j' slot is not increasing within rows"));
	    jmin = jk;
	    ++k;
	}
    }
    return ScalarLogical(1);
}

SEXP TsparseMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), m = pdim[0], n = pdim[1];
    
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    if (TYPEOF(i) != INTSXP)
	return mkString(_("'i' slot is not of type \"integer\""));
    SEXP j = GET_SLOT(obj, Matrix_jSym);
    if (TYPEOF(j) != INTSXP)
	return mkString(_("'j' slot is not of type \"integer\""));
    R_xlen_t nnz = XLENGTH(i);
    if (XLENGTH(j) != nnz)
	return mkString(_("'i' and 'j' slots do not have equal length"));
    if (nnz == 0)
	return ScalarLogical(1);
    if (m == 0 || n == 0)
	return mkString(_("'i' slot has nonzero length but prod(Dim) is 0"));
    int *pi = INTEGER(i), *pj = INTEGER(j);
    while ((nnz--) > 0) {
	if (*pi == NA_LOGICAL)
	    return mkString(_("'i' slot contains NA"));
	if (*pj == NA_LOGICAL)
	    return mkString(_("'j' slot contains NA"));
	if (*pi < 0 || *pi >= m)
	    return mkString(_("'i' slot has elements not in {0,...,Dim[1]-1}"));
	if (*pj < 0 || *pj >= n)
	    return mkString(_("'j' slot has elements not in {0,...,Dim[2]-1}"));
	++pi;
	++pj;
    }
    return ScalarLogical(1);
}

SEXP sCMatrix_validate(SEXP obj)
{
    SEXP p = GET_SLOT(obj, Matrix_pSym);
    int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
    if (pp[n] == 0)
	return ScalarLogical(1);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    int *pi = INTEGER(i), k = 0, kend, j;
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] > j)
		    return mkString(_("uplo=\"U\" but there are entries below the diagonal"));
		++k;
	    }
	}
    } else {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] < j)
		    return mkString(_("uplo=\"L\" but there are entries above the diagonal"));
		++k;
	    }
	}
    }
    return ScalarLogical(1);
}

SEXP tCMatrix_validate(SEXP obj)
{
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_diagSym), 0)) == 'N')
	return sCMatrix_validate(obj);
    SEXP p = GET_SLOT(obj, Matrix_pSym);
    int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
    if (pp[n] == 0)
	return ScalarLogical(1);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    int *pi = INTEGER(i), k = 0, kend, j;
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] >= j)
		    return mkString(_((pi[k] == j) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal"));
		++k;
	    }
	}
    } else {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] <= j)
		    return mkString(_((pi[k] == j) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal"));
		++k;
	    }
	}
    }
    return ScalarLogical(1);
}

SEXP sRMatrix_validate(SEXP obj)
{
    SEXP p = GET_SLOT(obj, Matrix_pSym);
    int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
    if (pp[m] == 0)
	return ScalarLogical(1);
    SEXP j = GET_SLOT(obj, Matrix_jSym);
    int *pj = INTEGER(j), k = 0, kend, i;
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] < i)
		    return mkString(_("uplo=\"U\" but there are entries below the diagonal"));
		++k;
	    }
	}
    } else {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] > i)
		    return mkString(_("uplo=\"L\" but there are entries above the diagonal"));
		++k;
	    }
	}
    }
    return ScalarLogical(1);
}

SEXP tRMatrix_validate(SEXP obj)
{
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_diagSym), 0)) == 'N')
	return sRMatrix_validate(obj);
    SEXP p = GET_SLOT(obj, Matrix_pSym);
    int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
    if (pp[m] == 0)
	return ScalarLogical(1);
    SEXP j = GET_SLOT(obj, Matrix_jSym);
    int *pj = INTEGER(j), k = 0, kend, i;
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] <= i)
		    return mkString(_((pj[k] == i) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal"));
		++k;
	    }
	}
    } else {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] >= i)
		    return mkString(_((pj[k] == i) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal"));
		++k;
	    }
	}
    }
    return ScalarLogical(1);
}

SEXP sTMatrix_validate(SEXP obj)
{
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    R_xlen_t nnz = XLENGTH(i);
    if (nnz == 0)
	return ScalarLogical(1);
    int *pi = INTEGER(i), *pj = INTEGER(GET_SLOT(obj, Matrix_jSym));
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	while ((nnz--) > 0)
	    if (*(pi++) > *(pj++))
		return mkString(_("uplo=\"U\" but there are entries below the diagonal"));
    } else {
	while ((nnz--) > 0)
	    if (*(pi++) < *(pj++))
		return mkString(_("uplo=\"L\" but there are entries above the diagonal"));
    }
    return ScalarLogical(1);
}

SEXP tTMatrix_validate(SEXP obj)
{
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_diagSym), 0)) == 'N')
	return sTMatrix_validate(obj);
    SEXP i = GET_SLOT(obj, Matrix_iSym);
    R_xlen_t nnz = XLENGTH(i);
    if (nnz == 0)
	return ScalarLogical(1);
    int *pi = INTEGER(i), *pj = INTEGER(GET_SLOT(obj, Matrix_jSym));
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	while ((nnz--) > 0) {
	    if (*pi >= *pj)
		return mkString(_((*pi == *pj) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal"));
	    ++pi;
	    ++pj;
	}
    } else {
	while ((nnz--) > 0) {
	    if (*pi <= *pj)
		return mkString(_((*pi == *pj) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal"));
	    ++pi;
	    ++pj;
	}
    }
    return ScalarLogical(1);
}

SEXP xgCMatrix_validate(SEXP obj)
{
    if (XLENGTH(GET_SLOT(obj, Matrix_xSym)) !=
	XLENGTH(GET_SLOT(obj, Matrix_iSym)))
	return mkString(_("'i' and 'x' slots do not have equal length"));
    else
	return ScalarLogical(1);
}

SEXP xsCMatrix_validate(SEXP obj)
{
    SEXP val = xgCMatrix_validate(obj);
    return (isString(val)) ? val : sCMatrix_validate(obj);
}

SEXP xtCMatrix_validate(SEXP obj)
{
    SEXP val = xgCMatrix_validate(obj);
    return (isString(val)) ? val : tCMatrix_validate(obj);
}

SEXP xgRMatrix_validate(SEXP obj)
{
    if (XLENGTH(GET_SLOT(obj, Matrix_xSym)) !=
	XLENGTH(GET_SLOT(obj, Matrix_jSym)))
	return mkString(_("'j' and 'x' slots do not have equal length"));
    else
	return ScalarLogical(1);
}

SEXP xsRMatrix_validate(SEXP obj)
{
    SEXP val = xgRMatrix_validate(obj);
    return (isString(val)) ? val : sRMatrix_validate(obj);
}

SEXP xtRMatrix_validate(SEXP obj)
{
    SEXP val = xgRMatrix_validate(obj);
    return (isString(val)) ? val : tRMatrix_validate(obj);
}

SEXP xgTMatrix_validate(SEXP obj)
{
    if (XLENGTH(GET_SLOT(obj, Matrix_xSym)) !=
	XLENGTH(GET_SLOT(obj, Matrix_iSym)))
	return mkString(_("'i' and 'x' slots do not have equal length"));
    else
	return ScalarLogical(1);
}

SEXP xsTMatrix_validate(SEXP obj)
{
    SEXP val = xgTMatrix_validate(obj);
    return (isString(val)) ? val : sTMatrix_validate(obj);
}

SEXP xtTMatrix_validate(SEXP obj)
{
    SEXP val = xgTMatrix_validate(obj);
    return (isString(val)) ? val : tTMatrix_validate(obj);
}

SEXP unpackedMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (XLENGTH(GET_SLOT(obj, Matrix_xSym)) != (R_xlen_t) pdim[0] * pdim[1])
	return mkString(_("'x' slot does not have length prod(Dim)"));
    else
	return ScalarLogical(1);
}

SEXP packedMatrix_validate(SEXP obj)
{
    R_xlen_t n = (R_xlen_t) INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    if (2 * XLENGTH(GET_SLOT(obj, Matrix_xSym)) != n * (n + 1))
        return mkString(_("'x' slot does not have length n*(n+1)/2, n=Dim[1]"));
    else
	return ScalarLogical(1);
}

SEXP dpoMatrix_validate(SEXP obj)
{
    double *x = REAL(GET_SLOT(obj, Matrix_xSym));
    int i, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t pos = 0, np1 = (R_xlen_t) n + 1;
    
    /* Non-negative diagonal elements are necessary _but not_ sufficient */
    for (i = 0; i < n; ++i, pos += np1)
	if (!ISNAN(x[pos]) && x[pos] < 0)
	    return mkString(_("matrix is not positive semidefinite"));
    return ScalarLogical(1);
}

SEXP dppMatrix_validate(SEXP obj)
{
    double *x = REAL(GET_SLOT(obj, Matrix_xSym));
    int i, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t pos = 0;

    /* Non-negative diagonal elements are necessary _but not_ sufficient */
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	for (i = 0; i < n; pos += (++i)+1)
	    if (x[pos] < 0)
		return mkString(_("matrix is not positive semidefinite"));
    } else {
	for (i = 0; i < n; pos += n-(i++))
	    if (x[pos] < 0)
		return mkString(_("matrix is not positive semidefinite"));
    }
    return ScalarLogical(1);
}

SEXP corMatrix_validate(SEXP obj)
{
    int i, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    SEXP sd = GET_SLOT(obj, Matrix_sdSym);
    if (XLENGTH(sd) != n)
	return mkString(_("'sd' slot does not have length n=Dim[1]"));
    double *psd = REAL(sd);
    for (i = 0; i < n; ++i) {
	if (!R_FINITE(psd[i]))
	    return mkString(_("'sd' slot has nonfinite elements"));
	if (psd[i] < 0)
	    return mkString(_("'sd' slot has negative elements"));
    }
    return ScalarLogical(1);
}

SEXP Cholesky_validate(SEXP obj)
{
    double *x = REAL(GET_SLOT(obj, Matrix_xSym));
    int i, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t pos = 0, np1 = (R_xlen_t) n + 1;
    
    /* Non-negative diagonal elements are necessary _and_ sufficient */
    for (i = 0; i < n; ++i, pos += np1)
	if (!ISNAN(x[pos]) && x[pos] < 0)
	    return mkString(_("matrix has negative diagonal elements"));
    return ScalarLogical(1);
}

SEXP pCholesky_validate(SEXP obj)
{
    double *x = REAL(GET_SLOT(obj, Matrix_xSym));
    int i, n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    R_xlen_t pos = 0;

    /* Non-negative diagonal elements are necessary _and_ sufficient */
    if (*CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0)) == 'U') {
	for (i = 0; i < n; pos += (++i)+1)
	    if (x[pos] < 0)
		return mkString(_("matrix has negative diagonal elements"));
    } else {
	for (i = 0; i < n; pos += n-(i++))
	    if (x[pos] < 0)
		return mkString(_("matrix has negative diagonal elements"));
    }
    return ScalarLogical(1);
}

SEXP BunchKaufman_validate(SEXP obj)
{
    SEXP perm = GET_SLOT(obj, Matrix_permSym);
    if (TYPEOF(perm) != INTSXP)
	return mkString(_("'perm' slot is not of type \"integer\""));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    if (LENGTH(perm) != n)
	return mkString(_("'perm' slot does not have length n=Dim[1]"));
    int n_ = n, *pperm = INTEGER(perm);
    while (n_ > 0) {
	if (*pperm == NA_INTEGER)
	    return mkString(_("'perm' slot contains NA"));
	if (*pperm < -n || *pperm == 0 || *pperm > n)
	    return mkString(_("'perm' slot has elements not in {-n,...,n}\\{0}, n=Dim[1]"));
	if (*pperm > 0) {
	    pperm += 1;
	    n_ -= 1;
	} else if (n_ > 1 && *(pperm + 1) == *pperm) {
	    pperm += 2;
	    n_ -= 2;
	} else {
	    return mkString(_("'perm' slot has an unpaired negative element"));
	}
    }
    return ScalarLogical(1);
}

SEXP pBunchKaufman_validate(SEXP obj)
{
    return BunchKaufman_validate(obj); /* since we only look at 'perm' */
}

SEXP Schur_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    SEXP Q = GET_SLOT(obj, Matrix_QSym), T = GET_SLOT(obj, Matrix_TSym);
    if (!(IS_S4_OBJECT(Q) &&
	  R_has_slot(Q, Matrix_DimSym) &&
	  R_has_slot(Q, Matrix_DimNamesSym) &&
	  !isString(Matrix_validate(Q))))
	return mkString(_("'Q' slot is not a valid Matrix"));
    if (!(IS_S4_OBJECT(T) &&
	  R_has_slot(T, Matrix_DimSym) &&
	  R_has_slot(T, Matrix_DimNamesSym) &&
	  !isString(Matrix_validate(T))))
	return mkString(_("'T' slot is not a valid Matrix"));
    pdim = INTEGER(GET_SLOT(Q, Matrix_DimSym));
    if (pdim[0] != n || pdim[1] != n)
	return mkString(_("dimensions of 'Q' are not identical to 'Dim'"));
    pdim = INTEGER(GET_SLOT(T, Matrix_DimSym));
    if (pdim[0] != n || pdim[1] != n)
	return mkString(_("dimensions of 'T' are not identical to 'Dim'"));
    SEXP e = GET_SLOT(obj, install("EValues"));
    SEXPTYPE te = TYPEOF(e);
    if (te != REALSXP && te != CPLXSXP)
	return mkString(_("'EValues' slot does not have type \"double\" "
			  "or type \"complex\""));
    if (LENGTH(e) != n)
	return mkString(_("'EValues' slot does not have length n=Dim[1]"));
    return ScalarLogical(1);
}

SEXP denseLU_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    SEXP val = DimNames_validate(GET_SLOT(obj, Matrix_DimNamesSym), pdim);
    if (isString(val))
	return val;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (TYPEOF(x) != REALSXP)
	return mkString(_("'x' slot is not of type \"double\""));
    int m = pdim[0], n = pdim[1];
    if (XLENGTH(x) != (double) m * n)
	return mkString(_("'x' slot does not have length prod(Dim)"));
    SEXP perm = GET_SLOT(obj, Matrix_permSym);
    if (TYPEOF(perm) != INTSXP)
	return mkString(_("'perm' slot is not of type \"integer\""));
    int r = (m < n) ? m : n;
    if (XLENGTH(perm) != r)
	return mkString(_("'perm' slot does not have length min(Dim)"));
    int *pperm = INTEGER(perm);
    while ((r--) > 0) {
	if (*pperm == NA_INTEGER)
	    return mkString(_("'perm' slot contains NA"));
	if (*pperm < 1 || *pperm > m)
	    return mkString(_("'perm' slot has elements not in {1,...,Dim[1]}"));
    }
    return ScalarLogical(1);
}

SEXP sparseLU_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    static const char *valid[] = { "dtCMatrix", "" };
    static SEXP validObject = NULL;
    if (!validObject)
	validObject = install("validObject");
    SEXP L = GET_SLOT(obj, Matrix_LSym), U = GET_SLOT(obj, Matrix_USym),
	call = PROTECT(lang3(validObject,
			     /* object = */ L, /* test = */ ScalarLogical(1)));
    if (R_check_class_etc(L, valid) != 0 || isString(eval(call, R_GlobalEnv))) {
	UNPROTECT(1);
	return mkString(_("'L' slot is not a valid dtCMatrix"));
    }
    SETCADR(call, U);
    if (R_check_class_etc(U, valid) != 0 || isString(eval(call, R_GlobalEnv))) {
	UNPROTECT(1);
	return mkString(_("'U' slot is not a valid dtCMatrix"));
    }
    UNPROTECT(1);
    pdim = INTEGER(GET_SLOT(L, Matrix_DimSym));
    if (pdim[0] != n || pdim[1] != n)
	return mkString(_("dimensions of 'L' slot are not identical to 'Dim'"));
    pdim = INTEGER(GET_SLOT(U, Matrix_DimSym));
    if (pdim[0] != n || pdim[1] != n)
	return mkString(_("dimensions of 'U' slot are not identical to 'Dim'"));
    if (*CHAR(STRING_ELT(GET_SLOT(L, Matrix_uploSym), 0)) == 'U')
	return mkString(_("'L' slot is upper (not lower) triangular"));
    if (*CHAR(STRING_ELT(GET_SLOT(U, Matrix_uploSym), 0)) != 'U')
	return mkString(_("'U' slot is lower (not upper) triangular"));
    SEXP p = GET_SLOT(obj, Matrix_pSym), q = GET_SLOT(obj, Matrix_qSym);
    if (TYPEOF(p) != INTSXP)
	return mkString(_("'p' slot is not of type \"integer\""));
    if (TYPEOF(q) != INTSXP)
	return mkString(_("'q' slot is not of type \"integer\""));
    if (LENGTH(p) != n)
	return mkString(_("'p' slot does not have length Dim[1]"));
    if (LENGTH(q) != n)
	return mkString(_("'q' slot does not have length Dim[1]"));
    int i, *pp = INTEGER(p), *pq = INTEGER(q);
    for (i = 0; i < n; ++i) {
	if (pp[i] == NA_INTEGER)
	    return mkString(_("'p' slot contains NA"));
	if (pq[i] == NA_INTEGER)
	    return mkString(_("'q' slot contains NA"));
	if (pp[i] < 0 || pp[i] >= n)
	    return mkString(_("'p' slot has elements not in {0,...,Dim[1]-1}"));
	if (pq[i] < 0 || pq[i] >= n)
	    return mkString(_("'q' slot has elements not in {0,...,Dim[1]-1}"));
    }
    int *work;
    Calloc_or_Alloca_TO(work, n, int);
    Memzero(work, n);
    for (i = 0; i < n; ++i) {
	if (work[*pp])
	    return mkString(_("'p' slot contains duplicates"));
	work[*pp] = 1;
	++pp;
    }
    Memzero(work, n);
    for (i = 0; i < n; ++i) {
	if (work[*pq])
	    return mkString(_("'q' slot contains duplicates"));
	work[*pq] = 1;
	++pq;
    }
    Free_FROM(work, n);
    return ScalarLogical(1);
}

SEXP sparseQR_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), m = pdim[0], n = pdim[1];
    if (m < n)
	return mkString(_("matrix has more columns than rows"));
    SEXP beta = GET_SLOT(obj, Matrix_betaSym);
    if (TYPEOF(beta) != REALSXP)
	return mkString(_("'beta' slot is not of type \"double\""));
    if (LENGTH(beta) != n)
	return mkString(_("'beta' slot does not have length Dim[2]"));
    static const char *valid[] = { "dgCMatrix", "" };
    static SEXP validObject = NULL;
    if (!validObject)
	validObject = install("validObject");
    SEXP V = GET_SLOT(obj, Matrix_VSym), R = GET_SLOT(obj, Matrix_RSym),
	call = PROTECT(lang3(validObject,
			     /* object = */ V, /* test = */ ScalarLogical(1)));
    if (R_check_class_etc(V, valid) != 0 || isString(eval(call, R_GlobalEnv))) {
	UNPROTECT(1);
	return mkString(_("'V' slot is not a valid dgCMatrix"));
    }
    SETCADR(call, R);
    if (R_check_class_etc(R, valid) != 0 || isString(eval(call, R_GlobalEnv))) {
	UNPROTECT(1);
	return mkString(_("'R' slot is not a valid dgCMatrix"));
    }
    UNPROTECT(1);
    int m2;
    pdim = INTEGER(GET_SLOT(V, Matrix_DimSym));
    if ((m2 = pdim[0]) < m)
	return mkString(_("'V' slot has fewer than Dim[1] rows"));
    if (pdim[1] != n)
	return mkString(_("'V' slot does not have Dim[2] columns"));
    pdim = INTEGER(GET_SLOT(R, Matrix_DimSym));
    if (pdim[0] != m2)
	return mkString(_("'R' slot does not have nrow(V) rows"));
    if (pdim[1] != n)
	return mkString(_("'R' slot does not have Dim[2] columns"));
    SEXP R_p = GET_SLOT(R, Matrix_pSym), R_i = GET_SLOT(R, Matrix_iSym);
    int *R_pp = INTEGER(R_p), *R_pi = INTEGER(R_i), k = 0, kend, j;
    for (j = 0; j < n; ++j) {
	kend = *(++R_pp);
	while (k < kend) {
	    if (R_pi[k] > j)
		return mkString(_("'R' slot has entries below the diagonal"));
	    ++k;
	}
    }
    SEXP p = GET_SLOT(obj, Matrix_pSym), q = GET_SLOT(obj, Matrix_qSym);
    if (TYPEOF(p) != INTSXP)
	return mkString(_("'p' slot is not of type \"integer\""));
    if (TYPEOF(q) != INTSXP)
	return mkString(_("'q' slot is not of type \"integer\""));
    if (LENGTH(p) != m2)
	return mkString(_("'p' slot does not have length nrow(V)"));
    if (LENGTH(q) != 0 && LENGTH(q) != n)
	return mkString(_("'q' slot does not have length Dim[2] or length 0"));
    int i, *pp = INTEGER(p), *pq = INTEGER(q);
    for (i = 0; i < m2; ++i) {
	if (pp[i] == NA_INTEGER)
	    return mkString(_("'p' slot contains NA"));
	if (pp[i] < 0 || pp[i] >= m2)
	    return mkString(_("'p' slot has elements not in {0,...,nrow(V)-1}"));
    }
    int *work;
    Calloc_or_Alloca_TO(work, m2, int); /* n <= m <= m2 */
    Memzero(work, m2);
    for (i = 0; i < m2; ++i) {
	if (work[*pp])
	    return mkString(_("'p' slot contains duplicates"));
	work[*pp] = 1;
	++pp;
    }
    if (LENGTH(q) == n) {
	for (i = 0; i < n; ++i) {
	    if (pq[i] == NA_INTEGER)
		return mkString(_("'q' slot contains NA"));
	    if (pq[i] < 0 || pq[i] >= n)
		return mkString(_("'q' slot has elements not in {0,...,Dim[2]-1}"));
	}
	Memzero(work, n);
	for (i = 0; i < n; ++i) {
	    if (work[*pq])
		return mkString(_("'q' slot contains duplicates"));
	    work[*pq] = 1;
	    ++pq;
	}
    }
    Free_FROM(work, m2);
    return ScalarLogical(1);
}

/* NB: below three should use tests in ./chm_common.c ... C-s validate */

SEXP CHMfactor_validate(SEXP obj) /* TODO */
{
    return ScalarLogical(1);
}

SEXP CHMsimpl_validate(SEXP obj) /* TODO */
{
    return ScalarLogical(1);
}

SEXP CHMsuper_validate(SEXP obj) /* TODO */
{
    return ScalarLogical(1);
}


