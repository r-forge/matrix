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
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	val = Dim_validate(dim, CHAR(STRING_ELT(domain, 0)));
    UNPROTECT(1); /* dim */
    return val;
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
    } else if (LENGTH(dimnames) != 2) {
	SPRINTF(buf, _("'Dimnames' slot does not have length 2"));
	return mkString(buf);
    }

    SEXP s;
    int i;
    R_xlen_t ns;

    /* Behave as do_matrix() from src/main/array.c:
       Dimnames[[i]] must be NULL or _coercible to_ character
       of length Dim[i] or 0 ... see R_Dimnames_fixup() below */
    for (i = 0; i < 2; ++i) {
	if (!isNull(s = VECTOR_ELT(dimnames, i))) {
	    if (!isVector(s)) {
		SPRINTF(buf, _("Dimnames[[%d]] is not NULL or a vector"), i+1);
		return mkString(buf);
	    } else if ((ns = XLENGTH(s)) != pdim[i] && ns != 0) {
		SPRINTF(buf, _("length of Dimnames[[%d]] (%lld) "
			       "is not equal to Dim[%d] (%d)"),
			i+1, (long long) ns, i+1, pdim[i]);
		return mkString(buf);
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
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	val = DimNames_validate(dimnames, INTEGER(dim));
    UNPROTECT(2); /* dimnames, dim */
    return val;
}

#endif

/**
 * @brief Sanitize user-supplied `[dD]imnames`.
 *
 * Replaces length-0 vectors with `NULL` and non-character vectors
 * with the result of coercing to character. Intended to emulate the
 * behaviour of `do_matrix()` from `src/main/array.c`.
 *
 * @param dn A list of length 2 passing `DimNames_validate()`.
 * 
 * @return A modified copy of `dn`, or `dn` if no modification is
 *    necessary.
 */
SEXP R_DimNames_fixup(SEXP dn)
{
    SEXP s;
    int i;
    Rboolean do_fixup = FALSE;
    for (i = 0; i < 2; ++i) {
	if (!isNull(s = VECTOR_ELT(dn, i)) &&
	    (LENGTH(s) == 0 || TYPEOF(s) != STRSXP)) {
	    do_fixup = TRUE;
	    break;
	}
    }
    if (do_fixup) {
	PROTECT(dn = duplicate(dn));
	for (i = 0; i < 2; ++i) {
	    if (isNull(s = VECTOR_ELT(dn, i)))
		continue;
	    if (LENGTH(s) == 0)
		SET_VECTOR_ELT(dn, i, R_NilValue);
	    else if (TYPEOF(s) != STRSXP) {
		if (inherits(s, "factor"))
		    PROTECT(s = asCharacterFactor(s));
		else {
		    PROTECT(s = coerceVector(s, STRSXP));
		    SET_ATTRIB(s, R_NilValue);
		    SET_OBJECT(s, 0);
		}
		SET_VECTOR_ELT(dn, i, s);
		UNPROTECT(1); /* s */
	    }
	}
	UNPROTECT(1); /* dn */
    }
    return dn;
}


/* Class validity methods ==============================================
   NB: These assume that validity methods for superclasses 
   have already been called via validObject() ...
*/

SEXP Matrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = Dim_validate(dim, "Matrix"), &pid);
    if (TYPEOF(val) != STRSXP) {
	SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	REPROTECT(val = DimNames_validate(dimnames, INTEGER(dim)), pid);
	UNPROTECT(1); /* dimnames */
    }
    UNPROTECT(2); /* val, dim */
    return val;
}

SEXP MatrixFactorization_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	val = Dim_validate(dim, "MatrixFactorization");
    UNPROTECT(1); /* dim */
    return val;
}

SEXP compMatrix_validate(SEXP obj)
{
    SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorSym)), val = NULL;
    if (TYPEOF(factors) != VECSXP)
	val = mkString(_("'factors' slot is not a list"));
    else if (LENGTH(factors) > 0) {
	SEXP nms = PROTECT(getAttrib(factors, R_NamesSymbol));
	if (isNull(nms))
	    val = mkString(_("'factors' slot has no 'names' attribute"));
	UNPROTECT(1); /* nms */
    }
    UNPROTECT(1); /* factors */
    return (val) ? val : ScalarLogical(1);
}

#define TYPEMATRIX_VALIDATE(_PREFIX_, _SEXPTYPE_, _T2C_SEXPTYPE_)	\
SEXP _PREFIX_ ## Matrix_validate(SEXP obj)				\
{									\
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)), val = NULL;		\
    if (TYPEOF(x) != _SEXPTYPE_)					\
	val =  mkString(_("'x' slot is not of type \"" #_T2C_SEXPTYPE_ "\"")); \
    UNPROTECT(1); /* x */						\
    return (val) ? val : ScalarLogical(1);				\
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
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	val = mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    UNPROTECT(1); /* dim */
    if (val) return val;
    
#ifdef ENFORCE_SYMMETRIC_DIMNAMES
    /* This check can be expensive when both rownames and colnames have 
       nonzero length, and even more so when coercions to character are 
       required ... Users can avoid the expense by setting at least one 
       of rownames and colnames to NULL or by ensuring that they are the 
       same object, as testing for pointer equality is fast ... */
    
# define ANY_TO_STRING(x)					\
    (TYPEOF(x) == STRSXP					\
     ? x							\
     : (inherits(x, "factor")					\
	? asCharacterFactor(x)					\
	: coerceVector(x, STRSXP)))
    
    SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
	ndn = PROTECT(getAttrib(dn, R_NamesSymbol));
    const char *ndn0, *ndn1;
    if (!isNull(ndn) &&
	*(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
	*(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
	strcmp(ndn0, ndn1) != 0)
	val = mkString(_("Dimnames[1] differs from Dimnames[2]"));
    else if (n > 0) {
	/* NB: It is already known that the length of 'dn[[i]]' is 0 or 'n' */ 
	SEXP rn, cn;
	if (!isNull(rn = VECTOR_ELT(dn, 0)) &&
	    !isNull(cn = VECTOR_ELT(dn, 1)) &&
	    LENGTH(rn) == n && LENGTH(cn) == n && rn != cn) {
	    PROTECT(rn = ANY_TO_STRING(rn));
	    PROTECT(cn = ANY_TO_STRING(cn));
	    if (!equal_string_vectors(rn, cn, n))
		val = mkString(_("Dimnames[1] differs from Dimnames[2]"));
	    UNPROTECT(2); /* cn, rn */
	}
    }
    UNPROTECT(2); /* ndn, dn */
    if (val) return val;
    
# undef ANY_TO_STRING
#endif

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    if (TYPEOF(uplo) != STRSXP)
	val = mkString(_("'uplo' slot is not of type \"character\""));
    else if (LENGTH(uplo) != 1)
	val = mkString(_("'uplo' slot does not have length 1"));
    else {
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
	    val = mkString(_("'uplo' slot is not \"U\" or \"L\""));
    }
    UNPROTECT(1); /* uplo */
    return (val) ? val : ScalarLogical(1);
}

SEXP triangularMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	val = mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    UNPROTECT(1); /* dim */
    if (val) return val;
	
    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    if (TYPEOF(uplo) != STRSXP)
	val = mkString(_("'uplo' slot is not of type \"character\""));
    else if (LENGTH(uplo) != 1)
	val = mkString(_("'uplo' slot does not have length 1"));
    else {
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
	    val = mkString(_("'uplo' slot is not \"U\" or \"L\""));
    }
    UNPROTECT(1); /* uplo */
    if (val) return val;
    
    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
    if (TYPEOF(diag) != STRSXP)
	val = mkString(_("'diag' slot is not of type \"character\""));
    else if (LENGTH(diag) != 1)
	val = mkString(_("'diag' slot does not have length 1"));
    else {
	const char *di = CHAR(STRING_ELT(diag, 0));
	if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
	    val = mkString(_("'diag' slot is not \"N\" or \"U\""));
    }
    UNPROTECT(1); /* diag */
    return (val) ? val : ScalarLogical(1);
}

SEXP diagonalMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	val = mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    UNPROTECT(1); /* dim */
    if (val) return val;
    
    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
    int nonunit = 0;
    if (TYPEOF(diag) != STRSXP)
	val = mkString(_("'diag' slot is not of type \"character\""));
    else if (LENGTH(diag) != 1)
	val = mkString(_("'diag' slot does not have length 1"));
    else {
	const char *di = CHAR(STRING_ELT(diag, 0));
	if (di[0] == '\0' || di[1] != '\0' || (di[0] != 'N' && di[0] != 'U'))
	    val = mkString(_("'diag' slot is not \"N\" or \"U\""));
	else
	    nonunit = di[0] == 'N';
    }
    UNPROTECT(1); /* diag */
    if (val) return val;

    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    if (nonunit) {
	if (XLENGTH(x) != n)
	    val = mkString(_("'diag' slot is \"N\" but "
			     "'x' slot does not have length n=Dim[1]"));
    } else {
	if (XLENGTH(x) != 0)
	    val = mkString(_("'diag' slot is \"U\" (identity matrix) but "
			     "'x' slot does not have length 0"));
    }
    UNPROTECT(1); /* x */
    return (val) ? val : ScalarLogical(1);
}

SEXP indMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (n == 0 && m > 0)
	val = mkString(_("m-by-0 indMatrix invalid for positive 'm'"));
    UNPROTECT(1); /* dim */
    if (val) return val;
    
    SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
    if (TYPEOF(perm) != INTSXP)
	val = mkString(_("'perm' slot is not of type \"integer\""));
    else if (XLENGTH(perm) != m)
	val = mkString(_("'perm' slot does not have length Dim[1]"));
    else {
	int *pperm = INTEGER(perm);
	while (m--) {
	    if (*pperm == NA_INTEGER) {
		val = mkString(_("'perm' slot contains NA"));
		break;
	    } else if (*pperm < 1 || *pperm > n) {
		val = mkString(_("'perm' slot has elements not in "
				 "{1,...,Dim[2]}"));
		break;
	    }
	    ++pperm;
	}
    }
    UNPROTECT(1); /* perm */
    return (val) ? val : ScalarLogical(1);
}

SEXP pMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	val = mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    else if (n > 1) {
	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	int i, *pperm = INTEGER(perm);
	char *work;
	Calloc_or_Alloca_TO(work, n, char);
	Memzero(work, n);
	--work;
	for (i = 0; i < n; ++i) {
	    if (work[*pperm]) {
		val = mkString(_("'perm' slot contains duplicates"));
		break;
	    }
	    work[*(pperm++)] = 1;
	}
	++work;
	Free_FROM(work, n);
	UNPROTECT(1); /* perm */
    }
    UNPROTECT(1); /* dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP CsparseMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1); /* dim */
    
    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
    int *pp = NULL;
    if (TYPEOF(p) != INTSXP)
	val = mkString(_("'p' slot is not of type \"integer\""));
    else if (XLENGTH(p) - 1 != n)
	val = mkString(_("'p' slot does not have length Dim[2]+1"));
    else {
	pp = INTEGER(p);
	if (pp[0] != 0)
	    val = mkString(_("first element of 'p' slot is not 0"));
	else {
	    int j;
	    for (j = 1; j <= n; ++j) {
		if (pp[j] == NA_INTEGER) {
		    val = mkString(_("'p' slot contains NA"));
		    break;
		} else if (pp[j] < pp[j-1]) {
		    val = mkString(_("'p' slot is not nondecreasing"));
		    break;
		} else if (pp[j] - pp[j-1] > m) {
		    val = mkString(_("first differences of 'p' slot exceed "
				     "Dim[1]"));
		    break;
		}
	    }
	}
    }
    if (val) {
	UNPROTECT(1); /* p */
	return val;
    }

    SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
    if (TYPEOF(i) != INTSXP)
	val = mkString(_("'i' slot is not of type \"integer\""));
    else if (XLENGTH(i) < pp[n])
	val = mkString(_("'i' slot has length less than p[length(p)]"));
    else {
	int *pi = INTEGER(i), j, k = 0, kend, ik, i0;
	for (j = 1; j <= n; ++j) {
	    kend = pp[j];
	    i0 = -1;
	    while (k < kend) {
		ik = pi[k];
		if (ik == NA_INTEGER) {
		    val = mkString(_("'i' slot contains NA"));
		    break;
		} else if (ik < 0 || ik >= m) {
		    val = mkString(_("'i' slot has elements not in "
				     "{0,...,Dim[1]-1}"));
		    break;
		} else if (ik <= i0) {
		    val = mkString(_("'i' slot is not increasing within "
				     "columns"));
		    break;
		}
		i0 = ik;
		++k;
	    }
	}
    }
    UNPROTECT(2); /* i, p */
    return (val) ? val : ScalarLogical(1);
}

SEXP RsparseMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1); /* dim */
    
    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
    int *pp = NULL;
    if (TYPEOF(p) != INTSXP)
	val = mkString(_("'p' slot is not of type \"integer\""));
    else if (XLENGTH(p) - 1 != m)
	val = mkString(_("'p' slot does not have length Dim[1]+1"));
    else {
	pp = INTEGER(p);
	if (pp[0] != 0)
	    val = mkString(_("first element of 'p' slot is not 0"));
	else {
	    int i;
	    for (i = 1; i <= m; ++i) {
		if (pp[i] == NA_INTEGER) {
		    val = mkString(_("'p' slot contains NA"));
		    break;
		} else if (pp[i] < pp[i-1]) {
		    val = mkString(_("'p' slot is not nondecreasing"));
		    break;
		} else if (pp[i] - pp[i-1] > n) {
		    val = mkString(_("first differences of 'p' slot exceed "
				     "Dim[2]"));
		    break;
		}
	    }
	}
    }
    if (val) {
	UNPROTECT(1); /* p */
	return val;
    }

    SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
    if (TYPEOF(j) != INTSXP)
	val = mkString(_("'j' slot is not of type \"integer\""));
    else if (XLENGTH(j) < pp[m])
	val = mkString(_("'j' slot has length less than p[length(p)]"));
    else {
	int *pj = INTEGER(j), i, k = 0, kend, jk, j0;
	for (i = 1; i <= m; ++i) {
	    kend = pp[i];
	    j0 = -1;
	    while (k < kend) {
		jk = pj[k];
		if (jk == NA_INTEGER) {
		    val = mkString(_("'j' slot contains NA"));
		    break;
		} else if (jk < 0 || jk >= n) {
		    val = mkString(_("'j' slot has elements not in "
				     "{0,...,Dim[2]-1}"));
		    break;
		} else if (jk <= j0) {
		    val = mkString(_("'j' slot is not increasing within "
				     "rows"));
		    break;
		}
		j0 = jk;
		++k;
	    }
	}
    }
    UNPROTECT(2); /* j, p */
    return (val) ? val : ScalarLogical(1);
}

SEXP TsparseMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    UNPROTECT(1); /* dim */
    
    SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
	j = PROTECT(GET_SLOT(obj, Matrix_jSym));
    R_xlen_t nnz;
    if (TYPEOF(i) != INTSXP)
	val = mkString(_("'i' slot is not of type \"integer\""));
    else if (TYPEOF(j) != INTSXP)
	val = mkString(_("'j' slot is not of type \"integer\""));
    else if ((nnz = XLENGTH(i)) != XLENGTH(j))
	val = mkString(_("'i' and 'j' slots do not have equal length"));
    else if (nnz == 0)
	val = ScalarLogical(1);
    else if (m == 0 || n == 0)
	val = mkString(_("'i' slot has nonzero length but prod(Dim) is 0"));
    else {
	int *pi = INTEGER(i), *pj = INTEGER(j);
	while (nnz--) {
	    if (*pi == NA_LOGICAL) {
		val = mkString(_("'i' slot contains NA"));
		break;
	    } else if (*pj == NA_LOGICAL) {
		val = mkString(_("'j' slot contains NA"));
		break;
	    } else if (*pi < 0 || *pi >= m) {
		val = mkString(_("'i' slot has elements not in "
				 "{0,...,Dim[1]-1}"));
		break;
	    } else if (*pj < 0 || *pj >= n) {
		val = mkString(_("'j' slot has elements not in "
				 "{0,...,Dim[2]-1}"));
		break;

	    }
	    ++pi;
	    ++pj;
	}
    }
    UNPROTECT(2); /* j, i */
    return (val) ? val : ScalarLogical(1);
}

SEXP sCMatrix_validate(SEXP obj)
{
    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
    int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
    if (pp[n] == 0) {
	UNPROTECT(1);
	return ScalarLogical(1);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
    SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
    int *pi = INTEGER(i), j, k = 0, kend;
    if (ul == 'U') {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] > j) {
		    UNPROTECT(2); /* i, p */
		    return mkString(_("uplo=\"U\" but there are entries below the diagonal"));
		}
		++k;
	    }
	}
    } else {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] < j) {
		    UNPROTECT(2); /* i, p */
		    return mkString(_("uplo=\"L\" but there are entries above the diagonal"));
		}
		++k;
	    }
	}
    }
    UNPROTECT(2); /* i, p */
    return ScalarLogical(1);
}

SEXP tCMatrix_validate(SEXP obj)
{
    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
    char di = *CHAR(STRING_ELT(diag, 0));
    UNPROTECT(1); /* diag */
    if (di == 'N')
	return sCMatrix_validate(obj);
    
    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
    int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
    if (pp[n] == 0) {
	UNPROTECT(1); /* p */
	return ScalarLogical(1);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
    SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
    int *pi = INTEGER(i), j, k = 0, kend;
    if (ul == 'U') {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] >= j) {
		    int s = pi[k] == j;
		    UNPROTECT(2); /* i, p */
		    return mkString(_((s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal"));
		}
		++k;
	    }
	}
    } else {
	for (j = 0; j < n; ++j) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pi[k] <= j) {
		    int s = pi[k] == j;
		    UNPROTECT(2); /* i, p */
		    return mkString(_((s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal"));
		}
		++k;
	    }
	}
    }
    UNPROTECT(2); /* i, p */
    return ScalarLogical(1);
}

SEXP sRMatrix_validate(SEXP obj)
{
    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
    int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
    if (pp[m] == 0) {
	UNPROTECT(1); /* p */
	return ScalarLogical(1);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
    SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
    int *pj = INTEGER(j), i, k = 0, kend;
    if (ul == 'U') {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] < i) {
		    UNPROTECT(2); /* j, p */
		    return mkString(_("uplo=\"U\" but there are entries below the diagonal"));
		}
		++k;
	    }
	}
    } else {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] > i) {
		    UNPROTECT(2); /* j, p */
		    return mkString(_("uplo=\"L\" but there are entries above the diagonal"));
		}
		++k;
	    }
	}
    }
    UNPROTECT(2); /* j, p */
    return ScalarLogical(1);
}

SEXP tRMatrix_validate(SEXP obj)
{
    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
    char di = *CHAR(STRING_ELT(diag, 0));
    UNPROTECT(1); /* diag */
    if (di == 'N')
	return sRMatrix_validate(obj);

    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
    int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
    if (pp[m] == 0) {
	UNPROTECT(1); /* p */
	return ScalarLogical(1);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
    SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
    int *pj = INTEGER(j), i, k = 0, kend;
    if (ul == 'U') {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] <= i) {
		    int s = pj[k] == i;
		    UNPROTECT(2); /* j, p */
		    return mkString(_((s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal"));
		}
		++k;
	    }
	}
    } else {
	for (i = 0; i < m; ++i) {
	    kend = *(++pp);
	    while (k < kend) {
		if (pj[k] >= i) {
		    int s = pj[k] == i;
		    UNPROTECT(2); /* j, p */
		    return mkString(_((s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal"));
		}
		++k;
	    }
	}
    }
    UNPROTECT(2); /* j, p */
    return ScalarLogical(1);
}

SEXP sTMatrix_validate(SEXP obj)
{
    SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
    R_xlen_t nnz = XLENGTH(i);
    if (nnz == 0) {
	UNPROTECT(1); /* i */
	return ScalarLogical(1);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */

    SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
    int *pi = INTEGER(i), *pj = INTEGER(j);
    if (ul == 'U') {
	while (nnz--)
	    if (*(pi++) > *(pj++)) {
		UNPROTECT(2); /* j, i */
		return mkString(_("uplo=\"U\" but there are entries below the diagonal"));
	    }
    } else {
	while (nnz--)
	    if (*(pi++) < *(pj++)) {
		UNPROTECT(2); /* j, i */
		return mkString(_("uplo=\"L\" but there are entries above the diagonal"));
	    }
    }
    UNPROTECT(2); /* j, i */
    return ScalarLogical(1);
}

SEXP tTMatrix_validate(SEXP obj)
{
    SEXP diag = PROTECT(GET_SLOT(obj, Matrix_diagSym));
    char di = *CHAR(STRING_ELT(diag, 0));
    UNPROTECT(1); /* diag */
    if (di == 'N')
	return sTMatrix_validate(obj);

    SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym));
    R_xlen_t nnz = XLENGTH(i);
    if (nnz == 0) {
	UNPROTECT(1); /* i */
	return ScalarLogical(1);
    }

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */

    SEXP j = PROTECT(GET_SLOT(obj, Matrix_jSym));
    int *pi = INTEGER(i), *pj = INTEGER(j);
    if (ul == 'U') {
	while (nnz--) {
	    if (*pi >= *pj) {
		int s = *pi == *pj;
		UNPROTECT(2); /* j, i */
		return mkString(_((s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"U\" but there are entries below the diagonal"));
	    }
	    ++pi;
	    ++pj;
	}
    } else {
	while (nnz--) {
	    if (*pi <= *pj) {
		int s = *pi == *pj;
		UNPROTECT(2); /* j, i */
		return mkString(_((s) ? "diag=\"U\" but there are entries on the diagonal" : "uplo=\"L\" but there are entries above the diagonal"));
	    }
	    ++pi;
	    ++pj;
	}
    }
    UNPROTECT(2); /* j, i */
    return ScalarLogical(1);
}

SEXP xgCMatrix_validate(SEXP obj)
{
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
	val = NULL;
    if (XLENGTH(x) != XLENGTH(i))
	val = mkString(_("'i' and 'x' slots do not have equal length"));
    UNPROTECT(2); /* i, x */
    return (val) ? val : ScalarLogical(1);
}

SEXP xsCMatrix_validate(SEXP obj)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = xgCMatrix_validate(obj), &pid);
    if (TYPEOF(val) != STRSXP)
	REPROTECT(val = sCMatrix_validate(obj), pid);
    UNPROTECT(1); /* val */
    return val;
}

SEXP xtCMatrix_validate(SEXP obj)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = xgCMatrix_validate(obj), &pid);
    if (TYPEOF(val) != STRSXP)
	REPROTECT(val = tCMatrix_validate(obj), pid);
    UNPROTECT(1); /* val */
    return val;
}

SEXP xgRMatrix_validate(SEXP obj)
{
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	j = PROTECT(GET_SLOT(obj, Matrix_jSym)),
	val = NULL;
    if (XLENGTH(x) != XLENGTH(j))
	val = mkString(_("'j' and 'x' slots do not have equal length"));
    UNPROTECT(2); /* j, x */
    return (val) ? val : ScalarLogical(1);
}

SEXP xsRMatrix_validate(SEXP obj)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = xgRMatrix_validate(obj), &pid);
    if (TYPEOF(val) != STRSXP)
	REPROTECT(val = sRMatrix_validate(obj), pid);
    UNPROTECT(1); /* val */
    return val;
}

SEXP xtRMatrix_validate(SEXP obj)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = xgRMatrix_validate(obj), &pid);
    if (TYPEOF(val) != STRSXP)
	REPROTECT(val = tRMatrix_validate(obj), pid);
    UNPROTECT(1); /* val */
    return val;
}

SEXP xgTMatrix_validate(SEXP obj)
{
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
	val = NULL;
    if (XLENGTH(x) != XLENGTH(i))
	val = mkString(_("'i' and 'x' slots do not have equal length"));
    UNPROTECT(2); /* i, x */
    return (val) ? val : ScalarLogical(1);
}

SEXP xsTMatrix_validate(SEXP obj)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = xgTMatrix_validate(obj), &pid);
    if (TYPEOF(val) != STRSXP)
	REPROTECT(val = sTMatrix_validate(obj), pid);
    UNPROTECT(1); /* val */
    return val;
}

SEXP xtTMatrix_validate(SEXP obj)
{
    SEXP val;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(val = xgTMatrix_validate(obj), &pid);
    if (TYPEOF(val) != STRSXP)
	REPROTECT(val = tTMatrix_validate(obj), pid);
    UNPROTECT(1); /* val */
    return val;
}

SEXP unpackedMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	val = NULL;
    int *pdim = INTEGER(dim);
    if (XLENGTH(x) != (double) pdim[0] * pdim[1])
	val = mkString(_("'x' slot does not have length prod(Dim)"));
    UNPROTECT(2); /* x, dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP packedMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	val = NULL;
    int *pdim = INTEGER(dim);
    if (XLENGTH(x) != 0.5 * pdim[0] * (pdim[0] + 1.0))
        val = mkString(_("'x' slot does not have length n*(n+1)/2, n=Dim[1]"));
    UNPROTECT(2); /* x, dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP dpoMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	val = NULL;
    double *px = REAL(x);
    int i, n = INTEGER(dim)[0];
    R_xlen_t np1 = (R_xlen_t) n + 1;
    
    /* Non-negative diagonal elements are necessary _but not_ sufficient */
    for (i = 0; i < n; ++i, px += np1) {
	if (!ISNAN(*px) && *px < 0.0) {
	    val = mkString(_("matrix is not positive semidefinite"));
	    break;
	}
    }
    UNPROTECT(2); /* x, dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP dppMatrix_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	val = NULL;
    double *px = REAL(x);
    int i, n = INTEGER(dim)[0];

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
    /* Non-negative diagonal elements are necessary _but not_ sufficient */
    if (ul == 'U') {
	for (i = 0; i < n; px += (++i)+1) {
	    if (!ISNAN(*px) && *px < 0.0) {
		val = mkString(_("matrix is not positive semidefinite"));
		break;
	    }
	}
    } else {
	for (i = 0; i < n; px += n-(i++)) {
	    if (!ISNAN(*px) && *px < 0.0) {
		val = mkString(_("matrix is not positive semidefinite"));
		break;
	    }
	}
    }
    UNPROTECT(2); /* x, dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP corMatrix_validate(SEXP obj)
{
    SEXP sd = PROTECT(GET_SLOT(obj, Matrix_sdSym)), val = NULL;
    if (TYPEOF(sd) != REALSXP)
	val = mkString(_("'sd' slot is not of type \"double\""));
    else {
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	
	if (XLENGTH(sd) != n)
	    val = mkString(_("'sd' slot does not have length n=Dim[1]"));
	else {
	    int i;
	    double *psd = REAL(sd);
	    for (i = 0; i < n; ++i) {
		if (!R_FINITE(psd[i])) {
		    val = mkString(_("'sd' slot has non-finite elements"));
		    break;
		} else if (psd[i] < 0.0) {
		    val = mkString(_("'sd' slot has negative elements"));
		    break;
		}
	    }
	}
    }
    UNPROTECT(1); /* sd */
    return (val) ? val : ScalarLogical(1);
}

SEXP Cholesky_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	val = NULL;
    double *px = REAL(x);
    int i, n = INTEGER(dim)[0];
    R_xlen_t np1 = (R_xlen_t) n + 1;
    
    /* Non-negative diagonal elements are necessary _and_ sufficient */
    for (i = 0; i < n; ++i, px += np1) {
	if (!ISNAN(*px) && *px < 0.0) {
	    val = mkString(_("matrix has negative diagonal elements"));
	    break;
	}
    }
    UNPROTECT(2); /* x, dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP pCholesky_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)),
	x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
	val = NULL;
    double *px = REAL(x);
    int i, n = INTEGER(dim)[0];

    SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
    char ul = *CHAR(STRING_ELT(uplo, 0));
    UNPROTECT(1); /* uplo */
    
    /* Non-negative diagonal elements are necessary _and_ sufficient */
    if (ul == 'U') {
	for (i = 0; i < n; px += (++i)+1) {
	    if (!ISNAN(*px) && *px < 0.0) {
		val = mkString(_("matrix has negative diagonal elements"));
		break;
	    }
	}
    } else {
	for (i = 0; i < n; px += n-(i++)) {
	    if (!ISNAN(*px) && *px < 0.0) {
		val = mkString(_("matrix has negative diagonal elements"));
		break;
	    }
	}
    }
    UNPROTECT(2); /* x, dim */
    return (val) ? val : ScalarLogical(1);
}

SEXP BunchKaufman_validate(SEXP obj)
{
    SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym)), val = NULL;
    if (TYPEOF(perm) != INTSXP)
	val = mkString(_("'perm' slot is not of type \"integer\""));
    else {
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	UNPROTECT(1); /* dim */
	
	if (XLENGTH(perm) != n)
	    val = mkString(_("'perm' slot does not have length n=Dim[1]"));
	else {
	    int n_ = n, *pperm = INTEGER(perm);
	    while (n_) {
		if (*pperm == NA_INTEGER) {
		    val = mkString(_("'perm' slot contains NA"));
		    break;
		} else if (*pperm < -n || *pperm == 0 || *pperm > n) {
		    val = mkString(_("'perm' slot has elements not in "
				     "{-n,...,n}\\{0}, n=Dim[1]"));
		    break;
		} else if (*pperm > 0) {
		    pperm += 1;
		    n_ -= 1;
		} else if (n_ > 1 && *(pperm + 1) == *pperm) {
		    pperm += 2;
		    n_ -= 2;
		} else {
		    val = mkString(_("'perm' slot has an unpaired "
				     "negative element"));
		    break;
		}
	    }
	}
    }
    UNPROTECT(1); /* perm */
    return (val) ? val : ScalarLogical(1);
}

SEXP pBunchKaufman_validate(SEXP obj)
{
    return BunchKaufman_validate(obj); /* since we only look at 'perm' */
}

SEXP Schur_validate(SEXP obj)
{
    /* MJ: assuming for simplicity that 'Q' and 'T' slots are formally valid */
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	val = mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    UNPROTECT(1); /* dim */
    if (val) return val;
    
    SEXP Q = PROTECT(GET_SLOT(obj, Matrix_QSym));
    PROTECT(dim = GET_SLOT(Q, Matrix_DimSym));
    pdim = INTEGER(dim);
    if (pdim[0] != n || pdim[1] != n)
	val = mkString(_("dimensions of 'Q' slot are not identical to 'Dim'"));
    UNPROTECT(2); /* dim, Q */
    if (val) return val;
    
    SEXP T = PROTECT(GET_SLOT(obj, Matrix_TSym));
    PROTECT(dim = GET_SLOT(T, Matrix_DimSym));
    pdim = INTEGER(dim);
    if (pdim[0] != n || pdim[1] != n)
	val = mkString(_("dimensions of 'T' slot are not identical to 'Dim'"));
    UNPROTECT(2); /* dim, T */
    if (val) return val;
    
    SEXP e = PROTECT(GET_SLOT(obj, install("EValues")));
    SEXPTYPE te = TYPEOF(e);
    if (te != REALSXP && te != CPLXSXP)
	val = mkString(_("'EValues' slot does not have type \"double\" "
			 "or type \"complex\""));
    else if (XLENGTH(e) != n)
	val = mkString(_("'EValues' slot does not have length n=Dim[1]"));
    UNPROTECT(1); /* e */
    return (val) ? val : ScalarLogical(1);
}

SEXP denseLU_validate(SEXP obj)
{
    /* MJ: assuming for simplicity that the 'Dimnames' slot is a valid list 
       partly because I'd like denseLU to formally extend dgeMatrix and in 
       that case checking 'Dimnames' here would be redundant */
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
    UNPROTECT(1); /* dim */
    
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    if (TYPEOF(x) != REALSXP)
	val = mkString(_("'x' slot is not of type \"double\""));
    else if (XLENGTH(x) != (double) m * n)
	val = mkString(_("'x' slot does not have length prod(Dim)"));
    UNPROTECT(1); /* x */
    if (val) return val;
    
    SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
    if (TYPEOF(perm) != INTSXP)
	val = mkString(_("'perm' slot is not of type \"integer\""));
    else if (XLENGTH(perm) != r)
	val = mkString(_("'perm' slot does not have length min(Dim)"));
    else {
	int *pperm = INTEGER(perm);
	while (r--) {
	    if (*pperm == NA_INTEGER) {
		val = mkString(_("'perm' slot contains NA"));
		break;
	    } else if (*pperm < 1 || *pperm > m) {
		val = mkString(_("'perm' slot has elements not in "
				 "{1,...,Dim[1]}"));
		break;
	    }
	    ++pperm;
	}
    }
    UNPROTECT(1); /* perm */
    return (val) ? val : ScalarLogical(1);
}

SEXP sparseLU_validate(SEXP obj)
{
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), n = pdim[0];
    if (pdim[1] != n)
	val = mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    UNPROTECT(1); /* dim */
    if (val) return val;
    
    SEXP L = PROTECT(GET_SLOT(obj, Matrix_LSym));
    PROTECT(dim = GET_SLOT(L, Matrix_DimSym));
    pdim = INTEGER(dim);
    if (pdim[0] != n || pdim[1] != n)
        val = mkString(_("dimensions of 'L' slot are not identical to 'Dim'"));
    else {
	SEXP uplo = PROTECT(GET_SLOT(L, Matrix_uploSym));
	if (*CHAR(STRING_ELT(uplo, 0)) == 'U')
	    val = mkString(_("'L' slot is upper (not lower) triangular"));
	UNPROTECT(1); /* uplo */
    }
    UNPROTECT(1); /* L */
    if (val) return val;
    
    SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym));
    PROTECT(dim = GET_SLOT(U, Matrix_DimSym));
    pdim = INTEGER(dim);
    if (pdim[0] != n || pdim[1] != n)
        val = mkString(_("dimensions of 'U' slot are not identical to 'Dim'"));
    else {
	SEXP uplo = PROTECT(GET_SLOT(U, Matrix_uploSym));
	if (*CHAR(STRING_ELT(uplo, 0)) != 'U')
	    val = mkString(_("'U' slot is lower (not upper) triangular"));
	UNPROTECT(1); /* uplo */
    }
    UNPROTECT(1); /* U */
    if (val) return val;
    
    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
	q = PROTECT(GET_SLOT(obj, Matrix_qSym));
    if (TYPEOF(p) != INTSXP)
	val = mkString(_("'p' slot is not of type \"integer\""));
    else if (TYPEOF(q) != INTSXP)
	val = mkString(_("'q' slot is not of type \"integer\""));
    else if (XLENGTH(p) != n)
	val = mkString(_("'p' slot does not have length Dim[1]"));
    else if (XLENGTH(q) != n)
	val = mkString(_("'q' slot does not have length Dim[1]"));
    else {
	int i, *pp = INTEGER(p), *pq = INTEGER(q);
	char *work;
	Calloc_or_Alloca_TO(work, n, char);
	Memzero(work, n);
	for (i = 0; i < n; ++i) {
	    if (*pp == NA_INTEGER) {
		val = mkString(_("'p' slot contains NA"));
		break;
	    } else if (*pp < 0 || *pp >= n) {
		val = mkString(_("'p' slot has elements not in "
				 "{0,...,Dim[1]-1}"));
		break;
	    } else if (work[*pp]  % 2) {
		val = mkString(_("'p' slot contains duplicates"));
		break;
	    } else if (*pq == NA_INTEGER) {
		val = mkString(_("'q' slot contains NA"));
		break;
	    } else if (*pq < 0 || *pq >= n) {
		val = mkString(_("'q' slot has elements not in "
				 "{0,...,Dim[1]-1}"));
		break;
	    } else if (work[*pq] >= 2) {
		val = mkString(_("'q' slot contains duplicates"));
		break;
	    }
	    work[*(pp++)] += 1;
	    work[*(pq++)] += 2;
	}
	Free_FROM(work, n);
    }
    UNPROTECT(2); /* q, p */
    return (val) ? val : ScalarLogical(1);
}

SEXP sparseQR_validate(SEXP obj)
{
    /* MJ: assuming for simplicity that 'V' and 'R' slots are formally valid */
    
    SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym)), val = NULL;
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (m < n)
	val = mkString(_("matrix has more columns than rows"));
    UNPROTECT(1); /* dim */
    if (val) return val;
    
    SEXP beta = PROTECT(GET_SLOT(obj, Matrix_betaSym));
    if (TYPEOF(beta) != REALSXP)
	val = mkString(_("'beta' slot is not of type \"double\""));
    else if (XLENGTH(beta) != n)
	val = mkString(_("'beta' slot does not have length Dim[2]"));
    UNPROTECT(1); /* beta */
    if (val) return val;

    int m2;
    SEXP V = PROTECT(GET_SLOT(obj, Matrix_VSym));
    PROTECT(dim = GET_SLOT(V, Matrix_DimSym));
    pdim = INTEGER(dim);
    if ((m2 = pdim[0]) < m)
	val = mkString(_("'V' slot has fewer than Dim[1] rows"));
    else if (pdim[1] != n)
	val = mkString(_("'V' slot does not have Dim[2] columns"));
    UNPROTECT(2); /* dim, V */
    if (val) return val;
    
    SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
    PROTECT(dim = GET_SLOT(R, Matrix_DimSym));
    pdim = INTEGER(dim);
    if (pdim[0] != m2)
	val = mkString(_("'R' slot does not have nrow(V) rows"));
    else if (pdim[1] != n)
	val = mkString(_("'R' slot does not have Dim[2] columns"));
    else {
	SEXP R_p = PROTECT(GET_SLOT(R, Matrix_pSym)),
	    R_i = PROTECT(GET_SLOT(R, Matrix_iSym));
	int *R_pp = INTEGER(R_p), *R_pi = INTEGER(R_i), j, k = 0, kend;
	for (j = 0; j < n; ++j) {
	    kend = *(++R_pp);
	    while (k < kend) {
		if (R_pi[k] > j) {
		    val = mkString(_("'R' slot must be upper trapezoidal "
				     "but has entries below the diagonal"));
		    break;
		}
		++k;
	    }
	}
	UNPROTECT(2); /* R_i, R_p */
    }
    UNPROTECT(2); /* dim, R */
    if (val) return val;

    SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
	q = PROTECT(GET_SLOT(obj, Matrix_qSym));
    if (TYPEOF(p) != INTSXP)
	val = mkString(_("'p' slot is not of type \"integer\""));
    else if (TYPEOF(q) != INTSXP)
	val = mkString(_("'q' slot is not of type \"integer\""));
    else if (XLENGTH(p) != m2)
	val = mkString(_("'p' slot does not have length nrow(V)"));
    else if (XLENGTH(q) != n && XLENGTH(q) != 0)
	val = mkString(_("'q' slot does not have length Dim[2] or length 0"));
    else {
	int i, *pp = INTEGER(p);
	char *work;
	Calloc_or_Alloca_TO(work, m2, char); /* n <= m <= m2 */
	Memzero(work, m2);
	for (i = 0; i < m2; ++i) {
	    if (*pp == NA_INTEGER) {
		val = mkString(_("'p' slot contains NA"));
		break;
	    } else if (*pp < 0 || *pp >= m2) {
		val = mkString(_("'p' slot has elements not in "
				 "{0,...,nrow(V)-1}"));
		break;
	    } else if (work[*pp]) {
		val = mkString(_("'p' slot contains duplicates"));
		break;
	    }
	    work[*(pp++)] = 1;
	}
	if (!val && LENGTH(q) == n) {
	    int *pq = INTEGER(q);
	    Memzero(work, n);
	    for (i = 0; i < n; ++i) {
		if (*pq == NA_INTEGER) {
		    val = mkString(_("'q' slot contains NA"));
		    break;
		} else if (*pq < 0 || *pq >= n) {
		    val = mkString(_("'q' slot has elements not in "
				     "{0,...,Dim[2]-1}"));
		    break;
		} else if (work[*pq]) {
		    val = mkString(_("'q' slot contains duplicates"));
		    break;
		}
		work[*(pq++)] = 1;
	    }
	}
	Free_FROM(work, m2);
    }
    UNPROTECT(2); /* q, p */
    return (val) ? val : ScalarLogical(1);
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


