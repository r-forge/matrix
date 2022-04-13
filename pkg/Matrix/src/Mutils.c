#include <limits.h>
#include <R_ext/Lapack.h>
#include "Mutils.h"

/* Slot validity methods ===============================================
   Called by various class validity methods (see below).
*/

/**
 * Test that `dim` is a length-2, non-negative integer vector.
 *
 * @param obj A `SEXP`, 
 *     typically the `Dim` slot of a (to be validated) `Matrix`.
 * @param domain A string specifying a domain for message translation.
 *
 * @return Either `TRUE` (indicating success) or a length-1 `STRSXP`
 *     containing an error message.
 */
SEXP Dim_validate(SEXP dim, const char* domain) {
    /* TODO? coerce from REALSXP to INTSXP?
       // if (TYPEOF(dim) != INTSXP && TYPEOF(dim) != REALSXP)
       //     return mkString(_("'Dim' slot is not numeric"));
       though above is not enough as we must prohibit Dim[i] > INT_MAX

       FIXME? Prohibit inherits(dim, "factor") and return a different
       error message in that case? What about S3-classed slot values,
       more generally?
    */
    if (TYPEOF(dim) != INTSXP)
	return mkString(_("'Dim' slot is not of type \"integer\""));
    if (LENGTH(dim) != 2)
	return mkString(_("'Dim' slot does not have length 2"));
    int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
    if (m < 0 || n < 0)
	return mkString(dngettext(domain,
				  "'Dim' slot contains negative value",
				  "'Dim' slot contains negative values",
				  (m < 0 && n < 0) ? 2 : 1));
    return ScalarLogical(1);
}

/**
 * Test that `dimnames` is a valid length-2 list.
 *
 * @param obj A `SEXP`,
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
    buf = Alloca(Matrix_Error_Buf_Size, char); R_CheckStack(); sprintf

    if (TYPEOF(dimnames) != VECSXP) {
	SPRINTF(buf, _("'Dimnames' slot is not a list"));
	return mkString(buf);
    }
    if (LENGTH(dimnames) != 2) {
	SPRINTF(buf, _("'Dimnames' slot does not have length 2"));
	return mkString(buf);
    }
    SEXP s;
    for (int j = 0; j < 2; ++j) {
	/* Behave as 'do_matrix()' from src/main/array.c:
	   Dimnames[[j]] must be NULL or _coercible to_ character
	   of length Dim[j] or 0 ... see 'R_Dimnames_fixup()' below
	*/
	s = VECTOR_ELT(dimnames, j);
	if (!isNull(s)) {
	    if (!isVector(s)) {
		SPRINTF(buf, _("Dimnames[[%d]] is not NULL or a vector"), j+1);
		return mkString(buf);
	    }
	    if (LENGTH(s) != pdim[j]) {
		if (LENGTH(s) != 0) {
		    SPRINTF(buf, _("length of Dimnames[[%d]] (%d) is not equal to Dim[%d] (%d)"), j+1, LENGTH(s), j+1, pdim[j]);
		    return mkString(buf);
		}
	    }
	}
    }
    return ScalarLogical(1);
}

/**
 * Test that an R object is a valid 1-character string,
 * given a set of allowed characters.
 *
 * @param s A `SEXP`, typically an S4 slot.
 * @param valid A string containing allowed characters.
 * @param nm A string naming the object being checked, for error messages.
 *
 * @return Either `TRUE` (indicating success) or a length-1 `STRSXP` 
 *     containing an error message.
 */
SEXP string_scalar_validate(SEXP s, char *valid, char *nm)
{
    char *buf;
    
    if (TYPEOF(s) != STRSXP) {
	SPRINTF(buf, _("%s is not of type \"character\""), nm);
    } else if (LENGTH(s) != 1) {
	SPRINTF(buf, _("%s does not have length 1"), nm);
    } else {
	const char *str = CHAR(STRING_ELT(s, 0));
	if (strlen(str) != 1) {
	    SPRINTF(buf, _("%s does not have string length 1"), nm);
	} else {
	    int nvalid = strlen(valid);
	    for (int i = 0; i < nvalid; ++i) {
		if (valid[i] == *str) {
		    return ScalarLogical(1);
		}
	    }
	    SPRINTF(buf, _("%s is not a character in \"%s\""), nm, valid);
	}
    }

#undef SPRINTF
    
    return mkString(buf);
}


/* Virtual class validity methods ======================================
   NB: These assume that validity methods for superclasses 
   have already been called via 'validObject()' ...
*/

SEXP Matrix_validate(SEXP obj)
{
    SEXP dim = GET_SLOT(obj, Matrix_DimSym), val = Dim_validate(dim, "Matrix");
    if (isString(val))
	return val;
    else
	return DimNames_validate(GET_SLOT(obj, Matrix_DimNamesSym), INTEGER(dim));
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

SEXP triangularMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (pdim[1] != pdim[0])
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    SEXP val = string_scalar_validate(GET_SLOT(obj, Matrix_uploSym),
				      "UL", "'uplo' slot");
    if (isString(val))
	return val;
    else
	return string_scalar_validate(GET_SLOT(obj, Matrix_diagSym),
				      "NU", "'diag' slot");
}

SEXP symmetricMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n)
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    
#undef ENFORCE_SYMMETRIC_DIMNAMES /* NOT YET */
#ifdef ENFORCE_SYMMETRIC_DIMNAMES
    /* This check is expensive when both rownames and colnames 
       have nonzero length, and even more so when a coercion to
       character is required ... Users can avoid the expense by
       setting at least one of rownames and colnames to NULL ...
     */
    
# define ANY_TO_STRING(x)					\
    (isString(x)						\
     ? x							\
     : (inherits(x, "factor")					\
	? asCharacterFactor(x)					\
	: coerceVector(x, STRSXP)))
    
    if (n > 0) {
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym), rn, cn;
	if (LENGTH(rn = VECTOR_ELT(dn, 0)) == n &&
	    LENGTH(cn = VECTOR_ELT(dn, 1)) == n &&
	    !equal_string_vectors(ANY_TO_STRING(rn), ANY_TO_STRING(cn)))
	    return mkString(_("Dimnames[[1]] differs from Dimnames[[2]]"));
    }
    
# undef ANY_TO_STRING
#endif

    return string_scalar_validate(GET_SLOT(obj, Matrix_uploSym),
				  "UL", "'uplo' slot");
}

SEXP diagonalMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym)), n = pdim[0];
    if (pdim[1] != n) {
	return mkString(_("Dim[1] != Dim[2] (matrix is not square)"));
    }
    SEXP diag = GET_SLOT(obj, Matrix_diagSym);
    SEXP val = string_scalar_validate(diag, "NU", "'diag' slot");
    if (isString(val))
	return val;
    if (*CHAR(asChar(diag)) == 'N') {
	if (LENGTH(GET_SLOT(obj, Matrix_xSym)) != n) {
	    return mkString(_("'diag' slot equal to \"N\" requires 'x' slot of length n=Dim[1]"));
	}
    } else {
	if (LENGTH(GET_SLOT(obj, Matrix_xSym)) != 0) {
	    return mkString(_("'diag' slot equal to \"U\" (identity matrix) requires 'x' slot of length 0"));
	}
    }
    return ScalarLogical(1);
}

SEXP unpackedMatrix_validate(SEXP obj)
{
    int *pdim = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    if (XLENGTH(GET_SLOT(obj, Matrix_xSym)) != pdim[0] * (R_xlen_t) pdim[1])
	return mkString(_("length of 'x' slot is not equal to prod(Dim)"));
    else
	return ScalarLogical(1);
}

SEXP packedMatrix_validate(SEXP obj)
{
    R_xlen_t n = (R_xlen_t) INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    if (2 * XLENGTH(GET_SLOT(obj, Matrix_xSym)) != n * (n + 1))
        return mkString(_("length of 'x' slot is not equal to n*(n+1)/2, n=Dim[1]"));
    else
	return ScalarLogical(1);
}

#define TYPEMATRIX_VALIDATE(_PREFIX_, _T2C_SEXPTYPE_, _SEXPTYPE_)	\
SEXP _PREFIX_ ## Matrix_validate(SEXP obj)				\
{									\
    if (TYPEOF(GET_SLOT(obj, Matrix_xSym)) != _SEXPTYPE_)		\
	return mkString(_("'x' slot is not of type \"" #_T2C_SEXPTYPE_ "\"")); \
    else								\
	return ScalarLogical(1);					\
}
/* dMatrix_validate() */
TYPEMATRIX_VALIDATE(     d,  double, REALSXP)
/* lMatrix_validate() */
TYPEMATRIX_VALIDATE(     l, logical,  LGLSXP)
/* ndenseMatrix_validate() */
/* NB: "nsparseMatrix" has no 'x' slot, only "ndenseMatrix" ... */
TYPEMATRIX_VALIDATE(ndense, logical,  LGLSXP)
/* iMatrix_validate() */
TYPEMATRIX_VALIDATE(     i, integer,  INTSXP)
/* zMatrix_validate() */
TYPEMATRIX_VALIDATE(     z, complex, CPLXSXP)
#undef TYPEMATRIX_VALIDATE


/* More for 'Dimnames' ============================================== */

/**
 * @brief Standardize user-supplied `Dimnames`.
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
	s = VECTOR_ELT(dn, i);
	if (!isNull(s) && (LENGTH(s) == 0 || !isString(s))) {
	    do_fixup = TRUE;
	    break;
	}
    }
    if (do_fixup) {
	PROTECT(dn = duplicate(dn));
	for (i = 0; i < 2; ++i) {
	    if (isNull(s = VECTOR_ELT(dn, i))) {
		continue;
	    }
	    if (LENGTH(s) == 0) {
		SET_VECTOR_ELT(dn, i, R_NilValue);
	    } else if (!isString(s)) {
		if (inherits(s, "factor")) {
		    SET_VECTOR_ELT(dn, i, asCharacterFactor(s));
		} else {
		    PROTECT(s = coerceVector(s, STRSXP));
		    SET_ATTRIB(s, R_NilValue);
		    SET_OBJECT(s, 0);
		    SET_VECTOR_ELT(dn, i, s);
		    UNPROTECT(1);
		}
	    }
	}
	UNPROTECT(1);
    }
    return dn;
}

/**
 * @brief Produce symmetric `Dimnames` from possibly asymmetric ones.
 * 
 * Called from `symmDN` in `../R/Auxiliaries.R`.
 *
 * @param dn A list of length 2, typically the `Dimnames` slot
 *     of a `symmetricMatrix`.
 * 
 * @return `rep(dn[1], 2)` if `is.null(dn[[2]])` and `!is.null(dn[[1]])`, 
 *     `rep(dn[2], 2)` otherwise.
 */
SEXP R_symmDN(SEXP dn) {
    
#define NON_TRIVIAL_DN							\
    !(isNull(VECTOR_ELT(dn, 1)) &&					\
      isNull(VECTOR_ELT(dn, 0)) &&					\
      isNull(getAttrib(dn, R_NamesSymbol)))
    
    /* Be fast (do nothing!) when dimnames = list(NULL, NULL) : */
    if (NON_TRIVIAL_DN) {
	PROTECT(dn = duplicate(dn));
	symmDN(dn, -1);
	UNPROTECT(1);
    }
    return dn;
}

void symmDN(SEXP dn, int J /* -1|0|1 */) {
    SEXP s;
    if (J < 0) {
	/* Use rownames if only rownames are non-NULL, otherwise use colnames */
	if (!isNull(s = VECTOR_ELT(dn, J = 1)) ||
	    !isNull(s = VECTOR_ELT(dn, J = 0))) {
	    SET_VECTOR_ELT(dn, !J, s);
	} else {
	    J = 1;
	}
    } else {
	/* Use rownames or colnames according to J */
    	SET_VECTOR_ELT(dn, !J, VECTOR_ELT(dn, J));	
    }
    /* names(dimnames(.)) */
    if (!isNull(s = getAttrib(dn, R_NamesSymbol))) {
	PROTECT(s);
	SET_STRING_ELT(s, !J, STRING_ELT(s, J));
	setAttrib(dn, R_NamesSymbol, s);
	UNPROTECT(1);
    }
    return;
}

/**
 * A convenience wrapper for `R_symmDN`, getting the `Dimnames` slot 
 * from a square `Matrix` and symmetrizing _if necessary_.
 *
 * @param from A square `Matrix`, typically a `symmetricMatrix`.
 *
 * @return A list `dn` of length 2 satisfying `identical(dn[1], dn[2])`,
 *     giving a suitable (symmetric) value for the `Dimnames` slot of `from`.
 */
SEXP get_symmetrized_DimNames(SEXP x) {
    return R_symmDN(GET_SLOT(x, Matrix_DimNamesSym));
}

/**
 * Symmetrize the (possibly asymmetric) `Dimnames` of a square `Matrix` 
 * and install them as the `Dimnames` of another.
 *
 * @param dest,src Square `Matrix` of equal size.
 */
void set_symmetrized_DimNames(SEXP dest, SEXP src) {
    SEXP dn = GET_SLOT(src, Matrix_DimNamesSym);
    /* Be fast (do nothing!) when dimnames = list(NULL, NULL) : */
    if (NON_TRIVIAL_DN) {
	PROTECT(dn = duplicate(dn));
	symmDN(dn, -1);
	SET_SLOT(dest, Matrix_DimNamesSym, dn);
	UNPROTECT(1);
    }
    return;
}

void set_DimNames(SEXP dest, SEXP src)
{
    SEXP dn = GET_SLOT(src, Matrix_DimNamesSym);
    /* Be fast (do nothing!) when dimnames = list(NULL, NULL) : */
    if (!(isNull(VECTOR_ELT(dn, 0)) && isNull(VECTOR_ELT(dn, 1))))
	SET_SLOT(dest, Matrix_DimNamesSym, duplicate(dn));
    return;
}


/* More for 'factors' =============================================== */

SEXP get_factor(SEXP obj, char *nm)
{
    SEXP fac = GET_SLOT(obj, Matrix_factorSym);
    R_xlen_t i = strmatch(nm, getAttrib(fac, R_NamesSymbol));
    return (i >= 0) ? VECTOR_ELT(fac, i) : R_NilValue;
}

void set_factor(SEXP obj, char *nm, SEXP val)
{
    PROTECT(val);
    SEXP fac = GET_SLOT(obj, Matrix_factorSym);
    R_xlen_t i = strmatch(nm, getAttrib(fac, R_NamesSymbol));
    if (i >= 0) {
	/* If there is already a 'nm' entry, then reset it : */
	PROTECT(fac);
	SET_VECTOR_ELT(fac, i, val);
	UNPROTECT(1);
    } else {
	/* Otherwise, install an augmented list with a 'nm' entry : */
	SET_SLOT(obj, Matrix_factorSym, append_to_named_list(fac, nm, val));
    }
    UNPROTECT(1);
    return;
}

/** 
 * @brief Subassign by name to the `factors` slot of a `compMatrix`.
 * 
 * Like `obj\@factors[[nm]] <- val`, but modifying `obj` (rather than a copy)
 * _even if_ `obj` is referenced elsewhere, supporting "automagic" caching of
 * factorizations by R functions taking `compMatrix` as an argument.
 * _Use with care!_
 * 
 * @param obj A `compMatrix`.
 * @param nm A length-1 `STRSXP` giving a factor name.
 * @param val A `SEXP`, usually a `MatrixFactorization`.
 * @param warn A length-1 `LGLSXP`. Warn if `obj` has no `factors` slot 
 *     (in which case `obj` is untouched)?
 *
 * @return `val`.
 */
SEXP R_set_factor(SEXP obj, SEXP nm, SEXP val, SEXP warn)
{
    if (R_has_slot(obj, Matrix_factorSym)) {
	PROTECT(obj);
	set_factor(obj, (char *) CHAR(asChar(nm)), val);
	UNPROTECT(1);
    } else {
	if (asLogical(warn)) {
	    warning(_("attempt to set factor on 'Matrix' without 'factors' slot"));
	}
    }
    return val;
}

/** 
 * @brief Empty the 'factors' slot of a 'compMatrix'.
 * 
 * Like `obj\@factors <- list()`, but modifying `obj` (rather than a copy)
 * _even if_ `obj` is referenced elsewhere, supporting "automagic" clearing
 * of the `factors` slot by R functions taking `compMatrix` as an argument. 
 * _Use with care!_
 * 
 * @param obj A `compMatrix`.
 * @param warn A length-1 LGLSXP. Warn if `obj` has no `factors` slot 
 *     (in which case `obj` is untouched)?
 *
 * @return `TRUE` if `obj` has a nonempty `factors` slot, `FALSE` otherwise.
 */
SEXP R_empty_factors(SEXP obj, SEXP warn)
{
    /* If there is a nonempty 'factors' slot, then replace it with list() */
    if (R_has_slot(obj, Matrix_factorSym)) {
	if (LENGTH(GET_SLOT(obj, Matrix_factorSym)) > 0) {
	    PROTECT(obj);
	    SET_SLOT(obj, Matrix_factorSym, allocVector(VECSXP, 0));
	    UNPROTECT(1);
	    return ScalarLogical(1); /* slot was reset */
	}
    } else {
	if (asLogical(warn)) {
	    warning(_("attempt to discard factors from 'Matrix' without 'factors' slot"));
	}
    }
    return ScalarLogical(0); /* no-op */
}


/* For coercions ==================================================== */

#define MAKE_TRIANGULAR_LOOP(_X_, _DIM_, _UPLO_, _DIAG_, _ZERO_, _ONE_)	\
    do {								\
        int i, j, m = _DIM_[0], n = _DIM_[1], r = (m < n) ? m : n;	\
	R_xlen_t pos;							\
									\
	if (_UPLO_ == 'U') {						\
	    pos = 1;							\
	    for (j = 0; j < r; pos += (++j)+1) {			\
		for (i = j+1; i < m; ++i) {				\
		    _X_[pos++] = _ZERO_;				\
		}							\
	    }								\
	} else {							\
	    pos = m;							\
	    for (j = 1; j < r; pos += m-(j++)) {			\
		for (i = 0; i < j; ++i) {				\
		    _X_[pos++] = _ZERO_;				\
		}							\
	    }								\
	    for (j = r; j < n; ++j) {					\
		for (i = 0; i < m; ++i) {				\
		    _X_[pos++] = _ZERO_;				\
		}							\
	    }								\
	}								\
	if (_DIAG_ == 'U') {						\
	    pos = 0;							\
	    R_xlen_t mp1 = (R_xlen_t) m + 1;				\
	    for (j = 0; j < r; ++j, pos += mp1) {			\
		_X_[pos] = _ONE_;					\
	    }								\
	}								\
    } while(0)

#define MAKE_TRIANGULAR(_PREFIX_, _TYPE_, _ZERO_, _ONE_)		\
void _PREFIX_ ## dense_unpacked_make_triangular(_TYPE_ *to, SEXP from)  \
{									\
    int *pdim = INTEGER(GET_SLOT(from, Matrix_DimSym));			\
    MAKE_TRIANGULAR_LOOP(to, pdim, *uplo_P(from), *diag_P(from),	\
			 _ZERO_, _ONE_);				\
}

/**
 * @brief Make triangular an `unpackedMatrix`.
 * 
 * Fills the "trivial" part of an `unpackedMatrix` with zeros to force 
 * triangularity, or with zeros _and_ ones to force unit triangularity.
 * The `unpackedMatrix` need not be square, though all `triangularMatrix` 
 * _are_.
 * 
 * @param to A pointer to the first element of a length-`m*n` array,
 *     usually the "data" of the `x` slot of an `unpackedMatrix`.
 * @param from A `Matrix` with `Dim`, `uplo` and `diag` slots.
 */
/* ddense_unpacked_make_triangular() */
MAKE_TRIANGULAR(d, double, 0.0, 1.0)
/* ldense_unpacked_make_triangular() */
MAKE_TRIANGULAR(l, int, 0, 1)

#undef MAKE_TRIANGULAR
#undef MAKE_TRIANGULAR_LOOP

#define MAKE_SYMMETRIC_LOOP(_X_, _DIM_, _UPLO_)				\
    do {								\
        int i, j, n = _DIM_[0];						\
	R_xlen_t upos = 0, lpos = 0;					\
									\
	if (_UPLO_ == 'U') {						\
	    for (j = 0; j < n; lpos += (++j)+1, upos = lpos) {		\
		for (i = j+1; i < n; ++i) {				\
		    _X_[++lpos] = _X_[upos += n];			\
		}							\
	    }								\
	} else {							\
	    for (j = 0; j < n; lpos += (++j)+1, upos = lpos) {		\
		for (i = j+1; i < n; ++i) {				\
		    _X_[upos += n] = _X_[++lpos];			\
		}							\
	    }								\
	}								\
    } while(0)

#define MAKE_SYMMETRIC(_PREFIX_, _TYPE_)				\
void _PREFIX_ ## dense_unpacked_make_symmetric(_TYPE_ *to, SEXP from)   \
{									\
    int *pdim = INTEGER(GET_SLOT(from, Matrix_DimSym));			\
    MAKE_SYMMETRIC_LOOP(to, pdim, *uplo_P(from));			\
}

/**
 * @brief Make symmetric a square `unpackedMatrix`.
 * 
 * Symmetrizes the elements of a square `unpackedMatrix`.
 * 
 * @param to A pointer to the first element of a length-`n*n` array,
 *     usually the "data" of the `x` slot of a square `unpackedMatrix`.
 * @param from A `Matrix` with `Dim` and `uplo` slots.
 */
/* ddense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(d, double)
/* ldense_unpacked_make_symmetric() */
MAKE_SYMMETRIC(l, int)

#undef MAKE_SYMMETRIC
#undef MAKE_SYMMETRIC_LOOP

#define MAKE_DIAGONAL(_PREFIX_, _TYPE_, _PTR_, _ONE_)		        \
void _PREFIX_ ## dense_unpacked_make_diagonal(_TYPE_ *to, SEXP from)    \
{									\
    int n = INTEGER(GET_SLOT(from, Matrix_DimSym))[0];			\
    R_xlen_t pos = 0, np1 = (R_xlen_t) n + 1;				\
    									\
    Memzero(to, n * (size_t) n);					\
    if (*diag_P(from) == 'U') {						\
	for (int j = 0; j < n; ++j, pos += np1) {			\
	    to[pos] = _ONE_;						\
	}								\
    } else {								\
	_TYPE_ *px = _PTR_(GET_SLOT(from, Matrix_xSym));		\
	for (int j = 0; j < n; ++j, pos += np1) {			\
	    to[pos] = px[j];						\
	}								\
    }									\
}

/**
 * @brief Make diagonal a square `unpackedMatrix`.
 *
 * Copy a length-`n` diagonal from a `[dl]diMatrix` to a length-`n*n` 
 * array (after zero-ing the array).
 * 
 * @param to A pointer to the first element of a length-`n*n` array,
 *     usually the "data" of the `x` slot of a square `unpackedMatrix`.
 * @param from A `[dl]diMatrix`.
 */
/* ddense_unpacked_make_diagonal() */
MAKE_DIAGONAL(d, double, REAL, 1.0)
/* ldense_unpacked_make_diagonal() */
MAKE_DIAGONAL(l, int, LOGICAL, 1)

#undef MAKE_DIAGONAL

#define PACK(_PREFIX_, _TYPE_, _ONE_)				        \
_TYPE_* _PREFIX_ ## dense_pack(_TYPE_ *dest, const _TYPE_ *src, int n,  \
			       enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag) \
{									\
    int i, j;								\
    R_xlen_t dpos = 0, spos = 0;					\
									\
    switch (uplo) {							\
    case UPP:								\
        for (j = 0; j < n; spos += n-(++j))				\
	    for (i = 0; i <= j; ++i)					\
		dest[dpos++] = src[spos++];				\
	if (diag == UNT) {						\
	    dpos = 0;							\
	    for (j = 0; j < n; dpos += (++j)+1)				\
		dest[dpos] = _ONE_;					\
	}								\
	break;								\
    case LOW:								\
    	for (j = 0; j < n; spos += (++j))				\
	    for (i = j; i <  n; ++i)					\
		dest[dpos++] = src[spos++];				\
	if (diag == UNT) {						\
	    dpos = 0;							\
	    for (j = 0; j < n; dpos += n-(j++))				\
		dest[dpos] = _ONE_;					\
	}								\
	break;								\
    default:								\
    	error(_("'uplo' must be UPP or LOW"));				\
	break;								\
    }									\
    return dest;							\
}

/**
 * @brief Pack a square `unpackedMatrix`.
 * 
 * Copies the upper or lower triangular part of `src` to `dest`,
 * where it is stored contiguously ("packed"). Optionally resets 
 * the diagonal elements to 1.
 * 
 * @param dest,src Pointers to the first elements of length-`n*(n+1)/2`
 *     and length-`n*n` (resp.) arrays, usually the "data" of the  
 *     `x` slot of a `.(sp|tp)Matrix` and `.(sy|tr)Matrix` (resp.).
 * @param n Size of matrix being packed.
 * @param uplo,diag `enum` constants specifying whether the 
 *     "nontrivial part" is upper (`UPP`) or lower (`LOW`) and
 *     whether to "force" a unit diagonal (`UNT`) or not (`NUN`).
 * 
 * @return `dest`
 */
/* ddense_pack() */
PACK(d, double, 1.0)
/* ldense_pack() */
PACK(l, int, 1)

#undef PACK

#define UNPACK(_PREFIX_, _TYPE_, _ONE_)					\
_TYPE_* _PREFIX_ ## dense_unpack(_TYPE_ *dest, const _TYPE_ *src, int n, \
				 enum CBLAS_UPLO uplo, enum CBLAS_DIAG diag) \
{									\
    int i, j;								\
    R_xlen_t dpos = 0, spos = 0;					\
    									\
    Memzero(dest, n * (size_t) n);					\
    switch (uplo) {							\
    case UPP:								\
        for (j = 0; j < n; dpos += n-(++j))				\
	    for (i = 0; i <= j; ++i)					\
		dest[dpos++] = src[spos++];				\
	break;								\
    case LOW:								\
    	for (j = 0; j < n; dpos += (++j))				\
	    for (i = j; i <  n; ++i)					\
		dest[dpos++] = src[spos++];				\
	break;								\
    default:								\
    	error(_("'uplo' must be UPP or LOW"));				\
	break;								\
    }									\
    if (diag == UNT) {							\
	dpos = 0;							\
	R_xlen_t np1 = (R_xlen_t) n + 1;				\
	for (j = 0; j < n; ++j, dpos += np1)				\
	    dest[dpos] = _ONE_;						\
    }									\
    return dest;							\
}

/**
 * @brief Unpack a `packedMatrix`.
 * 
 * Copies `src` to the upper or lower triangular part of `dest`
 * (after zero-ing `dest`), where it is stored _non_-contiguously 
 * ("unpacked"). Optionally resets the diagonal elements to 1.
 * 
 * @param dest,src Pointers to the first elements of length-`n*n` 
 *     and length-`n*(n+1)/2` (resp.) arrays, usually the "data" of the  
 *     `x` slot of a `.(sp|tp)Matrix` and `.(sy|tr)Matrix` (resp.).
 * @param n Size of matrix being unpacked.
 * @param uplo,diag `enum` constants specifying whether to copy `src`
 *     to the upper (`UPP`) or lower (`LOW`) triangle of `dest` and 
 *     whether to "force" a unit diagonal (`UNT`) or not (`NUN`).
 * 
 * @return `dest`
 */
/* ddense_unpack() */
UNPACK(d, double, 1.0)
/* ldense_unpack() */
UNPACK(l, int, 1)

#undef UNPACK

/** @brief Duplicate `.denseMatrix` (and others) as `.geMatrix`.
 *
 *  This utility supports the many `*_matrix_{prod,crossprod,tcrossprod,..}`
 *  functions that should work with both classed and unclassed matrices.
 *  It is used in many places for `.geMatrix` ("generalized") dispatch,
 *  where needed.
 *
 * @param A A `denseMatrix`, a `diagonalMatrix`, a numeric or logical `matrix`, 
 *     or a numeric or logical vector (to be handled as an `n`-by-1 `matrix`).
 * @param force A length-1 LGLSXP. Duplicate if `A` is already a `.geMatrix`?
 *
 * @param A `.geMatrix`.
 */
SEXP R_dup_mMatrix_as_geMatrix(SEXP A, SEXP force)
{
    return dup_mMatrix_as_geMatrix2(A, asLogical(force), FALSE);
}
SEXP dup_mMatrix_as_geMatrix(SEXP A, Rboolean force)
{
    return dup_mMatrix_as_geMatrix2(A, force, FALSE);
}
SEXP dup_mMatrix_as_geMatrix2(SEXP A, Rboolean force,
			      Rboolean transpose_if_vector)
{
    /* NOTA BENE: If you enlarge this list, then change '14' and '6' below!
     * ---------  '[dl]diMatrix' are no longer '[dl]denseMatrix' at R level,
     *            but are still dealt with here. 
     */
    static const char *valid[] = {
	"_NOT_A_CLASS_",
        /* defined in ./Mutils.h : */
	MATRIX_VALID_ddense, /* 14 */
	MATRIX_VALID_ldense, /*  6 */
	MATRIX_VALID_ndense, /*  5 */
	""};
    int ctype = R_check_class_etc(A, valid);

    /* Be fast if 'A' is already a '.geMatrix' */
    if (ctype == 1 | ctype == 1+14 || ctype == 1+14+6) {
	return force ? duplicate(A) : A;
    }
    
    SEXP ans, dim = R_NilValue, dimnames = R_NilValue; /* -Wall */
    int nprot = 0, *pdim;
    R_xlen_t len = 0;
    Rboolean symmetric = FALSE;
    enum dense_enum mtype = ddense;
    
    if (ctype > 0) { /* [dln]denseMatrix or [dl]diMatrix */
	
	mtype = (ctype <= 14) ? ddense : ((ctype <= 14+6) ? ldense : ndense);
	dim = GET_SLOT(A, Matrix_DimSym);
	dimnames = GET_SLOT(A, Matrix_DimNamesSym);
	
    } else { /* ctype <= 0 */

	/* We need a double or logical vector (including "matrix") */
	if (isReal(A)) {
	    mtype = ddense;
	} else if (isInteger(A)) { /* FALSE for "factor" */
	    A = PROTECT(coerceVector(A, REALSXP)); nprot++;
	    mtype = ddense;
	} else if (isLogical(A)) {
	    mtype = ldense;
	} else {
	    error(_("invalid class \"%s\" to 'dup_mMatrix_as_geMatrix()'"),
		  class_P(A));
	}

#define DUP_MMATRIX_NON_VALID(_TRANSPOSE_IF_VECTOR_)			\
	do {								\
	    ctype = 0;							\
	    if (isMatrix(A)) { /* "matrix" */				\
		dim = getAttrib(A, R_DimSymbol);			\
		dimnames = getAttrib(A, R_DimNamesSymbol);		\
		if (isNull(dimnames)) {					\
		    dimnames = PROTECT(allocVector(VECSXP, 2)); nprot++; \
		}							\
	    } else { /* non-"matrix" vector */				\
		len = XLENGTH(A);					\
		if (len > INT_MAX)					\
		    error(_("vector of length exceeding 2^31-1 to 'dup_mMatrix_as_geMatrix()'")); \
		dim = PROTECT(allocVector(INTSXP, 2)); nprot++;		\
		pdim = INTEGER(dim);					\
		if (_TRANSPOSE_IF_VECTOR_) {				\
		    pdim[0] = 1;					\
		    pdim[1] = (int) len;				\
		} else {						\
		    pdim[0] = (int) len;				\
		    pdim[1] = 1;					\
		}							\
		dimnames = PROTECT(allocVector(VECSXP, 2)); nprot++;	\
		SEXP nms = getAttrib(A, R_NamesSymbol);			\
		if (!isNull(nms)) {					\
		    SET_VECTOR_ELT(dimnames,				\
				   _TRANSPOSE_IF_VECTOR_ ? 1 : 0,	\
				   nms);				\
		}							\
	    }								\
        } while (0)
	
	DUP_MMATRIX_NON_VALID(transpose_if_vector);

    } /* ctype <= 0 */

    ans = PROTECT(NEW_OBJECT_OF_CLASS(mtype == ddense
				      ? "dgeMatrix"
				      : (mtype == ldense
					 ? "lgeMatrix"
					 : "ngeMatrix")));
    nprot++;
    pdim = INTEGER(dim);
    len = pdim[0] * (R_xlen_t) pdim[1];

    if (mtype == ddense) {
	
	/* ddenseMatrix, ddiMatrix -> dgeMatrix */
	double *px;

#define DUP_MMATRIX_CASES_ddense					\
	do {								\
	    px = REAL(ALLOC_SLOT(ans, Matrix_xSym, REALSXP, len));	\
	    switch (ctype) {						\
	    case 0: /* double vector, TYPEOF(A) == REALSXP */		\
		Memcpy(px, REAL(A), len);				\
		break;							\
	    case 1: /* dgeMatrix */					\
		Memcpy(px, REAL(GET_SLOT(A, Matrix_xSym)), len);	\
		break;							\
	    case 2: /* dtrMatrix and subclasses: */			\
	    case 9: case 10: case 11: /* Cholesky, LDL, BunchKaufman */	\
		Memcpy(px, REAL(GET_SLOT(A, Matrix_xSym)), len);	\
		ddense_unpacked_make_triangular(px, A);			\
		break;							\
	    case 3: /* dsyMatrix and subclasses: */			\
	    case 4: /* dpoMatrix and subclasses: */			\
	    case 14: /* corMatrix */					\
		Memcpy(px, REAL(GET_SLOT(A, Matrix_xSym)), len);	\
		ddense_unpacked_make_symmetric(px, A);			\
		symmetric = TRUE;					\
		break;							\
	    case 5: /* ddiMatrix */					\
		ddense_unpacked_make_diagonal(px, A);			\
		break;							\
	    case 6: /* dtpMatrix and subclasses: */			\
	    case 12: case 13: /* pCholesky, pBunchKaufman */		\
		ddense_unpack(px, REAL(GET_SLOT(A, Matrix_xSym)),	\
			      pdim[0],					\
			      *uplo_P(A) == 'U' ? UPP : LOW,		\
			      *diag_P(A) == 'N' ? NUN : UNT);		\
		ddense_unpacked_make_triangular(px, A);			\
		break;							\
	    case 7: /* dspMatrix and subclasses: */			\
	    case 8: /* dppMatrix */					\
		ddense_unpack(px, REAL(GET_SLOT(A, Matrix_xSym)),	\
			      pdim[0],					\
			      *uplo_P(A) == 'U' ? UPP : LOW,		\
			      NUN);					\
		ddense_unpacked_make_symmetric(px, A);			\
		symmetric = TRUE;					\
		break;							\
	    default:							\
		error(_("unexpected ctype=%d in 'dup_mMatrix_as_geMatrix()'"), \
		      ctype);						\
		break;							\
	    } /* switch (ctype) */					\
	} while(0)
	
	DUP_MMATRIX_CASES_ddense;

    } else { /* mtype == ldense || mtype == ndense */
	
	/* ldenseMatrix, ldiMatrix --> lgeMatrix */
	/* ndenseMatrix            --> ngeMatrix */
        int *px;

#define DUP_MMATRIX_CASES_ldense					\
	do {								\
	    px = LOGICAL(ALLOC_SLOT(ans, Matrix_xSym, LGLSXP, len));	\
	    switch (ctype) {						\
	    case 0: /* logical vector, TYPEOF(A) == LGLSXP */		\
		Memcpy(px, LOGICAL(A), len);				\
		break;							\
	    case 1+14: /* lgeMatrix */					\
	    case 1+14+6: /* ngeMatrix */				\
		Memcpy(px, LOGICAL(GET_SLOT(A, Matrix_xSym)), len);	\
		break;							\
	    case 2+14: /* ltrMatrix */					\
	    case 2+14+6: /* ntrMatrix */				\
		Memcpy(px, LOGICAL(GET_SLOT(A, Matrix_xSym)), len);	\
		ldense_unpacked_make_triangular(px, A);			\
		break;							\
	    case 3+14: /* lsyMatrix */					\
	    case 3+14+6: /* nsyMatrix */				\
		Memcpy(px, LOGICAL(GET_SLOT(A, Matrix_xSym)), len);	\
		ldense_unpacked_make_symmetric(px, A);			\
		symmetric = TRUE;					\
		break;							\
	    case 4+14: /* ldiMatrix */					\
		/* NB: ndiMatrix _does not exist_ */			\
		ldense_unpacked_make_diagonal(px, A);			\
		break;							\
	    case 5+14: /* ltpMatrix */					\
	    case 5+14+6-1: /* ntpMatrix */				\
		ldense_unpack(px, LOGICAL(GET_SLOT(A, Matrix_xSym)),	\
			      pdim[0],					\
			      *uplo_P(A) == 'U' ? UPP : LOW,		\
			      *diag_P(A) == 'N' ? NUN : UNT);		\
		ldense_unpacked_make_triangular(px, A);			\
		break;							\
	    case 6+14: /* lspMatrix */					\
	    case 6+14+6-1: /* nspMatrix */				\
		ldense_unpack(px, LOGICAL(GET_SLOT(A, Matrix_xSym)),	\
			      pdim[0],					\
			      *uplo_P(A) == 'U' ? UPP : LOW,		\
			      NUN);					\
		ldense_unpacked_make_symmetric(px, A);			\
		symmetric = TRUE;					\
		break;							\
	    default:							\
		error(_("unexpected ctype=%d in 'dup_mMatrix_as_geMatrix()'"), \
		      ctype);						\
		break;							\
	    } /* switch (ctype) */					\
	} while(0)

	DUP_MMATRIX_CASES_ldense;

    } /* mtype == ldense || mtype == ndense */

#define DUP_MMATRIX_FINALIZE(_SYMMETRIC_)				\
    do {								\
	SET_SLOT(ans, Matrix_DimSym, duplicate(dim));			\
	SET_SLOT(ans, Matrix_DimNamesSym,				\
		 _SYMMETRIC_ ? R_symmDN(dimnames) : duplicate(dimnames)); \
    } while (0)

    DUP_MMATRIX_FINALIZE(symmetric);
    UNPROTECT(nprot);
    return ans;
}

SEXP R_dup_mMatrix_as_dgeMatrix(SEXP A, SEXP force)
{
    return dup_mMatrix_as_dgeMatrix2(A, asLogical(force), FALSE);
}
SEXP dup_mMatrix_as_dgeMatrix(SEXP A, Rboolean force)
{
    return dup_mMatrix_as_dgeMatrix2(A, force, FALSE);
}
SEXP dup_mMatrix_as_dgeMatrix2(SEXP A, Rboolean force,
			       Rboolean transpose_if_vector)
{
    static const char *valid[] = {"_NOT_A_CLASS_", MATRIX_VALID_ddense, ""};
    int ctype = R_check_class_etc(A, valid);

    /* Be fast if 'A' is already a 'dgeMatrix' */
    if (ctype == 1) {
	return force ? duplicate(A) : A;
    }

    SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
	dim = R_NilValue, /* -Wall */
	dimnames = R_NilValue; /* -Wall */
    int nprot = 1, *pdim;
    R_xlen_t len = 0;
    Rboolean symmetric = FALSE;
    double *px;

    if (ctype > 0) { /* ddenseMatrix or ddiMatrix */
	dim = GET_SLOT(A, Matrix_DimSym);
	dimnames = GET_SLOT(A, Matrix_DimNamesSym);
    } else { /* ctype <= 0 */
	if (!isReal(A)) {
	    if (isInteger(A) || isLogical(A)) { /* FALSE for "factor" */
		A = PROTECT(coerceVector(A, REALSXP)); nprot++;
	    } else {
		error(_("invalid class \"%s\" to 'dup_mMatrix_as_dgeMatrix()'"),
		      class_P(A));
	    }
	}
	DUP_MMATRIX_NON_VALID(transpose_if_vector); /* updates 'dim(names)?' */
    }
    pdim = INTEGER(dim);
    len = pdim[0] * (R_xlen_t) pdim[1];
    DUP_MMATRIX_CASES_ddense; /* updates 'symmetric' */
    DUP_MMATRIX_FINALIZE(symmetric);
    UNPROTECT(nprot);
    return ans;
}


/* "General" purpose ================================================ */

/* That both 's1' and 's2' are STRSXP must be checked by the caller ...
   see symmetricMatrix_validate() above (currently the only use case)
*/
Rboolean equal_string_vectors(SEXP s1, SEXP s2)
{
    int n = LENGTH(s1);
    if (LENGTH(s2) != n) {
	return FALSE;
    }
    for (int i = 0; i < n; ++i) {
	/* Note that 'R_compute_identical()' in src/main/identical.c
	   is careful to distinguish between NA_STRING and "NA" in STRSXP, 
	   but we need not be here ...
	   
	   MJ: Why not?
	*/
	if (strcmp(CHAR(STRING_ELT(s1, i)), CHAR(STRING_ELT(s2, i)))) {
	    return FALSE;
	}
    }
    return TRUE;
}

/* That 'valid' is a STRSXP must be checked by the caller ...
   NOTE: ./Mutils.h has
         int Matrix_check_class_(char *x, const char **valid);
	 int Matrix_check_class(SEXP x, const char **valid);
	 ... and this is yet another variant ...
*/
R_xlen_t strmatch(char *x, SEXP valid)
{
    R_xlen_t len = xlength(valid);
    for (R_xlen_t i = 0; i < len; ++i)
	if (!strcmp(x, CHAR(STRING_ELT(valid, i)))) return i;
    return -1;
}

#define APPEND_TO_NAMED(_T2C_SEXPTYPE_, _SEXPTYPE_, _CTYPE_,		\
			_XPTR_, _YPTR_, _COPY_, _APPEND_)		\
SEXP append_to_named_ ## _T2C_SEXPTYPE_(SEXP x, char *nm, _CTYPE_ val)  \
{									\
    R_xlen_t len = XLENGTH(x);						\
    SEXP nx = getAttrib(x, R_NamesSymbol),				\
	y = PROTECT(allocVector(_SEXPTYPE_, len + 1)),			\
	ny = PROTECT(allocVector(STRSXP, len + 1));			\
    _XPTR_; _YPTR_;							\
    for (R_xlen_t i = 0; i < len; ++i) {				\
	_COPY_;								\
	SET_STRING_ELT(ny, i, STRING_ELT(nx, i));			\
    }									\
    _APPEND_;								\
    SET_STRING_ELT(ny, len, mkChar(nm));				\
    setAttrib(y, R_NamesSymbol, ny);					\
    UNPROTECT(2);							\
    return y;								\
}
/* append_to_named_list() */
APPEND_TO_NAMED(list, VECSXP, SEXP,
		,
		,
		SET_VECTOR_ELT(y, i, VECTOR_ELT(x, i)),
		SET_VECTOR_ELT(y, len, val))
#undef APPEND_TO_NAMED


/* "General purpose" ================================================ */
    
/* La_norm_type() and La_rcond_type() have been in src/include/R_ext/Lapack.h
   and later in src/modules/lapack/Lapack.c but have still not been available 
   to package writers ...
*/
char La_norm_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* aliases */
    else if (typup == 'E')
	typup = 'F';
    else if (typup != 'M' && typup != 'O' && typup != 'I' && typup != 'F')
	error(_("argument type[1]='%s' must be one of 'M','1','O','I','F', or 'E'"),
	      typstr);
    return typup;
}

char La_rcond_type(const char *typstr)
{
    char typup;

    if (strlen(typstr) != 1)
	error(_("argument type[1]='%s' must be a character string of string length 1"),
	      typstr);
    typup = (char) toupper(*typstr);
    if (typup == '1')
	typup = 'O'; /* alias */
    else if (typup != 'O' && typup != 'I')
	error(_("argument type[1]='%s' must be one of '1','O', or 'I'"),
	      typstr);
    return typup; /* 'O' or 'I' */
}

SEXP as_det_obj(double mod, int log, int sign)
{
    SEXP det = PROTECT(allocVector(VECSXP, 2)),
	nms = PROTECT(allocVector(STRSXP, 2)),
	val = PROTECT(ScalarReal(mod));

    setAttrib(det, R_NamesSymbol, nms);
    SET_STRING_ELT(nms, 0, mkChar("modulus"));
    SET_STRING_ELT(nms, 1, mkChar("sign"));
    setAttrib(val, install("logarithm"), ScalarLogical(log));
    SET_VECTOR_ELT(det, 0, val);
    SET_VECTOR_ELT(det, 1, ScalarInteger(sign));
    setAttrib(det, R_ClassSymbol, mkString("det"));
    UNPROTECT(3);
    return det;
}

/* MJ: No longer needed ... replacement in ./packedMatrix.c */
#if 0

/**
 * Copy the diagonal elements of the packed denseMatrix x to dest
 *
 * @param dest vector of length ncol(x)
 * @param x (pointer to) a "d?pMatrix" object
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

/** diag(x) <- D  for   x a  <dspMatrix>  or dppMatrix, ..etc
 */
SEXP d_packed_setDiag(double *diag, int l_d, SEXP x, int n)
{
#define SET_packed_setDiag				\
    SEXP ret = PROTECT(duplicate(x)),			\
	r_x = GET_SLOT(ret, Matrix_xSym);		\
    Rboolean d_full = (l_d == n);			\
    if (!d_full && l_d != 1)				\
	error(_("replacement diagonal has wrong length"))

#define END_packed_setDiag						\
    int j, pos = 0;							\
									\
    if (*uplo_P(x) == 'U') {						\
	if(d_full)							\
	    for(pos= 0, j=0; j < n; pos += 1+(++j))	 xx[pos] = diag[j]; \
	else /* l_d == 1 */						\
	    for(pos= 0, j=0; j < n; pos += 1+(++j))	 xx[pos] = *diag; \
    } else {								\
	if(d_full)							\
	    for(pos= 0, j=0; j < n; pos += (n - j), j++) xx[pos] = diag[j]; \
	else /* l_d == 1 */						\
	    for(pos= 0, j=0; j < n; pos += (n - j), j++) xx[pos] = *diag; \
    }									\
    UNPROTECT(1);							\
    return ret

    SET_packed_setDiag; double *xx = REAL(r_x);
    END_packed_setDiag;
}

SEXP l_packed_setDiag(int *diag, int l_d, SEXP x, int n)
{
    SET_packed_setDiag; int *xx = LOGICAL(r_x);
    END_packed_setDiag;
}

#define tr_END_packed_setDiag						\
    if (*diag_P(x) == 'U') { /* uni-triangular */			\
	/* after setting, typically is not uni-triangular anymore: */	\
	SEXP ch_N = PROTECT(mkChar("N"));				\
	SET_STRING_ELT(GET_SLOT(ret, Matrix_diagSym), 0, ch_N);		\
	UNPROTECT(1);							\
    }									\
    END_packed_setDiag

SEXP tr_d_packed_setDiag(double *diag, int l_d, SEXP x, int n)
{
    SET_packed_setDiag; double *xx = REAL(r_x);
    tr_END_packed_setDiag;
}

SEXP tr_l_packed_setDiag(int *diag, int l_d, SEXP x, int n)
{
    SET_packed_setDiag; int *xx = LOGICAL(r_x);
    tr_END_packed_setDiag;
}

#undef SET_packed_setDiag
#undef END_packed_setDiag
#undef tr_END_packed_setDiag

void tr_d_packed_getDiag(double *dest, SEXP x, int n)
{
    if (*diag_P(x) == 'U') {
	for (int j = 0; j < n; j++) dest[j] = 1.;
    } else {
	d_packed_getDiag(dest, x, n);
    }
    return;
}

void tr_l_packed_getDiag(   int *dest, SEXP x, int n)
{
    if (*diag_P(x) == 'U')
	for (int j = 0; j < n; j++) dest[j] = 1;
    else
	l_packed_getDiag(dest, x, n);
    return;
}

/* These two *_addDiag() were unused and not replaced */
SEXP d_packed_addDiag(double *diag, int l_d, SEXP x, int n)
{
    SEXP ret = PROTECT(duplicate(x)),
	r_x = GET_SLOT(ret, Matrix_xSym);
    double *xx = REAL(r_x);
    int j, pos = 0;

    if (*uplo_P(x) == 'U') {
	for(pos= 0, j=0; j < n; pos += 1+(++j))	     xx[pos] += diag[j];
    } else {
	for(pos= 0, j=0; j < n; pos += (n - j), j++) xx[pos] += diag[j];
    }
    UNPROTECT(1);
    return ret;
}

SEXP tr_d_packed_addDiag(double *diag, int l_d, SEXP x, int n)
{
    SEXP ret = PROTECT(d_packed_addDiag(diag, l_d, x, n));
    if (*diag_P(x) == 'U') { /* uni-triangular */
	SEXP ch_N = PROTECT(mkChar("N"));
	SET_STRING_ELT(GET_SLOT(ret, Matrix_diagSym), 0, ch_N);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return ret;
}

#endif /* MJ */

#if 0 /* unused */

double get_double_by_name(SEXP obj, char *nm)
{
    SEXP nms = PROTECT(getAttrib(obj, R_NamesSymbol));
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error(_("object must be a named, numeric vector"));
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    UNPROTECT(1);
	    return REAL(obj)[i];
	}
    }
    UNPROTECT(1);
    return R_NaReal;
}

SEXP set_double_by_name(SEXP obj, double val, char *nm)
{
    SEXP nms = PROTECT(getAttrib(obj, R_NamesSymbol));
    int i, len = length(obj);

    if ((!isReal(obj)) || (length(obj) > 0 && nms == R_NilValue))
	error(_("object must be a named, numeric vector"));
    // case 1:  replace existing entry named  <nm>
    for (i = 0; i < len; i++) {
	if (!strcmp(nm, CHAR(STRING_ELT(nms, i)))) {
	    REAL(obj)[i] = val;
	    UNPROTECT(1);
	    return obj;
	}
    }
    // case 2:  no such name --> add new entry with that name at end of vec
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
	UNPROTECT(2);
	return nx;
    }
}

/* useful for all the ..CMatrix classes (and ..R by [0] <-> [1]) */
SEXP CMatrix_set_Dim(SEXP x, int nrow)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));

    dims[0] = nrow;
    dims[1] = length(GET_SLOT(x, Matrix_pSym)) - 1;
    return x;
}

SEXP new_dgeMatrix(int nrow, int ncol)
{
    SEXP ans = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix")),
	 ad = PROTECT(allocVector(INTSXP, 2));

    INTEGER(ad)[0] = nrow;
    INTEGER(ad)[1] = ncol;
    SET_SLOT(ans, Matrix_DimSym, ad);
    SET_SLOT(ans, Matrix_DimNamesSym, allocVector(VECSXP, 2));
    ALLOC_SLOT(ans, Matrix_xSym, REALSXP, ((R_xlen_t) nrow) * ncol);

    UNPROTECT(2);
    return ans;
}

#endif /* unused */

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
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param ij: 2-column integer matrix
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd(SEXP ij, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
    SEXP ans;
    int *ij_di = NULL, n, nprot=1;
    Rboolean check_bounds = asLogical(chk_bnds), one_ind = asLogical(orig_1);

    if(TYPEOF(di) != INTSXP) {di = PROTECT(coerceVector(di, INTSXP)); nprot++; }
    if(TYPEOF(ij) != INTSXP) {ij = PROTECT(coerceVector(ij, INTSXP)); nprot++; }
    if(!isMatrix(ij) ||
       (ij_di = INTEGER(getAttrib(ij, R_DimSymbol)))[1] != 2)
	error(_("Argument ij must be 2-column integer matrix"));
    n = ij_di[0];
    int *Di = INTEGER(di), *IJ = INTEGER(ij),
	*j_ = IJ+n;/* pointer offset! */

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];
#define do_ii_FILL(_i_, _j_)						\
	int i;								\
	if(check_bounds) {						\
	    for(i=0; i < n; i++) {					\
		if(_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER)	\
		    ii[i] = NA_INTEGER;					\
		else {							\
		    register int i_i, j_i;				\
	            if(one_ind) { i_i = _i_[i]-1; j_i = _j_[i]-1; }	\
	            else        { i_i = _i_[i]  ; j_i = _j_[i]  ; }	\
		    if(i_i < 0 || i_i >= Di[0])				\
			error(_("subscript 'i' out of bounds in M[ij]")); \
		    if(j_i < 0 || j_i >= Di[1])				\
			error(_("subscript 'j' out of bounds in M[ij]")); \
		    ii[i] = i_i + j_i * nr;				\
		}							\
	    }								\
	} else {							\
	    for(i=0; i < n; i++)					\
		ii[i] = (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER)	\
		    ? NA_INTEGER					\
 	            : (one_ind ? ((_i_[i]-1) + (_j_[i]-1)*nr)		\
	                       :   _i_[i]    +  _j_[i]   *nr);		\
	}

	do_ii_FILL(IJ, j_);
    } else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(IJ, j_);
    }
    UNPROTECT(nprot);
    return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param i: integer vector
 * @param j: integer vector of same length as 'i'
 * @param orig_1: logical: if TRUE, "1-origin" otherwise "0-origin"
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
    SEXP ans;
    int n = LENGTH(i), nprot = 1;
    Rboolean check_bounds = asLogical(chk_bnds), one_ind = asLogical(orig_1);

    if(TYPEOF(di)!= INTSXP) {di = PROTECT(coerceVector(di,INTSXP)); nprot++; }
    if(TYPEOF(i) != INTSXP) { i = PROTECT(coerceVector(i, INTSXP)); nprot++; }
    if(TYPEOF(j) != INTSXP) { j = PROTECT(coerceVector(j, INTSXP)); nprot++; }
    if(LENGTH(j) != n)
	error(_("i and j must be integer vectors of the same length"));
    int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

    if((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
	ans = PROTECT(allocVector(REALSXP, n));
	double *ii = REAL(ans), nr = (double) Di[0];

	do_ii_FILL(i_, j_);
    } else {
	ans = PROTECT(allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(i_, j_);
    }
    UNPROTECT(nprot);
    return ans;
}
#undef do_ii_FILL

// Almost "Cut n Paste" from ...R../src/main/array.c  do_matrix() :
// used in ../R/Matrix.R as
//
// .External(Mmatrix,
//	     data, nrow, ncol, byrow, dimnames,
//	     missing(nrow), missing(ncol))
SEXP Mmatrix(SEXP args)
{
    SEXP vals, ans, snr, snc, dimnames;
    int nr = 1, nc = 1, byrow, miss_nr, miss_nc;
    R_xlen_t lendat;

    args = CDR(args); /* skip 'name' */
    vals = CAR(args); args = CDR(args);
    /* Supposedly as.vector() gave a vector type, but we check */
    switch(TYPEOF(vals)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
	case STRSXP:
	case RAWSXP:
	case EXPRSXP:
	case VECSXP:
	    break;
	default:
	    error(_("'data' must be of a vector type"));
    }
    lendat = XLENGTH(vals);
    snr = CAR(args); args = CDR(args);
    snc = CAR(args); args = CDR(args);
    byrow = asLogical(CAR(args)); args = CDR(args);
    if (byrow == NA_INTEGER)
	error(_("invalid '%s' argument"), "byrow");
    dimnames = CAR(args);
    args = CDR(args);
    miss_nr = asLogical(CAR(args)); args = CDR(args);
    miss_nc = asLogical(CAR(args));

    if (!miss_nr) {
	if (!isNumeric(snr)) error(_("non-numeric matrix extent"));
	nr = asInteger(snr);
	if (nr == NA_INTEGER)
	    error(_("invalid 'nrow' value (too large or NA)"));
	if (nr < 0)
	    error(_("invalid 'nrow' value (< 0)"));
    }
    if (!miss_nc) {
	if (!isNumeric(snc)) error(_("non-numeric matrix extent"));
	nc = asInteger(snc);
	if (nc == NA_INTEGER)
	    error(_("invalid 'ncol' value (too large or NA)"));
	if (nc < 0)
	    error(_("invalid 'ncol' value (< 0)"));
    }
    if (miss_nr && miss_nc) {
	if (lendat > INT_MAX) error("data is too long");
	nr = (int) lendat;
    } else if (miss_nr) {
	if (lendat > (double) nc * INT_MAX) error("data is too long");
	nr = (int) ceil((double) lendat / (double) nc);
    } else if (miss_nc) {
	if (lendat > (double) nr * INT_MAX) error("data is too long");
	nc = (int) ceil((double) lendat / (double) nr);
    }

    if(lendat > 0) {
	R_xlen_t nrc = (R_xlen_t) nr * nc;
	if (lendat > 1 && nrc % lendat != 0) {
	    if (((lendat > nr) && (lendat / nr) * nr != lendat) ||
		((lendat < nr) && (nr / lendat) * lendat != nr))
		warning(_("data length [%d] is not a sub-multiple or multiple of the number of rows [%d]"), lendat, nr);
	    else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
		     ((lendat < nc) && (nc / lendat) * lendat != nc))
		warning(_("data length [%d] is not a sub-multiple or multiple of the number of columns [%d]"), lendat, nc);
	}
	else if ((lendat > 1) && (nrc == 0)){
	    warning(_("data length exceeds size of matrix"));
	}
    }

#ifndef LONG_VECTOR_SUPPORT
   if ((double)nr * (double)nc > INT_MAX)
	error(_("too many elements specified"));
#endif

    PROTECT(ans = allocMatrix(TYPEOF(vals), nr, nc));
    if(lendat) {
	if (isVector(vals))
	    copyMatrix(ans, vals, byrow);
	else
	    copyListMatrix(ans, vals, byrow);
    } else if (isVector(vals)) { /* fill with NAs */
	R_xlen_t N = (R_xlen_t) nr * nc, i;
	switch(TYPEOF(vals)) {
	case STRSXP:
	    for (i = 0; i < N; i++)
		SET_STRING_ELT(ans, i, NA_STRING);
	    break;
	case LGLSXP:
	    for (i = 0; i < N; i++)
		LOGICAL(ans)[i] = NA_LOGICAL;
	    break;
	case INTSXP:
	    for (i = 0; i < N; i++)
		INTEGER(ans)[i] = NA_INTEGER;
	    break;
	case REALSXP:
	    for (i = 0; i < N; i++)
		REAL(ans)[i] = NA_REAL;
	    break;
	case CPLXSXP:
	    {
		Rcomplex na_cmplx;
		na_cmplx.r = NA_REAL;
		na_cmplx.i = 0;
		for (i = 0; i < N; i++)
		    COMPLEX(ans)[i] = na_cmplx;
	    }
	    break;
	case RAWSXP:
	    // FIXME:  N may overflow size_t !!
	    memset(RAW(ans), 0, N);
	    break;
	default:
	    /* don't fill with anything */
	    ;
	}
    }
    if(!isNull(dimnames)&& length(dimnames) > 0)
	ans = dimnamesgets(ans, dimnames);
    UNPROTECT(1);
    return ans;
}

/**
 * From the two 'x' slots of two dense matrices a and b,
 * compute the 'x' slot of rbind(a, b)
 *
 * Currently, an auxiliary only for setMethod rbind2(<denseMatrix>, <denseMatrix>)
 * in ../R/bind2.R
 *
 * @param a
 * @param b
 *
 * @return
 */SEXP R_rbind2_vector(SEXP a, SEXP b) {
    int *d_a = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*d_b = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	n1 = d_a[0], m = d_a[1],
	n2 = d_b[0];
    if(d_b[1] != m)
	error(_("the number of columns differ in R_rbind2_vector: %d != %d"),
	      m, d_b[1]);
    SEXP
	a_x = GET_SLOT(a, Matrix_xSym),
	b_x = GET_SLOT(b, Matrix_xSym);
    int nprot = 1;
    // Care: can have "ddenseMatrix" "ldenseMatrix" or "ndenseMatrix"
    if(TYPEOF(a_x) != TYPEOF(b_x)) { // choose the "common type"
	// Now know: either LGLSXP or REALSXP. FIXME for iMatrix, zMatrix,..
	if(TYPEOF(a_x) != REALSXP) {
	    a_x = PROTECT(duplicate(coerceVector(a_x, REALSXP))); nprot++;
	} else if(TYPEOF(b_x) != REALSXP) {
	    b_x = PROTECT(duplicate(coerceVector(b_x, REALSXP))); nprot++;
	}
    }

    SEXP ans = PROTECT(allocVector(TYPEOF(a_x), m * (n1 + n2)));
    int ii = 0;
    switch(TYPEOF(a_x)) {
    case LGLSXP: {
	int
	    *r = LOGICAL(ans),
	    *ax= LOGICAL(a_x),
	    *bx= LOGICAL(b_x);

#define COPY_a_AND_b_j					\
	for(int j=0; j < m; j++) {			\
	    Memcpy(r+ii, ax+ j*n1, n1); ii += n1;	\
	    Memcpy(r+ii, bx+ j*n2, n2); ii += n2;	\
	} ; break

	COPY_a_AND_b_j;
    }
    case REALSXP: {
	double
	    *r = REAL(ans),
	    *ax= REAL(a_x),
	    *bx= REAL(b_x);

	COPY_a_AND_b_j;
    }
    } // switch
    UNPROTECT(nprot);
    return ans;
}

#define TRUE_  ScalarLogical(1)
#define FALSE_ ScalarLogical(0)

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// all0     <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
// allFalse <- function(x) !any(x) && !any(is.na(x)) ## ~= all0
SEXP R_all0(SEXP x) {
    if (!isVectorAtomic(x)) {
	if(length(x) == 0) return TRUE_;
	// Typically S4.  TODO: Call the R code above, instead!
	error(_("Argument must be numeric-like atomic vector"));
    }
    R_xlen_t i, n = XLENGTH(x);
    if(n == 0) return TRUE_;

    switch(TYPEOF(x)) {
    case LGLSXP: {
	int *xx = LOGICAL(x);
	for(i=0; i < n; i++)
	    if(xx[i] == NA_LOGICAL || xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    case INTSXP: {
	int *xx = INTEGER(x);
	for(i=0; i < n; i++)
	    if(xx[i] == NA_INTEGER || xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    case REALSXP: {
	double *xx = REAL(x);
	for(i=0; i < n; i++)
	    if(ISNAN(xx[i]) || xx[i] != 0.) return FALSE_;
	return TRUE_;
    }
    case RAWSXP: {
	unsigned char *xx = RAW(x);
	for(i=0; i < n; i++)
	    if(xx[i] != 0) return FALSE_;
	return TRUE_;
    }
    }
    error(_("Argument must be numeric-like atomic vector"));
    return R_NilValue; // -Wall
}

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// any0 <- function(x) isTRUE(any(x == 0)) ## ~= anyFalse
// anyFalse <- function(x) isTRUE(any(!x)) ## ~= any0
SEXP R_any0(SEXP x) {
    if (!isVectorAtomic(x)) {
	if(length(x) == 0) return FALSE_;
	// Typically S4.  TODO: Call the R code above, instead!
	error(_("Argument must be numeric-like atomic vector"));
    }
    R_xlen_t i, n = XLENGTH(x);
    if(n == 0) return FALSE_;

    switch(TYPEOF(x)) {
    case LGLSXP: {
	int *xx = LOGICAL(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    case INTSXP: {
	int *xx = INTEGER(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    case REALSXP: {
	double *xx = REAL(x);
	for(i=0; i < n; i++) if(xx[i] == 0.) return TRUE_;
	return FALSE_;
    }
    case RAWSXP: {
	unsigned char *xx = RAW(x);
	for(i=0; i < n; i++) if(xx[i] == 0) return TRUE_;
	return FALSE_;
    }
    }
    error(_("Argument must be numeric-like atomic vector"));
    return R_NilValue; // -Wall
}

#undef TRUE_
#undef FALSE_

/**
 * A safe NEW_OBJECT(MAKE_CLASS(what)), where the caller must protect the
 * return value of this function
 *
 * @param an R character string specifying the name of a known S4 class
 */
SEXP NEW_OBJECT_OF_CLASS(const char* what)
{
    SEXP ans = NEW_OBJECT(PROTECT(MAKE_CLASS(what)));
    UNPROTECT(1);
    return ans;
}
