#include "Mdefines.h"

#define     MK(_FORMAT_     )    Rf_mkString(_FORMAT_             )
#define     MS(_FORMAT_, ...) Matrix_sprintf(_FORMAT_, __VA_ARGS__)

#define    RMK(_FORMAT_     ) \
	return MK(   _FORMAT_              )
#define    RMS(_FORMAT_, ...) \
	return    MS(_FORMAT_, __VA_ARGS__)
#define  RMKMS(_FORMAT_, ...) \
	return MK(MS(_FORMAT_, __VA_ARGS__))

#define   FRMK(_FORMAT_     ) \
	do { \
		Matrix_Free(work, lwork); \
		RMK  (_FORMAT_             ); \
	} while (0)
#define   FRMS(_FORMAT_, ...) \
	do { \
		Matrix_Free(work, lwork); \
		RMS  (_FORMAT_, __VA_ARGS__); \
	} while (0)
#define FRMKMS(_FORMAT_, ...) \
	do { \
		Matrix_Free(work, lwork); \
		RMKMS(_FORMAT_, __VA_ARGS__); \
	} while (0)


/* Slot validity methods ===============================================
   Called by various class validity methods (see below).
*/

static
char *valid_slot_Dim(SEXP dim)
{
	if (TYPEOF(dim) != INTSXP)
		RMS(_("'%s' slot is not of type \"%s\""), "Dim", "integer");
	if (XLENGTH(dim) != 2)
		RMS(_("'%s' slot does not have length %d"), "Dim", 2);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m == NA_INTEGER || n == NA_INTEGER)
		RMS(_("'%s' slot contains NA"), "Dim");
	if (m < 0 || n < 0)
		RMS(_("'%s' slot has negative elements"), "Dim");

	return NULL;
}

static
char *valid_slot_Dimnames(SEXP dimnames, SEXP dim)
{
	if (TYPEOF(dimnames) != VECSXP)
		RMS(_("'%s' slot is not a list"), "Dimnames");
	if (XLENGTH(dimnames) != 2)
		RMS(_("'%s' slot does not have length %d"), "Dimnames", 2);

	/* Behave as do_matrix() from src/main/array.c:
	   Dimnames[[i]] must be NULL or _coercible to_ character
	   of length Dim[i] or 0 ... see R_DimNames_fixup [./attrib.c]
	*/

	SEXP s;
	int i, *pdim = INTEGER(dim);
	R_xlen_t ns;

	for (i = 0; i < 2; ++i) {
		s = VECTOR_ELT(dimnames, i);
		if (s == R_NilValue)
			continue;
		if (!Rf_isVector(s))
			RMS(_("%s[[%d]] is not NULL or a vector"), "Dimnames", i + 1);
		ns = XLENGTH(s);
		if (ns != pdim[i] && ns != 0)
			RMS(_("length of %s[[%d]] (%lld) is not equal to %s[%d] (%d)"),
			    "Dimnames", i + 1, (long long) ns,
			    "Dim"     , i + 1,        pdim[i]);
	}

	return NULL;
}

SEXP R_valid_slot_Dim(SEXP dim)
{
	char *msg = valid_slot_Dim(dim);
	return (msg) ? Rf_mkString(msg) : Rf_ScalarLogical(1);
}

SEXP R_valid_slot_Dimnames(SEXP dimnames, SEXP dim)
{
	char *msg = valid_slot_Dimnames(dimnames, dim);
	return (msg) ? Rf_mkString(msg) : Rf_ScalarLogical(1);
}


/* Class validity methods ==============================================
   NB: These assume that validity methods for superclasses
   have already been called via validObject() ...
*/

SEXP R_valid_Matrix(SEXP obj)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	char *msg = valid_slot_Dim(dim);
	if (!msg) {
		SEXP dimnames = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
		msg = valid_slot_Dimnames(dimnames, dim);
		UNPROTECT(1); /* dimnames */
	}
	UNPROTECT(1); /* dim */
	return (msg) ? Rf_mkString(msg) : Rf_ScalarLogical(1);
}

#define TEMPLATE(_CLASS_, _SEXPTYPE_) \
SEXP R_valid_ ## _CLASS_(SEXP obj) \
{ \
	SEXP x = GET_SLOT(obj, Matrix_xSym); \
	if (TYPEOF(x) != _SEXPTYPE_) \
		RMKMS(_("'%s' slot is not of type \"%s\""), \
		      "x", Rf_type2char(_SEXPTYPE_)); \
	return Rf_ScalarLogical(1); \
}
TEMPLATE(nMatrix,  LGLSXP)
TEMPLATE(lMatrix,  LGLSXP)
TEMPLATE(iMatrix,  INTSXP)
TEMPLATE(dMatrix, REALSXP)
TEMPLATE(zMatrix, CPLXSXP)
#undef TEMPLATE

SEXP R_valid_generalMatrix(SEXP obj)
{
	SEXP factors = GET_SLOT(obj, Matrix_factorsSym);
	if (TYPEOF(factors) != VECSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "factors", "list");
	if (XLENGTH(factors) > 0) {
		PROTECT(factors);
		SEXP nms = Rf_getAttrib(factors, R_NamesSymbol);
		UNPROTECT(1); /* factors */
		if (nms == R_NilValue)
			RMKMS(_("'%s' slot has no '%s' attribute"), "factors", "names");
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_symmetricMatrix(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

#ifdef ENFORCE_SYMMETRIC_DIMNAMES

	/* This check can be expensive when both rownames and colnames have
	   nonzero length, and even more so when coercions to character are
	   required ... Users can avoid the expense by setting at least one
	   of rownames and colnames to NULL or by ensuring that they are the
	   same object, as testing for pointer equality is fast ...
	*/

# define ANY_TO_STRING(x) \
	((TYPEOF(x) == INTSXP && Rf_inherits(x, "factor")) \
	 ? Rf_asCharacterFactor(x) \
	 : Rf_coerceVector(x, STRSXP))

	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
		ndn = Rf_getAttrib(dn, R_NamesSymbol);
	UNPROTECT(1); /* dn */

	const char *ndn0, *ndn1;
	if (ndn != R_NilValue &&
	    *(ndn0 = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
	    *(ndn1 = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
	    strcmp(ndn0, ndn1) != 0)
		RMKMS(_("%s[1] differs from %s[2]"), "Dimnames", "Dimnames");
	if (n > 0) {
		/* NB: It is already known that the length of 'dn[[i]]' is 0 or 'n' */
		SEXP rn, cn;
		if ((rn = VECTOR_ELT(dn, 0)) != R_NilValue &&
		    (cn = VECTOR_ELT(dn, 1)) != R_NilValue &&
		    LENGTH(rn) == n && LENGTH(cn) == n && rn != cn) {
			PROTECT(rn);
			PROTECT(cn);
			PROTECT(rn = ANY_TO_STRING(rn));
			PROTECT(cn = ANY_TO_STRING(cn));
			UNPROTECT(4); /* cn, rn */
			if (!equalString(rn, cn, n))
				RMKMS(_("%s[1] differs from %s[2]"), "Dimnames", "Dimnames");
		}
	}

# undef ANY_TO_STRING

#endif

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	if (TYPEOF(uplo) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "uplo", "character");
	if (XLENGTH(uplo) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "uplo", 1);
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "uplo", "U", "L");

	if (HAS_SLOT(obj, Matrix_transSym)) {
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	if (TYPEOF(trans) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "trans", "character");
	if (XLENGTH(trans) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "trans", 1);
	const char *ct = CHAR(STRING_ELT(trans, 0));
	if (ct[0] == '\0' || ct[1] != '\0' || (ct[0] != 'C' && ct[0] != 'T'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "trans", "C", "T");
	}

	return R_valid_generalMatrix(obj);
}

SEXP R_valid_triangularMatrix(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	if (TYPEOF(uplo) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "uplo", "character");
	if (XLENGTH(uplo) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "uplo", 1);
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "uplo", "U", "L");

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	if (TYPEOF(diag) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "diag", "character");
	if (XLENGTH(diag) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "diag", 1);
	const char *nu = CHAR(STRING_ELT(diag, 0));
	if (nu[0] == '\0' || nu[1] != '\0' || (nu[0] != 'N' && nu[0] != 'U'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "diag", "N", "U");

	return Rf_ScalarLogical(1);
}

SEXP R_valid_unpackedMatrix(SEXP obj)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (XLENGTH(x) != (int_fast64_t) m * n)
		RMKMS(_("'%s' slot does not have length %s"), "x", "prod(Dim)");
	return Rf_ScalarLogical(1);
}

SEXP R_valid_packedMatrix(SEXP obj)
{
	int n = DIM(obj)[1];
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (XLENGTH(x) != n + ((int_fast64_t) n * (n - 1)) / 2)
		RMKMS(_("'%s' slot does not have length %s"), "x", "Dim[1]*(Dim[1]+1)/2");
	return Rf_ScalarLogical(1);
}

SEXP R_valid_CsparseMatrix(SEXP obj)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(2); /* i, p */

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[2]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "p");
	int j;
	for (j = 1; j <= n; ++j) {
		if (pp[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[j] < pp[j - 1])
			RMKMS(_("'%s' slot is not nondecreasing"), "p");
		if (pp[j] - pp[j - 1] > m)
			RMKMS(_("first differences of '%s' slot exceed %s"), "p", "Dim[1]");
	}

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (XLENGTH(i) < pp[n])
		RMKMS(_("'%s' slot has length less than %s"), "i", "p[length(p)]");
	int *pi = INTEGER(i), k, kend, ik, i0;
	for (j = 1, k = 0; j <= n; ++j) {
		kend = pp[j];
		i0 = -1;
		while (k < kend) {
			ik = pi[k];
			if (ik == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (ik < 0 || ik >= m)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (ik <= i0)
				RMKMS(_("'%s' slot is not increasing within columns"), "i");
			i0 = ik;
			++k;
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_RsparseMatrix(SEXP obj)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	UNPROTECT(2); /* j, p */

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != m)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[1]+1");
	int *pp = INTEGER(p);
	if (pp[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "p");
	int i;
	for (i = 1; i <= m; ++i) {
		if (pp[i] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[i] < pp[i - 1])
			RMKMS(_("'%s' slot is not nondecreasing"), "p");
		if (pp[i] - pp[i - 1] > n)
			RMKMS(_("first differences of '%s' slot exceed %s"), "p", "Dim[2]");
	}

	if (TYPEOF(j) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "j", "integer");
	if (XLENGTH(j) < pp[m])
		RMKMS(_("'%s' slot has length less than %s"), "j", "p[length(p)]");
	int *pj = INTEGER(j), k, kend, jk, j0;
	for (i = 1, k = 0; i <= m; ++i) {
		kend = pp[i];
		j0 = -1;
		while (k < kend) {
			jk = pj[k];
			if (jk == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "j");
			if (jk < 0 || jk >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "j", "0,...,Dim[2]-1");
			if (jk <= j0)
				RMKMS(_("'%s' slot is not increasing within rows"), "j");
			j0 = jk;
			++k;
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_TsparseMatrix(SEXP obj)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	UNPROTECT(2); /* j, i */

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (TYPEOF(j) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "j", "integer");
	R_xlen_t nnz = XLENGTH(i);
	if (XLENGTH(j) != nnz)
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "j");
	if (nnz > 0) {
		if (m == 0 || n == 0)
			RMKMS(_("'%s' slot has nonzero length but %s is 0"), "i", "prod(Dim)");
		int *pi = INTEGER(i), *pj = INTEGER(j);
		while (nnz--) {
			if (*pi == NA_LOGICAL)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (*pj == NA_LOGICAL)
				RMKMS(_("'%s' slot contains NA"), "j");
			if (*pi < 0 || *pi >= m)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (*pj < 0 || *pj >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "j", "0,...,Dim[2]-1");
			++pi;
			++pj;
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_diagonalMatrix(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	if (TYPEOF(diag) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "diag", "character");
	if (XLENGTH(diag) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "diag", 1);
	const char *nu = CHAR(STRING_ELT(diag, 0));
	if (nu[0] == '\0' || nu[1] != '\0' || (nu[0] != 'N' && nu[0] != 'U'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "diag", "N", "U");
	int nonunit = nu[0] == 'N';

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (XLENGTH(x) != ((nonunit) ? n : 0))
		RMKMS(_("'%s' slot is \"%s\" but '%s' slot does not have length %s"),
		      "diag", (nonunit) ? "N" : "U", "x", (nonunit) ? "Dim[1]" : "0");

	return Rf_ScalarLogical(1);
}

SEXP R_valid_indMatrix(SEXP obj)
{
	SEXP margin = GET_SLOT(obj, Matrix_marginSym);
	if (TYPEOF(margin) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "margin", "integer");
	if (XLENGTH(margin) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "margin", 1);
	int mg = INTEGER(margin)[0] - 1;
	if (mg != 0 && mg != 1)
		RMKMS(_("'%s' slot is not %d or %d"), "margin", 1, 2);

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];
	if (m > 0 && n == 0)
		RMKMS(_("%s-by-%s %s invalid for positive '%s' when %s=%d"),
		      (mg == 0) ? "m" : "0", (mg == 0) ? "0" : "n", "indMatrix",
		      (mg == 0) ? "m" : "n", "margin", (mg == 0) ? 1 : 2);

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != m)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[margin]");
	int *pperm = INTEGER(perm);
	while (m--) {
		if (*pperm == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < 1 || *pperm > n)
			RMKMS(_("'%s' slot has elements not in {%s}"),
			      "perm", "1,...,Dim[1+margin%%2]");
		++pperm;
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_pMatrix(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	if (n > 1) {
		SEXP perm = GET_SLOT(obj, Matrix_permSym);
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int j, *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (work[*pperm - 1])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++) - 1] = 1;
		}
		Matrix_Free(work, lwork);
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sCMatrix(SEXP obj)
{
	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (pp[n] > 0) {
		PROTECT(p);

		char ul = UPLO(obj);

		SEXP i = GET_SLOT(obj, Matrix_iSym);
		int *pi = INTEGER(i), j, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] > j)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					++k;
				}
			}
		} else {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] < j)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					++k;
				}
			}
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_tCMatrix(SEXP obj)
{
	if (DIAG(obj) == 'N')
		return R_valid_sCMatrix(obj);

	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), n = (int) (XLENGTH(p) - 1);
	if (pp[n] > 0) {
		PROTECT(p);

		char ul = UPLO(obj);

		SEXP i = GET_SLOT(obj, Matrix_iSym);
		int *pi = INTEGER(i), j, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] >  j)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					if (pi[k] == j)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		} else {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp[j + 1];
				while (k < kend) {
					if (pi[k] <  j)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					if (pi[k] == j)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sRMatrix(SEXP obj)
{
	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
	if (pp[m] > 0) {
		PROTECT(p);

		char ul = UPLO(obj);

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pj = INTEGER(j), i, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] < i)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					++k;
				}
			}
		} else {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] > i)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					++k;
				}
			}
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_tRMatrix(SEXP obj)
{
	if (DIAG(obj) == 'N')
		return R_valid_sRMatrix(obj);

	SEXP p = GET_SLOT(obj, Matrix_pSym);
	int *pp = INTEGER(p), m = (int) (XLENGTH(p) - 1);
	if (pp[m] > 0) {
		PROTECT(p);

		char ul = UPLO(obj);

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pj = INTEGER(j), i, k, kend;

		UNPROTECT(1); /* p */
		if (ul == 'U') {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] <  i)
						RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
						      "uplo", "U");
					if (pj[k] == i)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		} else {
			for (i = 0, k = 0; i < m; ++i) {
				kend = pp[i + 1];
				while (k < kend) {
					if (pj[k] >  i)
						RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
						      "uplo", "L");
					if (pj[k] == i)
						RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
						      "diag", "U");
					++k;
				}
			}
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sTMatrix(SEXP obj)
{
	SEXP i = GET_SLOT(obj, Matrix_iSym);
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > 0) {
		PROTECT(i);

		char ul = UPLO(obj);

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pi = INTEGER(i), *pj = INTEGER(j);

		UNPROTECT(1); /* i */
		if (ul == 'U') {
			while (nnz--)
				if (*(pi++) > *(pj++))
					RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
					      "uplo", "U");
		} else {
			while (nnz--)
				if (*(pi++) < *(pj++))
					RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
					      "uplo", "L");
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_tTMatrix(SEXP obj)
{
	if (DIAG(obj) == 'N')
		return R_valid_sTMatrix(obj);

	SEXP i = GET_SLOT(obj, Matrix_iSym);
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > 0) {
		PROTECT(i);

		char ul = UPLO(obj);

		SEXP j = GET_SLOT(obj, Matrix_jSym);
		int *pi = INTEGER(i), *pj = INTEGER(j);

		UNPROTECT(1); /* i */
		if (ul == 'U') {
			while (nnz--) {
				if (*pi >  *pj)
					RMKMS(_("%s=\"%s\" but there are entries below the diagonal"),
					      "uplo", "U");
				if (*pi == *pj)
					RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
					      "diag", "U");
				++pi;
				++pj;
			}
		} else {
			while (nnz--) {
				if (*pi <  *pj)
					RMKMS(_("%s=\"%s\" but there are entries above the diagonal"),
					      "uplo", "L");
				if (*pi == *pj)
					RMKMS(_("%s=\"%s\" but there are entries on the diagonal"),
					      "diag", "U");
				++pi;
				++pj;
			}
		}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_xgCMatrix(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(2); /* i, x */
	if (XLENGTH(x) != XLENGTH(i))
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "x");
	return Rf_ScalarLogical(1);
}

SEXP R_valid_xsCMatrix(SEXP obj)
{
	SEXP val = R_valid_xgCMatrix(obj);
	if (TYPEOF(val) != STRSXP)
		val = R_valid_sCMatrix(obj);
	return val;
}

SEXP R_valid_xtCMatrix(SEXP obj)
{
	SEXP val = R_valid_xgCMatrix(obj);
	if (TYPEOF(val) != STRSXP)
		val = R_valid_tCMatrix(obj);
	return val;
}

SEXP R_valid_xgRMatrix(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym));
	UNPROTECT(2); /* j, x */
	if (XLENGTH(x) != XLENGTH(j))
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "j", "x");
	return Rf_ScalarLogical(1);
}

SEXP R_valid_xsRMatrix(SEXP obj)
{
	SEXP val = R_valid_xgRMatrix(obj);
	if (TYPEOF(val) != STRSXP)
		val = R_valid_sRMatrix(obj);
	return val;
}

SEXP R_valid_xtRMatrix(SEXP obj)
{
	SEXP val = R_valid_xgRMatrix(obj);
	if (TYPEOF(val) != STRSXP)
		val = R_valid_tRMatrix(obj);
	return val;
}

SEXP R_valid_xgTMatrix(SEXP obj)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym));
	UNPROTECT(2); /* i, x */
	if (XLENGTH(x) != XLENGTH(i))
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "i", "x");
	return Rf_ScalarLogical(1);
}

SEXP R_valid_xsTMatrix(SEXP obj)
{
	SEXP val = R_valid_xgTMatrix(obj);
	if (TYPEOF(val) != STRSXP)
		val = R_valid_sTMatrix(obj);
	return val;
}

SEXP R_valid_xtTMatrix(SEXP obj)
{
	SEXP val = R_valid_xgTMatrix(obj);
	if (TYPEOF(val) != STRSXP)
		val = R_valid_tTMatrix(obj);
	return val;
}

/* NB: Non-finite entries are "valid" because we consider
   crossprod(x) and tcrossprod(x) to be positive semidefinite
   even if 'x' contains non-finite entries (for speed) ...
*/

SEXP R_valid_xpoMatrix(SEXP obj)
{
	int n = DIM(obj)[1], j;
	R_xlen_t n1a = (R_xlen_t) n + 1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += n1a)
		if (!ISNAN(*px) && *px < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	} else {
	if (TRANS(obj) != 'C')
		RMKMS(_("'%s' slot is not \"%s\""), "trans", "C");
	Rcomplex *px = COMPLEX(x);
	for (j = 0; j < n; ++j, px += n1a)
		if (!ISNAN((*px).r) && (*px).r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_xppMatrix(SEXP obj)
{
	int n = DIM(obj)[1], j;
	char ul = UPLO(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	if (ul == 'U') {
	for (j = 0; j < n; px += (++j)+1)
		if (!ISNAN(*px) && *px < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	} else {
	for (j = 0; j < n; px += n-(j++))
		if (!ISNAN(*px) && *px < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	}
	} else {
	if (TRANS(obj) != 'C')
		RMKMS(_("'%s' slot is not \"%s\""), "trans", "C");
	Rcomplex *px = COMPLEX(x);
	if (ul == 'U') {
	for (j = 0; j < n; px += (++j)+1)
		if (!ISNAN((*px).r) && (*px).r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	} else {
	for (j = 0; j < n; px += n-(j++))
		if (!ISNAN((*px).r) && (*px).r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_xpCMatrix(SEXP obj)
{
	int n = DIM(obj)[1], j;
	char ul = UPLO(obj);

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	UNPROTECT(3); /* x, i, p */
	int *pp = INTEGER(p), *pi = INTEGER(i);
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	if (ul == 'U') {
	for (j = 0; j < n; ++j)
		if (pp[j + 1] - pp[j] > 0 && pi[pp[j + 1] - 1] == j &&
		    !ISNAN(px[pp[j + 1] - 1]) && px[pp[j + 1] - 1] < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	} else {
	for (j = 0; j < n; ++j)
		if (pp[j + 1] - pp[j] > 0 && pi[pp[j]] == j &&
		    !ISNAN(px[pp[j]]) && px[pp[j]] < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	}
	} else {
	if (TRANS(obj) != 'C')
		RMKMS(_("'%s' slot is not \"%s\""), "trans", "C");
	Rcomplex *px = COMPLEX(x);
	if (ul == 'U') {
	for (j = 0; j < n; ++j)
		if (pp[j + 1] - pp[j] > 0 && pi[pp[j + 1] - 1] == j &&
		    !ISNAN(px[pp[j + 1] - 1].r) && px[pp[j + 1] - 1].r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	} else {
	for (j = 0; j < n; ++j)
		if (pp[j + 1] - pp[j] > 0 && pi[pp[j]] == j &&
		    !ISNAN(px[pp[j]].r) && px[pp[j]].r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_xpRMatrix(SEXP obj)
{
	int m = DIM(obj)[0], i;
	char ul = UPLO(obj);

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	UNPROTECT(3); /* x, j, p */
	int *pp = INTEGER(p), *pj = INTEGER(j);
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	if (ul == 'U') {
	for (i = 0; i < m; ++i)
		if (pp[i + 1] - pp[i] > 0 && pj[pp[i]] == i &&
		    !ISNAN(px[pp[i]]) && px[pp[i]] < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	} else {
	for (i = 0; i < m; ++i)
		if (pp[i + 1] - pp[i] > 0 && pj[pp[i + 1] - 1] == i &&
		    !ISNAN(px[pp[i + 1] - 1]) && px[pp[i + 1] - 1] < 0.0)
			RMK(_("matrix has negative diagonal elements"));
	}
	} else {
	if (TRANS(obj) != 'C')
		RMKMS(_("'%s' slot is not \"%s\""), "trans", "C");
	Rcomplex *px = COMPLEX(x);
	if (ul == 'U') {
	for (i = 0; i < m; ++i)
		if (pp[i + 1] - pp[i] > 0 && pj[pp[i]] == i &&
		    !ISNAN(px[pp[i]].r) && px[pp[i]].r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	} else {
	for (i = 0; i < m; ++i)
		if (pp[i + 1] - pp[i] > 0 && pj[pp[i + 1] - 1] == i &&
		    !ISNAN(px[pp[i + 1] - 1].r) && px[pp[i + 1] - 1].r < 0.0)
			RMK(_("matrix has diagonal elements with negative real part"));
	}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_xpTMatrix(SEXP obj)
{
	int n = DIM(obj)[1];

	SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		j = PROTECT(GET_SLOT(obj, Matrix_jSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));
	UNPROTECT(3); /* x, j, i */
	int *pi = INTEGER(i), *pj = INTEGER(j);
	R_xlen_t nnz = XLENGTH(x);

	double *work;
	int lwork = n;
	Matrix_Calloc(work, lwork, double);
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	while (nnz--) {
		if (*pi == *pj)
			work[*pi] += *px;
		++pi; ++pj; ++px;
	}
	while (n--)
		if (!ISNAN(work[n]) && work[n] < 0.0)
			FRMK(_("matrix has negative diagonal elements"));
	} else {
	if (TRANS(obj) != 'C')
		RMKMS(_("'%s' slot is not \"%s\""), "trans", "C");
	Rcomplex *px = COMPLEX(x);
	while (nnz--) {
		if (*pi == *pj)
			work[*pi] += (*px).r;
		++pi; ++pj; ++px;
	}
	while (n--)
		if (!ISNAN(work[n]) && work[n] < 0.0)
			FRMK(_("matrix has diagonal elements with negative real part"));
	}
	Matrix_Free(work, lwork);

	return Rf_ScalarLogical(1);
}

SEXP R_valid_corMatrix(SEXP obj)
{
	int n = DIM(obj)[1], j;
	R_xlen_t n1a = (R_xlen_t) n + 1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	for (j = 0; j < n; ++j, px += n1a)
		if (ISNAN(*px) || *px != 1.0)
			RMK(_("matrix has nonunit diagonal elements"));

	SEXP sd = GET_SLOT(obj, Matrix_sdSym);
	if (TYPEOF(sd) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "sd", "double");
	if (XLENGTH(sd) != n)
		RMKMS(_("'%s' slot does not have length %s"), "sd", "Dim[1]");
	double *psd = REAL(sd);
	for (j = 0; j < n; ++j)
		if (!ISNAN(psd[j]) && psd[j] < 0.0)
			RMKMS(_("'%s' slot has negative elements"), "sd");

	return Rf_ScalarLogical(1);
}

SEXP R_valid_copMatrix(SEXP obj)
{
	int n = DIM(obj)[1], j;
	char ul = UPLO(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	double *px = REAL(x);
	if (ul == 'U') {
	for (j = 0; j < n; px += (++j)+1)
		if (ISNAN(*px) || *px != 1.0)
			RMK(_("matrix has nonunit diagonal elements"));
	} else {
	for (j = 0; j < n; px += n-(j++))
		if (ISNAN(*px) || *px != 1.0)
			RMK(_("matrix has nonunit diagonal elements"));
	}

	SEXP sd = GET_SLOT(obj, Matrix_sdSym);
	if (TYPEOF(sd) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "sd", "double");
	if (XLENGTH(sd) != n)
		RMKMS(_("'%s' slot does not have length %s"), "sd", "Dim[1]");
	double *psd = REAL(sd);
	for (j = 0; j < n; ++j)
		if (!ISNAN(psd[j]) && psd[j] < 0.0)
			RMKMS(_("'%s' slot has negative elements"), "sd");

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sparseVector(SEXP obj)
{
	SEXP length = GET_SLOT(obj, Matrix_lengthSym);
	if (TYPEOF(length) != INTSXP && TYPEOF(length) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "length", "integer", "double");
	if (XLENGTH(length) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "length", 1);
	int_fast64_t n;
	if (TYPEOF(length) == INTSXP) {
		int n_ = INTEGER(length)[0];
		if (n_ == NA_INTEGER)
			RMKMS(_("'%s' slot is NA"), "length");
		if (n_ < 0)
			RMKMS(_("'%s' slot is negative"), "length");
		n = (int_fast64_t) n_;
	} else {
		double n_ = REAL(length)[0];
		if (ISNAN(n_))
			RMKMS(_("'%s' slot is NA"), "length");
		if (n_ < 0.0)
			RMKMS(_("'%s' slot is negative"), "length");
		if (n_ > 0x1.0p+53)
			RMKMS(_("'%s' slot exceeds %s"), "length", "2^53");
		n = (int_fast64_t) n_;
	}

	SEXP i = GET_SLOT(obj, Matrix_iSym);
	if (TYPEOF(i) != INTSXP && TYPEOF(i) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "i", "integer", "double");
	R_xlen_t nnz = XLENGTH(i);
	if (nnz > n)
		RMKMS(_("'%s' slot has length greater than '%s' slot"), "i", "length");
	if (TYPEOF(i) == INTSXP) {
		int *pi = INTEGER(i), max = (n > INT_MAX) ? INT_MAX : (int) n, last = 0;
		while (nnz--) {
			if (*pi == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (*pi < 1 || *pi > max)
				RMKMS(_("'%s' slot has elements not in {%s}"),
					  "i", "1,...,length");
			if (*pi <= last)
				RMKMS(_("'%s' slot is not increasing"), "i");
			last = *(pi++);
		}
	} else {
		double *pi = REAL(i), max = (double) n, last = 0.0, tmp;
		while (nnz--) {
			if (ISNAN(*pi))
				RMKMS(_("'%s' slot contains NA"), "i");
			tmp = trunc(*(pi++));
			if (tmp < 1.0 || tmp > max)
				RMKMS(_("'%s' slot has elements not in {%s} after truncation towards zero"),
					  "i", "1,...,length");
			if (tmp <= last)
				RMKMS(_("'%s' slot is not increasing after truncation towards zero"), "i");
			last = tmp;
		}
	}

	return Rf_ScalarLogical(1);
}

#define TEMPLATE(_CLASS_, _SEXPTYPE_) \
SEXP R_valid_ ## _CLASS_(SEXP obj) \
{ \
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)); \
	UNPROTECT(2); /* i, x */ \
	if (TYPEOF(x) != _SEXPTYPE_) \
		RMKMS(_("'%s' slot is not of type \"%s\""), \
		      "x", Rf_type2char(_SEXPTYPE_)); \
	if (XLENGTH(x) != XLENGTH(i)) \
		RMKMS(_("'%s' and '%s' slots do not have equal length"), \
		      "i", "x"); \
	return Rf_ScalarLogical(1); \
}
TEMPLATE(lsparseVector,  LGLSXP)
TEMPLATE(isparseVector,  INTSXP)
TEMPLATE(dsparseVector, REALSXP)
TEMPLATE(zsparseVector, CPLXSXP)
#undef TEMPLATE

SEXP R_valid_MatrixFactorization(SEXP obj)
{
	return R_valid_Matrix(obj);
}

SEXP R_valid_denseSchur(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");
	int_fast64_t nn = (int_fast64_t) n * n;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	if (XLENGTH(x) != nn && XLENGTH(x) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "x", "prod(Dim)", "0");
	int normal = n > 0 && XLENGTH(x) == 0;

	SEXP vectors = GET_SLOT(obj, Matrix_vectorsSym);
	if (TYPEOF(vectors) != TYPEOF(x))
		RMKMS(_("'%s' and '%s' slots do not have the same type"),
		      "x", "vectors");
	if (XLENGTH(vectors) != nn && XLENGTH(vectors) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "vectors", "prod(Dim)", "0");

	SEXP values = GET_SLOT(obj, Matrix_valuesSym);
	if (TYPEOF(values) != REALSXP) {
		if (normal)
		RMKMS(_("'%s' slot is not of type \"%s\""),
		      "values", "double");
		else if (TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "values", "double", "complex");
	}
	if (XLENGTH(values) != n)
		RMKMS(_("'%s' slot does not have length %s"), "values", "Dim[1]");

	return Rf_ScalarLogical(1);
}

SEXP R_valid_denseQR(SEXP obj)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	if (XLENGTH(x) != (int_fast64_t) m * n)
		RMKMS(_("'%s' slot does not have length %s"), "x", "prod(Dim)");

	SEXP beta = GET_SLOT(obj, Matrix_betaSym);
	if (TYPEOF(beta) != REALSXP && TYPEOF(beta) != CPLXSXP)
		RMKMS(_("'%s' and '%s' slots do not have the same type"),
		      "x", "beta");
	if (XLENGTH(beta) != r)
		RMKMS(_("'%s' slot does not have length %s"), "beta", "min(Dim)");

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != n)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[2]");
	char *work;
	int lwork = n;
	Matrix_Calloc(work, lwork, char);
	int j, *pperm = INTEGER(perm);
	for (j = 0; j < n; ++j) {
		if (*pperm == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < 1 || *pperm > n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			      "perm", "1,...,Dim[2]");
		if (work[*pperm - 1])
			FRMKMS(_("'%s' slot contains duplicates"), "perm");
		work[*(pperm++)] = 1;
	}
	Matrix_Free(work, lwork);

	return Rf_ScalarLogical(1);
}

SEXP R_valid_denseLU(SEXP obj)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	if (XLENGTH(x) != (int_fast64_t) m * n)
		RMKMS(_("'%s' slot does not have length %s"), "x", "prod(Dim)");

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != r)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "min(Dim)");
	int *pperm = INTEGER(perm);
	while (r--) {
		if (*pperm == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < 1 || *pperm > m)
			RMKMS(_("'%s' slot has elements not in {%s}"),
			      "perm", "1,...,Dim[1]");
		++pperm;
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_denseBunchKaufman(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	if (TYPEOF(uplo) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "uplo", "character");
	if (XLENGTH(uplo) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "uplo", 1);
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "uplo", "U", "L");

	if (HAS_SLOT(obj, Matrix_transSym)) {
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	if (TYPEOF(trans) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "trans", "character");
	if (XLENGTH(trans) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "trans", 1);
	const char *ct = CHAR(STRING_ELT(trans, 0));
	if (ct[0] == '\0' || ct[1] != '\0' || (ct[0] != 'C' && ct[0] != 'T'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "trans", "C", "T");
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	int packed = XLENGTH(x) != (int_fast64_t) n * n;
	if (packed && XLENGTH(x) != n + ((int_fast64_t) n * (n - 1)) / 2)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "x", "prod(Dim)", "Dim[1]*(Dim[1]+1)/2");

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != n)
		RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[1]");
	int n_ = n, *pperm = INTEGER(perm);
	while (n_) {
		if (*pperm == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "perm");
		if (*pperm < -n || *pperm == 0 || *pperm > n)
			RMKMS(_("'%s' slot has elements not in {%s}\\{%s}"),
			      "perm", "-Dim[1],...,Dim[1]", "0");
		if (*pperm > 0) {
			pperm += 1;
			n_ -= 1;
		} else if (n_ > 1 && *(pperm + 1) == *pperm) {
			pperm += 2;
			n_ -= 2;
		} else
			RMKMS(_("'%s' slot has unpaired negative elements"), "perm");
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_denseCholesky(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
	if (TYPEOF(uplo) != STRSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "uplo", "character");
	if (XLENGTH(uplo) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "uplo", 1);
	const char *ul = CHAR(STRING_ELT(uplo, 0));
	if (ul[0] == '\0' || ul[1] != '\0' || (ul[0] != 'U' && ul[0] != 'L'))
		RMKMS(_("'%s' slot is not \"%s\" or \"%s\""), "uplo", "U", "L");

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	int j, packed = XLENGTH(x) != (int_fast64_t) n * n;
	if (packed && XLENGTH(x) != n + ((int_fast64_t) n * (n - 1)) / 2)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "x", "prod(Dim)", "Dim[1]*(Dim[1]+1)/2");

	/* Real, non-negative diagonal elements are necessary and sufficient */
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	if (!packed) {
		R_xlen_t n1a = (R_xlen_t) n + 1;
		for (j = 0; j < n; ++j, px += n1a)
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("Cholesky factor has negative diagonal elements"));
	} else if (*ul == 'U') {
		for (j = 0; j < n; ++j, px += (++j)+1)
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("Cholesky factor has negative diagonal elements"));
	} else {
		for (j = 0; j < n; ++j, px += n-(j++))
			if (!ISNAN(*px) && *px < 0.0)
				RMK(_("Cholesky factor has negative diagonal elements"));
	}
	} else {
	Rcomplex *px = COMPLEX(x);
	if (!packed) {
		R_xlen_t n1a = (R_xlen_t) n + 1;
		for (j = 0; j < n; ++j, px += n1a)
			if (!ISNAN((*px).r) && (*px).r < 0.0)
				RMK(_("Cholesky factor has diagonal elements with negative real part"));
	} else if (*ul == 'U') {
		for (j = 0; j < n; ++j, px += (++j)+1)
			if (!ISNAN((*px).r) && (*px).r < 0.0)
				RMK(_("Cholesky factor has diagonal elements with negative real part"));
	} else {
		for (j = 0; j < n; ++j, px += n-(j++))
			if (!ISNAN((*px).r) && (*px).r < 0.0)
				RMK(_("Cholesky factor has diagonal elements with negative real part"));
	}
	}

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (XLENGTH(perm) != n && XLENGTH(perm) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "perm", "Dim[1]", "0");
	if (LENGTH(perm) == n) {
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (*pperm == NA_INTEGER)
				FRMKMS(_("'%s' slot contains NA"), "perm");
			if (*pperm < 0 || *pperm >= n)
				FRMKMS(_("'%s' slot has elements not in {%s}"),
				       "perm", "0,...,Dim[1]-1");
			if (work[*pperm])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, lwork);
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sparseQR(SEXP obj)
{
	/* MJ: assuming for simplicity that 'V' and 'R' slots are formally valid */

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];
	if (m < n)
		RMK(_("matrix has more columns than rows"));

	SEXP beta = GET_SLOT(obj, Matrix_betaSym);
	if (TYPEOF(beta) != REALSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "beta", "double");
	if (XLENGTH(beta) != n)
		RMKMS(_("'%s' slot does not have length %s"), "beta", "Dim[2]");

	SEXP p, i, q;
	int *pp, *pi, *pq, j, k, kend, m0;

	SEXP V = PROTECT(GET_SLOT(obj, Matrix_VSym));
	PROTECT(p = GET_SLOT(V, Matrix_pSym));
	PROTECT(i = GET_SLOT(V, Matrix_iSym));
	pdim = DIM(V); m0 = pdim[0];
	UNPROTECT(3); /* i, p, V */
	if (m0 < m)
		RMKMS(_("'%s' slot has fewer than %s rows"), "V", "Dim[1]");
	if (m0 > m + n)
		RMKMS(_("'%s' slot has more than %s rows"), "V", "Dim[1]+Dim[2]");
	if (pdim[1] != n)
		RMKMS(_("'%s' slot does not have %s columns"), "V", "Dim[2]");
	pp = INTEGER(p);
	pi = INTEGER(i);
	for (j = 0, k = 0; j < n; ++j) {
		kend = pp[j + 1];
		if (k < kend) {
			if (pi[k] < j)
				RMKMS(_("'%s' slot must be lower trapezoidal but has entries above the diagonal"), "V");
		}
		k = kend;
	}

	SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
	PROTECT(p = GET_SLOT(R, Matrix_pSym));
	PROTECT(i = GET_SLOT(R, Matrix_iSym));
	pdim = DIM(R);
	UNPROTECT(3); /* i, p, R */
	if (pdim[0] != m0)
		RMKMS(_("'%s' slot does not have %s row"), "R", "nrow(V)");
	if (pdim[1] != n)
		RMKMS(_("'%s' slot does not have %s columns"), "R", "Dim[2]");
	pp = INTEGER(p);
	pi = INTEGER(i);
	for (j = 0, k = 0; j < n; ++j) {
		kend = pp[j + 1];
		if (k < kend && pi[kend - 1] > j)
			RMKMS(_("'%s' slot must be upper trapezoidal but has entries below the diagonal"), "R");
		k = kend;
	}

	PROTECT(p = GET_SLOT(obj, Matrix_pSym));
	PROTECT(q = GET_SLOT(obj, Matrix_qSym));
	UNPROTECT(2); /* q, p */
	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (TYPEOF(q) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "q", "integer");
	if (XLENGTH(p) != m0)
		RMKMS(_("'%s' slot does not have length %s"), "p", "nrow(V)");
	if (XLENGTH(q) != n && XLENGTH(q) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "q", "Dim[2]", "0");
	char *work;
	int lwork = m0; /* n <= m <= m0 */
	Matrix_Calloc(work, lwork, char);
	pp = INTEGER(p);
	for (j = 0; j < m0; ++j) {
		if (*pp == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "p");
		if (*pp < 0 || *pp >= m0)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "p", "0,...,nrow(V)-1");
		if (work[*pp])
			FRMKMS(_("'%s' slot contains duplicates"), "p");
		work[*(pp++)] = 1;
	}
	if (LENGTH(q) == n) {
	pq = INTEGER(q);
	for (j = 0; j < n; ++j) {
		if (*pq == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "q");
		if (*pq < 0 || *pq >= n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "q", "0,...,Dim[2]-1");
		if (!work[*pq])
			FRMKMS(_("'%s' slot contains duplicates"), "q");
		work[*(pq++)] = 0;
	}
	}
	Matrix_Free(work, lwork);

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sparseLU(SEXP obj)
{
	/* MJ: assuming for simplicity that 'L' and 'U' slots are formally valid */

	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");
	char ul, nu;

	SEXP L = PROTECT(GET_SLOT(obj, Matrix_LSym));
	ul = UPLO(L);
	nu = DIAG(L);
	pdim = DIM(L);
	UNPROTECT(1); /* L */
	if (pdim[0] != n || pdim[1] != n)
		RMKMS(_("dimensions of '%s' slot are not identical to '%s'"), "L", "Dim");
	if (ul == 'U')
		RMKMS(_("'%s' slot is upper (not lower) triangular"), "L");
	if (nu == 'N') {
		PROTECT(L);
		SEXP p = PROTECT(GET_SLOT(L, Matrix_pSym)),
			i = PROTECT(GET_SLOT(L, Matrix_iSym)),
			x = PROTECT(GET_SLOT(L, Matrix_xSym));
		UNPROTECT(4); /* x, i, p, L */

		int *pp = INTEGER(p), *pi = INTEGER(i), j, k, kend;
		if (TYPEOF(x) == REALSXP) {
		double *px = REAL(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j + 1];
			if (kend == k || pi[k] != j || px[k] != 1.0)
				RMKMS(_("'%s' slot has nonunit diagonal elements"), "L");
			k = kend;
		}
		} else {
		Rcomplex *px = COMPLEX(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j + 1];
			if (kend == k || pi[k] != j || px[k].r != 1.0 || px[k].i != 0.0)
				RMKMS(_("'%s' slot has nonunit diagonal elements"), "L");
			k = kend;
		}
		}
	}

	SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym));
	ul = UPLO(U);
	pdim = DIM(U);
	UNPROTECT(1); /* U */
	if (pdim[0] != n || pdim[1] != n)
		RMKMS(_("dimensions of '%s' slot are not identical to '%s'"), "U", "Dim");
	if (ul != 'U')
		RMKMS(_("'%s' slot is lower (not upper) triangular"), "U");

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		q = PROTECT(GET_SLOT(obj, Matrix_qSym));
	UNPROTECT(2); /* q, p */
	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (TYPEOF(q) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "q", "integer");
	if (XLENGTH(p) != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[1]");
	if (XLENGTH(q) != n && XLENGTH(q) != 0)
		RMKMS(_("'%s' slot does not have length %s or length %s"),
		      "q", "Dim[2]", "0");
	char *work;
	int lwork = n;
	Matrix_Calloc(work, lwork, char);
	int j, *pp = INTEGER(p);
	for (j = 0; j < n; ++j) {
		if (*pp == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "p");
		if (*pp < 0 || *pp >= n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "p", "0,...,Dim[1]-1");
		if (work[*pp])
			FRMKMS(_("'%s' slot contains duplicates"), "p");
		work[*(pp++)] = 1;
	}
	if (LENGTH(q) == n) {
	int *pq = INTEGER(q);
	for (j = 0; j < n; ++j) {
		if (*pq == NA_INTEGER)
			FRMKMS(_("'%s' slot contains NA"), "q");
		if (*pq < 0 || *pq >= n)
			FRMKMS(_("'%s' slot has elements not in {%s}"),
			       "q", "0,...,Dim[2]-1");
		if (!work[*pq])
			FRMKMS(_("'%s' slot contains duplicates"), "q");
		work[*(pq++)] = 0;
	}
	}
	Matrix_Free(work, lwork);

	return Rf_ScalarLogical(1);
}

SEXP R_valid_sparseCholesky(SEXP obj)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		RMKMS(_("%s[1] != %s[2] (matrix is not square)"), "Dim", "Dim");

	SEXP ordering = GET_SLOT(obj, Matrix_orderingSym);
	if (TYPEOF(ordering) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "ordering", "integer");
	if (XLENGTH(ordering) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "ordering", 1);
	int or = INTEGER(ordering)[0];
	if (or < 0 || or > 6)
		RMKMS(_("'%s' is not in %d:%d"), "ordering", 0, 6);

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	if (TYPEOF(perm) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "perm", "integer");
	if (or == 0) {
		if (XLENGTH(perm) != 0)
			RMKMS(_("'%s' slot does not have length %d"), "perm", 0);
	} else {
		if (XLENGTH(perm) != n)
			RMKMS(_("'%s' slot does not have length %s"), "perm", "Dim[1]");
		char *work;
		int lwork = n;
		Matrix_Calloc(work, lwork, char);
		int j, *pperm = INTEGER(perm);
		for (j = 0; j < n; ++j) {
			if (*pperm == NA_INTEGER)
				FRMKMS(_("'%s' slot contains NA"), "perm");
			if (*pperm < 0 || *pperm >= n)
				FRMKMS(_("'%s' slot has elements not in {%s}"),
				       "perm", "0,...,Dim[1]-1");
			if (work[*pperm])
				FRMKMS(_("'%s' slot contains duplicates"), "perm");
			work[*(pperm++)] = 1;
		}
		Matrix_Free(work, lwork);
	}

	SEXP colcount = GET_SLOT(obj, Matrix_colcountSym);
	if (TYPEOF(colcount) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "colcount", "integer");
	if (XLENGTH(colcount) != n)
		RMKMS(_("'%s' slot does not have length %s"), "colcount", "Dim[2]");
	int j, *pcolcount = INTEGER(colcount);
	for (j = 0; j < n; ++j) {
		if (pcolcount[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "colcount");
		if (pcolcount[j] < 0 || pcolcount[j] > n - j)
			RMKMS(_("%s is not in {%s}"), "colcount[j]", "0,...,Dim[2]-j+1");
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_simplicialCholesky(SEXP obj)
{
	int pattern = !HAS_SLOT(obj, Matrix_minorSym);
	if (pattern)
		return Rf_ScalarLogical(1);

	int n = DIM(obj)[1];
	if (n == INT_MAX)
		RMKMS(_("%s is not representable as \"%s\""), "Dim[2]+1", "integer");

	SEXP minor = GET_SLOT(obj, Matrix_minorSym);
	int mr = INTEGER(minor)[0];
	if (mr < 0 || mr > n)
		RMKMS(_("'%s' slot is not in {%s}"), "minor", "{0,...,Dim[2]}");

	SEXP is_ll = GET_SLOT(obj, Matrix_isllSym);
	if (TYPEOF(is_ll) != LGLSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "is_ll", "logical");
	if (XLENGTH(is_ll) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "is_ll", 1);
	int ll = LOGICAL(is_ll)[0];
	if (ll == NA_LOGICAL)
		RMKMS(_("'%s' slot is NA"), "is_ll");

	SEXP is_monotonic = GET_SLOT(obj, Matrix_ismtSym);
	if (TYPEOF(is_monotonic) != LGLSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "is_monotonic", "logical");
	if (XLENGTH(is_monotonic) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "is_monotonic", 1);
	int mt = LOGICAL(is_monotonic)[0];
	if (mt == NA_LOGICAL)
		RMKMS(_("'%s' slot is NA"), "is_monotonic");

	SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		nz = PROTECT(GET_SLOT(obj, Matrix_nzSym)),
		next = PROTECT(GET_SLOT(obj, Matrix_nextSym)),
		prev = PROTECT(GET_SLOT(obj, Matrix_prevSym));
	UNPROTECT(6); /* prev, next, nz, x, i, p */

	if (TYPEOF(next) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "next", "integer");
	if (TYPEOF(prev) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "prev", "integer");
	if (XLENGTH(next) - 2 != n)
		RMKMS(_("'%s' slot does not have length %s"), "next", "Dim[2]+2");
	if (XLENGTH(prev) - 2 != n)
		RMKMS(_("'%s' slot does not have length %s"), "prev", "Dim[2]+2");
	int *pnext = INTEGER(next), *pprev = INTEGER(prev),
		j1 = pnext[n + 1], j2 = pprev[n], count = n + 1;
	while (count--) {
		if (j1 < 0 || j1 > n)
			RMKMS(_("%s has elements not in {%s}"),
			      "`next`[-(Dim[2]+1)]", "0,...,Dim[2]");
		if (j2 < 0 || j2 > n + 1 || j2 == n)
			RMKMS(_("%s has elements not in {%s}\\{%s}"),
			      "`prev`[-(Dim[2]+2)]", "0,...,Dim[2]+1", "Dim[2]");
		if ((count >  1) && mt && (pnext[j1] != j1 + 1 || pprev[j2] != j2 - 1))
			RMKMS(_("'%s' slot is %s but columns are not stored in increasing order"),
			      "is_monotonic", "TRUE");
		if ((count >= 1) ? j1 == n : j1 != n)
			RMKMS(_("traversal of '%s' slot does not complete in exactly %s steps"),
			      "next", "length(`next`)");
		if ((count >= 1) ? j2 == n + 1 : j2 != n + 1)
			RMKMS(_("traversal of '%s' slot does not complete in exactly %s steps"),
			      "prev", "length(`prev`)");
		j1 = pnext[j1];
		j2 = pprev[j2];
	}
	if (j1 != -1)
		RMKMS(_("%s is not %d"), "`next`[Dim[2]+1]", -1);
	if (j2 != -1)
		RMKMS(_("%s is not %d"), "`prev`[Dim[2]+2]", -1);

	if (TYPEOF(nz) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "nz", "integer");
	if (XLENGTH(nz) != n)
		RMKMS(_("'%s' slot does not have length %s"), "nz", "Dim[2]");
	int j, *pnz = INTEGER(nz);
	for (j = 0; j < n; ++j) {
		if (pnz[j] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "nz");
		if (pnz[j] < 1 || pnz[j] > n - j)
			RMKMS(_("%s is not in {%s}"), "nz[j]", "1,...,Dim[1]-j+1");
	}

	if (TYPEOF(p) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "p", "integer");
	if (XLENGTH(p) - 1 != n)
		RMKMS(_("'%s' slot does not have length %s"), "p", "Dim[2]+1");
	j1 = pnext[n + 1];
	int *pp = INTEGER(p);
	if (pp[j1] != 0)
		RMKMS(_("column '%s' is stored first but %s is not 0"), "j", "p[j]");
	for (j = 0; j < n; ++j) {
		j2 = pnext[j1];
		if (pp[j2] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "p");
		if (pp[j2] < pp[j1])
			RMKMS(_("'%s' slot is not increasing when traversed in stored column order"), "p");
		if (pp[j2] - pp[j1] < pnz[j1])
			RMKMS(_("'%s' slot allocates fewer than %s elements for column '%s'"),
			      "i", "nz[j]", "j");
		if (pp[j2] - pp[j1] > n - j1)
			RMKMS(_("'%s' slot allocates more than %s elements for column '%s'"),
			      "i", "Dim[2]-j+1", "j");
		j1 = j2;
	}

	if (TYPEOF(i) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "i", "integer");
	if (XLENGTH(i) != pp[n])
		RMKMS(_("'%s' slot does not have length %s"), "i", "p[length(p)]");
	int *pi = INTEGER(i), *pi_, k;
	j1 = pnext[n + 1];
	for (j = 0; j < n; ++j) {
		pi_ = pi + pp[j1];
		if (pi_[0] != j1)
			RMKMS(_("first entry in column '%s' does not have row index '%s'"),
			      "j", "j");
		for (k = 1; k < pnz[j1]; ++k) {
			if (pi_[k] == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "i");
			if (pi_[k] < 0 || pi_[k] >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "i", "0,...,Dim[1]-1");
			if (pi_[k] <= pi_[k - 1])
				RMKMS(_("'%s' slot is not increasing within columns"), "i");
		}
		j1 = pnext[j1];
	}

	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	if (XLENGTH(x) != pp[n])
		RMKMS(_("'%s' slot does not have length %s"), "x", "p[length(p)]");

	if (ll) {
	/* Real, non-negative diagonal elements are necessary and sufficient */
	if (TYPEOF(x) == REALSXP) {
	double *px = REAL(x);
	for (j = 0; j < n; ++j)
		if (!ISNAN(px[pp[j]]) && px[pp[j]] < 0.0)
			RMK(_("Cholesky factor has negative diagonal elements"));
	} else {
	Rcomplex *px = COMPLEX(x);
	for (j = 0; j < n; ++j)
		if (!ISNAN(px[pp[j]].r) && px[pp[j]].r < 0.0)
			RMK(_("Cholesky factor has diagonal elements with negative real part"));
	}
	}

	return Rf_ScalarLogical(1);
}

SEXP R_valid_supernodalCholesky(SEXP obj)
{
	int pattern = !HAS_SLOT(obj, Matrix_minorSym), n = DIM(obj)[1];

	if (!pattern) {
	SEXP minor = GET_SLOT(obj, Matrix_minorSym);
	int mr = INTEGER(minor)[0];
	if (mr < 0 || mr > n)
		RMKMS(_("'%s' slot is not in {%s}"), "minor", "{0,...,Dim[2]}");
	}

	SEXP maxcsize = PROTECT(GET_SLOT(obj, Matrix_maxcsizeSym)),
		maxesize = PROTECT(GET_SLOT(obj, Matrix_maxesizeSym)),
		super = PROTECT(GET_SLOT(obj, Matrix_superSym)),
		pi = PROTECT(GET_SLOT(obj, Matrix_piSym)),
		px = PROTECT(GET_SLOT(obj, Matrix_pxSym)),
		s = PROTECT(GET_SLOT(obj, Matrix_sSym)),
		x = (HAS_SLOT(obj, Matrix_xSym)) ? GET_SLOT(obj, Matrix_xSym) : NULL;
	UNPROTECT(6); /* s, px, pi, super, maxesize, maxcsize */

	if (TYPEOF(super) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "super", "integer");
	R_xlen_t nsuper1a = XLENGTH(super);
	if (nsuper1a - 1 < ((n > 0) ? 1 : 0))
		RMKMS(_("'%s' slot has length less than %d"), "super", 2);
	if (nsuper1a - 1 > n)
		RMKMS(_("'%s' slot has length greater than %s"), "super", "Dim[2]+1");
	int k, nsuper = (int) (nsuper1a - 1), *psuper = INTEGER(super);
	if (psuper[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "super");
	if (psuper[nsuper] != n)
		RMKMS(_("last element of '%s' slot is not %s"), "super", "Dim[2]");
	for (k = 1; k <= nsuper; ++k) {
		if (psuper[k] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "super");
		if (psuper[k] <= psuper[k - 1])
			RMKMS(_("'%s' slot is not increasing"), "super");
	}

	if (TYPEOF(pi) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "pi", "integer");
	if (TYPEOF(px) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "px", "integer");
	if (XLENGTH(pi) != nsuper1a)
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "pi", "super");
	if (XLENGTH(px) != nsuper1a)
		RMKMS(_("'%s' and '%s' slots do not have equal length"), "px", "super");
	int *ppi = INTEGER(pi), *ppx = INTEGER(px), nr, nc, l, ml = 0;
	if (ppi[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "pi");
	if (ppx[0] != 0)
		RMKMS(_("first element of '%s' slot is not 0"), "px");
	for (k = 1; k <= nsuper; ++k) {
		if (ppi[k] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "pi");
		if (ppx[k] == NA_INTEGER)
			RMKMS(_("'%s' slot contains NA"), "px");
		if (ppi[k] <= ppi[k - 1])
			RMKMS(_("'%s' slot is not increasing"), "pi");
		if (ppx[k] <= ppx[k - 1])
			RMKMS(_("'%s' slot is not increasing"), "px");
		nr = ppi[k] - ppi[k - 1];
		nc = psuper[k] - psuper[k - 1];
		if (nr < nc)
			RMKMS(_("first differences of '%s' slot are less than those of '%s' slot"),
			      "pi", "super");
		if ((int_fast64_t) nr * nc > INT_MAX)
			RMKMS(_("supernode lengths exceed %s"), "2^31-1");
		l = nr * nc;
		if (ppx[k] - ppx[k - 1] != l)
			RMKMS(_("first differences of '%s' slot are not equal to supernode lengths"),
			      "px");
		if (l > ml)
			ml = l;
	}

	if (TYPEOF(s) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "s", "integer");
	if (XLENGTH(s) != ppi[nsuper])
		RMKMS(_("'%s' slot does not have length %s"), "s", "pi[length(pi)]");
	int i, j, *ps = INTEGER(s);
	for (k = 1; k <= nsuper; ++k) {
		nr = ppi[k] - ppi[k-1];
		nc = psuper[k] - (j = psuper[k-1]);
		for (i = 0; i < nr; ++i) {
			if (ps[i] == NA_INTEGER)
				RMKMS(_("'%s' slot contains NA"), "s");
			if (ps[i] < 0 || ps[i] >= n)
				RMKMS(_("'%s' slot has elements not in {%s}"),
				      "s", "0,...,Dim[1]-1");
			if (i < nc) {
				if (ps[i] != j + i)
					RMKMS(_("'%s' slot is wrong within diagonal blocks (row and column indices do not coincide)"), "s");
			} else {
				if (ps[i] <= ps[i-1])
					RMKMS(_("'%s' slot is not increasing within supernodes"), "s");
			}
		}
		ps += nr;
	}

	if (TYPEOF(maxcsize) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "maxesize", "integer");
	if (XLENGTH(maxcsize) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "maxcsize", 1);
	int mc = INTEGER(maxcsize)[0];
	if (mc < 0 || mc > ml)
		RMKMS(_("'%s' slot is negative or exceeds maximum supernode length"), "maxcsize");

	if (TYPEOF(maxesize) != INTSXP)
		RMKMS(_("'%s' slot is not of type \"%s\""), "maxesize", "integer");
	if (XLENGTH(maxesize) != 1)
		RMKMS(_("'%s' slot does not have length %d"), "maxesize", 1);
	int me = INTEGER(maxesize)[0];
	if (me < 0 || me > n)
		RMKMS(_("'%s' slot is negative or exceeds %s"), "maxesize", "Dim[1]");

	if (pattern)
		return Rf_ScalarLogical(1);

	if (TYPEOF(x) != REALSXP && TYPEOF(x) != CPLXSXP)
		RMKMS(_("'%s' slot is not of type \"%s\" or \"%s\""),
		      "x", "double", "complex");
	if (XLENGTH(x) != ppx[nsuper])
		RMKMS(_("'%s' slot does not have length %s"), "x", "px[length(px)]");

	/* Real, non-negative diagonal elements are necessary and sufficient */
	if (TYPEOF(x) == REALSXP) {
	double *pu = REAL(x), *pv;
	for (k = 0; k < nsuper; ++k) {
		nr = ppi[k+1] - ppi[k];
		nc = psuper[k+1] - psuper[k];
		pv = pu + ppx[k];
		for (j = 0; j < nc; ++j) {
			if (!ISNAN(*pv) && *pv < 0.0)
				RMK(_("Cholesky factor has diagonal elements with negative real part"));
			pv += (R_xlen_t) nr + 1;
		}
	}
	} else {
	Rcomplex *pu = COMPLEX(x), *pv;
	for (k = 0; k < nsuper; ++k) {
		nr = ppi[k+1] - ppi[k];
		nc = psuper[k+1] - psuper[k];
		pv = pu + ppx[k];
		for (j = 0; j < nc; ++j) {
			if (!ISNAN((*pv).r) && (*pv).r < 0.0)
				RMK(_("Cholesky factor has diagonal elements with negative real part"));
			pv += (R_xlen_t) nr + 1;
		}
	}
	}

	return Rf_ScalarLogical(1);
}

/* where 'cl' must be an element of 'valid_matrix' */
void validObject(SEXP obj, const char *cl)
{
#ifndef MATRIX_DISABLE_VALIDITY

	SEXP status;

# define IS_VALID(_CLASS_) \
	do { \
		status = R_valid_ ## _CLASS_(obj); \
		if (TYPEOF(status) == STRSXP) \
			Rf_error(_("invalid class \"%s\" object: %s"), \
			         cl, CHAR(STRING_ELT(status, 0))); \
	} while (0)

#define IS_VALID_SPARSE(_C_) \
	do { \
		IS_VALID(_C_ ## sparseMatrix); \
		if (cl[0] == 'n') { \
			if (cl[1] == 's') \
				IS_VALID(s ## _C_ ## Matrix); \
			else if (cl[1] == 't') \
				IS_VALID(t ## _C_ ## Matrix); \
		} else { \
			if (cl[1] == 'g') \
				IS_VALID(xg ## _C_ ## Matrix); \
			else if (cl[1] == 's' || cl[1] == 'p') \
				IS_VALID(xs ## _C_ ## Matrix); \
			else if (cl[1] == 't') \
				IS_VALID(xt ## _C_ ## Matrix); \
			if (cl[1] == 'p') \
				IS_VALID(xp ## _C_ ## Matrix); \
		} \
	} while (0)

	IS_VALID(Matrix);

	const char *cl_ = cl;
	if (cl[0] == 'c')
		cl = (cl[2] != 'p') ? "dpoMatrix" : "dppMatrix";
	else if (cl[0] == 'p')
		cl = "indMatrix";

	if (cl[0] == 'i' && cl[1] == 'n' && cl[2] == 'd') {
		IS_VALID(indMatrix);
		if (cl_[0] == 'p')
			IS_VALID(pMatrix);
		return;
	}

	if (cl[0] == 'n' && cl[2] != 'C' && cl[2] != 'R' && cl[2] != 'T')
		IS_VALID(nMatrix);
	else if (cl[0] == 'l')
		IS_VALID(lMatrix);
	else if (cl[0] == 'i')
		IS_VALID(iMatrix);
	else if (cl[0] == 'd')
		IS_VALID(dMatrix);
	else if (cl[0] == 'z')
		IS_VALID(zMatrix);

	if (cl[1] == 'g')
		IS_VALID(generalMatrix);
	else if (cl[1] == 's' || cl[1] == 'p')
		IS_VALID(symmetricMatrix);
	else if (cl[1] == 't')
		IS_VALID(triangularMatrix);
	else if (cl[1] == 'd') {
		IS_VALID(diagonalMatrix);
		return;
	}

	if (cl[2] == 'C')
		IS_VALID_SPARSE(C);
	else if (cl[2] == 'R')
		IS_VALID_SPARSE(R);
	else if (cl[2] == 'T')
		IS_VALID_SPARSE(T);
	else if (cl[2] != 'p') {
		IS_VALID(unpackedMatrix);
		if (cl[1] == 'p') {
			IS_VALID(xpoMatrix);
			if (cl_[0] == 'c')
				IS_VALID(corMatrix);
		}
	} else {
		IS_VALID(packedMatrix);
		if (cl[1] == 'p') {
			IS_VALID(xppMatrix);
			if (cl_[0] == 'c')
				IS_VALID(copMatrix);
		}
	}

# undef IS_VALID_SPARSE
# undef IS_VALID

#endif

	return;
}
