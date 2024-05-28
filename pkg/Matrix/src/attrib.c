#include "Mdefines.h"
#include "attrib.h"

/* .... Dimnames .................................................... */

int DimNames_is_trivial(SEXP dn)
{
	return
		isNull(VECTOR_ELT(dn, 0)) &&
		isNull(VECTOR_ELT(dn, 1)) &&
		isNull(getAttrib(dn, R_NamesSymbol));
}

int DimNames_is_symmetric(SEXP dn)
{
	SEXP rn, cn, ndn;
	const char *nrn, *ncn;
	int n;

	return
		!((!isNull(rn = VECTOR_ELT(dn, 0)) &&
		   !isNull(cn = VECTOR_ELT(dn, 1)) &&
		   rn != cn &&
		   ((n = LENGTH(rn)) != LENGTH(cn) ||
		    !equal_character_vectors(rn, cn, n))) ||
		  ((!isNull(ndn = getAttrib(dn, R_NamesSymbol)) &&
		    *(nrn = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
		    *(ncn = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
		    strcmp(nrn, ncn) != 0)));
}

SEXP R_DimNames_is_symmetric(SEXP s_dn)
{
	return ScalarLogical(DimNames_is_symmetric(s_dn));
}

void symDN(SEXP dest, SEXP src, int J /* -1|0|1 */)
{
	SEXP s;
	if (J < 0) {
		if (!isNull(s = VECTOR_ELT(src, J = 1)) ||
		    !isNull(s = VECTOR_ELT(src, J = 0))) {
			SET_VECTOR_ELT(dest, 0, s);
			SET_VECTOR_ELT(dest, 1, s);
		} else {
			J = 1;
		}
	} else {
		if (!isNull(s = VECTOR_ELT(src, J))) {
			SET_VECTOR_ELT(dest, 0, s);
			SET_VECTOR_ELT(dest, 1, s);
		}
	}
	PROTECT(s = getAttrib(src, R_NamesSymbol));
	if (!isNull(s)) {
		SEXP destnms = PROTECT(allocVector(STRSXP, 2));
		if (CHAR(s = STRING_ELT(s, J))[0] != '\0') {
			SET_STRING_ELT(destnms, 0, s);
			SET_STRING_ELT(destnms, 1, s);
		}
		setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

void revDN(SEXP dest, SEXP src) {
	SEXP s;
	if (!isNull(s = VECTOR_ELT(src, 0)))
		SET_VECTOR_ELT(dest, 1, s);
	if (!isNull(s = VECTOR_ELT(src, 1)))
		SET_VECTOR_ELT(dest, 0, s);
	PROTECT(s = getAttrib(src, R_NamesSymbol));
	if (!isNull(s)) {
		SEXP srcnms = s, destnms = PROTECT(allocVector(STRSXP, 2));
		if (CHAR(s = STRING_ELT(srcnms, 0))[0] != '\0')
			SET_STRING_ELT(destnms, 1, s);
		if (CHAR(s = STRING_ELT(srcnms, 1))[0] != '\0')
			SET_STRING_ELT(destnms, 0, s);
		setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

SEXP R_symDN(SEXP s_dn)
{
	if (DimNames_is_trivial(s_dn))
		return s_dn;
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symDN(newdn, s_dn, -1);
	UNPROTECT(1);
	return newdn;
}

SEXP R_revDN(SEXP s_dn)
{
	if (DimNames_is_trivial(s_dn))
		return s_dn;
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	revDN(newdn, s_dn);
	UNPROTECT(1);
	return newdn;
}

SEXP get_symmetrized_DimNames(SEXP obj, int J) {
	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	if (DimNames_is_trivial(dn)) {
		UNPROTECT(1);
		return dn;
	}
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	symDN(newdn, dn, J);
	UNPROTECT(2);
	return newdn;
}

SEXP get_reversed_DimNames(SEXP obj) {
	SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
	if (DimNames_is_trivial(dn)) {
		UNPROTECT(1);
		return dn;
	}
	SEXP newdn = PROTECT(allocVector(VECSXP, 2));
	revDN(newdn, dn);
	UNPROTECT(2);
	return newdn;
}

void set_symmetrized_DimNames(SEXP obj, SEXP dn, int J) {
	if (!DimNames_is_trivial(dn)) {
		SEXP newdn = PROTECT(allocVector(VECSXP, 2));
		symDN(newdn, dn, J);
		SET_SLOT(obj, Matrix_DimNamesSym, newdn);
		UNPROTECT(1);
	}
	return;
}

void set_reversed_DimNames(SEXP obj, SEXP dn) {
	if (!DimNames_is_trivial(dn)) {
		SEXP newdn = PROTECT(allocVector(VECSXP, 2));
		revDN(newdn, dn);
		SET_SLOT(obj, Matrix_DimNamesSym, newdn);
		UNPROTECT(1);
	}
	return;
}


/* .... factors ..................................................... */

static
R_xlen_t strmatch(const char *s, SEXP nms)
{
	if (TYPEOF(nms) == STRSXP) {
		for (R_xlen_t i = 0, n = XLENGTH(nms); i < n; ++i)
			if (strcmp(s, CHAR(STRING_ELT(nms, i))) == 0)
				return i;
	}
	return (R_xlen_t) -1;
}

SEXP get_factor(SEXP obj, const char *nm)
{
	SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorsSym)),
		nms = PROTECT(getAttrib(factors, R_NamesSymbol)),
		val = R_NilValue;
	R_xlen_t i = strmatch(nm, nms);
	if (i >= 0)
		val = VECTOR_ELT(factors, i);
	UNPROTECT(2);
	return val;
}

void set_factor(SEXP obj, const char *nm, SEXP val)
{
	PROTECT(obj);
	PROTECT(val);
	SEXP factors = PROTECT(GET_SLOT(obj, Matrix_factorsSym)),
		nms = PROTECT(getAttrib(factors, R_NamesSymbol));
	R_xlen_t i = strmatch(nm, nms);
	if (i >= 0) {
		SET_VECTOR_ELT(factors, i, val);
		UNPROTECT(4);
		return;
	}
	R_xlen_t n = XLENGTH(factors);
	SEXP factors1 = PROTECT(allocVector(VECSXP, n + 1)),
		nms1 = PROTECT(allocVector(STRSXP, n + 1));
	for (i = 0; i < n; ++i) {
		SET_VECTOR_ELT(factors1, i, VECTOR_ELT(factors, i));
		if (nms != R_NilValue)
		SET_STRING_ELT(    nms1, i, STRING_ELT(    nms, i));
	}
	SET_VECTOR_ELT(factors1, n,                 val);
	SET_STRING_ELT(    nms1, n, PROTECT(mkChar(nm)));
	setAttrib(factors1, R_NamesSymbol, nms1);
	SET_SLOT(obj, Matrix_factorsSym, factors1);
	UNPROTECT(7);
	return;
}

SEXP R_set_factor(SEXP s_obj, SEXP s_nm, SEXP s_val, SEXP s_warn)
{
	if (TYPEOF(s_nm) != STRSXP || LENGTH(s_nm) < 1 ||
	    (s_nm = STRING_ELT(s_nm, 0)) == NA_STRING)
		error(_("invalid factor name"));
	else if (TYPEOF(getAttrib(s_obj, Matrix_factorsSym)) == VECSXP)
		set_factor(s_obj, CHAR(s_nm), s_val);
	else if (asLogical(s_warn) != 0)
		warning(_("attempt to set factor on %s without '%s' slot"),
		        "Matrix", "factors");
	return s_val;
}


/* .... x ........................................................... */

SEXP get_forced_x(SEXP obj, const char *class)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	if (class[0] != 'z' || (class[1] != 's' && class[1] != 'p'))
		return x;
	
	SEXP trans = GET_SLOT(obj, Matrix_transSym);
	char ct = CHAR(STRING_ELT(trans, 0))[0];
	if (ct != 'C')
		return x;

	PROTECT(x);
	Rcomplex *px = COMPLEX(x);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0], status = 0;
	
	switch (class[2]) {
	case 'y':
	case 'o':
	{
		R_xlen_t k = 0, n1a = (R_xlen_t) n + 1;
		for (int j = 0; j < n; ++j) {
			if (px[k].i != 0.0) {
				status = 1;
				break;
			}
			k += n1a;
		}
		if (status) {
			SEXP y = allocVector(CPLXSXP, XLENGTH(x));
			Rcomplex *py = COMPLEX(y);
			memcpy(py, px, (size_t) XLENGTH(x) * sizeof(Rcomplex));
			k = 0;
			for (int j = 0; j < n; ++j) {
				py[k].i = 0.0;
				k += n1a;
			}
			x = y;
		}
		break;
	}
	case 'p':
	{
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		int upper = ul == 'U';

		R_xlen_t k = 0;
		for (int j = 0; j < n; ++j) {
			if (px[k].i != 0.0) {
				status = 1;
				break;
			}
			k += (upper) ? j + 2 : n - j;
		}
		if (status) {
			SEXP y = allocVector(CPLXSXP, XLENGTH(x));
			Rcomplex *py = COMPLEX(y);
			memcpy(py, px, (size_t) XLENGTH(x) * sizeof(Rcomplex));
			k = 0;
			for (int j = 0; j < n; ++j) {
				py[k].i = 0.0;
				k += (upper) ? j + 2 : n - j;
			}
			x = y;
		}
		break;
	}
	case 'C':
	case 'R':
	{
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		int upper = (ul == 'U') == (class[2] == 'C');

		SEXP p = PROTECT(GET_SLOT(obj, Matrix_pSym));
		int *pp = INTEGER(p) + 1, k, kend;

		for (int j = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && px[(upper) ? kend - 1 : k].i != 0.0) {
				status = 1;
				break;
			}
			k = kend;
		}
		if (status) {
			SEXP y = allocVector(CPLXSXP, XLENGTH(x));
			Rcomplex *py = COMPLEX(y);
			memcpy(py, px, (size_t) XLENGTH(x) * sizeof(Rcomplex));
			for (int j = 0; j < n; ++j) {
				kend = pp[j];
				if (k < kend)
					py[(upper) ? kend - 1 : k].i = 0.0;
				k = kend;
			}
			x = y;
		}
		UNPROTECT(1);
		break;
	}
	case 'T':
	{
		SEXP i = PROTECT(GET_SLOT(obj, Matrix_iSym)),
			j = PROTECT(GET_SLOT(obj, Matrix_jSym));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		R_xlen_t nnz = XLENGTH(x);
		
		for (R_xlen_t k = 0; k < nnz; ++k) {
			if (pi[k] == pj[k] && px[k].i != 0.0) {
				status = 1;
				break;
			}
		}
		if (status) {
			SEXP y = allocVector(CPLXSXP, nnz);
			Rcomplex *py = COMPLEX(y);
			for (R_xlen_t k = 0; k < nnz; ++k) {
				py[k].r =                          px[k].r;
				py[k].i = (pi[k] == pj[k]) ? 0.0 : px[k].i;
			}
			x = y;
		}
		UNPROTECT(2);
	}
	default:
		break;
	}
	UNPROTECT(1);
	return x;
}
