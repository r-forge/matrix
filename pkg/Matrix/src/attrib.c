#include "Mdefines.h"

/* .... Dim ......................................................... */

SEXP R_Dim_prod(SEXP dim)
{
	SEXP ans;
	int m = INTEGER(dim)[0], n = INTEGER(dim)[1];
	int_fast64_t mn = (int_fast64_t) m * n;
	if (mn <= INT_MAX) {
		ans = Rf_allocVector(INTSXP, 1);
		INTEGER(ans)[0] = (int) mn;
	} else {
		int_fast64_t mn_ = (int_fast64_t) (double) mn;
		if (mn_ > mn)
			mn_ = (int_fast64_t) nextafter((double) mn, 0.0);
		ans = Rf_allocVector(REALSXP, 1);
		REAL(ans)[0] = (double) mn_;
		if (mn_ != mn) {
			SEXP off;
			PROTECT(ans);
			PROTECT(off = Rf_allocVector(REALSXP, 1));
			REAL(off)[0] = (double) (mn - mn_);
			Rf_setAttrib(ans, Matrix_offSym, off);
#if 0
			Rf_warning(_("true length %llu truncated to %llu"),
			           mn, mn_);
#endif
			UNPROTECT(2);
		}
	}
	return ans;
}


/* .... Dimnames .................................................... */

int DimNames_is_trivial(SEXP dn)
{
	return
		VECTOR_ELT(dn, 0) == R_NilValue &&
		VECTOR_ELT(dn, 1) == R_NilValue &&
		Rf_getAttrib(dn, R_NamesSymbol) == R_NilValue;
}

int DimNames_is_symmetric(SEXP dn)
{
	SEXP rn, cn, ndn;
	const char *nrn, *ncn;
	int n;

	return
		!(((rn = VECTOR_ELT(dn, 0)) != R_NilValue &&
		   (cn = VECTOR_ELT(dn, 1)) != R_NilValue &&
		   rn != cn &&
		   ((n = LENGTH(rn)) != LENGTH(cn) ||
		    !equalString(rn, cn, n))) ||
		  (((ndn = Rf_getAttrib(dn, R_NamesSymbol)) != R_NilValue &&
		    *(nrn = CHAR(STRING_ELT(ndn, 0))) != '\0' &&
		    *(ncn = CHAR(STRING_ELT(ndn, 1))) != '\0' &&
		    strcmp(nrn, ncn) != 0)));
}

SEXP R_DimNames_is_symmetric(SEXP dn)
{
	return Rf_ScalarLogical(DimNames_is_symmetric(dn));
}

void symDN(SEXP dest, SEXP src, int J /* -1|0|1 */)
{
	J = (J < 0) ? -1 : (J != 0);
	SEXP s;
	if (J < 0) {
		if ((s = VECTOR_ELT(src, J = 1)) != R_NilValue ||
		    (s = VECTOR_ELT(src, J = 0)) != R_NilValue) {
			SET_VECTOR_ELT(dest, 0, s);
			SET_VECTOR_ELT(dest, 1, s);
		} else
			J = 1;
	} else {
		if ((s = VECTOR_ELT(src, J)) != R_NilValue) {
			SET_VECTOR_ELT(dest, 0, s);
			SET_VECTOR_ELT(dest, 1, s);
		}
	}
	PROTECT(s = Rf_getAttrib(src, R_NamesSymbol));
	if (s != R_NilValue) {
		SEXP destnms = PROTECT(Rf_allocVector(STRSXP, 2));
		if (CHAR(s = STRING_ELT(s, J))[0] != '\0') {
			SET_STRING_ELT(destnms, 0, s);
			SET_STRING_ELT(destnms, 1, s);
		}
		Rf_setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

SEXP R_symDN(SEXP dn)
{
	if (DimNames_is_trivial(dn))
		return dn;
	SEXP value = PROTECT(Rf_allocVector(VECSXP, 2));
	symDN(value, dn, -1);
	UNPROTECT(1);
	return value;
}

void cpyDN(SEXP dest, SEXP src, int J /* 0|1 */)
{
	J = J != 0;
	SEXP s;
	if ((s = VECTOR_ELT(src, 0)) != R_NilValue)
		SET_VECTOR_ELT(dest,  J, s);
	if ((s = VECTOR_ELT(src, 1)) != R_NilValue)
		SET_VECTOR_ELT(dest, !J, s);
	PROTECT(s = Rf_getAttrib(src, R_NamesSymbol));
	if (s != R_NilValue) {
		SEXP srcnms = s, destnms = PROTECT(Rf_allocVector(STRSXP, 2));
		if (CHAR(s = STRING_ELT(srcnms, 0))[0] != '\0')
			SET_STRING_ELT(destnms,  J, s);
		if (CHAR(s = STRING_ELT(srcnms, 1))[0] != '\0')
			SET_STRING_ELT(destnms, !J, s);
		Rf_setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return;
}

SEXP (DIMNAMES)(SEXP obj, int mode)
{
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym);
	if (mode != 0) {
		PROTECT(dn);
		if (!DimNames_is_trivial(dn)) {
			SEXP value = PROTECT(Rf_allocVector(VECSXP, 2));
			if (mode < 0)
				symDN(value, dn, mode);
			else
				cpyDN(value, dn, mode);
			UNPROTECT(1);
			dn = value;
		}
		UNPROTECT(1);
	}
	return dn;
}

void (SET_DIMNAMES)(SEXP obj, int mode, SEXP value)
{
	PROTECT(value);
	if (!DimNames_is_trivial(value)) {
		SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym));
		if (mode < 0)
			symDN(dn, value, mode);
		else
			cpyDN(dn, value, mode);
		UNPROTECT(1);
	}
	UNPROTECT(1);
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
		nms = PROTECT(Rf_getAttrib(factors, R_NamesSymbol)),
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
		nms = PROTECT(Rf_getAttrib(factors, R_NamesSymbol));
	R_xlen_t i = strmatch(nm, nms);
	if (i >= 0) {
		SET_VECTOR_ELT(factors, i, val);
		UNPROTECT(4);
		return;
	}
	R_xlen_t n = XLENGTH(factors);
	SEXP factors1 = PROTECT(Rf_allocVector(VECSXP, n + 1)),
		nms1 = PROTECT(Rf_allocVector(STRSXP, n + 1));
	for (i = 0; i < n; ++i) {
		SET_VECTOR_ELT(factors1, i, VECTOR_ELT(factors, i));
		if (nms != R_NilValue)
		SET_STRING_ELT(    nms1, i, STRING_ELT(    nms, i));
	}
	SET_VECTOR_ELT(factors1, n,           val);
	SET_STRING_ELT(    nms1, n, Rf_mkChar(nm));
	Rf_setAttrib(factors1, R_NamesSymbol, nms1);
	SET_SLOT(obj, Matrix_factorsSym, factors1);
	UNPROTECT(6);
	return;
}

SEXP R_set_factor(SEXP s_obj, SEXP s_nm, SEXP s_val, SEXP s_warn)
{
	if (TYPEOF(s_nm) != STRSXP || LENGTH(s_nm) < 1 ||
	    (s_nm = STRING_ELT(s_nm, 0)) == NA_STRING)
		Rf_error(_("invalid factor name"));
	else if (TYPEOF(Rf_getAttrib(s_obj, Matrix_factorsSym)) == VECSXP)
		set_factor(s_obj, CHAR(s_nm), s_val);
	else if (Rf_asLogical(s_warn))
		Rf_warning(_("attempt to set factor on %s without '%s' slot"),
		           "Matrix", "factors");
	return s_val;
}
