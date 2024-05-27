#include "Mdefines.h"
#include "idz.h"
#include "dense.h"

SEXP dense_band(SEXP from, const char *class, int a, int b)
{
	int packed = class[2] == 'p';

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
	if (a <= 1 - m && b >= n - 1 &&
	    (class[1] == 't' || m != n || m > 1 || n > 1))
		return from;

	int ge = 0, sy = 0, tr = 0;
	ge = m != n || !((tr = a >= 0 || b <= 0 || class[1] == 't') ||
	                 (sy = a == -b && class[1] == 's'));

#define BAND_CASES(_DO_) \
	do { \
		switch (class[0]) { \
		case 'n': \
		case 'l': \
			_DO_(i, int, LOGICAL); \
			break; \
		case 'i': \
			_DO_(i, int, INTEGER); \
			break; \
		case 'd': \
			_DO_(d, double, REAL); \
			break; \
		case 'z': \
			_DO_(z, Rcomplex, COMPLEX); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define BAND(_PREFIX_, _CTYPE_, _PTR_) \
	_PREFIX_ ## band2(_PTR_(x1), m, n, a, b)

	if (class[1] != 'g' && ge) {
		/* defined in ./coerce.c : */
		SEXP dense_as_general(SEXP, const char *, int);
		PROTECT(from = dense_as_general(from, class, 1));
		SEXP x1 = PROTECT(GET_SLOT(from, Matrix_xSym));
		BAND_CASES(BAND);
		UNPROTECT(2); /* x1, from */
		return from;
	}

	char ul0 = 'U', ul1 = 'U', ct = 'C', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul0 = CHAR(STRING_ELT(uplo, 0))[0];

		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = GET_SLOT(from, Matrix_transSym);
			ct = CHAR(STRING_ELT(trans, 0))[0];
		}

		if (class[1] == 't') {
			/* Be fast if band contains entire triangle */
			if ((ul0 == 'U') ? (a <= 0 && b >= n - 1) : (b >= 0 && a <= 1 - m))
				return from;
			SEXP diag = GET_SLOT(from, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (ge) ? 'g' :                   ((sy) ? 's' : 't') ;
	cl[2] = (ge) ? 'e' : ((packed) ? 'p' : ((sy) ? 'y' : 'r'));
	SEXP to = PROTECT(newObject(cl));

	dim = GET_SLOT(to, Matrix_DimSym);
	pdim = INTEGER(dim);
	pdim[0] = m;
	pdim[1] = n;

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's' && !sy)
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1;

	if (class[1] == 'g' && ge) {
		PROTECT(x1 = duplicate(x0));
		if (ATTRIB(x1) != R_NilValue) {
			SET_ATTRIB(x1, R_NilValue);
			if (OBJECT(x1))
				SET_OBJECT(x1, 0);
		}
		SET_SLOT(to, Matrix_xSym, x1);
		BAND_CASES(BAND);
		UNPROTECT(3); /* x1, x0, to */
		return to;
	}

#undef BAND

#define BAND(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px1 = _PTR_(x1); \
		if (!packed) \
			_PREFIX_ ## band2(px1, m, n, a, b); \
		else \
			_PREFIX_ ## band1(px1,    n, a, b, ul1); \
	} while (0)

#define DCOPY(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		Matrix_memset(px1, 0, XLENGTH(x1), sizeof(_CTYPE_)); \
		if (di == 'N' && a <= 0 && b >= 0) { \
		if (!packed) \
			_PREFIX_ ## dcopy2(px1, px0, n, XLENGTH(x1),      ul0, di); \
		else \
			_PREFIX_ ## dcopy1(px1, px0, n, XLENGTH(x1), ul1, ul0, di); \
		} \
	} while (0)

#define TRANS(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (!packed) \
			_PREFIX_ ## trans2(px1, px0, m, n,      ct); \
		else \
			_PREFIX_ ## trans1(px1, px0,    n, ul0, ct); \
	} while (0)

	/* Returning .(sy|sp|tr|tp)Matrix ... */

	ul1 = (class[1] == 't' || sy) ? ul0 : ((a >= 0) ? 'U' : 'L');
	if (ul1 != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (ct != 'C' && sy) {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (di != 'N' && tr && a <= 0 && b >= 0) {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	if (class[1] == 't') {
		if ((ul0 == 'U') ? (b <= 0) : (a >= 0)) {
			/* Result is either a diagonal matrix or a zero matrix : */
			PROTECT(x1 = allocVector(TYPEOF(x0), XLENGTH(x0)));
			BAND_CASES(DCOPY);
		} else {
			PROTECT(x1 = duplicate(x0));
			BAND_CASES(BAND);
		}
	} else {
		if (sy || (tr && (class[1] == 'g' || ul0 == ul1 || n <= 1))) {
			PROTECT(x1 = duplicate(x0));
			if (ATTRIB(x1) != R_NilValue) {
				SET_ATTRIB(x1, R_NilValue);
				if (OBJECT(x1))
					SET_OBJECT(x1, 0);
			}
		} else {
			/* Band is "opposite" the stored triangle : */
			PROTECT(x1 = allocVector(TYPEOF(x0), XLENGTH(x0)));
			BAND_CASES(TRANS);
		}
		BAND_CASES(BAND);
		if (class[1] == 's' && class[0] == 'z' && ct == 'C' &&
		    !sy && (a == 0 || b == 0)) {
			if (!packed)
				zdreal2(COMPLEX(x1), n);
			else
				zdreal1(COMPLEX(x1), n, ul1);
		}
	}
	SET_SLOT(to, Matrix_xSym, x1);

#undef BAND_CASES
#undef BAND
#undef DCOPY
#undef TRANS

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* band(<denseMatrix>, k1, k2), tri[ul](<denseMatrix>, k) */
/* NB: argument validation more or less copied by R_sparse_band() */
SEXP R_dense_band(SEXP s_from, SEXP s_a, SEXP s_b)
{
	if (TYPEOF(s_from) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_from = matrix_as_dense(s_from, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_from);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_from, __func__);

	SEXP dim = GET_SLOT(s_from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	int a, b;
	if (s_a == R_NilValue)
		a = -m;
	else if ((a = asInteger(s_a)) == NA_INTEGER || a < -m || a > n)
		error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		      "k1", a, "-Dim[1]", -m, "Dim[2]", n);
	if (s_b == R_NilValue)
		b = n;
	else if ((b = asInteger(s_b)) == NA_INTEGER || b < -m || b > n)
		error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		      "k2", b, "-Dim[1]", -m, "Dim[2]", n);
	else if (b < a)
		error(_("'%s' (%d) must be less than or equal to '%s' (%d)"),
		      "k1", a, "k2", b);

	s_from = dense_band(s_from, valid[ivalid], a, b);
	UNPROTECT(1);
	return s_from;
}

SEXP dense_diag_get(SEXP obj, const char *class, int names)
{
	int packed = class[2] == 'p';

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n, j;

	char ul = 'U', ct = 'C', di = 'N';
	if (class[1] != 'g') {
		if (packed) {
			SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
			ul = CHAR(STRING_ELT(uplo, 0))[0];
		}
		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = GET_SLOT(obj, Matrix_transSym);
			ct = CHAR(STRING_ELT(trans, 0))[0];
		}
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		ans = PROTECT(allocVector(TYPEOF(x), r));

#define DG_LOOP(_CTYPE_, _PTR_, _ONE_) \
	do { \
		_CTYPE_ *pans = _PTR_(ans), *px = _PTR_(x); \
		if (di != 'N') \
			for (j = 0; j < r; ++j) \
				*(pans++) = _ONE_; \
		else if (!packed) { \
			R_xlen_t m1a = (R_xlen_t) m + 1; \
			for (j = 0; j < r; ++j, px += m1a) \
				*(pans++) = *px; \
		} \
		else if (ul == 'U') \
			for (j = 0; j < n; px += (++j) + 1) \
				*(pans++) = *px; \
		else \
			for (j = 0; j < n; px += n - (j++)) \
				*(pans++) = *px; \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		DG_LOOP(int, LOGICAL, 1);
		break;
	case 'i':
		DG_LOOP(int, INTEGER, 1);
		break;
	case 'd':
		DG_LOOP(double, REAL, 1.0);
		break;
	case 'z':
		DG_LOOP(Rcomplex, COMPLEX, Matrix_zone);
		if (class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), XLENGTH(ans));
		break;
	default:
		break;
	}

	if (names) {
		/* NB: The logic here must be adjusted once the validity method
		   for 'symmetricMatrix' enforce symmetric 'Dimnames'
		*/
		SEXP dn = PROTECT(GET_SLOT(obj, Matrix_DimNamesSym)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (cn == R_NilValue) {
			if (class[1] == 's' && rn != R_NilValue)
				setAttrib(ans, R_NamesSymbol, rn);
		} else {
			if (class[1] == 's')
				setAttrib(ans, R_NamesSymbol, cn);
			else if (rn != R_NilValue &&
			         (rn == cn || equal_character_vectors(rn, cn, r)))
				setAttrib(ans, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

#undef DG_LOOP

	UNPROTECT(2); /* x, ans */
	return ans;
}

/* diag(<denseMatrix>, names=) */
SEXP R_dense_diag_get(SEXP s_obj, SEXP s_names)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	int names;
	if (TYPEOF(s_names) != LGLSXP || LENGTH(s_names) < 1 ||
	    (names = LOGICAL(s_names)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "names", "TRUE", "FALSE");

	return dense_diag_get(s_obj, valid[ivalid], names);
}

SEXP dense_diag_set(SEXP from, const char *class, SEXP value, int new)
{
	int packed = class[2] == 'p';

	SEXP to = PROTECT(newObject(class));
	int v = LENGTH(value) != 1;

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n, j;
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (new) {
		x = duplicate(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
	}
	SET_SLOT(to, Matrix_xSym, x);

#define DS_LOOP(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pvalue = _PTR_(value); \
		if (!packed) { \
			R_xlen_t m1a = (R_xlen_t) m + 1; \
			if (v) \
				for (j = 0; j < r; ++j, px += m1a) \
					*px = *(pvalue++); \
			else \
				for (j = 0; j < r; ++j, px += m1a) \
					*px = *pvalue; \
		} else if (ul == 'U') { \
			if (v) \
				for (j = 0; j < n; px += (++j) + 1) \
					*px = *(pvalue++); \
			else \
				for (j = 0; j < n; px += (++j) + 1) \
					*px = *pvalue; \
		} else { \
			if (v) \
				for (j = 0; j < n; px += n - (j++)) \
					*px = *(pvalue++); \
			else \
				for (j = 0; j < n; px += n - (j++)) \
					*px = *pvalue; \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		DS_LOOP(int, LOGICAL);
		break;
	case 'i':
		DS_LOOP(int, INTEGER);
		break;
	case 'd':
		DS_LOOP(double, REAL);
		break;
	case 'z':
		DS_LOOP(Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

#undef DS_LOOP

	UNPROTECT(2); /* x, to */
	return to;
}

/* diag(<denseMatrix>) <- value */
SEXP R_dense_diag_set(SEXP s_from, SEXP s_value)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_from, __func__);
	const char *class = valid[ivalid];

	SEXPTYPE tx = kindToType(class[0]), tv = TYPEOF(s_value);

	switch (tv) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		error(_("replacement diagonal has incompatible type \"%s\""),
		      type2char(tv));
		break;
	}

	SEXP dim = GET_SLOT(s_from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	R_xlen_t len = XLENGTH(s_value);
	if (len != 1 && len != r)
		error(_("replacement diagonal has wrong length"));

	int new = 1;
	if (tv <= tx) {
		/* defined in ./coerce.c : */
		SEXP dense_as_general(SEXP, const char *, int);
		if (class[1] == 's' && class[0] == 'z' && tv == tx) {
			SEXP trans = GET_SLOT(s_from, Matrix_transSym);
			int ct = CHAR(STRING_ELT(trans, 0))[0];
			if (ct == 'C') {
				PROTECT(s_from = dense_as_general(s_from, class, 1));
				class = valid[R_check_class_etc(s_from, valid)];
				new = 0;
				UNPROTECT(1);
			}
		}
		PROTECT(s_from);
		PROTECT(s_value = coerceVector(s_value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
#ifndef MATRIX_ENABLE_IMATRIX
		if (tv == INTSXP) {
		PROTECT(s_from = dense_as_kind(s_from, class, 'd', 0));
		PROTECT(s_value = coerceVector(s_value, REALSXP));
		} else {
#endif
		PROTECT(s_from = dense_as_kind(s_from, class, typeToKind(tv), 0));
		PROTECT(s_value);
#ifndef MATRIX_ENABLE_IMATRIX
		}
#endif
		class = valid[R_check_class_etc(s_from, valid)];
		new = 0;
	}

	s_from = dense_diag_set(s_from, class, s_value, new);
	UNPROTECT(2);
	return s_from;
}

SEXP dense_transpose(SEXP from, const char *class, char ct)
{
	int packed = class[2] == 'p';

	SEXP to = PROTECT(newObject(class));

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n) {
		dim = GET_SLOT(to, Matrix_DimSym);
		pdim = INTEGER(dim);
		pdim[0] = n;
		pdim[1] = m;
	} else if (n > 0) {
		PROTECT(dim);
		SET_SLOT(to, Matrix_DimSym, dim);
		UNPROTECT(1); /* dim */
	}

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's' || class[1] == 'p' || class[1] == 'o')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_reversed_DimNames(to, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul == 'U') {
			PROTECT(uplo = mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
			char ct = CHAR(STRING_ELT(trans, 0))[0];
			if (ct != 'C')
				SET_SLOT(to, Matrix_transSym, trans);
			UNPROTECT(1); /* trans */
		}
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = CHAR(STRING_ELT(diag, 0))[0];
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorsSym, factors);
			UNPROTECT(1); /* factors */
			if (class[1] == 'o' && n > 0) {
				SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
				SET_SLOT(to, Matrix_sdSym, sd);
				UNPROTECT(1); /* sd */
			}
		}
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
	SET_SLOT(to, Matrix_xSym, x1);

#define TRANS(_PREFIX_, _CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (!packed) \
			_PREFIX_ ## trans2(px1, px0, m, n,     ct); \
		else \
			_PREFIX_ ## trans1(px1, px0,    n, ul, ct); \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		TRANS(i, int, LOGICAL);
		break;
	case 'i':
		TRANS(i, int, INTEGER);
		break;
	case 'c':
	case 'd':
		TRANS(d, double, REAL);
		break;
	case 'z':
		TRANS(z, Rcomplex, COMPLEX);
		break;
	default:
		break;
	}

#undef TRANS

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* t(<denseMatrix>) */
SEXP R_dense_transpose(SEXP s_from, SEXP s_trans)
{
	static const char *valid[] = {
		"corMatrix", "copMatrix",
		"dpoMatrix", "dppMatrix",
		"zpoMatrix", "zppMatrix",
		VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_from, __func__);

	char ct = 'C';
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	return dense_transpose(s_from, valid[ivalid], ct);
}

SEXP dense_force_symmetric(SEXP from, const char *class, char ul, char ct)
{
	char ul0 = 'U', ul1 = 'U', ct0 = 'C', ct1 = 'C', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul0 = ul1 = CHAR(STRING_ELT(uplo, 0))[0];
		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = GET_SLOT(from, Matrix_transSym);
			ct0 = ct1 = CHAR(STRING_ELT(trans, 0))[0];
		}
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(from, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}
	if (ul != '\0')
		ul1 = ul;
	if (ct != '\0' && class[0] == 'z')
		ct1 = ct;

	if (class[1] == 's') {
		if (ul0 != ul1)
			from = dense_transpose(from, class, ct0);
		if (ct0 == ct1 || class[0] != 'z')
			return from;
	}
	PROTECT(from);

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = 's';
	cl[2] = (packed) ? 'p' : 'y';
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to symmetrize a non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	if (ul1 != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	if (ct1 != 'C' && class[0] == 'z') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	PROTECT_INDEX pid;
	SEXP x0 = GET_SLOT(from, Matrix_xSym);
	PROTECT_WITH_INDEX(x0, &pid);
	if (class[1] == 's' && class[0] == 'z') {
		if (ul0 == ul1)
			REPROTECT(x0 = duplicate(x0), pid);
		if (!packed)
			zdreal2(COMPLEX(x0), n);
		else
			zdreal1(COMPLEX(x0), n, ul1);
	}

	if (class[1] == 'g' || class[1] == 's' || (ul0 == ul1 && di == 'N'))
		SET_SLOT(to, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
		SET_SLOT(to, Matrix_xSym, x1);

#define FORCE(_PREFIX_, _CTYPE_, _PTR_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			if (ul0 == ul1) \
				Matrix_memcpy(px1, px0, XLENGTH(x1), sizeof(_CTYPE_)); \
			else \
				Matrix_memset(px1,   0, XLENGTH(x1), sizeof(_CTYPE_)); \
			if (!packed) \
				_PREFIX_ ## dcopy2(px1, px0, n, XLENGTH(x0),      ul0, di); \
			else \
				_PREFIX_ ## dcopy1(px1, px0, n, XLENGTH(x0), ul1, ul0, di); \
		} while (0)

		switch (class[0]) {
		case 'n':
		case 'l':
			FORCE(i, int, LOGICAL);
			break;
		case 'i':
			FORCE(i, int, INTEGER);
			break;
		case 'd':
			FORCE(d, double, REAL);
			break;
		case 'z':
			FORCE(z, Rcomplex, COMPLEX);
			break;
		default:
			break;
		}

#undef FORCE

		UNPROTECT(1); /* x1 */
	}

	UNPROTECT(3); /* x0, to, from */
	return to;
}

/* forceSymmetric(<denseMatrix>, uplo, trans) */
SEXP R_dense_force_symmetric(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_from, __func__);

	char ul = '\0', ct = '\0';
	if (s_uplo != R_NilValue) {
		if (TYPEOF(s_uplo) != STRSXP || LENGTH(s_uplo) < 1 ||
		    (s_uplo = STRING_ELT(s_uplo, 0)) == NA_STRING ||
		    ((ul = CHAR(s_uplo)[0]) != 'U' && ul != 'L'))
			error(_("invalid '%s' to '%s'"), "uplo", __func__);
	}
	if (s_trans != R_NilValue) {
		if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
		    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
		    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
			error(_("invalid '%s' to '%s'"), "trans", __func__);
	}

	return dense_force_symmetric(s_from, valid[ivalid], ul, ct);
}

SEXP dense_symmpart(SEXP from, const char *class, char ct)
{
	char ct0 = ct, ct1 = ct;
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct0 = CHAR(STRING_ELT(trans, 0))[0];
	}

	if (class[0] != 'z' && class[0] != 'd') {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
		from = dense_as_kind(from, class, 'd', 0);
	}
	if (class[1] == 's' && ct0 == ct1)
		return from;
	PROTECT(from);

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = 's';
	cl[2] = (packed) ? 'p' : 'y';
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error((ct == 'C')
		      ? _("attempt to get Hermitian part of non-square matrix")
		      : _("attempt to get symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	if (class[0] == 'z' && ct1 != 'C') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	char di = 'N';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	PROTECT_INDEX pid;
	SEXP x = GET_SLOT(from, Matrix_xSym);
	PROTECT_WITH_INDEX(x, &pid);
	if (class[0] == 'z' || class[0] == 'd')
		REPROTECT(x = duplicate(x), pid);
	SET_SLOT(to, Matrix_xSym, x);

	if (class[1] == 's' && ct0 != ct1) {
		/* Symmetric part of Hermitian matrix is real part */
		/* Hermitian part of symmetric matrix is real part */
		zvreal(COMPLEX(x), XLENGTH(x));
		UNPROTECT(3); /* x, to, from */
		return to;
	}

	int i, j;

#define SP_LOOP(_PREFIX_, _CTYPE_, _PTR_, \
	            _INCREMENT_ID_, _INCREMENT_CJ_, _SCALE1_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		if (class[1] == 'g') { \
			_CTYPE_ *py = px; \
			for (j = 0; j < n; ++j) { \
				for (i = j + 1; i < n; ++i) { \
					px += n; \
					py += 1; \
					if (ct1 == 'C') \
						_INCREMENT_CJ_((*px), (*py)); \
					else \
						_INCREMENT_ID_((*px), (*py)); \
					_SCALE1_((*px), 0.5); \
				} \
				px = (py += j + 2); \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					_SCALE1_((*px), 0.5); \
					px += 1; \
				} \
				px += (!packed) ? n - j : 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				px += (!packed) ? j + 1 : 1; \
				for (i = j + 1; i < n; ++i) { \
					_SCALE1_((*px), 0.5); \
					px += 1; \
				} \
			} \
		} \
		if (di != 'N') { \
		if (!packed) \
			_PREFIX_ ## dcopy2(_PTR_(x), NULL, n, -1,     'U', di); \
		else \
			_PREFIX_ ## dcopy1(_PTR_(x), NULL, n, -1, ul, 'U', di); \
		} \
	} while (0)

	if (class[0] == 'z')
		SP_LOOP(z, Rcomplex, COMPLEX,
		        INCREMENT_COMPLEX_ID, INCREMENT_COMPLEX_CJ, SCALE1_COMPLEX);
	else
		SP_LOOP(d, double, REAL,
		        INCREMENT_REAL, INCREMENT_REAL, SCALE1_REAL);

#undef SP_LOOP

	UNPROTECT(3); /* x, to, from */
	return to;
}

/* symmpart(<denseMatrix>, trans) */
SEXP R_dense_symmpart(SEXP s_from, SEXP s_trans)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_from, __func__);

	char ct = 'C';
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	return dense_symmpart(s_from, valid[ivalid], ct);
}

SEXP dense_skewpart(SEXP from, const char *class, char ct)
{
	if (class[0] != 'z' && class[0] != 'd') {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
		from = dense_as_kind(from, class, 'd', 0);
	}
	PROTECT(from);

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 's') ? class[1] : 'g';
	cl[2] = (class[1] == 's') ? class[2] : 'e';
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error((ct == 'C')
		      ? _("attempt to get skew-Hermitian part of non-square matrix")
		      : _("attempt to get skew-symmetric part of non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U' && class[1] == 's')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	char ct0 = ct, ct1 = ct;
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		ct0 = CHAR(STRING_ELT(trans, 0))[0];
		if (ct0 != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1);
	}

	char di = 'N';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1 = x0;

	if (class[1] == 's') {
		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		/* Skew-Hermitian part of symmetric matrix is imaginary part */
		R_xlen_t len = XLENGTH(x1);
		PROTECT(x1 = allocVector(TYPEOF(x0), len));
		SET_SLOT(to, Matrix_xSym, x1);
		if (class[0] == 'z') {
			Rcomplex *px1 = COMPLEX(x1);
			Matrix_memset(px1, 0, len, sizeof(Rcomplex));
			if (ct0 != ct1) {
				Rcomplex *px0 = COMPLEX(x1);
				while (len--)
					(*(px1++)).i = (*(px0++)).i;
			}
		} else {
			double *px1 = REAL(x1);
			Matrix_memset(px1, 0, len, sizeof(double));
		}
		UNPROTECT(4); /* x1, x0, to, from */
		return to;
	}

	if (class[0] == 'z' || class[0] == 'd' || packed) {
		if ((int_fast64_t) n * n > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");
		x1 = allocVector(TYPEOF(x0), (R_xlen_t) n * n);
	}
	PROTECT(x1);
	SET_SLOT(to, Matrix_xSym, x1);

	int i, j;
	R_xlen_t upos = 0, lpos = 0;

#define SP_LOOP(_CTYPE_, _PTR_, _ZERO_, \
	            _ASSIGN2_ID_, _ASSIGN2_IM_, _INCREMENT_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				lpos = j; \
				for (i = 0; i < j; ++i) { \
					_ASSIGN2_ID_(px1[upos], 0.5 * px0[upos]); \
					_INCREMENT_(px1[upos], -0.5 * px0[lpos]); \
					_ASSIGN2_ID_(px1[lpos], -px1[upos]); \
					upos += 1; \
					lpos += n; \
				} \
				if (ct1 == 'C') \
					_ASSIGN2_IM_(px1[upos], px0[upos]); \
				else \
					px1[upos] = _ZERO_; \
				upos += n - j; \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				lpos = j; \
				for (i = 0; i < j; ++i) { \
					_ASSIGN2_ID_(px1[upos], 0.5 * (*px0)); \
					_ASSIGN2_ID_(px1[lpos], -px1[upos]); \
					px0 += 1; \
					upos += 1; \
					lpos += n; \
				} \
				if (ct1 == 'C' && di == 'N') \
					_ASSIGN2_IM_(px1[upos], (*px0)); \
				else \
					px1[upos] = _ZERO_; \
				px0 += (!packed) ? n - j : 1; \
				upos += n - j; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				upos = lpos; \
				if (ct1 == 'C' && di == 'N') \
					_ASSIGN2_IM_(px1[lpos], (*px0)); \
				else \
					px1[lpos] = _ZERO_; \
				for (i = j + 1; i < n; ++i) { \
					px0 += 1; \
					upos += n; \
					lpos += 1; \
					_ASSIGN2_ID_(px1[lpos], 0.5 * (*px0)); \
					_ASSIGN2_ID_(px1[upos], -px1[lpos]); \
				} \
				px0 += (!packed) ? j + 2 : 1; \
				lpos += j + 2; \
			} \
		} \
	} while (0)

	if (class[0] == 'z')
		SP_LOOP(Rcomplex, COMPLEX, Matrix_zzero,
		        ASSIGN2_COMPLEX_ID, ASSIGN2_COMPLEX_CJ, INCREMENT_COMPLEX_ID);
	else
		SP_LOOP(double, REAL, 0.0,
		        ASSIGN2_REAL_ID, ASSIGN2_REAL_IM, INCREMENT_REAL);

#undef SP_LOOP

	UNPROTECT(4); /* x1, x0, to, from */
	return to;
}

/* skewpart(<denseMatrix>, trans) */
SEXP R_dense_skewpart(SEXP s_from, SEXP s_trans)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_from, __func__);

	char ct = 'C';
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	return dense_skewpart(s_from, valid[ivalid], ct);
}

int dense_is_symmetric(SEXP obj, const char *class,
                       int exact, char ct, int checkDN)
{
	exact = exact || class[0] == 'n' || class[0] == 'l' || class[0] == 'i';
	
	if (class[1] == 's') {
		if (class[0] != 'z')
			return 1;
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		if (CHAR(STRING_ELT(trans, 0))[0] == ct)
			return 1;
		checkDN = 0;
	}

	if (checkDN) {
		SEXP dimnames = GET_SLOT(obj, Matrix_DimNamesSym);
		if (!DimNames_is_symmetric(dimnames))
			return 0;
	}

	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (class[0] == 'n' || class[0] == 'l' || class[0] == 'i' ||
		    (exact && (class[0] != 'z' || ct != 'C' || di != 'N')))
			return dense_is_diagonal(obj, class);
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n == 0 || (n == 1 && (class[0] != 'z' || ct != 'C')))
		return 1;
	if (!exact)
		return NA_LOGICAL; /* => do inexact numerical test in R */

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p';

#define IS_LOOP(_CTYPE_, _PTR_, _NOTEQUAL_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *py = px; \
		for (j = 0; j < n; px = (py += (++j) + 1)) { \
			for (i = j + 1; i < n; ++i) { \
				px += n; \
				py += 1; \
				if (_NOTEQUAL_((*px), (*py))) \
					return 0; \
			} \
		} \
		return 1; \
	} while (0)

	switch (class[0]) {
	case 'n':
		IS_LOOP(int, LOGICAL, NOTEQUAL_PATTERN);
		break;
	case 'l':
		IS_LOOP(int, LOGICAL, NOTEQUAL_LOGICAL);
		break;
	case 'i':
		IS_LOOP(int, INTEGER, NOTEQUAL_INTEGER);
		break;
	case 'd':
		IS_LOOP(double, REAL, NOTEQUAL_REAL);
		break;
	case 'z':
		if (class[1] == 'g' && ct != 'C')
		IS_LOOP(Rcomplex, COMPLEX, NOTEQUAL_COMPLEX);
		else {
		Rcomplex *px = COMPLEX(x), *py = px;
		if (class[1] == 'g') {
			for (j = 0; j < n; px = (py += (++j) + 1)) {
				if (NOTREAL_COMPLEX((*px)))
					return 0;
				for (i = j + 1; i < n; ++i) {
					px += 1;
					py += n;
					if (NOTCONJ_COMPLEX((*px), (*py)))
						return 0;
				}
			}
		} else {
			if (ul == 'U') {
			for (j = 0; j < n; ++j) {
				for (i = 0; i < j; ++i) {
					if (NOTREAL_COMPLEX((*px)))
						return 0;
					px += 1;
				}
				if (ct == 'C' && NOTREAL_COMPLEX((*px)))
					return 0;
				px += 1;
				if (!packed)
					px += n - j - 1;
			}
			} else {
				for (j = 0; j < n; ++j) {
					if (!packed)
						px += j;
					if (ct == 'C' && NOTREAL_COMPLEX((*px)))
						return 0;
					px += 1;
					for (i = j + 1; i < n; ++i) {
						if (NOTREAL_COMPLEX((*px)))
							return 0;
						px += 1;
					}
				}
			}
		}
		}
		break;
	default:
		break;
	}

#undef IS_LOOP

	return 0;
}

/* isSymmetric(<denseMatrix>, tol, tol1, trans, checkDN) */
SEXP R_dense_is_symmetric(SEXP s_obj,
                          SEXP s_exact, SEXP s_trans, SEXP s_checkDN)
{
	if (TYPEOF(s_obj) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_obj = matrix_as_dense(s_obj, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_obj);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	int exact;
	if (TYPEOF(s_exact) != LGLSXP || LENGTH(s_exact) < 1 ||
	    (exact = LOGICAL(s_exact)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "exact", "TRUE", "FALSE");

	char ct = 'C';
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	int checkDN;
	if (TYPEOF(s_checkDN) != LGLSXP || LENGTH(s_checkDN) < 1 ||
	    (checkDN = LOGICAL(s_checkDN)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "checkDN", "TRUE", "FALSE");

	int ans_ = dense_is_symmetric(s_obj, valid[ivalid], exact, ct, checkDN);
	SEXP ans = ScalarLogical(ans_);
	UNPROTECT(1);
	return ans;
}

int dense_is_triangular(SEXP obj, const char *class, int upper)
{
	if (class[1] == 't') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (upper == NA_LOGICAL || (upper != 0) == (ul == 'U'))
			return (ul == 'U') ? 1 : -1;
		else if (dense_is_diagonal(obj, class))
			return (ul == 'U') ? -1 : 1;
		else
			return 0;
	}

	if (class[1] == 's') {
		if (!dense_is_diagonal(obj, class))
			return 0;
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (upper == NA_LOGICAL)
			return (ul == 'U') ? 1 : -1;
		else
			return (upper != 0) ? 1 : -1;
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return (upper != 0) ? 1 : -1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

#define IT_LOOP(_CTYPE_, _PTR_, _NOTZERO_) \
	do { \
		_CTYPE_ *px; \
		if (upper == NA_LOGICAL) { \
			px = _PTR_(x); \
			for (j = 0; j < n; px += (++j)) { \
				px += 1; \
				for (i = j + 1; i < n; ++i, px += 1) { \
					if (_NOTZERO_(*px)) { \
						j = n; \
						break; \
					} \
				} \
			} \
			if (j == n) \
				return  1; \
			px = _PTR_(x); \
			for (j = 0; j < n; px += n - (++j)) { \
				for (i = 0; i < j; ++i, px += 1) { \
					if (_NOTZERO_(*px)) { \
						j = n; \
						break; \
					} \
				} \
				px += 1; \
			} \
			if (j == n) \
				return -1; \
			return 0; \
		} else if (upper != 0) { \
			px = _PTR_(x); \
			for (j = 0; j < n; px += (++j)) { \
				px += 1; \
				for (i = j + 1; i < n; ++i, px += 1) \
					if (_NOTZERO_(*px)) \
						return 0; \
			} \
			return  1; \
		} else { \
			px = _PTR_(x); \
			for (j = 0; j < n; px += n - (++j)) { \
				for (i = 0; i < j; ++i, px += 1) \
					if (_NOTZERO_(*px)) \
						return 0; \
				px += 1; \
			} \
			return -1; \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
		IT_LOOP(int, LOGICAL, NOTZERO_PATTERN);
		break;
	case 'l':
		IT_LOOP(int, LOGICAL, NOTZERO_LOGICAL);
		break;
	case 'i':
		IT_LOOP(int, INTEGER, NOTZERO_INTEGER);
		break;
	case 'd':
		IT_LOOP(double, REAL, NOTZERO_REAL);
		break;
	case 'z':
		IT_LOOP(Rcomplex, COMPLEX, NOTZERO_COMPLEX);
		break;
	default:
		break;
	}

#undef IT_LOOP

	return 0;
}

/* isTriangular(<denseMatrix>, upper) */
SEXP R_dense_is_triangular(SEXP s_obj, SEXP s_upper)
{
	if (TYPEOF(s_obj) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_obj = matrix_as_dense(s_obj, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_obj);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	if (TYPEOF(s_upper) != LGLSXP || LENGTH(s_upper) < 1)
		error(_("'%s' must be %s or %s or %s"), "upper", "TRUE", "FALSE", "NA");
	int upper = LOGICAL(s_upper)[0];

	int ans_ = dense_is_triangular(s_obj, valid[ivalid], upper);
	SEXP ans = allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = ans_ != 0;
	if (upper == NA_LOGICAL && ans_ != 0) {
		PROTECT(ans);
		static
		SEXP kindSym = NULL;
		SEXP kindVal = PROTECT(mkString((ans_ > 0) ? "U" : "L"));
		if (!kindSym)
			kindSym = install("kind");
		setAttrib(ans, kindSym, kindVal);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

int dense_is_diagonal(SEXP obj, const char *class)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		return 0;
	if (n <= 1)
		return 1;

	char ul = 'U';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p';

#define ID_LOOP(_CTYPE_, _PTR_, _NOTZERO_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (_NOTZERO_(*px)) \
						return 0; \
					px += 1; \
				} \
				px += 1; \
				for (i = j + 1; i < n; ++i) { \
					if (_NOTZERO_(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (_NOTZERO_(*px)) \
						return 0; \
					px += 1; \
				} \
				px += (!packed) ? n - j : 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				px += (!packed) ? j + 1 : 1; \
				for (i = j + 1; i < n; ++i) { \
					if (_NOTZERO_(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} \
		return 1; \
	} while (0)

	switch (class[0]) {
	case 'n':
		ID_LOOP(int, LOGICAL, NOTZERO_PATTERN);
		break;
	case 'l':
		ID_LOOP(int, LOGICAL, NOTZERO_LOGICAL);
		break;
	case 'i':
		ID_LOOP(int, INTEGER, NOTZERO_INTEGER);
		break;
	case 'd':
		ID_LOOP(double, REAL, NOTZERO_REAL);
		break;
	case 'z':
		ID_LOOP(Rcomplex, COMPLEX, NOTZERO_COMPLEX);
		break;
	default:
		break;
	}

#undef ID_LOOP

	return 0;
}

/* isDiagonal(<denseMatrix>) */
SEXP R_dense_is_diagonal(SEXP s_obj)
{
	if (TYPEOF(s_obj) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_obj = matrix_as_dense(s_obj, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_obj);
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	int ans_ = dense_is_diagonal(s_obj, valid[ivalid]);
	SEXP ans = ScalarLogical(ans_ != 0);
	UNPROTECT(1);
	return ans;
}

#define CAST_PATTERN(_X_) (_X_ != 0)
#define CAST_LOGICAL(_X_) (_X_ != 0)
#define CAST_INTEGER(_X_)  _X_
#define CAST_REAL(_X_)     _X_
#define CAST_COMPLEX(_X_)  _X_

#define SUM_CASES \
do { \
	switch (class[0]) { \
	case 'n': \
		if (mean) \
		SUM_LOOP(int, LOGICAL, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_PATTERN, CAST_PATTERN, \
		         INCREMENT_REAL, INCREMENT_REAL, SCALE2_REAL); \
		else \
		SUM_LOOP(int, LOGICAL, int, INTEGER, \
		         0, 1, NA_INTEGER, ISNA_PATTERN, CAST_PATTERN, \
		         INCREMENT_INTEGER, INCREMENT_INTEGER, SCALE2_REAL); \
		break; \
	case 'l': \
		if (mean) \
		SUM_LOOP(int, LOGICAL, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_LOGICAL, CAST_LOGICAL, \
		         INCREMENT_REAL, INCREMENT_REAL, SCALE2_REAL); \
		else \
		SUM_LOOP(int, LOGICAL, int, INTEGER, \
		         0, 1, NA_INTEGER, ISNA_LOGICAL, CAST_LOGICAL, \
		         INCREMENT_INTEGER, INCREMENT_INTEGER, SCALE2_REAL); \
		break; \
	case 'i': \
		SUM_LOOP(int, INTEGER, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_INTEGER, CAST_INTEGER, \
		         INCREMENT_REAL, INCREMENT_REAL, SCALE2_REAL); \
		break; \
	case 'd': \
		SUM_LOOP(double, REAL, double, REAL, \
		         0.0, 1.0, NA_REAL, ISNA_REAL, CAST_REAL, \
		         INCREMENT_REAL, INCREMENT_REAL, SCALE2_REAL); \
		break; \
	case 'z': \
		SUM_LOOP(Rcomplex, COMPLEX, Rcomplex, COMPLEX, \
		         Matrix_zzero, Matrix_zone, Matrix_zna, ISNA_COMPLEX, CAST_COMPLEX, \
		         INCREMENT_COMPLEX_ID, INCREMENT_COMPLEX_CJ, SCALE2_COMPLEX); \
		break; \
	default: \
		break; \
	} \
} while (0)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void dense_colsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char di, int narm, int mean,
                  SEXP ans)
{
	int i, j, count = -1, narm_ = narm && mean && class[0] != 'n',
		packed = class[2] == 'p';

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, \
	             _ZERO_, _ONE_, _NA_, _ISNA_, \
	             _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, _SCALE2_) \
	do { \
		_CTYPE0_ *px0 = _PTR0_(x); \
		_CTYPE1_ *px1 = _PTR1_(ans), tmp; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				*px1 = _ZERO_; \
				SUM_KERNEL(for (i = 0; i < m; ++i), _NA_, _ISNA_, \
				           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, \
				           _SCALE2_); \
				px1 += 1; \
			} \
		} else if (di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					*px1 = _ZERO_; \
					SUM_KERNEL(for (i = 0; i <= j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, \
					           _SCALE2_); \
					if (!packed) \
						px0 += n - j - 1; \
					px1 += 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px0 += j; \
					*px1 = _ZERO_; \
					SUM_KERNEL(for (i = j; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, \
					           _SCALE2_); \
					px1 += 1; \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					*px1 = _ONE_; \
					SUM_KERNEL(for (i = 0; i < j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, \
					           _SCALE2_); \
					px0 += 1; \
					if (!packed) \
						px0 += n - j - 1; \
					px1 += 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px0 += j; \
					px0 += 1; \
					*px1 = _ONE_; \
					SUM_KERNEL(for (i = j + 1; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, \
					           _SCALE2_); \
					px1 += 1; \
				} \
			} \
		} \
	} while (0)

#define SUM_KERNEL(_FOR_, _NA_, _ISNA_, \
	               _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, _SCALE2_) \
	do { \
		if (mean) \
			count = m; \
		_FOR_ { \
			if (_ISNA_(*px0)) { \
				if (!narm) \
					*px1 = _NA_; \
				else if (narm_) \
					--count; \
			} else { \
				tmp = _CAST_(*px0); \
				_INCREMENT_ID_((*px1), tmp); \
			} \
			px0 += 1; \
		} \
		if (mean) \
			_SCALE2_((*px1), count); \
	} while (0)

	SUM_CASES;

#undef SUM_LOOP
#undef SUM_KERNEL

	return;
}

static
void dense_rowsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char di, int narm, int mean,
                  SEXP ans)
{
	int i, j, *count = NULL, narm_ = narm && mean && class[0] != 'n',
		packed = class[2] == 'p', sy = class[1] == 's', he = sy && ct == 'C';

	if (narm_) {
		Matrix_Calloc(count, m, int);
		for (i = 0; i < m; ++i)
			count[i] = n;
	}

#define SUM_LOOP(_CTYPE0_, _PTR0_, _CTYPE1_, _PTR1_, \
	             _ZERO_, _ONE_, _NA_, _ISNA_, \
	             _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_, _SCALE2_) \
	do { \
		_CTYPE0_ *px0 = _PTR0_(x); \
		_CTYPE1_ *px1 = _PTR1_(ans), tmp = (di == 'N') ? _ZERO_ : _ONE_; \
		for (i = 0; i < m; ++i) \
			px1[i] = tmp; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i), _NA_, _ISNA_, \
				           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_); \
		} else if (class[1] == 's' || di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i <= j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_); \
					if (!packed) \
						px0 += n - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px0 += j; \
					SUM_KERNEL(for (i = j; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_); \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i < j; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_); \
					px0 += 1; \
					if (!packed) \
						px0 += n - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px0 += j; \
					px0 += 1; \
					SUM_KERNEL(for (i = j + 1; i < n; ++i), _NA_, _ISNA_, \
					           _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_); \
				} \
			} \
		} \
		if (mean) { \
			if (narm_) \
				for (i = 0; i < m; ++i) \
					_SCALE2_(px1[i], count[i]); \
			else \
				for (i = 0; i < m; ++i) \
					_SCALE2_(px1[i], n); \
		} \
	} while (0)

#define SUM_KERNEL(_FOR_, _NA_, _ISNA_, \
	               _CAST_, _INCREMENT_ID_, _INCREMENT_CJ_) \
	do { \
		_FOR_ { \
			if (_ISNA_(*px0)) { \
				if (!narm) { \
					px1[i] = _NA_; \
					if (sy && i != j) \
					px1[j] = _NA_; \
				} else if (narm_) { \
					--count[i]; \
					if (sy && i != j) \
					--count[j]; \
				} \
			} else { \
				tmp = _CAST_(*px0); \
				_INCREMENT_ID_(px1[i], tmp); \
				if (sy && i != j) { \
				if (he) \
				_INCREMENT_CJ_(px1[j], tmp); \
				else \
				_INCREMENT_ID_(px1[j], tmp); \
				} \
			} \
			px0 += 1; \
		} \
	} while (0)

	SUM_CASES;

#undef SUM_LOOP
#undef SUM_KERNEL

	if (narm_)
		Matrix_Free(count, m);
	return;
}

SEXP dense_marginsum(SEXP obj, const char *class, int mg, int narm, int mean)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (mg == 0) ? m : n;

	SEXP ans = PROTECT(allocVector(SUM_TYPEOF(class[0]), r)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));

	SEXP dimnames = (class[1] == 's')
		? get_symmetrized_DimNames(obj, -1)
		: GET_SLOT(obj, Matrix_DimNamesSym),
		marnames = VECTOR_ELT(dimnames, mg);
	if (marnames != R_NilValue) {
		PROTECT(marnames);
		setAttrib(ans, R_NamesSymbol, marnames);
		UNPROTECT(1); /* marnames */
	}

	char ul = 'U', ct = 'C', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = GET_SLOT(obj, Matrix_transSym);
			ct = CHAR(STRING_ELT(trans, 0))[0];
		}
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}

	if (mg == 0 || class[1] == 's')
		dense_rowsum(x, class, m, n, ul, ct, di, narm, mean, ans);
	else
		dense_colsum(x, class, m, n, ul, ct, di, narm, mean, ans);

	UNPROTECT(2); /* x, ans */
	return ans;
}

/* (row|col)(Sums|Means)(<denseMatrix>, na.rm=) */
SEXP R_dense_marginsum(SEXP s_obj, SEXP s_margin, SEXP s_narm, SEXP s_mean)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	int mg;
	if (TYPEOF(s_margin) != INTSXP || LENGTH(s_margin) < 1 ||
	    ((mg = INTEGER(s_margin)[0]) != 0 && mg != 1))
		error(_("'%s' must be %d or %d"), "margin", 0, 1);

	int narm;
	if (TYPEOF(s_narm) != LGLSXP || LENGTH(s_narm) < 1 ||
	    (narm = LOGICAL(s_narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	int mean;
	if (TYPEOF(s_mean) != LGLSXP || LENGTH(s_mean) < 1 ||
	    (mean = LOGICAL(s_mean)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "mean", "TRUE", "FALSE");

	return dense_marginsum(s_obj, valid[ivalid], mg, narm, mean);
}

#undef SUM_CASES
#undef SUM_TYPEOF

#define TRY_INCREMENT(_LABEL_) \
	do { \
		if ((s >= 0) \
			? ( t <= INT_FAST64_MAX - s) \
			: (-t <= s - INT_FAST64_MIN)) { \
			s += t; \
			t = 0; \
			count = 0; \
		} else { \
			over = 1; \
			goto _LABEL_; \
		} \
	} while (0)

#define LONGDOUBLE_AS_DOUBLE(v) \
	(v > DBL_MAX) ? R_PosInf : ((v < -DBL_MAX) ? R_NegInf : (double) v);

SEXP dense_sum(SEXP obj, const char *class, int narm)
{
	SEXP ans;

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = 'U', ct = 'C', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = GET_SLOT(obj, Matrix_transSym);
			ct = CHAR(STRING_ELT(trans, 0))[0];
		}
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j,
		packed = class[2] == 'p', sy = class[1] == 's', he = sy && ct == 'C';

#define SUM_LOOP \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (sy || di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i <= j; ++i)); \
					if (!packed) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px += j; \
					SUM_KERNEL(for (i = j; i < m; ++i)); \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					SUM_KERNEL(for (i = 0; i < j; ++i)); \
					px += 1; \
					if (!packed) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px += j; \
					px += 1; \
					SUM_KERNEL(for (i = j + 1; i < m; ++i)); \
				} \
			} \
		} \
	} while (0)

	if (class[0] == 'n') {
		int *px = LOGICAL(x);
		int_fast64_t s = (di == 'N') ? 0LL : n;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (*px != 0) \
					s += (sy && i != j) ? 2 : 1; \
				px += 1; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		if (s <= INT_MAX) {
			ans = allocVector(INTSXP, 1);
			INTEGER(ans)[0] = (int) s;
		} else {
			ans = allocVector(REALSXP, 1);
			REAL(ans)[0] = (double) s;
		}
		return ans;
	}

	if (!narm && (class[0] == 'l' || class[0] == 'i')) {
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (*px == NA_INTEGER) { \
					ans = allocVector(INTSXP, 1); \
					INTEGER(ans)[0] = NA_INTEGER; \
					return ans; \
				} \
				px += 1; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

	}

	if (class[0] == 'z') {
		Rcomplex *px = COMPLEX(x);
		long double zr = (di == 'N') ? 0.0L : n, zi = 0.0L;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
					zr += (sy && i != j) ?                2.0L * (*px).r  : (*px).r; \
					zi += (sy && i != j) ? ((he) ? 0.0L : 2.0L * (*px).i) : (*px).i; \
				} \
				px += 1; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		ans = allocVector(CPLXSXP, 1);
		COMPLEX(ans)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
		COMPLEX(ans)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
	} else if (class[0] == 'd') {
		double *px = REAL(x);
		long double zr = (di == 'N') ? 0.0L : n;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && ISNAN(*px))) \
					zr += (sy && i != j) ? 2.0L * *px : *px; \
				px += 1; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		ans = allocVector(REALSXP, 1);
		REAL(ans)[0] = LONGDOUBLE_AS_DOUBLE(zr);
	} else {
		int *px = (class[0] == 'i') ? INTEGER(x) : LOGICAL(x);
		int_fast64_t s = (di == 'N') ? 0LL : n, t = 0LL;
		unsigned int count = 0;
		int over = 0;

#define SUM_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && *px == NA_INTEGER)) { \
					int d = (sy && i != j) ? 2 : 1; \
					if (count > UINT_MAX - d) \
						TRY_INCREMENT(ifover); \
					t += (d == 2) ? 2LL * *px : *px; \
					count += d; \
				} \
				px += 1; \
			} \
		} while (0)

		SUM_LOOP;

#undef SUM_KERNEL

		TRY_INCREMENT(ifover);
	ifover:
		if (over) {
			long double zr = (di == 'N') ? 0.0L : n; /* FIXME: wasteful */
			px = (class[0] == 'i') ? INTEGER(x) : LOGICAL(x);

#define SUM_KERNEL(_FOR_) \
			do { \
				_FOR_ { \
					if (!(narm && *px == NA_INTEGER)) \
						zr += (sy && i != j) ? 2.0L * *px : *px; \
					px += 1; \
				} \
			} while (0)

			SUM_LOOP;

#undef SUM_KERNEL

			ans = allocVector(REALSXP, 1);
			REAL(ans)[0] = LONGDOUBLE_AS_DOUBLE(zr);
		} else if (s > INT_MIN && s <= INT_MAX) {
			ans = allocVector(INTSXP, 1);
			INTEGER(ans)[0] = (int) s;
		} else {
			ans = allocVector(REALSXP, 1);
			REAL(ans)[0] = (double) s;
		}
	}

#undef SUM_LOOP

	return ans;
}

/* sum(<denseMatrix>, na.rm=) */
SEXP R_dense_sum(SEXP s_obj, SEXP s_narm)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	int narm;
	if (TYPEOF(s_narm) != LGLSXP || LENGTH(s_narm) < 1 ||
	    (narm = LOGICAL(s_narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return dense_sum(s_obj, valid[ivalid], narm);
}

SEXP dense_prod(SEXP obj, const char *class, int narm)
{
	SEXP ans = PROTECT(allocVector((class[0] == 'z') ? CPLXSXP : REALSXP, 1));

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = 'U', ct = 'C', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (class[1] == 's' && class[0] == 'z') {
			SEXP trans = GET_SLOT(obj, Matrix_transSym);
			ct = CHAR(STRING_ELT(trans, 0))[0];
		}
		if (class[1] == 't') {
			SEXP diag = GET_SLOT(obj, Matrix_diagSym);
			di = CHAR(STRING_ELT(diag, 0))[0];
		}
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j,
		packed = class[2] == 'p', sy = class[1] == 's', he = sy && ct == 'C';
	long double zr = 1.0L, zi = 0.0L;

#define PROD_LOOP \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				PROD_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					PROD_KERNEL(for (i = 0; i <= j; ++i)); \
					if (!packed) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (!packed) \
						px += j; \
					PROD_KERNEL(for (i = j; i < m; ++i)); \
				} \
			} \
		} else if (di == 'N') { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					PROD_KERNEL(for (i = 0; i <= j; ++i)); \
					if (!packed) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					if (!packed) \
						px += j; \
					PROD_KERNEL(for (i = j; i < m; ++i)); \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					PROD_KERNEL(for (i = 0; i < j; ++i)); \
					px += 1; \
					if (!packed) \
						px += m - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					if (j == 1) { zr *= 0.0L; zi *= 0.0L; } \
					if (!packed) \
						px += j; \
					px += 1; \
					PROD_KERNEL(for (i = j + 1; i < m; ++i)); \
				} \
			} \
		} \
	} while (0)

	if (class[0] == 'n') {
		int *px = LOGICAL(x);
		if (class[1] == 't')
			REAL(ans)[0] = (n > 1 || (n == 1 && *px == 0)) ? 0.0 : 1.0;
		else {

#define PROD_KERNEL(_FOR_) \
			do { \
				_FOR_ { \
					if (*px == 0) { \
						REAL(ans)[0] = 0.0; \
						UNPROTECT(1); /* ans */ \
						return ans; \
					} \
					px += 1; \
				} \
			} while (0)

			PROD_LOOP;

#undef PROD_KERNEL

			REAL(ans)[0] = 1.0;
		}
		UNPROTECT(1); /* ans */
		return ans;
	}

	if (class[0] == 'z') {
		Rcomplex *px = COMPLEX(x);
		long double zr0, zi0;

#define PROD_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
					zr0 = zr; zi0 = zi; \
					zr = zr0 * (*px).r - zi0 * (*px).i; \
					zi = zi0 * (*px).r + zr0 * (*px).i; \
					if (sy && i != j) { \
					zr0 = zr; zi0 = zi; \
					if (he) { \
					zr = zr0 * (*px).r + zi0 * (*px).i; \
					zi = zi0 * (*px).r - zr0 * (*px).i; \
					} else { \
					zr = zr0 * (*px).r - zi0 * (*px).i; \
					zi = zi0 * (*px).r + zr0 * (*px).i; \
					} \
					} \
				} \
				px += 1; \
			} \
		} while (0)

		PROD_LOOP;

#undef PROD_KERNEL

	} else if (class[0] == 'd') {
		double *px = REAL(x);

#define PROD_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (!(narm && ISNAN(*px))) \
					zr *= (sy && i != j) ? (long double) *px * *px : *px; \
				px += 1; \
			} \
		} while (0)

		PROD_LOOP;

#undef PROD_KERNEL

	} else {
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);

#define PROD_KERNEL(_FOR_) \
		do { \
			_FOR_ { \
				if (*px != NA_INTEGER) \
					zr *= (sy && i != j) ? (long double) *px * *px : *px; \
				else if (!narm) \
					zr *= NA_REAL; \
				px += 1; \
			} \
		} while (0)

		PROD_LOOP;

#undef PROD_KERNEL

	}

#undef PROD_LOOP

	if (class[0] == 'z') {
		COMPLEX(ans)[0].r = LONGDOUBLE_AS_DOUBLE(zr);
		COMPLEX(ans)[0].i = LONGDOUBLE_AS_DOUBLE(zi);
	} else
		   REAL(ans)[0]   = LONGDOUBLE_AS_DOUBLE(zr);
	UNPROTECT(1); /* ans */
	return ans;
}

/* prod(<denseMatrix>, na.rm=) */
SEXP R_dense_prod(SEXP s_obj, SEXP s_narm)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int ivalid = R_check_class_etc(s_obj, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(s_obj, __func__);

	int narm;
	if (TYPEOF(s_narm) != LGLSXP || LENGTH(s_narm) < 1 ||
	    (narm = LOGICAL(s_narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return dense_prod(s_obj, valid[ivalid], narm);
}

#undef TRY_INCREMENT
#undef LONGDOUBLE_AS_DOUBLE
