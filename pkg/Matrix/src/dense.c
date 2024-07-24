#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "dense.h"

SEXP dense_band(SEXP from, const char *class, int a, int b)
{
	int packed = class[2] == 'p';

	SEXP dim = GET_SLOT(from, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
	if ((m == 0 || n == 0 || (a <= 1 - m && b >= n - 1)) &&
	    (m != n || n > 1 || class[1] == 't'))
		return from;

	int ge, sy, tr;
	tr = class[1] == 't' || (m == n && (a >= 0 || b <= 0));
	sy = !tr && class[1] == 's' && a == -b;
	ge = !tr && !sy;

	char ul0 = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul0 = CHAR(STRING_ELT(uplo, 0))[0];
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		/* Be fast if band contains entire triangle */
		if ((ul0 == 'U') ? (a <= 0 && b >= n - 1) : (a <= 1 - m && b >= 0))
			return from;
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = CHAR(STRING_ELT(diag, 0))[0];
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

	char ul1 = (tr && class[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ul0;
	if (ul1 != '\0' && ul1 != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (ct  != '\0' && ct  != 'C' && sy) {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (di  != '\0' && di  != 'N' && tr && a <= 0 && b >= 0) {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (!packed || ge) ? (R_xlen_t) m * n : (R_xlen_t) PACKED_LENGTH((size_t) n)));
	SET_SLOT(to, Matrix_xSym, x1);

	size_t
		m_ = (size_t) m,
		n_ = (size_t) n,
		a_ = (size_t) ((int_fast64_t) m + a),
		b_ = (size_t) ((int_fast64_t) m + b);

#define BAND(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (ge && class[1] != 'g') { \
		if (!packed) \
			c##NAME(force2)(px1, px0, n_, ul0, ct, di); \
		else \
			c##NAME( pack1)(px1, px0, n_, ul0, ct, di); \
		px0 = NULL; \
		packed = 0; \
		} else if (tr && class[1] == 's' && ul0 != ul1) { \
		if (!packed) \
			c##NAME(trans2)(px1, px0, m_, n_ , ct); \
		else \
			c##NAME(trans1)(px1, px0, n_, ul0, ct); \
		px0 = NULL; \
		ul0 = ul1; \
		} \
		if (!packed) \
			c##NAME( band2)(px1, px0, m_, n_ , a_, b_); \
		else \
			c##NAME( band1)(px1, px0, n_, ul0, a_, b_); \
	} while (0)

	SWITCH4(class[0], BAND);

#undef BAND

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
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

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

	s_from = dense_band(s_from, class, a, b);
	UNPROTECT(1);
	return s_from;
}

SEXP dense_diag_get(SEXP obj, const char *class, int names)
{
	int packed = class[2] == 'p';

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
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

	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		ans = PROTECT(allocVector(TYPEOF(x), r));

	size_t m_ = (size_t) m, n_ = (size_t) n, r_ = (size_t) r;

#define DIAG(c) \
	do { \
		c##TYPE *pa = c##PTR(ans), *px = c##PTR(x); \
		if (di != '\0' && di != 'N') \
			for (int j = 0; j < r; ++j) \
				*(pa++) = c##UNIT; \
		else if (!packed) \
			c##NAME(copy2)(r_, pa, 1, px, m_ + 1); \
		else if (ul == 'U') \
			c##NAME(copy1)(n_, pa, 1, 0, 0, px, 2 , 1, 0); \
		else \
			c##NAME(copy1)(n_, pa, 1, 0, 0, px, n_, 1, 1); \
	} while (0)

	SWITCH4(class[0], DIAG);

#undef DIAG

	if (class[0] == 'n')
		naToUnit(ans);
	if (class[0] == 'z' && class[1] == 's' && ct == 'C')
		zvreal(COMPLEX(ans), r_);

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
			         (rn == cn || equalString(rn, cn, r)))
				setAttrib(ans, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(2); /* x, ans */
	return ans;
}

/* diag(<denseMatrix>, names=) */
SEXP R_dense_diag_get(SEXP s_obj, SEXP s_names)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int names;
	if (TYPEOF(s_names) != LGLSXP || LENGTH(s_names) < 1 ||
	    (names = LOGICAL(s_names)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "names", "TRUE", "FALSE");

	return dense_diag_get(s_obj, class, names);
}

SEXP dense_diag_set(SEXP from, const char *class, SEXP value, int new)
{
	int packed = class[2] == 'p';

	SEXP to = PROTECT(newObject(class));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (m < n) ? m : n;
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = '\0';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (new) {
		x = duplicateVector(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
	}
	SET_SLOT(to, Matrix_xSym, x);

	size_t m_ = (size_t) m, n_ = (size_t) n, r_ = (size_t) r,
		v_ = (LENGTH(value) == r) ? 1 : 0;

#define DIAG(c) \
	do { \
		c##TYPE *px = c##PTR(x), *pvalue = c##PTR(value); \
		if (!packed) \
			c##NAME(copy2)(r_, px, m_ + 1, pvalue, v_); \
		else if (ul == 'U') \
			c##NAME(copy1)(n_, px, 2 , 1, 0, pvalue, v_, 0, 0); \
		else  \
			c##NAME(copy1)(n_, px, n_, 1, 1, pvalue, v_, 0, 0); \
	} while (0)

	SWITCH4(class[0], DIAG);

#undef DIAG

	UNPROTECT(2); /* x, to */
	return to;
}

/* diag(<denseMatrix>) <- value */
SEXP R_dense_diag_set(SEXP s_from, SEXP s_value)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);
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
	if (XLENGTH(s_value) != 1 && XLENGTH(s_value) != r)
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
				class = Matrix_class(s_from, valid_dense, 6, __func__);
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
		class = Matrix_class(s_from, valid_dense, 6, __func__);
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

	char ul = '\0';
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

	size_t m_ = (size_t) m, n_ = (size_t) n;

#define TRANS(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed) \
			c##NAME(trans2)(px1, px0, m_, n_, ct); \
		else \
			c##NAME(trans1)(px1, px0, n_, ul, ct); \
	} while (0)

	SWITCH4((class[0] == 'c') ? 'd' : class[0], TRANS);

#undef TRANS

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* t(<denseMatrix>) */
SEXP R_dense_transpose(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);

	char ct;
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	return dense_transpose(s_from, class, ct);
}

SEXP dense_force_symmetric(SEXP from, const char *class, char ul, char ct)
{
	char ul0 = '\0', ul1 = 'U', ct0 = '\0', ct1 = (class[0] == 'z') ? 'C' : '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(from, Matrix_uploSym);
		ul0 = ul1 = CHAR(STRING_ELT(uplo, 0))[0];
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct0 = ct1 = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}
	if (ul != '\0')
		ul1 = ul;
	if (ct != '\0' && class[0] == 'z')
		ct1 = ct;

	if (class[1] == 's' && ul0 == ul1 && (ct0 == ct1 || ct1 == 'C'))
		return from;

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

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));

	if (class[1] == 'g' || (class[0] == 't' && ul0 == ul1 && di == 'N'))
		SET_SLOT(to, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
		SET_SLOT(to, Matrix_xSym, x1);

		size_t n_ = (size_t) n;

#define FORCE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			if (!packed) \
				c##NAME(force2)(px1,  px0, n_, ul0, ct0, di); \
			else if (ul0 == ul1) \
				c##NAME(force1)(px1,  px0, n_, ul0, ct0, di); \
			else { \
				c##NAME(trans1)(px1,  px0, n_, ul0, ct0); \
				c##NAME(force1)(px1, NULL, n_, ul1, ct0, di); \
			} \
		} while (0)

		SWITCH4(class[0], FORCE);

#undef FORCE

		UNPROTECT(1); /* x1 */
	}

	UNPROTECT(2); /* x0, to */
	return to;
}

/* forceSymmetric(<denseMatrix>, uplo, trans) */
SEXP R_dense_force_symmetric(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

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

	return dense_force_symmetric(s_from, class, ul, ct);
}

SEXP dense_symmpart(SEXP from, const char *class, char ct)
{
	/* defined in ./coerce.c : */
	SEXP dense_as_kind(SEXP, const char *, char, int);
	PROTECT(from = dense_as_kind(from, class, ',', 0));

	char ct0 = ct, ct1 = ct;
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct0 = CHAR(STRING_ELT(trans, 0))[0];
	}

	if (class[1] == 's' && ct0 == ct1) {
		UNPROTECT(1);
		return from;
	}

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = 's';
	cl[2] = (packed) ? 'p' : 'y';
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error((ct1 == 'C')
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

	char ul = '\0';
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

	char di = '\0';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
	SET_SLOT(to, Matrix_xSym, x1);

	if (class[1] == 's') {

		/* Symmetric part of Hermitian matrix is real part */
		/* Hermitian part of symmetric matrix is real part */

		zvreal(COMPLEX(x1), (size_t) XLENGTH(x1));

	} else {

		int i, j;

#define SPART(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *pu0 = px0, *pl0 = px0; \
			c##TYPE *px1 = c##PTR(x1), *pu1 = px1, *pl1 = px1; \
			if (!packed) \
				memset(px1, 0, sizeof(c##TYPE) * (size_t) XLENGTH(x1)); \
			if (class[1] == 'g') { \
				if (ct1 == 'C') \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##INCREMENT_CONJ(*pu1, *pl0); \
							c##MULTIPLY(*pu1, 0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl0 += n; \
						} \
						c##ASSIGN_PROJ_REAL(*pu1, *pu0); \
						pu0 += n - j; \
						pu1 += n - j; \
						pl0 = px0 + j + 1; \
					} \
				else \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##INCREMENT_IDEN(*pu1, *pl0); \
							c##MULTIPLY(*pu1, 0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl0 += n; \
						} \
						c##ASSIGN_IDEN(*pu1, *pu0); \
						pu0 += n - j; \
						pu1 += n - j; \
						pl0 = px0 + j + 1; \
					} \
			} else { \
				if (ul == 'U') \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##MULTIPLY(*pu1, 0.5); \
							pu0 += 1; \
							pu1 += 1; \
						} \
						if (di != 'N') \
							c##SET_UNIT(*pu1); \
						else if (ct1 == 'C') \
							c##ASSIGN_PROJ_REAL(*pu1, *pu0); \
						else \
							c##ASSIGN_IDEN(*pu1, *pu0); \
						pu0 += 1; \
						pu1 += 1; \
						if (!packed) { \
							pu0 += n - j - 1; \
							pu1 += n - j - 1; \
						} \
					} \
				else \
					for (j = 0; j < n; ++j) { \
						if (!packed) { \
							pl0 += j; \
							pl1 += j; \
						} \
						if (di != 'N') \
							c##SET_UNIT(*pl1); \
						else if (ct1 == 'C') \
							c##ASSIGN_PROJ_REAL(*pl1, *pl0); \
						else \
							c##ASSIGN_IDEN(*pl1, *pl0); \
						pl0 += 1; \
						pl1 += 1; \
						for (i = j + 1; i < n; ++i) { \
							c##ASSIGN_IDEN(*pl1, *pl0); \
							c##MULTIPLY(*pl1, 0.5); \
							pl0 += 1; \
							pl1 += 1; \
						} \
					} \
			} \
		} while (0)

		SWITCH2(class[0], SPART);

#undef SPART

	}

	UNPROTECT(4); /* x1, x0, to, from */
	return to;
}

/* symmpart(<denseMatrix>, trans) */
SEXP R_dense_symmpart(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ct;
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	return dense_symmpart(s_from, class, ct);
}

SEXP dense_skewpart(SEXP from, const char *class, char ct)
{
	/* defined in ./coerce.c : */
	SEXP dense_as_kind(SEXP, const char *, char, int);
	PROTECT(from = dense_as_kind(from, class, ',', 0));

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 's') ? 's' : 'g';
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

	char ul = '\0';
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

	char di = '\0';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1;

	if (class[1] == 's') {

		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		/* Skew-Hermitian part of symmetric matrix is imaginary part */

		PROTECT(x1 = allocVector(TYPEOF(x0), XLENGTH(x0)));
		if (class[0] == 'z') {
			Rcomplex *px1 = COMPLEX(x1);
			if (ct0 == ct1)
			memset(px1, 0, sizeof(Rcomplex) * (size_t) XLENGTH(x1));
			else {
			Rcomplex *px0 = COMPLEX(x0);
			for (R_xlen_t k = 0, kend = XLENGTH(x1); k < kend; ++k) {
				zASSIGN_PROJ_IMAG(*px1, *px0);
				px0 += 1;
				px1 += 1;
			}
			}
		} else {
			double *px1 = REAL(x1);
			memset(px1, 0, sizeof(double) * (size_t) XLENGTH(x1));
		}

	} else {

		if (class[0] == 't' && (int_fast64_t) n * n > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");
		PROTECT(x1 = allocVector(TYPEOF(x0), (R_xlen_t) n * n));

		int i, j;

#define SPART(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *pu0 = px0, *pl0 = px0; \
			c##TYPE *px1 = c##PTR(x1), *pu1 = px1, *pl1 = px1; \
			if (class[1] == 'g') { \
				if (ct1 == 'C') \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##DECREMENT_CONJ(*pu1, *pl0); \
							c##ASSIGN_CONJ(*pl1, *pu1); \
							c##MULTIPLY(*pu1,  0.5); \
							c##MULTIPLY(*pl1, -0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl0 += n; \
							pl1 += n; \
						} \
						c##ASSIGN_PROJ_IMAG(*pu1, *pu0); \
						pu0 += n - j; \
						pu1 += n - j; \
						pl0 = px0 + j + 1; \
						pl1 = px1 + j + 1; \
					} \
				else \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##DECREMENT_IDEN(*pu1, *pl0); \
							c##ASSIGN_IDEN(*pl1, *pu1); \
							c##MULTIPLY(*pu1,  0.5); \
							c##MULTIPLY(*pl1, -0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl0 += n; \
							pl1 += n; \
						} \
						c##SET_ZERO(*pu1); \
						pu0 += n - j; \
						pu1 += n - j; \
						pl0 = px0 + j + 1; \
						pl1 = px1 + j + 1; \
					} \
			} else { \
				if (ul == 'U') \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##ASSIGN_IDEN(*pl1, *pu0); \
							c##MULTIPLY(*pu1,  0.5); \
							c##MULTIPLY(*pl1, -0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl1 += n; \
						} \
						if (di != 'N' || ct1 != 'C') \
							c##SET_ZERO(*pu1); \
						else \
							c##ASSIGN_PROJ_IMAG(*pu1, *pu0); \
						pu0 += 1; \
						pu1 += 1; \
						pl1 += n; \
						if (!packed) \
						pu0 += n - j - 1; \
						pu1 += n - j - 1; \
						pl1 = px1 + j + 1; \
					} \
				else \
					for (j = 0; j < n; ++j) { \
						if (!packed) \
						pl0 += j; \
						pl1 += j; \
						pu1 = pl1; \
						if (di != 'N' || ct1 != 'C') \
							c##SET_ZERO(*pl1); \
						else \
							c##ASSIGN_PROJ_IMAG(*pl1, *pl0); \
						pl0 += 1; \
						pl1 += 1; \
						pu1 += n; \
						for (i = j + 1; i < n; ++i) { \
							c##ASSIGN_IDEN(*pl1, *pl0); \
							c##ASSIGN_IDEN(*pu1, *pl0); \
							c##MULTIPLY(*pl1,  0.5); \
							c##MULTIPLY(*pu1, -0.5); \
							pl0 += 1; \
							pl1 += 1; \
							pu1 += n; \
						} \
					} \
			} \
		} while (0)

		SWITCH2(class[0], SPART);

#undef SPART

	}

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(4); /* x1, x0, to, from */
	return to;
}

/* skewpart(<denseMatrix>, trans) */
SEXP R_dense_skewpart(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ct = '\0';
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	return dense_skewpart(s_from, class, ct);
}

int dense_is_symmetric(SEXP obj, const char *class,
                       char ct, int exact, int checkDN)
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
		if (exact && (class[0] != 'z' || ct != 'C' || di != 'N'))
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

	char ul = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

	if (class[1] == 'g' && (class[0] != 'z' || ct != 'C')) {

#define ISS(c) \
	do { \
		c##TYPE *px = c##PTR(x), *pu = px, *pl = px; \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) { \
				if (c##NOT_IDEN(*pu, *pl)) \
					return 0; \
				pu += 1; \
				pl += n; \
			} \
			pu += n - j; \
			pl = px + j + 1; \
		} \
	} while (0)

	SWITCH5(class[0], ISS);

#undef ISS

	} else {

	Rcomplex *px = COMPLEX(x), *pu = px, *pl = px;
	int packed = class[2] == 'p';
	if (class[1] == 'g') {
		for (j = 0; j < n; ++j) {
			for (i = 0; i < j; ++i) {
				if (zNOT_CONJ(*pu, *pl))
					return 0;
				pu += 1;
				pl += n;
			}
			if (zNOT_REAL(*pu))
				return 0;
			pu += n - j;
			pl = px + j + 1;
		}
	} else if (class[1] == 's') {
		if (ul == 'U')
		for (j = 0; j < n; ++j) {
			for (i = 0; i < j; ++i) {
				if (zNOT_REAL(*pu))
					return 0;
				pu += 1;
			}
			if (ct == 'C' && zNOT_REAL(*pu))
				return 0;
			pu += 1;
			if (!packed)
			pu += n - j - 1;
		}
		else
		for (j = 0; j < n; ++j) {
			if (!packed)
			pl += j;
			if (ct == 'C' && zNOT_REAL(*pl))
				return 0;
			pl += 1;
			for (i = j + 1; i < n; ++i) {
				if (zNOT_REAL(*pl))
					return 0;
				pl += 1;
			}
		}
	} else {
		if (ul == 'U')
		for (j = 0; j < n; ++j) {
			for (i = 0; i < j; ++i) {
				if (zNOT_ZERO(*pu))
					return 0;
				pu += 1;
			}
			if (zNOT_REAL(*pu))
				return 0;
			pu += 1;
			if (!packed)
			pu += n - j - 1;
		}
		else
		for (j = 0; j < n; ++j) {
			if (!packed)
			pl += j;
			if (zNOT_REAL(*pl))
				return 0;
			pl += 1;
			for (i = j + 1; i < n; ++i) {
				if (zNOT_ZERO(*pl))
					return 0;
				pl += 1;
			}
		}
	}

	}

	return 1;
}

/* isSymmetric(<denseMatrix>, tol, tol1, trans, checkDN) */
SEXP R_dense_is_symmetric(SEXP s_obj,
                          SEXP s_trans, SEXP s_exact, SEXP s_checkDN)
{
	if (TYPEOF(s_obj) != OBJSXP) {
		/* defined in ./coerce.c : */
		SEXP matrix_as_dense(SEXP, const char *, char, char, char, int, int);
		s_obj = matrix_as_dense(s_obj, ".ge", '\0', '\0', '\0', 1, 0);
	}
	PROTECT(s_obj);
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	char ct;
	if (TYPEOF(s_trans) != STRSXP || LENGTH(s_trans) < 1 ||
	    (s_trans = STRING_ELT(s_trans, 0)) == NA_STRING ||
	    ((ct = CHAR(s_trans)[0]) != 'C' && ct != 'T'))
		error(_("invalid '%s' to '%s'"), "trans", __func__);

	int exact;
	if (TYPEOF(s_exact) != LGLSXP || LENGTH(s_exact) < 1 ||
	    (exact = LOGICAL(s_exact)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "exact", "TRUE", "FALSE");

	int checkDN;
	if (TYPEOF(s_checkDN) != LGLSXP || LENGTH(s_checkDN) < 1 ||
	    (checkDN = LOGICAL(s_checkDN)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "checkDN", "TRUE", "FALSE");

	int ans_ = dense_is_symmetric(s_obj, class, ct, exact, checkDN);
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
			return (upper != 0) ? 1 : -1;
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

#define IST(c) \
	do { \
		c##TYPE *px = c##PTR(x), *pu = px, *pl = px; \
		if (upper == NA_LOGICAL) { \
			for (j = 0; j < n; ++j) { \
				pl += j + 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*pl)) {	\
						j = n; \
						break; \
					} \
					pl += 1; \
				} \
			} \
			if (j == n) \
				return  1; \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (c##NOT_ZERO(*pu)) { \
						j = n; \
						break; \
					} \
					pu += 1; \
				} \
				pu += n - j; \
			} \
			if (j == n) \
				return -1; \
			return 0; \
		} else if (upper != 0) { \
			for (j = 0; j < n; ++j) { \
				pl += j + 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*pl)) \
						return 0; \
					pl += 1; \
				} \
			} \
			return  1; \
		} else { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (c##NOT_ZERO(*pu)) \
						return 0; \
					pu += 1; \
				} \
				pu += n - j; \
			} \
			return -1; \
		} \
	} while (0)

	SWITCH4(class[0], IST);

#undef IST

	return (upper != 0) ? 1 : -1;
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
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	if (TYPEOF(s_upper) != LGLSXP || LENGTH(s_upper) < 1)
		error(_("'%s' must be %s or %s or %s"), "upper", "TRUE", "FALSE", "NA");
	int upper = LOGICAL(s_upper)[0];

	int ans_ = dense_is_triangular(s_obj, class, upper);
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

	char ul = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p';

#define ISD(c) \
	do { \
		c##TYPE *px = c##PTR(x); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
				px += 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
				px += 1; \
				if (!packed) \
				px += n - j - 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
				px += j; \
				px += 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*px)) \
						return 0; \
					px += 1; \
				} \
			} \
		} \
		return 1; \
	} while (0)

	SWITCH4(class[0], ISD);

#undef ISD

	return 1;
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
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int ans_ = dense_is_diagonal(s_obj, class);
	SEXP ans = ScalarLogical(ans_ != 0);
	UNPROTECT(1);
	return ans;
}

#define nCAST(x) (x != 0)
#define lCAST(x) (x != 0)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

static
void dense_colsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char di, int narm, int mean,
                  SEXP ans)
{
	int i, j, count = -1, packed = class[2] == 'p';

#define SUM(c, d) \
	do { \
		c##TYPE *px = c##PTR(x); \
		d##TYPE *pa = d##PTR(ans), tmp; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				*pa = d##ZERO; \
				SUM_KERNEL(c, d, for (i = 0; i < m; ++i)); \
				pa += 1; \
			} \
		} else if (di == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				*pa = d##ZERO; \
				SUM_KERNEL(c, d, for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
				pa += 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				*pa = d##ZERO; \
				if (!packed) \
					px += j; \
				SUM_KERNEL(c, d, for (i = j; i < n; ++i)); \
				pa += 1; \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				*pa = d##UNIT; \
				SUM_KERNEL(c, d, for (i = 0; i < j; ++i)); \
				px += 1; \
				if (!packed) \
					px += n - j - 1; \
				pa += 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				*pa = d##UNIT; \
				if (!packed) \
					px += j; \
				px += 1; \
				SUM_KERNEL(c, d, for (i = j + 1; i < n; ++i)); \
				pa += 1; \
			} \
		} \
	} while (0)

#define SUM_KERNEL(c, d, __for__) \
	do { \
		if (mean) \
			count = m; \
		__for__ { \
			if (c##NOT_NA(*px)) { \
				tmp = c##CAST(*px); \
				d##INCREMENT_IDEN(*pa, tmp); \
			} \
			else if (!narm)	\
				*pa = d##NA; \
			else if (mean) \
				--count; \
			px += 1; \
		} \
		if (mean) \
			d##DIVIDE(*pa, count); \
	} while (0)

	switch (class[0]) {
	case 'n': if (mean) SUM(n, d); else SUM(n, i); break;
	case 'l': if (mean) SUM(l, d); else SUM(l, i); break;
	case 'i':                           SUM(i, d); break;
	case 'd':                           SUM(d, d); break;
	case 'z':                           SUM(z, z); break;
	default:                                       break;
	}

#undef SUM
#undef SUM_KERNEL

	return;
}

static
void dense_rowsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char di, int narm, int mean,
                  SEXP ans)
{
	int i, j, *count = NULL, packed = XLENGTH(x) != (int_fast64_t) m * n,
		sy = class[1] == 's', he = sy && ct == 'C';

	if (mean && narm) {
		Matrix_Calloc(count, m, int);
		for (i = 0; i < m; ++i)
			count[i] = n;
	}

#define SUM(c, d) \
	do { \
		c##TYPE *px = c##PTR(x); \
		d##TYPE *pa = d##PTR(ans), \
			tmp = (class[1] != 't' || di == 'N') ? c##ZERO : c##UNIT; \
		for (i = 0; i < m; ++i) \
			pa[i] = tmp; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(c, d, for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || di == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(c, d, for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				SUM_KERNEL(c, d, for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(c, d, for (i = 0; i < j; ++i)); \
				px += 1; \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				px += 1; \
				SUM_KERNEL(c, d, for (i = j + 1; i < n; ++i)); \
			} \
		} \
		if (mean) { \
			if (narm) \
				for (i = 0; i < m; ++i) \
					d##DIVIDE(pa[i], count[i]); \
			else \
				for (i = 0; i < m; ++i) \
					d##DIVIDE(pa[i], n); \
		} \
	} while (0)

#define SUM_KERNEL(c, d, __for__) \
	do { \
		__for__ { \
			if (c##NOT_NA(*px)) { \
				tmp = c##CAST(*px); \
				d##INCREMENT_IDEN(pa[i], tmp); \
				if (sy && i != j) { \
				if (he) \
				d##INCREMENT_CONJ(pa[j], tmp); \
				else \
				d##INCREMENT_IDEN(pa[j], tmp); \
				} \
			} else if (!narm) { \
				pa[i] = d##NA; \
				if (sy && i != j) \
				pa[j] = d##NA; \
			} else if (mean) { \
				--count[i]; \
				if (sy && i != j) \
				--count[j]; \
			} \
			px += 1; \
		} \
	} while (0)

	switch (class[0]) {
	case 'n': if (mean) SUM(n, d); else SUM(n, i); break;
	case 'l': if (mean) SUM(l, d); else SUM(l, i); break;
	case 'i':                           SUM(i, d); break;
	case 'd':                           SUM(d, d); break;
	case 'z':                           SUM(z, z); break;
	default:                                       break;
	}

#undef SUM
#undef SUM_KERNEL

	if (mean && narm)
		Matrix_Free(count, m);
	return;
}

SEXP dense_marginsum(SEXP obj, const char *class, int mg, int narm, int mean)
{
	narm = narm && class[0] != 'n';

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1], r = (mg == 0) ? m : n;

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

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

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
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
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int mg;
	if (TYPEOF(s_margin) != INTSXP || LENGTH(s_margin) < 1 ||
	    ((mg = INTEGER(s_margin)[0] - 1) != 0 && mg != 1))
		error(_("'%s' must be %d or %d"), "margin", 1, 2);

	int narm;
	if (TYPEOF(s_narm) != LGLSXP || LENGTH(s_narm) < 1 ||
	    (narm = LOGICAL(s_narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	int mean;
	if (TYPEOF(s_mean) != LGLSXP || LENGTH(s_mean) < 1 ||
	    (mean = LOGICAL(s_mean)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "mean", "TRUE", "FALSE");

	return dense_marginsum(s_obj, class, mg, narm, mean);
}

SEXP dense_sum(SEXP obj, const char *class, int narm)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
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

	SEXP ans = R_NilValue, x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p',
		sy = class[1] == 's', he = sy && ct == 'C';

#define SUM \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || di == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				SUM_KERNEL(for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(for (i = 0; i < j; ++i)); \
				px += 1; \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				px += 1; \
				SUM_KERNEL(for (i = j + 1; i < n; ++i)); \
			} \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	{
		int *px = LOGICAL(x);
		int_fast64_t s = (class[1] != 't' || di == 'N') ? 0 : n;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != 0) \
					s += (sy && i != j) ? 2 : 1; \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		if (s <= INT_MAX)
			ans = ScalarInteger((int) s);
		else
			ans = ScalarReal((double) s);
		break;
	}
	case 'l':
	case 'i':
	{
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		int_fast64_t s = (class[0] != 't' || di == 'N') ? 0 : n, t = 0;
		unsigned int count = 0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) { \
					unsigned int d = (sy && i != j) ? 2 : 1; \
					if (count > UINT_MAX - d) \
						TRY_INCREMENT(s, t); \
					t += (sy && i != j) ? 2LL * *px : *px; \
					count += d; \
				} else if (!narm) \
					return ScalarInteger(NA_INTEGER); \
				px += 1; \
			} \
		} while (0)

#define TRY_INCREMENT(x, y) \
		do { \
			if ((x >= 0) \
				? ( y <= INT_FAST64_MAX - x) \
				: (-y <= x - INT_FAST64_MIN)) { \
				x += y; \
				y = 0; \
				count = 0; \
			} else \
				goto over; \
		} while (0)

		SUM;

#undef SUM_KERNEL

		TRY_INCREMENT(s, t);

		if (s > INT_MIN && s <= INT_MAX)
			ans = ScalarInteger((int) s);
		else
			ans = ScalarReal((double) s);
		break;

over:
		px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		long double lr = (class[1] != 't' || di == 'N') ? 0.0 : (double) n;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) \
					lr += (sy && i != j) ? 2.0L * *px : *px; \
				else if (!narm) \
					return ScalarInteger(NA_INTEGER); \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		ans = ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
		break;
	}
	case 'd':
	{
		double *px = REAL(x);
		long double lr = (class[1] != 't' || di == 'N') ? 0.0 : (double) n;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && ISNAN(*px))) \
					lr += (sy && i != j) ? 2.0L * *px : *px; \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		ans = ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
		break;
	}
	case 'z':
	{
		Rcomplex *px = COMPLEX(x);
		long double lr = (class[1] != 't' || di == 'N') ? 0.0 : (double) n;
		long double li = 0.0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
					lr += (sy && i != j) ?                2.0L * (*px).r  : (*px).r; \
					li += (sy && i != j) ? ((he) ? 0.0L : 2.0L * (*px).i) : (*px).i; \
				} \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		Rcomplex tmp;
		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = ScalarComplex(tmp);
		break;
	}
	default:
		break;
	}

#undef SUM

	return ans;
}

/* sum(<denseMatrix>, na.rm=) */
SEXP R_dense_sum(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	if (TYPEOF(s_narm) != LGLSXP || LENGTH(s_narm) < 1 ||
	    (narm = LOGICAL(s_narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return dense_sum(s_obj, class, narm);
}

SEXP dense_prod(SEXP obj, const char *class, int narm)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
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

	SEXP ans, x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p',
		sy = class[1] == 's', he = sy && ct == 'C';
	long double lr = 1.0, li = 0.0;

#define PROD \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				PROD_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				PROD_KERNEL(for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				PROD_KERNEL(for (i = j; i < n; ++i)); \
			} \
		} else if (di == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				PROD_KERNEL(for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				if (!packed) \
					px += j; \
				PROD_KERNEL(for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				PROD_KERNEL(for (i = 0; i < j; ++i)); \
				px += 1; \
				if (!packed) \
					px += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (j == 1) { lr *= 0.0L; li *= 0.0L; } \
				if (!packed) \
					px += j; \
				px += 1; \
				PROD_KERNEL(for (i = j + 1; i < n; ++i)); \
			} \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	{
		int *px = LOGICAL(x);
		if (class[1] == 't') {
			if (n > 1 || (n == 1 && *px == 0))
				return ScalarReal(0.0);
			break;
		}

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px == 0) \
					return ScalarReal(0.0); \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	case 'l':
	case 'i':
	{
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) \
					lr *= (sy && i != j) ? (long double) *px * *px : *px; \
				else if (!narm) \
					lr *= NA_REAL; \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	case 'd':
	{
		double *px = REAL(x);

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && ISNAN(*px))) \
					lr *= (sy && i != j) ? (long double) *px * *px : *px; \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	case 'z':
	{
		Rcomplex *px = COMPLEX(x);
		long double lr0, li0;

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && (ISNAN((*px).r) || ISNAN((*px).i)))) { \
					lr0 = lr; li0 = li; \
					lr = lr0 * (*px).r - li0 * (*px).i; \
					li = li0 * (*px).r + lr0 * (*px).i; \
					if (sy && i != j) { \
					lr0 = lr; li0 = li; \
					if (he) { \
					lr = lr0 * (*px).r + li0 * (*px).i; \
					li = li0 * (*px).r - lr0 * (*px).i; \
					} else { \
					lr = lr0 * (*px).r - li0 * (*px).i; \
					li = li0 * (*px).r + lr0 * (*px).i; \
					} \
					} \
				} \
				px += 1; \
			} \
		} while (0)

		PROD;

#undef PROD_KERNEL

		break;
	}
	default:
		break;
	}

#undef PROD

	if (class[0] != 'z')
		ans = ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
	else {
		Rcomplex tmp;
		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = ScalarComplex(tmp);
	}
	return ans;
}

/* prod(<denseMatrix>, na.rm=) */
SEXP R_dense_prod(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	if (TYPEOF(s_narm) != LGLSXP || LENGTH(s_narm) < 1 ||
	    (narm = LOGICAL(s_narm)[0]) == NA_LOGICAL)
		error(_("'%s' must be %s or %s"), "narm", "TRUE", "FALSE");

	return dense_prod(s_obj, class, narm);
}
