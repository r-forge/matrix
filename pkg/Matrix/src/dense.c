#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "dense.h"

SEXP dense_band(SEXP from, const char *class, int a, int b)
{
	int packed = class[2] == 'p';

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) to be triangularMatrix */
	if ((m == 0 || n == 0 || (a <= 1 - m && b >= n - 1)) &&
	    (m != n || n > 1 || class[1] == 't'))
		return from;

	int ge, sy, tr;
	tr = class[1] == 't' || (m == n && (a >= 0 || b <= 0));
	sy = !tr && class[1] == 's' && a == -b;
	ge = !tr && !sy;

	char ul0 = '\0', ct0 = '\0', nu0 = '\0';
	if (class[1] != 'g')
		ul0 = UPLO(from);
	if (class[1] == 's' && class[0] == 'z')
		ct0 = TRANS(from);
	if (class[1] == 't') {
		/* Be fast if band contains entire triangle */
		if ((ul0 == 'U') ? (a <= 0 && b >= n - 1) : (a <= 1 - m && b >= 0))
			return from;
		nu0 = DIAG(from);
	}

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (ge) ? 'g' :                   ((sy) ? 's' : 't') ;
	cl[2] = (ge) ? 'e' : ((packed) ? 'p' : ((sy) ? 'y' : 'r'));
	SEXP to = PROTECT(newObject(cl));

	SET_DIM(to, m, n);
	SET_DIMNAMES(to, -(class[1] == 's' && !sy), DIMNAMES(from, 0));

	char ul1 = (tr && class[1] != 't') ? ((a >= 0) ? 'U' : 'L') : ul0;
	if (ul1 != '\0' && ul1 != 'U')
		SET_UPLO(to);
	if (ct0 != '\0' && ct0 != 'C' && sy)
		SET_TRANS(to);
	if (nu0 != '\0' && nu0 != 'N' && tr && a <= 0 && b >= 0)
		SET_DIAG(to);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), (!packed || ge) ? (R_xlen_t) m * n : (R_xlen_t) PACKED_LENGTH((size_t) n)));

	size_t
		m_ = (size_t) m,
		n_ = (size_t) n,
		a_ = (size_t) ((int_fast64_t) m + a),
		b_ = (size_t) ((int_fast64_t) m + b);

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (ge && class[1] != 'g') { \
		if (!packed) \
			c##NAME(force2)(px1, px0, n_, ul0, ct0, nu0); \
		else \
			c##NAME( pack1)(px1, px0, n_, ul0, ct0, nu0); \
		px0 = NULL; \
		packed = 0; \
		} else if (tr && class[1] == 's' && ul0 != ul1) { \
		if (!packed) \
			c##NAME(trans2)(px1, px0, m_, n_ , ct0); \
		else \
			c##NAME(trans1)(px1, px0, n_, ul0, ct0); \
		px0 = NULL; \
		ul0 = ul1; \
		} \
		if (!packed) \
			c##NAME( band2)(px1, px0, m_, n_ , a_, b_); \
		else \
			c##NAME( band1)(px1, px0, n_, ul0, a_, b_); \
	} while (0)

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

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

	int *pdim = DIM(s_from), m = pdim[0], n = pdim[1], a, b;
	if (s_a == R_NilValue)
		a = -m;
	else if ((a = Rf_asInteger(s_a)) == NA_INTEGER || a < -m || a > n)
		Rf_error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		         "k1", a, "-Dim[1]", -m, "Dim[2]", n);
	if (s_b == R_NilValue)
		b = n;
	else if ((b = Rf_asInteger(s_b)) == NA_INTEGER || b < -m || b > n)
		Rf_error(_("'%s' (%d) must be an integer from %s (%d) to %s (%d)"),
		         "k2", b, "-Dim[1]", -m, "Dim[2]", n);
	else if (b < a)
		Rf_error(_("'%s' (%d) must be less than or equal to '%s' (%d)"),
		         "k1", a, "k2", b);

	s_from = dense_band(s_from, class, a, b);
	UNPROTECT(1); /* s_from */
	return s_from;
}

SEXP dense_diag_get(SEXP obj, const char *class, int names)
{
	int packed = class[2] == 'p';

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP ans;

	if (nu != '\0' && nu != 'N') {

		PROTECT(ans = allocUnit  (kindToType(class[0]), r));

	} else {

		PROTECT(ans = Rf_allocVector(kindToType(class[0]), r));

		SEXP x = GET_SLOT(obj, Matrix_xSym);
		size_t m_ = (size_t) m, n_ = (size_t) n, r_ = (size_t) r;

#define TEMPLATE(c) \
		do { \
			c##TYPE *pa = c##PTR(ans), *px = c##PTR(x); \
			if (!packed) \
				c##NAME(copy2)(r_, pa, 1, px, m_ + 1); \
			else if (ul == 'U') \
				c##NAME(copy1)(n_, pa, 1, 0, 0, px, 2 , 1, 0); \
			else \
				c##NAME(copy1)(n_, pa, 1, 0, 0, px, n_, 1, 1); \
		} while (0)

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		if (class[0] == 'n')
			naToUnit(ans);
		if (class[0] == 'z' && class[1] == 's' && ct == 'C')
			zvreal(COMPLEX(ans), NULL, r_);

	}

	if (names) {
		/* NB: The logic here must be adjusted once the validity method
		   for 'symmetricMatrix' enforces symmetric 'Dimnames'
		*/
		SEXP dn = PROTECT(DIMNAMES(obj, 0)),
			rn = VECTOR_ELT(dn, 0),
			cn = VECTOR_ELT(dn, 1);
		if (cn == R_NilValue) {
			if (class[1] == 's' && rn != R_NilValue)
				Rf_setAttrib(ans, R_NamesSymbol, rn);
		} else {
			if (class[1] == 's')
				Rf_setAttrib(ans, R_NamesSymbol, cn);
			else if (rn != R_NilValue &&
			         (rn == cn || equalString(rn, cn, r)))
				Rf_setAttrib(ans, R_NamesSymbol, (r == m) ? rn : cn);
		}
		UNPROTECT(1); /* dn */
	}

	UNPROTECT(1); /* ans */
	return ans;
}

/* diag(<denseMatrix>, names=) */
SEXP R_dense_diag_get(SEXP s_obj, SEXP s_names)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int names;
	VALID_LOGIC2(s_names, names);

	return dense_diag_get(s_obj, class, names);
}

SEXP dense_diag_set(SEXP from, const char *class, SEXP value, int new)
{
	SEXP to = PROTECT(newObject(class));
	int packed = class[2] == 'p';

	int *pdim = DIM(from), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	char ul = '\0', ct = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);

	SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (new) {
		x = duplicateVector(x);
		UNPROTECT(1); /* x */
		PROTECT(x);
	}

	size_t m_ = (size_t) m, n_ = (size_t) n, r_ = (size_t) r,
		d_ = (LENGTH(value) == r) ? 1 : 0;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x), *pv = c##PTR(value); \
		if (!packed) \
			c##NAME(copy2)(r_, px, m_ + 1, pv, d_); \
		else if (ul == 'U') \
			c##NAME(copy1)(n_, px, 2 , 1, 0, pv, d_, 0, 0); \
		else \
			c##NAME(copy1)(n_, px, n_, 1, 1, pv, d_, 0, 0); \
	} while (0)

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x);

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
		Rf_error(_("replacement diagonal has incompatible type \"%s\""),
		         Rf_type2char(tv));
		break;
	}

	int *pdim = DIM(s_from), m = pdim[0], n = pdim[1],
		r = (m < n) ? m : n;
	if (XLENGTH(s_value) != 1 && XLENGTH(s_value) != r)
		Rf_error(_("replacement diagonal has wrong length"));

	int new = 1;
	if (tv <= tx) {
		/* defined in ./coerce.c : */
		SEXP dense_as_general(SEXP, const char *, int);
		if (class[1] == 's' && class[0] == 'z' && tv == tx &&
		    TRANS(s_from) == 'C') {
			PROTECT(s_from = dense_as_general(s_from, class, 1));
			class = Matrix_class(s_from, valid_dense, 6, __func__);
			UNPROTECT(1); /* s_from */
			new = 0;
		}
		PROTECT(s_from);
		PROTECT(s_value = Rf_coerceVector(s_value, tx));
	} else {
		/* defined in ./coerce.c : */
		SEXP dense_as_kind(SEXP, const char *, char, int);
#ifndef MATRIX_ENABLE_IMATRIX
		if (tv == INTSXP) {
		PROTECT(s_from = dense_as_kind(s_from, class, 'd', 0));
		PROTECT(s_value = Rf_coerceVector(s_value, REALSXP));
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
	UNPROTECT(2); /* s_value, s_from */
	return s_from;
}

SEXP dense_transpose(SEXP from, const char *class, char op_ct)
{
	if (class[0] != 'z')
		op_ct = '\0';

	SEXP to = PROTECT(newObject(class));
	int packed = class[2] == 'p';

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, n, m);
	SET_DIMNAMES(to, class[1] != 's' && class[1] != 'p' && class[1] != 'o', DIMNAMES(from, 0));

	char ul = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) == 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] == 'o')
		COPY_SLOT(to, from, Matrix_sdSym);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));
	size_t m_ = (size_t) m, n_ = (size_t) n;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed) \
			c##NAME(trans2)(px1, px0, m_, n_, op_ct); \
		else \
			c##NAME(trans1)(px1, px0, n_, ul, op_ct); \
	} while (0)

	SWITCH4((class[0] == 'c') ? 'd' : class[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* t(<denseMatrix>) */
SEXP R_dense_transpose(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	return dense_transpose(s_from, class, ct);
}

SEXP dense_force_symmetric(SEXP from, const char *class, char op_ul, char op_ct)
{
	if (class[0] != 'z')
		op_ct = '\0';

	char
		ul0 = '\0', ul1 = 'U',
		ct0 = '\0', ct1 = (class[0] == 'z') ? 'C' : '\0',
		nu0 = '\0';
	if (class[1] != 'g')
		ul0 = ul1 = UPLO(from);
	if (class[1] == 's' && class[0] == 'z')
		ct0 = ct1 = TRANS(from);
	if (class[1] == 't')
		nu0 = DIAG(from);
	if (op_ul != '\0')
		ul1 = op_ul;
	if (op_ct != '\0')
		ct1 = op_ct;

	if (class[1] == 's' && ul0 == ul1 && ct0 == ct1)
		return from;

	int packed = class[2] == 'p';

	char cl[] = ".s.Matrix";
	cl[0] = class[0];
	cl[2] = (packed) ? 'p' : 'y';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("attempt to symmetrize a non-square matrix"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));
	if (ul1 != 'U')
		SET_UPLO(to);
	if (ct1 != 'C' && ct1 != '\0')
		SET_TRANS(to);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));

	if ((class[1] == 'g') ||
	    (class[1] == 's' && ul0 == ul1 && ct0 != 'C') ||
	    (class[1] == 't' && ul0 == ul1 && nu0 == 'N'))
		SET_SLOT(to, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));
		size_t n_ = (size_t) n;

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			if (!packed) \
				c##NAME(force2)(px1,  px0, n_, ul0, ct0, nu0); \
			else if (ul0 == ul1) \
				c##NAME(force1)(px1,  px0, n_, ul0, ct0, nu0); \
			else { \
				c##NAME(trans1)(px1,  px0, n_, ul0, ct0); \
				c##NAME(force1)(px1, NULL, n_, ul1, ct0, nu0); \
			} \
		} while (0)

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		SET_SLOT(to, Matrix_xSym, x1);
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
	if (s_uplo != R_NilValue)
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return dense_force_symmetric(s_from, class, ul, ct);
}

SEXP dense_symmpart(SEXP from, const char *class, char op_ul, char op_ct)
{
	/* defined in ./coerce.c : */
	SEXP dense_as_kind(SEXP, const char *, char, int);
	PROTECT(from = dense_as_kind(from, class, ',', 0));

	if (class[1] != 'g')
		op_ul = '\0';
	if (class[0] != 'z')
		op_ct = '\0';

	char ct = '\0';
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(from);

	if (class[1] == 's' && op_ct == ct) {
		UNPROTECT(1); /* from */
		return from;
	}

	int packed = class[2] == 'p';

	char cl[] = ".s.Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[2] = (packed) ? 'p' : 'y';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error((op_ct == 'C')
		         ? _("attempt to get Hermitian part of non-square matrix")
		         : _("attempt to get symmetric part of non-square matrix"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));

	char ul = '\0', nu = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U')
		SET_UPLO(to);
	if (class[1] == 't' && (nu = DIAG(from)) != 'N')
		;

	if (op_ul != '\0' && op_ul != 'U')
		SET_UPLO(to);
	if (op_ct != '\0' && op_ct != 'C')
		SET_TRANS(to);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));

	if (class[1] == 's') {

		/* Symmetric part of Hermitian matrix is real part */
		/* Hermitian part of symmetric matrix is real part */
		zvreal(COMPLEX(x1), COMPLEX(x0), (size_t) XLENGTH(x0));

	} else {

		int i, j;

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *pu0 = px0, *pl0 = px0; \
			c##TYPE *px1 = c##PTR(x1), *pu1 = px1, *pl1 = px1; \
			if (!packed) \
				memset(px1, 0, sizeof(c##TYPE) * (size_t) XLENGTH(x1)); \
			if (class[1] == 'g') { \
				if (op_ul == 'U') { \
					if (op_ct == 'C') \
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
						pu0 += 1; \
						pu1 += 1; \
						pu0 += n - j - 1; \
						pu1 += n - j - 1; \
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
						pu0 += 1; \
						pu1 += 1; \
						pu0 += n - j - 1; \
						pu1 += n - j - 1; \
						pl0 = px0 + j + 1; \
					} \
				} else { \
					if (op_ct == 'C') \
					for (j = 0; j < n; ++j) { \
						pl0 += j; \
						pl1 += j; \
						pu0 = pl0 + n; \
						c##ASSIGN_PROJ_REAL(*pl1, *pl0); \
						pl0 += 1; \
						pl1 += 1; \
						for (i = j + 1; i < n; ++i) { \
							c##ASSIGN_IDEN(*pl1, *pl0); \
							c##INCREMENT_CONJ(*pl1, *pu0); \
							c##MULTIPLY(*pl1, 0.5); \
							pl0 += 1; \
							pl1 += 1; \
							pu0 += n; \
						} \
					} \
					else \
					for (j = 0; j < n; ++j) { \
						pl0 += j; \
						pl1 += j; \
						pu0 = pl0 + n; \
						c##ASSIGN_IDEN(*pl1, *pl0); \
						pl0 += 1; \
						pl1 += 1; \
						for (i = j + 1; i < n; ++i) { \
							c##ASSIGN_IDEN(*pl1, *pl0); \
							c##INCREMENT_IDEN(*pl1, *pu0); \
							c##MULTIPLY(*pl1, 0.5); \
							pl0 += 1; \
							pl1 += 1; \
							pu0 += n; \
						} \
					} \
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
						if (nu != 'N') \
							c##SET_UNIT(*pu1); \
						else if (op_ct == 'C') \
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
						if (nu != 'N') \
							c##SET_UNIT(*pl1); \
						else if (op_ct == 'C') \
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

		SWITCH2(class[0], TEMPLATE);

#undef TEMPLATE

	}

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(4); /* x1, x0, to, from */
	return to;
}

/* symmpart(<denseMatrix>, uplo, trans) */
SEXP R_dense_symmpart(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ul, ct;
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return dense_symmpart(s_from, class, ul, ct);
}

SEXP dense_skewpart(SEXP from, const char *class, char op_ct)
{
	/* defined in ./coerce.c : */
	SEXP dense_as_kind(SEXP, const char *, char, int);
	PROTECT(from = dense_as_kind(from, class, ',', 0));

	if (class[0] != 'z')
		op_ct = '\0';

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 's') ? 's' : 'g';
	cl[2] = (class[1] == 's') ? class[2] : 'e';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error((op_ct == 'C')
		         ? _("attempt to get skew-Hermitian part of non-square matrix")
		         : _("attempt to get skew-symmetric part of non-square matrix"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U' && class[1] == 's')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && (nu = DIAG(from)) != 'N')
		;

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1;

	if (class[1] == 's') {

		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		/* Skew-Hermitian part of symmetric matrix is imaginary part */
		if (op_ct == ct)
		PROTECT(x1 = allocZero(TYPEOF(x0), XLENGTH(x0)));
		else {
		PROTECT(x1 = Rf_allocVector(CPLXSXP, XLENGTH(x0)));
		zvimag(COMPLEX(x1), COMPLEX(x0), (size_t) XLENGTH(x0));
		}

	} else {

		if ((int_fast64_t) n * n > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		PROTECT(x1 = Rf_allocVector(TYPEOF(x0), (R_xlen_t) n * n));

		int i, j;

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *pu0 = px0, *pl0 = px0; \
			c##TYPE *px1 = c##PTR(x1), *pu1 = px1, *pl1 = px1; \
			if (class[1] == 'g') { \
				if (op_ct == 'C') \
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
						if (nu != 'N' || op_ct != 'C') \
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
						if (nu != 'N' || op_ct != 'C') \
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

		SWITCH2(class[0], TEMPLATE);

#undef TEMPLATE

	}

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(4); /* x1, x0, to, from */
	return to;
}

/* skewpart(<denseMatrix>, trans) */
SEXP R_dense_skewpart(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	return dense_skewpart(s_from, class, ct);
}

int dense_is_symmetric(SEXP obj, const char *class,
                       char op_ct, int exact, int checkDN)
{
	if (class[0] != 'z') {
		op_ct = '\0';
		if (class[0] != 'd')
			exact = 1;
	}

	char ct = '\0';
	if (class[1] == 's') {
		if (class[0] != 'z')
			return 1;
		ct = TRANS(obj);
		if (op_ct == ct)
			return 1;
		checkDN = 0;
	}

	if (checkDN && !DimNames_is_symmetric(DIMNAMES(obj, 0)))
		return 0;

	char nu = '\0';
	if (class[1] == 't') {
		nu = DIAG(obj);
		if (exact && (nu != 'N' || op_ct != 'C'))
			return dense_is_diagonal(obj, class);
	}

	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		return 0;
	if (n == 0 || (n == 1 && op_ct != 'C'))
		return 1;
	if (!exact)
		return NA_LOGICAL; /* do inexact numerical test in R */

	char ul = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

	if (class[1] == 'g' && op_ct != 'C') {

#define TEMPLATE(c) \
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

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

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
			if (zNOT_ZERO_IMAG(*pu))
				return 0;
			pu += n - j;
			pl = px + j + 1;
		}
	} else if (class[1] == 's') {
		/* Testing if Hermitian matrix is symmetric */
		/*      or if symmetric matrix is Hermitian */
		/* <=====> if matrix is real                */
		if (ul == 'U')
		for (j = 0; j < n; ++j) {
			for (i = 0; i < j; ++i) {
				if (zNOT_ZERO_IMAG(*pu))
					return 0;
				pu += 1;
			}
			if (op_ct == 'C' && zNOT_ZERO_IMAG(*pu))
				return 0;
			pu += 1;
			if (!packed)
			pu += n - j - 1;
		}
		else
		for (j = 0; j < n; ++j) {
			if (!packed)
			pl += j;
			if (op_ct == 'C' && zNOT_ZERO_IMAG(*pl))
				return 0;
			pl += 1;
			for (i = j + 1; i < n; ++i) {
				if (zNOT_ZERO_IMAG(*pl))
					return 0;
				pl += 1;
			}
		}
	} else {
		/* Testing if non-unit triangular matrix is Hermitian */
		/* <=====> if matrix is real and diagonal             */
		if (ul == 'U')
		for (j = 0; j < n; ++j) {
			for (i = 0; i < j; ++i) {
				if (zNOT_ZERO(*pu))
					return 0;
				pu += 1;
			}
			if (zNOT_ZERO_IMAG(*pu))
				return 0;
			pu += 1;
			if (!packed)
			pu += n - j - 1;
		}
		else
		for (j = 0; j < n; ++j) {
			if (!packed)
			pl += j;
			if (zNOT_ZERO_IMAG(*pl))
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
	VALID_TRANS(s_trans, ct);

	int exact, checkDN;
	VALID_LOGIC2(s_exact  , exact  );
	VALID_LOGIC2(s_checkDN, checkDN);

	int ans_ = dense_is_symmetric(s_obj, class, ct, exact, checkDN);
	SEXP ans = Rf_ScalarLogical(ans_);
	UNPROTECT(1); /* s_obj */
	return ans;
}

int dense_is_triangular(SEXP obj, const char *class, char op_ul)
{
	if (class[1] == 't') {
		char ul = UPLO(obj);
		if (op_ul == '\0' || op_ul == ul)
			return (   ul == 'U') ? 1 : -1;
		else if (dense_is_diagonal(obj, class))
			return (op_ul == 'U') ? 1 : -1;
		else
			return 0;
	}

	if (class[1] == 's') {
		if (!dense_is_diagonal(obj, class))
			return 0;
		else if (op_ul != '\0')
			return (op_ul == 'U') ? 1 : -1;
		else {
			char ul = UPLO(obj);
			return (   ul == 'U') ? 1 : -1;
		}
	}

	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		return 0;
	if (n <= 1)
		return (op_ul == '\0' || op_ul == 'U') ? 1 : -1;

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x), *pu = px, *pl = px; \
		if (op_ul == '\0') { \
			for (j = 0; j < n; ++j) { \
				pl += j + 1; \
				for (i = j + 1; i < n; ++i) { \
					if (c##NOT_ZERO(*pl)) { \
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
		} else if (op_ul == 'U') { \
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

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

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
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int up;
	VALID_LOGIC3(s_upper, up);

	int ans_ = dense_is_triangular(s_obj, class,
		(up == NA_LOGICAL) ? '\0' : ((up != 0) ? 'U' : 'L'));
	SEXP ans = Rf_allocVector(LGLSXP, 1);
	LOGICAL(ans)[0] = ans_ != 0;
	if (up == NA_LOGICAL && ans_ != 0) {
		PROTECT(ans);
		static
		SEXP kindSym = NULL;
		SEXP kindVal = PROTECT(Rf_mkString((ans_ > 0) ? "U" : "L"));
		if (!kindSym)
			kindSym = Rf_install("kind");
		Rf_setAttrib(ans, kindSym, kindVal);
		UNPROTECT(2); /* kindVal, ans */
	}
	UNPROTECT(1); /* s_obj */
	return ans;
}

int dense_is_diagonal(SEXP obj, const char *class)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		return 0;
	if (n <= 1)
		return 1;

	char ul = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p';

#define TEMPLATE(c) \
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

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

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
	SEXP ans = Rf_ScalarLogical(ans_ != 0);
	UNPROTECT(1); /* s_obj */
	return ans;
}

#define nCAST(x) (x != 0)
#define lCAST(x) (x)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#define SUM_TYPEOF(c) (c == 'z') ? CPLXSXP : ((mean || c == 'd' || c == 'i') ? REALSXP : INTSXP)

static
void dense_colsum(SEXP x, const char *class,
                  int m, int n, char ul, char ct, char nu, int narm, int mean,
                  SEXP ans)
{
	int i, j, count = -1, packed = class[2] == 'p';

#define SUM(c0, c1) \
	do { \
		c0##TYPE *px0 = c0##PTR(  x); \
		c1##TYPE *px1 = c1##PTR(ans); \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##ZERO; \
				SUM_KERNEL(c0, c1, for (i = 0; i < m; ++i)); \
				px1 += 1; \
			} \
		} else if (nu == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##ZERO; \
				SUM_KERNEL(c0, c1, for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px0 += n - j - 1; \
				px1 += 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##ZERO; \
				if (!packed) \
					px0 += j; \
				SUM_KERNEL(c0, c1, for (i = j; i < n; ++i)); \
				px1 += 1; \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##UNIT; \
				SUM_KERNEL(c0, c1, for (i = 0; i < j; ++i)); \
				px0 += 1; \
				if (!packed) \
					px0 += n - j - 1; \
				px1 += 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				*px1 = c1##UNIT; \
				if (!packed) \
					px0 += j; \
				px0 += 1; \
				SUM_KERNEL(c0, c1, for (i = j + 1; i < n; ++i)); \
				px1 += 1; \
			} \
		} \
	} while (0)

#define SUM_KERNEL(c0, c1, __for__) \
	do { \
		if (mean) \
			count = m; \
		__for__ { \
			if (c0##NOT_NA(*px0)) \
				c1##INCREMENT_IDEN(*px1, c0##CAST(*px0)); \
			else if (!narm) \
				*px1 = c1##NA; \
			else if (mean) \
				--count; \
			px0 += 1; \
		} \
		if (mean) \
			c1##DIVIDE(*px1, count); \
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
                  int m, int n, char ul, char ct, char nu, int narm, int mean,
                  SEXP ans)
{
	int i, j, *count = NULL, packed = XLENGTH(x) != (int_fast64_t) m * n,
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';

	if (mean && narm) {
		Matrix_Calloc(count, m, int);
		for (i = 0; i < m; ++i)
			count[i] = n;
	}

#define SUM(c0, c1) \
	do { \
		c0##TYPE *px0 = c0##PTR(  x), tmp0; \
		c1##TYPE *px1 = c1##PTR(ans), tmp1 = (un) ? c0##UNIT : c0##ZERO; \
		for (i = 0; i < m; ++i) \
			px1[i] = tmp1; \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(c0, c1, for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || nu == 'N') { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(c0, c1, for (i = 0; i <= j; ++i)); \
				if (!packed) \
					px0 += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px0 += j; \
				SUM_KERNEL(c0, c1, for (i = j; i < n; ++i)); \
			} \
		} else { \
			if (ul == 'U') \
			for (j = 0; j < n; ++j) { \
				SUM_KERNEL(c0, c1, for (i = 0; i < j; ++i)); \
				px0 += 1; \
				if (!packed) \
					px0 += n - j - 1; \
			} \
			else \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px0 += j; \
				px0 += 1; \
				SUM_KERNEL(c0, c1, for (i = j + 1; i < n; ++i)); \
			} \
		} \
		if (mean) { \
			if (!narm) \
				for (i = 0; i < m; ++i) \
					c1##DIVIDE(px1[i], n); \
			else \
				for (i = 0; i < m; ++i) \
					c1##DIVIDE(px1[i], count[i]); \
		} \
	} while (0)

#define SUM_KERNEL(c0, c1, __for__) \
	do { \
		__for__ { \
			if (he && i == j) \
			c0##ASSIGN_PROJ_REAL(tmp0, *px0); \
			else \
			c0##ASSIGN_IDEN     (tmp0, *px0); \
			if (c0##NOT_NA(tmp0)) { \
				c1##INCREMENT_IDEN(px1[i], c0##CAST(tmp0)); \
				if (sy && i != j) { \
				if (he) \
				c1##INCREMENT_CONJ(px1[j], c0##CAST(tmp0)); \
				else \
				c1##INCREMENT_IDEN(px1[j], c0##CAST(tmp0)); \
				} \
			} else if (!narm) { \
				px1[i] = c1##NA; \
				if (sy && i != j) \
				px1[j] = c1##NA; \
			} else if (mean) { \
				--count[i]; \
				if (sy && i != j) \
				--count[j]; \
			} \
			px0 += 1; \
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

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		r = (mg == 0) ? m : n;

	SEXP ans = PROTECT(Rf_allocVector(SUM_TYPEOF(class[0]), r)),
		x = PROTECT(GET_SLOT(obj, Matrix_xSym));

	SEXP dimnames = DIMNAMES(obj, -(class[1] == 's')),
		marnames = VECTOR_ELT(dimnames, mg);
	if (marnames != R_NilValue) {
		PROTECT(marnames);
		Rf_setAttrib(ans, R_NamesSymbol, marnames);
		UNPROTECT(1); /* marnames */
	}

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	if (mg == 0 || class[1] == 's')
		dense_rowsum(x, class, m, n, ul, ct, nu, narm, mean, ans);
	else
		dense_colsum(x, class, m, n, ul, ct, nu, narm, mean, ans);

	UNPROTECT(2); /* x, ans */
	return ans;
}

/* (row|col)(Sums|Means)(<denseMatrix>, na.rm=) */
SEXP R_dense_marginsum(SEXP s_obj, SEXP s_margin, SEXP s_narm, SEXP s_mean)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int mg;
	VALID_MARGIN(s_margin, mg);

	int narm, mean;
	VALID_LOGIC2(s_narm, narm);
	VALID_LOGIC2(s_mean, mean);

	return dense_marginsum(s_obj, class, mg, narm, mean);
}

SEXP dense_sum(SEXP obj, const char *class, int narm)
{
	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym),
		ans = R_NilValue;
	int i, j, packed = class[2] == 'p',
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';

#define SUM \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) \
				SUM_KERNEL(for (i = 0; i < m; ++i)); \
		} else if (class[1] == 's' || nu == 'N') { \
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
		int_fast64_t s = (un) ? n : 0;

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
			ans = Rf_ScalarInteger((int) s);
		else
			ans = Rf_ScalarReal((double) s);
		break;
	}
	case 'l':
	case 'i':
	{
		int *px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		int_fast64_t s = (un) ? n : 0, t = 0;
		unsigned int count = 0;

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

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) { \
					unsigned int d = (sy && i != j) ? 2 : 1; \
					if (count > UINT_MAX - d) \
						TRY_INCREMENT(s, t); \
					t += (sy && i != j) ? 2LL * *px : *px; \
					count += d; \
				} \
				else if (!narm) \
					return Rf_ScalarInteger(NA_INTEGER); \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		TRY_INCREMENT(s, t);

		if (s > INT_MIN && s <= INT_MAX)
			ans = Rf_ScalarInteger((int) s);
		else
			ans = Rf_ScalarReal((double) s);
		break;
over:
		px = (class[0] == 'l') ? LOGICAL(x) : INTEGER(x);
		long double lr = (un) ? (long double) n : 0.0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px != NA_INTEGER) \
					lr += (sy && i != j) ? 2.0L * *px : *px; \
				else if (!narm) \
					return Rf_ScalarInteger(NA_INTEGER); \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
		break;
	}
	case 'd':
	{
		double *px = REAL(x);
		long double lr = (un) ? (long double) n : 0.0;

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

		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
		break;
	}
	case 'z':
	{
		Rcomplex *px = COMPLEX(x), tmp;
		long double lr = (un) ? (long double) n : 0.0;
		long double li = 0.0;

#define SUM_KERNEL(__for__) \
		do { \
			__for__ { \
				if (!(narm && (ISNAN((*px).r) || (!he && ISNAN((*px).i))))) { \
					lr += (sy && i != j) ? 2.0L * (*px).r : (*px).r; \
					if (!he) \
					li += (sy && i != j) ? 2.0L * (*px).i : (*px).i; \
				} \
				px += 1; \
			} \
		} while (0)

		SUM;

#undef SUM_KERNEL

		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = Rf_ScalarComplex(tmp);
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
	VALID_LOGIC2(s_narm, narm);

	return dense_sum(s_obj, class, narm);
}

SEXP dense_prod(SEXP obj, const char *class, int narm)
{
	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = class[2] == 'p',
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N';
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
		} else if (nu == 'N') { \
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
			if (n > 1 || (n == 1 && !un && *px == 0))
				return Rf_ScalarReal(0.0);
			break;
		}

#define PROD_KERNEL(__for__) \
		do { \
			__for__ { \
				if (*px == 0) \
					return Rf_ScalarReal(0.0); \
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
				if (!(narm && (ISNAN((*px).r) || (!(he && i == j) && ISNAN((*px).i))))) { \
					if (he) { \
						lr0 = (*px).r; \
						if (i != j) { \
						li0 = (*px).i; \
						lr0 = lr0 * lr0 + li0 * li0; \
						} \
						lr *= lr0; \
						li *= lr0; \
					} else { \
						lr0 = lr; li0 = li; \
						lr = lr0 * (*px).r - li0 * (*px).i; \
						li = li0 * (*px).r + lr0 * (*px).i; \
						if (i != j) { \
						lr0 = lr; li0 = li; \
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

	SEXP ans;
	if (class[0] != 'z')
		ans = Rf_ScalarReal(LONGDOUBLE_AS_DOUBLE(lr));
	else {
		Rcomplex tmp;
		tmp.r = LONGDOUBLE_AS_DOUBLE(lr);
		tmp.i = LONGDOUBLE_AS_DOUBLE(li);
		ans = Rf_ScalarComplex(tmp);
	}
	return ans;
}

/* prod(<denseMatrix>, na.rm=) */
SEXP R_dense_prod(SEXP s_obj, SEXP s_narm)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);

	int narm;
	VALID_LOGIC2(s_narm, narm);

	return dense_prod(s_obj, class, narm);
}
