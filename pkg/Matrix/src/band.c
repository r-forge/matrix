/* C implementation of methods for band */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

SEXP dense_band(SEXP from, const char *class, int a, int b)
{
	int packed = class[2] == 'p';

	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) */
	/* to be triangularMatrix                       */
	/* Don't return 'from' if from@x has attributes */
	/* (notably 'class' or 'dim'), which can happen */
	/* if from = matrix_as_dense(..., new=0).       */
	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	if ((m == 0 || n == 0 || (a <= 1 - m && b >= n - 1)) &&
	    (m != n || n > 1 || class[1] == 't') &&
	    !ANY_ATTRIB(GET_SLOT(from, Matrix_xSym)))
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

SEXP sparse_band(SEXP from, const char *class, int a, int b)
{
	/* Need tri[ul](<0-by-0>) and tri[ul](<1-by-1>) */
	/* to be triangularMatrix                       */
	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
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
	cl[1] = (ge) ? 'g' : ((sy) ? 's' : 't');
	cl[2] = class[2];
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

	if (class[2] != 'T') {

		if (class[2] == 'R') {
			SWAP(m, n, int, );
			SWAP(a, b, int, -);
		}

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0)));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			*pp1 = INTEGER(p1), *iwork = NULL,
			j, k, kend, nnz0 = pp0[n], nnz1 = 0, d;
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] != 's' || sy) {
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++nnz1;
					++k;
				}
				pp1[j] = nnz1;
			}
			if (nnz0 == nnz1) {
				SET_SLOT(to, iSym, i0);
				if (class[0] != 'n') {
					SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
					SET_SLOT(to, Matrix_xSym, x0);
					UNPROTECT(1); /* x0 */
				}
				UNPROTECT(4); /* p1, i0, p0, to */
				return to;
			}
		} else {
			memset(pp1, 0, sizeof(int) * (size_t) n);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if ((d = j - pi0[k]) >= a && d <= b)
						++pp1[j];
					if (d != 0 && -d >= a && -d <= b)
						++pp1[pi0[k]];
					++k;
				}
			}
			Matrix_Calloc(iwork, n, int);
			for (j = 0; j < n; ++j) {
				iwork[j] = nnz1;
				nnz1 += pp1[j];
				pp1[j] = nnz1;
			}
		}

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define TEMPLATE(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz1)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] != 's' || sy) \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							*(pi1++) = pi0[k]; \
							c##IF_NPATTERN( \
							*(px1++) = px0[k]; \
							); \
						} \
						++k; \
					} \
				} \
			else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if ((d = j - pi0[k]) >= a && d <= b) { \
							c##IF_NPATTERN( \
							if (ct0 == 'C' && d == 0) \
								c##ASSIGN_PROJ_REAL(px1[iwork[j]], px0[k]); \
							else \
								c##ASSIGN_IDEN     (px1[iwork[j]], px0[k]); \
							); \
							pi1[iwork[j]++] = pi0[k]; \
						} \
						if (d != 0 && -d >= a && -d <= b) { \
							c##IF_NPATTERN( \
							if (ct0 == 'C') \
								c##ASSIGN_CONJ(px1[iwork[pi0[k]]], px0[k]); \
							else \
								c##ASSIGN_IDEN(px1[iwork[pi0[k]]], px0[k]); \
							); \
							pi1[iwork[pi0[k]]++] = j; \
						} \
						++k; \
					} \
				} \
				Matrix_Free(iwork, n); \
			} \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(4); /* i1, p1, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), d;
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1 = 0;

		if (class[1] != 's' || sy) {
			for (k = 0; k < nnz0; ++k)
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
			if (nnz0 == nnz1) {
				SET_SLOT(to, Matrix_iSym, i0);
				SET_SLOT(to, Matrix_jSym, j0);
				if (class[0] != 'n') {
					SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
					SET_SLOT(to, Matrix_xSym, x0);
					UNPROTECT(1); /* x0 */
				}
				UNPROTECT(3); /* j0, i0, to */
				return to;
			}
		} else
			for (k = 0; k < nnz0; ++k) {
				if ((d = pj0[k] - pi0[k]) >= a && d <= b)
					++nnz1;
				if (d != 0 && -d >= a && -d <= b)
					++nnz1;
			}

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
			j1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#define TEMPLATE(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz1)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			if (class[1] != 's' || sy) \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						c##IF_NPATTERN( \
						*(px1++) = px0[k]; \
						); \
					} \
				} \
			else \
				for (k = 0; k < nnz0; ++k) { \
					if ((d = pj0[k] - pi0[k]) >= a && d <= b) { \
						*(pi1++) = pi0[k]; \
						*(pj1++) = pj0[k]; \
						c##IF_NPATTERN( \
						if (ct0 == 'C' && d == 0) \
							c##ASSIGN_PROJ_REAL(*px1, px0[k]); \
						else \
							c##ASSIGN_IDEN     (*px1, px0[k]); \
						px1++; \
						); \
					} \
					if (d != 0 && -d >= a && -d <= b) { \
						*(pi1++) = pj0[k]; \
						*(pj1++) = pi0[k]; \
						c##IF_NPATTERN( \
						if (ct0 == 'C') \
							c##ASSIGN_CONJ(*px1, px0[k]); \
						else \
							c##ASSIGN_IDEN(*px1, px0[k]); \
						px1++; \
						); \
					} \
				} \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

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
	UNPROTECT(1);
	return s_from;
}

SEXP R_sparse_band(SEXP s_from, SEXP s_a, SEXP s_b)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

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

	return sparse_band(s_from, class, a, b);
}
