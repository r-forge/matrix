/* C implementation of methods for coerce */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "coerce.h"

SEXP vector_as_dense(SEXP from, const char *zzz,
                     char ul, char ct, char nu,
                     int m, int n, int byrow, SEXP dimnames)
{
	SEXPTYPE tf = TYPEOF(from);
	char cl[] = "...Matrix";
	cl[0] = (zzz[0] == '.') ? typeToKind(tf) : ((zzz[0] == ',') ? ((tf == CPLXSXP) ? 'z' : 'd') : zzz[0]);
	cl[1] = zzz[1];
	cl[2] = zzz[2];
#ifndef MATRIX_ENABLE_IMATRIX
	if (cl[0] == 'i')
		cl[0] = 'd';
#endif
	SEXPTYPE tt = kindToType(cl[0]);
	int packed = cl[2] == 'p';
	PROTECT(from = Rf_coerceVector(from, tt));

	if (cl[1] != 'g' && m != n)
		Rf_error(_("attempt to construct non-square %s"),
		         (cl[1] == 's' || cl[1] == 'p') ? "symmetricMatrix" : "triangularMatrix");

	int_fast64_t mn = (int_fast64_t) m * n,
		xlen = (!packed) ? mn : n + (mn - n) / 2;
	if (xlen > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	SEXP to = PROTECT(newObject(cl));
	SET_DIM(to, m, n);
	if (dimnames != R_NilValue)
	SET_DIMNAMES(to, -(cl[1] == 's' || cl[1] == 'p'), dimnames);
	if (cl[1] != 'g' && ul != 'U')
		SET_UPLO(to);
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z')
		SET_TRANS(to);
	if (cl[1] == 't' && nu != 'N')
		SET_DIAG(to);

	/* FIXME: add argument 'new' and conditionally avoid allocation */
	SEXP x = PROTECT(Rf_allocVector(tt, (R_xlen_t) xlen));
	R_xlen_t k, r = XLENGTH(from);
	int i, j, recycle = r < mn;

#define TEMPLATE(c) \
	do { \
		c##TYPE *dest = c##PTR(x), *src = c##PTR(from); \
		if (r == 0) { \
			while (mn-- > 0) \
				*(dest++) = c##NA; \
		} else if (r == 1) { \
			while (mn-- > 0) \
				*(dest++) = *src; \
		} else if (!packed) { \
			if (!recycle) \
				c##NAME(trans2)(dest, src, (size_t) n, (size_t) m, \
				                (!byrow) ? 'N' : 'T'); \
			else { \
				if (!byrow) { \
					k = 0; \
					while (mn-- > 0) { \
						k %= r; \
						*(dest++) = src[k++]; \
					} \
				} else { \
					for (j = 0; j < n; ++j) { \
						k = j; \
						for (i = 0; i < m; ++i) { \
							k %= r; \
							*(dest++) = src[k]; \
							k += n; \
						} \
					} \
				} \
			} \
		} else if (ul == 'U') { \
			if (!byrow) { \
				k = 0; \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i <= j; ++i) { \
						if (recycle) k %= r; \
						*(dest++) = src[k++]; \
					} \
					k += n - j - 1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					k = j; \
					for (i = 0; i <= j; ++i) { \
						if (recycle) k %= r; \
						*(dest++) = src[k]; \
						k += n; \
					} \
				} \
			} \
		} else { \
			if (!byrow) { \
				k = 0; \
				for (j = 0; j < n; ++j) { \
					for (i = j; i < n; ++i) { \
						if (recycle) k %= r; \
						*(dest++) = src[k++]; \
					} \
					k += j + 1; \
				} \
			} else { \
				R_xlen_t d = 0; \
				for (j = 0; j < n; ++j) { \
					k = j + d; \
					for (i = 0; i <= j; ++i) { \
						if (recycle) k %= r; \
						*(dest++) = src[k]; \
						k += n; \
					} \
					d += n; \
				} \
			} \
		} \
	} while (0)

	SWITCH4(cl[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(3); /* x, to, from */
	return to;
}

SEXP R_vector_as_dense(SEXP s_from, SEXP s_zzz,
                       SEXP s_uplo, SEXP s_trans, SEXP s_diag,
                       SEXP s_m, SEXP s_n, SEXP s_byrow,
                       SEXP s_dimnames)
{
	switch (TYPEOF(s_from)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		ERROR_INVALID_TYPE(s_from, __func__);
		break;
	}

	const char *zzz;
	if (TYPEOF(s_zzz) != STRSXP || LENGTH(s_zzz) < 1 ||
	    (s_zzz = STRING_ELT(s_zzz, 0)) == NA_STRING ||
	    (zzz = CHAR(s_zzz))[0] == '\0' ||
	    (zzz              )[1] == '\0' ||
	    !((zzz[1] == 'g' && (zzz[2] == 'e'                  )) ||
	      (zzz[1] == 's' && (zzz[2] == 'y' || zzz[2] == 'p')) ||
	      (zzz[1] == 'p' && (zzz[2] == 'o' || zzz[2] == 'p')) ||
	      (zzz[1] == 't' && (zzz[2] == 'r' || zzz[2] == 'p'))))
		Rf_error(_("second argument of '%s' does not specify a subclass of %s"),
		         __func__, "denseMatrix");

	char ul = '\0', ct = '\0', nu = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , nu);

	int m = -1;
	if (s_m != R_NilValue) {
		if (TYPEOF(s_m) == INTSXP) {
			int tmp;
			if (LENGTH(s_m) >= 1 && (tmp = INTEGER(s_m)[0]) != NA_INTEGER &&
			    tmp >= 0)
				m = tmp;
		} else if (TYPEOF(s_m) == REALSXP) {
			double tmp;
			if (LENGTH(s_m) >= 1 && !ISNAN(tmp = REAL(s_m)[0]) &&
			    tmp >= 0.0) {
				if (trunc(tmp) > INT_MAX)
					Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
				m = (int) tmp;
			}
		}
		if (m < 0)
			Rf_error(_("invalid '%s' to '%s'"), "m", __func__);
	}

	int n = -1;
	if (s_n != R_NilValue) {
		if (TYPEOF(s_n) == INTSXP) {
			int tmp;
			if (LENGTH(s_n) >= 1 && (tmp = INTEGER(s_n)[0]) != NA_INTEGER &&
			    tmp >= 0)
				n = tmp;
		} else if (TYPEOF(s_n) == REALSXP) {
			double tmp;
			if (LENGTH(s_n) >= 1 && !ISNAN(tmp = REAL(s_n)[0]) &&
			    tmp >= 0.0) {
				if (trunc(tmp) > INT_MAX)
					Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
				n = (int) tmp;
			}
		}
		if (n < 0)
			Rf_error(_("invalid '%s' to '%s'"), "n", __func__);
	}

	int byrow;
	VALID_LOGIC2(s_byrow, byrow);

	if (s_dimnames != R_NilValue)
		if (TYPEOF(s_dimnames) != VECSXP || LENGTH(s_dimnames) != 2)
			Rf_error(_("invalid '%s' to '%s'"), "dimnames", __func__);

	R_xlen_t vlen = XLENGTH(s_from);
	if (zzz[1] != 'g' && (m < 0) != (n < 0)) {
		if (m < 0)
			m = n;
		else
			n = m;
	} else if (m < 0 && n < 0) {
		if (vlen > INT_MAX)
			Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		m = (int) vlen;
		n = 1;
	} else if (m < 0) {
		if (vlen > (int_fast64_t) INT_MAX * n) {
			if (n == 0)
				Rf_error(_("nonempty vector supplied for empty matrix"));
			else
				Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		m = (int) ((n == 0) ? 0 : vlen / n + (vlen % n != 0));
	} else if (n < 0) {
		if (vlen > (int_fast64_t) m * INT_MAX) {
			if (m == 0)
				Rf_error(_("nonempty vector supplied for empty matrix"));
			else
				Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		n = (int) ((m == 0) ? 0 : vlen / m + (vlen % m != 0));
	}

	int_fast64_t mlen = (int_fast64_t) m * n;
	if (vlen <= 1)
		/* do nothing */ ;
	else if (mlen == 0)
		Rf_warning(_("nonempty vector supplied for empty matrix"));
	else if (vlen > mlen)
		Rf_warning(_("vector length (%lld) exceeds matrix length (%d * %d)"),
		           (long long) vlen, m, n);
	else if (mlen % vlen != 0)
		Rf_warning(_("matrix length (%d * %d) is not a multiple of vector length (%lld)"),
		           m, n, (long long) vlen);

	return
	vector_as_dense(s_from, zzz, ul, ct, nu, m, n, byrow, s_dimnames);
}

SEXP matrix_as_dense(SEXP from, const char *zzz,
                     char ul, char ct, char nu, int mg, int new)
{
	SEXPTYPE tf = TYPEOF(from);
	char cl[] = "...Matrix";
	cl[0] = (zzz[0] == '.') ? typeToKind(tf) : ((zzz[0] == ',') ? ((tf == CPLXSXP) ? 'z' : 'd') : zzz[0]);
	cl[1] = zzz[1];
	cl[2] = zzz[2];
#ifndef MATRIX_ENABLE_IMATRIX
	if (cl[0] == 'i')
		cl[0] = 'd';
#endif
	SEXPTYPE tt = kindToType(cl[0]);
	int packed = cl[2] == 'p';
	PROTECT(from = Rf_coerceVector(from, tt));

	SEXP to = PROTECT(newObject(cl));
	int nprotect = 2;

	SEXP dim = Rf_getAttrib(from, R_DimSymbol), dimnames;
	int *pdim, isM, m, n, doDN;
	R_xlen_t mn = XLENGTH(from);

	isM = TYPEOF(dim) == INTSXP && LENGTH(dim) == 2;
	if (isM) {

		pdim = INTEGER(dim);
		m = pdim[0];
		n = pdim[1];
		pdim = DIM(to);
		pdim[0] = m;
		pdim[1] = n;

		PROTECT(dimnames = Rf_getAttrib(from, R_DimNamesSymbol));
		++nprotect;
		doDN = dimnames != R_NilValue;

	} else {

		if (mn > INT_MAX)
			Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		pdim = DIM(to);
		if (mg == 0) {
			pdim[0] = m = 1;
			pdim[1] = n = (int) mn;
		} else {
			pdim[0] = m = (int) mn;
			pdim[1] = n = 1;
		}

		SEXP nms = PROTECT(Rf_getAttrib(from, R_NamesSymbol));
		++nprotect;
		doDN = nms != R_NilValue;
		if (doDN) {
			PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
			++nprotect;
			SET_VECTOR_ELT(dimnames, mg == 0, nms);
		}

	}

	if (cl[1] != 'g' && m != n)
		Rf_error(_("attempt to construct non-square %s"),
		         (cl[1] == 's' || cl[1] == 'p') ? "symmetricMatrix" : "triangularMatrix");

	if (doDN)
		SET_DIMNAMES(to, -(cl[1] == 's' || cl[1] == 'p'), dimnames);
	if (cl[1] != 'g' && ul != 'U')
		SET_UPLO(to);
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z')
		SET_TRANS(to);
	if (cl[1] == 't' && nu != 'N')
		SET_DIAG(to);

	SEXP x;

	if (!packed) {

	if (new <= 0 || (new <= 1 && !ANY_ATTRIB(from)) ||
	    !MAYBE_REFERENCED(from)) {

		if (new > 0 && ANY_ATTRIB(from))
			/* 'from' has attributes and no references : */
			CLEAR_ATTRIB(from);
		x = from;

	} else {

		PROTECT(x = duplicateVector(from));
		++nprotect;

	}

	} else {

		PROTECT(x = Rf_allocVector(tt, n + (mn - n) / 2));
		++nprotect;

#define TEMPLATE(c) \
		c##NAME(pack2)(c##PTR(x), c##PTR(from), (size_t) n, \
		               ul, '\0', 'N')

		SWITCH4(cl[0], TEMPLATE);

#undef TEMPLATE

	}

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(nprotect);
	return to;
}

/* as(<matrix>, ".(ge|sy|sp|po|pp|tr|tp)Matrix") */
SEXP R_matrix_as_dense(SEXP s_from, SEXP s_zzz,
                       SEXP s_uplo, SEXP s_trans, SEXP s_diag, SEXP s_margin)
{
	switch (TYPEOF(s_from)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		ERROR_INVALID_TYPE(s_from, __func__);
		break;
	}

	const char *zzz;
	if (TYPEOF(s_zzz) != STRSXP || LENGTH(s_zzz) < 1 ||
	    (s_zzz = STRING_ELT(s_zzz, 0)) == NA_STRING ||
	    (zzz = CHAR(s_zzz))[0] == '\0' ||
	    (zzz              )[1] == '\0' ||
	    !((zzz[1] == 'g' && (zzz[2] == 'e'                  )) ||
	      (zzz[1] == 's' && (zzz[2] == 'y' || zzz[2] == 'p')) ||
	      (zzz[1] == 'p' && (zzz[2] == 'o' || zzz[2] == 'p')) ||
	      (zzz[1] == 't' && (zzz[2] == 'r' || zzz[2] == 'p'))))
		Rf_error(_("second argument of '%s' does not specify a subclass of %s"),
		         __func__, "denseMatrix");

	char ul = '\0', ct = '\0', nu = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , nu);

	int mg;
	VALID_MARGIN(s_margin, mg);

	return matrix_as_dense(s_from, zzz, ul, ct, nu, mg, 1);
}

SEXP sparse_as_dense(SEXP from, const char *class, int packed)
{
	packed = packed && class[1] != 'g';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = (packed) ? 'p' : ((class[1] == 'g') ? 'e' : ((class[1] == 's') ? 'y' : ((class[1] == 'p') ? 'o' : 'r')));
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n,
		xlen = (!packed) ? mn : n + (mn - n) / 2;
	if (xlen > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");
	if (class[2] != 'C' && packed && mn > R_XLEN_T_MAX)
		Rf_error(_("coercing n-by-n %s to %s is not supported for n*n exceeding %s"),
		         "[RT]sparseMatrix", "packedMatrix", "R_XLEN_T_MAX");
	double bytes = (double) xlen * kindToSize(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		Rf_warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		           0x1.0p-30 * bytes);
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && (nu = DIAG(from)) != 'N')
		SET_DIAG(to);

	/* It remains to fill 'x' ... */

	SEXP p0, i0, j0;
	int *pp0, *pi0, *pj0, nprotect = 1;
	p0 = i0 = j0 = NULL;
	pp0 = pi0 = pj0 = NULL;
	if (class[2] != 'T') {
		PROTECT(p0 = GET_SLOT(from, Matrix_pSym));
		++nprotect;
		pp0 = INTEGER(p0) + 1;
	}
	if (class[2] != 'R') {
		PROTECT(i0 = GET_SLOT(from, Matrix_iSym));
		++nprotect;
		pi0 = INTEGER(i0);
	}
	if (class[2] != 'C') {
		PROTECT(j0 = GET_SLOT(from, Matrix_jSym));
		++nprotect;
		pj0 = INTEGER(j0);
	}

#define TEMPLATE(c) \
	do { \
		c##IF_NPATTERN( \
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
		c##TYPE *px0 = c##PTR(x0); \
		); \
		SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, (R_xlen_t) xlen)); \
		c##TYPE *px1 = c##PTR(x1); \
		memset(px1, 0, sizeof(c##TYPE) * (size_t) xlen); \
		switch (class[2]) { \
		case 'C': \
		{ \
			int j, k, kend; \
			if (!packed) \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						px1[*pi0] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pi0; c##IF_NPATTERN(++px0); \
					} \
					px1 += m; \
				} \
			else if (ul == 'U') \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						px1[*pi0] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pi0; c##IF_NPATTERN(++px0); \
					} \
					px1 += j + 1; \
				} \
			else \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						px1[*pi0 - j] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pi0; c##IF_NPATTERN(++px0); \
					} \
					px1 += n - j; \
				} \
			break; \
		} \
		case 'R': \
		{ \
			int i, k, kend; \
			int_fast64_t index, m1 = (int_fast64_t) m; \
			if (!packed) \
				for (i = 0, k = 0; i < m; ++i) { \
					kend = pp0[i]; \
					while (k < kend) { \
						index = m1 * *pj0; \
						px1[index] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pj0; c##IF_NPATTERN(++px0); \
					} \
					px1 += 1; \
				} \
			else if (ul == 'U') \
				for (i = 0, k = 0; i < m; ++i) { \
					kend = pp0[i]; \
					while (k < kend) { \
						index = DENSE_INDEX_U((int_fast64_t) i, (int_fast64_t) *pj0, m1); \
						px1[index] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pj0; c##IF_NPATTERN(++px0); \
					} \
				} \
			else \
				for (i = 0, k = 0; i < m; ++i) { \
					kend = pp0[i]; \
					while (k < kend) { \
						index = DENSE_INDEX_L((int_fast64_t) i, (int_fast64_t) *pj0, m1); \
						px1[index] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pj0; c##IF_NPATTERN(++px0); \
					} \
				} \
			break; \
		} \
		case 'T': \
		{ \
			R_xlen_t k, kend = XLENGTH(i0); \
			int_fast64_t index, m1 = (int_fast64_t) m; \
			if (!packed) \
				for (k = 0; k < kend; ++k) { \
					index = m1 * *pj0 + *pi0; \
					c##INCREMENT_IDEN(px1[index], *px0); \
					++pi0; ++pj0; c##IF_NPATTERN(++px0); \
				} \
			else if (ul == 'U') \
				for (k = 0; k < kend; ++k) { \
					index = DENSE_INDEX_U((int_fast64_t) *pi0, (int_fast64_t) *pj0, m1); \
					c##INCREMENT_IDEN(px1[index], *px0); \
					++pi0; ++pj0; c##IF_NPATTERN(++px0); \
				} \
			else \
				for (k = 0; k < kend; ++k) { \
					index = DENSE_INDEX_L((int_fast64_t) *pi0, (int_fast64_t) *pj0, m1); \
					c##INCREMENT_IDEN(px1[index], *px0); \
					++pi0; ++pj0; c##IF_NPATTERN(++px0); \
				} \
			break; \
		} \
		default: \
			break; \
		} \
		SET_SLOT(to, Matrix_xSym, x1); \
		UNPROTECT(1); /* x1 */ \
		c##IF_NPATTERN( \
		UNPROTECT(1); /* x0 */ \
		); \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

	UNPROTECT(nprotect);
	return to;
}

/* as(<[CRT]sparseMatrix>, "(un)?packedMatrix") */
SEXP R_sparse_as_dense(SEXP s_from, SEXP s_packed)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);

	int packed = 0;
	if (class[1] != 'g') VALID_LOGIC2(s_packed, packed);

	return sparse_as_dense(s_from, class, packed);
}

SEXP diagonal_as_dense(SEXP from, const char *class,
                       char kind, char shape, int packed,
                       char ul, char ct)
{
	packed = packed && shape != 'g';

	char cl[] = "...Matrix";
	cl[0] = (kind == '.') ? class[0] : ((kind == ',') ? ((class[0] == 'z') ? 'z' : 'd') : kind);
	cl[1] = shape;
	cl[2] = (packed) ? 'p' : ((shape == 'g') ? 'e' : ((shape == 's') ? 'y' : ((shape == 'p') ? 'o' : 'r')));
	SEXP to = PROTECT(newObject(cl));

	int n = DIM(from)[0];
	int_fast64_t nn = (int_fast64_t) n * n,
		xlen = (!packed) ? nn : n + (nn - n) / 2;
	if (xlen > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");
	double bytes = (double) xlen * kindToSize(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		Rf_warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		           0x1.0p-30 * bytes);
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(cl[1] == 's' || cl[1] == 'p'), DIMNAMES(from, 0));
	if (cl[1] != 'g' && ul != 'U')
		SET_UPLO(to);
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z')
		SET_TRANS(to);

	char nu = DIAG(from);
	if (cl[1] == 't' && nu != 'N')
		SET_DIAG(to);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (cl[0] != class[0]) {
		if (class[0] == 'n' && cl[0] == 'l')
			x0 = duplicateVector(x0);
		else
			x0 = Rf_coerceVector(x0, kindToType(cl[0]));
		if (class[0] == 'n')
			naToUnit(x0);
		UNPROTECT(1); /* x0 */
		PROTECT(x0);
	}

	SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), (R_xlen_t) xlen));

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed) \
			c##NAME(force2)(px1, px0, (size_t) n, ul, 1, -nu); \
		else \
			c##NAME(force1)(px1, px0, (size_t) n, ul, 1, -nu); \
	} while (0)

	SWITCH4(cl[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<diagonalMatrix>, ".(ge|sy|sp|po|pp|tr|tp)Matrix") */
SEXP R_diagonal_as_dense(SEXP s_from,
                         SEXP s_kind, SEXP s_shape, SEXP s_packed,
                         SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_diagonal, 0, __func__);

	char kind, shape;
	VALID_KIND (s_kind , kind );
	VALID_SHAPE(s_shape, shape);

	int packed = 0;
	if (shape != 'g') VALID_LOGIC2(s_packed, packed);

	char ul = '\0', ct = '\0';
	if (shape != 'g') VALID_UPLO (s_uplo , ul);
	if (shape == 's') VALID_TRANS(s_trans, ct);

	return diagonal_as_dense(s_from, class, kind, shape, packed, ul, ct);
}

SEXP index_as_dense(SEXP from, const char *class, char kind)
{
	int mg = MARGIN(from);

	char cl[] = ".geMatrix";
	cl[0] = (kind == '.') ? 'n' : ((kind == ',') ? 'd' : kind);
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	int_fast64_t xlen = (int_fast64_t) m * n;
	if (xlen > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");
	double bytes = (double) xlen * kindToSize(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		Rf_warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		           0x1.0p-30 * bytes);
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	SEXP perm = PROTECT(GET_SLOT(from, Matrix_permSym));
	int *pperm = INTEGER(perm);

	SEXP x = PROTECT(Rf_allocVector(kindToType(cl[0]), (R_xlen_t) xlen));

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x); \
		memset(px, 0, sizeof(c##TYPE) * (size_t) xlen); \
		if (mg == 0) { \
			int_fast64_t m1 = (int_fast64_t) m; \
			for (int i = 0; i < m; ++i) { \
				px[m1 * (*(pperm++) - 1)] = c##UNIT; \
				px += 1; \
			} \
		} else { \
			for (int j = 0; j < n; ++j) { \
				px[      *(pperm++) - 1 ] = c##UNIT; \
				px += m; \
			} \
		} \
	} while (0)

	SWITCH4(cl[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(3); /* x, perm, to */
	return to;
}

/* as(<indMatrix>, ".geMatrix") */
SEXP R_index_as_dense(SEXP s_from, SEXP s_kind)
{
	const char *class = Matrix_class(s_from, valid_index, 0, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	return index_as_dense(s_from, class, kind);
}

SEXP Vector_as_sparse(SEXP from, const char *zzz,
                      char ul, char ct, char nu,
                      int m, int n, int byrow, SEXP dimnames)
{
	SEXP length0 = GET_SLOT(from, Matrix_lengthSym);
	int_fast64_t r = (int_fast64_t)
		((TYPEOF(length0) == INTSXP) ? INTEGER(length0)[0] : REAL(length0)[0]);

	SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
		x0 = PROTECT(Rf_getAttrib(from, Matrix_xSym));

	SEXPTYPE tf = TYPEOF(x0);
	char cl[] = "...Matrix";
	cl[0] = (zzz[0] == '.')
		? ((x0 == R_NilValue) ? 'n' : typeToKind(tf))
		: ((zzz[0] == ',') ? ((tf == CPLXSXP) ? 'z' : 'd') : zzz[0]);
	cl[1] = zzz[1];
	cl[2] = (byrow) ? 'R' : 'C';
#ifndef MATRIX_ENABLE_IMATRIX
	if (cl[0] == 'i')
		cl[0] = 'd';
#endif
	SEXPTYPE tt = kindToType(cl[0]);
	if (x0 != R_NilValue && cl[0] != 'n') {
		x0 = Rf_coerceVector(x0, tt);
		UNPROTECT(1); /* x0 */
		PROTECT(x0);
	}

	if (cl[1] != 'g' && m != n)
		Rf_error(_("attempt to construct non-square %s"),
		         (cl[1] == 's' || cl[1] == 'p') ? "symmetricMatrix" : "triangularMatrix");

	SEXP to = PROTECT(newObject(cl));
	SET_DIM(to, m, n);
	if (dimnames != R_NilValue)
	SET_DIMNAMES(to, -(cl[1] == 's' || cl[1] == 'p'), dimnames);
	if (cl[1] != 'g' && ul != 'U')
		SET_UPLO(to);
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z')
		SET_TRANS(to);
	if (cl[1] == 't' && nu != 'N')
		SET_DIAG(to);

	int_fast64_t pos, mn = (int_fast64_t) m * n, nnz1 = 0;
	R_xlen_t k = 0, nnz0 = XLENGTH(i0);

#define TEMPLATE__(d) \
	do { \
		d##TYPE *pi0 = d##PTR(i0); \
		if (nnz0 == 0) \
			/* do nothing */ ; \
		else if (cl[1] == 'g') { \
			if (r == 0) \
				nnz1 = mn; \
			else if (r == mn) \
				nnz1 = nnz0; \
			else if (r > mn) \
				while (k < nnz0 && (int_fast64_t) pi0[k++] <= mn) \
					++nnz1; \
			else { \
				int_fast64_t mn_mod_r = mn % r; \
				nnz1 = nnz0 * (mn / r); \
				while (k < nnz0 && (int_fast64_t) pi0[k++] <= mn_mod_r) \
					++nnz1; \
			} \
		} \
		else if (cl[1] == 's' || cl[1] == 'p' || nu == 'N') { \
			if (r == 0) \
				nnz1 = n + (mn - n) / 2; \
			else if (r >= mn) { \
				if ((ul == 'U') == !byrow) { \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n <= pos / n) \
						++nnz1; \
				} else { \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n >= pos / n) \
						++nnz1; \
				} \
			} \
			else { \
				int_fast64_t a = 0; \
				if ((ul == 'U') == !byrow) \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n <= pos / n) \
						++nnz1; \
				a += r; \
				} \
				else \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n >= pos / n) \
						++nnz1; \
				a += r; \
				} \
			} \
		} \
		else { \
			if (r == 0) \
				nnz1 = (mn - n) / 2; \
			else if (r >= mn) { \
				if ((ul == 'U') == !byrow) { \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n < pos / n) \
						++nnz1; \
				} else { \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n > pos / n) \
						++nnz1; \
				} \
			} \
			else { \
				int_fast64_t a = 0; \
				if ((ul == 'U') == !byrow) \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n < pos / n) \
						++nnz1; \
				a += r; \
				} \
				else \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k++] - 1) < mn) \
					if (pos % n > pos / n) \
						++nnz1; \
				a += r; \
				} \
			} \
		} \
	} while (0)

	if (TYPEOF(i0) == INTSXP)
		TEMPLATE__(i);
	else
		TEMPLATE__(d);

#undef TEMPLATE__

	if (nnz1 > INT_MAX)
		Rf_error(_("attempt to construct %s with more than %s nonzero entries"),
		         "sparseMatrix", "2^31-1");

	if (byrow)
		SWAP(m, n, int, );

	SEXP iSym = (!byrow) ? Matrix_iSym : Matrix_jSym,
		p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) n + 1)),
		i1 = PROTECT(Rf_allocVector(INTSXP, nnz1));
	int *pp1 = INTEGER(p1) + 1, *pi1 = INTEGER(i1), i, j;
	memset(pp1 - 1, 0, sizeof(int) * ((size_t) n + 1));
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to,        iSym, i1);

	k = 0;

#define TEMPLATE__(d, c0, c1) \
	do { \
		d ##TYPE *pi0 = d ##PTR(i0); \
		c0##IF_NPATTERN( \
		c0##TYPE *px0 = c0##PTR(x0); \
		); \
		c1##IF_NPATTERN( \
		SEXP x1 = PROTECT(Rf_allocVector(c1##TYPESXP, nnz1)); \
		c1##TYPE *px1 = c1##PTR(x1); \
		); \
		if (nnz1 == 0) \
			/* do nothing */ ; \
		else if (cl[1] == 'g') { \
			if (r == 0) \
				for (j = 0; j < n; ++j) { \
					pp1[j] = m; \
					for (i = 0; i < m; ++i) { \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c1##NA; \
						); \
					} \
				} \
			else if (r >= mn) \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k] - 1) < mn) { \
					++pp1[pos / m]; \
					*(pi1++) = pos % m; \
					c1##IF_NPATTERN( \
					*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
					); \
					++k; \
				} \
			else { \
				int_fast64_t a = 0; \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k] - 1) < mn) { \
					++pp1[pos / m]; \
					*(pi1++) = pos % m; \
					c1##IF_NPATTERN( \
					*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
					); \
					++k; \
				} \
				a += r; \
				} \
			} \
		} \
		else if (cl[1] == 's' || cl[1] == 'p' || nu == 'N') { \
			if (r == 0) { \
				if ((ul == 'U') == !byrow) \
				for (j = 0; j < n; ++j) { \
					pp1[j] = j + 1; \
					for (i = 0; i <= j; ++i) { \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c1##NA; \
						); \
					} \
				} \
				else \
				for (j = 0; j < n; ++j) { \
					pp1[j] = n - j; \
					for (i = j; i < n; ++i) { \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c1##NA; \
						); \
					} \
				} \
			} \
			else if (r >= mn) { \
				if ((ul == 'U') == !byrow) \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) <= (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
				else \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) >= (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
			} \
			else { \
				int_fast64_t a = 0; \
				if ((ul == 'U') == !byrow) \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) <= (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
				a += r; \
				} \
				else \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) >= (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
				a += r; \
				} \
			} \
		} \
		else { \
			if (r == 0) { \
				if ((ul == 'U') == !byrow) \
				for (j = 0; j < n; ++j) { \
					pp1[j] = j; \
					for (i = 0; i < j; ++i) { \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c1##NA; \
						); \
					} \
				} \
				else \
				for (j = 0; j < n; ++j) { \
					pp1[j] = n - j - 1; \
					for (i = j + 1; i < n; ++i) { \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c1##NA; \
						); \
					} \
				} \
			} \
			else if (r >= mn) { \
				if ((ul == 'U') == !byrow) \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) < (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
				else \
				while (k < nnz0 && (pos =     (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) > (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
			} \
			else { \
				int_fast64_t a = 0; \
				if ((ul == 'U') == !byrow) \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) < (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
				a += r; \
				} \
				else \
				while (a < mn) { \
				k = 0; \
				while (k < nnz0 && (pos = a + (int_fast64_t) pi0[k] - 1) < mn) { \
					if ((i = (int) pos % n) > (j = (int) pos / n)) { \
						++pp1[j]; \
						*(pi1++) = i; \
						c1##IF_NPATTERN( \
						*(px1++) = c0##IFELSE_NPATTERN(px0[k], c1##UNIT); \
						); \
					} \
					++k; \
				} \
				a += r; \
				} \
			} \
		} \
		c1##IF_NPATTERN( \
		SET_SLOT(to, Matrix_xSym, x1); \
		UNPROTECT(1); /* x1 */ \
		); \
	} while (0)

	if (x0 == R_NilValue) {

#define TEMPLATE(c) \
		do { \
			if (TYPEOF(i0) == INTSXP) \
				TEMPLATE__(i, n, c); \
			else \
				TEMPLATE__(d, n, c); \
		} while (0)

		SWITCH5(cl[0], TEMPLATE);

#undef TEMPLATE

	} else {

#define TEMPLATE(c) \
		do { \
			if (TYPEOF(i0) == INTSXP) \
				TEMPLATE__(i, c, c); \
			else \
				TEMPLATE__(d, c, c); \
		} while (0)

		SWITCH5(cl[0], TEMPLATE);

#undef TEMPLATE

	}

#undef TEMPLATE__

	for (j = 0; j < n; ++j)
		pp1[j] += pp1[j - 1];

	switch (zzz[2]) {
	case 'C':
		to = sparse_as_Csparse(to, cl);
		break;
	case 'R':
		to = sparse_as_Rsparse(to, cl);
		break;
	case 'T':
		to = sparse_as_Tsparse(to, cl);
		break;
	default:
		break;
	}

	UNPROTECT(5); /* i1, p1, to, x0, i0 */
	return to;
}

SEXP R_Vector_as_sparse(SEXP s_from, SEXP s_zzz,
                        SEXP s_uplo, SEXP s_trans, SEXP s_diag,
                        SEXP s_m, SEXP s_n, SEXP s_byrow,
                        SEXP s_dimnames)
{
	Matrix_class(s_from, valid_vector, 0, __func__);

	const char *zzz;
	if (TYPEOF(s_zzz) != STRSXP || LENGTH(s_zzz) < 1 ||
	    (s_zzz = STRING_ELT(s_zzz, 0)) == NA_STRING ||
	    (zzz = CHAR(s_zzz))[0] == '\0' ||
	    (zzz[1] != 'g' && zzz[1] != 's' && zzz[1] != 'p' && zzz[1] != 't') ||
	    (zzz[2] != 'C' && zzz[2] != 'R' && zzz[2] != 'T'))
		Rf_error(_("second argument of '%s' does not specify a subclass of %s"),
		         __func__, "[CRT]sparseMatrix");

	char ul = '\0', ct = '\0', nu = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , nu);

	int m = -1;
	if (s_m != R_NilValue) {
		if (TYPEOF(s_m) == INTSXP) {
			int tmp;
			if (LENGTH(s_m) >= 1 && (tmp = INTEGER(s_m)[0]) != NA_INTEGER &&
			    tmp >= 0)
				m = tmp;
		} else if (TYPEOF(s_m) == REALSXP) {
			double tmp;
			if (LENGTH(s_m) >= 1 && !ISNAN(tmp = REAL(s_m)[0]) &&
			    tmp >= 0.0) {
				if (trunc(tmp) > INT_MAX)
					Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
				m = (int) tmp;
			}
		}
		if (m < 0)
			Rf_error(_("invalid '%s' to '%s'"), "m", __func__);
	}

	int n = -1;
	if (s_n != R_NilValue) {
		if (TYPEOF(s_n) == INTSXP) {
			int tmp;
			if (LENGTH(s_n) >= 1 && (tmp = INTEGER(s_n)[0]) != NA_INTEGER &&
			    tmp >= 0)
				n = tmp;
		} else if (TYPEOF(s_n) == REALSXP) {
			double tmp;
			if (LENGTH(s_n) >= 1 && !ISNAN(tmp = REAL(s_n)[0]) &&
			    tmp >= 0.0) {
				if (trunc(tmp) > INT_MAX)
					Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
				n = (int) tmp;
			}
		}
		if (n < 0)
			Rf_error(_("invalid '%s' to '%s'"), "n", __func__);
	}

	int byrow;
	VALID_LOGIC2(s_byrow, byrow);

	if (s_dimnames != R_NilValue)
		if (TYPEOF(s_dimnames) != VECSXP || LENGTH(s_dimnames) != 2)
			Rf_error(_("invalid '%s' to '%s'"), "dimnames", __func__);

	SEXP tmp = GET_SLOT(s_from, Matrix_lengthSym);
	int_fast64_t vlen = (int_fast64_t)
		((TYPEOF(tmp) == INTSXP) ? INTEGER(tmp)[0] : REAL(tmp)[0]);
	if (zzz[1] != 'g' && (m < 0) != (n < 0)) {
		if (m < 0)
			m = n;
		else
			n = m;
	} else if (m < 0 && n < 0) {
		if (vlen > INT_MAX)
			Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		m = (int) vlen;
		n = 1;
	} else if (m < 0) {
		if (vlen > (int_fast64_t) INT_MAX * n) {
			if (n == 0)
				Rf_error(_("nonempty vector supplied for empty matrix"));
			else
				Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		m = (int) ((n == 0) ? 0 : vlen / n + (vlen % n != 0));
	} else if (n < 0) {
		if (vlen > (int_fast64_t) m * INT_MAX) {
			if (m == 0)
				Rf_error(_("nonempty vector supplied for empty matrix"));
			else
				Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		n = (int) ((m == 0) ? 0 : vlen / m + (vlen % m != 0));
	}

	int_fast64_t mlen = (int_fast64_t) m * n;
	if (vlen <= 1)
		/* do nothing */ ;
	else if (mlen == 0)
		Rf_warning(_("nonempty vector supplied for empty matrix"));
	else if (vlen > mlen)
		Rf_warning(_("vector length (%lld) exceeds matrix length (%d * %d)"),
		           (long long) vlen, m, n);
	else if (mlen % vlen != 0)
		Rf_warning(_("matrix length (%d * %d) is not a multiple of vector length (%lld)"),
		           m, n, (long long) vlen);

	return
	Vector_as_sparse(s_from, zzz, ul, ct, nu, m, n, byrow, s_dimnames);
}

SEXP matrix_as_sparse(SEXP from, const char *zzz,
                      char ul, char ct, char nu, int mg)
{
	char cl[] = "...Matrix";
	cl[0] = typeToKind(TYPEOF(from));
	cl[1] = zzz[1];
	cl[2] = (zzz[1] == 'g') ? 'e' : ((zzz[1] == 's') ? 'y' : ((zzz[1] == 'p') ? 'o' : 'r'));
#ifndef MATRIX_ENABLE_IMATRIX
	if (cl[0] == 'i')
		cl[0] = 'd';
#endif
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);
	REPROTECT(from = matrix_as_dense(from, cl, ul, ct, nu, mg, 0), pid);
	REPROTECT(from = dense_as_sparse(from, cl, zzz[2]), pid);
	cl[2] = zzz[2];
	REPROTECT(from = sparse_as_kind(from, cl, zzz[0]), pid);
	UNPROTECT(1); /* from */
	return from;
}

/* as(<matrix>, ".[gspt][CRT]Matrix") */
SEXP R_matrix_as_sparse(SEXP s_from, SEXP s_zzz,
                        SEXP s_uplo, SEXP s_trans, SEXP s_diag, SEXP s_margin)
{
	switch (TYPEOF(s_from)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		ERROR_INVALID_TYPE(s_from, __func__);
		break;
	}

	const char *zzz;
	if (TYPEOF(s_zzz) != STRSXP || LENGTH(s_zzz) < 1 ||
	    (s_zzz = STRING_ELT(s_zzz, 0)) == NA_STRING ||
	    (zzz = CHAR(s_zzz))[0] == '\0' ||
	    (zzz[1] != 'g' && zzz[1] != 's' && zzz[1] != 'p' && zzz[1] != 't') ||
	    (zzz[2] != 'C' && zzz[2] != 'R' && zzz[2] != 'T'))
		Rf_error(_("second argument of '%s' does not specify a subclass of %s"),
		         __func__, "[CRT]sparseMatrix");

	char ul = '\0', ct = '\0', nu = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , nu);

	int mg;
	VALID_MARGIN(s_margin, mg);

	return matrix_as_sparse(s_from, zzz, ul, ct, nu, mg);
}

SEXP dense_as_sparse(SEXP from, const char *class, char repr)
{
	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = repr;
#ifndef MATRIX_ENABLE_POSDEF
	if (cl[1] == 'p')
		cl[1] = 's';
#endif
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && (nu = DIAG(from)) != 'N')
		SET_DIAG(to);

	SEXP p1, i1, j1;
	int i, j, *pp1, *pi1, *pj1;
	p1 = i1 = j1 = NULL;
	pp1 = pi1 = pj1 = NULL;

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	R_xlen_t nnz = 0;
	int nprotect = 2;

	if (cl[2] != 'T') {
		int r = (cl[2] == 'C') ? n : m;
		PROTECT(p1 = Rf_allocVector(INTSXP, (R_xlen_t) r + 1));
		++nprotect;
		pp1 = INTEGER(p1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);
	}

#define DAS_CHECK \
	do { \
		if (nnz > INT_MAX) \
			Rf_error(_("attempt to construct %s with more than %s nonzero entries"), \
			         "sparseMatrix", "2^31-1"); \
		*(pp1++) = (int) nnz; \
	} while (0)

#define DAS_BYCOL(c, kernel, hook) \
	do { \
		if (class[1] == 'g') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < m; ++i) { \
					if (c##NOT_ZERO(*px0)) \
						kernel; \
					px0 += 1; \
				} \
				hook; \
			} \
		} else if (!packed && (nu == '\0' || nu == 'N')) { \
			if (ul == 'U') \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					for (; i <= j; ++i) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					px0 += m - j - 1; \
					hook; \
				} \
			else \
				for (j = 0; j < n; ++j) { \
					px0 += j; \
					for (i = j; i <= j; ++i) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					for (; i < m; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					hook; \
				} \
		} else if (!packed) { \
			if (ul == 'U') \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					px0 += m - j; \
					hook; \
				} \
			else \
				for (j = 0; j < n; ++j) { \
					px0 += j + 1; \
					for (i = j + 1; i < m; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					hook; \
				} \
		} else if (nu == '\0' || nu == 'N') { \
			if (ul == 'U') \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					for (; i <= j; ++i) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					hook; \
				} \
			else \
				for (j = 0; j < n; ++j) { \
					for (i = j; i <= j; ++i) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					for (; i < m; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					hook; \
				} \
		} else { \
			if (ul == 'U') \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					px0 += 1; \
					hook; \
				} \
			else \
				for (j = 0; j < n; ++j) { \
					px0 += 1; \
					for (i = j + 1; i < m; ++i) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += 1; \
					} \
					hook; \
				} \
		} \
	} while (0)

#define DAS_BYROW(c, kernel, hook) \
	do { \
		int_fast64_t d; \
		if (class[1] == 'g') { \
			d = (int_fast64_t) m * n - 1; \
			for (i = 0; i < m; ++i) { \
				for (j = 0; j < n; ++j) { \
					if (c##NOT_ZERO(*px0)) \
						kernel; \
					px0 += m; \
				} \
				px0 -= d; \
				hook; \
			} \
		} else if (!packed && (nu == '\0' || nu == 'N')) { \
			if (ul == 'U') { \
				d = (int_fast64_t) m * n - 1; \
				for (i = 0; i < m; ++i) { \
					for (j = i; j <= i; ++j) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m; \
					} \
					for (; j < n; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m; \
					} \
					px0 -= (d -= m); \
					hook; \
				} \
			} else { \
				d = -1; \
				for (i = 0; i < m; ++i) { \
					for (j = 0; j < i; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m; \
					} \
					for (; j <= i; ++j) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m; \
					} \
					px0 -= (d += m); \
					hook; \
				} \
			} \
		} else if (!packed) { \
			if (ul == 'U') { \
				d = (int_fast64_t) m * n - 1; \
				for (i = 0; i < m; ++i) { \
					px0 += m; \
					for (j = i + 1; j < n; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m; \
					} \
					px0 -= (d -= m); \
					hook; \
				} \
			} else { \
				d = -1; \
				for (i = 0; i < m; ++i) { \
					for (j = 0; j < i; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m; \
					} \
					px0 += m; \
					px0 -= (d += m); \
					hook; \
				} \
			} \
		} else if (nu == '\0' || nu == 'N')	{ \
			if (ul == 'U') { \
				d = PACKED_LENGTH((int_fast64_t) n) - 1; \
				for (i = 0; i < m; ++i) { \
					for (j = i; j <= i; ++j) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += j + 1; \
					} \
					for (; j < n; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += j + 1; \
					} \
					px0 -= (d -= i + 1); \
					hook; \
				} \
			} else { \
				d = -1; \
				for (i = 0; i < m; ++i) { \
					for (j = 0; j < i; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m - j - 1; \
					} \
					for (; j <= i; ++j) { \
						if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m - j - 1; \
					} \
					px0 -= (d += m - i - 1); \
					hook; \
				} \
			} \
		} else { \
			if (ul == 'U') { \
				d = PACKED_LENGTH((int_fast64_t) n) - 1; \
				for (i = 0; i < m; ++i) { \
					px0 += i + 1; \
					for (j = i + 1; j < n; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += j + 1; \
					} \
					px0 -= (d -= i + 1); \
					hook; \
				} \
			} else { \
				d = -1; \
				for (i = 0; i < m; ++i) { \
					for (j = 0; j < i; ++j) { \
						if (c##NOT_ZERO(*px0)) \
							kernel; \
						px0 += m - j - 1; \
					} \
					px0 += m - i - 1; \
					px0 -= (d += m - i - 1); \
					hook; \
				} \
			} \
		} \
	} while (0)

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0); \
		switch (cl[2]) { \
		case 'C': \
			DAS_BYCOL(c, ++nnz, DAS_CHECK); \
			break; \
		case 'R': \
			DAS_BYROW(c, ++nnz, DAS_CHECK); \
			break; \
		case 'T': \
			DAS_BYCOL(c, ++nnz, ); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

	if (cl[2] != 'R') {
		PROTECT(i1 = Rf_allocVector(INTSXP, nnz));
		++nprotect;
		pi1 = INTEGER(i1);
		SET_SLOT(to, Matrix_iSym, i1);
	}
	if (cl[2] != 'C') {
		PROTECT(j1 = Rf_allocVector(INTSXP, nnz));
		++nprotect;
		pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_jSym, j1);
	}

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
		c##TYPE *px1 = c##PTR(x1); \
		); \
		switch (cl[2]) { \
		case 'C': \
			DAS_BYCOL(c, \
			          do { \
			          	*(pi1++) = i; \
			          	c##IF_NPATTERN(*(px1++) = *px0); \
			          } while (0), ); \
			break; \
		case 'R': \
			DAS_BYROW(c, \
			          do { \
			          	*(pj1++) = j; \
			          	c##IF_NPATTERN(*(px1++) = *px0); \
			          } while (0), ); \
			break; \
		case 'T': \
			DAS_BYCOL(c, \
			          do { \
			          	*(pi1++) = i; \
			          	*(pj1++) = j; \
			          	c##IF_NPATTERN(*(px1++) = *px0); \
			          } while (0), ); \
			break; \
		default: \
			break; \
		} \
		c##IF_NPATTERN( \
		SET_SLOT(to, Matrix_xSym, x1); \
		UNPROTECT(1); /* x1 */ \
		); \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

	UNPROTECT(nprotect);
	return to;
}

/* as(<denseMatrix>, "[CRT]sparseMatrix") */
SEXP R_dense_as_sparse(SEXP s_from, SEXP s_repr)
{
	const char *class = Matrix_class(s_from, valid_dense, 2, __func__);

	char repr;
	VALID_REPR(s_repr, repr, 0);

	return dense_as_sparse(s_from, class, repr);
}

SEXP diagonal_as_sparse(SEXP from, const char *class,
                        char kind, char shape, char repr,
                        char ul, char ct)
{
	char cl[] = "...Matrix";
	cl[0] = (kind == '.') ? class[0] : ((kind == ',') ? ((class[0] == 'z') ? 'z' : 'd') : kind);
	cl[1] = shape;
	cl[2] = repr;
	SEXP to = PROTECT(newObject(cl));

	int n = DIM(from)[1];
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(cl[1] == 's' || cl[1] == 'p'), DIMNAMES(from, 0));

	char nu = DIAG(from);
	if (cl[1] != 'g' && ul != 'U')
		SET_UPLO(to);
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z')
		SET_TRANS(to);
	if (cl[1] == 't' && nu != 'N')
		SET_DIAG(to);

	if (nu != 'N') {
		if (cl[1] == 't') {
			if (cl[2] != 'T') {
				SEXP p1 = PROTECT(allocZero(INTSXP, (R_xlen_t) n + 1));
				SET_SLOT(to, Matrix_pSym, p1);
				UNPROTECT(1); /* p1 */
			}
		} else {
			if (cl[2] != 'T') {
				SEXP p1 = PROTECT(allocSeqInt(0, (R_xlen_t) n + 1));
				SET_SLOT(to, Matrix_pSym, p1);
				UNPROTECT(1); /* p1 */
			}
			SEXP i1 = PROTECT(allocSeqInt(0, n));
			if (cl[2] != 'R')
				SET_SLOT(to, Matrix_iSym, i1);
			if (cl[2] != 'C')
				SET_SLOT(to, Matrix_jSym, i1);
			UNPROTECT(1); /* i1 */
			if (cl[0] != 'n') {
				SEXP x1 = PROTECT(allocUnit(kindToType(cl[0]), n));
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
		}
		UNPROTECT(1); /* to */
		return to;
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (cl[0] != class[0]) {
		if (class[0] == 'n' && cl[0] == 'l')
			x0 = duplicateVector(x0);
		else
			x0 = Rf_coerceVector(x0, kindToType(cl[0]));
		if (class[0] == 'n')
			naToUnit(x0);
		UNPROTECT(1); /* x0 */
		PROTECT(x0);
	}

	int j, nnz = 0;
	if (cl[2] == 'T') {

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			for (j = 0; j < n; ++j) { \
				if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
					++nnz; \
				px0 += 1; \
			} \
		} while (0)

		SWITCH4(cl[0], TEMPLATE);

#undef TEMPLATE

	} else {

		SEXP p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) n + 1));
		int *pp1 = INTEGER(p1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			for (j = 0; j < n; ++j) { \
				if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
					++nnz; \
				px0 += 1; \
				*(pp1++) = nnz; \
			} \
		} while (0)

		SWITCH4(cl[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(1); /* p1 */

	}

	SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz));
	int *pi1 = INTEGER(i1);
	if (cl[2] != 'R')
		SET_SLOT(to, Matrix_iSym, i1);
	if (cl[2] != 'C')
		SET_SLOT(to, Matrix_jSym, i1);

	if (nnz == n) {
		for (j = 0; j < n; ++j)
			*(pi1++) = j;
		if (cl[0] != 'n')
			SET_SLOT(to, Matrix_xSym, x0);
	} else {

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			c##IF_NPATTERN( \
			SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
			c##TYPE *px1 = c##PTR(x1); \
			); \
			for (j = 0; j < n; ++j) { \
				if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) { \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = *px0; \
					); \
				} \
				px0 += 1; \
			} \
			c##IF_NPATTERN( \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(1); /* x1 */ \
			); \
		} while (0)

		SWITCH5(cl[0], TEMPLATE);

#undef TEMPLATE

	}

	UNPROTECT(3); /* i1, x0, to */
	return to;
}

/* as(<diagonalMatrix>, ".[gspt][CRT]Matrix") */
SEXP R_diagonal_as_sparse(SEXP s_from,
                          SEXP s_kind, SEXP s_shape, SEXP s_repr,
                          SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_diagonal, 0, __func__);

	char kind, shape, repr;
	VALID_KIND (s_kind , kind );
	VALID_SHAPE(s_shape, shape);
	VALID_REPR (s_repr , repr , 0);

	char ul = '\0', ct = '\0';
	if (shape != 'g') VALID_UPLO (s_uplo , ul);
	if (shape == 's') VALID_TRANS(s_trans, ct);

	return diagonal_as_sparse(s_from, class, kind, shape, repr, ul, ct);
}

SEXP index_as_sparse(SEXP from, const char *class, char kind, char repr)
{
	int mg = MARGIN(from);

	char cl[] = ".g.Matrix";
	cl[0] = (kind == '.') ? 'n' : ((kind == ',') ? 'd' : kind);
	cl[2] = (repr == '.') ? ((mg == 0) ? 'R' : 'C') : repr;
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1],
		r = (mg == 0) ? m : n, s = (mg == 0) ? n : m;
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	SEXP perm = PROTECT(GET_SLOT(from, Matrix_permSym));
	int *pperm = INTEGER(perm), k;

	if (cl[2] == 'T') {
		SEXP i = PROTECT(Rf_allocVector(INTSXP, r)),
			j = PROTECT(Rf_allocVector(INTSXP, r));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		for (k = 0; k < r; ++k) {
			*(pi++) = k;
			*(pj++) = *(pperm++) - 1;
		}
		SET_SLOT(to, Matrix_iSym, (mg == 0) ? i : j);
		SET_SLOT(to, Matrix_jSym, (mg == 0) ? j : i);
		UNPROTECT(2); /* j, i */
	} else if ((cl[2] == 'C') == (mg != 0)) {
		SEXP p = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) r + 1)),
			i = PROTECT(Rf_allocVector(INTSXP, r));
		int *pp = INTEGER(p), *pi = INTEGER(i);
		for (k = 0; k < r; ++k) {
			*(pp++) = k;
			*(pi++) = *(pperm++) - 1;
		}
		*pp = r;
		SET_SLOT(to, Matrix_pSym, p);
		SET_SLOT(to, (mg != 0) ? Matrix_iSym : Matrix_jSym, i);
		UNPROTECT(2); /* i, p */
	} else {
		SEXP p = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) s + 1));
		int *pp = INTEGER(p);
		memset(pp, 0, sizeof(int) * ((size_t) s + 1));
		for (k = 0; k < r; ++k)
			++pp[pperm[k]];
		for (k = 0; k < s; ++k)
			pp[k + 1] += pp[k];
		SEXP j = PROTECT(Rf_allocVector(INTSXP, r));
		int *pj = INTEGER(j), *work;
		Matrix_Calloc(work, s, int);
		memcpy(work, pp, sizeof(int) * (size_t) s);
		--work;
		for (k = 0; k < r; ++k)
			pj[work[pperm[k]]++] = k;
		++work;
		Matrix_Free(work, s);
		SET_SLOT(to, Matrix_pSym, p);
		SET_SLOT(to, (mg != 0) ? Matrix_jSym : Matrix_iSym, j);
		UNPROTECT(2); /* j, p */
	}

	if (cl[0] != 'n') {
		SEXP x = PROTECT(allocUnit(kindToType(cl[0]), r));
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(2); /* perm, to */
	return to;
}

/* as(<indMatrix>, ".g[CRT]Matrix") */
SEXP R_index_as_sparse(SEXP s_from, SEXP s_kind, SEXP s_repr)
{
	const char *class = Matrix_class(s_from, valid_index, 0, __func__);

	char kind, repr;
	VALID_KIND(s_kind, kind);
	VALID_REPR(s_repr, repr, 1);

	return index_as_sparse(s_from, class, kind, repr);
}

SEXP dense_as_kind(SEXP from, const char *class, char kind, int new)
{
	if (kind == '.')
		kind = class[0];
	else if (kind == ',')
		kind = (class[0] == 'z') ? 'z' : 'd';
	if (kind == class[0])
		return from;
	SEXPTYPE tt = kindToType(kind);

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = kind;
	cl[1] = class[1];
	cl[2] = class[2];
	if (cl[1] == 'p' && cl[0] != 'z' && cl[0] != 'd') {
		cl[1] = 's';
		if (!packed)
			cl[2] = 'y';
	}
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);

	PROTECT_INDEX pid;
	SEXP x;
	PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);

	if (TYPEOF(x) != tt) {
		REPROTECT(x = Rf_coerceVector(x, tt), pid);
		if (class[0] == 'n')
			/* n->[idz] */
			naToUnit(x);
	} else if (new) {
		REPROTECT(x = duplicateVector(x), pid);
		if (class[0] == 'n')
			/* n->l */
			naToUnit(x);
	} else {
		if (class[0] == 'n') {
			/* n->l */
			int *px = LOGICAL(x);
			for (R_xlen_t k = 0, kend = XLENGTH(x); k < kend; ++k) {
				if (*(px++) == NA_LOGICAL) {
					REPROTECT(x = duplicateVector(x), pid);
					naToUnit(x);
					break;
				}
			}
		}
	}

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(2); /* x, to */
	return to;
}

/* as(<denseMatrix>, "[nlidz]Matrix") */
SEXP R_dense_as_kind(SEXP s_from, SEXP s_kind)
{
	const char *class = Matrix_class(s_from, valid_dense, 2, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	return dense_as_kind(s_from, class, kind, 0);
}

SEXP sparse_as_kind(SEXP from, const char *class, char kind)
{
	if (kind == '.')
		kind = class[0];
	else if (kind == ',')
		kind = (class[0] == 'z') ? 'z' : 'd';
	if (kind == class[0])
		return from;
	SEXPTYPE tt = kindToType(kind);

	if (class[2] == 'T' && (class[0] == 'n' || class[0] == 'l') &&
	    kind != 'n' && kind != 'l') {
		/* defined in ./aggregate.c : */
		SEXP sparse_aggregate(SEXP, const char *);
		from = sparse_aggregate(from, class);
	}
	PROTECT(from);

	char cl[] = "...Matrix";
	cl[0] = kind;
	cl[1] = class[1];
	cl[2] = class[2];
	if (cl[1] == 'p' && cl[0] != 'z' && cl[0] != 'd')
		cl[1] = 's';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);

	R_xlen_t nnz = -1;
	if (class[2] != 'T') {
		SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
		if (nnz < 0)
			nnz = INTEGER(p)[XLENGTH(p) - 1];
		SET_SLOT(to, Matrix_pSym, p);
		UNPROTECT(1); /* p */
	}
	if (class[2] != 'R') {
		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
		if (nnz < 0)
			nnz = XLENGTH(i);
		SET_SLOT(to, Matrix_iSym, i);
		UNPROTECT(1); /* i */
	}
	if (class[2] != 'C') {
		SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
		if (nnz < 0)
			nnz = XLENGTH(j);
		SET_SLOT(to, Matrix_jSym, j);
		UNPROTECT(1); /* j */
	}
	if (class[0] == 'n') {
		SEXP x = PROTECT(allocUnit(tt, nnz));
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	} else if (kind != 'n') {
		PROTECT_INDEX pid;
		SEXP x;
		PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
		REPROTECT(x = Rf_coerceVector(x, tt), pid);
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(2); /* to, from */
	return to;
}

/* as(<denseMatrix>, "[nlidz]Matrix") */
SEXP R_sparse_as_kind(SEXP s_from, SEXP s_kind)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	return sparse_as_kind(s_from, class, kind);
}

SEXP diagonal_as_kind(SEXP from, const char *class, char kind)
{
	if (kind == '.')
		kind = class[0];
	else if (kind == ',')
		kind = (class[0] == 'z') ? 'z' : 'd';
	if (kind == class[0])
		return from;
	SEXPTYPE tt = kindToType(kind);

	char cl[] = ".diMatrix";
	cl[0] = kind;
	SEXP to = PROTECT(newObject(cl));

	int n = DIM(from)[1];
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	if (DIAG(from) != 'N')
		SET_DIAG(to);
	else {
		PROTECT_INDEX pid;
		SEXP x;
		PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
		if (TYPEOF(x) == tt) {
			if (class[0] == 'n') {
				/* n->l */
				int *px = LOGICAL(x);
				for (int j = 0; j < n; ++j) {
					if (*(px++) == NA_LOGICAL) {
						REPROTECT(x = duplicateVector(x), pid);
						naToUnit(x);
						break;
					}
				}
			}
		} else {
			REPROTECT(x = Rf_coerceVector(x, tt), pid);
			if (class[0] == 'n')
				/* n->[idz] */
				naToUnit(x);
		}
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<diagonalMatrix>, "[nlidz]Matrix") */
SEXP R_diagonal_as_kind(SEXP s_from, SEXP s_kind)
{
	const char *class = Matrix_class(s_from, valid_diagonal, 0, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	return diagonal_as_kind(s_from, class, kind);
}

SEXP index_as_kind(SEXP from, const char *class, char kind)
{
	return index_as_sparse(from, class, kind, '.');
}

/* as(<indMatrix>, "[nlidz]Matrix") */
SEXP R_index_as_kind(SEXP s_from, SEXP s_kind)
{
	const char *class = Matrix_class(s_from, valid_index, 0, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	return index_as_kind(s_from, class, kind);
}

SEXP dense_as_general(SEXP from, const char *class, int new)
{
	if (class[1] == 'g')
		return from;

	int packed = class[2] == 'p';

	char cl[] = ".geMatrix";
	cl[0] = class[0];
	SEXP to = PROTECT(newObject(cl));

	int n = DIM(from)[1];
	if ((int_fast64_t) n * n > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] == 's' || class[1] == 'p'), DIMNAMES(from, 0));

	char ul = UPLO(from), ct = '\0', nu = '\0';
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(from);
	if (class[1] == 't')
		nu = DIAG(from);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1 = x0;
	int nprotect = 2;
	if (packed || new) {
		PROTECT(x1 = Rf_allocVector(TYPEOF(x0), (R_xlen_t) n * n));
		++nprotect;
	}

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (packed) \
			c##NAME( pack1)(px1,  px0, (size_t) n, ul, ct, nu); \
		else if (new) \
			c##NAME(force2)(px1,  px0, (size_t) n, ul, ct, nu); \
		else \
			c##NAME(force2)(px1, NULL, (size_t) n, ul, ct, nu); \
	} while (0)

	SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(nprotect);
	return to;
}

/* as(<denseMatrix>, "generalMatrix") */
SEXP R_dense_as_general(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_dense, 2, __func__);
	return dense_as_general(s_from, class, 1);
}

SEXP sparse_as_general(SEXP from, const char *class)
{
	if (class[1] == 'g')
		return from;

	char cl[] = ".g.Matrix";
	cl[0] = class[0];
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

	int n = DIM(from)[1];
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] == 's' || class[1] == 'p'), DIMNAMES(from, 0));

	if (class[1] == 't' && DIAG(from) == 'N') {
		if (class[2] != 'T')
			COPY_SLOT(to, from, Matrix_pSym);
		if (class[2] != 'R')
			COPY_SLOT(to, from, Matrix_iSym);
		if (class[2] != 'C')
			COPY_SLOT(to, from, Matrix_jSym);
		if (class[0] != 'n')
			COPY_SLOT(to, from, Matrix_xSym);
		UNPROTECT(1); /* to */
		return to;
	}

	char ul = UPLO(from), ct = '\0';
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(from);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0)));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), *pp1 = INTEGER(p1),
			j, k, kend;
		pp0++; *(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

		if (class[1] == 's' || class[1] == 'p') {
			memset(pp1, 0, sizeof(int) * (size_t) n);
			for (j = 0, k = 0; j < n; ++j) {
				kend = pp0[j];
				while (k < kend) {
					if (pi0[k] != j)
						++pp1[pi0[k]];
					++k;
				}
			}
			for (j = 1; j < n; ++j)
				pp1[j] += pp1[j - 1];
			if (pp1[n - 1] > INT_MAX - pp0[n - 1])
				Rf_error(_("attempt to construct %s with more than %s nonzero entries"),
				         "sparseMatrix", "2^31-1");
			for (j = 0; j < n; ++j)
				pp1[j] += pp0[j];
		} else {
			if (n > INT_MAX - pp0[n - 1])
				Rf_error(_("attempt to construct %s with more than %s nonzero entries"),
				         "sparseMatrix", "2^31-1");
			for (j = 0; j < n; ++j)
				pp1[j] = pp0[j] + j + 1;
		}

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, pp1[n - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define TEMPLATE(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, pp1[n - 1])); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			); \
			if (class[1] == 's' || class[1] == 'p') { \
				int *pp1_; \
				Matrix_Calloc(pp1_, n, int); \
				memcpy(pp1_, pp1 - 1, sizeof(int) * (size_t) n); \
				if (ct == 'C') { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] == j) { \
							pi1[pp1_[j]] = pi0[k]; \
							c##IF_NPATTERN( \
							c##ASSIGN_PROJ_REAL(px1[pp1_[j]], px0[k]); \
							); \
							++pp1_[j]; \
						} else { \
							pi1[pp1_[j]] = pi0[k]; \
							c##IF_NPATTERN( \
							c##ASSIGN_IDEN(px1[pp1_[j]], px0[k]); \
							); \
							++pp1_[j]; \
							pi1[pp1_[pi0[k]]] = j; \
							c##IF_NPATTERN( \
							c##ASSIGN_CONJ(px1[pp1_[pi0[k]]], px0[k]); \
							); \
							++pp1_[pi0[k]]; \
						} \
						++k; \
					} \
				} \
				} else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] == j) { \
							pi1[pp1_[j]] = pi0[k]; \
							c##IF_NPATTERN( \
							c##ASSIGN_IDEN(px1[pp1_[j]], px0[k]); \
							); \
							++pp1_[j]; \
						} else { \
							pi1[pp1_[j]] = pi0[k]; \
							c##IF_NPATTERN( \
							c##ASSIGN_IDEN(px1[pp1_[j]], px0[k]); \
							); \
							++pp1_[j]; \
							pi1[pp1_[pi0[k]]] = j; \
							c##IF_NPATTERN( \
							c##ASSIGN_IDEN(px1[pp1_[pi0[k]]], px0[k]); \
							); \
							++pp1_[pi0[k]]; \
						} \
						++k; \
					} \
				} \
				} \
				Matrix_Free(pp1_, n); \
			} else if ((class[2] == 'C') == (ul == 'U')) { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						*(pi1++) = *(pi0++); \
						c##IF_NPATTERN( \
						*(px1++) = *(px0++); \
						); \
						++k; \
					} \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
				} \
			} else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					*(pi1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
					while (k < kend) { \
						*(pi1++) = *(pi0++); \
						c##IF_NPATTERN( \
						*(px1++) = *(px0++); \
						); \
						++k; \
					} \
				} \
			} \
			c##IF_NPATTERN( \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(4); /* i1, p1, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), j;
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1;

		if (class[1] == 's' || class[1] == 'p') {
			nnz1 = nnz0;
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					--nnz1;
		} else
			nnz1 = n;
		if (nnz1 > R_XLEN_T_MAX - nnz0)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		nnz1 += nnz0;

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
			); \
			if (class[1] == 's' || class[1] == 'p') { \
				if (ct == 'C') { \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 == *pj0) { \
						*(pi1++) = *pi0; \
						*(pj1++) = *pj0; \
						c##IF_NPATTERN( \
						c##ASSIGN_PROJ_REAL(*px1, *px0); \
						px1++; \
						); \
					} else { \
						*(pi1++) = *pi0; \
						*(pj1++) = *pj0; \
						*(pi1++) = *pj0; \
						*(pj1++) = *pi0; \
						c##IF_NPATTERN( \
						c##ASSIGN_IDEN(*px1, *px0); \
						px1++; \
						c##ASSIGN_CONJ(*px1, *px0); \
						px1++; \
						); \
					} \
					++pi0; ++pj0; c##IF_NPATTERN(++px0); \
				} \
				} else { \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 == *pj0) { \
						*(pi1++) = *pi0; \
						*(pj1++) = *pj0; \
						c##IF_NPATTERN( \
						c##ASSIGN_IDEN(*px1, *px0); \
						px1++; \
						); \
					} else { \
						*(pi1++) = *pi0; \
						*(pj1++) = *pj0; \
						*(pi1++) = *pj0; \
						*(pj1++) = *pi0; \
						c##IF_NPATTERN( \
						c##ASSIGN_IDEN(*px1, *px0); \
						px1++; \
						c##ASSIGN_IDEN(*px1, *px0); \
						px1++; \
						); \
					} \
					++pi0; ++pj0; c##IF_NPATTERN(++px0); \
				} \
				} \
			} else { \
				memcpy(pi1, pi0, sizeof(    int) * (size_t) nnz0); \
				memcpy(pj1, pj0, sizeof(    int) * (size_t) nnz0); \
				pi1 += nnz0; \
				pj1 += nnz0; \
				c##IF_NPATTERN( \
				memcpy(px1, px0, sizeof(c##TYPE) * (size_t) nnz0); \
				px1 += nnz0; \
				); \
				for (j = 0; j < n; ++j) { \
					*(pi1++) = *(pj1++) = j; \
					c##IF_NPATTERN( \
					*(px1++) = c##UNIT; \
					); \
				} \
			} \
			c##IF_NPATTERN( \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		UNPROTECT(4); /* j1, i1, j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "generalMatrix") */
SEXP R_sparse_as_general(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_as_general(s_from, class);
}

SEXP dense_as_unpacked(SEXP from, const char *class)
{
	if (class[2] != 'p')
		return from;

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = (class[1] == 's') ? 'y' : ((class[1] == 'p') ? 'o' : 'r');
	SEXP to = PROTECT(newObject(cl));

	int n = DIM(from)[0];
	if ((int_fast64_t) n * n > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	char ul = UPLO(from);
	if (ul != 'U')
		SET_UPLO(to);
	if (cl[1] == 's' && cl[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (cl[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (cl[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);
	if (cl[1] == 'o')
		COPY_SLOT(to, from, Matrix_sdSym);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), (R_xlen_t) n * n));

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		c##NAME(pack1)(px1, px0, (size_t) n, ul, '\0', 'N'); \
	} while (0)

	SWITCH4((cl[0] == 'c') ? 'd' : cl[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<denseMatrix>, "unpackedMatrix") */
SEXP R_dense_as_unpacked(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);
	return dense_as_unpacked(s_from, class);
}

SEXP dense_as_packed(SEXP from, const char *class, char ul, char ct, char nu)
{
	if (class[2] == 'p')
		return from;

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = 'p';
	if (class[1] == 'g')
		cl[1] = (nu == '\0') ? 's' : 't';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("attempt to pack non-square matrix"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));

	if (class[1] == 'g') {
		if (ul == '\0')
			ul = 'U';
		if (ul != 'U')
			SET_UPLO(to);
		if (ct == '\0')
			ct = 'C';
		if (cl[1] == 's' && ct != 'C' && cl[0] == 'z')
			SET_TRANS(to);
		if (cl[1] == 't' && nu != 'N')
			SET_DIAG(to);
	} else {
		if ((ul = UPLO(from)) != 'U')
			SET_UPLO(to);
		if (cl[1] == 's' && cl[0] == 'z' && (ct = TRANS(from)) != 'C')
			SET_TRANS(to);
		if (cl[1] == 't' && (nu = DIAG(from)) != 'N')
			SET_DIAG(to);
		if (cl[1] != 't')
			COPY_SLOT(to, from, Matrix_factorsSym);
		if (cl[1] == 'o')
			COPY_SLOT(to, from, Matrix_sdSym);
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), PACKED_LENGTH((R_xlen_t) n)));

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		c##NAME(pack2)(px1, px0, (size_t) n, ul, '\0', 'N'); \
	} while (0)

	SWITCH4((cl[0] == 'c') ? 'd' : cl[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<denseMatrix>, "packedMatrix") */
SEXP R_dense_as_packed(SEXP s_from, SEXP s_uplo, SEXP s_trans, SEXP s_diag)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] == 'g') {
	if (s_uplo  != R_NilValue) VALID_UPLO (s_uplo , ul);
	if (s_trans != R_NilValue) VALID_TRANS(s_trans, ct);
	if (s_diag  != R_NilValue) VALID_DIAG (s_diag , nu);
	}

	return dense_as_packed(s_from, class, ul, ct, nu);
}

SEXP sparse_as_Csparse(SEXP from, const char *class)
{
	if (class[2] == 'C')
		return from;

	char cl[] = "..CMatrix";
	cl[0] = class[0];
	cl[1] = class[1];
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);

	if (class[2] == 'R') {
		SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) n + 1)),
			i1 = PROTECT(Rf_allocVector(INTSXP, INTEGER(p0)[m]));
		int *pp0 = INTEGER(p0), *pj0 = INTEGER(j0),
			*pp1 = INTEGER(p1), *pi1 = INTEGER(i1), *iwork = NULL;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to, Matrix_iSym, i1);
		Matrix_Calloc(iwork, n, int);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, INTEGER(p0)[m])); \
			px0 = c##PTR(x0); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			c##csptrans(pp1, pi1, px1, pp0, pj0, px0, n, m, 'T', iwork); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		Matrix_Free(iwork, n);
		UNPROTECT(4); /* i1, p1, j0, p0 */
	} else {
		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym));
		if (XLENGTH(i0) > INT_MAX)
			Rf_error(_("number of triplets to be aggregated exceeds %s"),
			         "2^31-1");
		int *pi0 = INTEGER(i0), nnz = (int) XLENGTH(i0);

		SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym)), p1, i1;
		int *pj0 = INTEGER(j0), *pp1 = NULL, *pi1 = NULL, *iwork = NULL;
		size_t liwork = (size_t) ((int_fast64_t) m + 1 + m + n + nnz),
			lwork = (size_t) nnz;
		Matrix_Calloc(iwork, liwork, int);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL, *work = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
			px0 = c##PTR(x0); \
			Matrix_Calloc(work, lwork, c##TYPE); \
			); \
			c##tspsort(pp1, pi1, px1, pi0, pj0, px0, m, n, &nnz, iwork, work); \
			PROTECT(p1 = Rf_allocVector(INTSXP, (R_xlen_t) n + 1)), \
			PROTECT(i1 = Rf_allocVector(INTSXP, nnz)); \
			pp1 = INTEGER(p1); \
			pi1 = INTEGER(i1); \
			SET_SLOT(to, Matrix_pSym, p1); \
			SET_SLOT(to, Matrix_iSym, i1); \
			c##IF_NPATTERN( \
			SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			); \
			c##tspsort(pp1, pi1, px1, pi0, pj0, px0, m, n, &nnz, iwork, work); \
			c##IF_NPATTERN( \
			UNPROTECT(1); /* x1 */ \
			); \
			UNPROTECT(2); /* i1, p1 */ \
			c##IF_NPATTERN( \
			Matrix_Free(work, lwork); \
			UNPROTECT(1); /* x0 */ \
			); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		Matrix_Free(iwork, liwork);
		UNPROTECT(2); /* j0, i0 */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "CsparseMatrix") */
SEXP R_sparse_as_Csparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_as_Csparse(s_from, class);
}

SEXP sparse_as_Rsparse(SEXP from, const char *class)
{
	if (class[2] == 'R')
		return from;

	char cl[] = "..RMatrix";
	cl[0] = class[0];
	cl[1] = class[1];
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);

	if (class[2] == 'C') {
		SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) m + 1)),
			j1 = PROTECT(Rf_allocVector(INTSXP, INTEGER(p0)[n]));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			*pp1 = INTEGER(p1), *pj1 = INTEGER(j1), *iwork = NULL;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to, Matrix_jSym, j1);
		Matrix_Calloc(iwork, m, int);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, INTEGER(p0)[n])); \
			px0 = c##PTR(x0); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			c##csptrans(pp1, pj1, px1, pp0, pi0, px0, m, n, 'T', iwork); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		Matrix_Free(iwork, m);
		UNPROTECT(4); /* j1, p1, i0, p0 */
	} else {
		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym));
		if (XLENGTH(i0) > INT_MAX)
			Rf_error(_("number of triplets to be aggregated exceeds %s"),
			         "2^31-1");
		int *pi0 = INTEGER(i0), nnz = (int) XLENGTH(i0);

		SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym)), p1, j1;
		int *pj0 = INTEGER(j0), *pp1 = NULL, *pj1 = NULL, *iwork = NULL;
		size_t liwork = (size_t) ((int_fast64_t) n + 1 + n + m + nnz),
			lwork = (size_t) nnz;
		Matrix_Calloc(iwork, liwork, int);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL, *work = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
			px0 = c##PTR(x0); \
			Matrix_Calloc(work, nnz, c##TYPE); \
			); \
			c##tspsort(pp1, pj1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
			PROTECT(p1 = Rf_allocVector(INTSXP, (R_xlen_t) m + 1)), \
			PROTECT(j1 = Rf_allocVector(INTSXP, nnz)); \
			pp1 = INTEGER(p1); \
			pj1 = INTEGER(j1); \
			SET_SLOT(to, Matrix_pSym, p1); \
			SET_SLOT(to, Matrix_jSym, j1); \
			c##IF_NPATTERN( \
			SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			); \
			c##tspsort(pp1, pj1, px1, pj0, pi0, px0, n, m, &nnz, iwork, work); \
			c##IF_NPATTERN( \
			UNPROTECT(1); /* x1 */ \
			); \
			UNPROTECT(2); /* j1, p1 */ \
			c##IF_NPATTERN( \
			Matrix_Free(work, lwork); \
			UNPROTECT(1); /* x0 */ \
			); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		Matrix_Free(iwork, liwork);
		UNPROTECT(2); /* j0, i0 */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "RsparseMatrix") */
SEXP R_sparse_as_Rsparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_as_Rsparse(s_from, class);
}

SEXP sparse_as_Tsparse(SEXP from, const char *class)
{
	if (class[2] == 'T')
		return from;

	char cl[] = "..TMatrix";
	cl[0] = class[0];
	cl[1] = class[1];
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, m, n);
	SET_DIMNAMES(to, 0, DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) != 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] != 't')
		COPY_SLOT(to, from, Matrix_factorsSym);

	if (class[2] == 'R')
		SWAP(m, n, int, );

	SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
	int *pp0 = INTEGER(p0), nnz = pp0[n];
	pp0++;

	SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		i0 = PROTECT(GET_SLOT(from, iSym));
	if (XLENGTH(i0) == nnz)
		SET_SLOT(to, iSym, i0);
	else {
		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz));
		memcpy(INTEGER(i1), INTEGER(i0), sizeof(int) * (size_t) nnz);
		SET_SLOT(to, iSym, i1);
		UNPROTECT(1); /* i1 */
	}

	SEXP jSym = (class[2] == 'C') ? Matrix_jSym : Matrix_iSym,
		j1 = PROTECT(Rf_allocVector(INTSXP, nnz));
	int *pj1 = INTEGER(j1), j, k, kend;
	SET_SLOT(to, jSym, j1);
	for (j = 0, k = 0; j < n; ++j) {
		kend = pp0[j];
		while (k < kend)
			pj1[k++] = j;
	}

	if (class[0] != 'n') {
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		if (XLENGTH(x0) == nnz)
			SET_SLOT(to, Matrix_xSym, x0);
		else {
			SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz));

#define TEMPLATE(c) \
			memcpy(c##PTR(x1), c##PTR(x0), sizeof(c##TYPE) * (size_t) nnz)

			SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */
		}
		UNPROTECT(1); /* x0 */
	}

	UNPROTECT(4); /* j1, i0, p0, to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "TsparseMatrix") */
SEXP R_sparse_as_Tsparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_as_Tsparse(s_from, class);
}

SEXP vector_as_Vector(SEXP from, char kind)
{
	R_xlen_t vlen = XLENGTH(from);
	if (vlen > 0x1.0p+53)
		Rf_error(_("%s length cannot exceed %s"), "sparseVector", "2^53");

	SEXPTYPE tf = TYPEOF(from);
	char cl[] = ".sparseVector";
	cl[0] = (kind == '.') ? typeToKind(tf) : ((kind == ',') ? ((tf == CPLXSXP) ? 'z' : 'd') : kind);
	SEXPTYPE tt = kindToType(cl[0]);
	SEXP to = newObject(cl);
	if (vlen == 0)
		return to;
	PROTECT(to);

	SEXP length = PROTECT(Rf_allocVector((vlen <= INT_MAX) ? INTSXP : REALSXP, 1));
	if (TYPEOF(length) == INTSXP)
	INTEGER(length)[0] = (int) vlen;
	else
	REAL(length)[0] = (double) vlen;
	SET_SLOT(to, Matrix_lengthSym, length);
	UNPROTECT(1); /* length */

	SEXP x0 = from;
	R_xlen_t k, nnz = 0;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0); \
		for (k = 0; k < vlen; ++k) \
			if (c##NOT_ZERO(px0[k])) \
				++nnz; \
	} while (0)

	SWITCH4(typeToKind(tf), TEMPLATE);

#undef TEMPLATE

	SEXP i1 = PROTECT(Rf_allocVector(TYPEOF(length), nnz));
	SET_SLOT(to, Matrix_iSym, i1);

#define TEMPLATE(c) \
	do { \
		if (TYPEOF(i1) == INTSXP) \
			TEMPLATE__(i, c); \
		else \
			TEMPLATE__(d, c); \
	} while (0)

#define TEMPLATE__(d, c) \
	do { \
		d##TYPE *pi1 = d##PTR(i1); \
		c##TYPE *px0 = c##PTR(x0); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
		c##TYPE *px1 = c##PTR(x1); \
		); \
		for (k = 0; k < vlen; ++k) \
			if (c##NOT_ZERO(px0[k])) { \
				*(pi1++) = (d##TYPE) (k + 1); \
				c##IF_NPATTERN( \
				*(px1++) = px0[k]; \
				); \
			} \
		c##IF_NPATTERN( \
		PROTECT(x1 = Rf_coerceVector(x1, tt)); \
		SET_SLOT(to, Matrix_xSym, x1); \
		UNPROTECT(2); /* x1 */ \
		); \
	} while (0)

	SWITCH5((kind == 'n') ? 'n' : typeToKind(tf), TEMPLATE);

#undef TEMPLATE__
#undef TEMPLATE

	UNPROTECT(2); /* i1, to */
	return to;
}

SEXP R_vector_as_Vector(SEXP s_from, SEXP s_kind)
{
	switch (TYPEOF(s_from)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
		break;
	default:
		ERROR_INVALID_TYPE(s_from, __func__);
		break;
	}

	char kind;
	VALID_KIND(s_kind, kind);

	return vector_as_Vector(s_from, kind);
}

SEXP sparse_as_Vector(SEXP from, const char *class)
{
	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;
	if (mn > 0x1.0p+53)
		Rf_error(_("%s length cannot exceed %s"), "sparseVector", "2^53");

	char cl[] = ".sparseVector";
	cl[0] = class[0];
	SEXP to = newObject(cl);
	if (mn == 0)
		return to;
	PROTECT(to);

	SEXP length = PROTECT(Rf_allocVector((mn <= INT_MAX) ? INTSXP : REALSXP, 1));
	if (TYPEOF(length) == INTSXP)
	INTEGER(length)[0] = (int) mn;
	else
	REAL(length)[0] = (double) mn;
	SET_SLOT(to, Matrix_lengthSym, length);
	UNPROTECT(1); /* length */

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);
	if (class[2] == 'T') {
	REPROTECT(from = sparse_as_Csparse(from, class), pid);
	class = Matrix_class(from, valid_sparse, 2, __func__);
	}
	REPROTECT(from = sparse_as_general(from, class), pid);

	SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
	int *pp0 = INTEGER(p0), nnz = pp0[(class[2] == 'C') ? n : m];
	pp0++;
	if (nnz == 0) {
		UNPROTECT(3); /* p0, from, to */
		return to;
	}

	SEXP i1 = PROTECT(Rf_allocVector(TYPEOF(length), nnz));
	SET_SLOT(to, Matrix_iSym, i1);

#define TEMPLATE(c) \
	do { \
		if (TYPEOF(i1) == INTSXP) \
			TEMPLATE__(i, c); \
		else \
			TEMPLATE__(d, c); \
	} while (0)

	if (class[2] == 'C') {
		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym));
		int *pi0 = INTEGER(i0), j, k, kend;

#define TEMPLATE__(d, c) \
		do { \
			d##TYPE *pi1 = d##PTR(i1), l = (d##TYPE) 1, dl = (d##TYPE) m; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			for (j = 0, k = 0; j < n; ++j) { \
				kend = pp0[j]; \
				while (k < kend) { \
					c##IF_NPATTERN( \
					*(px1++) = px0[k]; \
					); \
					*(pi1++) = pi0[k] + l; \
					++k; \
				} \
				l += dl; \
			} \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE__

		UNPROTECT(1); /* i0 */
	} else {
		SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pj0 = INTEGER(j0), i, j, k, kend, *iwork = NULL;
		size_t liwork = (size_t) n + 1;
		Matrix_Calloc(iwork, liwork, int);

		++iwork;
		for (k = 0; k < nnz; ++k)
			++iwork[pj0[k]];
		for (j = 0; j < n; ++j)
			iwork[j] += iwork[j - 1];
		--iwork;

#define TEMPLATE__(d, c) \
		do { \
			d##TYPE *pi1 = d##PTR(i1); \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			for (i = 0, k = 0; i < m; ++i) { \
				kend = pp0[i]; \
				while (k < kend) { \
					j = pj0[k]; \
					c##IF_NPATTERN( \
					px1[iwork[j]  ] = px0[k]; \
					); \
					pi1[iwork[j]++] = i + 1 + j * (d##TYPE) m; \
					++k; \
				} \
			} \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE__

		Matrix_Free(iwork, liwork);
		UNPROTECT(1); /* j0 */
	}

#undef TEMPLATE

	UNPROTECT(4); /* i1, p0, from, to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "sparseVector") */
SEXP R_sparse_as_Vector(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_as_Vector(s_from, class);
}

SEXP diagonal_as_Vector(SEXP from, const char *class)
{
	int n = DIM(from)[0];
	int_fast64_t nn = (int_fast64_t) n * n;
	if (nn > 0x1.0p+53)
		Rf_error(_("%s length cannot exceed %s"), "sparseVector", "2^53");

	char cl[] = ".sparseVector";
	cl[0] = class[0];
	SEXP to = newObject(cl);
	if (nn == 0)
		return to;
	PROTECT(to);

	SEXP length = PROTECT(Rf_allocVector((nn <= INT_MAX) ? INTSXP : REALSXP, 1));
	if (TYPEOF(length) == INTSXP)
	INTEGER(length)[0] = (int) nn;
	else
	REAL(length)[0] = (double) nn;
	SET_SLOT(to, Matrix_lengthSym, length);
	UNPROTECT(1); /* length */

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	int j, nnz1;
	char nu = DIAG(from);

	if (nu != 'N')
		nnz1 = n;
	else {
		nnz1 = 0;

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			for (j = 0; j < n; ++j) \
				if (c##NOT_ZERO(px0[j])) \
					++nnz1; \
			} while (0)

		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

	}

	SEXP i1 = PROTECT(Rf_allocVector(TYPEOF(length), nnz1));
	SET_SLOT(to, Matrix_iSym, i1);

#define TEMPLATE(c) \
	do { \
		if (TYPEOF(i1) == INTSXP) \
			TEMPLATE__(i, c); \
		else \
			TEMPLATE__(d, c); \
	} while (0)

#define TEMPLATE__(d, c) \
	do { \
		d##TYPE *pi1 = d##PTR(i1), l = (d##TYPE) 1, dl = (d##TYPE) (n + 1); \
		c##TYPE *px0 = c##PTR(x0); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz1)); \
		c##TYPE *px1 = c##PTR(x1); \
		SET_SLOT(to, Matrix_xSym, x1); \
		UNPROTECT(1); /* x1 */ \
		); \
		for (j = 0; j < n; ++j) { \
			if (nu != 'N' || c##NOT_ZERO(px0[j])) { \
				*(pi1++) = l; \
				c##IF_NPATTERN( \
				*(px1++) = (nu != 'N') ? c##UNIT : px0[j]; \
				); \
			} \
			l += dl; \
		} \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE__
#undef TEMPLATE

	UNPROTECT(3); /* i1, x0, to */
	return to;
}

/* as(<diagonalMatrix>, "sparseVector") */
SEXP R_diagonal_as_Vector(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_diagonal, 0, __func__);
	return diagonal_as_Vector(s_from, class);
}

SEXP index_as_Vector(SEXP from, const char *class)
{
	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;
	if (mn > 0x1.0p+53)
		Rf_error(_("%s length cannot exceed %s"), "sparseVector", "2^53");

	SEXP to = newObject("nsparseVector");
	if (mn == 0)
		return to;
	PROTECT(to);

	SEXP length = PROTECT(Rf_allocVector((mn <= INT_MAX) ? INTSXP : REALSXP, 1));
	if (TYPEOF(length) == INTSXP)
	INTEGER(length)[0] = (int) mn;
	else
	REAL(length)[0] = (double) mn;
	SET_SLOT(to, Matrix_lengthSym, length);
	UNPROTECT(1); /* length */

	SEXP perm = PROTECT(GET_SLOT(from, Matrix_permSym));
	int *pperm = INTEGER(perm), mg = MARGIN(from);

	SEXP i1 = PROTECT(Rf_allocVector(TYPEOF(length), (mg == 0) ? m : n));
	SET_SLOT(to, Matrix_iSym, i1);

	if (mg == 0) {
		int i, j, *iwork = NULL;
		size_t liwork = (size_t) n + 1;
		Matrix_Calloc(iwork, liwork, int);
		for (i = 0; i < m; ++i)
			++iwork[pperm[i]];
		for (j = 1; j <= n; ++j)
			iwork[j] += iwork[j - 1];
		if (TYPEOF(i1) == INTSXP) {
			int *pi1 = INTEGER(i1);
			for (i = 0; i < m; ++i) {
				j = pperm[i] - 1;
				pi1[iwork[j]++] = i + 1 + j * m;
			}
		} else {
			double *pi1 = REAL(i1);
			for (i = 0; i < m; ++i) {
				j = pperm[i] - 1;
				pi1[iwork[j]++] = i + 1 + j * (double) m;
			}
		}
		Matrix_Free(iwork, liwork);
	} else {
		int j;
		if (TYPEOF(i1) == INTSXP) {
			int *pi1 = INTEGER(i1);
			for (j = 0; j < n; ++j)
				pi1[j] = pperm[j] + j * m;
		} else {
			double *pi1 = REAL(i1);
			for (j = 0; j < n; ++j)
				pi1[j] = pperm[j] + j * (double) m;
		}
	}

	UNPROTECT(3); /* i1, perm, to */
	return to;
}

/* as(<indMatrix>, "sparseVector") */
SEXP R_index_as_Vector(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_index, 0, __func__);
	return index_as_Vector(s_from, class);
}

/* as(<Matrix>, "vector") */
SEXP R_Matrix_as_vector(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 7, __func__);

	SEXP to = R_NilValue;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(s_from, &pid);

	switch (class[2]) {
	case 'e':
		to = GET_SLOT(s_from, Matrix_xSym);
		if (class[0] == 'n') {
			int *px = LOGICAL(to);
			for (R_xlen_t k = 0, kend = XLENGTH(to); k < kend; ++k) {
				if (*(px++) == NA_LOGICAL) {
					PROTECT(to);
					to = duplicateVector(to);
					UNPROTECT(1); /* to */
					break;
				}
			}
		}
		break;
	case 'y':
	case 'r':
	case 'p':
		REPROTECT(s_from = dense_as_general(s_from, class, 1), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	case 'C':
	case 'R':
	case 'T':
		REPROTECT(s_from = sparse_as_dense(s_from, class, 0), pid);
		REPROTECT(s_from = dense_as_general(s_from, class, 0), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	case 'i':
		REPROTECT(s_from = diagonal_as_dense(s_from, class, '.', 'g', 0, '\0', '\0'), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	case 'd':
		REPROTECT(s_from = index_as_dense(s_from, class, 'n'), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	default:
		break;
	}

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
	case 'i':
		if (class[0] == 'n') {
			PROTECT(to);
			naToUnit(to);
			UNPROTECT(1); /* to */
		}
		break;
	default:
		break;
	}

	UNPROTECT(1); /* s_from */
	return to;
}

/* as(<Matrix>, "matrix") */
SEXP R_Matrix_as_matrix(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 7, __func__);

	SEXP to = R_NilValue;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(s_from, &pid);

	switch (class[2]) {
	case 'e':
		PROTECT(to = GET_SLOT(s_from, Matrix_xSym));
		to = duplicateVector(to);
		UNPROTECT(1); /* to */
		break;
	case 'y':
	case 'r':
	case 'p':
		REPROTECT(s_from = dense_as_general(s_from, class, 1), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	case 'C':
	case 'R':
	case 'T':
		REPROTECT(s_from = sparse_as_dense(s_from, class, 0), pid);
		REPROTECT(s_from = dense_as_general(s_from, class, 0), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	case 'i':
		REPROTECT(s_from = diagonal_as_dense(s_from, class, '.', 'g', 0, '\0', '\0'), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	case 'd':
		REPROTECT(s_from = index_as_dense(s_from, class, 'n'), pid);
		to = GET_SLOT(s_from, Matrix_xSym);
		break;
	default:
		break;
	}
	PROTECT(to);

	SEXP dim = PROTECT(GET_SLOT(s_from, Matrix_DimSym));
	Rf_setAttrib(to, R_DimSymbol, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(DIMNAMES(s_from, 0));
	if (!DimNames_is_trivial(dimnames))
		Rf_setAttrib(to, R_DimNamesSymbol, dimnames);
	UNPROTECT(1); /* dimnames */

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
	case 'i':
		if (class[0] == 'n')
			naToUnit(to);
		break;
	default:
		break;
	}

	UNPROTECT(2); /* to, s_from */
	return to;
}

/* as(<Matrix>, "unpackedMatrix") */
SEXP R_Matrix_as_unpacked(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 1, __func__);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'o':
	case 'r':
		return s_from;
	case 'p':
		return dense_as_unpacked(s_from, class);
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_dense(s_from, class, 0);
	case 'i':
		return diagonal_as_dense(s_from, class, '.', 't', 0, 'U', '\0');
	case 'd':
		return index_as_dense(s_from, class, 'n');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "packedMatrix") */
SEXP R_Matrix_as_packed(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 1, __func__);

	switch (class[2]) {
	case 'e':
		Rf_error(_("attempt to pack a %s"), "generalMatrix");
		return R_NilValue;
	case 'y':
	case 'o':
	case 'r':
		return dense_as_packed(s_from, class, '\0', '\0', '\0');
	case 'p':
		return s_from;
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_dense(s_from, class, 1);
	case 'i':
		return diagonal_as_dense(s_from, class, '.', 't', 1, 'U', '\0');
	case 'd':
		Rf_error(_("attempt to pack an %s"), "indMatrix");
		return R_NilValue;
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "CsparseMatrix") */
SEXP R_Matrix_as_Csparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 3, __func__);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'o':
	case 'r':
	case 'p':
		return dense_as_sparse(s_from, class, 'C');
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Csparse(s_from, class);
	case 'i':
		return diagonal_as_sparse(s_from, class, '.', 't', 'C', 'U', '\0');
	case 'd':
		return index_as_sparse(s_from, class, 'n', 'C');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "RsparseMatrix") */
SEXP R_Matrix_as_Rsparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 3, __func__);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'o':
	case 'r':
	case 'p':
		return dense_as_sparse(s_from, class, 'R');
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Rsparse(s_from, class);
	case 'i':
		return diagonal_as_sparse(s_from, class, '.', 't', 'R', 'U', '\0');
	case 'd':
		return index_as_sparse(s_from, class, 'n', 'R');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "TsparseMatrix") */
SEXP R_Matrix_as_Tsparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 3, __func__);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'o':
	case 'r':
	case 'p':
		return dense_as_sparse(s_from, class, 'T');
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Tsparse(s_from, class);
	case 'i':
		return diagonal_as_sparse(s_from, class, '.', 't', 'T', 'U', '\0');
	case 'd':
		return index_as_sparse(s_from, class, 'n', 'T');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "sparseVector") */
SEXP R_Matrix_as_Vector(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 7, __func__);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		PROTECT(s_from = dense_as_sparse(s_from, class, 'C'));
		class = Matrix_class(s_from, valid_matrix, 7, __func__);
		s_from = sparse_as_Vector(s_from, class);
		UNPROTECT(1); /* s_from */
		return s_from;
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Vector(s_from, class);
	case 'i':
		return diagonal_as_Vector(s_from, class);
	case 'd':
		return index_as_Vector(s_from, class);
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "[nldiz]Matrix") */
SEXP R_Matrix_as_kind(SEXP s_from, SEXP s_kind, SEXP s_sparse)
{
	const char *class = Matrix_class(s_from, valid_matrix, 3, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	int sparse;
	VALID_LOGIC3(s_sparse, sparse);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'o':
	case 'r':
	case 'p':
		if (sparse == NA_LOGICAL || sparse == 0)
			s_from = dense_as_kind(s_from, class, kind, 0);
		else {
			s_from = dense_as_sparse(s_from, class, 'C');
			PROTECT(s_from);
			class = Matrix_class(s_from, valid_matrix, 3, __func__);
			s_from = sparse_as_kind(s_from, class, kind);
			UNPROTECT(1); /* s_from */
		}
		return s_from;
	case 'C':
	case 'R':
	case 'T':
		s_from = sparse_as_kind(s_from, class, kind);
		if (sparse == 0) {
			PROTECT(s_from);
			class = Matrix_class(s_from, valid_matrix, 3, __func__);
			s_from = sparse_as_dense(s_from, class, 0);
			UNPROTECT(1); /* s_from */
		}
		return s_from;
	case 'i':
		if (sparse == NA_LOGICAL)
			s_from = diagonal_as_kind(s_from, class, kind);
		else if (sparse != 0)
			s_from = diagonal_as_sparse(s_from, class, kind, 't', 'C', 'U', '\0');
		else
			s_from = diagonal_as_dense(s_from, class, kind, 't', 0, 'U', '\0');
		return s_from;
	case 'd':
		if (sparse == NA_LOGICAL || sparse != 0)
			s_from = index_as_sparse(s_from, class, kind, '.');
		else
			s_from = index_as_dense(s_from, class, kind);
		return s_from;
	default:
		return R_NilValue;
	}
}

/* as(as(<Matrix>, "[nlidz]Matrix"), "generalMatrix") */
SEXP R_Matrix_as_general(SEXP s_from, SEXP s_kind)
{
	const char *class = Matrix_class(s_from, valid_matrix, 3, __func__);

	char kind;
	VALID_KIND(s_kind, kind);

	switch (class[2]) {
	case 'e':
	case 'y':
	case 'o':
	case 'r':
	case 'p':
	{
		char z = class[0];
		s_from = dense_as_kind(s_from, class, kind, 1);
		PROTECT(s_from);
		class = Matrix_class(s_from, valid_matrix, 3, __func__);
		s_from = dense_as_general(s_from, class, kindToType(class[0]) == kindToType(z));
		UNPROTECT(1); /* s_from */
		return s_from;
	}
	case 'C':
	case 'R':
	case 'T':
		s_from = sparse_as_kind(s_from, class, kind);
		PROTECT(s_from);
		class = Matrix_class(s_from, valid_matrix, 3, __func__);
		s_from = sparse_as_general(s_from, class);
		UNPROTECT(1); /* s_from */
		return s_from;
	case 'i':
		return diagonal_as_sparse(s_from, class, kind, 'g', 'C', '\0', '\0');
	case 'd':
		return index_as_sparse(s_from, class, kind, '.');
	default:
		return R_NilValue;
	}
}
