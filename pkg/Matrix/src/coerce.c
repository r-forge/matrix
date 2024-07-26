#include <math.h> /* trunc */
#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "coerce.h"

SEXP vector_as_dense(SEXP from, const char *zzz,
                     char ul, char ct, char di,
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
	PROTECT(from = coerceVector(from, tt));

	if (cl[1] != 'g' && m != n)
		error(_("attempt to construct non-square %s"),
		      (cl[1] == 's' || cl[1] == 'p') ? "symmetricMatrix" : "triangularMatrix");

	int_fast64_t mn = (int_fast64_t) m * n,
		xlen = (!packed) ? mn : n + (mn - n) / 2;
	if (xlen > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP to = PROTECT(newObject(cl));

	SEXP dim = GET_SLOT(to, Matrix_DimSym);
	int *pdim = INTEGER(dim);
	pdim[0] = m;
	pdim[1] = n;

	if (cl[1] == 's' || cl[1] == 'p')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (cl[1] == 't' && di != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	/* FIXME: add argument 'new' and conditionally avoid allocation */
	SEXP x = PROTECT(allocVector(tt, (R_xlen_t) xlen));
	R_xlen_t k, r = XLENGTH(from);
	int i, j, recycle = r < mn;

#define VAD(c) \
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
						if (k == r) k = 0; \
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

	SWITCH4(cl[0], VAD);

#undef VAD

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(3); /* x, to, from */
	return to;
}

SEXP R_vector_as_dense(SEXP s_from, SEXP s_zzz,
                       SEXP s_uplo, SEXP s_trans, SEXP s_diag,
                       SEXP s_m, SEXP s_n, SEXP s_byrow, SEXP s_dimnames)
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
		error(_("second argument of '%s' does not specify a subclass of %s"),
		      __func__, "denseMatrix");

	char ul = '\0', ct = '\0', di = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , di);

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
					error(_("dimensions cannot exceed %s"), "2^31-1");
				m = (int) tmp;
			}
		}
		if (m < 0)
			error(_("invalid '%s' to '%s'"), "m", __func__);
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
					error(_("dimensions cannot exceed %s"), "2^31-1");
				n = (int) tmp;
			}
		}
		if (n < 0)
			error(_("invalid '%s' to '%s'"), "n", __func__);
	}

	int byrow;
	VALID_LOGIC2(s_byrow, byrow);

	if (s_dimnames != R_NilValue)
		if (TYPEOF(s_dimnames) != VECSXP || LENGTH(s_dimnames) != 2)
			error(_("invalid '%s' to '%s'"), "dimnames", __func__);

	R_xlen_t vlen = XLENGTH(s_from);
	if (zzz[1] != 'g' && (m < 0) != (n < 0)) {
		if (m < 0)
			m = n;
		else
			n = m;
	} else if (m < 0 && n < 0) {
		if (vlen > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
		m = (int) vlen;
		n = 1;
	} else if (m < 0) {
		if (vlen > (int_fast64_t) INT_MAX * n) {
			if (n == 0)
				error(_("nonempty vector supplied for empty matrix"));
			else
				error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		m = (int) ((n == 0) ? 0 : vlen / n + (vlen % n != 0));
	} else if (n < 0) {
		if (vlen > (int_fast64_t) m * INT_MAX) {
			if (m == 0)
				error(_("nonempty vector supplied for empty matrix"));
			else
				error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		n = (int) ((m == 0) ? 0 : vlen / m + (vlen % m != 0));
	}

	int_fast64_t mlen = (int_fast64_t) m * n;
	if (vlen <= 1)
		/* do nothing */ ;
	else if (mlen == 0)
		warning(_("nonempty vector supplied for empty matrix"));
	else if (vlen > mlen)
		warning(_("vector length (%lld) exceeds matrix length (%d * %d)"),
		        (long long) vlen, m, n);
	else if (mlen % vlen != 0)
		warning(_("matrix length (%d * %d) is not a multiple of vector length (%lld)"),
		        m, n, (long long) vlen);

	return
	vector_as_dense(s_from, zzz, ul, ct, di, m, n, byrow, s_dimnames);
}

SEXP matrix_as_dense(SEXP from, const char *zzz,
                     char ul, char ct, char di, int mg, int new)
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
	PROTECT(from = coerceVector(from, tt));

	SEXP to = PROTECT(newObject(cl));
	int nprotect = 2;

	SEXP dim = getAttrib(from, R_DimSymbol), dimnames;
	int *pdim, isM, m, n, doDN;
	R_xlen_t mn = XLENGTH(from);

	isM = TYPEOF(dim) == INTSXP && LENGTH(dim) == 2;
	if (isM) {

		pdim = INTEGER(dim);
		m = pdim[0];
		n = pdim[1];
		if (m != n || n > 0) {
			dim = GET_SLOT(to, Matrix_DimSym);
			pdim = INTEGER(dim);
			pdim[0] = m;
			pdim[1] = n;
		}

		PROTECT(dimnames = getAttrib(from, R_DimNamesSymbol));
		++nprotect;
		doDN = dimnames != R_NilValue;

	} else {

		if (mn > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
		dim = GET_SLOT(to, Matrix_DimSym);
		pdim = INTEGER(dim);
		if (mg == 0) {
			pdim[0] = m = 1;
			pdim[1] = n = (int) mn;
		} else {
			pdim[0] = m = (int) mn;
			pdim[1] = n = 1;
		}

		SEXP nms = PROTECT(getAttrib(from, R_NamesSymbol));
		++nprotect;
		doDN = nms != R_NilValue;
		if (doDN) {
			PROTECT(dimnames = allocVector(VECSXP, 2));
			++nprotect;
			SET_VECTOR_ELT(dimnames, (mg == 0) ? 1 : 0, nms);
		}

	}

	if (cl[1] != 'g' && m != n)
		error(_("attempt to construct non-square %s"),
		      (cl[1] == 's' || cl[1] == 'p') ? "symmetricMatrix" : "triangularMatrix");

	if (doDN) {
		if (cl[1] == 's' || cl[1] == 'p')
			set_symmetrized_DimNames(to, dimnames, -1);
		else
			SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	}

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (cl[1] == 't' && di != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	SEXP x;

	if (!packed) {

	if (new <= 0 || (new <= 1 && ATTRIB(from) == R_NilValue) ||
	    !MAYBE_REFERENCED(from)) {

		if (ATTRIB(from) != R_NilValue && new >= 1) {
			/* 'from' has attributes and no references : */
			SET_ATTRIB(from, R_NilValue);
			if (OBJECT(from))
				SET_OBJECT(from, 0);
		}
		x = from;

	} else {

		PROTECT(x = duplicateVector(from));
		++nprotect;

	}

	} else {

		PROTECT(x = allocVector(tt, n + (mn - n) / 2));
		++nprotect;

#define MAD(c) c##NAME(pack2)(c##PTR(x), c##PTR(from), (size_t) n, \
                              ul, '\0', 'N')

		SWITCH4(cl[0], MAD);

#undef MAD

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
		error(_("second argument of '%s' does not specify a subclass of %s"),
		      __func__, "denseMatrix");

	char ul = '\0', ct = '\0', di = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , di);

	int mg;
	VALID_MARGIN(s_margin, mg);

	return matrix_as_dense(s_from, zzz, ul, ct, di, mg, 1);
}

SEXP sparse_as_dense(SEXP from, const char *class, int packed)
{
	packed = packed && class[1] != 'g';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = (packed) ? 'p' : ((class[1] == 'g') ? 'e' : ((class[1] == 's') ? 'y' : ((class[1] == 'p') ? 'o' : 'r')));
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n,
		xlen = (!packed) ? mn : n + (mn - n) / 2;
	if (xlen > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	if (class[2] != 'C' && packed && mn > R_XLEN_T_MAX)
		error(_("coercing n-by-n %s to %s is not supported for n*n exceeding %s"),
		      "[RT]sparseMatrix", "packedMatrix", "R_XLEN_T_MAX");
	double bytes = (double) xlen * kindToSize(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * bytes);
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

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

#define SAD(c) \
	do { \
		c##IF_NPATTERN( \
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)); \
		c##TYPE *px0 = c##PTR(x0); \
		); \
		SEXP x1 = PROTECT(allocVector(c##TYPESXP, (R_xlen_t) xlen)); \
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
						index = PACKED_AR21_UP((int_fast64_t) i, (int_fast64_t) *pj0); \
						px1[index] = c##IFELSE_NPATTERN(*px0, 1); \
						++k; ++pj0; c##IF_NPATTERN(++px0); \
					} \
				} \
			else \
				for (i = 0, k = 0; i < m; ++i) { \
					kend = pp0[i]; \
					while (k < kend) { \
						index = PACKED_AR21_LO((int_fast64_t) i, (int_fast64_t) *pj0, m1); \
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
					index = PACKED_AR21_UP((int_fast64_t) *pi0, (int_fast64_t) *pj0); \
					c##INCREMENT_IDEN(px1[index], *px0); \
					++pi0; ++pj0; c##IF_NPATTERN(++px0); \
				} \
			else \
				for (k = 0; k < kend; ++k) { \
					index = PACKED_AR21_LO((int_fast64_t) *pi0, (int_fast64_t) *pj0, m1); \
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

	SWITCH5(class[0], SAD);

#undef SAD

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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	int_fast64_t nn = (int_fast64_t) n * n,
		xlen = (!packed) ? nn : n + (nn - n) / 2;
	if (xlen > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	double bytes = (double) xlen * kindToSize(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * bytes);
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (cl[1] == 's' || cl[1] == 'p')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = CHAR(STRING_ELT(diag, 0))[0];
	if (cl[1] == 't' && di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	if (cl[0] != class[0]) {
		if (class[0] == 'n' && cl[0] == 'l')
			x0 = duplicateVector(x0);
		else
			x0 = coerceVector(x0, kindToType(cl[0]));
		if (class[0] == 'n')
			naToUnit(x0);
		UNPROTECT(1); /* x0 */
		PROTECT(x0);
	}

	SEXP x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) xlen));

#define DAD(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed) \
			c##NAME(force2)(px1, px0, (size_t) n, ul, 1, -di); \
		else \
			c##NAME(force1)(px1, px0, (size_t) n, ul, 1, -di); \
	} while (0)

	SWITCH4(cl[0], DAD);

#undef DAD

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
	SEXP margin = GET_SLOT(from, Matrix_marginSym);
	int mg = INTEGER(margin)[0] - 1;

	char cl[] = ".geMatrix";
	cl[0] = (kind == '.') ? 'n' : ((kind == ',') ? 'd' : kind);
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	int_fast64_t xlen = (int_fast64_t) m * n;
	if (xlen > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	double bytes = (double) xlen * kindToSize(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * bytes);
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP perm = PROTECT(GET_SLOT(from, Matrix_permSym));
	int *pperm = INTEGER(perm);

	SEXP x = PROTECT(allocVector(kindToType(cl[0]), (R_xlen_t) xlen));

#define IAD(c) \
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

	SWITCH4(cl[0], IAD);

#undef IAD

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

SEXP vector_as_sparse(SEXP from, const char *zzz,
                      char ul, char ct, char di,
                      int m, int n, int byrow, SEXP dimnames)
{
	SEXP length0 = GET_SLOT(from, Matrix_lengthSym);
	int_fast64_t r = (int_fast64_t)
		((TYPEOF(length0) == INTSXP) ? INTEGER(length0)[0] : REAL(length0)[0]);

	SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
		x0 = PROTECT(getAttrib(from, Matrix_xSym));

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
		x0 = coerceVector(x0, tt);
		UNPROTECT(1); /* x0 */
		PROTECT(x0);
	}

	if (cl[1] != 'g' && m != n)
		error(_("attempt to construct non-square %s"),
		      (cl[1] == 's' || cl[1] == 'p') ? "symmetricMatrix" : "triangularMatrix");

	SEXP to = PROTECT(newObject(cl));

	SEXP dim = GET_SLOT(to, Matrix_DimSym);
	int *pdim = INTEGER(dim);
	pdim[0] = m;
	pdim[1] = n;

	if (cl[1] == 's' || cl[1] == 'p')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (cl[1] == 't' && di != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	int_fast64_t pos, mn = (int_fast64_t) m * n, nnz1 = 0;
	R_xlen_t k = 0, nnz0 = XLENGTH(i0);

#define VAS__(d) \
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
		else if (cl[1] == 's' || cl[1] == 'p' || di == 'N') { \
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
		VAS__(i);
	else
		VAS__(d);

#undef VAS__

	if (nnz1 > INT_MAX)
		error(_("attempt to construct %s with more than %s nonzero entries"),
		      "sparseMatrix", "2^31-1");

	if (byrow) {
		int tmp;
		tmp = m; m = n; n = tmp;
	}
	
	SEXP iSym = (byrow) ? Matrix_jSym : Matrix_iSym,
		p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
		i1 = PROTECT(allocVector(INTSXP, nnz1));
	int *pp1 = INTEGER(p1) + 1, *pi1 = INTEGER(i1), i, j;
	memset(pp1 - 1, 0, sizeof(int) * ((size_t) n + 1));
	SET_SLOT(to, Matrix_pSym, p1);
	SET_SLOT(to,        iSym, i1);

	k = 0;

#define VAS__(d, c0, c1) \
	do { \
		d ##TYPE *pi0 = d ##PTR(i0); \
		c0##IF_NPATTERN( \
		c0##TYPE *px0 = c0##PTR(x0); \
		); \
		c1##IF_NPATTERN( \
		SEXP x1 = PROTECT(allocVector(c1##TYPESXP, nnz1)); \
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
		else if (cl[1] == 's' || cl[1] == 'p' || di == 'N') { \
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

#define VAS(c) \
		do { \
			if (TYPEOF(i0) == INTSXP) \
				VAS__(i, n, c); \
			else \
				VAS__(d, n, c); \
		} while (0)

		SWITCH5(cl[0], VAS);

#undef VAS

	} else {

#define VAS(c) \
		do { \
			if (TYPEOF(i0) == INTSXP) \
				VAS__(i, c, c); \
			else \
				VAS__(d, c, c); \
		} while (0)

		SWITCH5(cl[0], VAS);

#undef VAS

	}

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

SEXP R_vector_as_sparse(SEXP s_from, SEXP s_zzz,
                        SEXP s_uplo, SEXP s_trans, SEXP s_diag,
                        SEXP s_m, SEXP s_n, SEXP s_byrow, SEXP s_dimnames)
{
	Matrix_class(s_from, valid_vector, 0, __func__);

	const char *zzz;
	if (TYPEOF(s_zzz) != STRSXP || LENGTH(s_zzz) < 1 ||
	    (s_zzz = STRING_ELT(s_zzz, 0)) == NA_STRING ||
	    (zzz = CHAR(s_zzz))[0] == '\0' ||
	    (zzz[1] != 'g' && zzz[1] != 's' && zzz[1] != 'p' && zzz[1] != 't') ||
	    (zzz[2] != 'C' && zzz[2] != 'R' && zzz[2] != 'T'))
		error(_("second argument of '%s' does not specify a subclass of %s"),
		      __func__, "[CRT]sparseMatrix");

	char ul = '\0', ct = '\0', di = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , di);

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
					error(_("dimensions cannot exceed %s"), "2^31-1");
				m = (int) tmp;
			}
		}
		if (m < 0)
			error(_("invalid '%s' to '%s'"), "m", __func__);
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
					error(_("dimensions cannot exceed %s"), "2^31-1");
				n = (int) tmp;
			}
		}
		if (n < 0)
			error(_("invalid '%s' to '%s'"), "n", __func__);
	}

	int byrow;
	VALID_LOGIC2(s_byrow, byrow);

	if (s_dimnames != R_NilValue)
		if (TYPEOF(s_dimnames) != VECSXP || LENGTH(s_dimnames) != 2)
			error(_("invalid '%s' to '%s'"), "dimnames", __func__);

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
			error(_("dimensions cannot exceed %s"), "2^31-1");
		m = (int) vlen;
		n = 1;
	} else if (m < 0) {
		if (vlen > (int_fast64_t) INT_MAX * n) {
			if (n == 0)
				error(_("nonempty vector supplied for empty matrix"));
			else
				error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		m = (int) ((n == 0) ? 0 : vlen / n + (vlen % n != 0));
	} else if (n < 0) {
		if (vlen > (int_fast64_t) m * INT_MAX) {
			if (m == 0)
				error(_("nonempty vector supplied for empty matrix"));
			else
				error(_("dimensions cannot exceed %s"), "2^31-1");
		}
		n = (int) ((m == 0) ? 0 : vlen / m + (vlen % m != 0));
	}

	int_fast64_t mlen = (int_fast64_t) m * n;
	if (vlen <= 1)
		/* do nothing */ ;
	else if (mlen == 0)
		warning(_("nonempty vector supplied for empty matrix"));
	else if (vlen > mlen)
		warning(_("vector length (%lld) exceeds matrix length (%d * %d)"),
		        (long long) vlen, m, n);
	else if (mlen % vlen != 0)
		warning(_("matrix length (%d * %d) is not a multiple of vector length (%lld)"),
		        m, n, (long long) vlen);

	return
	vector_as_sparse(s_from, zzz, ul, ct, di, m, n, byrow, s_dimnames);
}

SEXP matrix_as_sparse(SEXP from, const char *zzz,
                      char ul, char ct, char di, int mg)
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
	REPROTECT(from = matrix_as_dense(from, cl, ul, ct, di, mg, 0), pid);
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
		error(_("second argument of '%s' does not specify a subclass of %s"),
		      __func__, "[CRT]sparseMatrix");

	char ul = '\0', ct = '\0', di = '\0';
	if (zzz[1] != 'g') VALID_UPLO (s_uplo , ul);
	if (zzz[1] == 's') VALID_TRANS(s_trans, ct);
	if (zzz[1] == 't') VALID_DIAG (s_diag , di);

	int mg;
	VALID_MARGIN(s_margin, mg);

	return matrix_as_sparse(s_from, zzz, ul, ct, di, mg);
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	SEXP p1, i1, j1;
	int i, j, *pp1, *pi1, *pj1;
	p1 = i1 = j1 = NULL;
	pp1 = pi1 = pj1 = NULL;

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	R_xlen_t nnz = 0;
	int nprotect = 2;

	if (cl[2] != 'T') {
		int r = (cl[2] == 'C') ? n : m;
		PROTECT(p1 = allocVector(INTSXP, (R_xlen_t) r + 1));
		++nprotect;
		pp1 = INTEGER(p1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);
	}

#define DAS_CHECK \
	do { \
		if (nnz > INT_MAX) \
			error(_("attempt to construct %s with more than %s nonzero entries"), \
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
		} else if (!packed && (di == '\0' || di == 'N')) { \
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
		} else if (di == '\0' || di == 'N') { \
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
		} else if (!packed && (di == '\0' || di == 'N')) { \
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
		} else if (di == '\0' || di == 'N')	{ \
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

#define DAS(c) \
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

	SWITCH5(class[0], DAS);

#undef DAS

	if (cl[2] != 'R') {
		PROTECT(i1 = allocVector(INTSXP, nnz));
		++nprotect;
		pi1 = INTEGER(i1);
		SET_SLOT(to, Matrix_iSym, i1);
	}
	if (cl[2] != 'C') {
		PROTECT(j1 = allocVector(INTSXP, nnz));
		++nprotect;
		pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_jSym, j1);
	}

#define DAS(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0); \
		c##IF_NPATTERN( \
		SEXP x1 = PROTECT(allocVector(c##TYPESXP, nnz)); \
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

	SWITCH5(class[0], DAS);

#undef DAS

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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (cl[1] == 's' || cl[1] == 'p')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	if (cl[1] == 's' && ct != 'C' && cl[0] == 'z') {
		SEXP trans = PROTECT(mkString("T"));
		SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = CHAR(STRING_ELT(diag, 0))[0];
	if (cl[1] == 't' && di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	if (di != 'N') {
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
			x0 = coerceVector(x0, kindToType(cl[0]));
		if (class[0] == 'n')
			naToUnit(x0);
		UNPROTECT(1); /* x0 */
		PROTECT(x0);
	}

	int j, nnz = 0;
	if (cl[2] == 'T') {

#define DAS(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			for (j = 0; j < n; ++j) { \
				if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
					++nnz; \
				px0 += 1; \
			} \
		} while (0)

		SWITCH4(cl[0], DAS);

#undef DAS

	} else {

		SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
		int *pp1 = INTEGER(p1);
		*(pp1++) = 0;
		SET_SLOT(to, Matrix_pSym, p1);

#define DAS(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			for (j = 0; j < n; ++j) { \
				if ((ct == 'C') ? c##NOT_ZERO_REAL(*px0) : c##NOT_ZERO(*px0)) \
					++nnz; \
				px0 += 1; \
				*(pp1++) = nnz; \
			} \
		} while (0)

		SWITCH4(cl[0], DAS);

#undef DAS

		UNPROTECT(1); /* p1 */

	}

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz));
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

#define DAS(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			c##IF_NPATTERN( \
			SEXP x1 = PROTECT(allocVector(c##TYPESXP, nnz)); \
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

		SWITCH5(cl[0], DAS);

#undef DAS

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
	SEXP margin = GET_SLOT(from, Matrix_marginSym);
	int mg = INTEGER(margin)[0] - 1;

	char cl[] = ".g.Matrix";
	cl[0] = (kind == '.') ? 'n' : ((kind == ',') ? 'd' : kind);
	cl[2] = (repr == '.') ? ((mg == 0) ? 'R' : 'C') : repr;
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1],
		r = (mg == 0) ? m : n, s = (mg == 0) ? n : m;
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP perm = PROTECT(GET_SLOT(from, Matrix_permSym));
	int k, *pperm = INTEGER(perm);

	if (cl[2] == 'T') {
		SEXP i = PROTECT(allocVector(INTSXP, r)),
			j = PROTECT(allocVector(INTSXP, r));
		int *pi = INTEGER(i), *pj = INTEGER(j);
		for (k = 0; k < r; ++k) {
			*(pi++) = k;
			*(pj++) = *(pperm++) - 1;
		}
		SET_SLOT(to, Matrix_iSym, (mg == 0) ? i : j);
		SET_SLOT(to, Matrix_jSym, (mg == 0) ? j : i);
		UNPROTECT(2); /* j, i */
	} else if ((cl[2] == 'C') == (mg != 0)) {
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) r + 1)),
			i = PROTECT(allocVector(INTSXP, r));
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
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) s + 1));
		int *pp = INTEGER(p);
		memset(pp, 0, sizeof(int) * ((size_t) s + 1));
		for (k = 0; k < r; ++k)
			++pp[pperm[k]];
		for (k = 0; k < s; ++k)
			pp[k + 1] += pp[k];
		SEXP j = PROTECT(allocVector(INTSXP, r));
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	PROTECT_INDEX pid;
	SEXP x;
	PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);

	if (TYPEOF(x) != tt) {
		REPROTECT(x = coerceVector(x, tt), pid);
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
		/* defined in ./sparse.c : */
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

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
		REPROTECT(x = coerceVector(x, tt), pid);
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = CHAR(STRING_ELT(diag, 0))[0];
	if (di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	if (di == 'N') {
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
			REPROTECT(x = coerceVector(x, tt), pid);
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if ((int_fast64_t) n * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's' || class[1] == 'p')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP uplo = GET_SLOT(from, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	char ct = '\0';
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}

	char di = '\0';
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		di = CHAR(STRING_ELT(diag, 0))[0];
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1 = x0;
	int nprotect = 2;
	if (packed || new) {
		PROTECT(x1 = allocVector(TYPEOF(x0), (R_xlen_t) n * n));
		++nprotect;
	}
#define DAG(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (packed) \
			c##NAME( pack1)(px1,  px0, (size_t) n, ul, ct, di); \
		else if (new) \
			c##NAME(force2)(px1,  px0, (size_t) n, ul, ct, di); \
		else \
			c##NAME(force2)(px1, NULL, (size_t) n, ul, ct, di); \
	} while (0)

	SWITCH4(class[0], DAG);

#undef DAG

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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's' || class[1] == 'p')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] == 't') {
		SEXP diag = GET_SLOT(from, Matrix_diagSym);
		char di = CHAR(STRING_ELT(diag, 0))[0];
		if (di == 'N') {
			if (class[2] != 'T') {
				SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym));
				SET_SLOT(to, Matrix_pSym, p0);
				UNPROTECT(1); /* p0 */
			}
			if (class[2] != 'R') {
				SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym));
				SET_SLOT(to, Matrix_iSym, i0);
				UNPROTECT(1); /* i0 */
			}
			if (class[2] != 'C') {
				SEXP j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
				SET_SLOT(to, Matrix_jSym, j0);
				UNPROTECT(1); /* j0 */
			}
			if (class[0] != 'n') {
				SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
				SET_SLOT(to, Matrix_xSym, x0);
				UNPROTECT(1); /* x0 */
			}
			UNPROTECT(1); /* to */
			return to;
		}
	} else {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */
	}

	SEXP uplo = GET_SLOT(from, Matrix_uploSym);
	char ul = CHAR(STRING_ELT(uplo, 0))[0];

	char ct = '\0';
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(from, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			p1 = PROTECT(allocVector(INTSXP, XLENGTH(p0)));
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
				error(_("attempt to construct %s with more than %s nonzero entries"),
				      "sparseMatrix", "2^31-1");
			for (j = 0; j < n; ++j)
				pp1[j] += pp0[j];
		} else {
			if (n > INT_MAX - pp0[n - 1])
				error(_("attempt to construct %s with more than %s nonzero entries"),
				      "sparseMatrix", "2^31-1");
			for (j = 0; j < n; ++j)
				pp1[j] = pp0[j] + j + 1;
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[n - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(to, iSym, i1);

#define SAG(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, pp1[n - 1])); \
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
			} else if (ul == ((class[2] == 'C') ? 'U' : 'L')) { \
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

		SWITCH5(class[0], SAG);

#undef SAG

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
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");
		nnz1 += nnz0;

		SEXP i1 = PROTECT(allocVector(INTSXP, nnz1)),
			j1 = PROTECT(allocVector(INTSXP, nnz1));
		int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, j1);

#define SAG(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, nnz1)); \
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

		SWITCH5(class[0], SAG);

#undef SAG

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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if ((int_fast64_t) n * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = CHAR(STRING_ELT(uplo, 0))[0];
	if (ul != 'U')
		SET_SLOT(to, Matrix_uploSym, uplo);
	UNPROTECT(1); /* uplo */

	if (cl[1] == 's' && cl[0] == 'z') {
		SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
		char ct = CHAR(STRING_ELT(trans, 0))[0];
		if (ct != 'C')
			SET_SLOT(to, Matrix_transSym, trans);
		UNPROTECT(1); /* trans */
	}

	if (cl[1] == 't') {
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
	}

	if (cl[1] == 'o') {
		SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
		if (LENGTH(sd) > 0)
			SET_SLOT(to, Matrix_sdSym, sd);
		UNPROTECT(1); /* sd */
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) n * n));

#define UNPACK(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		c##NAME(pack1)(px1, px0, (size_t) n, ul, '\0', 'N'); \
	} while (0)

	SWITCH4((cl[0] == 'c') ? 'd' : cl[0], UNPACK);

#undef UNPACK

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

SEXP dense_as_packed(SEXP from, const char *class, char ul, char ct, char di)
{
	if (class[2] == 'p')
		return from;

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = 'p';
	if (class[1] == 'g')
		cl[1] = (di == '\0') ? 's' : 't';
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), n = pdim[0];
	if (pdim[1] != n)
		error(_("attempt to pack non-square matrix"));
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] == 'g') {
		if (ul == '\0')
			ul = 'U';
		if (ul != 'U') {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}

		if (ct == '\0')
			ct = 'C';
		if (cl[1] == 's' && ct != 'C' && cl[0] == 'z') {
			SEXP trans = PROTECT(mkString("T"));
			SET_SLOT(to, Matrix_transSym, trans);
			UNPROTECT(1); /* trans */
		}

		if (cl[1] == 't' && di != 'N') {
			SEXP diag = PROTECT(mkString("U"));
			SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	} else {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (cl[1] == 's' && cl[0] == 'z') {
			SEXP trans = PROTECT(GET_SLOT(from, Matrix_transSym));
			ct = CHAR(STRING_ELT(trans, 0))[0];
			if (ct != 'C')
				SET_SLOT(to, Matrix_transSym, trans);
			UNPROTECT(1); /* trans */
		}

		if (cl[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = CHAR(STRING_ELT(diag, 0))[0];
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorsSym, factors);
			UNPROTECT(1); /* factors */
		}

		if (cl[1] == 'o') {
			SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
			if (LENGTH(sd) > 0)
				SET_SLOT(to, Matrix_sdSym, sd);
			UNPROTECT(1); /* sd */
		}
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), PACKED_LENGTH((R_xlen_t) n)));

#define PACK(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		c##NAME(pack2)(px1, px0, (size_t) n, ul, '\0', 'N'); \
	} while (0)

	SWITCH4((cl[0] == 'c') ? 'd' : cl[0], PACK);

#undef PACK

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<denseMatrix>, "packedMatrix") */
SEXP R_dense_as_packed(SEXP s_from, SEXP s_uplo, SEXP s_trans, SEXP s_diag)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);

	char ul = '\0', ct = '\0', di = '\0';
	if (class[1] == 'g') {
	if (s_uplo  != R_NilValue) VALID_UPLO (s_uplo , ul);
	if (s_trans != R_NilValue) VALID_TRANS(s_trans, ct);
	if (s_diag  != R_NilValue) VALID_DIAG (s_diag , di);
	}

	return dense_as_packed(s_from, class, ul, ct, di);
}

void trans(SEXP p0, SEXP i0, SEXP x0, SEXP p1, SEXP i1, SEXP x1,
           int m, int n)
{
	int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1),
		*pi0 = INTEGER(i0), *pi1 = INTEGER(i1),
		i, j, k, kend, nnz = pp0[n];
	pp0++; *(pp1++) = 0;
	memset(pp1, 0, sizeof(int) * (size_t) m);
	for (k = 0; k < nnz; ++k)
		++pp1[pi0[k]];
	for (i = 0; i < m; ++i)
		pp1[i] += pp1[i - 1];

	int *pp1_;
	Matrix_Calloc(pp1_, m, int);
	memcpy(pp1_, pp1 - 1, sizeof(int) * (size_t) m);

#define TRANS(c) \
	do { \
		c##IF_NPATTERN( \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		); \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			while (k < kend) { \
				i = pi0[k]; \
				pi1[pp1_[i]] = j; \
				c##IF_NPATTERN( \
				px1[pp1_[i]] = px0[k]; \
				); \
				++pp1_[i]; \
				++k; \
			} \
		} \
	} while (0)

	SWITCH5((x0 && x1 && TYPEOF(x0) == TYPEOF(x1)) ? typeToKind(TYPEOF(x0)) : 'n', TRANS);

#undef TRANS

	Matrix_Free(pp1_, m);
	return;
}

void tsort(SEXP i0, SEXP j0, SEXP x0, SEXP *p1, SEXP *i1, SEXP *x1,
           int m, int n)
{
	if (XLENGTH(i0) > INT_MAX)
		error(_("unable to aggregate %s with '%s' and '%s' slots of length exceeding %s"),
		      "TsparseMatrix", "i", "j", "2^31-1");
	int nnz0 = (int) XLENGTH(i0), nnz1 = 0;

	PROTECT(*p1 = allocVector(INTSXP, (R_xlen_t) n + 1));
	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), *pp1 = INTEGER(*p1), *pi1,
		i, j, r = (m < n) ? n : m, k, kstart, kend, kend_;
	*(pp1++) = 0;

	int *workA, *workB, *workC, *pj_;
	int_fast64_t lwork = (int_fast64_t) m + r + m + nnz0;
	Matrix_Calloc(workA, lwork, int);
	workB = workA + m;
	workC = workB + r;
	pj_   = workC + m;

#define TSORT(c) \
	do { \
		c##IF_NPATTERN( \
		c##TYPE *px0 = c##PTR(x0), *px1, *px_; \
		Matrix_Calloc(px_, nnz0, c##TYPE); \
		); \
		 \
		/* 1. Tabulate column indices in workA[i]                         */ \
		 \
		for (k = 0; k < nnz0; ++k) \
			++workA[pi0[k]]; \
		 \
		/* 2. Compute cumulative sum in workA[i], copying to workB[i]     */ \
		/*                                                                */ \
		/* workA[i]: number of column indices listed for row i,           */ \
		/*           incl. duplicates                                     */ \
		 \
		for (i = 1; i < m; ++i) \
			workA[i] += (workB[i] = workA[i - 1]); \
		 \
		/* 3. Group column indices and data by row in pj_[k], px_[k]      */ \
		/*                                                                */ \
		/* workA[i]: number of column indices listed for row <= i,        */ \
		/*           incl. duplicates                                     */ \
		/* workB[i]: number of column indices listed for row <  i,        */ \
		/*           incl. duplicates                                     */ \
		 \
		for (k = 0; k < nnz0; ++k) { \
			pj_[workB[pi0[k]]] = pj0[k]; \
			c##IF_NPATTERN( \
			px_[workB[pi0[k]]] = px0[k]; \
			); \
			++workB[pi0[k]]; \
		} \
		 \
		/* 4. Gather _unique_ column indices at the front of each group,  */ \
		/*    aggregating data accordingly; record in workC[i] where      */ \
		/*    the unique column indices stop and the duplicates begin     */ \
		/*                                                                */ \
		/* workB[.]: unused                                               */ \
		/*   pj_[k]: column indices grouped by row,                       */ \
		/*           incl. duplicates, unsorted                           */ \
		/*   px_[k]: corresponding data                                   */ \
		 \
		k = 0; \
		for (j = 0; j < n; ++j) \
			workB[j] = -1; \
		for (i = 0; i < m; ++i) { \
			kstart = kend_ = k; \
			kend = workA[i]; \
			while (k < kend) { \
				if (workB[pj_[k]] < kstart) { \
					/* Have not yet seen this column index */ \
					workB[pj_[k]] = kend_; \
					pj_[kend_] = pj_[k]; \
					c##IF_NPATTERN( \
					px_[kend_] = px_[k]; \
					); \
					++kend_; \
				} else { \
					/* Have already seen this column index */ \
					c##IF_NPATTERN( \
					c##INCREMENT_IDEN(px_[workB[pj_[k]]], px_[k]); \
					); \
				} \
				++k; \
			} \
			workC[i] = kend_; \
			nnz1 += kend_ - kstart; \
		} \
		 \
		/* 5. Tabulate _unique_ column indices in workB[j]                */ \
		/*                                                                */ \
		/* workC[i]: pointer to first non-unique column index in row i    */ \
		/*   pi_[k]: column indices grouped by row,                       */ \
		/*           with unique indices in front,                        */ \
		/*           i.e., in positions workA[i - 1] <= k < workC[i]      */ \
		/*   px_[k]: corresponding data, "cumulated" appropriately        */ \
		 \
		k = 0; \
		memset(workB, 0, sizeof(int) * (size_t) n); \
		for (i = 0; i < m; ++i) { \
			kend_ = workC[i]; \
			while (k < kend_) { \
				++workB[pj_[k]]; \
				++k; \
			} \
			k = workA[i]; \
		} \
		 \
		/* 6. Compute cumulative sum in pp1[j], copying to workB[j]       */ \
		/*                                                                */ \
		/* workB[j]: number of nonzero elements in column j               */ \
		 \
		for (j = 0; j < n; ++j) { \
			pp1[j] = pp1[j - 1] + workB[j]; \
			workB[j] = pp1[j - 1]; \
		} \
		 \
		PROTECT(*i1 = allocVector(    INTSXP, nnz1)); \
		pi1 = INTEGER(*i1); \
		c##IF_NPATTERN( \
		PROTECT(*x1 = allocVector(c##TYPESXP, nnz1)); \
		px1 =  c##PTR(*x1); \
		); \
		 \
		/* 7. Pop unique (i,j) pairs from the unsorted stacks 0 <= i < m  */ \
		/*    onto new stacks 0 <= j < n, which will be sorted            */ \
		/*                                                                */ \
		/* workB[j]: number of nonzero elements in columns <  j           */ \
		/*   pp1[j]: number of nonzero elements in columns <= j           */ \
		 \
		k = 0; \
		for (i = 0; i < m; ++i) { \
			kend_ = workC[i]; \
			while (k < kend_) { \
				pi1[workB[pj_[k]]] = i; \
				c##IF_NPATTERN( \
				px1[workB[pj_[k]]] = px_[k]; \
				); \
				++workB[pj_[k]]; \
				++k; \
			} \
			k = workA[i]; \
		} \
		 \
		c##IF_NPATTERN( \
		Matrix_Free(px_, nnz0); \
		UNPROTECT(1); /* *x1 */ \
		); \
		UNPROTECT(1); /* *i1 */ \
	} while (0)

	SWITCH5((x0) ? typeToKind(TYPEOF(x0)) : 'n', TSORT);

#undef TSORT

	Matrix_Free(workA, lwork);
	UNPROTECT(1); /* *p1 */
	return;
}

void taggr(SEXP i0, SEXP j0, SEXP x0, SEXP *i1, SEXP *j1, SEXP *x1,
           int m, int n)
{
	if (XLENGTH(i0) > INT_MAX)
		error(_("unable to aggregate %s with '%s' and '%s' slots of length exceeding %s"),
		      "TsparseMatrix", "i", "j", "2^31-1");
	int nnz0 = (int) XLENGTH(i0), nnz1 = 0;

	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), *pi1, *pj1,
		i, j, r = (m < n) ? n : m, k, kstart, kend, kend_;

	int *workA, *workB, *workC, *pj_;
	int_fast64_t lwork = (int_fast64_t) m + r + m + nnz0;
	Matrix_Calloc(workA, lwork, int);
	workB = workA + m;
	workC = workB + r;
	pj_   = workC + m;

#define TAGGR(c) \
	do { \
		c##IF_NPATTERN( \
		c##TYPE *px0 = c##PTR(x0), *px1, *px_; \
		Matrix_Calloc(px_, nnz0, c##TYPE); \
		); \
		 \
		/* 1. Tabulate column indices in workA[i]                         */ \
		 \
		for (k = 0; k < nnz0; ++k) \
			++workA[pi0[k]]; \
		 \
		/* 2. Compute cumulative sum in workA[i], copying to workB[i]     */ \
		/*                                                                */ \
		/* workA[i]: number of column indices listed for row i,           */ \
		/*           incl. duplicates                                     */ \
		 \
		for (i = 1; i < m; ++i) \
			workA[i] += (workB[i] = workA[i - 1]); \
		 \
		/* 3. Group column indices and data by row in pj_[k], px_[k]      */ \
		/*                                                                */ \
		/* workA[i]: number of column indices listed for row <= i,        */ \
		/*           incl. duplicates                                     */ \
		/* workB[i]: number of column indices listed for row <  i,        */ \
		/*           incl. duplicates                                     */ \
		 \
		for (k = 0; k < nnz0; ++k) { \
			pj_[workB[pi0[k]]] = pj0[k]; \
			c##IF_NPATTERN( \
			px_[workB[pi0[k]]] = px0[k]; \
			); \
			++workB[pi0[k]]; \
		} \
		 \
		/* 4. Gather _unique_ column indices at the front of each group,  */ \
		/*    aggregating data accordingly; record in workC[i] where      */ \
		/*    the unique column indices stop and the duplicates begin     */ \
		/*                                                                */ \
		/* workB[.]: unused                                               */ \
		/*   pj_[k]: column indices grouped by row,                       */ \
		/*           incl. duplicates, unsorted                           */ \
		/*   px_[k]: corresponding data                                   */ \
		 \
		k = 0; \
		for (j = 0; j < n; ++j) \
			workB[j] = -1; \
		for (i = 0; i < m; ++i) { \
			kstart = kend_ = k; \
			kend = workA[i]; \
			while (k < kend) { \
				if (workB[pj_[k]] < kstart) { \
					/* Have not yet seen this column index */ \
					workB[pj_[k]] = kend_; \
					pj_[kend_] = pj_[k]; \
					c##IF_NPATTERN( \
					px_[kend_] = px_[k]; \
					); \
					++kend_; \
				} else { \
					/* Have already seen this column index */ \
					c##IF_NPATTERN( \
					c##INCREMENT_IDEN(px_[workB[pj_[k]]], px_[k]); \
					); \
				} \
				++k; \
			} \
			workC[i] = kend_; \
			nnz1 += kend_ - kstart; \
		} \
		if (nnz1 != nnz0) { \
			PROTECT(*i1 = allocVector(    INTSXP, nnz1)); \
			PROTECT(*j1 = allocVector(    INTSXP, nnz1)); \
			pi1 = INTEGER(*i1); \
			pj1 = INTEGER(*j1); \
			c##IF_NPATTERN( \
			PROTECT(*x1 = allocVector(c##TYPESXP, nnz1)); \
			px1 =  c##PTR(*x1); \
			); \
			 \
			k = 0; \
			for (i = 0; i < m; ++i) { \
				kend_ = workC[i]; \
				while (k < kend_) { \
					*(pi1++) =      i ; \
					*(pj1++) = pj_[k] ; \
					c##IF_NPATTERN( \
					*(px1++) = px_[k]; \
					); \
					++k; \
				} \
				k = workA[i]; \
			} \
			 \
			c##IF_NPATTERN( \
			UNPROTECT(1); /* x1 */ \
			); \
			UNPROTECT(2); /* j1, i1 */ \
		} \
		c##IF_NPATTERN( \
		Matrix_Free(px_, nnz0); \
		); \
	} while (0)

	SWITCH5((x0) ? typeToKind(TYPEOF(x0)) : 'n', TAGGR);

#undef TAGGR

	Matrix_Free(workA, lwork);
	return;
}

SEXP sparse_as_Csparse(SEXP from, const char *class)
{
	if (class[2] == 'C')
		return from;

	char cl[] = "..CMatrix";
	cl[0] = class[0];
	cl[1] = class[1];
	SEXP to = PROTECT(newObject(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
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
	}

	if (class[2] == 'R') {
		SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i1 = PROTECT(allocVector(INTSXP, INTEGER(p0)[m]));
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to, Matrix_iSym, i1);
		if (class[0] == 'n')
			trans(p0, j0, NULL, p1, i1, NULL, n, m);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), INTEGER(p0)[m]));
			SET_SLOT(to, Matrix_xSym, x1);
			trans(p0, j0, x0, p1, i1, x1, n, m);
			UNPROTECT(2); /* x1, x0 */
		}
		UNPROTECT(4); /* i1, p1, j0, p0 */
	} else {
		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			p1 = NULL, i1 = NULL;
		if (class[0] == 'n') {
			tsort(i0, j0, NULL, &p1, &i1, NULL, m, n);
			PROTECT(p1);
			PROTECT(i1);
			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to, Matrix_iSym, i1);
			UNPROTECT(2); /* i1, p1 */
		} else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = NULL;
			tsort(i0, j0, x0, &p1, &i1, &x1, m, n);
			PROTECT(p1);
			PROTECT(i1);
			PROTECT(x1);
			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(4); /* x1, i1, p1, x0 */
		}
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
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
	}

	if (class[2] == 'C') {
		SEXP p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) m + 1)),
			j1 = PROTECT(allocVector(INTSXP, INTEGER(p0)[n]));
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to, Matrix_jSym, j1);
		if (class[0] == 'n')
			trans(p0, i0, NULL, p1, j1, NULL, m, n);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), INTEGER(p0)[n]));
			SET_SLOT(to, Matrix_xSym, x1);
			trans(p0, i0, x0, p1, j1, x1, m, n);
			UNPROTECT(2); /* x1, x0 */
		}
		UNPROTECT(4); /* j1, p1, i0, p0 */
	} else {
		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			p1 = NULL, j1 = NULL;
		if (class[0] == 'n') {
			tsort(j0, i0, NULL, &p1, &j1, NULL, n, m);
			PROTECT(p1);
			PROTECT(j1);
			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to, Matrix_jSym, j1);
			UNPROTECT(2); /* j1, p1 */
		} else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = NULL;
			tsort(j0, i0, x0, &p1, &j1, &x1, n, m);
			PROTECT(p1);
			PROTECT(j1);
			PROTECT(x1);
			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(4); /* x1, j1, p1, x0 */
		}
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

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		char ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (ul != 'U')
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
	}

	SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(from, iSym));
	int *pp0 = INTEGER(p0) + 1, r = (class[2] == 'C') ? n : m, nnz = pp0[r - 1];
	if (XLENGTH(i0) == nnz)
		SET_SLOT(to, iSym, i0);
	else {
		SEXP i1 = PROTECT(allocVector(INTSXP, nnz));
		memcpy(INTEGER(i1), INTEGER(i0), sizeof(int) * (size_t) nnz);
		SET_SLOT(to, iSym, i1);
		UNPROTECT(1); /* i1 */
	}
	SEXP jSym = (class[2] == 'C') ? Matrix_jSym : Matrix_iSym,
		j1 = PROTECT(allocVector(INTSXP, nnz));
	int *pj1 = INTEGER(j1), j, k, kend;
	SET_SLOT(to, jSym, j1);
	for (j = 0, k = 0; j < r; ++j) {
		kend = pp0[j];
		while (k < kend)
			pj1[k++] = j;
	}
	UNPROTECT(3); /* j1, i0, p0 */

	if (class[0] != 'n') {
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		if (XLENGTH(x0) == nnz)
			SET_SLOT(to, Matrix_xSym, x0);
		else {
			SEXP x1 = PROTECT(allocVector(TYPEOF(x0), nnz));

#define SAT(c) memcpy(c##PTR(x1), c##PTR(x0), sizeof(c##TYPE) * (size_t) nnz)

			SWITCH4(class[0], SAT);

#undef SAT

			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */
		}
		UNPROTECT(1); /* x0 */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "TsparseMatrix") */
SEXP R_sparse_as_Tsparse(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_sparse, 2, __func__);
	return sparse_as_Tsparse(s_from, class);
}

/* as(<Matrix>, "vector") */
SEXP R_Matrix_as_vector(SEXP s_from)
{
	const char *class = Matrix_class(s_from, valid_matrix, 7, __func__);

	SEXP to = NULL;
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

	SEXP to = NULL;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(s_from, &pid);

	switch (class[2]) {
	case 'e':
		PROTECT(to = GET_SLOT(s_from, Matrix_xSym));
		to = duplicate(to);
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
	setAttrib(to, R_DimSymbol, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(s_from, Matrix_DimNamesSym));
	if (!DimNames_is_trivial(dimnames))
		setAttrib(to, R_DimNamesSymbol, dimnames);
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
		error(_("attempt to pack a %s"), "generalMatrix");
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
		error(_("attempt to pack an %s"), "indMatrix");
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
