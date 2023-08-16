#include "coerce.h"

SEXP matrix_as_dense(SEXP from, const char *zzz, char ul, char di,
                     int transpose_if_vector, int new)
{
	SEXPTYPE tf = TYPEOF(from);
	char cl[] = "...Matrix";
	cl[0] = (zzz[0] == '.') ? type2kind(tf) : zzz[0];
	cl[1] = zzz[1];
	cl[2] = zzz[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl)),
		dim = getAttrib(from, R_DimSymbol),
		dimnames;
	int *pdim, m, n, isM, doDN, nprotect = 1;
	R_xlen_t len = XLENGTH(from);

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

		if (len > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
		dim = GET_SLOT(to, Matrix_DimSym);
		pdim = INTEGER(dim);
		if (transpose_if_vector) {
			pdim[0] = m = 1;
			pdim[1] = n = (int) len;
		} else {
			pdim[0] = m = (int) len;
			pdim[1] = n = 1;
		}

		SEXP nms = PROTECT(getAttrib(from, R_NamesSymbol));
		++nprotect;
		doDN = nms != R_NilValue;
		if (doDN) {
			PROTECT(dimnames = allocVector(VECSXP, 2));
			++nprotect;
			SET_VECTOR_ELT(dimnames, transpose_if_vector ? 1 : 0, nms);
		}

	}

	if (cl[1] != 'g' && m != n)
		error(_("attempt to construct %s or %s from non-square matrix"),
		      "symmetricMatrix", "triangularMatrix");

	if (doDN) {
		if (cl[1] != 's')
			SET_SLOT(to, Matrix_DimNamesSym, dimnames);
		else
			set_symmetrized_DimNames(to, dimnames, -1);
	}

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	if (cl[1] == 't' && di != 'N') {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	SEXPTYPE tt = kind2type(cl[0]);
	if (tf != tt) {
		PROTECT(from = coerceVector(from, tt));
		++nprotect;
	}

	SEXP x;

	if (cl[2] != 'p') {

		if (!new || tf != tt || !MAYBE_REFERENCED(from)) {
			x = from;
			if (new && ATTRIB(x) != R_NilValue) {
				SET_ATTRIB(x, R_NilValue);
				if (OBJECT(x))
					SET_OBJECT(x, 0);
			}
		} else {
			PROTECT(x = allocVector(tt, len));
			++nprotect;
			switch (tt) {
			case LGLSXP:
				Matrix_memcpy(LOGICAL(x), LOGICAL(from), len, sizeof(int));
				break;
			case INTSXP:
				Matrix_memcpy(INTEGER(x), INTEGER(from), len, sizeof(int));
				break;
			case REALSXP:
				Matrix_memcpy(REAL(x), REAL(from), len, sizeof(double));
				break;
			case CPLXSXP:
				Matrix_memcpy(COMPLEX(x), COMPLEX(from), len, sizeof(Rcomplex));
				break;
			default:
				break;
			}
		}

	} else {

		PROTECT(x = allocVector(tt, PM_LENGTH(n)));
		++nprotect;
		switch (tt) {
		case LGLSXP:
			idense_pack(LOGICAL(x), LOGICAL(from), n, ul, di);
			break;
		case INTSXP:
			idense_pack(INTEGER(x), INTEGER(from), n, ul, di);
			break;
		case REALSXP:
			ddense_pack(REAL(x), REAL(from), n, ul, di);
			break;
		case CPLXSXP:
			zdense_pack(COMPLEX(x), COMPLEX(from), n, ul, di);
			break;
		default:
			break;
		}

	}

	SET_SLOT(to, Matrix_xSym, x);

	UNPROTECT(nprotect);
	return to;
}

/* as(<matrix|vector>, ".(ge|sy|sp|tr|tp)Matrix") */
SEXP R_matrix_as_dense(SEXP from, SEXP class, SEXP uplo, SEXP diag)
{
	switch (TYPEOF(from)) {
	case LGLSXP:
#ifdef MATRIX_ENABLE_IMATRIX
	case INTSXP:
#endif
	case REALSXP:
#ifdef MATRIX_ENABLE_ZMATRIX
	case CPLXSXP:
#endif
		break;
#ifndef MATRIX_ENABLE_IMATRIX
	case INTSXP:
		if (!inherits(from, "factor"))
			break;
#endif
	default:
		ERROR_INVALID_CLASS(from, __func__);
		break;
	}

	const char *zzz;
	if (TYPEOF(class) != STRSXP || LENGTH(class) < 1 ||
	    (class = STRING_ELT(class, 0)) == NA_STRING ||
	    (zzz = CHAR(class))[0] == '\0' ||
	    (zzz              )[1] == '\0' ||
	    !((zzz[1] == 'g' && (zzz[2] == 'e'                 )) ||
	      (zzz[1] == 's' && (zzz[2] == 'y' || zzz[2] == 'p')) ||
	      (zzz[1] == 't' && (zzz[2] == 'r' || zzz[2] == 'p'))))
		error(_("invalid '%s' to %s()"), "class", __func__);

	char ul = 'U', di = 'N';
	if (zzz[1] != 'g') {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid '%s' to %s()"), "uplo", __func__);
		if (zzz[1] == 't') {
			if (TYPEOF(diag) != STRSXP || LENGTH(diag) < 1 ||
			    (diag = STRING_ELT(diag, 0)) == NA_STRING ||
			    ((di = *CHAR(diag)) != 'N' && di != 'U'))
				error(_("invalid '%s' to %s()"), "diag", __func__);
		}
	}

	return matrix_as_dense(from, zzz, ul, di, 0, 1);
}

SEXP sparse_as_dense(SEXP from, const char *class, int packed)
{
	packed = packed && class[1] != 'g';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = (packed) ? 'p' :
		((class[1] == 'g') ? 'e' : ((class[1] == 's') ? 'y' : 'r'));
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	Matrix_int_fast64_t len = (Matrix_int_fast64_t) m * n;
	if (packed)
		len = (len + n) / 2;
	if (len > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	if (class[2] != 'C' && packed && len > R_XLEN_T_MAX)
		error(_("coercing n-by-n %s to %s is not supported for n*n exceeding %s"),
		      "[RT]sparseMatrix", "packedMatrix", "R_XLEN_T_MAX");
	double bytes = (double) len * kind2size(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * bytes);
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	}

	/* It remains to fill 'x' ... */

	SEXP x1 = PROTECT(allocVector(kind2type(class[0]), (R_xlen_t) len)),
		p0, i0, j0;
	int *pp, *pi, *pj, nprotect = 2;
	p0 = i0 = j0 = NULL;
	pp = pi = pj = NULL;
	SET_SLOT(to, Matrix_xSym, x1);

	if (class[2] != 'T') {
		PROTECT(p0 = GET_SLOT(from, Matrix_pSym));
		++nprotect;
		pp = INTEGER(p0) + 1;
	}
	if (class[2] != 'R') {
		PROTECT(i0 = GET_SLOT(from, Matrix_iSym));
		++nprotect;
		pi = INTEGER(i0);
	}
	if (class[2] != 'C') {
		PROTECT(j0 = GET_SLOT(from, Matrix_jSym));
		++nprotect;
		pj = INTEGER(j0);
	}

#define SAD_CASES \
	do { \
		switch (class[0]) { \
		case 'l': \
			SAD_SUBCASES(int, LOGICAL, SHOW, FIRSTOF, INCREMENT_LOGICAL, 1); \
			break; \
		case 'i': \
			SAD_SUBCASES(int, INTEGER, SHOW, FIRSTOF, INCREMENT_INTEGER, 1); \
			break; \
		case 'd': \
			SAD_SUBCASES(double, REAL, SHOW, FIRSTOF, INCREMENT_REAL, 1.0); \
			break; \
		case 'z': \
			SAD_SUBCASES(Rcomplex, COMPLEX, SHOW, FIRSTOF, INCREMENT_COMPLEX, Matrix_zone); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define SAD_SUBCASES(_CTYPE_, _PTR_, _MASK_, _REPLACE_, _INCREMENT_, _ONE_) \
	do { \
		_MASK_(_CTYPE_ *px0 = _PTR_(x0)); \
		       _CTYPE_ *px1 = _PTR_(x1) ; \
		Matrix_memset(px1, 0, len, sizeof(_CTYPE_)); \
		if (!packed) { \
			/* .(ge|sy|tr)Matrix */ \
			SAD_SUBSUBCASES(SAD_LOOP_C2NP, SAD_LOOP_R2NP, SAD_LOOP_T2NP, \
			                _MASK_, _REPLACE_, _INCREMENT_); \
			if (class[1] == 't' && di != 'N') { \
				px1 = _PTR_(x1); \
				R_xlen_t n1a = (R_xlen_t) n + 1; \
				for (int d = 0; d < n; ++d, px1 += n1a) \
					*px1 = _ONE_; \
			} \
		} else if (ul == 'U') { \
			/* upper triangular .(sp|tp)Matrix */ \
			SAD_SUBSUBCASES(SAD_LOOP_C2UP, SAD_LOOP_R2UP, SAD_LOOP_T2UP, \
			                _MASK_, _REPLACE_, _INCREMENT_); \
			if (class[1] == 't' && di != 'N') { \
				px1 = _PTR_(x1); \
				for (int d = 0; d < n; px1 += (++d) + 1) \
					*px1 = _ONE_; \
			} \
		} else { \
			/* lower triangular .(sp|tp)Matrix */ \
			SAD_SUBSUBCASES(SAD_LOOP_C2LP, SAD_LOOP_R2LP, SAD_LOOP_T2LP, \
			                _MASK_, _REPLACE_, _INCREMENT_); \
			if (class[1] == 't' && di != 'N') { \
				px1 = _PTR_(x1); \
				for (int d = 0; d < n; px1 += n - (d++)) \
					*px1 = _ONE_; \
			} \
		} \
	} while (0)

#define SAD_SUBSUBCASES(_LOOP_C_, _LOOP_R_, _LOOP_T_, _MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		switch (class[2]) { \
		case 'C': \
		{ \
			int j, k, kend; \
			_LOOP_C_(_MASK_, _REPLACE_, _INCREMENT_); \
			break; \
		} \
		case 'R': \
		{ \
			int i, k, kend; \
			_LOOP_R_(_MASK_, _REPLACE_, _INCREMENT_); \
			break; \
		} \
		case 'T': \
		{ \
			R_xlen_t index, k, nnz = XLENGTH(i0); \
			_LOOP_T_(_MASK_, _REPLACE_, _INCREMENT_); \
			break; \
		} \
		default: \
			break; \
		} \
	} while (0)

#define SAD_LOOP_C2NP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		for (j = 0, k = 0; j < n; ++j, px1 += m) { \
			kend = pp[j]; \
			while (k < kend) { \
				px1[*pi] = _REPLACE_(*px0, 1); \
				++k; ++pi; _MASK_(++px0); \
			} \
		} \
	} while (0)

#define SAD_LOOP_C2UP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		for (j = 0, k = 0; j < n; px1 += (++j)) { \
			kend = pp[j]; \
			while (k < kend) { \
				px1[*pi] = _REPLACE_(*px0, 1); \
				++k; ++pi; _MASK_(++px0); \
			} \
		} \
	} while (0)

#define SAD_LOOP_C2LP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		for (j = 0, k = 0; j < n; px1 += n - (j++)) { \
			kend = pp[j]; \
			while (k < kend) { \
				px1[*pi - j] = _REPLACE_(*px0, 1); \
				++k; ++pi; _MASK_(++px0); \
			} \
		} \
	} while (0)

#define SAD_LOOP_R2NP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		R_xlen_t m1 = (R_xlen_t) m; \
		for (i = 0, k = 0; i < m; ++i, ++px1) { \
			kend = pp[i]; \
			while (k < kend) { \
				px1[m1 * *pj] = _REPLACE_(*px0, 1); \
				++k; ++pj; _MASK_(++px0); \
			} \
		} \
	} while (0)

#define SAD_LOOP_R2UP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		for (i = 0, k = 0; i < n; ++i) { \
			kend = pp[i]; \
			while (k < kend) { \
				px1[PM_AR21_UP(i, *pj)] = _REPLACE_(*px0, 1); \
				++k; ++pj; _MASK_(++px0); \
			} \
		} \
	} while (0)

#define SAD_LOOP_R2LP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		R_xlen_t n2 = (R_xlen_t) n * 2; \
		for (i = 0, k = 0; i < n; ++i) { \
			kend = pp[i]; \
			while (k < kend) { \
				px1[PM_AR21_LO(i, *pj, n2)] = _REPLACE_(*px0, 1); \
				++k; ++pj; _MASK_(++px0); \
			} \
		} \
	} while (0)

#define SAD_LOOP_T2NP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		R_xlen_t m1 = (R_xlen_t) m; \
		for (k = 0; k < nnz; ++k) { \
			index = m1 * *pj + *pi; \
			_INCREMENT_(px1[index], (*px0)); \
			++pi; ++pj; _MASK_(++px0); \
		} \
	} while (0)

#define SAD_LOOP_T2UP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		for (k = 0; k < nnz; ++k) { \
			index = PM_AR21_UP(*pi, *pj); \
			_INCREMENT_(px1[index], (*px0)); \
			++pi; ++pj; _MASK_(++px0); \
		} \
	} while (0)

#define SAD_LOOP_T2LP(_MASK_, _REPLACE_, _INCREMENT_) \
	do { \
		R_xlen_t n2 = (R_xlen_t) n * 2; \
		for (k = 0; k < nnz; ++k) { \
			index = PM_AR21_LO(*pi, *pj, n2); \
			_INCREMENT_(px1[index], (*px0)); \
			++pi; ++pj; _MASK_(++px0); \
		} \
	} while (0)

	if (class[0] == 'n')
		SAD_SUBCASES(int, LOGICAL, HIDE, SECONDOF, INCREMENT_PATTERN, 1);
	else {
		SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		SAD_CASES;
		UNPROTECT(1); /* x0 */
	}

#undef SAD_CASES
#undef SAD_SUBCASES
#undef SAD_SUBSUBCASES
#undef SAD_LOOP_C2NP
#undef SAD_LOOP_C2UP
#undef SAD_LOOP_C2LP
#undef SAD_LOOP_R2NP
#undef SAD_LOOP_R2UP
#undef SAD_LOOP_R2LP
#undef SAD_LOOP_T2NP
#undef SAD_LOOP_T2UP
#undef SAD_LOOP_T2LP

	UNPROTECT(nprotect);
	return to;
}

/* as(<[CRT]sparseMatrix>, "(un)?packedMatrix") */
SEXP R_sparse_as_dense(SEXP from, SEXP packed)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	int packed_;
	if (TYPEOF(packed) != LGLSXP || LENGTH(packed) < 1 ||
	    (packed_ = LOGICAL(packed)[0]) == NA_LOGICAL)
		error(_("invalid '%s' to %s()"), "packed", __func__);

	return sparse_as_dense(from, valid[ivalid], packed_);
}

SEXP diagonal_as_dense(SEXP from, const char *class,
                       char shape, int packed, char ul)
{
	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = shape;
	cl[2] = (cl[1] == 'g') ? 'e' : ((packed) ? 'p' : ((cl[1] == 's') ? 'y' : 'r'));
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	Matrix_int_fast64_t len = (Matrix_int_fast64_t) n * n;
	if (len > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	double bytes = (double) len * kind2size(cl[0]);
	if (bytes > 0x1.0p+30 /* 1 GiB */)
		warning(_("sparse->dense coercion: allocating vector of size %0.1f GiB"),
		        0x1.0p-30 * bytes);
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (cl[1] == 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (cl[1] == 't' && di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) len));
	SET_SLOT(to, Matrix_xSym, x1);

#define DAD_SUBCASES(_CTYPE_, _PTR_, _PREFIX_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		Matrix_memset(px1, 0, (R_xlen_t) len, sizeof(_CTYPE_)); \
		if (di == 'N' || cl[1] != 't') { \
			if (cl[2] != 'p') \
				_PREFIX_ ## dense_unpacked_copy_diagonal( \
					px1, px0, n, n,     ul, di); \
			else \
				_PREFIX_ ## dense_packed_copy_diagonal( \
					px1, px0, n, n, ul, ul, di); \
		} \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		DAD_SUBCASES(int, LOGICAL, i);
		break;
	case 'i':
		DAD_SUBCASES(int, INTEGER, i);
		break;
	case 'd':
		DAD_SUBCASES(double, REAL, d);
		break;
	case 'z':
		DAD_SUBCASES(Rcomplex, COMPLEX, z);
		break;
	default:
		break;
	}

#undef DAD_CASES
#undef DAD_SUBCASES

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<diagonalMatrix>, ".(ge|sy|sp|tr|tp)Matrix") */
SEXP R_diagonal_as_dense(SEXP from, SEXP shape, SEXP packed, SEXP uplo)
{
	static const char *valid[] = { VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char shape_;
	if (TYPEOF(shape) != STRSXP || LENGTH(shape) < 1 ||
	    (shape = STRING_ELT(shape, 0)) == NA_STRING ||
	    ((shape_ = CHAR(shape)[0]) != 'g' && shape_ != 's' && shape_ != 't'))
		error(_("invalid '%s' to %s()"), "shape", __func__);

	int packed_ = 0;
	if (shape_ != 'g') {
		if (TYPEOF(packed) != LGLSXP || LENGTH(packed) < 1 ||
		    (packed_ = LOGICAL(packed)[0]) == NA_LOGICAL)
			error(_("invalid '%s' to %s()"), "packed", __func__);
	}

	char ul = 'U';
	if (shape_ != 'g') {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid '%s' to %s()"), "uplo", __func__);
	}

	return diagonal_as_dense(from, valid[ivalid], shape_, packed_, ul);
}

SEXP index_as_dense(SEXP from, const char *class, char kind)
{
	SEXP margin = PROTECT(GET_SLOT(from, Matrix_marginSym));
	int mg = INTEGER(margin)[0] - 1;
	UNPROTECT(1); /* margin */

	char cl[] = ".geMatrix";
	cl[0] = (kind != '.') ? kind : 'n';
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	Matrix_int_fast64_t len = (Matrix_int_fast64_t) m * n;
	if (len > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	double bytes = (double) len * kind2size(cl[0]);
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

	SEXP x = PROTECT(allocVector(kind2type(cl[0]), (R_xlen_t) len));
	SET_SLOT(to, Matrix_xSym, x);

#define IAD_SUBCASES(_CTYPE_, _PTR_, _ONE_) \
	do { \
		_CTYPE_ *px = _PTR_(x); \
		Matrix_memset(px, 0, (R_xlen_t) len, sizeof(_CTYPE_)); \
		if (mg == 0) { \
			R_xlen_t m1 = (R_xlen_t) m; \
			for (int i = 0; i < m; ++i) \
				px[i + m1 * (*(pperm++) - 1)] = _ONE_; \
		} else { \
			for (int j = 0; j < n; ++j, px += m) \
				px[          *(pperm++) - 1 ] = _ONE_; \
		} \
	} while (0)

	switch (cl[0]) {
	case 'n':
	case 'l':
		IAD_SUBCASES(int, LOGICAL, 1);
		break;
	case 'i':
		IAD_SUBCASES(int, INTEGER, 1);
		break;
	case 'd':
		IAD_SUBCASES(double, REAL, 1.0);
		break;
	case 'z':
		IAD_SUBCASES(Rcomplex, COMPLEX, Matrix_zone);
		break;
	default:
		break;
	}

#undef IAD_SUBCASES

	UNPROTECT(3); /* x, perm, to */
	return to;
}

/* as(<indMatrix>, ".geMatrix") */
SEXP R_index_as_dense(SEXP from, SEXP kind)
{
	static const char *valid[] = { "indMatrix", "pMatrix" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	return index_as_dense(from, valid[ivalid], kind_);
}

SEXP matrix_as_sparse(SEXP from, const char *zzz, char ul, char di,
                      int transpose_if_vector)
{
	char cl[] = "...Matrix";
	cl[0] = type2kind(TYPEOF(from));
	cl[1] = zzz[1];
	cl[2] = (zzz[1] == 'g') ? 'e' : ((zzz[1] == 's') ? 'y' : 'r');
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);
	REPROTECT(from = matrix_as_dense(from, cl, ul, di, transpose_if_vector, 1), pid);
	REPROTECT(from = dense_as_sparse(from, cl, zzz[2]), pid);
	cl[2] = zzz[2];
	REPROTECT(from = sparse_as_kind(from, cl, zzz[0]), pid);
	UNPROTECT(1);
	return from;
}

/* as(<matrix|vector>, ".[gst][CRT]Matrix") */
SEXP R_matrix_as_sparse(SEXP from, SEXP class, SEXP uplo, SEXP diag)
{
	switch (TYPEOF(from)) {
	case LGLSXP:
#ifdef MATRIX_ENABLE_IMATRIX
	case INTSXP:
#endif
	case REALSXP:
#ifdef MATRIX_ENABLE_ZMATRIX
	case CPLXSXP:
#endif
		break;
#ifndef MATRIX_ENABLE_IMATRIX
	case INTSXP:
		if (!inherits(from, "factor"))
			break;
#endif
	default:
		ERROR_INVALID_CLASS(from, __func__);
		break;
	}

	const char *zzz;
	if (TYPEOF(class) != STRSXP || LENGTH(class) < 1 ||
	    (class = STRING_ELT(class, 0)) == NA_STRING ||
	    (zzz = CHAR(class))[0] == '\0' ||
	    (zzz[1] != 'g' && zzz[1] != 's' && zzz[1] != 't') ||
	    (zzz[2] != 'C' && zzz[2] != 'R' && zzz[2] != 'T'))
		error(_("invalid '%s' to %s()"), "class", __func__);

	char ul = 'U', di = 'N';
	if (zzz[1] != 'g') {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid '%s' to %s()"), "uplo", __func__);
		if (zzz[1] == 't') {
			if (TYPEOF(diag) != STRSXP || LENGTH(diag) < 1 ||
			    (diag = STRING_ELT(diag, 0)) == NA_STRING ||
			    ((di = *CHAR(diag)) != 'N' && di != 'U'))
				error(_("invalid '%s' to %s()"), "diag", __func__);
		}
	}

	return matrix_as_sparse(from, zzz, ul, di, 0);
}

SEXP dense_as_sparse(SEXP from, const char *class, char repr)
{
	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = class[1];
	cl[2] = repr;
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	if (m != n || n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	char ul = 'U', di = 'N';
	if (class[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		p1, i1, j1;
	int i, j, *pp, *pi, *pj;
	R_xlen_t nnz = 0;
	p1 = i1 = j1 = NULL;
	pp = pi = pj = NULL;

#define DAS_CASES(_MASK_) \
	do { \
		switch (class[0]) { \
		case 'l': \
			DAS_SUBCASES(int, LOGICAL, _MASK_, ISNZ_LOGICAL); \
			break; \
		case 'i': \
			DAS_SUBCASES(int, INTEGER, _MASK_, ISNZ_INTEGER); \
			break; \
		case 'd': \
			DAS_SUBCASES(double, REAL, _MASK_, ISNZ_REAL); \
			break; \
		case 'z': \
			DAS_SUBCASES(Rcomplex, COMPLEX, _MASK_, ISNZ_COMPLEX); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define DAS_SUBCASES(_CTYPE_, _PTR_, _MASK_, _ISNZ_) \
	do { \
		       _CTYPE_ *px0 = _PTR_(x0) ; \
		_MASK_(_CTYPE_ *px1 = _PTR_(x1)); \
		if (class[1] == 'g') \
			/* .geMatrix */ \
			DAS_SUBSUBCASES(DAS_LOOP_GEN2C, DAS_LOOP_GEN2R, DAS_LOOP_GEN2C, \
			                _MASK_, _ISNZ_);	 \
		else if (class[2] != 'p' && di == 'N') \
			/* .syMatrix, non-unit diagonal .trMatrix */ \
			DAS_SUBSUBCASES(DAS_LOOP_TRN2C, DAS_LOOP_TRN2R, DAS_LOOP_TRN2C, \
			                _MASK_, _ISNZ_); \
		else if (class[2] != 'p') \
			/* unit diagonal .trMatrix */ \
			DAS_SUBSUBCASES(DAS_LOOP_TRU2C, DAS_LOOP_TRU2R, DAS_LOOP_TRU2C, \
			                _MASK_, _ISNZ_); \
		else if (di == 'N') \
			/* .spMatrix, non-unit diagonal .tpMatrix */ \
			DAS_SUBSUBCASES(DAS_LOOP_TPN2C, DAS_LOOP_TPN2R, DAS_LOOP_TPN2C, \
			                _MASK_, _ISNZ_); \
		else \
			/* unit diagonal .tpMatrix */ \
			DAS_SUBSUBCASES(DAS_LOOP_TPU2C, DAS_LOOP_TPU2R, DAS_LOOP_TPU2C, \
			                _MASK_, _ISNZ_); \
	} while (0)

#undef DAS_SUBSUBCASES
#define DAS_SUBSUBCASES(_LOOP_C_, _LOOP_R_, _LOOP_T_, _MASK_, _ISNZ_) \
	do { \
		switch (cl[2]) { \
		case 'C': \
			_LOOP_C_(_ISNZ_, ++nnz, DAS_VALID2C); \
			break; \
		case 'R': \
			_LOOP_R_(_ISNZ_, ++nnz, DAS_VALID2R); \
			break; \
		case 'T': \
			_LOOP_T_(_ISNZ_, ++nnz, DAS_VALID2T); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define DAS_LOOP_GEN2C(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < m; ++i, ++px0) \
				if (_ISNZ_(*px0)) _DO_INNER_; \
			_DO_OUTER_; \
		} \
	} while (0)

#define DAS_LOOP_GEN2R(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		R_xlen_t mn1s = (R_xlen_t) m * n - 1; \
		for (i = 0; i < m; ++i, px0 -= mn1s) { \
			for (j = 0; j < n; ++j, px0 += m) \
				if (_ISNZ_(*px0)) _DO_INNER_; \
			_DO_OUTER_; \
		} \
	} while (0)

#define DAS_LOOP_TRN2C(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		if (ul == 'U') { \
			for (j = 0; j < n; px0 += n - (++j)) { \
				for (i = 0; i <= j; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} else { \
			for (j = 0; j < n; px0 += (++j)) { \
				for (i = j; i < n; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} \
	} while (0)

#define DAS_LOOP_TRN2R(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		R_xlen_t d; \
		if (ul == 'U') { \
			d = (R_xlen_t) n * n - 1; \
			for (i = 0; i < n; ++i, px0 -= (d -= n)) { \
				for (j = i; j < n; ++j, px0 += n) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} else { \
			d = -1; \
			for (i = 0; i < n; ++i, px0 -= (d += n)) { \
				for (j = 0; j <= i; ++j, px0 += n) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} \
	} while (0)

#define DAS_LOOP_TRU2C(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		if (ul == 'U') { \
			px0 += n; \
			for (j = 1; j < n; ++j) { \
				for (i = 0; i < j; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
				px0 += n - j; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				px0 += j + 1; \
				for (i = j + 1; i < n; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} \
	} while (0)

#define DAS_LOOP_TRU2R(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		R_xlen_t d; \
		if (ul == 'U') { \
			d = (R_xlen_t) n * (n - 1) - 1; \
			for (i = 0; i < n; ++i) { \
				for (j = i + 1; j < n; ++j) { \
					px0 += n; \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				} \
				_DO_OUTER_; \
				px0 -= (d -= n); \
			} \
		} else { \
			++px0; \
			d = -1; \
			for (i = 1; i < n; ++i) { \
				for (j = 0; j < i; ++j) { \
					if (_ISNZ_(*px0)) _DO_INNER_; \
					px0 += n; \
				} \
				_DO_OUTER_; \
				px0 -= (d += n); \
			} \
		} \
	} while (0)

#define DAS_LOOP_TPN2C(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i <= j; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				for (i = j; i < n; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} \
	} while (0)

#define DAS_LOOP_TPN2R(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		R_xlen_t d; \
		if (ul == 'U') { \
			d = PM_LENGTH(n) - 1; \
			for (i = 0; i < n; px0 -= (d -= (++i))) { \
				for (j = i; j < n; px0 += (++j)) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} else { \
			d = -1; \
			for (i = 0; i < n; px0 -= (d += n - (++i))) { \
				for (j = 0; j <= i; px0 += n - (++j)) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} \
	} while (0)

#define DAS_LOOP_TPU2C(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		if (ul == 'U') { \
			for (j = 1; j < n; ++j) { \
				++px0; \
				for (i = 0; i < j; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				++px0; \
				for (i = j + 1; i < n; ++i, ++px0) \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				_DO_OUTER_; \
			} \
		} \
	} while (0)

#define DAS_LOOP_TPU2R(_ISNZ_, _DO_INNER_, _DO_OUTER_) \
	do { \
		R_xlen_t d; \
		if (ul == 'U') { \
			d = PM_LENGTH(n - 1) - 1; \
			for (i = 0; i < n; ++i) { \
				for (j = i + 1; j < n; ++j) { \
					px0 += j; \
					if (_ISNZ_(*px0)) _DO_INNER_; \
				} \
				_DO_OUTER_; \
				px0 -= (d -= i + 1); \
			} \
		} else { \
			++px0; \
			d = -1; \
			for (i = 1; i < n; ++i) { \
				for (j = 0; j < i; ++j) { \
					if (_ISNZ_(*px0)) _DO_INNER_; \
					px0 += n - j - 1; \
				} \
				_DO_OUTER_; \
				px0 -= (d += n - i); \
			} \
		} \
	} while (0)

#define DAS_VALID2T \
	if (nnz > INT_MAX) \
		error(_("attempt to construct %s with more than %s nonzero entries"), \
		      "sparseMatrix", "2^31-1")

#define DAS_VALID2C \
	do { DAS_VALID2T; else *(pp++) = (int) nnz; } while (0)

#define DAS_VALID2R DAS_VALID2C

	/* First we loop over the _nontrivial part_ of the denseMatrix 'from',
	   by row ('R' case) or by column ('C' and 'T' cases), counting the
	   nonzero entries and filling the 'p' slot of the result accordingly
	   ('C' and 'R' cases) ...
	*/

	int nprotect = 2;

	if (cl[2] != 'T') {
		int r = (cl[2] == 'C') ? n : m;
		PROTECT(p1 = allocVector(INTSXP, (R_xlen_t) r + 1));
		++nprotect;
		SET_SLOT(to, Matrix_pSym, p1);
		pp = INTEGER(p1);
		*(pp++) = 0;
		if (r > 0 && di != 'N' && ul == ((cl[2] == 'C') ? 'U' : 'L'))
			*(pp++) = 0; /* first row or column skipped in these loops */
	}
	if (class[0] == 'n')
		DAS_SUBCASES(int, LOGICAL, HIDE, ISNZ_LOGICAL);
	else
		DAS_CASES(HIDE);
	if (cl[2] != 'R') {
		PROTECT(i1 = allocVector(INTSXP, nnz));
		++nprotect;
		SET_SLOT(to, Matrix_iSym, i1);
		pi = INTEGER(i1);
	}
	if (cl[2] != 'C') {
		PROTECT(j1 = allocVector(INTSXP, nnz));
		++nprotect;
		SET_SLOT(to, Matrix_jSym, j1);
		pj = INTEGER(j1);
	}

#undef DAS_SUBSUBCASES
#define DAS_SUBSUBCASES(_LOOP_C_, _LOOP_R_, _LOOP_T_, _MASK_, _ISNZ_) \
	do { \
		switch (repr) { \
		case 'C': \
			_LOOP_C_(_ISNZ_, \
			         do { \
			             *(pi++) = i; \
			             _MASK_(*(px1++) = *px0); \
			         } while (0), ); \
			break; \
		case 'R': \
			_LOOP_R_(_ISNZ_, \
			         do { \
			             *(pj++) = j; \
			             _MASK_(*(px1++) = *px0); \
			         } while (0), ); \
			break; \
		case 'T': \
			_LOOP_T_(_ISNZ_, \
			         do { \
			             *(pi++) = i; \
			             *(pj++) = j; \
			             _MASK_(*(px1++) = *px0); \
			         } while (0), ); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	/* Then we loop back over the same entries in order to fill
	   the 'i', 'j', and 'x' slots of the result (whichever exist) ...
	*/

	if (class[0] == 'n')
		DAS_SUBCASES(int, LOGICAL, HIDE, ISNZ_LOGICAL);
	else if (cl[2] != 'R' && nnz == XLENGTH(x0)) {
		SET_SLOT(to, Matrix_xSym, x0);
		DAS_CASES(HIDE);
	} else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), nnz));
		SET_SLOT(to, Matrix_xSym, x1);
		DAS_CASES(SHOW);
		UNPROTECT(1); /* x1 */
	}

#undef DAS_CASES
#undef DAS_SUBCASES
#undef DAS_SUBSUBCASES
#undef DAS_VALID2C
#undef DAS_VALID2R
#undef DAS_VALID2T
#undef DAS_LOOP_GEN2C
#undef DAS_LOOP_GEN2R
#undef DAS_LOOP_TRN2C
#undef DAS_LOOP_TRN2R
#undef DAS_LOOP_TRU2C
#undef DAS_LOOP_TRU2R
#undef DAS_LOOP_TPN2C
#undef DAS_LOOP_TPN2R
#undef DAS_LOOP_TPU2C
#undef DAS_LOOP_TPU2R

	UNPROTECT(nprotect);
	return to;
}

/* as(<denseMatrix>, "[CRT]sparseMatrix") */
SEXP R_dense_as_sparse(SEXP from, SEXP repr)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char repr_;
	if (TYPEOF(repr) != STRSXP || LENGTH(repr) < 1 ||
	    (repr = STRING_ELT(repr, 0)) == NA_STRING ||
	    ((repr_ = CHAR(repr)[0]) != 'C' && repr_ != 'R' && repr_ != 'T'))
		error(_("invalid '%s' to %s()"), "repr", __func__);

	return dense_as_sparse(from, valid[ivalid], repr_);
}

SEXP diagonal_as_sparse(SEXP from, const char *class,
                           char shape, char repr, char ul)
{
	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = shape;
	cl[2] = repr;
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (cl[1] == 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (cl[1] != 'g' && ul != 'U') {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (cl[1] == 't' && di != 'N')
		SET_SLOT(to, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	if (cl[1] == 't' && di != 'N') {
		if (cl[2] != 'T') {
			SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
			SET_SLOT(to, Matrix_pSym, p);
			Matrix_memset(INTEGER(p), 0, (R_xlen_t) n + 1, sizeof(int));
			UNPROTECT(1); /* p */
		}
		UNPROTECT(1); /* to */
		return to;
	}

#define DAS_CASES(_MASK_) \
	do { \
		switch (class[0]) { \
		case 'l': \
			DAS_LOOP(int, LOGICAL, _MASK_, ISNZ_LOGICAL, 1); \
			break; \
		case 'i': \
			DAS_LOOP(int, INTEGER, _MASK_, ISNZ_INTEGER, 1); \
			break; \
		case 'd': \
			DAS_LOOP(double, REAL, _MASK_, ISNZ_REAL, 1.0); \
			break; \
		case 'z': \
			DAS_LOOP(Rcomplex, COMPLEX, _MASK_, ISNZ_COMPLEX, Matrix_zone); \
			break; \
		default: \
			break; \
		} \
	} while (0)

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
	int d, nnz;
	if (cl[2] != 'T') {
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1));
		SET_SLOT(to, Matrix_pSym, p);
		int *pp = INTEGER(p);
		*(pp++) = 0;
		if (di == 'N') {
			nnz = 0;

#undef DAS_LOOP
#define DAS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_, _ONE_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0); \
				for (d = 0; d < n; ++d) { \
					if (_ISNZ_(*px0)) \
						++nnz; \
					*(pp++) = nnz; \
					++px0; \
				} \
			} while (0)

			if (class[0] == 'n')
				DAS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL, 1);
			else
				DAS_CASES(SHOW);
		} else {
			nnz = n;
			for (d = 1; d <= n; ++d)
				*(pp++) = d;
		}
		UNPROTECT(1); /* p */
	} else {
		if (di == 'N') {
			nnz = 0;

#undef DAS_LOOP
#define DAS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_, _ONE_) \
			do { \
				_CTYPE_ *px0 = _PTR_(x0); \
				for (d = 0; d < n; ++d) { \
					if (_ISNZ_(*px0)) \
						++nnz; \
					++px0; \
				} \
			} while (0)

			if (class[0] == 'n')
				DAS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL, 1);
			else
				DAS_CASES(SHOW);
		} else
			nnz = n;
	}

	SEXP i1 = PROTECT(allocVector(INTSXP, nnz));
	if (cl[2] != 'T')
		SET_SLOT(to, (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym, i1);
	else {
		SET_SLOT(to, Matrix_iSym, i1);
		SET_SLOT(to, Matrix_jSym, i1);
	}
	int *pi1 = INTEGER(i1);

#undef DAS_LOOP
#define DAS_LOOP(_CTYPE_, _PTR_, _MASK_, _ISNZ_, _ONE_) \
	do { \
		_MASK_(_CTYPE_ *px1 = _PTR_(x1)); \
		if (di == 'N') { \
			_CTYPE_ *px0 = _PTR_(x0); \
			for (d = 0; d < n; ++d) { \
				if (_ISNZ_(*px0)) { \
					*(pi1++) = d; \
					_MASK_(*(px1++) = *px0); \
				} \
				++px0; \
			} \
		} else { \
			for (d = 0; d < n; ++d) { \
				*(pi1++) = d; \
				_MASK_(*(px1++) = _ONE_); \
			} \
		} \
	} while (0)

	if (class[0] == 'n')
		DAS_LOOP(int, LOGICAL, HIDE, ISNZ_LOGICAL, 1);
	else if (di == 'N' && nnz == n) {
		SET_SLOT(to, Matrix_xSym, x0);
		DAS_CASES(HIDE);
	} else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), nnz));
		SET_SLOT(to, Matrix_xSym, x1);
		DAS_CASES(SHOW);
		UNPROTECT(1); /* x1 */
	}

	UNPROTECT(3); /* i1, x0, to */
	return to;
}

/* as(<diagonalMatrix>, ".[gst][CRT]Matrix") */
SEXP R_diagonal_as_sparse(SEXP from, SEXP shape, SEXP repr, SEXP uplo)
{
	static const char *valid[] = { VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char shape_;
	if (TYPEOF(shape) != STRSXP || LENGTH(shape) < 1 ||
	    (shape = STRING_ELT(shape, 0)) == NA_STRING ||
	    ((shape_ = CHAR(shape)[0]) != 'g' && shape_ != 's' && shape_ != 't'))
		error(_("invalid '%s' to %s()"), "shape", __func__);

	char repr_;
	if (TYPEOF(repr) != STRSXP || LENGTH(repr) < 1 ||
	    (repr = STRING_ELT(repr, 0)) == NA_STRING ||
	    ((repr_ = CHAR(repr)[0]) != 'C' && repr_ != 'R' && repr_ != 'T'))
		error(_("invalid '%s' to %s()"), "repr", __func__);

	char ul = 'U';
	if (shape_ != 'g') {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid '%s' to %s()"), "uplo", __func__);
	}

	return diagonal_as_sparse(from, valid[ivalid], shape_, repr_, ul);
}

SEXP index_as_sparse(SEXP from, const char *class, char kind, char repr)
{
	SEXP margin = PROTECT(GET_SLOT(from, Matrix_marginSym));
	int mg = INTEGER(margin)[0] - 1;
	UNPROTECT(1); /* margin */

	char cl[] = ".g.Matrix";
	cl[0] = (kind != '.') ? kind : 'n';
	cl[2] = (repr != '.') ? repr : ((mg == 0) ? 'R' : 'C');
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
	int *pperm = INTEGER(perm);

	if (cl[2] == 'T') {
		SEXP i = PROTECT(allocVector(INTSXP, r)),
			j = PROTECT(allocVector(INTSXP, r));
		int k, *pi = INTEGER(i), *pj = INTEGER(j);
		for (k = 0; k < r; ++k) {
			*(pi++) = k;
			*(pj++) = *(pperm++) - 1;
		}
		SET_SLOT(to, Matrix_iSym, (mg == 0) ? i : j);
		SET_SLOT(to, Matrix_jSym, (mg == 0) ? j : i);
	} else if ((cl[2] == 'C') == (mg != 0)) {
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) r + 1)),
			i = PROTECT(allocVector(INTSXP, r));
		int k, *pp = INTEGER(p), *pi = INTEGER(i);
		for (k = 0; k < r; ++k) {
			*(pp++) = k;
			*(pi++) = *(pperm++) - 1;
		}
		*pp = r;
		SET_SLOT(to, Matrix_pSym, p);
		SET_SLOT(to, (mg != 0) ? Matrix_iSym : Matrix_jSym, i);
	} else {
		SEXP p = PROTECT(allocVector(INTSXP, (R_xlen_t) s + 1));
		int k, *pp = INTEGER(p);
		Matrix_memset(pp, 0, (R_xlen_t) s + 1, sizeof(int));
		for (k = 0; k < r; ++k)
			++pp[pperm[k]];
		for (k = 0; k < s; ++k)
			pp[k + 1] += pp[k];
		SEXP j = PROTECT(allocVector(INTSXP, r));
		int *pj = INTEGER(j), *work;
		Matrix_Calloc(work, s, int);
		Matrix_memcpy(work, pp, s, sizeof(int));
		--work;
		for (k = 0; k < r; ++k)
			pj[work[pperm[k]]++] = k;
		++work;
		Matrix_Free(work, s);
		SET_SLOT(to, Matrix_pSym, p);
		SET_SLOT(to, (mg != 0) ? Matrix_jSym : Matrix_iSym, j);
	}
	UNPROTECT(2);

	if (cl[0] != 'n') {
		SEXP x = PROTECT(allocVector(kind2type(cl[0]), r));
		SET_SLOT(to, Matrix_xSym, x);

#define IAS_SUBCASES(_CTYPE_, _PTR_, _ONE_) \
		do { \
			_CTYPE_ *px = _PTR_(x); \
			for (int k = 0; k < r; ++k) \
				*(px++) = _ONE_; \
		} while (0)

		switch (cl[0]) {
		case 'l':
			IAS_SUBCASES(int, LOGICAL, 1);
			break;
		case 'i':
			IAS_SUBCASES(int, INTEGER, 1);
			break;
		case 'd':
			IAS_SUBCASES(double, REAL, 1.0);
			break;
		case 'z':
			IAS_SUBCASES(Rcomplex, COMPLEX, Matrix_zone);
			break;
		default:
			break;
		}
		UNPROTECT(1); /* x */
	}

#undef IAS_SUBCASES

	UNPROTECT(2); /* perm, to */
	return to;
}

/* as(<indMatrix>, ".g[CRT]Matrix") */
SEXP R_index_as_sparse(SEXP from, SEXP kind, SEXP repr)
{
	static const char *valid[] = { "indMatrix", "pMatrix" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	char repr_;
	if (TYPEOF(repr) != STRSXP || LENGTH(repr) < 1 ||
	    (repr = STRING_ELT(repr, 0)) == NA_STRING ||
	    ((repr_ = CHAR(repr)[0]) != '.' &&
	     repr_ != 'C' && repr_ != 'R' && repr_ != 'T'))
		error(_("invalid '%s' to %s()"), "repr", __func__);

	return index_as_sparse(from, valid[ivalid], kind_, repr_);
}

SEXP dense_as_kind(SEXP from, const char *class, char kind)
{
	if (kind == '.' || kind == class[0])
		return from;
	SEXPTYPE tt = kind2type(kind);

	char cl[] = "...Matrix";
	cl[0] = kind;
	cl[1] = class[1];
	cl[2] = class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	}

	PROTECT_INDEX pid;
	SEXP x;
	PROTECT_WITH_INDEX(x = GET_SLOT(from, Matrix_xSym), &pid);
	if (TYPEOF(x) == tt) {
		if (class[0] == 'n') {
			/* n->l */
			R_xlen_t i, len = XLENGTH(x);
			int *px = LOGICAL(x);
			for (i = 0; i < len; ++i, ++px) {
				if (*px == NA_LOGICAL) {
					REPROTECT(x = duplicate(x), pid);
					na2one(x);
					break;
				}
			}
		}
	} else {
		REPROTECT(x = coerceVector(x, tt), pid);
		if (class[0] == 'n')
			/* n->[idz] */
			na2one(x);
	}
	SET_SLOT(to, Matrix_xSym, x);
	UNPROTECT(2); /* x, to */
	return to;
}

/* as(<denseMatrix>, "[nlidz]Matrix") */
SEXP R_dense_as_kind(SEXP from, SEXP kind)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	return dense_as_kind(from, valid[ivalid], kind_);
}

SEXP sparse_as_kind(SEXP from, const char *class, char kind)
{
	if (kind == '.' || kind == class[0])
		return from;
	SEXPTYPE tt = kind2type(kind);

	if (class[2] == 'T' && (class[0] == 'n' || class[0] == 'l') &&
	    kind != 'n' && kind != 'l') {
		/* defined in ./sparse.c : */
		SEXP Tsparse_aggregate(SEXP);
		from = Tsparse_aggregate(from);
	}
	PROTECT(from);

	char cl[] = "...Matrix";
	cl[0] = kind;
	cl[1] = class[1];
	cl[2] = class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
		if (class[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			char di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	}

	R_xlen_t nnz = -1;
	if (class[2] != 'T') {
		SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
		SET_SLOT(to, Matrix_pSym, p);
		nnz = INTEGER(p)[XLENGTH(p) - 1];
		UNPROTECT(1); /* p */
	}
	if (class[2] != 'R') {
		SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
		SET_SLOT(to, Matrix_iSym, i);
		if (nnz < 1)
		nnz = XLENGTH(i);
		UNPROTECT(1); /* i */
	}
	if (class[2] != 'C') {
		SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_jSym, j);
		UNPROTECT(1); /* j */
	}
	if (class[0] == 'n') {
		SEXP x = PROTECT(allocVector(tt, nnz));
		switch (tt) {
		case LGLSXP:
		{
			int *px = LOGICAL(x);
			while (nnz--)
				*(px++) = 1;
			break;
		}
		case INTSXP:
		{
			int *px = INTEGER(x);
			while (nnz--)
				*(px++) = 1;
			break;
		}
		case REALSXP:
		{
			double *px = REAL(x);
			while (nnz--)
				*(px++) = 1.0;
			break;
		}
		case CPLXSXP:
		{
			Rcomplex *px = COMPLEX(x);
			while (nnz--)
				*(px++) = Matrix_zone;
			break;
		}
		default:
			break;
		}
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
SEXP R_sparse_as_kind(SEXP from, SEXP kind)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	return sparse_as_kind(from, valid[ivalid], kind_);
}

SEXP diagonal_as_kind(SEXP from, const char *class, char kind)
{
	if (kind == '.' || kind == class[0])
		return from;
	SEXPTYPE tt = kind2type(kind);

	char cl[] = ".diMatrix";
	cl[0] = kind;
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
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
				R_xlen_t i, len = XLENGTH(x);
				int *px = LOGICAL(x);
				for (i = 0; i < len; ++i, ++px) {
					if (*px == NA_LOGICAL) {
						REPROTECT(x = duplicate(x), pid);
						na2one(x);
						break;
					}
				}
			}
		} else {
			REPROTECT(x = coerceVector(x, tt), pid);
			if (class[0] == 'n')
				/* n->[idz] */
				na2one(x);
		}
		SET_SLOT(to, Matrix_xSym, x);
		UNPROTECT(1); /* x */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<diagonalMatrix>, "[nlidz]Matrix") */
SEXP R_diagonal_as_kind(SEXP from, SEXP kind)
{
	static const char *valid[] = { VALID_DIAGONAL, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	return diagonal_as_kind(from, valid[ivalid], kind_);
}

SEXP index_as_kind(SEXP from, const char *class, char kind)
{
	return index_as_sparse(from, class, kind, '.');
}

/* as(<indMatrix>, "[nlidz]Matrix") */
SEXP R_index_as_kind(SEXP from, SEXP kind)
{
	static const char *valid[] = { "indMatrix", "pMatrix" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	return index_as_kind(from, valid[ivalid], kind_);
}

SEXP dense_as_general(SEXP from, const char *class, int new)
{
	if (class[1] == 'g')
		return from;

	char cl[] = ".geMatrix";
	cl[0] = class[0];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] != 's')
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	else
		set_symmetrized_DimNames(to, dimnames, -1);
	UNPROTECT(1); /* dimnames */

	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0)), di = 'N';
	UNPROTECT(1); /* uplo */

	if (class[1] == 's') {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */
	} else {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */
	}

	if ((Matrix_int_fast64_t) n * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1 = x0;
	int nprotect = 2;
	if (class[2] == 'p' || new > 0) {
		PROTECT(x1 = allocVector(TYPEOF(x0), (R_xlen_t) n * n));
		++nprotect;
	}
	SET_SLOT(to, Matrix_xSym, x1);

#define DAG_SUBCASES(_CTYPE_, _PTR_, _PREFIX_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (class[2] == 'p') \
			_PREFIX_ ## dense_unpack(px1, px0, n, ul, di); \
		else if (new > 0) \
			Matrix_memcpy(px1, px0, (R_xlen_t) n * n, sizeof(_CTYPE_)); \
		if (class[1] == 's') \
			_PREFIX_ ## dense_unpacked_make_symmetric(px1, n, ul); \
		else \
			_PREFIX_ ## dense_unpacked_make_triangular(px1, n, n, ul, di); \
	} while (0)

	switch (class[0]) {
	case 'n':
	case 'l':
		DAG_SUBCASES(int, LOGICAL, i);
		break;
	case 'i':
		DAG_SUBCASES(int, INTEGER, i);
		break;
	case 'd':
		DAG_SUBCASES(double, REAL, d);
		break;
	case 'z':
		DAG_SUBCASES(Rcomplex, COMPLEX, z);
		break;
	default:
		break;
	}

#undef DAG_CASES
#undef DAG_SUBCASES

	UNPROTECT(nprotect);
	return to;
}

/* as(<denseMatrix>, "generaMatrix") */
SEXP R_dense_as_general(SEXP from)
{
	static const char *valid[] = {
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return dense_as_general(from, valid[ivalid], 1);
}

SEXP sparse_as_general(SEXP from, const char *class)
{
	if (class[1] == 'g')
		return from;

	char cl[] = ".g.Matrix";
	cl[0] = class[0];
	cl[2] = class[2];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (class[1] == 's')
		set_symmetrized_DimNames(to, dimnames, -1);
	else
		SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	if (class[1] == 's') {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */
	} else {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		UNPROTECT(1); /* diag */

		if (di == 'N') {
			if (class[2] != 'T') {
				SEXP p = PROTECT(GET_SLOT(from, Matrix_pSym));
				SET_SLOT(to, Matrix_pSym, p);
				UNPROTECT(1); /* p */
			}
			if (class[2] != 'R') {
				SEXP i = PROTECT(GET_SLOT(from, Matrix_iSym));
				SET_SLOT(to, Matrix_iSym, i);
				UNPROTECT(1); /* i */
			}
			if (class[2] != 'C') {
				SEXP j = PROTECT(GET_SLOT(from, Matrix_jSym));
				SET_SLOT(to, Matrix_jSym, j);
				UNPROTECT(1); /* j */
			}
			if (class[0] != 'n') {
				SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
				SET_SLOT(to, Matrix_xSym, x);
				UNPROTECT(1); /* x */
			}
			UNPROTECT(1); /* to */
			return to;
		}
	}

	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	UNPROTECT(1); /* uplo */

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) n + 1)),
			i0 = PROTECT(GET_SLOT(from, iSym));
		int j, k, kend,
			*pp0 = INTEGER(p0), *pp1 = INTEGER(p1), *pi0 = INTEGER(i0);
		SET_SLOT(to, Matrix_pSym, p1);
		pp0++; *(pp1++) = 0;

		if (class[1] == 's') {
			Matrix_memset(pp1, 0, n, sizeof(int));
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

#define SAG_CASES \
		do { \
			switch (class[0]) { \
			case 'l': \
				SAG_SUBCASES(int, LOGICAL, SHOW, 1); \
				break; \
			case 'i': \
				SAG_SUBCASES(int, INTEGER, SHOW, 1); \
				break; \
			case 'd': \
				SAG_SUBCASES(double, REAL, SHOW, 1.0); \
				break; \
			case 'z': \
				SAG_SUBCASES(Rcomplex, COMPLEX, SHOW, Matrix_zone); \
				break; \
			default: \
				break; \
			} \
		} while (0)

#undef SAG_SUBCASES
#define SAG_SUBCASES(_CTYPE_, _PTR_, _MASK_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			if (class[1] == 's') { \
				int *pp1_; \
				Matrix_Calloc(pp1_, n, int); \
				Matrix_memcpy(pp1_, pp1 - 1, n, sizeof(int)); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						pi1[pp1_[j]] = pi0[k]; \
						_MASK_(px1[pp1_[j]] = px0[k]); \
						++pp1_[j]; \
						if (pi0[k] != j) { \
							pi1[pp1_[pi0[k]]] = j; \
							_MASK_(px1[pp1_[pi0[k]]] = px0[k]); \
							++pp1_[pi0[k]]; \
						} \
						++k; \
					} \
				} \
				Matrix_Free(pp1_, n); \
			} else if (ul == ((class[2] == 'C') ? 'U' : 'L')) { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						*(pi1++) = *(pi0++); \
						_MASK_(*(px1++) = *(px0++)); \
						++k; \
					} \
					*(pi1++) = j; \
					_MASK_(*(px1++) = _ONE_); \
				} \
			} else { \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					*(pi1++) = j; \
					_MASK_(*(px1++) = _ONE_); \
					while (k < kend) { \
						*(pi1++) = *(pi0++); \
						_MASK_(*(px1++) = *(px0++)); \
						++k; \
					} \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			SAG_SUBCASES(int, LOGICAL, HIDE, 1);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), pp1[n - 1]));
			SET_SLOT(to, Matrix_xSym, x1);
			SAG_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t nnz0 = XLENGTH(i0), nnz1;

		if (class[1] == 's') {
			nnz1 = nnz0;
			for (R_xlen_t k = 0; k < nnz0; ++k)
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
		Matrix_memcpy(pi1, pi0, nnz0, sizeof(int));
		Matrix_memcpy(pj1, pj0, nnz0, sizeof(int));
		pi1 += nnz0;
		pj1 += nnz0;

#undef SAG_SUBCASES
#define SAG_SUBCASES(_CTYPE_, _PTR_, _MASK_, _ONE_) \
		do { \
			_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
			_MASK_(Matrix_memcpy(px1, px0, nnz0, sizeof(_CTYPE_))); \
			_MASK_(px1 += nnz0); \
			if (class[1] == 's') { \
				for (R_xlen_t k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						*(pi1++) = *pj0; \
						*(pj1++) = *pi0; \
						_MASK_(*(px1++) = *px0) ; \
					} \
					++pi0; ++pj0; _MASK_(++px0); \
				} \
			} else { \
				for (int d = 0; d < n; ++d) { \
					*(pi1++) = *(pj1++) = d; \
					_MASK_(*(px1++) = _ONE_); \
				} \
			} \
		} while (0)

		if (class[0] == 'n')
			SAG_SUBCASES(int, LOGICAL, HIDE, 1);
		else {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
				x1 = PROTECT(allocVector(TYPEOF(x0), nnz1));
			SET_SLOT(to, Matrix_xSym, x1);
			SAG_CASES;
			UNPROTECT(2); /* x1, x0 */
		}

	}

#undef SAG_CASES
#undef SAG_SUBCASES

	UNPROTECT(5);
	return to;
}

/* as(<[CRT]sparseMatrix>, "generalMatrix") */
SEXP R_sparse_as_general(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_as_general(from, valid[ivalid]);
}

SEXP dense_as_unpacked(SEXP from, const char *class)
{
	if (class[0] != 'p' && class[2] != 'p')
		return from;

	char cl[] = "...Matrix";
	if (class[0] == 'p') {
		cl[0] = 'c'; cl[1] = 'o'; cl[2] = 'r';
	} else if (class[1] == 'p') {
		cl[0] = 'd'; cl[1] = 'p'; cl[2] = 'o';
	} else {
		cl[0] = class[0];
		cl[1] = class[1];
		cl[2] = (class[1] == 's') ? 'y' : 'r';
	}
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	int n = INTEGER(dim)[0];
	if ((Matrix_int_fast64_t) n * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	if (n > 0)
		SET_SLOT(to, Matrix_DimSym, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	SET_SLOT(to, Matrix_DimNamesSym, dimnames);
	UNPROTECT(1); /* dimnames */

	SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	if (ul != 'U')
		SET_SLOT(to, Matrix_uploSym, uplo);
	UNPROTECT(1); /* uplo */

	if (cl[1] != 't') {
		SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
		if (LENGTH(factors) > 0)
			SET_SLOT(to, Matrix_factorsSym, factors);
		UNPROTECT(1); /* factors */

		if (cl[0] == 'c') {
			SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
			if (LENGTH(sd) > 0)
				SET_SLOT(to, Matrix_sdSym, sd);
			UNPROTECT(1); /* sd */
		}
	} else {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N')
			SET_SLOT(to, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) n * n));
	SET_SLOT(to, Matrix_xSym, x1);

#define UNPACK(_CTYPE_, _PTR_, _PREFIX_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		Matrix_memset(px1, 0, (R_xlen_t) n * n, sizeof(_CTYPE_)); \
		_PREFIX_ ## dense_unpack(px1, px0, n, ul, 'N'); \
	} while (0)

	switch (cl[0]) {
	case 'n':
	case 'l':
		UNPACK(int, LOGICAL, i);
		break;
	case 'i':
		UNPACK(int, INTEGER, i);
		break;
	case 'c':
	case 'd':
		UNPACK(double, REAL, d);
		break;
	case 'z':
		UNPACK(Rcomplex, COMPLEX, z);
		break;
	default:
		break;
	}

#undef UNPACK

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<denseMatrix>, "unpackedMatrix") */
SEXP R_dense_as_unpacked(SEXP from)
{
	static const char *valid[] = {
		"dpoMatrix", "dppMatrix", "corMatrix", "pcorMatrix",
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return dense_as_unpacked(from, valid[ivalid]);
}

SEXP dense_as_packed(SEXP from, const char *class, char ul, char di)
{
	if (class[0] == 'p' || class[2] == 'p')
		return from;

	char cl_[] = "p...Matrix", *cl = ((char *) &cl_) + 1;
	int ge = 0;
	if (class[0] == 'c') {
		cl[0] = 'c'; cl[1] = 'o'; cl[2] = 'r';
	} else if (class[1] == 'p') {
		cl[0] = 'd'; cl[1] = 'p'; cl[2] = 'p';
	} else {
		ge = class[1] == 'g';
		cl[0] = class[0];
		cl[1] = (!ge) ? class[1] : ((di == '\0') ? 's' : 't');
		cl[2] = 'p';
	}
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS((class[0] == 'c') ? cl - 1 : cl));

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

	if (ge) {
		if (ul != 'U') {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(to, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
		if (cl[1] == 't' && di != 'N') {
			SEXP diag = PROTECT(mkString("U"));
			SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}
	} else {
		SEXP uplo = PROTECT(GET_SLOT(from, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */

		if (cl[1] == 't') {
			SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
			di = *CHAR(STRING_ELT(diag, 0));
			if (di != 'N')
				SET_SLOT(to, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		} else {
			SEXP factors = PROTECT(GET_SLOT(from, Matrix_factorsSym));
			if (LENGTH(factors) > 0)
				SET_SLOT(to, Matrix_factorsSym, factors);
			UNPROTECT(1); /* factors */

			if (cl[0] == 'c') {
				SEXP sd = PROTECT(GET_SLOT(from, Matrix_sdSym));
				if (LENGTH(sd) > 0)
					SET_SLOT(to, Matrix_sdSym, sd);
				UNPROTECT(1); /* sd */
			}
		}
	}

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), PM_LENGTH(n)));
	SET_SLOT(to, Matrix_xSym, x1);

#define PACK(_CTYPE_, _PTR_, _PREFIX_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		_PREFIX_ ## dense_pack(px1, px0, n, ul, 'N'); \
	} while (0)

	switch (cl[0]) {
	case 'n':
	case 'l':
		PACK(int, LOGICAL, i);
		break;
	case 'i':
		PACK(int, INTEGER, i);
		break;
	case 'c':
	case 'd':
		PACK(double, REAL, d);
		break;
	case 'z':
		PACK(Rcomplex, COMPLEX, z);
		break;
	default:
		break;
	}

#undef PACK

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

/* as(<denseMatrix>, "packedMatrix") */
SEXP R_dense_as_packed(SEXP from, SEXP uplo, SEXP diag)
{
	static const char *valid[] = {
		"dpoMatrix", "dppMatrix", "corMatrix", "pcorMatrix",
		VALID_DDENSE, VALID_LDENSE, VALID_NDENSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	char ul = 'U', di = '\0';
	if (valid[ivalid][1] == 'g') {
		if (TYPEOF(uplo) != STRSXP || LENGTH(uplo) < 1 ||
		    (uplo = STRING_ELT(uplo, 0)) == NA_STRING ||
		    ((ul = *CHAR(uplo)) != 'U' && ul != 'L'))
			error(_("invalid '%s' to %s()"), "uplo", __func__);
		if (TYPEOF(diag) != STRSXP || LENGTH(diag) < 1 ||
		    ((diag = STRING_ELT(diag, 0)) != NA_STRING &&
		     (di = *CHAR(diag)) != '\0' && di != 'N' && di != 'U'))
			error(_("invalid '%s' to %s()"), "diag", __func__);
	}

	return dense_as_packed(from, valid[ivalid], ul, di);
}

void trans(SEXP p0, SEXP i0, SEXP x0, SEXP p1, SEXP i1, SEXP x1,
           int m, int n)
{
	int *pp0 = INTEGER(p0), *pp1 = INTEGER(p1),
		*pi0 = INTEGER(i0), *pi1 = INTEGER(i1),
		i, j, k, kend, nnz = pp0[n];
	Matrix_memset(pp1, 0, (R_xlen_t) m + 1, sizeof(int));
	++pp0; ++pp1;
	for (k = 0; k < nnz; ++k)
		++pp1[pi0[k]];
	for (i = 0; i < m; ++i)
		pp1[i] += pp1[i - 1];

	int *pp1_;
	Matrix_Calloc(pp1_, m, int);
	Matrix_memcpy(pp1_, pp1 - 1, m, sizeof(int));

#define TRANS_LOOP(_CTYPE_, _PTR_, _MASK_) \
	do { \
		_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1)); \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp0[j]; \
			while (k < kend) { \
				i = pi0[k]; \
				pi1[pp1_[i]] = j; \
				_MASK_(px1[pp1_[i]] = px0[k]); \
				++pp1_[i]; \
				++k; \
			} \
		} \
	} while (0)

	if (!(x0 && x1) || TYPEOF(x0) != TYPEOF(x1))
		TRANS_LOOP(int, LOGICAL, HIDE);
	else {
		switch (TYPEOF(x0)) {
		case LGLSXP:
			TRANS_LOOP(int, LOGICAL, SHOW);
			break;
		case INTSXP:
			TRANS_LOOP(int, INTEGER, SHOW);
			break;
		case REALSXP:
			TRANS_LOOP(double, REAL, SHOW);
			break;
		case CPLXSXP:
			TRANS_LOOP(Rcomplex, COMPLEX, SHOW);
			break;
		default:
			break;
		}
	}

#undef TRANS_LOOP

	Matrix_Free(pp1_, m);
	return;
}

void tsort(SEXP i0, SEXP j0, SEXP x0, SEXP *p1, SEXP *i1, SEXP *x1,
           int m, int n)
{
	R_xlen_t nnz0 = XLENGTH(i0), nnz1 = 0;
	if (nnz0 > INT_MAX)
		error(_("unable to aggregate %s with '%s' and '%s' slots of length exceeding %s"),
		      "TsparseMatrix", "i", "j", "2^31-1");

	/* FIXME: test for overflow and throw error only in that case */

	PROTECT(*p1 = allocVector(INTSXP, (R_xlen_t) n + 1));
	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), *pp1 = INTEGER(*p1), *pi1,
		i, j, r = (m < n) ? n : m;
	R_xlen_t k, kstart, kend, kend_;
	*(pp1++) = 0;

	int *workA, *workB, *workC, *pj_;
	size_t lwork = (size_t) m + r + m + nnz0;
	Matrix_Calloc(workA, lwork, int);
	workB = workA + m;
	workC = workB + r;
	pj_   = workC + m;

#define TSORT_LOOP(_CTYPE_, _PTR_, _MASK_, _INCREMENT_) \
	do { \
		_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1, *px_); \
		_MASK_(Matrix_Calloc(px_, nnz0, _CTYPE_)); \
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
			       pj_[workB[pi0[k]]] = pj0[k] ; \
			_MASK_(px_[workB[pi0[k]]] = px0[k]); \
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
					       pj_[kend_] = pj_[k] ; \
					_MASK_(px_[kend_] = px_[k]); \
					++kend_; \
				} else { \
					/* Have already seen this column index */ \
					_MASK_(_INCREMENT_(px_[workB[pj_[k]]], px_[k])); \
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
		Matrix_memset(workB, 0, n, sizeof(int)); \
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
		       PROTECT(*i1 = allocVector(    INTSXP, nnz1)) ; \
		_MASK_(PROTECT(*x1 = allocVector(TYPEOF(x0), nnz1))); \
		       pi1 = INTEGER(*i1) ; \
		_MASK_(px1 =   _PTR_(*x1)); \
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
				_MASK_(px1[workB[pj_[k]]] = px_[k]); \
				++workB[pj_[k]]; \
				++k; \
			} \
			k = workA[i]; \
		} \
		 \
		_MASK_(Matrix_Free(px_, nnz0)); \
		_MASK_(UNPROTECT(1)); /* *px1 */ \
		       UNPROTECT(1) ; /* *pi1 */ \
	} while (0)

	if (!x0)
		TSORT_LOOP(int, LOGICAL, HIDE, INCREMENT_PATTERN);
	else {
		switch (TYPEOF(x0)) {
		case LGLSXP:
			TSORT_LOOP(int, LOGICAL, SHOW, INCREMENT_LOGICAL);
			break;
		case INTSXP:
			TSORT_LOOP(int, INTEGER, SHOW, INCREMENT_INTEGER);
			break;
		case REALSXP:
			TSORT_LOOP(double, REAL, SHOW, INCREMENT_REAL);
			break;
		case CPLXSXP:
			TSORT_LOOP(Rcomplex, COMPLEX, SHOW, INCREMENT_COMPLEX);
			break;
		default:
			break;
		}
	}

#undef TSORT_LOOP

	Matrix_Free(workA, lwork);
	UNPROTECT(1); /* *p1 */
	return;
}

void taggr(SEXP i0, SEXP j0, SEXP x0, SEXP *i1, SEXP *j1, SEXP *x1,
           int m, int n)
{
	R_xlen_t nnz0 = XLENGTH(i0), nnz1 = 0;
	if (nnz0 > INT_MAX)
		error(_("unable to aggregate %s with '%s' and '%s' slots of length exceeding %s"),
		      "TsparseMatrix", "i", "j", "2^31-1");

	/* FIXME: test for overflow and throw error only in that case */

	int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0), *pi1, *pj1,
		i, j, r = (m < n) ? n : m;
	R_xlen_t k, kstart, kend, kend_;

	int *workA, *workB, *workC, *pj_;
	size_t lwork = (size_t) m + r + m + nnz0;
	Matrix_Calloc(workA, lwork, int);
	workB = workA + m;
	workC = workB + r;
	pj_   = workC + m;

#define TAGGR_LOOP(_CTYPE_, _PTR_, _MASK_, _INCREMENT_) \
	do { \
		_MASK_(_CTYPE_ *px0 = _PTR_(x0), *px1, *px_); \
		_MASK_(Matrix_Calloc(px_, nnz0, _CTYPE_)); \
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
			       pj_[workB[pi0[k]]] = pj0[k] ; \
			_MASK_(px_[workB[pi0[k]]] = px0[k]); \
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
					       pj_[kend_] = pj_[k] ; \
					_MASK_(px_[kend_] = px_[k]); \
					++kend_; \
				} else { \
					/* Have already seen this column index */ \
					_MASK_(_INCREMENT_(px_[k], px_[workB[pj_[k]]])); \
				} \
				++k; \
			} \
			workC[i] = kend_; \
			nnz1 += kend_ - kstart; \
		} \
		if (nnz1 != nnz0) { \
			       PROTECT(*i1 = allocVector(    INTSXP, nnz1)) ; \
			       PROTECT(*j1 = allocVector(    INTSXP, nnz1)) ; \
			_MASK_(PROTECT(*x1 = allocVector(TYPEOF(x0), nnz1))); \
			       pi1 = INTEGER(*i1) ; \
			       pj1 = INTEGER(*j1) ; \
			_MASK_(px1 =   _PTR_(*x1)); \
			 \
			k = 0; \
			for (i = 0; i < m; ++i) { \
				kend_ = workC[i]; \
				while (k < kend_) { \
					       *(pi1++) =      i ; \
					       *(pj1++) = pj_[k] ; \
					_MASK_(*(px1++) = px_[k]); \
					++k; \
				} \
				k = workA[i]; \
			} \
			 \
			_MASK_(UNPROTECT(1)); /* *px1 */ \
			       UNPROTECT(2) ; /* *pj1, *px1 */ \
		} \
		_MASK_(Matrix_Free(px_, nnz0)); \
	} while (0)

	if (!x0)
		TAGGR_LOOP(int, LOGICAL, HIDE, );
	else {
		switch (TYPEOF(x0)) {
		case LGLSXP:
			TAGGR_LOOP(int, LOGICAL, SHOW, INCREMENT_LOGICAL);
			break;
		case INTSXP:
			TAGGR_LOOP(int, INTEGER, SHOW, INCREMENT_INTEGER);
			break;
		case REALSXP:
			TAGGR_LOOP(double, REAL, SHOW, INCREMENT_REAL);
			break;
		case CPLXSXP:
			TAGGR_LOOP(Rcomplex, COMPLEX, SHOW, INCREMENT_COMPLEX);
			break;
		default:
			break;
		}
	}

#undef TAGGR_LOOP

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
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
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
SEXP R_sparse_as_Csparse(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_as_Csparse(from, valid[ivalid]);
}

SEXP sparse_as_Rsparse(SEXP from, const char *class)
{
	if (class[2] == 'R')
		return from;

	char cl[] = "..RMatrix";
	cl[0] = class[0];
	cl[1] = class[1];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
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
SEXP R_sparse_as_Rsparse(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_as_Rsparse(from, valid[ivalid]);
}

SEXP sparse_as_Tsparse(SEXP from, const char *class)
{
	if (class[2] == 'T')
		return from;

	char cl[] = "..TMatrix";
	cl[0] = class[0];
	cl[1] = class[1];
	SEXP to = PROTECT(NEW_OBJECT_OF_CLASS(cl));

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
		char ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(to, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (class[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT(from, Matrix_diagSym));
		char di = *CHAR(STRING_ELT(diag, 0));
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
		jSym = (class[2] == 'C') ? Matrix_jSym : Matrix_iSym,
		p = PROTECT(GET_SLOT(from, Matrix_pSym)),
		i = PROTECT(GET_SLOT(from, iSym));
	int *pp = INTEGER(p), *pi, r = (class[2] == 'C') ? n : m, nnz = pp[r];
	if (XLENGTH(i) == nnz) {
		SET_SLOT(to, iSym, i);
		UNPROTECT(1); /* i */
	} else {
		SEXP i_ = PROTECT(allocVector(INTSXP, nnz));
		Matrix_memcpy(INTEGER(i_), INTEGER(i), nnz, sizeof(int));
		SET_SLOT(to, iSym, i_);
		UNPROTECT(2); /* i_, i */
	}
	PROTECT(i = allocVector(INTSXP, nnz));
	SET_SLOT(to, jSym, i);
	pi = INTEGER(i); ++pp;
	int j, k, kend;
	for (j = 0, k = 0; j < r; ++j) {
		kend = pp[j];
		while (k < kend)
			pi[k++] = j;
	}
	UNPROTECT(2); /* i, p */

	if (class[0] != 'n') {
		SEXP x = PROTECT(GET_SLOT(from, Matrix_xSym));
		if (XLENGTH(x) == nnz)
			SET_SLOT(to, Matrix_xSym, x);
		else {
			SEXP x_ = PROTECT(allocVector(TYPEOF(x), nnz));
			SET_SLOT(to, Matrix_xSym, x_);
			switch (class[0]) {
			case 'l':
				Matrix_memcpy(LOGICAL(x_), LOGICAL(x), nnz, sizeof(int));
				break;
			case 'i':
				Matrix_memcpy(INTEGER(x_), INTEGER(x), nnz, sizeof(int));
				break;
			case 'd':
				Matrix_memcpy(REAL(x_), REAL(x), nnz, sizeof(double));
				break;
			case 'z':
				Matrix_memcpy(COMPLEX(x_), COMPLEX(x), nnz, sizeof(Rcomplex));
				break;
			default:
				break;
			}
			UNPROTECT(1); /* x_ */
		}
		UNPROTECT(1); /* x */
	}

	UNPROTECT(1); /* to */
	return to;
}

/* as(<[CRT]sparseMatrix>, "TsparseMatrix") */
SEXP R_sparse_as_Tsparse(SEXP from)
{
	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);

	return sparse_as_Tsparse(from, valid[ivalid]);
}

/* as(<Matrix>, "vector") */
SEXP R_Matrix_as_vector(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	SEXP to = NULL;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);

	switch (cl[2]) {
	case 'e':
		to = GET_SLOT(from, Matrix_xSym);
		if (cl[0] == 'n') {
			R_xlen_t len = XLENGTH(to);
			int *px = LOGICAL(to);
			while (len--)
				if (*(px++) == NA_LOGICAL) {
					PROTECT(to);
					to = duplicate(to);
					UNPROTECT(1);
					break;
				}
		}
		break;
	case 'y':
	case 'r':
	case 'p':
		REPROTECT(from = dense_as_general(from, cl, 1), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	case 'C':
	case 'R':
	case 'T':
		REPROTECT(from = sparse_as_dense(from, cl, 0), pid);
		REPROTECT(from = dense_as_general(from, cl, 0), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	case 'i':
		REPROTECT(from = diagonal_as_dense(from, cl, 'g', 0, '\0'), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	case 'd':
		REPROTECT(from = index_as_dense(from, cl, 'n'), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	default:
		break;
	}

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
	case 'i':
		if (cl[0] == 'n') {
			PROTECT(to);
			na2one(to);
			UNPROTECT(1);
		}
		break;
	default:
		break;
	}

	UNPROTECT(1); /* from */
	return to;
}

/* as(<Matrix>, "matrix") */
SEXP R_Matrix_as_matrix(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	SEXP to = NULL;
	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(from, &pid);

	switch (cl[2]) {
	case 'e':
		PROTECT(to = GET_SLOT(from, Matrix_xSym));
		to = duplicate(to);
		UNPROTECT(1);
		break;
	case 'y':
	case 'r':
	case 'p':
		REPROTECT(from = dense_as_general(from, cl, 1), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	case 'C':
	case 'R':
	case 'T':
		REPROTECT(from = sparse_as_dense(from, cl, 0), pid);
		REPROTECT(from = dense_as_general(from, cl, 0), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	case 'i':
		REPROTECT(from = diagonal_as_dense(from, cl, 'g', 0, '\0'), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	case 'd':
		REPROTECT(from = index_as_dense(from, cl, 'n'), pid);
		to = GET_SLOT(from, Matrix_xSym);
		break;
	default:
		break;
	}
	PROTECT(to);

	SEXP dim = PROTECT(GET_SLOT(from, Matrix_DimSym));
	setAttrib(to, R_DimSymbol, dim);
	UNPROTECT(1); /* dim */

	SEXP dimnames = PROTECT(GET_SLOT(from, Matrix_DimNamesSym));
	if (!DimNames_is_trivial(dimnames))
		setAttrib(to, R_DimNamesSymbol, dimnames);
	UNPROTECT(1); /* dimnames */

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
	case 'i':
		if (cl[0] == 'n')
			na2one(to);
		break;
	default:
		break;
	}

	UNPROTECT(2); /* to, from */
	return to;
}

/* as(<Matrix>, "unpackedMatrix") */
SEXP R_Matrix_as_unpacked(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
		return from;
	case 'p':
		return dense_as_unpacked(from, valid[ivalid]);
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_dense(from, cl, 0);
	case 'i':
		return diagonal_as_dense(from, cl, 't', 0, 'U');
	case 'd':
		return index_as_dense(from, cl, 'n');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "packedMatrix") */
SEXP R_Matrix_as_packed(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	if (cl[1] == 'g' || cl[2] == 'd')
		error(_("attempt to pack a %s"), "generalMatrix");

	switch (cl[2]) {
	case 'y':
	case 'r':
		return dense_as_packed(from, valid[ivalid], '\0', '\0');
	case 'p':
		return from;
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_dense(from, cl, 1);
	case 'i':
		return diagonal_as_dense(from, cl, 't', 1, 'U');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "CsparseMatrix") */
SEXP R_Matrix_as_Csparse(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		return dense_as_sparse(from, cl, 'C');
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Csparse(from, cl);
	case 'i':
		return diagonal_as_sparse(from, cl, 't', 'C', 'U');
	case 'd':
		return index_as_sparse(from, cl, 'n', 'C');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "RsparseMatrix") */
SEXP R_Matrix_as_Rsparse(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		return dense_as_sparse(from, cl, 'R');
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Rsparse(from, cl);
	case 'i':
		return diagonal_as_sparse(from, cl, 't', 'R', 'U');
	case 'd':
		return index_as_sparse(from, cl, 'n', 'R');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "TsparseMatrix") */
SEXP R_Matrix_as_Tsparse(SEXP from)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		return dense_as_sparse(from, cl, 'T');
	case 'C':
	case 'R':
	case 'T':
		return sparse_as_Tsparse(from, cl);
	case 'i':
		return diagonal_as_sparse(from, cl, 't', 'T', 'U');
	case 'd':
		return index_as_sparse(from, cl, 'n', 'T');
	default:
		return R_NilValue;
	}
}

/* as(<Matrix>, "[nldiz]Matrix") */
SEXP R_Matrix_as_kind(SEXP from, SEXP kind, SEXP sparse)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	if (TYPEOF(sparse) != LGLSXP || LENGTH(sparse) < 1)
		error(_("invalid '%s' to %s()"), "sparse", __func__);
	int sparse_ = LOGICAL(sparse)[0];

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		if (sparse_ == NA_LOGICAL || !sparse_)
			from = dense_as_kind(from, cl, kind_);
		else {
			from = dense_as_sparse(from, cl, 'C');
			PROTECT(from);
			char cl_[] = "..CMatrix";
			cl_[0] = cl[0];
			cl_[1] = cl[1];
			from = sparse_as_kind(from, cl_, kind_);
			UNPROTECT(1);
		}
		return from;
	case 'C':
	case 'R':
	case 'T':
		from = sparse_as_kind(from, cl, kind_);
		if (sparse_ != NA_LOGICAL && !sparse_) {
			PROTECT(from);
			char cl_[] = "...Matrix";
			cl_[0] = (kind_ == '.') ? cl[0] : kind_;
			cl_[1] = cl[1];
			cl_[2] = cl[2];
			from = sparse_as_dense(from, cl_, 0);
			UNPROTECT(1);
		}
		return from;
	case 'i':
		/* Ugly because there is no ndiMatrix,
		   and because [ld]diMatrix does not extend [ld]sparseMatrix
		*/
		if (kind_ != 'n') {
			from = diagonal_as_kind(from, cl, kind_);
			if (sparse_ != NA_LOGICAL) {
				PROTECT(from);
				char cl_[] = ".diMatrix";
				cl_[0] = (kind_ == '.') ? cl[0] : kind_;
				if (sparse_)
					from = diagonal_as_sparse(from, cl_, 't', 'C', 'U');
				else
					from = diagonal_as_dense(from, cl_, 't', 0, 'U');
				UNPROTECT(1);
			}
		} else {
			from = diagonal_as_sparse(from, cl, 't', 'C', 'U');
			PROTECT(from);
			char cl_[] = ".tCMatrix";
			cl_[0] = cl[0];
			from = sparse_as_kind(from, cl_, 'n');
			UNPROTECT(1);
			if (sparse_ != NA_LOGICAL && !sparse_) {
				PROTECT(from);
				cl_[0] = 'n';
				from = sparse_as_dense(from, cl_, 0);
				UNPROTECT(1);
			}
		}
		return from;
	case 'd':
		if (sparse_ == NA_LOGICAL || sparse_)
			from = index_as_sparse(from, cl, kind_, '.');
		else
			from = index_as_dense(from, cl, kind_);
		return from;
	default:
		return R_NilValue;
	}
}

/* as(as(<Matrix>, "[nlidz]Matrix"), "generalMatrix") */
SEXP R_Matrix_as_general(SEXP from, SEXP kind)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	int ivalid = R_check_class_etc(from, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(from, __func__);
	const char *cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

	char kind_;
	if (TYPEOF(kind) != STRSXP || LENGTH(kind) < 1 ||
	    (kind = STRING_ELT(kind, 0)) == NA_STRING ||
	    (kind_ = CHAR(kind)[0]) == '\0')
		error(_("invalid '%s' to %s()"), "kind", __func__);

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
	{
		char cl_[] = "...Matrix";
		cl_[0] = (kind_ == '.') ? cl[0] : kind_;
		cl_[1] = cl[1];
		cl_[2] = cl[2];
		int new = (cl[0] == cl_[0]) || ((cl[0] == 'n' && cl_[0] == 'l') ||
		                                (cl[0] == 'l' && cl_[0] == 'n'));
		from = dense_as_kind(from, cl, kind_);
		PROTECT(from);
		from = dense_as_general(from, cl_, new);
		UNPROTECT(1);
		return from;
	}
	case 'C':
	case 'R':
	case 'T':
	{
		char cl_[] = "...Matrix";
		cl_[0] = (kind_ == '.') ? cl[0] : kind_;
		cl_[1] = cl[1];
		cl_[2] = cl[2];
		from = sparse_as_kind(from, cl, kind_);
		PROTECT(from);
		from = sparse_as_general(from, cl_);
		UNPROTECT(1);
		return from;
	}
	case 'i':
	{
		/* As there is no ndiMatrix, we cannot use diagonal_as_kind here: */
		char cl_[] = ".gCMatrix";
		cl_[0] = cl[0];
		from = diagonal_as_sparse(from, cl, 'g', 'C', '\0');
		PROTECT(from);
		from = sparse_as_kind(from, cl_, kind_);
		UNPROTECT(1);
		return from;
	}
	case 'd':
		/* indMatrix extends generalMatrix, but we typically do want this: */
		return index_as_sparse(from, cl, kind_, '.');
	default:
		return R_NilValue;
	}
}
