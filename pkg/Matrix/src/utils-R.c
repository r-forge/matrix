#include "Mdefines.h"
#include "M5.h"

// ../R/Auxiliaries.R  -- indTri(n, packed = FALSE, upper = TRUE, diag = FALSE)
SEXP R_index_triangle(SEXP s_n, SEXP s_packed, SEXP s_upper, SEXP s_diag)
{
	SEXP r;
	int i, j, n = Rf_asInteger(s_n);
	Rboolean packed = Rf_asRboolean(s_packed), upper = Rf_asRboolean(s_upper),
	    diag = Rf_asRboolean(s_diag);
	int_fast64_t
		nn = (int_fast64_t) n * n,
		nx = (packed) ? n + (nn - n) / 2 : nn,
		nr = (diag) ? n + (nn - n) / 2 : (nn - n) / 2;
	if (nx > (1LL << 53))
		Rf_error(_("maximum index would exceed %s"), "2^53");
	if (nr > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

#define DO_INDEX \
	do { \
		if (packed) { \
			if (diag) { \
				while (k <= nr_) \
					*(pr++) = k++; \
			} else if (upper) { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) \
						*(pr++) = k++; \
					k++; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					k++; \
					for (i = j+1; i < n; ++i) \
						*(pr++) = k++; \
				} \
			} \
		} else if (diag) { \
			if (upper) { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i <= j; ++i) \
						*(pr++) = k++; \
					k += n-j-1; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					k += j; \
					for (i = j; i < n; ++i) \
						*(pr++) = k++; \
				} \
			} \
		} else { \
			if (upper) { \
				for (j = 0; j < n; ++j) { \
					for (i = 0; i < j; ++i) \
						*(pr++) = k++; \
					k += n-j; \
				} \
			} else { \
				for (j = 0; j < n; ++j) { \
					k += j+1; \
					for (i = j+1; i < n; ++i) \
						*(pr++) = k++; \
				} \
			} \
		} \
	} while (0)

	if (nx > INT_MAX) {

		PROTECT(r = Rf_allocVector(REALSXP, (R_xlen_t) nr));
		double k = 1.0, nr_ = (double) nr, *pr = REAL(r);

		DO_INDEX;

	} else {

		PROTECT(r = Rf_allocVector(INTSXP, (R_xlen_t) nr));
		int k = 1, nr_ = (int) nr, *pr = INTEGER(r);

		DO_INDEX;

	}

#undef DO_INDEX

	UNPROTECT(1);
	return r;
}

// ../R/Auxiliaries.R  -- indDiag(n, packed = FALSE, upper = TRUE)
SEXP R_index_diagonal(SEXP s_n, SEXP s_packed, SEXP s_upper)
{
	SEXP r;
	int j, n = Rf_asInteger(s_n);
	Rboolean
	    packed = Rf_asRboolean(s_packed),
	    upper  = Rf_asRboolean(s_upper);
	int_fast64_t
		nn = (int_fast64_t) n * n,
		nx = (packed) ? n + (nn - n) / 2 : nn;
	if (nx > (1LL << 53))
		Rf_error(_("maximum index would exceed %s"), "2^53");

#define DO_INDEX \
	do { \
		if (!packed) { \
			for (j = 0; j < n; ++j) { \
				*(pr++) = k++; \
				k += n; \
			} \
		} else if (upper) { \
			for (j = 0; j < n; ++j) { \
				*(pr++) = k; \
				k += j+2; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				*(pr++) = k; \
				k += n-j; \
			} \
		} \
	} while (0)

	if (nx > INT_MAX) {

		PROTECT(r = Rf_allocVector(REALSXP, n));
		double k = 1.0, *pr = REAL(r);

		DO_INDEX;

	} else {

		PROTECT(r = Rf_allocVector(INTSXP, n));
		int k = 1, *pr = INTEGER(r);
		DO_INDEX;

	}

#undef DO_INDEX

	UNPROTECT(1);
	return r;
}

SEXP R_nnz(SEXP s_x, SEXP s_countNA, SEXP s_nnzmax)
{
	int countNA = Rf_asLogical(s_countNA);
	R_xlen_t n = XLENGTH(s_x), nnz = 0;
	double nnzmax = Rf_asReal(s_nnzmax);
	if (!ISNAN(nnzmax) && nnzmax >= 0.0 && nnzmax < (double) n)
		n = (R_xlen_t) nnzmax;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(s_x); \
		if (countNA == NA_LOGICAL) { \
			while (n-- > 0) { \
				if (!c##NOT_NA(*px)) \
					return Rf_ScalarInteger(NA_INTEGER); \
				if (c##NOT_ZERO(*px)) \
					++nnz; \
				++px; \
			} \
		} else if (countNA != 0) { \
			while (n-- > 0) { \
				if (c##NOT_ZERO(*px)) \
					++nnz; \
				++px; \
			} \
		} else { \
			while (n-- > 0) { \
				if (c##NOT_NA(*px) && c##NOT_ZERO(*px)) \
					++nnz; \
				++px; \
			} \
		} \
	} while (0)

	SWITCH4(typeToKind(TYPEOF(s_x)), TEMPLATE);

#undef TEMPLATE

	return (nnz <= INT_MAX)
		? Rf_ScalarInteger((int) nnz) : Rf_ScalarReal((double) nnz);
}


/* ================================================================== */
/* ================================================================== */


#define TRUE_  Rf_ScalarLogical(1)
#define FALSE_ Rf_ScalarLogical(0)

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// all0	 <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
// allFalse <- function(x) !any(x) && !any(is.na(x)) ## ~= all0
SEXP R_all0(SEXP x) {
	if (!Rf_isVectorAtomic(x)) {
		if (Rf_length(x) == 0) return TRUE_;
		// Typically S4.  TODO: Call the R code above, instead!
		Rf_error(_("Argument must be numeric-like atomic vector"));
	}
	R_xlen_t i, n = XLENGTH(x);
	if (n == 0) return TRUE_;

	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *xx = LOGICAL(x);
		for (i = 0; i < n; i++)
			if (xx[i] == NA_LOGICAL || xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	case INTSXP:
	{
		int *xx = INTEGER(x);
		for (i = 0; i < n; i++)
			if (xx[i] == NA_INTEGER || xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	case REALSXP:
	{
		double *xx = REAL(x);
		for (i = 0; i < n; i++)
			if (ISNAN(xx[i]) || xx[i] != 0.) return FALSE_;
		return TRUE_;
	}
	case CPLXSXP:
	{
		Rcomplex *xx = COMPLEX(x);
		for (i = 0; i < n; i++)
			if (ISNAN(xx[i].r) || xx[i].r != 0. ||
			    ISNAN(xx[i].i) || xx[i].i != 0.) return FALSE_;
		return TRUE_;
	}
	case RAWSXP:
	{
		unsigned char *xx = RAW(x);
		for (i = 0; i < n; i++)
			if (xx[i] != 0) return FALSE_;
		return TRUE_;
	}
	}
	Rf_error(_("Argument must be numeric-like atomic vector"));
	return R_NilValue; // -Wall
}

// Fast implementation of [ originally in  ../R/Auxiliaries.R ]
// any0 <- function(x) isTRUE(any(x == 0)) ## ~= anyFalse
// anyFalse <- function(x) isTRUE(any(!x)) ## ~= any0
SEXP R_any0(SEXP x) {
	if (!Rf_isVectorAtomic(x)) {
		if (Rf_length(x) == 0) return FALSE_;
		// Typically S4.  TODO: Call the R code above, instead!
		Rf_error(_("Argument must be numeric-like atomic vector"));
	}
	R_xlen_t i, n = XLENGTH(x);
	if (n == 0) return FALSE_;

	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *xx = LOGICAL(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	case INTSXP:
	{
		int *xx = INTEGER(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	case REALSXP:
	{
		double *xx = REAL(x);
		for (i = 0; i < n; i++) if (xx[i] == 0.) return TRUE_;
		return FALSE_;
	}
	case CPLXSXP:
	{
		Rcomplex *xx = COMPLEX(x);
		for (i = 0; i < n; i++) if (xx[i].r == 0. && xx[i].i == 0.) return TRUE_;
		return FALSE_;
	}
	case RAWSXP:
	{
		unsigned char *xx = RAW(x);
		for (i = 0; i < n; i++) if (xx[i] == 0) return TRUE_;
		return FALSE_;
	}
	}
	Rf_error(_("Argument must be numeric-like atomic vector"));
	return R_NilValue; // -Wall
}

#undef TRUE_
#undef FALSE_

/* used in ../R/construct.R as

    data <- .External(Mmatrix,
                     data, nrow, ncol, byrow, dimnames, mnrow, mncol)

Almost "Cut n Paste" from ...R../src/main/array.c  do_matrix()
*/
SEXP Mmatrix(SEXP args)
{
	SEXP vals, ans, snr, snc, dimnames;
	int nr = 1, nc = 1;

	args = CDR(args); /* skip 'name' */
	vals = CAR(args); args = CDR(args);
	/* Supposedly as.vector() gave a vector type, but we check */
	switch (TYPEOF(vals)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
	case CPLXSXP:
	case STRSXP:
	case RAWSXP:
	case EXPRSXP:
	case VECSXP:
		break;
	default:
		Rf_error(_("'data' must be of a vector type"));
	}
	R_xlen_t lendat = XLENGTH(vals);
	snr = CAR(args); args = CDR(args);
	snc = CAR(args); args = CDR(args);
	Rboolean byrow = Rf_asRboolean(CAR(args)); args = CDR(args);
	dimnames = CAR(args);
	args = CDR(args);
	Rboolean miss_nr = Rf_asRboolean(CAR(args)); args = CDR(args);
	Rboolean miss_nc = Rf_asRboolean(CAR(args));

	if (!miss_nr) {
		if (!Rf_isNumeric(snr)) Rf_error(_("non-numeric matrix extent"));
		nr = Rf_asInteger(snr);
		if (nr == NA_INTEGER)
			Rf_error(_("invalid 'nrow' value (too large or NA)"));
		if (nr < 0)
			Rf_error(_("invalid 'nrow' value (< 0)"));
	}
	if (!miss_nc) {
		if (!Rf_isNumeric(snc)) Rf_error(_("non-numeric matrix extent"));
		nc = Rf_asInteger(snc);
		if (nc == NA_INTEGER)
			Rf_error(_("invalid 'ncol' value (too large or NA)"));
		if (nc < 0)
			Rf_error(_("invalid 'ncol' value (< 0)"));
	}
	if (miss_nr && miss_nc) {
		if (lendat > INT_MAX) Rf_error("data is too long");
		nr = (int) lendat;
	} else if (miss_nr) {
		if (lendat > (double) nc * INT_MAX) Rf_error("data is too long");
		nr = (int) ceil((double) lendat / (double) nc);
	} else if (miss_nc) {
		if (lendat > (double) nr * INT_MAX) Rf_error("data is too long");
		nc = (int) ceil((double) lendat / (double) nr);
	}

	if (lendat > 0) {
		R_xlen_t nrc = (R_xlen_t) nr * nc;
		if (lendat > 1 && nrc % lendat != 0) {
			if ((lendat > nr && (lendat / nr) * nr != lendat) ||
			    (lendat < nr && (nr / lendat) * lendat != nr))
				Rf_warning(_("data length [%lld] is not a sub-multiple "
				             "or multiple of the number of rows [%d]"),
				           (long long) lendat, nr);
			else if ((lendat > nc && (lendat / nc) * nc != lendat) ||
				 (lendat < nc && (nc / lendat) * lendat != nc))
				Rf_warning(_("data length [%lld] is not a sub-multiple "
				             "or multiple of the number of columns [%d]"),
				           (long long) lendat, nc);
		} else if (lendat > 1 && nrc == 0)
			Rf_warning(_("data length exceeds size of matrix"));
	}

#ifndef LONG_VECTOR_SUPPORT
	if ((double) nr * (double) nc > INT_MAX)
		Rf_error(_("too many elements specified"));
#endif

	PROTECT(ans = Rf_allocMatrix(TYPEOF(vals), nr, nc));
	if (Rf_isVector(vals)) {
	    if(lendat)
		Rf_copyMatrix(ans, vals, (Rboolean) byrow);
	    else { /* fill with NAs */
		R_xlen_t N = (R_xlen_t) nr * nc, i;
		switch (TYPEOF(vals)) {
		case STRSXP:
			for (i = 0; i < N; i++)
				SET_STRING_ELT(ans, i, NA_STRING);
			break;
		case LGLSXP:
			for (i = 0; i < N; i++)
				LOGICAL(ans)[i] = NA_LOGICAL;
			break;
		case INTSXP:
			for (i = 0; i < N; i++)
				INTEGER(ans)[i] = NA_INTEGER;
			break;
		case REALSXP:
			for (i = 0; i < N; i++)
				REAL(ans)[i] = NA_REAL;
			break;
		case CPLXSXP:
		{
			/* Initialization must work whether Rcomplex is typedef-ed
			   to a struct { R < 4.3.0 } or to a union { R >= 4.3.0 }
			*/
			Rcomplex zna = { .r = NA_REAL, .i = 0.0 };
			for (i = 0; i < N; i++)
				COMPLEX(ans)[i] = zna;
			break;
		}
		case RAWSXP:
			// FIXME:  N may overflow size_t !!
			memset(RAW(ans), 0, N);
			break;
		default:
			/* don't fill with anything */
			;
		}
	    }
	}
	if (dimnames != R_NilValue && Rf_length(dimnames) > 0)
		ans = Rf_dimnamesgets(ans, dimnames);
	UNPROTECT(1);
	return ans;
}

/**
 * Expand compressed pointers in the array mp into a full set of indices
 * in the array mj.
 *
 * @param ncol number of columns (or rows)
 * @param mp column pointer vector of length ncol + 1
 * @param mj vector of length mp[ncol] to hold the result
 *
 * @return mj
 */
static
int *expand_cmprPt(int ncol, const int mp[], int mj[])
{
	for (int j = 0; j < ncol; j++) {
		int j2 = mp[j+1], jj;
		for (jj = mp[j]; jj < j2; jj++)
			mj[jj] = j;
	}
	return mj;
}

/** From a [CR]sparseMatrix, return a 2 column matrix  `cbind(i, j)` of 0-origin index vectors
 *  (i,j) which entirely correspond to the (i,j) slots of as(x, "TsparseMatrix") :
 * @param x a Csparse- or RsparseMatrix
 * @param col TRUE     /  FALSE
 */
SEXP compressed_non_0_ij(SEXP x, SEXP colP)
{
    Rboolean col = Rf_asRboolean(colP); /* TRUE iff "C"olumn compressed;  FALSE if "R"ow */
    SEXP ans, indSym = col ? Matrix_iSym : Matrix_jSym;
    SEXP indP = PROTECT(GET_SLOT(x, indSym)),
	 pP   = PROTECT(GET_SLOT(x, Matrix_pSym));
    int i, *ij;
    int nouter = INTEGER(GET_SLOT(x, Matrix_DimSym))[col ? 1 : 0],
	n_el   = INTEGER(pP)[nouter]; /* is only == length(indP), if the
				     inner slot is not over-allocated */

    ij = INTEGER(ans = PROTECT(Rf_allocMatrix(INTSXP, n_el, 2)));
    /* expand the compressed margin to 'i' or 'j' : */
    expand_cmprPt(nouter, INTEGER(pP), &ij[col ? n_el : 0]);
    /* and copy the other one: */
    if (col)
	for(i = 0; i < n_el; i++)
	    ij[i] = INTEGER(indP)[i];
    else /* row compressed */
	for(i = 0; i < n_el; i++)
	    ij[i + n_el] = INTEGER(indP)[i];

    UNPROTECT(3);
    return ans;
}

SEXP Matrix_expand_pointers(SEXP pP)
{
	int n = Rf_length(pP) - 1;
	int *p = INTEGER(pP);
	SEXP ans = PROTECT(Rf_allocVector(INTSXP, p[n]));

	expand_cmprPt(n, p, INTEGER(ans));
	UNPROTECT(1);
	return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param ij: 2-column integer matrix
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd(SEXP ij, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
	SEXP ans;
	int *ij_di = NULL, n, nprot=1;
	int check_bounds = Rf_asLogical(chk_bnds), one_ind = Rf_asLogical(orig_1);

	if (TYPEOF(di) != INTSXP) {
		di = PROTECT(Rf_coerceVector(di, INTSXP));
		nprot++;
	}
	if (TYPEOF(ij) != INTSXP) {
		ij = PROTECT(Rf_coerceVector(ij, INTSXP));
		nprot++;
	}
	if (!Rf_isMatrix(ij) ||
	    (ij_di = INTEGER(Rf_getAttrib(ij, R_DimSymbol)))[1] != 2)
		Rf_error(_("Argument ij must be 2-column integer matrix"));
	n = ij_di[0];
	int *Di = INTEGER(di), *IJ = INTEGER(ij),
		*j_ = IJ+n;/* pointer offset! */

	if ((Di[0] * (double) Di[1]) >= 1 + (double)INT_MAX) { /* need double */
		ans = PROTECT(Rf_allocVector(REALSXP, n));
		double *ii = REAL(ans), nr = (double) Di[0];

#define do_ii_FILL(_i_, _j_) \
		int i; \
		if (check_bounds) { \
			for (i = 0; i < n; i++) { \
				if (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER) \
					ii[i] = NA_INTEGER; \
				else { \
					register int i_i, j_i; \
					if (one_ind) { \
						i_i = _i_[i]-1; \
						j_i = _j_[i]-1; \
					} else { \
						i_i = _i_[i]; \
						j_i = _j_[i]; \
					} \
					if (i_i < 0 || i_i >= Di[0]) \
						Rf_error(_("subscript 'i' out of bounds in M[ij]")); \
					if (j_i < 0 || j_i >= Di[1]) \
						Rf_error(_("subscript 'j' out of bounds in M[ij]")); \
					ii[i] = i_i + j_i * nr; \
				} \
			} \
		} else { \
			for (i = 0; i < n; i++) \
				ii[i] = (_i_[i] == NA_INTEGER || _j_[i] == NA_INTEGER) \
					? NA_INTEGER \
					: ((one_ind) \
					   ? ((_i_[i]-1) + (_j_[i]-1) * nr) \
					   :   _i_[i] + _j_[i] * nr); \
		}

		do_ii_FILL(IJ, j_);
	} else {
	ans = PROTECT(Rf_allocVector(INTSXP, n));
	int *ii = INTEGER(ans), nr = Di[0];

	do_ii_FILL(IJ, j_);
	}
	UNPROTECT(nprot);
	return ans;
}

/**
 * Encode Matrix index (i,j)  |-->  i + j * nrow   {i,j : 0-origin}
 *
 * @param i: integer vector
 * @param j: integer vector of same length as 'i'
 * @param orig_1: logical: if TRUE, "1-origin" otherwise "0-origin"
 * @param di: dim(.), i.e. length 2 integer vector
 * @param chk_bnds: logical indicating  0 <= ij[,k] < di[k]  need to be checked.
 *
 * @return encoded index; integer if prod(dim) is small; double otherwise
 */
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds)
{
	SEXP ans;
	int n = LENGTH(i), nprot = 1;
	int check_bounds = Rf_asLogical(chk_bnds), one_ind = Rf_asLogical(orig_1);

	if (TYPEOF(di)!= INTSXP) {
		di = PROTECT(Rf_coerceVector(di,INTSXP));
		nprot++;
	}
	if (TYPEOF(i) != INTSXP) {
		i = PROTECT(Rf_coerceVector(i, INTSXP));
		nprot++;
	}
	if (TYPEOF(j) != INTSXP) {
		j = PROTECT(Rf_coerceVector(j, INTSXP));
		nprot++;
	}
	if (LENGTH(j) != n)
		Rf_error(_("i and j must be integer vectors of the same length"));

	int *Di = INTEGER(di), *i_ = INTEGER(i), *j_ = INTEGER(j);

	if ((Di[0] * (double) Di[1]) >= 1 + (double) INT_MAX) { /* need double */
		ans = PROTECT(Rf_allocVector(REALSXP, n));
		double *ii = REAL(ans), nr = (double) Di[0];

		do_ii_FILL(i_, j_);
	} else {
		ans = PROTECT(Rf_allocVector(INTSXP, n));
		int *ii = INTEGER(ans), nr = Di[0];

		do_ii_FILL(i_, j_);
	}
	UNPROTECT(nprot);
	return ans;
}
#undef do_ii_FILL

#define _rle_d_
#include "t_rle.c"
#undef _rle_d_

#define _rle_i_
#include "t_rle.c"
#undef _rle_i_
