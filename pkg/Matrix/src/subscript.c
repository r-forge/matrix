/* C implementation of methods for [, [[ */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

static
int anyNA(int *p, int n)
{
	if (p && n > 0)
		while (n--)
			if (*(p++) == NA_INTEGER)
				return 1;
	return 0;
}

#define SUB1(c) \
do { \
	if (TYPEOF(s) == INTSXP) \
		SUB1__(i, c); \
	else \
		SUB1__(d, c); \
} while (0)

#define iOOB(s, t, mn) (s == NA_INTEGER || (t = s - 1) >= mn)
#define dOOB(s, t, mn) (ISNAN(s) || s >= 0x1.0p+63 || (t = (int_fast64_t) s - 1) >= mn)

#define lCAST(x) (x)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP dense_subscript_1ary(SEXP obj, const char *class, SEXP s)
{
	R_xlen_t slen = XLENGTH(s);
	SEXP ans = Rf_allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;

	int packed = class[2] == 'p';
	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	R_xlen_t l;
	int_fast64_t b, i, j, k;

	int ge = class[1] == 'g', sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N', up = ul == 'U';

#define SUB1__(d, c) \
	do { \
		d##TYPE *ps = d##PTR(s); \
		c##TYPE *pa = c##PTR(ans); \
		c##TYPE *px = c##PTR(x); \
		for (l = 0; l < slen; ++l) { \
			if (d##OOB(ps[l], k, mn)) { \
				pa[l] = c##NA; \
				continue; \
			} \
			if (ge) { \
				c##ASSIGN_IDEN(pa[l], c##CAST(px[k])); \
				continue; \
			} \
			i = k % m; \
			j = k / m; \
			b = (up) ? j - i : i - j; \
			if (!sy) { \
			if (b < 0) { \
				pa[l] = c##ZERO; \
				continue; \
			} \
			if (un && b == 0) { \
				pa[l] = c##UNIT; \
				continue; \
			} \
			} \
			if (!packed) { \
				if (b >= 0) \
					/* do nothing */; \
				else \
					k = i * m + j; \
			} else if (up) { \
				if (b >= 0) \
					k = DENSE_INDEX_U(i, j, m); \
				else \
					k = DENSE_INDEX_U(j, i, m); \
			} else { \
				if (b >= 0) \
					k = DENSE_INDEX_L(i, j, m); \
				else \
					k = DENSE_INDEX_L(j, i, m); \
			} \
			if (!he || b > 0) \
				c##ASSIGN_IDEN     (pa[l], c##CAST(px[k])); \
			else if (b == 0) \
				c##ASSIGN_PROJ_REAL(pa[l], c##CAST(px[k])); \
			else \
				c##ASSIGN_CONJ     (pa[l], c##CAST(px[k])); \
		} \
	} while (0)

	SWITCH5(class[0], SUB1);

#undef SUB1__

	UNPROTECT(1); /* ans */
	return ans;
}

#undef nCAST
#define nCAST(x) (1)

static
SEXP sparse_subscript_1ary(SEXP obj, const char *class, SEXP s, SEXP o)
{
	R_xlen_t slen = XLENGTH(s);
	SEXP ans = Rf_allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		obj = sparse_as_Csparse(obj, class);
	}
	PROTECT(obj);

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p_ = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i_ = PROTECT(GET_SLOT(obj,        iSym));
	int *pp = INTEGER(p_), *pi = INTEGER(i_), j_, k_, kend_;
	pp++;

	int *po_i = (TYPEOF(o) == INTSXP) ? INTEGER(o) : NULL;
	double *po_d = (TYPEOF(o) == REALSXP) ? REAL(o) : NULL;

	R_xlen_t l = 0, l_;
	int b = 0, i, j;
	int_fast64_t k;

	int byrow = class[2] == 'R',
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N', up = (!byrow) == (ul == 'U');

#define MAP(l) \
	(po_i) ? (R_xlen_t) po_i[l] - 1 : (po_d) ? (R_xlen_t) po_d[l] - 1 : l

#define ADV(d) \
	do { \
		if (l >= slen) \
			j = -1; \
		else { \
			l_ = MAP(l); \
			if (d##OOB(ps[l_], k, mn)) \
				j = -1; \
			else { \
				if (byrow) { \
					i = (int) (k / m); \
					j = (int) (k % m); \
				} \
				else { \
					i = (int) (k % m); \
					j = (int) (k / m); \
				} \
				if (sy && (b = (up) ? j - i : i - j) < 0) \
					SWAP(i, j, int, ); \
				++l; \
			} \
		} \
	} while (0)

#define SUB1__(d, c) \
	do { \
		d##TYPE *ps = d##PTR(s); \
		c##TYPE *pa = c##PTR(ans); \
		c##IF_NPATTERN( \
		SEXP x = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px = c##PTR(x); \
		); \
		ADV(d); \
		while (j >= 0) { \
			j_ = j; \
			k_ = pp[j_ - 1]; \
			kend_ = pp[j_]; \
			while (k_ < kend_ && j == j_) { \
				if (pi[k_] < i) \
					++k_; \
				else { \
					if (pi[k_] > i) \
						pa[l_] = (un && i == j) ? c##UNIT : c##ZERO; \
					else if (!he || b > 0) \
						c##ASSIGN_IDEN     (pa[l_], c##CAST(px[k_])); \
					else if (b == 0) \
						c##ASSIGN_PROJ_REAL(pa[l_], c##CAST(px[k_])); \
					else \
						c##ASSIGN_CONJ     (pa[l_], c##CAST(px[k_])); \
					ADV(d); \
				} \
			} \
			while (j == j_) { \
				pa[l_] = (un && i == j) ? c##UNIT : c##ZERO; \
				ADV(d); \
			} \
		} \
		while (l < slen) { \
			l_ = MAP(l); \
			pa[l_] = c##NA; \
			++l; \
		} \
	} while (0)

	SWITCH5(class[0], SUB1);

#undef MAP
#undef ADV
#undef SUB1__

	UNPROTECT(4); /* i, p, obj, ans */
	return ans;
}

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP diagonal_subscript_1ary(SEXP obj, const char *class, SEXP s)
{
	R_xlen_t slen = XLENGTH(s);
	SEXP ans = Rf_allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);

	int n = DIM(obj)[1];
	int_fast64_t nn = (int_fast64_t) n * n, n1a = (int_fast64_t) n + 1;

	int un = DIAG(obj) != 'N';

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	R_xlen_t l;
	int_fast64_t k;


#define SUB1__(d, c) \
	do { \
		d##TYPE *ps = d##PTR(s); \
		c##TYPE *pa = c##PTR(ans); \
		c##TYPE *px = c##PTR(x); \
		for (l = 0; l < slen; ++l) { \
			if (d##OOB(ps[l], k, nn)) \
				pa[l] = c##NA; \
			else if (k % n1a) \
				pa[l] = c##ZERO; \
			else if (un) \
				pa[l] = c##UNIT; \
			else \
				pa[l] = c##CAST(px[k / n1a]); \
		} \
	} while (0)

	SWITCH5(class[0], SUB1);

#undef SUB1__

	UNPROTECT(1); /* ans */
	return ans;
}

static
SEXP index_subscript_1ary(SEXP obj, const char *class, SEXP s)
{
	R_xlen_t slen = XLENGTH(s);
	SEXP ans = Rf_allocVector(LGLSXP, slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *pa = LOGICAL(ans);

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;

	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	int *pperm = INTEGER(perm), mg = MARGIN(obj);

	R_xlen_t l;
	int_fast64_t k;

#define SUB1__(d, c) \
	do { \
		d##TYPE *ps = d##PTR(s); \
		for (l = 0; l < slen; ++l) { \
			if (d##OOB(ps[l], k, mn)) \
				pa[l] = NA_LOGICAL; \
			else if (mg == 0) \
				pa[l] = k / m == pperm[k % m] - 1; \
			else \
				pa[l] = k % m == pperm[k / m] - 1; \
		} \
	} while (0)

	SUB1();

#undef SUB1__

	UNPROTECT(2); /* perm, ans */
	return ans;
}

/* x[s] where 's' is a vector of type "integer" or "double"      */
/* with elements greater than or equal to 1 (or NA)              */
/* NB: for [CRT], the caller must provide a permutation 'o' of   */
/*     the processed subscript op(s) sorting the latter by row   */
/*     then by column [R] or by column then by row [CT].         */
/*     For general and triangular matrices, 'op' is an identity. */
/*     For symmetric matrices:                                   */
/*                                                               */
/*                     / i+j*m if (i,j) is " on-side"            */
/*         op(i+j*m) =                                           */
/*                     \ j+i*m if (i,j) is "off-side"            */
/*                                                               */
/*     o=NULL <==> o=seq_along(s)                                */
SEXP R_subscript_1ary(SEXP s_obj, SEXP s_s, SEXP s_o)
{
	const char *class = Matrix_class(s_obj, valid_matrix, 7, __func__);
	validObject(s_obj, class);
	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		return    dense_subscript_1ary(s_obj, class, s_s);
	case 'C':
	case 'R':
	case 'T':
		return   sparse_subscript_1ary(s_obj, class, s_s, s_o);
	case 'i':
		return diagonal_subscript_1ary(s_obj, class, s_s);
	case 'd':
		return    index_subscript_1ary(s_obj, class, s_s);
	default:
		return R_NilValue;
	}
}

#undef SUB1

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP dense_subscript_1ary_2col(SEXP obj, const char *class, SEXP s)
{
	int slen = (int) (XLENGTH(s) / 2);
	SEXP ans = Rf_allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int m = DIM(obj)[0], *ps0 = INTEGER(s), *ps1 = ps0 + slen;

	int packed = class[2] == 'p';
	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int l;
	int_fast64_t b, i, j, k;

	int ge = class[1] == 'g', sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N', up = ul == 'U';

#define SUB1(c) \
	do { \
		c##TYPE *pa = c##PTR(ans); \
		c##TYPE *px = c##PTR(x); \
		for (l = 0; l < slen; ++l) { \
			if (ps0[l] == NA_INTEGER || ps1[l] == NA_INTEGER) { \
				pa[l] = c##NA; \
				continue; \
			} \
			i = (int_fast64_t) ps0[l] - 1; \
			j = (int_fast64_t) ps1[l] - 1; \
			if (ge) { \
				c##ASSIGN_IDEN(pa[l], c##CAST(px[DENSE_INDEX_N(i, j, m)])); \
				continue; \
			} \
			b = (up) ? j - i : i - j; \
			if (!sy) { \
			if (b < 0) { \
				pa[l] = c##ZERO; \
				continue; \
			} \
			if (un && b == 0) { \
				pa[l] = c##UNIT; \
				continue; \
			} \
			} \
			if (!packed) { \
				if (b >= 0) \
					k = j * m + i; \
				else \
					k = i * m + j; \
			} else if (up) { \
				if (b >= 0) \
					k = DENSE_INDEX_U(i, j, m); \
				else \
					k = DENSE_INDEX_U(j, i, m); \
			} else { \
				if (b >= 0) \
					k = DENSE_INDEX_L(i, j, m); \
				else \
					k = DENSE_INDEX_L(j, i, m); \
			} \
			if (!he || b > 0) \
				c##ASSIGN_IDEN     (pa[l], c##CAST(px[k])); \
			else if (b == 0) \
				c##ASSIGN_PROJ_REAL(pa[l], c##CAST(px[k])); \
			else \
				c##ASSIGN_CONJ     (pa[l], c##CAST(px[k])); \
		} \
	} while (0)

	SWITCH5(class[0], SUB1);

#undef SUB1

	UNPROTECT(1); /* ans */
	return ans;
}

#undef nCAST
#define nCAST(x) (1)

static
SEXP sparse_subscript_1ary_2col(SEXP obj, const char *class, SEXP s, SEXP o)
{
	int slen = (int) (XLENGTH(s) / 2);
	SEXP ans = Rf_allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *ps0 = INTEGER(s), *ps1 = ps0 + slen;

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		obj = sparse_as_Csparse(obj, class);
	}
	PROTECT(obj);

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g')
		ul = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct = TRANS(obj);
	if (class[1] == 't')
		nu = DIAG(obj);

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p_ = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i_ = PROTECT(GET_SLOT(obj,        iSym));
	int *pp = INTEGER(p_), *pi = INTEGER(i_), j_, k_, kend_;
	pp++;

	int *po_i = (TYPEOF(o) == INTSXP) ? INTEGER(o) : NULL;
	double *po_d = (TYPEOF(o) == REALSXP) ? REAL(o) : NULL;

	int l = 0, l_;
	int b = 0, i, j;

	int byrow = class[2] == 'R',
		sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N', up = (!byrow) == (ul == 'U');

#define MAP(l) \
	(po_i) ? (int) po_i[l] - 1 : (po_d) ? (int) po_d[l] - 1 : l

#define ADV \
	do { \
		if (l >= slen) \
			j = -1; \
		else { \
			l_ = MAP(l); \
			if (ps0[l_] == NA_INTEGER || ps1[l_] == NA_INTEGER) \
				j = -1; \
			else { \
				if (byrow) { \
					i = ps1[l_] - 1; \
					j = ps0[l_] - 1; \
				} \
				else { \
					i = ps0[l_] - 1; \
					j = ps1[l_] - 1; \
				} \
				if (sy && (b = (up) ? j - i : i - j) < 0) \
					SWAP(i, j, int, ); \
				++l; \
			} \
		} \
	} while (0)

#define SUB1(c) \
	do { \
		c##TYPE *pa = c##PTR(ans); \
		c##IF_NPATTERN( \
		SEXP x = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px = c##PTR(x); \
		); \
		ADV; \
		while (j >= 0) { \
			j_ = j; \
			k_ = pp[j_ - 1]; \
			kend_ = pp[j_]; \
			while (k_ < kend_ && j == j_) { \
				if (pi[k_] < i) \
					++k_; \
				else { \
					if (pi[k_] > i) \
						pa[l_] = (un && i == j) ? c##UNIT : c##ZERO; \
					else if (!he || b > 0) \
						c##ASSIGN_IDEN     (pa[l_], c##CAST(px[k_])); \
					else if (b == 0) \
						c##ASSIGN_PROJ_REAL(pa[l_], c##CAST(px[k_])); \
					else \
						c##ASSIGN_CONJ     (pa[l_], c##CAST(px[k_])); \
					ADV; \
				} \
			} \
			while (j == j_) { \
				pa[l_] = (un && i == j) ? c##UNIT : c##ZERO; \
				ADV; \
			} \
		} \
		while (l < slen) { \
			l_ = MAP(l); \
			pa[l_] = c##NA; \
			++l; \
		} \
	} while (0)

	SWITCH5(class[0], SUB1);

#undef MAP
#undef ADV
#undef SUB1

	UNPROTECT(4); /* i, p, obj, ans */
	return ans;
}

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP diagonal_subscript_1ary_2col(SEXP obj, const char *class, SEXP s)
{
	int slen = (int) (XLENGTH(s) / 2);
	SEXP ans = Rf_allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *ps0 = INTEGER(s), *ps1 = ps0 + slen;

	int un = DIAG(obj) != 'N';

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int l;

#define SUB1(c) \
	do { \
		c##TYPE *pa = c##PTR(ans); \
		c##TYPE *px = c##PTR(x); \
		for (l = 0; l < slen; ++l) { \
			if (ps0[l] == NA_INTEGER || ps1[l] == NA_INTEGER) \
				pa[l] = c##NA; \
			else if (ps0[l] != ps1[l]) \
				pa[l] = c##ZERO; \
			else if (un) \
				pa[l] = c##UNIT; \
			else \
				pa[l] = c##CAST(px[ps0[l] - 1]); \
		} \
	} while (0)

	SWITCH5(class[0], SUB1);

#undef SUB1

	UNPROTECT(1); /* ans */
	return ans;
}

static
SEXP index_subscript_1ary_2col(SEXP obj, const char *class, SEXP s)
{
	int slen = (int) (XLENGTH(s) / 2);
	SEXP ans = Rf_allocVector(LGLSXP, slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *ps0 = INTEGER(s), *ps1 = ps0 + slen, *pa = LOGICAL(ans);

	SEXP perm = PROTECT(GET_SLOT(obj, Matrix_permSym));
	int *pperm = INTEGER(perm), mg = MARGIN(obj);

	int l;

	for (l = 0; l < slen; ++l) {
		if (ps0[l] == NA_INTEGER || ps1[l] == NA_INTEGER)
			pa[l] = NA_LOGICAL;
		else if (mg == 0)
			pa[l] = ps1[l] == pperm[ps0[l]] - 1;
		else
			pa[l] = ps0[l] == pperm[ps1[l]] - 1;
	}

	UNPROTECT(2); /* perm, ans */
	return ans;
}

/* x[s] where 's' is a 2-column matrix of type "integer"         */
/* with s[, 1L] in 1:m (or NA) and s[, 2L] in 1:n (or NA)        */
/* NB: for [CRT], the caller must provide a permutation 'o' of   */
/*     the processed subscript op(s) sorting the latter by row   */
/*     then by column [R] or by column then by row [CT].         */
/*     For general and triangular matrices, 'op' is an identity. */
/*     For symmetric matrices:                                   */
/*                                                               */
/*                     / (i,j) if (i,j) is " on-side"            */
/*         op((i,j)) =                                           */
/*                     \ (j,i) if (i,j) is "off-side"            */
/*                                                               */
/*     o=NULL <==> o=seq_len(nrow(s))                            */
SEXP R_subscript_1ary_2col(SEXP s_obj, SEXP s_s, SEXP s_o)
{
	const char *class = Matrix_class(s_obj, valid_matrix, 7, __func__);
	validObject(s_obj, class);
	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		return    dense_subscript_1ary_2col(s_obj, class, s_s);
	case 'C':
	case 'R':
	case 'T':
		return   sparse_subscript_1ary_2col(s_obj, class, s_s, s_o);
	case 'i':
		return diagonal_subscript_1ary_2col(s_obj, class, s_s);
	case 'd':
		return    index_subscript_1ary_2col(s_obj, class, s_s);
	default:
		return R_NilValue;
	}
}

static
int stay_sy(int *pi, int ni, int *pj, int nj, int n, char uplo, int checkNA)
{
	if (!pi || !pj || ni != nj ||
	    memcmp(pi, pj, sizeof(int) * (size_t) ni) != 0)
		return 0;
	int ki, r = (uplo == 'U') ? 1 : -1;
	if (checkNA && anyNA(pi, ni))
		return r;
	if (ni >= 2) {
		/* triangular iff monotone */
		if (pi[0] == pi[1])
			return r;
		else if (pi[0] < pi[1]) {
			for (ki = 2; ki < ni; ++ki)
				if (pi[ki - 1] >= pi[ki])
					return r;
			/* up->up, lo->lo */
		} else {
			for (ki = 2; ki < ni; ++ki)
				if (pi[ki - 1] <= pi[ki])
					return r;
			/* up->lo, lo->up */
			r = -r;
		}
	}
	return 2 * r;
}

static
int stay_tr(int *pi, int ni, int *pj, int nj, int n, char uplo, int checkNA)
{
	if (!pi || !pj || ni != nj)
		return 0;
	int ki, kj, iden = memcmp(pi, pj, sizeof(int) * (size_t) ni) == 0,
		r = (uplo == 'U') ? 1 : -1;
	if (checkNA && (anyNA(pi, ni) || (!iden && anyNA(pj, ni))))
		return 0;
	if (iden) {
		/* triangular iff monotone; unit diagonal is preserved */
		if (ni >= 2) {
			if (pi[0] == pi[1])
				return 0;
			else if (pi[0] < pi[1]) {
				for (ki = 2; ki < ni; ++ki)
					if (pi[ki - 1] >= pi[ki])
						return 0;
				/* up->up, lo->lo */
			} else {
				for (ki = 2; ki < ni; ++ki)
					if (pi[ki - 1] <= pi[ki])
						return 0;
				/* up->lo, lo->up */
				r = -r;
			}
		}
		return r * 2;
	} else {
		/* brute force ... */
		int j;
		if (uplo == 'U') {
			for (kj = 0; kj < nj; ++kj)
				for (ki = kj + 1, j = pj[kj]; ki < ni; ++ki)
					if (pi[ki] <= j)
						goto LO;
			/* up->up */
			return  r;
		LO:
			for (kj = 0; kj < nj; ++kj)
				for (ki = 0, j = pj[kj]; ki < kj; ++ki)
					if (pi[ki] <= j)
						return 0;
			/* up->lo */
			return -r;
		} else {
			for (kj = 0; kj < nj; ++kj)
				for (ki = 0, j = pj[kj]; ki < kj; ++ki)
					if (pi[ki] >= j)
						goto UP;
			/* lo->lo */
			return  r;
		UP:
			for (kj = 0; kj < nj; ++kj)
				for (ki = kj + 1, j = pj[kj]; ki < ni; ++ki)
					if (pi[ki] >= j)
						return 0;
			/* lo->up */
			return -r;
		}
	}
}

static
int stay_di(int *pi, int ni, int *pj, int nj, int n, int checkNA)
{
	if (!pi || !pj || ni != nj)
		return 0;
	int ki, kj, iden = memcmp(pi, pj, sizeof(int) * (size_t) ni) == 0;
	if (checkNA && (anyNA(pi, ni) || (!iden && anyNA(pj, ni))))
		return 0;
	if (iden) {
		/* diagonal iff no duplicates; unit diagonal is preserved */
		char *iwork;
		Matrix_Calloc(iwork, n, char);
		--iwork;
		for (ki = 0; ki < ni; ++ki) {
			if (iwork[pi[ki]])
				return 0;
			iwork[pi[ki]] = 1;
		}
		++iwork;
		Matrix_Free(iwork, n);
		return 2;
	} else {
		/* brute force ... */
		int j;
		for (kj = 0; kj < nj; ++kj) {
			j = pj[kj];
			for (ki = 0; ki < kj; ++ki)
				if (pi[ki] == j)
					return 0;
			for (ki = kj + 1; ki < ni; ++ki)
				if (pi[ki] == j)
					return 0;
		}
		return 1;
	}
}

static
SEXP dense_subscript_2ary(SEXP obj, const char *class, SEXP si, SEXP sj)
{
	if (si == R_NilValue && sj == R_NilValue)
		return obj;

	int *pdim = DIM(obj), m = pdim[0], n = pdim[1],
		ki, mi = si == R_NilValue, ni = (mi) ? m : LENGTH(si),
		kj, mj = sj == R_NilValue, nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);

	int packed = class[2] == 'p';
	char ul0 = '\0', ct0 = '\0', nu0 = '\0';
	if (class[1] != 'g')
		ul0 = UPLO(obj);
	if (class[1] == 's' && class[0] == 'z')
		ct0 = TRANS(obj);
	if (class[1] == 't')
		nu0 = DIAG(obj);

	int stay = (class[1] == 'g') ? 0 : (class[1] == 's')
		? stay_sy(pi, ni, pj, nj, n, ul0, 1)
		: stay_tr(pi, ni, pj, nj, n, ul0, 1);
	int_fast64_t ninj = (int_fast64_t) ni * nj,
		xlen = (!packed || stay == 0) ? ninj : ni + (ninj - ni) / 2;
	if (xlen > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	char ul1 = '\0', ct1 = '\0', nu1 = '\0';
	if (stay != 0) {
	if (class[1] != 'g')
		ul1 = (stay > 0) ? 'U' : 'L';
	if (class[1] == 's' && class[0] == 'z')
		ct1 = ct0;
	if (class[1] == 't')
		nu1 = (ABS(stay) == 1) ? 'N' : nu0;
	}

	int he = class[1] == 's' && ct0 == 'C' && ct1 != 'C',
		un = class[1] == 't' && nu0 != 'N';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (stay == 0) ? 'g' : class[1];
	cl[2] = (stay == 0) ? 'e' : class[2];
	SEXP ans = PROTECT(newObject(cl));

	SET_DIM(ans, ni, nj);
	if (cl[1] != 'g' && ul1 != 'U')
		SET_UPLO(ans);
	if (cl[1] == 's' && ct1 != 'C' && cl[0] == 'z')
		SET_TRANS(ans);
	if (cl[1] == 't' && nu1 != 'N')
		SET_DIAG(ans);

	SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), (R_xlen_t) xlen));
	SET_SLOT(ans, Matrix_xSym, x1);

	int_fast64_t i, j;

#define SUB2(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed && cl[1] != 'g') \
			memset(px1, 0, sizeof(c##TYPE) * (size_t) XLENGTH(x1)); \
		if (class[1] == 'g') { \
			SUB2__(c, ASSIGN_GE, N, N, for (ki =  0; ki <  ni; ++ki), , ); \
		} else if (class[1] == 's' && !packed) { \
			if (ul1 == '\0') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_SY, U, N, for (ki =  0; ki <  ni; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_SY, L, N, for (ki =  0; ki <  ni; ++ki), , ); \
			} else if (ul1 == 'U') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_SY, U, N, for (ki =  0; ki <= kj; ++ki), \
			       , px1 += ni - kj - 1); \
			else \
			SUB2__(c, ASSIGN_SY, L, N, for (ki =  0; ki <= kj; ++ki), \
			       , px1 += ni - kj - 1); \
			} else { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_SY, U, N, for (ki = kj; ki <  ni; ++ki), \
			       px1 += kj, ); \
			else \
			SUB2__(c, ASSIGN_SY, L, N, for (ki = kj; ki <  ni; ++ki), \
			       px1 += kj, ); \
			} \
		} else if (class[1] == 's') { \
			if (ul1 == '\0') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_SY, U, U, for (ki =  0; ki <  ni; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_SY, L, L, for (ki =  0; ki <  ni; ++ki), , ); \
			} else if (ul1 == 'U') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_SY, U, U, for (ki =  0; ki <= kj; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_SY, L, L, for (ki =  0; ki <= kj; ++ki), , ); \
			} else { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_SY, U, U, for (ki = kj; ki <  ni; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_SY, L, L, for (ki = kj; ki <  ni; ++ki), , ); \
			} \
		} else if (!packed) { \
			if (ul1 == '\0') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_TR, U, N, for (ki =  0; ki <  ni; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_TR, L, N, for (ki =  0; ki <  ni; ++ki), , ); \
			} else if (ul1 == 'U') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_TR, U, N, for (ki =  0; ki <= kj; ++ki), \
			       , px1 += ni - kj - 1); \
			else \
			SUB2__(c, ASSIGN_TR, L, N, for (ki =  0; ki <= kj; ++ki), \
			       , px1 += ni - kj - 1); \
			} else { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_TR, U, N, for (ki = kj; ki <  ni; ++ki), \
			       px1 += kj, ); \
			else \
			SUB2__(c, ASSIGN_TR, L, N, for (ki = kj; ki <  ni; ++ki), \
			       px1 += kj, ); \
			} \
		} else { \
			if (ul1 == '\0') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_TR, U, U, for (ki =  0; ki <  ni; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_TR, L, L, for (ki =  0; ki <  ni; ++ki), , ); \
			} else if (ul1 == 'U') { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_TR, U, U, for (ki =  0; ki <= kj; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_TR, L, L, for (ki =  0; ki <= kj; ++ki), , ); \
			} else { \
			if (ul0 == 'U') \
			SUB2__(c, ASSIGN_TR, U, U, for (ki = kj; ki <  ni; ++ki), , ); \
			else \
			SUB2__(c, ASSIGN_TR, L, L, for (ki = kj; ki <  ni; ++ki), , ); \
			} \
		} \
	} while (0)

#define SUB2__(c, assign, s, t, __for__, jump0, jump1) \
	do { \
		for (kj = 0; kj < nj; ++kj) { \
			if (mj) \
				j = kj; \
			else if (pj[kj] != NA_INTEGER) \
				j = pj[kj] - 1; \
			else { \
				jump0; \
				__for__ { \
					*(px1++) = c##NA; \
				} \
				jump1; \
				continue; \
			} \
			jump0; \
			__for__ { \
				if (mi) \
					i = ki; \
				else if (pi[ki] != NA_INTEGER) \
					i = pi[ki] - 1; \
				else { \
					*(px1++) = c##NA; \
					continue; \
				} \
				assign(c, px0, px1, i, j, m, s, t); \
				px1++; \
			} \
			jump1; \
		} \
	} while (0)

#define ASSIGN_GE(c, x, y, i, j, m, s, t) \
	do { \
		c##ASSIGN_IDEN(*y, x[DENSE_INDEX_N(i, j, m)]); \
	} while (0)

#define ASSIGN_SY(c, x, y, i, j, m, s, t) \
	do { \
		if (he && i == j) \
			c##ASSIGN_PROJ_REAL(*y, x[DENSE_INDEX_##t(i, j, m)]); \
		else if (CMP_##s(i, j)) \
			c##ASSIGN_IDEN     (*y, x[DENSE_INDEX_##t(i, j, m)]); \
		else if (he) \
			c##ASSIGN_CONJ     (*y, x[DENSE_INDEX_##t(j, i, m)]); \
		else \
			c##ASSIGN_IDEN     (*y, x[DENSE_INDEX_##t(j, i, m)]); \
	} while (0)

#define ASSIGN_TR(c, x, y, i, j, m, s, t) \
	do { \
		if (un && i == j) \
			c##ASSIGN_IDEN(*y, c##UNIT); \
		else if (CMP_##s(i, j)) \
			c##ASSIGN_IDEN(*y, x[DENSE_INDEX_##t(i, j, m)]); \
		else \
			c##ASSIGN_IDEN(*y, c##ZERO); \
	} while (0)

#define CMP_U(i, j) i <= j
#define CMP_L(i, j) i >= j

	SWITCH4(class[0], SUB2);

#undef SUB2__
#undef SUB2

	UNPROTECT(3); /* x1, x0, ans */
	return ans;
}

static
SEXP sparse_subscript_2ary(SEXP obj, const char *class, SEXP si, SEXP sj)
{
	if (si == R_NilValue && sj == R_NilValue)
		return obj;

	int mg = (class[2] == 'R') ? 0 : 1;
	if (mg == 0)
		SWAP(si, sj, SEXP, );

	int *pdim = DIM(obj), m = pdim[mg == 0], n = pdim[mg != 0],
		ki, mi = si == R_NilValue, ni = (mi) ? m : LENGTH(si),
		kj, mj = sj == R_NilValue, nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);
	if (anyNA(pi, ni) || anyNA(pj, nj))
		Rf_error(_("NA subscripts in %s not supported for '%s' inheriting from %s"),
		         "x[i, j]", "x", "sparseMatrix");

	char ul0 = '\0', ct0 = '\0', nu0 = '\0';
	if (class[1] != 'g') {
		ul0 = UPLO(obj);
		if (mg == 0)
			ul0 = (ul0 == 'U') ? 'L' : 'U';
	}
	if (class[1] == 's' && class[0] == 'z')
		ct0 = TRANS(obj);
	if (class[1] == 't')
		nu0 = DIAG(obj);

	int stay = (class[1] == 'g') ? 0 : (class[1] == 's')
		? stay_sy(pi, ni, pj, nj, n, ul0, 0)
		: stay_tr(pi, ni, pj, nj, n, ul0, 0);

	char ul1 = '\0', ct1 = '\0', nu1 = '\0';
	if (stay != 0) {
	if (class[1] != 'g')
		ul1 = (stay > 0) ? 'U' : 'L';
	if (class[1] == 's' && class[0] == 'z')
		ct1 = ct0;
	if (class[1] == 't')
		nu1 = (ABS(stay) == 1) ? 'N' : nu0;
	}

	char class__[] = "...Matrix";
	class__[0] = class[0];
	class__[1] = class[1];
	class__[2] = class[2];
	PROTECT_INDEX pid_obj;
	PROTECT_WITH_INDEX(obj, &pid_obj);
	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		REPROTECT(obj = sparse_as_Csparse(obj, class__), pid_obj);
		class__[2] = 'C';
	}
	if (class[1] != 'g' && ABS(stay) <= 1) {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		REPROTECT(obj = sparse_as_general(obj, class__), pid_obj);
		class__[1] = 'g';
	}

	char cl[] = "...Matrix";
	cl[0] = class  [0];
	cl[1] = class  [1];
	cl[2] = class__[2];
	if (ABS(stay) <= ((class[1] != 's') ? 0 : 1))
		cl[1] = 'g';
	PROTECT_INDEX pid_ans;
	SEXP ans;
	PROTECT_WITH_INDEX(ans = newObject(cl), &pid_ans);

	SET_DIM(ans, (mg == 0) ? nj : ni, (mg == 0) ? ni : nj);
	if (cl[1] != 'g' && (mg == (ul1 != 'U')))
		SET_UPLO(ans);
	if (cl[1] == 's' && ct1 != 'C' && cl[0] == 'z')
		SET_TRANS(ans);
	if (cl[1] == 't' && nu1 != 'N')
		SET_DIAG(ans);

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(obj,        iSym)),
		p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nj + 1));
	int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), *pp1 = INTEGER(p1),
		i, j, k, kend;
	pp0++; *(pp1++) = 0;
	SET_SLOT(ans, Matrix_pSym, p1);

	if (mi) {

		int_fast64_t nnz = 0;
		for (kj = 0; kj < nj; ++kj) {
			j = pj[kj] - 1;
			nnz += pp0[j] - pp0[j - 1];
			pp1[kj] = (int) nnz;
		}
		if (nnz > INT_MAX)
			Rf_error(_("%s too dense for %s; would have more than %s nonzero entries"),
			         "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, (int) nnz));
		int *pi1 = INTEGER(i1);
		SET_SLOT(ans, iSym, i1);

#define SUB2(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, (int) nnz)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(ans, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			for (kj = 0; kj < nj; ++kj) { \
				j = pj[kj] - 1; \
				k = pp0[j - 1]; \
				kend = pp0[j]; \
				while (k < kend) { \
					*(pi1++) = pi0[k]; \
					c##IF_NPATTERN( \
					*(px1++) = px0[k]; \
					); \
					++k; \
				} \
			} \
		} while (0)

		SWITCH5(class[0], SUB2);

#undef SUB2

		UNPROTECT(1); /* i1 */

	} else {

		int *iwork;
		size_t liwork = (size_t) ((int_fast64_t) m + m + ni);
		Matrix_Calloc(iwork, liwork, int);

		int *iworkA = iwork, *iworkB = iworkA + m, *iworkC = iworkB + m,
			prev = m, sort = 0;

		/* workA[ i]:    size of set { ki  : pi[ki] - 1 == i }             */
		/* workB[ i]: minimum of set { ki  : pi[ki] - 1 == i }             */
		/* workC[ki]: minimum of set { ki' : ki' > ki, pi[ki'] == pi[ki] } */

		for (ki = ni - 1; ki >= 0; --ki) {
			i = pi[ki] - 1;
			++iworkA[i];
			iworkC[ki] = iworkB[i];
			iworkB[i] = ki;
			if (i > prev)
				sort = 1;
			prev = i;
		}

		int_fast64_t nnz = 0;
		for (kj = 0; kj < nj; ++kj) {
			j = (mj) ? kj : pj[kj] - 1;
			k = pp0[j - 1];
			kend = pp0[j];
			while (k < kend) {
				nnz += iworkA[pi0[k]];
				++k;
			}
			pp1[kj] = (int) nnz;
		}
		if (nnz > INT_MAX)
			Rf_error(_("%s too dense for %s; would have more than %s nonzero entries"),
			         "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, (int) nnz));
		int *pi1 = INTEGER(i1), d;
		SET_SLOT(ans, iSym, i1);

#define SUB2(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, (int) nnz)); \
			c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
			SET_SLOT(ans, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			for (kj = 0; kj < nj; ++kj) { \
				j = (mj) ? kj : pj[kj] - 1; \
				k = pp0[j - 1]; \
				kend = pp0[j]; \
				while (k < kend) { \
					i = pi0[k]; \
					d = iworkA[i]; \
					ki = iworkB[i]; \
					while (d--) { \
						*(pi1++) = ki; \
						c##IF_NPATTERN( \
						*(px1++) = px0[k]; \
						); \
						ki = iworkC[ki]; \
					} \
					++k; \
				} \
			} \
		} while (0)

		SWITCH5(class[0], SUB2);

#undef SUB2

		Matrix_Free(iwork, liwork);

		if (sort) {
			liwork = (size_t) ((int_fast64_t) ni + 1 + ((ni < nj) ? nj : ni) + nnz);
			Matrix_Calloc(iwork, liwork, int);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px1 = NULL, *work = NULL; \
				c##IF_NPATTERN( \
				SEXP x1 = PROTECT(GET_SLOT(ans, Matrix_xSym)); \
				px1 = c##PTR(x1); \
				Matrix_Calloc(work, nnz, c##TYPE); \
				); \
				c##cspsort(INTEGER(p1), INTEGER(i1), px1, ni, nj, iwork, work); \
				c##IF_NPATTERN( \
				Matrix_Free(work, nnz); \
				UNPROTECT(1); /* x1 */ \
				); \
			} while (0)

			SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

			Matrix_Free(iwork, liwork);
		}

		UNPROTECT(1); /* i1 */

	}

	if (class[1] == 's' && ABS(stay) == 1) {
		/* defined in ./forceSymmetric.c : */
		SEXP sparse_force_symmetric(SEXP, const char *, char, char);
		REPROTECT(ans = sparse_force_symmetric(ans, cl, (mg == (ul1 == 'U')) ? 'U' : 'L', ct1), pid_ans);
		cl[1] = 's';
	}
	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Tsparse(SEXP, const char *);
		REPROTECT(ans = sparse_as_Tsparse(ans, cl), pid_ans);
		cl[2] = 'T';
	}

	UNPROTECT(5); /* p1, i0, p0, ans, obj */
	return ans;
}

static
SEXP diagonal_subscript_2ary(SEXP obj, const char *class, SEXP si, SEXP sj)
{
	if (si == R_NilValue && sj == R_NilValue)
		return obj;

	int n = DIM(obj)[1],
		ki, mi = si == R_NilValue, ni = (mi) ? n : LENGTH(si),
		kj, mj = sj == R_NilValue, nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);
	if (anyNA(pi, ni) || anyNA(pj, nj))
		Rf_error(_("NA subscripts in %s not supported for '%s' inheriting from %s"),
		         "x[i, j]", "x", "sparseMatrix");
	int stay = stay_di(pi, ni, pj, nj, n, 0);

	char nu0 = DIAG(obj),
		nu1 = (stay <= 0) ? '\0' : (stay <= 1) ? 'N' : nu0;

	int un = nu0 != 'N';

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (stay <= 0) ? 'g' : 'd';
	cl[2] = (stay <= 0) ? 'C' : 'i';
	SEXP ans = PROTECT(newObject(cl));

	SET_DIM(ans, ni, nj);

	if (nu1 != '\0' && nu1 != 'N') {

		SET_DIAG(ans);

	} else if (nu1 != '\0') {

		SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			x1 = PROTECT(Rf_allocVector(TYPEOF(x0), ni));
		SET_SLOT(ans, Matrix_xSym, x1);

		int j;

#define SUB2(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0) - 1, *px1 = c##PTR(x1); \
			for (ki = 0; ki < ni; ++ki) \
				if (*(pi++) != (j = *(pj++))) \
					*(px1++) = c##ZERO; \
				else \
					*(px1++) = (un) ? c##UNIT : px0[j]; \
		} while (0)

		SWITCH4(class[0], SUB2);

#undef SUB2

		UNPROTECT(2); /* x1, x0 */

	} else {

		SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym));

		SEXP p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nj + 1));
		int *pp1 = INTEGER(p1), j;
		*(pp1++) = 0;
		SET_SLOT(ans, Matrix_pSym, p1);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			for (kj = 0; kj < nj; ++kj) { \
				j = (mj) ? kj : pj[kj] - 1; \
				pp1[kj] = 0; \
				if (un || c##NOT_ZERO(px0[j])) { \
					if (mi) { \
						for (ki = 0; ki < ni; ++ki) \
							if (ki == j) \
								++pp1[kj]; \
					} else { \
						for (ki = 0; ki < ni; ++ki) \
							if (pi[ki] - 1 == j) \
								++pp1[kj]; \
					} \
					if (pp1[kj] > INT_MAX - pp1[kj - 1]) \
						break; \
				} \
				pp1[kj] += pp1[kj - 1]; \
			} \
		} while (0)

		kj = -1;
		SWITCH4(class[0], TEMPLATE);

#undef TEMPLATE

		if (kj < nj)
			Rf_error(_("%s too dense for %s; would have more than %s nonzero entries"),
			         "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, pp1[nj - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(ans, Matrix_iSym, i1);

#define SUB2(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			c##IF_NPATTERN( \
			SEXP x1 = PROTECT(Rf_allocVector(c##TYPESXP, pp1[nj - 1])); \
			c##TYPE *px1 = c##PTR(x1); \
			SET_SLOT(ans, Matrix_xSym, x1); \
			UNPROTECT(1); /* x1 */ \
			); \
			for (kj = 0; kj < nj; ++kj) { \
				j = (mj) ? kj : pj[kj] - 1; \
				if (un || c##NOT_ZERO(px0[j])) { \
					if (mi) { \
						for (ki = 0; ki < ni; ++ki) \
							if (ki == j) { \
								*(pi1++) = ki; \
								c##IF_NPATTERN( \
								*(px1++) = (un) ? c##UNIT : px0[j]; \
								); \
							} \
					} else { \
						for (ki = 0; ki < ni; ++ki) \
							if (pi[ki] - 1 == j) { \
								*(pi1++) = ki; \
								c##IF_NPATTERN( \
								*(px1++) = (un) ? c##UNIT : px0[j]; \
								); \
							} \
					} \
				} \
			} \
		} while (0)

		SWITCH4(class[0], SUB2);

#undef SUB2

		UNPROTECT(3); /* i1, p1, x0 */

	}

	UNPROTECT(1); /* ans */
	return ans;
}

static
SEXP index_subscript_2ary(SEXP obj, const char *class, SEXP si, SEXP sj)
{
	if (si == R_NilValue && sj == R_NilValue)
		return obj;

	/* defined in ./perm.c : */
	int isPerm(const int *, int, int);
	void invertPerm(const int *, int *, int, int, int);

	int mg = MARGIN(obj);
	if (mg != 0)
		SWAP(si, sj, SEXP, );

	int *pdim = DIM(obj), m = pdim[mg != 0], n = pdim[mg == 0],
		ki, mi = si == R_NilValue, ni = (mi) ? m : LENGTH(si),
		kj, mj = sj == R_NilValue, nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);
	if (anyNA(pi, ni) || anyNA(pj, nj))
		Rf_error(_("NA subscripts in %s not supported for '%s' inheriting from %s"),
		         "x[i, j]", "x", "sparseMatrix");
	int stay = class[0] == 'p';

	PROTECT_INDEX pid_perm;
	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	int *pperm = INTEGER(perm);
	PROTECT_WITH_INDEX(perm, &pid_perm);

	PROTECT_INDEX pid_ans;
	SEXP ans = obj;
	PROTECT_WITH_INDEX(ans, &pid_ans);

	if (!mi) {
		stay = stay && ni == m && isPerm(pi, m, 1);
		REPROTECT(ans = newObject((stay) ? "pMatrix" : "indMatrix"), pid_ans);

		SET_DIM(ans, (mg == 0) ? ni : n, (mg == 0) ? n : ni);
		SET_MARGIN(ans, mg);

		int *tmp = pperm;
		REPROTECT(perm = Rf_allocVector(INTSXP, ni), pid_perm);
		pperm = INTEGER(perm);

		--tmp; /* 1-indexed */
		for (ki = 0; ki < ni; ++ki)
			pperm[ki] = tmp[pi[ki]];
		++tmp; /* 0-indexed */

		SET_SLOT(ans, Matrix_permSym, perm);
		m = ni;
	}

	if (!mj) {
		stay = stay && nj == n && isPerm(pj, n, 1);
		REPROTECT(ans = newObject((stay) ? "pMatrix" : (mg == 0) ? "ngCMatrix" : "ngRMatrix"), pid_ans);

		SET_DIM(ans, (mg == 0) ? m : nj, (mg == 0) ? nj : m);

		int *iwork;
		size_t liwork = (size_t) ((stay) ? nj : (int_fast64_t) n + 1 + n + m);
		Matrix_Calloc(iwork, liwork, int);

		if (stay) {

		invertPerm(pj, iwork, nj, 1, 1);

		int *tmp = pperm;
		REPROTECT(perm = Rf_allocVector(INTSXP, nj), pid_perm);
		pperm = INTEGER(perm);

		--iwork; /* 1-indexed */
		for (kj = 0; kj < nj; ++kj)
			pperm[kj] = iwork[tmp[kj]];
		++iwork; /* 0-indexed */

		SET_SLOT(ans, Matrix_permSym, perm);

		} else {

		int *pp0 = iwork, *pp_ = pp0 + n + 1, *pi0 = pp_ + n;

		SEXP p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) nj + 1));
		int *pp1 = INTEGER(p1), i, j, k, kend, k_;
		*(pp1++) = 0;
		SET_SLOT(ans, Matrix_pSym, p1);

		/* 1. Compute old column counts */
		for (i = 0; i < m; ++i)
			++pp0[pperm[i]];

		/* 2. Compute new column pointers */
		int_fast64_t nnz = 0;
		for (kj = 0; kj < nj; ++kj) {
			nnz += pp0[pj[kj]];
			pp1[kj] = (int) nnz;
		}
		if (nnz > INT_MAX)
			Rf_error(_("%s too dense for %s; would have more than %s nonzero entries"),
			         "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(Rf_allocVector(INTSXP, pp1[nj - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(ans, (mg == 0) ? Matrix_iSym : Matrix_jSym, i1);

		/* 3. Compute old column pointers */
		for (j = 0; j < n; ++j)
			pp0[j + 1] += (pp_[j] = pp0[j]);

		/* 4. Sort old row indices by column */
		--pp0;
		for (i = 0; i < m; ++i)
			pi0[pp0[pperm[i]]++] = i;

		/* 5. Copy old row indices */
		--pp_;
		for (kj = 0, k = 0; kj < nj; ++kj) {
			k_ = pp_[pj[kj]];
			kend = pp1[kj];
			while (k < kend)
				pi1[k++] = pi0[k_++];
		}

		UNPROTECT(2); /* i1, p1 */

		}

		Matrix_Free(iwork, liwork);
		n = nj;
	}

	UNPROTECT(2); /* perm, ans */
	return ans;
}

/* x[i, j, drop = FALSE] where 'i' and 'j' are vectors of type "integer"   */
/* of length at most 2^31-1 with 'i' in 1:m (or NA) and 'j' in 1:n (or NA) */
/* ... _not_ handling 'Dimnames' here                                      */
SEXP R_subscript_2ary(SEXP s_obj, SEXP s_si, SEXP s_sj)
{
	const char *class = Matrix_class(s_obj, valid_matrix, 6, __func__);
	validObject(s_obj, class);
	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
	case 'p':
		return    dense_subscript_2ary(s_obj, class, s_si, s_sj);
	case 'C':
	case 'R':
	case 'T':
		return   sparse_subscript_2ary(s_obj, class, s_si, s_sj);
	case 'i':
		return diagonal_subscript_2ary(s_obj, class, s_si, s_sj);
	case 'd':
	case 'a':
		return    index_subscript_2ary(s_obj, class, s_si, s_sj);
	default:
		return R_NilValue;
	}
}
