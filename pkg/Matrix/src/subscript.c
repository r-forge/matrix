#include "Mdefines.h"
#include "M5.h"
#include "subscript.h"

#define AR21_UP(i, j, m) i + j + (        j * (    j - 1)) / 2
#define AR21_LO(i, j, m) i +     (j * m + j * (m - j - 1)) / 2

#define SHOW(...) __VA_ARGS__
#define HIDE(...)

#define lCAST(x) (x)
#define iCAST(x) (x)
#define dCAST(x) (x)
#define zCAST(x) (x)

#define iOOB(s, t, mn) (s == NA_INTEGER || (t = s - 1) >= mn)
#define dOOB(s, t, mn) (ISNAN(s) || s >= 0x1.0p+63 || (t = (int_fast64_t) s - 1) >= mn)

#define SUB1(c) \
do { \
	if (TYPEOF(s) == INTSXP) \
		SUB1__(i, c); \
	else \
		SUB1__(d, c); \
} while (0)

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP dense_subscript_1ary(SEXP obj, const char *class, SEXP s)
{
	R_xlen_t slen = XLENGTH(s);
	SEXP ans = allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;

	char ul = '\0', ct = '\0', nu = '\0';
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
		nu = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	R_xlen_t l;
	int_fast64_t b, i, j, k;

	int packed = class[2] == 'p',
		ge = class[1] == 'g', sy = class[1] == 's', he = sy && ct == 'C',
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
					k = PACKED_AR21_UP(i, j); \
				else \
					k = PACKED_AR21_UP(j, i); \
			} else { \
				if (b >= 0) \
					k = PACKED_AR21_LO(i, j, m); \
				else \
					k = PACKED_AR21_LO(j, i, m); \
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
	SEXP ans = allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		obj = sparse_as_Csparse(obj, class);
	}
	PROTECT(obj);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;

	char ul = '\0', ct = '\0', nu = '\0';
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
		nu = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p_ = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i_ = PROTECT(GET_SLOT(obj,        iSym));
	int *pp = INTEGER(p_), *pi = INTEGER(i_), j_, k_, kend_;
	pp++;

	int *po_i = (TYPEOF(o) == INTSXP) ? INTEGER(o) : NULL;
	double *po_d = (TYPEOF(o) == REALSXP) ? REAL(o) : NULL;

	R_xlen_t l = 0, l_;
	int b, i, j, tmp;
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
				if (sy && (b = (up) ? j - i : i - j) < 0) { \
					tmp =   i; \
					i   =   j; \
					j   = tmp; \
				} \
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

	UNPROTECT(4); /* i, p, ans, obj */
	return ans;
}

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP diagonal_subscript_1ary(SEXP obj, const char *class, SEXP s)
{
	R_xlen_t slen = XLENGTH(s);
	SEXP ans = allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];
	int_fast64_t nn = (int_fast64_t) n * n, n1a = (int_fast64_t) n + 1;

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char nu = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	R_xlen_t l;
	int_fast64_t k;

	int un = nu != 'N';

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
	SEXP ans = allocVector(LGLSXP, slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *pa = LOGICAL(ans);

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];
	int_fast64_t mn = (int_fast64_t) m * n;

	SEXP margin = GET_SLOT(obj, Matrix_marginSym);
	int mg = INTEGER(margin)[0] - 1;

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	int *pperm = INTEGER(perm);

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

	UNPROTECT(1); /* ans */
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
	SEXP ans = allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *psi = INTEGER(s), *psj = psi + slen;

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int m = INTEGER(dim)[0];

	char ul = '\0', ct = '\0', nu = '\0';
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
		nu = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int l;
	int_fast64_t b, i, j, k;

	int packed = class[2] == 'p',
		ge = class[1] == 'g', sy = class[1] == 's', he = sy && ct == 'C',
		un = class[1] == 't' && nu != 'N', up = ul == 'U';

#define SUB1(c) \
	do { \
		c##TYPE *pa = c##PTR(ans); \
		c##TYPE *px = c##PTR(x); \
		for (l = 0; l < slen; ++l) { \
			if (psi[l] == NA_INTEGER || psj[l] == NA_INTEGER) { \
				pa[l] = c##NA; \
				continue; \
			} \
			i = (int_fast64_t) psi[l] - 1; \
			j = (int_fast64_t) psj[l] - 1; \
			if (ge) { \
				c##ASSIGN_IDEN(pa[l], c##CAST(px[j * m + i])); \
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
					k = PACKED_AR21_UP(i, j); \
				else \
					k = PACKED_AR21_UP(j, i); \
			} else { \
				if (b >= 0) \
					k = PACKED_AR21_LO(i, j, m); \
				else \
					k = PACKED_AR21_LO(j, i, m); \
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
	SEXP ans = allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *psi = INTEGER(s), *psj = psi + slen;

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		obj = sparse_as_Csparse(obj, class);
	}
	PROTECT(obj);

	char ul = '\0', ct = '\0', nu = '\0';
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
		nu = CHAR(STRING_ELT(diag, 0))[0];
	}

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p_ = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i_ = PROTECT(GET_SLOT(obj,        iSym));
	int *pp = INTEGER(p_), *pi = INTEGER(i_), j_, k_, kend_;
	pp++;

	int *po_i = (TYPEOF(o) == INTSXP) ? INTEGER(o) : NULL;
	double *po_d = (TYPEOF(o) == REALSXP) ? REAL(o) : NULL;

	int l = 0, l_;
	int b, i, j, tmp;

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
			if (psi[l_] == NA_INTEGER || psj[l_] == NA_INTEGER) \
				j = -1; \
			else { \
				if (byrow) { \
					i = psj[l_] - 1; \
					j = psi[l_] - 1; \
				} \
				else { \
					i = psi[l_] - 1; \
					j = psj[l_] - 1; \
				} \
				if (sy && (b = (up) ? j - i : i - j) < 0) { \
					tmp =   i; \
					i   =   j; \
					j   = tmp; \
				} \
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

	UNPROTECT(4); /* i, p, ans, obj */
	return ans;
}

#undef nCAST
#define nCAST(x) (x != 0)

static
SEXP diagonal_subscript_1ary_2col(SEXP obj, const char *class, SEXP s)
{
	int slen = (int) (XLENGTH(s) / 2);
	SEXP ans = allocVector(kindToType(class[0]), slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *psi = INTEGER(s), *psj = psi + slen;

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	char nu = CHAR(STRING_ELT(diag, 0))[0];

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int l;

	int un = nu != 'N';

#define SUB1(c) \
	do { \
		c##TYPE *pa = c##PTR(ans); \
		c##TYPE *px = c##PTR(x); \
		for (l = 0; l < slen; ++l) { \
			if (psi[l] == NA_INTEGER || psj[l] == NA_INTEGER) \
				pa[l] = c##NA; \
			else if (psi[l] != psj[l]) \
				pa[l] = c##ZERO; \
			else if (un) \
				pa[l] = c##UNIT; \
			else \
				pa[l] = c##CAST(px[psi[l] - 1]); \
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
	SEXP ans = allocVector(LGLSXP, slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *psi = INTEGER(s), *psj = psi + slen, *pa = LOGICAL(ans);

	SEXP margin = GET_SLOT(obj, Matrix_marginSym);
	int mg = INTEGER(margin)[0] - 1;

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	int *pperm = INTEGER(perm);

	int l;

	for (l = 0; l < slen; ++l) {
		if (psi[l] == NA_INTEGER || psj[l] == NA_INTEGER)
			pa[l] = NA_LOGICAL;
		else if (mg == 0)
			pa[l] = psj[l] == pperm[psi[l]] - 1;
		else
			pa[l] = psi[l] == pperm[psj[l]] - 1;
	}

	UNPROTECT(1); /* ans */
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
int keep_tr(int *pi, int *pj, int n, int upper, int nonunit, int checkNA)
{
	int k, ident = memcmp(pi, pj, sizeof(int) * (size_t) n) == 0;
	if (checkNA) {
		if (ident) {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER)
					return 0;
		} else {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER || pj[k] == NA_INTEGER)
					return 0;
		}
	}
	int r = (upper) ? 1 : -1;
	if (ident) {
		/* triangular iff monotone; unit diagonal is preserved */
		if (n >= 2) {
			if (pi[0] == pi[1])
				return 0;
			else if (pi[0] < pi[1]) {
				for (k = 2; k < n; ++k)
					if (pi[k - 1] >= pi[k])
						return 0;
				/* up->up, lo->lo */
			} else {
				for (k = 2; k < n; ++k)
					if (pi[k - 1] <= pi[k])
						return 0;
				/* up->lo, lo->up */
				r = -r;
			}
		}
		if (!nonunit)
			r *= 2;
		return r;
	} else {
		/* brute force ... */
		int ki, kj, j_;
		if (upper) {
			for (kj = 0; kj < n; ++kj)
				for (ki = kj + 1, j_ = pj[kj]; ki < n; ++ki)
					if (pi[ki] <= j_)
						goto LO;
			/* up->up */
			return  r;
		LO:
			for (kj = 0; kj < n; ++kj)
				for (ki = 0, j_ = pj[kj]; ki < kj; ++ki)
					if (pi[ki] <= j_)
						return 0;
			/* up->lo */
			return -r;
		} else {
			for (kj = 0; kj < n; ++kj)
				for (ki = 0, j_ = pj[kj]; ki < kj; ++ki)
					if (pi[ki] >= j_)
						goto UP;
			/* lo->lo */
			return  r;
		UP:
			for (kj = 0; kj < n; ++kj)
				for (ki = kj + 1, j_ = pj[kj]; ki < n; ++ki)
					if (pi[ki] >= j_)
						return 0;
			/* lo->up */
			return -r;
		}
	}
}

static
int keep_sy(int *pi, int *pj, int n, int upper, int checkNA)
{
	if (memcmp(pi, pj, sizeof(int) * (size_t) n) != 0)
		return 0;
	int k, r = (upper) ? 1 : -1;
	if (checkNA) {
		for (k = 0; k < n; ++k)
			if (pi[k] == NA_INTEGER)
				return r;
	}
	if (n >= 2) {
		/* triangular iff monotone */
		if (pi[0] == pi[1])
			return r;
		else if (pi[0] < pi[1]) {
			for (k = 2; k < n; ++k)
				if (pi[k - 1] >= pi[k])
					return r;
			/* up->up, lo->lo */
		} else {
			for (k = 2; k < n; ++k)
				if (pi[k - 1] <= pi[k])
					return r;
			/* up->lo, lo->up */
			r = -r;
		}
	}
	return 2 * r;
}

static
int keep_di(int *pi, int *pj, int n, int nonunit, int checkNA, int lwork)
{
	int k, ident = memcmp(pi, pj, sizeof(int) * (size_t) n) == 0;
	if (checkNA) {
		if (ident) {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER)
					return 0;
		} else {
			for (k = 0; k < n; ++k)
				if (pi[k] == NA_INTEGER || pj[k] == NA_INTEGER)
					return 0;
		}
	}
	if (ident) {
		/* diagonal iff no duplicates; unit diagonal is preserved */
		char *work;
		Matrix_Calloc(work, lwork, char);
		--work;
		for (k = 0; k < n; ++k) {
			if (work[pi[k]])
				return 0;
			work[pi[k]] = 1;
		}
		++work;
		Matrix_Free(work, lwork);
		return (nonunit) ? 1 : 2;
	} else {
		/* brute force ... */
		int ki, kj, j_;
		for (kj = 0; kj < n; ++kj) {
			j_ = pj[kj];
			for (ki = 0; ki < kj; ++ki)
				if (pi[ki] == j_)
					return 0;
			for (ki = kj + 1; ki < n; ++ki)
				if (pi[ki] == j_)
					return 0;
		}
		return 1;
	}
}

static
void sort_cr(SEXP obj, const char *cl)
{
	SEXP dim = PROTECT(GET_SLOT(obj, Matrix_DimSym));
	int *pdim = INTEGER(dim),
		m = pdim[(cl[2] == 'C') ? 0 : 1],
		n = pdim[(cl[2] == 'C') ? 1 : 0],
		r = (m < n) ? n : m;
	UNPROTECT(1); /* dim */

	SEXP iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	int *pp = INTEGER(p), *pi = INTEGER(i);

	int i_, j_, k, kend, nnz = pp[n], *workA, *workB, *workC;
	int_fast64_t lwork = (int_fast64_t) m + 1 + r + nnz;
	Matrix_Calloc(workA, lwork, int);
	workB = workA + m + 1;
	workC = workB + r;

#define SORT_LOOP(_MASK_) \
	do { \
		++workA; \
		for (k = 0; k < nnz; ++k) \
			++workA[pi[k]]; \
		--workA; \
		 \
		for (i_ = 1; i_ < m; ++i_) \
			workA[i_] = workB[i_] = workB[i_ - 1] + workA[i_]; \
		workA[m] = nnz; \
		 \
		++pp; \
		k = 0; \
		for (j_ = 0; j_ < n; ++j_) { \
			kend = pp[j_]; \
			while (k < kend) { \
				i_ = pi[k]; \
				workC[workB[i_]] = j_; \
				_MASK_(workD[workB[i_]] = px[k]); \
				++workB[i_]; \
				++k; \
			} \
		} \
		--pp; \
		 \
		for (j_ = 0; j_ < n; ++j_) \
			workB[j_] = pp[j_]; \
		 \
		++workA; \
		k = 0; \
		for (i_ = 0; i_ < m; ++i_) { \
			kend = workA[i_]; \
			while (k < kend) { \
				j_ = workC[k]; \
				pi[workB[j_]] = i_; \
				_MASK_(px[workB[j_]] = workD[k]); \
				++workB[j_]; \
				++k; \
			} \
		} \
		--workA; \
	} while (0)

#define SORT(_CTYPE_, _PTR_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *workD; \
		Matrix_Calloc(workD, nnz,	_CTYPE_); \
		SORT_LOOP(SHOW); \
		Matrix_Free(workD, nnz); \
	} while (0)

	if (cl[0] == 'n')
		SORT_LOOP(HIDE);
	else {
		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		switch (cl[0]) {
		case 'l':
			SORT(int, LOGICAL);
			break;
		case 'i':
			SORT(int, INTEGER);
			break;
		case 'd':
			SORT(double, REAL);
			break;
		case 'z':
			SORT(Rcomplex, COMPLEX);
			break;
		default:
			break;
		}
		UNPROTECT(1); /* x */
	}

#undef SORT_LOOP
#undef SORT

	Matrix_Free(workA, lwork);
	UNPROTECT(2); /* i, p */
	return;
}

#define XIJ_GE(    _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	*(_X_ + _J_ * _M_ + _I_)

#define XIJ_TR_U_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : _ZERO_)

#define XIJ_TR_U_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ < _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TR_L_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : _ZERO_)

#define XIJ_TR_L_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ > _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TP_U_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? *(_X_ + AR21_UP(_I_, _J_, _M_)) \
	 : _ZERO_)

#define XIJ_TP_U_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ < _J_) \
	 ? *(_X_ + AR21_UP(_I_, _J_, _M_)) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_TP_L_N(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? *(_X_ + AR21_LO(_I_, _J_, _M_)) \
	 : _ZERO_)

#define XIJ_TP_L_U(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ > _J_) \
	 ? *(_X_ + AR21_LO(_I_, _J_, _M_)) \
	 : ((_I_ == _J_) ? _ONE_ : _ZERO_))

#define XIJ_SY_U(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : XIJ_GE(_X_, _J_, _I_, _M_, _ZERO_, _ONE_))

#define XIJ_SY_L(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? XIJ_GE(_X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	 : XIJ_GE(_X_, _J_, _I_, _M_, _ZERO_, _ONE_))

#define XIJ_SP_U(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ <= _J_) \
	 ? *(_X_ + AR21_UP(_I_, _J_, _M_)) \
	 : *(_X_ + AR21_UP(_J_, _I_, _M_)))

#define XIJ_SP_L(  _X_, _I_, _J_, _M_, _ZERO_, _ONE_) \
	((_I_ >= _J_) \
	 ? *(_X_ + AR21_LO(_I_, _J_, _M_)) \
	 : *(_X_ + AR21_LO(_J_, _I_, _M_)))

static
SEXP unpackedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, const char *cl)
{

#define SUB2_START \
	SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym)); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	UNPROTECT(1); /* dim */ \
	 \
	int ki, kj, \
		mi = i == R_NilValue, \
		mj = j == R_NilValue, \
		ni = (mi) ? m : LENGTH(i), \
		nj = (mj) ? n : LENGTH(j), \
		*pi = (mi) ? NULL : INTEGER(i), \
		*pj = (mj) ? NULL : INTEGER(j);

#define SUB2_START_EXTRA(_E_, _R_, _Y_, _DENSE_) \
	SUB2_START; \
	 \
	int upper = 1, nonunit = 1, keep = 0; \
	SEXP uplo, diag; \
	if (cl[1] != 'g') { \
		PROTECT(uplo = GET_SLOT(x, Matrix_uploSym)); \
		upper = CHAR(STRING_ELT(uplo, 0))[0] == 'U'; \
		UNPROTECT(1); /* uplo */ \
		if (cl[1] == 't') { \
			PROTECT(diag = GET_SLOT(x, Matrix_diagSym)); \
			nonunit = CHAR(STRING_ELT(diag, 0))[0] == 'N'; \
			UNPROTECT(1); /* diag */ \
		} \
	} \
	 \
	char cl_[] = "...Matrix"; \
	cl_[0] = cl[0]; \
	cl_[1] = 'g'; \
	cl_[2] = _E_; \
	if (cl[1] != 'g' && !(mi || mj) && ni == nj) { \
		if (cl[1] == 't') { \
			keep = keep_tr(pi, pj, ni, upper, nonunit, _DENSE_); \
			if (keep != 0) { \
				cl_[1] = 't'; \
				cl_[2] = _R_; \
			} \
		} else { \
			keep = keep_sy(pi, pj, ni, upper, 0); \
			if ((_DENSE_) ? keep != 0 : keep < -1 || keep > 1) { \
				cl_[1] = 's'; \
				cl_[2] = _Y_; \
			} \
		} \
	} \
	SEXP ans = PROTECT(newObject(cl_)); \
	 \
	PROTECT(dim = GET_SLOT(ans, Matrix_DimSym)); \
	pdim = INTEGER(dim); \
	pdim[0] = ni; \
	pdim[1] = nj; \
	UNPROTECT(1); /* dim */ \
	 \
	if ((cl[1] != 's') ? keep < 0 : keep < -1) { \
		PROTECT(uplo = GET_SLOT(ans, Matrix_uploSym)); \
		SEXP uplo_ = PROTECT(mkChar("L")); \
		SET_STRING_ELT(uplo, 0, uplo_); \
		UNPROTECT(2); /* uplo_, uplo */ \
	} \
	if (cl[1] == 't' && (keep < -1 || keep > 1)) { \
		PROTECT(diag = GET_SLOT(ans, Matrix_diagSym)); \
		SEXP diag_ = PROTECT(mkChar("U")); \
		SET_STRING_ELT(diag, 0, diag_); \
		UNPROTECT(2); /* diag_, diag */ \
	}

	SUB2_START_EXTRA('e', 'r', 'y', 1);

	double ninj = (double) ni * nj;
	if (ninj > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) ninj));

	int i_, j_;
	int_fast64_t m_ = m;

#define SUB2_CASES(_SUB2_) \
	do { \
		switch (cl[0]) { \
		case 'n': \
		case 'l': \
			_SUB2_(int, LOGICAL, NA_LOGICAL, 0, 1); \
			break; \
		case 'i': \
			_SUB2_(int, INTEGER, NA_INTEGER, 0, 1); \
			break; \
		case 'd': \
			_SUB2_(double, REAL, NA_REAL, 0.0, 1.0); \
			break; \
		case 'z': \
			_SUB2_(Rcomplex, COMPLEX, \
			       Matrix_zna, Matrix_zzero, Matrix_zunit); \
			break; \
		default: \
			break; \
		} \
	} while (0)

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (cl_[1] == 'g') { \
			if (cl[1] == 'g') \
				SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
				          XIJ_GE, , , _NA_, _ZERO_, _ONE_); \
			else if (cl[1] == 't') { \
				if (upper) { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_U_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_U_U, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_L_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TR_L_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (upper) \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SY_U, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SY_L, , , _NA_, _ZERO_, _ONE_); \
			} \
		} else if (cl_[1] == 't') { \
			memset(px1, 0, sizeof(_CTYPE_) * (size_t) XLENGTH(x1)); \
			if (upper) { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_U_N, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_U_N, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_U_U, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_U_U, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_L_N, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_L_N, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TR_L_U, , px1 += ni - kj - 1, \
						          _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TR_L_U, px1 += kj, , \
						          _NA_, _ZERO_, _ONE_); \
				} \
			} \
		} else { \
			memset(px1, 0, sizeof(_CTYPE_) * (size_t) XLENGTH(x1)); \
			if (upper) { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SY_U, , px1 += ni - kj - 1, \
					          _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SY_U, px1 += kj, , \
					          _NA_, _ZERO_, _ONE_); \
			} else { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SY_L, , px1 += ni - kj - 1, \
					          _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SY_L, px1 += kj, , \
					          _NA_, _ZERO_, _ONE_); \
			} \
		} \
	} while (0)

#define SUB2_LOOP(_FOR_, _XIJ_, _JUMP1_, _JUMP2_, \
                  _NA_, _ZERO_, _ONE_) \
	do { \
		for (kj = 0; kj < nj; ++kj) { \
			if (mj) \
				j_ = kj; \
			else { \
				j_ = pj[kj]; \
				if (j_ != NA_INTEGER) \
					j_ -= 1; \
				else { \
					_JUMP1_; \
					_FOR_ { \
						*(px1++) = _NA_; \
					} \
					_JUMP2_; \
					continue; \
				} \
			} \
			_JUMP1_; \
			_FOR_ { \
				if (mi) \
					i_ = ki; \
				else { \
					i_ = pi[ki]; \
					if (i_ != NA_INTEGER) \
						i_ -= 1; \
					else { \
						*(px1++) = _NA_; \
						continue; \
					} \
				} \
				*(px1++) = _XIJ_(px0, i_, j_, m_, _ZERO_, _ONE_); \
			} \
			_JUMP2_; \
		} \
	} while (0)

	SUB2_CASES(SUB2);

#undef SUB2

	SET_SLOT(ans, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, ans */
	return ans;
}

static
SEXP packedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, const char *cl)
{
	SUB2_START_EXTRA('e', 'p', 'p', 1);

	double ninj = (double) ni * nj,
	ninj_ = (keep) ? 0.5 * (ninj + ni) : ninj;
	if (ninj_ > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) ninj_));

	int i_, j_;
	int_fast64_t m_ = m;

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		if (cl_[1] == 'g') { \
			if (cl[1] == 't') { \
				if (upper) { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (nonunit) \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
						          XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (upper) \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SP_U, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = 0; ki < ni; ++ki), \
					          XIJ_SP_L, , , _NA_, _ZERO_, _ONE_); \
			} \
		} else if (cl_[1] == 't') { \
			if (upper) { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_U_N, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_U_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} else { \
				if (nonunit) { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_L_N, , , _NA_, _ZERO_, _ONE_); \
				} else { \
					if (keep > 0) \
						SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
						          XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_); \
					else \
						SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
						          XIJ_TP_L_U, , , _NA_, _ZERO_, _ONE_); \
				} \
			} \
		} else { \
			if (upper) { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SP_U, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SP_U, , , _NA_, _ZERO_, _ONE_); \
			} else { \
				if (keep > 0) \
					SUB2_LOOP(for (ki = 0; ki <= kj; ++ki), \
					          XIJ_SP_L, , , _NA_, _ZERO_, _ONE_); \
				else \
					SUB2_LOOP(for (ki = kj; ki < ni; ++ki), \
					          XIJ_SP_L, , , _NA_, _ZERO_, _ONE_); \
			} \
		} \
	} while (0)

	SUB2_CASES(SUB2);

#undef SUB2_LOOP
#undef SUB2

	SET_SLOT(ans, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, ans */
	return ans;
}

static
SEXP CsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, const char *cl)
{
	SUB2_START_EXTRA('C', 'C', 'C', 0);

	if (cl[1] != 'g' && cl_[1] == 'g') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		x = sparse_as_general(x, cl);
	}
	PROTECT(x);

	SEXP p0 = PROTECT(GET_SLOT(x, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(x, Matrix_iSym)),
		p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1)),
		i1 = NULL;
	int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
		*pp1 = INTEGER(p1), *pi1 = NULL,
		d, k, kend, doSort = 0;
	int_fast64_t nnz = 0;
	*(pp1++) = 0;

#define SUB2_FINISH \
	if (nnz > INT_MAX) \
		error(_("%s too dense for %s; would have more than %s nonzero entries"), \
		      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \
	 \
	PROTECT(i1 = allocVector(INTSXP, (R_xlen_t) nnz)); \
	pi1 = INTEGER(i1); \
	 \
	if (cl[0] == 'n') \
		SUB2_LOOP(HIDE); \
	else { \
		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)), \
			x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) nnz)); \
		SUB2_CASES(SUB2); \
		SET_SLOT(ans, Matrix_xSym, x1); \
		UNPROTECT(2); /* x1, x0 */ \
	}

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
	do { \
		_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
		SUB2_LOOP(SHOW); \
	} while (0)

	if (mi) {

		for (kj = 0; kj < nj; ++kj) {
			nnz += pp0[pj[kj]] - pp0[pj[kj] - 1];
			pp1[kj] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (kj = 0; kj < nj; ++kj) { \
				k = pp0[pj[kj] - 1]; \
				kend = pp0[pj[kj]]; \
				d = kend - k; \
				if (d) { \
					memcpy(pi1, pi0 + k, sizeof(int) * (size_t) d); \
					pi1 += d; \
					_MASK_(memcpy(px1, px0 + k, sizeof(*px1) * (size_t) d)); \
					_MASK_(px1 += d); \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

	} else {

		int *workA, *workB, *workC;
		int_fast64_t lwork = (int_fast64_t) m + m + ni;
		Matrix_Calloc(workA, lwork, int);
		workB = workA + m;
		workC = workB + m;

		/* workA[ i] : size of the set { ki : pi[ki] - 1 == i }
		   workB[ i] : smallest ki such that pi[ki] - 1 == i
		   workC[ki] : smallest ki' > ki such that pi[ki'] == pi[ki]
		*/

		int i_, j_, i_prev = m;

		for (ki = ni - 1; ki >= 0; --ki) {
			i_ = pi[ki] - 1;
			++workA[i_];
			workC[ki] = workB[i_];
			workB[i_] = ki;
			if (i_ > i_prev)
				doSort = 1;
			i_prev = i_;
		}

		for (kj = 0; kj < nj; ++kj) {
			j_ = (mj) ? kj : pj[kj] - 1;
			k = pp0[j_];
			kend = pp0[j_ + 1];
			while (k < kend) {
				nnz += workA[pi0[k]];
				++k;
			}
			pp1[kj] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (kj = 0; kj < nj; ++kj) { \
				j_ = (mj) ? kj : pj[kj] - 1; \
				k = pp0[j_]; \
				kend = pp0[j_ + 1]; \
				while (k < kend) { \
					i_ = pi0[k]; \
					d = workA[i_]; \
					ki = workB[i_]; \
					while (d--) { \
						*(pi1++) = ki; \
						_MASK_(*(px1++) = px0[k]); \
						ki = workC[ki]; \
					} \
					++k; \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

		Matrix_Free(workA, lwork);

	}

#undef SUB2_FINISH

	SET_SLOT(ans, Matrix_pSym, p1);
	SET_SLOT(ans, Matrix_iSym, i1);
	UNPROTECT(4); /* i1, p1, i0, p0 */

	if (doSort)
		sort_cr(ans, cl);
	if (cl[1] == 's' && (keep == -1 || keep == 1)) {
		/* defined in ./sparse.c : */
		SEXP sparse_force_symmetric(SEXP, const char *, char);
		ans = sparse_force_symmetric(ans, cl_, (keep == 1) ? 'U' : 'L');
	}
	UNPROTECT(2); /* x, ans */
	return ans;
}

static
SEXP RsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, const char *cl)
{
	SUB2_START_EXTRA('R', 'R', 'R', 0);

	if (cl[1] != 'g' && cl_[1] == 'g') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		x = sparse_as_general(x, cl);
	}
	PROTECT(x);

	SEXP p0 = PROTECT(GET_SLOT(x, Matrix_pSym)),
		j0 = PROTECT(GET_SLOT(x, Matrix_jSym)),
		p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) ni + 1)),
		j1 = NULL;
	int *pp0 = INTEGER(p0), *pj0 = INTEGER(j0),
		*pp1 = INTEGER(p1), *pj1 = NULL,
		d, k, kend, doSort = 0;
	int_fast64_t nnz = 0;
	*(pp1++) = 0;

#define SUB2_FINISH \
	if (nnz > INT_MAX) \
		error(_("%s too dense for %s; would have more than %s nonzero entries"), \
		      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \
	 \
	PROTECT(j1 = allocVector(INTSXP, (R_xlen_t) nnz)); \
	pj1 = INTEGER(j1); \
	 \
	if (cl[0] == 'n') \
		SUB2_LOOP(HIDE); \
	else { \
		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)), \
			x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) nnz)); \
		SUB2_CASES(SUB2); \
		SET_SLOT(ans, Matrix_xSym, x1); \
		UNPROTECT(2); /* x1, x0 */ \
	}

	if (mj) {

		for (ki = 0; ki < ni; ++ki) {
			nnz += pp0[pi[ki]] - pp0[pi[ki] - 1];
			pp1[ki] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (ki = 0; ki < ni; ++ki) { \
				k = pp0[pi[ki] - 1]; \
				kend = pp0[pi[ki]]; \
				d = kend - k; \
				if (d) { \
					memcpy(pj1, pj0 + k, sizeof(int) * (size_t) d); \
					pj1 += d; \
					_MASK_(memcpy(px1, px0 + k, sizeof(*px1) * (size_t) d)); \
					_MASK_(px1 += d); \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

	} else {

		int *workA, *workB, *workC;
		int_fast64_t lwork = (int_fast64_t) n + n + nj;
		Matrix_Calloc(workA, lwork, int);
		workB = workA + n;
		workC = workB + n;

		/* workA[ j] : size of the set { kj : pj[kj] - 1 == j }
		   workB[ j] : smallest ki such that pj[kj] - 1 == j
		   workC[kj] : smallest kj' > kj such that pj[kj'] == pj[kj]
		*/

		int i_, j_, j_prev = n;

		for (kj = nj - 1; kj >= 0; --kj) {
			j_ = pj[kj] - 1;
			++workA[j_];
			workC[kj] = workB[j_];
			workB[j_] = kj;
			if (j_ > j_prev)
				doSort = 1;
			j_prev = j_;
		}

		for (ki = 0; ki < ni; ++ki) {
			i_ = (mi) ? ki : pi[ki] - 1;
			k = pp0[i_];
			kend = pp0[i_ + 1];
			while (k < kend) {
				nnz += workA[pj0[k]];
				++k;
			}
			pp1[ki] = (int) nnz;
		}

#define SUB2_LOOP(_MASK_) \
		do { \
			for (ki = 0; ki < ni; ++ki) { \
				i_ = (mi) ? ki : pi[ki] - 1; \
				k = pp0[i_]; \
				kend = pp0[i_ + 1]; \
				while (k < kend) { \
					j_ = pj0[k]; \
					d = workA[j_]; \
					kj = workB[j_]; \
					while (d--) { \
						*(pj1++) = kj; \
						_MASK_(*(px1++) = px0[k]); \
						kj = workC[kj]; \
					} \
					++k; \
				} \
			} \
		} while (0)

		SUB2_FINISH;

#undef SUB2_LOOP

		Matrix_Free(workA, lwork);

	}

#undef SUB2_FINISH
#undef SUB2

	SET_SLOT(ans, Matrix_pSym, p1);
	SET_SLOT(ans, Matrix_jSym, j1);
	UNPROTECT(4); /* j1, p1, j0, p0 */

	if (doSort)
		sort_cr(ans, cl);
	if (cl[1] == 's' && (keep == -1 || keep == 1)) {
		/* defined in ./sparse.c : */
		SEXP sparse_force_symmetric(SEXP, const char *, char);
		ans = sparse_force_symmetric(ans, cl_, (keep == 1) ? 'U' : 'L');
	}
	UNPROTECT(2); /* x, ans */
	return ans;
}

static
SEXP diagonalMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, const char *cl)
{
	SUB2_START;

	int nonunit = 1, keep = 0;
	SEXP diag = PROTECT(GET_SLOT(x, Matrix_diagSym));
	nonunit = CHAR(STRING_ELT(diag, 0))[0] == 'N';
	UNPROTECT(1); /* diag */

	char cl_[] = ".gCMatrix";
	cl_[0] = cl[0];
	if (!(mi || mj) && ni == nj) {
		keep = keep_di(pi, pj, ni, nonunit, 0, m);
		if (keep > 0) {
			cl_[1] = 'd';
			cl_[2] = 'i';
		}
	}

	SEXP ans = PROTECT(newObject(cl_));

	PROTECT(dim = GET_SLOT(ans, Matrix_DimSym));
	pdim = INTEGER(dim);
	pdim[0] = (int) ni;
	pdim[1] = (int) nj;
	UNPROTECT(1); /* dim */

	if (keep > 1) {

		PROTECT(diag = GET_SLOT(ans, Matrix_diagSym));
		SEXP diag_ = PROTECT(mkChar("U"));
		SET_STRING_ELT(diag, 0, diag_);
		UNPROTECT(2); /* diag_, diag */

	} else if (keep > 0) {

		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym)),
			x1 = PROTECT(allocVector(TYPEOF(x0), ni));
		int j_;

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			while (ni--) \
				*(px1++) = \
					(*(pi++) != (j_ = *(pj++))) \
					? _ZERO_ \
					: ((nonunit) ? px0[j_ - 1] : _ONE_); \
		} while (0)

		SUB2_CASES(SUB2);

#undef SUB2

		SET_SLOT(ans, Matrix_xSym, x1);
		UNPROTECT(2); /* x0, x1 */

	} else {

		SEXP x0 = PROTECT(GET_SLOT(x, Matrix_xSym));
		char *work = NULL;
		int j_;

		if (nonunit) {

			Matrix_Calloc(work, n, char);

#define SUB2_WORK(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0); \
				for (j_ = 0; j_ < n; ++j_) \
					work[j_] = c##NOT_ZERO(px0[j_]); \
			} while (0)

			SWITCH5(cl[0], SUB2_WORK);

#undef SUB2_WORK

		}

		SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
		int *pp1 = INTEGER(p1);
		*(pp1++) = 0;

		for (kj = 0; kj < nj; ++kj) {
			pp1[kj] = 0;
			j_ = (mj) ? kj : pj[kj] - 1;
			if (!nonunit || work[j_]) {
				if (mi) {
					for (ki = 0; ki < ni; ++ki)
						if (ki == j_)
							++pp1[kj];
				} else {
					for (ki = 0; ki < ni; ++ki)
						if (pi[ki] - 1 == j_)
							++pp1[kj];
				}
				if (pp1[kj] > INT_MAX - pp1[kj - 1]) {
					if (nonunit)
						Matrix_Free(work, n);
					error(_("%s too dense for %s; would have more than %s nonzero entries"),
					      "x[i,j]", "[CR]sparseMatrix", "2^31-1");
				}
			}
			pp1[kj] += pp1[kj-1];
		}

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[nj - 1])),
			x1 = PROTECT(allocVector(TYPEOF(x0), pp1[nj - 1]));
		int *pi1 = INTEGER(i1);

#define SUB2(_CTYPE_, _PTR_, _NA_, _ZERO_, _ONE_) \
		do { \
			_CTYPE_ *px0 = _PTR_(x0), *px1 = _PTR_(x1); \
			for (kj = 0; kj < nj; ++kj) { \
				j_ = (mj) ? kj : pj[kj] - 1; \
				if (!nonunit || work[j_]) { \
					for (ki = 0; ki < ni; ++ki) { \
						if (((mi) ? ki : pi[ki] - 1) == j_) { \
							*(pi1++) = ki; \
							*(px1++) = (nonunit) ? px0[j_] : _ONE_; \
						} \
					} \
				} \
			} \
		} while (0)

		SUB2_CASES(SUB2);

#undef SUB2

		SET_SLOT(ans, Matrix_pSym, p1);
		SET_SLOT(ans, Matrix_iSym, i1);
		SET_SLOT(ans, Matrix_xSym, x1);
		UNPROTECT(4); /* x1, x0, i1, p1 */

		if (nonunit)
			Matrix_Free(work, n);

	}

	UNPROTECT(1); /* ans */
	return ans;
}

static
SEXP indMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, const char *cl)
{
	PROTECT_INDEX pidA;
	PROTECT_WITH_INDEX(x, &pidA);

	PROTECT_INDEX pidB;
	SEXP perm0 = GET_SLOT(x, Matrix_permSym);
	int *pperm0 = INTEGER(perm0);
	PROTECT_WITH_INDEX(perm0, &pidB);

	SEXP margin = PROTECT(GET_SLOT(x, Matrix_marginSym));
	int mg = INTEGER(margin)[0] - 1;

	SEXP dim = PROTECT(GET_SLOT(x, Matrix_DimSym));
	int *pdim = INTEGER(dim), m = pdim[mg], n = pdim[!mg];
	UNPROTECT(1); /* dim */

	if (mg) {
		SEXP i_tmp = i;
		i = j;
		j = i_tmp;
	}

	int ki, kj,
		mi = i == R_NilValue,
		mj = j == R_NilValue,
		ni = (mi) ? m : LENGTH(i),
		nj = (mj) ? n : LENGTH(j),
		*pi = (mi) ? NULL : INTEGER(i),
		*pj = (mj) ? NULL : INTEGER(j),
		isP = cl[0] == 'p';

	if (!mi) {
		isP = isP && ni == m;
		if (isP) {
			char *work;
			Matrix_Calloc(work, m, char);
			--work; /* now 1-indexed */
			for (ki = 0; ki < ni; ++ki) {
				if (work[pi[ki]]) {
					isP = 0;
					break;
				}
				work[pi[ki]] = 1;
			}
			++work; /* now 0-indexed */
			Matrix_Free(work, m);
		}

		x = newObject((isP) ? "pMatrix" : "indMatrix");
		REPROTECT(x, pidA);

		PROTECT(dim = GET_SLOT(x, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[ mg] = ni;
		pdim[!mg] = n;
		UNPROTECT(1); /* dim */

		if (mg)
			SET_SLOT(x, Matrix_marginSym, margin);

		SEXP perm1 = PROTECT(allocVector(INTSXP, ni));
		int *pperm1 = INTEGER(perm1);
		--pperm0; /* now 1-indexed */
		for (ki = 0; ki < ni; ++ki)
			pperm1[ki] = pperm0[pi[ki]];
		SET_SLOT(x, Matrix_permSym, perm1);
		UNPROTECT(1); /* perm1 */

		perm0 = perm1;
		pperm0 = pperm1;
		REPROTECT(perm0, pidB);

		m = ni;
	}

	if (!mj) {
		isP = isP && nj == n;
		if (isP) {
			char *work;
			Matrix_Calloc(work, nj, char);
			--work; /* now 1-indexed */
			for (kj = 0; kj < nj; ++kj) {
				if (work[pj[kj]]) {
					isP = 0;
					break;
				}
				work[pj[kj]] = 1;
			}
			++work; /* now 0-indexed */
			Matrix_Free(work, nj);
		}

		x = newObject((isP)
		                        ? "pMatrix"
		                        : ((!mg) ? "ngCMatrix" : "ngRMatrix"));
		REPROTECT(x, pidA);

		PROTECT(dim = GET_SLOT(x, Matrix_DimSym));
		pdim = INTEGER(dim);
		pdim[ mg] = m;
		pdim[!mg] = nj;
		UNPROTECT(1); /* dim */

		if (isP) {
			SEXP perm1 = PROTECT(allocVector(INTSXP, nj));
			int *pperm1 = INTEGER(perm1), *work;
			Matrix_Calloc(work, nj, int);
			--work; /* now 1-indexed */
			for (kj = 0; kj < nj; ++kj)
				work[pj[kj]] = kj + 1;
			for (kj = 0; kj < nj; ++kj)
				pperm1[kj] = work[pperm0[kj]];
			++work; /* now 0-indexed */
			Matrix_Free(work, nj);
			SET_SLOT(x, Matrix_permSym, perm1);
			UNPROTECT(1); /* perm1 */
		} else {
			int *workA, *workB, *workC;
			int_fast64_t lwork = (int_fast64_t) n + n + m;
			Matrix_Calloc(workA, lwork, int);
			workB = workA + n;
			workC = workB + n;
			--workA; /* now 1-indexed */
			--workB; /* now 1-indexed */

			SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
			int *pp1 = INTEGER(p1), k, kend;

			/* 1. Compute old column counts in 'workA' */
			for (k = 0; k < m; ++k)
				++workA[pperm0[k]];

			/* 2. Compute new column pointers in 'pp1' */
			*(pp1++) = 0;
			for (kj = 0; kj < nj; ++kj) {
				pp1[kj] = workA[pj[kj]];
				if (pp1[kj] > INT_MAX - pp1[kj - 1])
					error(_("%s too dense for %s; would have more than %s nonzero entries"), \
					      "x[i,j]", "[CR]sparseMatrix", "2^31-1"); \

				pp1[kj] += pp1[kj - 1];
			}

			/* 3. Compute old column pointers in 'workB' and copy to 'workA' */
			workB[1] = 0;
			for (k = 1; k < n; ++k) {
				workB[k + 1] = workB[k] + workA[k];
				workA[k] = workB[k];
			}
			workA[n] = workB[n];

			/* 4. Sort old row indices into 'workC' */
			for (k = 0; k < m; ++k)
				workC[workA[pperm0[k]]++] = k;

			SEXP i1 = PROTECT(allocVector(INTSXP, pp1[nj - 1]));
			int *pi1 = INTEGER(i1), pos;

			/* 5. Copy row indices from 'workC' to 'pi1' */
			k = 0;
			for (kj = 0; kj < nj; ++kj) {
				kend = pp1[kj];
				pos = workB[pj[kj]];
				while (k < kend)
					pi1[k++] = workC[pos++];
			}

			++workA; /* now 0-indexed */
			Matrix_Free(workA, lwork);
			SET_SLOT(x, Matrix_pSym, p1);
			SET_SLOT(x, (!mg) ? Matrix_iSym : Matrix_jSym, i1);
			UNPROTECT(2); /* i1, p1 */
		}

		n = nj;
	}

#undef SUB2_CASES
#undef SUB2_START
#undef SUB2_START_EXTRA

	UNPROTECT(3); /* margin, perm0, x */
	return x;
}

/* x[i,j,drop=FALSE] with 'i' and 'j' of type "integer" and length
   not exceeding 2^31-1 {'i' in 1:m or NA, 'j' in 1:n or NA} ...
   but _not_ handling 'Dimnames'
*/
SEXP R_subscript_2ary(SEXP s_x, SEXP s_i, SEXP s_j)
{
	if (s_i == R_NilValue && s_j == R_NilValue)
		return s_x;

	const char *cl = Matrix_class(s_x, valid_matrix, 7, __func__);
	validObject(s_x, cl);

	switch (cl[2]) {
	case 'e':
	case 'y':
	case 'r':
		return unpackedMatrix_subscript_2ary(s_x, s_i, s_j, cl);
	case 'p':
		return   packedMatrix_subscript_2ary(s_x, s_i, s_j, cl);
	default:
		break;
	}

	static SEXP anyNA = NULL;
	if (!anyNA)
		anyNA = install("anyNA");
	SEXP call = PROTECT(lang2(anyNA, R_NilValue)), value;

#define ERROR_IF_ANYNA(_I_) \
	do { \
		if ((_I_) != R_NilValue) { \
			SETCADR(call, _I_); \
			PROTECT(value = eval(call, R_BaseEnv)); \
			if (asLogical(value)) \
				error(_("NA subscripts in %s not supported for '%s' inheriting from %s"), \
					  "x[i,j]", "x", "sparseMatrix"); \
			UNPROTECT(1); \
		} \
	} while (0)

	ERROR_IF_ANYNA(s_i);
	ERROR_IF_ANYNA(s_j);

#undef ERROR_IF_ANYNA

	UNPROTECT(1);

	switch (cl[2]) {
	case 'C':
		return  CsparseMatrix_subscript_2ary(s_x, s_i, s_j, cl);
	case 'R':
		return  RsparseMatrix_subscript_2ary(s_x, s_i, s_j, cl);
	case 'T':
	{
		char cl_[] = "..CMatrix";
		cl_[0] = cl[0];
		cl_[1] = cl[1];

		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		SEXP sparse_as_Tsparse(SEXP, const char *);

		s_x = sparse_as_Csparse(s_x, cl);
		PROTECT(s_x);
		s_x = CsparseMatrix_subscript_2ary(s_x, s_i, s_j, cl_);
		UNPROTECT(1);
		PROTECT(s_x);
		s_x = sparse_as_Tsparse(s_x, Matrix_class(s_x, valid_matrix, 7, __func__));
		UNPROTECT(1);
		return s_x;
	}
	case 'i':
		return diagonalMatrix_subscript_2ary(s_x, s_i, s_j, cl);
	default:
		return      indMatrix_subscript_2ary(s_x, s_i, s_j, cl);
	}
}
