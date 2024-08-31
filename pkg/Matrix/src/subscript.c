#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "subscript.h"

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

	UNPROTECT(4); /* i, p, obj, ans */
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
	int un = CHAR(STRING_ELT(diag, 0))[0] != 'N';

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
	int *ps0 = INTEGER(s), *ps1 = ps0 + slen;

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
	SEXP ans = allocVector(kindToType(class[0]), slen);
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

	UNPROTECT(4); /* i, p, obj, ans */
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
	int *ps0 = INTEGER(s), *ps1 = ps0 + slen;

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	int un = CHAR(STRING_ELT(diag, 0))[0] != 'N';

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
	SEXP ans = allocVector(LGLSXP, slen);
	if (slen == 0)
		return ans;
	PROTECT(ans);
	int *ps0 = INTEGER(s), *ps1 = ps0 + slen, *pa = LOGICAL(ans);

	SEXP margin = GET_SLOT(obj, Matrix_marginSym);
	int mg = INTEGER(margin)[0] - 1;

	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	int *pperm = INTEGER(perm);

	int l;

	for (l = 0; l < slen; ++l) {
		if (ps0[l] == NA_INTEGER || ps1[l] == NA_INTEGER)
			pa[l] = NA_LOGICAL;
		else if (mg == 0)
			pa[l] = ps1[l] == pperm[ps0[l]] - 1;
		else
			pa[l] = ps0[l] == pperm[ps1[l]] - 1;
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

#define XIJ_GE(c, x, i, j, m) \
	*(x + DENSE_INDEX_N(i, j, m))

#define XIJ_SY_U(c, x, i, j, m) \
	((i <= j) \
	 ? XIJ_GE(c, x, i, j, m) \
	 : XIJ_GE(c, x, j, i, m))

#define XIJ_SY_L(c, x, i, j, m) \
	((i >= j) \
	 ? XIJ_GE(c, x, i, j, m) \
	 : XIJ_GE(c, x, j, i, m))

#define XIJ_SP_U(c, x, i, j, m) \
	((i <= j) \
	 ? *(x + DENSE_INDEX_U(i, j, m)) \
	 : *(x + DENSE_INDEX_U(j, i, m)))

#define XIJ_SP_L(c, x, i, j, m) \
	((i >= j) \
	 ? *(x + DENSE_INDEX_L(i, j, m)) \
	 : *(x + DENSE_INDEX_L(j, i, m)))

#define XIJ_TR_U_N(c, x, i, j, m) \
	((i <= j) \
	 ? XIJ_GE(c, x, i, j, m) \
	 : c##ZERO)

#define XIJ_TR_L_N(c, x, i, j, m) \
	((i >= j) \
	 ? XIJ_GE(c, x, i, j, m) \
	 : c##ZERO)

#define XIJ_TP_U_N(c, x, i, j, m) \
	((i <= j) \
	 ? *(x + DENSE_INDEX_U(i, j, m)) \
	 : c##ZERO)

#define XIJ_TP_L_N(c, x, i, j, m) \
	((i >= j) \
	 ? *(x + DENSE_INDEX_L(i, j, m)) \
	 : c##ZERO)

#define XIJ_TR_U_U(c, x, i, j, m) \
	((i < j) \
	 ? XIJ_GE(c, x, i, j, m) \
	 : ((i == j) ? c##UNIT : c##ZERO))

#define XIJ_TR_L_U(c, x, i, j, m) \
	((i > j) \
	 ? XIJ_GE(c, x, i, j, m) \
	 : ((i == j) ? c##UNIT : c##ZERO))

#define XIJ_TP_U_U(c, x, i, j, m) \
	((i < j) \
	 ? *(x + DENSE_INDEX_U(i, j, m)) \
	 : ((i == j) ? c##UNIT : c##ZERO))

#define XIJ_TP_L_U(c, x, i, j, m) \
	((i > j) \
	 ? *(x + DENSE_INDEX_L(i, j, m)) \
	 : ((i == j) ? c##UNIT : c##ZERO))

static
SEXP dense_subscript_2ary(SEXP obj, const char *class, SEXP si, SEXP sj)
{
	if (si == R_NilValue && sj == R_NilValue)
		return obj;

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1];

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

	int ki, kj,
		mi = si == R_NilValue,
		mj = sj == R_NilValue,
		ni = (mi) ? m : LENGTH(si),
		nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);

	int stay = (class[1] == 'g') ? 0 : (class[1] == 's') ? stay_sy(pi, ni, pj, nj, n, ul, 1) : stay_tr(pi, ni, pj, nj, n, ul, 1),
		packed = class[2] == 'p';
	int_fast64_t ninj = (int_fast64_t) ni * nj,
		xlen = (!packed || stay == 0) ? ninj : ni + (ninj - ni) / 2;
	if (xlen > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (stay == 0) ? 'g' : class[1];
	cl[2] = (stay == 0) ? 'e' : class[2];
	SEXP ans = PROTECT(newObject(cl));

	dim = GET_SLOT(ans, Matrix_DimSym);
	INTEGER(dim)[0] = ni;
	INTEGER(dim)[1] = nj;

	if (cl[1] != 'g' && stay < 0) {
		SEXP uplo = GET_SLOT(ans, Matrix_uploSym);
		SET_STRING_ELT(uplo, 0, mkChar("L"));
	}
	if (cl[1] == 's' && cl[0] == 'z' && ct != 'C') {
		SEXP trans = GET_SLOT(ans, Matrix_transSym);
		SET_STRING_ELT(trans, 0, mkChar("T"));
	}
	if (cl[1] == 't' && nu != 'N' && (stay < -1 || stay > 1)) {
		SEXP diag = GET_SLOT(ans, Matrix_diagSym);
		SET_STRING_ELT(diag, 0, mkChar("U"));
	}

	SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		x1 = PROTECT(allocVector(TYPEOF(x0), (R_xlen_t) xlen));
	SET_SLOT(ans, Matrix_xSym, x1);

	int_fast64_t i, j;

#define SUB2(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed && stay != 0) \
			memset(px1, 0, sizeof(c##TYPE) * (size_t) XLENGTH(x1)); \
		if (class[1] == 'g') { \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_GE, , ); \
		} else if (class[1] == 's' && !packed) { \
			if (stay == 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_SY_U, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_SY_L, , ); \
			} else if (stay > 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_SY_U, \
			       , px1 += ni - kj - 1); \
			else \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_SY_L, \
			       , px1 += ni - kj - 1); \
			} else { \
			if (ul == 'U') \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_SY_U, \
			       px1 += kj, ); \
			else \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_SY_L, \
			       px1 += kj, ); \
			} \
		} else if (class[1] == 's') { \
			if (stay == 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_SP_U, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_SP_L, , ); \
			} else if (stay > 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_SP_U, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_SP_L, , ); \
			} else { \
			if (ul == 'U') \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_SP_U, , ); \
			else \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_SP_L, , ); \
			} \
		} else if (nu == 'N' && !packed) { \
			if (stay == 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TR_U_N, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TR_L_N, , ); \
			} else if (stay > 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TR_U_N, \
			       , px1 += ni - kj - 1); \
			else \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TR_L_N, \
			       , px1 += ni - kj - 1); \
			} else { \
			if (ul == 'U') \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TR_U_N, \
			       px1 += kj, ); \
			else \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TR_L_N, \
			       px1 += kj, ); \
			} \
		} else if (nu == 'N') { \
			if (stay == 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TP_U_N, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TP_L_N, , ); \
			} else if (stay > 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TP_U_N, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TP_L_N, , ); \
			} else { \
			if (ul == 'U') \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TP_U_N, , ); \
			else \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TP_L_N, , ); \
			} \
		} else if (!packed) { \
			if (stay == 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TR_U_U, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TR_L_U, , ); \
			} else if (stay > 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TR_U_U, \
			       , px1 += ni - kj - 1); \
			else \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TR_L_U, \
			       , px1 += ni - kj - 1); \
			} else { \
			if (ul == 'U') \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TR_U_U, \
			       px1 += kj, ); \
			else \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TR_L_U, \
			       px1 += kj, ); \
			} \
		} else { \
			if (stay == 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TP_U_U, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <  ni; ++ki), XIJ_TP_L_U, , ); \
			} else if (stay > 0) { \
			if (ul == 'U') \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TP_U_U, , ); \
			else \
			SUB2__(c, for (ki =  0; ki <= kj; ++ki), XIJ_TP_L_U, , ); \
			} else { \
			if (ul == 'U') \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TP_U_U, , ); \
			else \
			SUB2__(c, for (ki = kj; ki <  ni; ++ki), XIJ_TP_L_U, , ); \
			} \
		} \
	} while (0)

#define SUB2__(c, __for__, index, jump0, jump1) \
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
				*(px1++) = index(c, px0, i, j, m); \
			} \
			jump1; \
		} \
	} while (0)

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
	if (!mg) {
		SEXP tmp;
		tmp = si; si = sj; sj = tmp;
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[!mg], n = pdim[mg];

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g') {
		SEXP uplo = GET_SLOT(obj, Matrix_uploSym);
		ul = CHAR(STRING_ELT(uplo, 0))[0];
		if (!mg)
			ul = (ul == 'U') ? 'L' : 'U';
	}
	if (class[1] == 's' && class[0] == 'z') {
		SEXP trans = GET_SLOT(obj, Matrix_transSym);
		ct = CHAR(STRING_ELT(trans, 0))[0];
	}
	if (class[1] == 't') {
		SEXP diag = GET_SLOT(obj, Matrix_diagSym);
		nu = CHAR(STRING_ELT(diag, 0))[0];
	}

	int ki, kj,
		mi = si == R_NilValue,
		mj = sj == R_NilValue,
		ni = (mi) ? m : LENGTH(si),
		nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);

	if (anyNA(pi, ni) || anyNA(pj, nj))
		error(_("NA subscripts in %s not supported for '%s' inheriting from %s"),
		      "x[i, j]", "x", "sparseMatrix");
	int stay = (class[1] == 'g') ? 0 : (class[1] == 's') ? stay_sy(pi, ni, pj, nj, n, ul, 0) : stay_tr(pi, ni, pj, nj, n, ul, 0);
	if (!mg)
		stay = -stay;

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
	if (class[1] != 'g' && stay >= -1 && stay <= 1) {
		/* defined in ./coerce.c : */
		SEXP sparse_as_general(SEXP, const char *);
		REPROTECT(obj = sparse_as_general(obj, class__), pid_obj);
		class__[1] = 'g';
	}

	char cl[] = "...Matrix";
	cl[0] = class  [0];
	cl[1] = class  [1];
	cl[2] = class__[2];
	if ((class[1] == 's' && stay >= -1 && stay <= 1) ||
	    (class[1] == 't' && stay == 0))
		cl[1] = 'g';
	PROTECT_INDEX pid_ans;
	SEXP ans;
	PROTECT_WITH_INDEX(ans = newObject(cl), &pid_ans);

	dim = GET_SLOT(ans, Matrix_DimSym);
	INTEGER(dim)[!mg] = ni;
	INTEGER(dim)[ mg] = nj;

	if (cl[1] != 'g' && stay < 0) {
		SEXP uplo = GET_SLOT(ans, Matrix_uploSym);
		SET_STRING_ELT(uplo, 0, mkChar("L"));
	}
	if (cl[1] == 's' && cl[0] == 'z' && ct != 'C') {
		SEXP trans = GET_SLOT(ans, Matrix_transSym);
		SET_STRING_ELT(trans, 0, mkChar("T"));
	}
	if (cl[1] == 't' && nu != 'N' && (stay < -1 || stay > 1)) {
		SEXP diag = GET_SLOT(ans, Matrix_diagSym);
		SET_STRING_ELT(diag, 0, mkChar("U"));
	}

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p0 = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i0 = PROTECT(GET_SLOT(obj,        iSym)),
		p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
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
			error(_("%s too dense for %s; would have more than %s nonzero entries"),
			      "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, (int) nnz));
		int *pi1 = INTEGER(i1);
		SET_SLOT(ans, iSym, i1);

#define SUB2(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, (int) nnz)); \
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
			error(_("%s too dense for %s; would have more than %s nonzero entries"),
			      "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, (int) nnz));
		int *pi1 = INTEGER(i1), d;
		SET_SLOT(ans, iSym, i1);

#define SUB2(c) \
		do { \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)), \
				x1 = PROTECT(allocVector(c##TYPESXP, (int) nnz)); \
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

#define SORT(c) \
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

			SWITCH5(class[0], SORT);

#undef SORT

			Matrix_Free(iwork, liwork);
		}

		UNPROTECT(1); /* i1 */

	}

	if (class[1] == 's' && cl[1] == 'g' && stay != 0) {
		/* defined in ./sparse.c : */
		SEXP sparse_force_symmetric(SEXP, const char *, char, char);
		REPROTECT(ans = sparse_force_symmetric(ans, cl, (stay > 0) ? 'U' : 'L', ct), pid_ans);
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

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int n = INTEGER(dim)[0];

	SEXP diag = GET_SLOT(obj, Matrix_diagSym);
	int un = CHAR(STRING_ELT(diag, 0))[0] != 'N';

	int ki, kj,
		mi = si == R_NilValue,
		mj = sj == R_NilValue,
		ni = (mi) ? n : LENGTH(si),
		nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);

	if (anyNA(pi, ni) || anyNA(pj, nj))
		error(_("NA subscripts in %s not supported for '%s' inheriting from %s"),
		      "x[i, j]", "x", "sparseMatrix");
	int stay = stay_di(pi, ni, pj, nj, n, 0);

	char cl[] = "...Matrix";
	cl[0] = class[0];
	cl[1] = (stay > 0) ? 'd' : 'g';
	cl[2] = (stay > 0) ? 'i' : 'C';
	SEXP ans = PROTECT(newObject(cl));

	dim = GET_SLOT(ans, Matrix_DimSym);
	INTEGER(dim)[0] = ni;
	INTEGER(dim)[1] = nj;

	if (un && stay > 1) {

		diag = GET_SLOT(ans, Matrix_diagSym);
		SET_STRING_ELT(diag, 0, mkChar("U"));

	} else if (stay > 0) {

		SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)),
			x1 = PROTECT(allocVector(TYPEOF(x0), ni));
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

		SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
		int *pp1 = INTEGER(p1), j;
		*(pp1++) = 0;
		SET_SLOT(ans, Matrix_pSym, p1);

#define COUNT(c) \
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
		SWITCH4(class[0], COUNT);

#undef COUNT

		if (kj < nj)
			error(_("%s too dense for %s; would have more than %s nonzero entries"),
			      "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[nj - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(ans, Matrix_iSym, i1);

#define SUB2(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0); \
			c##IF_NPATTERN( \
			SEXP x1 = PROTECT(allocVector(c##TYPESXP, pp1[nj - 1])); \
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

	PROTECT_INDEX pid_perm;
	SEXP perm = GET_SLOT(obj, Matrix_permSym);
	int *pperm = INTEGER(perm);
	PROTECT_WITH_INDEX(perm, &pid_perm);

	SEXP margin = GET_SLOT(obj, Matrix_marginSym);
	int mg = INTEGER(margin)[0] - 1;
	if (mg) {
		SEXP tmp;
		tmp = si; si = sj; sj = tmp;
	}

	SEXP dim = GET_SLOT(obj, Matrix_DimSym);
	int *pdim = INTEGER(dim), m = pdim[mg], n = pdim[!mg];

	int ki, kj,
		mi = si == R_NilValue,
		mj = sj == R_NilValue,
		ni = (mi) ? m : LENGTH(si),
		nj = (mj) ? n : LENGTH(sj),
		*pi = (mi) ? NULL : INTEGER(si),
		*pj = (mj) ? NULL : INTEGER(sj);

	if (anyNA(pi, ni) || anyNA(pj, nj))
		error(_("NA subscripts in %s not supported for '%s' inheriting from %s"),
		      "x[i, j]", "x", "sparseMatrix");
	int stay = class[0] == 'p';

	PROTECT_INDEX pid_ans;
	SEXP ans = obj;
	PROTECT_WITH_INDEX(ans, &pid_ans);

	if (!mi) {
		stay = stay && ni == m && isPerm(pi, m, 1);
		REPROTECT(ans = newObject((stay) ? "pMatrix" : "indMatrix"), pid_ans);

		dim = GET_SLOT(ans, Matrix_DimSym);
		INTEGER(dim)[ mg] = ni;
		INTEGER(dim)[!mg] = n ;

		margin = GET_SLOT(ans, Matrix_marginSym);
		INTEGER(margin)[0] = mg + 1;

		int *tmp = pperm;
		REPROTECT(perm = allocVector(INTSXP, ni), pid_perm);
		pperm = INTEGER(perm);
		SET_SLOT(ans, Matrix_permSym, perm);

		--tmp; /* 1-indexed */
		for (ki = 0; ki < ni; ++ki)
			pperm[ki] = tmp[pi[ki]];
		++tmp; /* 0-indexed */

		m = ni;
	}

	if (!mj) {
		stay = stay && nj == n && isPerm(pj, n, 1);
		REPROTECT(ans = newObject((stay) ? "pMatrix" : (!mg) ? "ngCMatrix" : "ngRMatrix"), pid_ans);

		dim = GET_SLOT(ans, Matrix_DimSym);
		INTEGER(dim)[ mg] = m ;
		INTEGER(dim)[!mg] = nj;

		int *iwork;
		size_t liwork = (size_t) ((stay) ? nj : (int_fast64_t) n + 1 + n + m);
		Matrix_Calloc(iwork, liwork, int);

		if (stay) {

		invertPerm(pj, iwork, nj, 1, 1);

		int *tmp = pperm;
		REPROTECT(perm = allocVector(INTSXP, nj), pid_perm);
		pperm = INTEGER(perm);
		SET_SLOT(ans, Matrix_permSym, perm);

		--iwork; /* 1-indexed */
		for (kj = 0; kj < nj; ++kj)
			pperm[kj] = iwork[tmp[kj]];
		++iwork; /* 0-indexed */

		} else {

		int *pp0 = iwork, *pp_ = pp0 + n + 1, *pi0 = pp_ + n;

		SEXP p1 = PROTECT(allocVector(INTSXP, (R_xlen_t) nj + 1));
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
			error(_("%s too dense for %s; would have more than %s nonzero entries"),
			      "x[i, j]", "[CR]sparseMatrix", "2^31-1");

		SEXP i1 = PROTECT(allocVector(INTSXP, pp1[nj - 1]));
		int *pi1 = INTEGER(i1);
		SET_SLOT(ans, (!mg) ? Matrix_iSym : Matrix_jSym, i1);

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
