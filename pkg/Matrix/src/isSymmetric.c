/* C implementation of methods for isSymmetric */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

/* defined in ./isDiagonal.c :*/
int dense_is_diagonal(SEXP, const char *);
int sparse_is_diagonal(SEXP, const char *);

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

	if (class[1] == 'g') {

#define TEMPLATE(c) \
	do { \
		c##TYPE *px = c##PTR(x); \
		if (c##NAME(test2)(px, (size_t) n, '\0', op_ct, '\0')) \
			return 0; \
	} while (0)

	SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

	} else {

	Rcomplex *px = COMPLEX(x), *pu = px, *pl = px;
	int packed = class[2] == 'p';
	if (class[1] == 's') {
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

int sparse_is_symmetric(SEXP obj, const char *class,
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
		if (!sparse_is_diagonal(obj, class))
			return 0;
		nu = DIAG(obj);
		if (nu != 'N' || op_ct != 'C')
			return 1;
	}

	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		return 0;
	if (n == 0 || (n == 1 && op_ct != 'C'))
		return 1;
	if (!exact && class[0] != 'g')
		return NA_LOGICAL; /* do inexact numerical test in R */

	if (class[2] == 'T') {
		/* defined in ./coerce.c : */
		SEXP sparse_as_Csparse(SEXP, const char *);
		if (exact)
		obj = sparse_as_Csparse(obj, class);
		else {
		/* Testing for symmetric nonzero pattern, hence : */
		char cl[] = "n.TMatrix";
		cl[1] = class[1];
		obj = sparse_as_Csparse(obj, cl);
		}
	}
	PROTECT(obj);

	SEXP iSym = (class[2] != 'R') ? Matrix_iSym : Matrix_jSym,
		p_ = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i_ = PROTECT(GET_SLOT(obj,        iSym));
	int *pp = INTEGER(p_), *pi = INTEGER(i_), i, j, k, kend,
		*iwork = NULL;
	if (class[1] == 'g') {
	Matrix_Calloc(iwork, n, int);
	memcpy(iwork, pp, sizeof(int) * (size_t) n);
	}
	pp++;

	int ans = 0;

	if (class[1] == 'g' && !(exact && op_ct == 'C')) {

#define TEMPLATE(c) \
	do { \
		c##IF_NPATTERN( \
		SEXP x = GET_SLOT(obj, Matrix_xSym); \
		c##TYPE *px = c##PTR(x); \
		); \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				i = pi[k]; \
				if (i < j) { \
					if (iwork[i] == pp[i] || pi[iwork[i]] != j) \
						goto done; \
					c##IF_NPATTERN( \
					if (c##NOT_IDEN(px[k], px[iwork[i]])) \
						goto done; \
					); \
					++iwork[i]; \
					++iwork[j]; \
					++k; \
				} else { \
					if (pi[k] == j) \
						++iwork[j]; \
					k = kend; \
				} \
			} \
		} \
	} while (0)

	SWITCH5((exact) ? class[0] : 'n', TEMPLATE);

#undef TEMPLATE

	} else {

	SEXP x = GET_SLOT(obj, Matrix_xSym);
	Rcomplex *px = COMPLEX(x);
	if (class[1] == 'g')
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			while (k < kend) {
				i = pi[k];
				if (i < j) {
					if (iwork[i] == pp[i] || pi[iwork[i]] != j) \
						goto done;
					if (zNOT_CONJ(px[k], px[iwork[i]]))
						goto done;
					++iwork[i];
					++iwork[j];
					++k;
				} else {
					if (i == j) {
						if (zNOT_ZERO_IMAG(px[k]))
							goto done;
						++iwork[j];
					}
					k = kend;
				}
			}
		}
	else if (class[1] == 's')
		/* Testing if Hermitian matrix is symmetric */
		/*      or if symmetric matrix is Hermitian */
		/* <=====> if matrix is real                */
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			while (k < kend) {
				if (zNOT_ZERO_IMAG(px[k]) && (pi[k] != j || op_ct == 'C'))
				    goto done;
				++k;
			}
		}
	else
		/* Testing if non-unit triangular matrix is Hermitian */
		/* <=====> if matrix is real and diagonal             */
		/* NB: diagonal nonzero pattern is already known      */
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && zNOT_ZERO_IMAG(px[k]))
				goto done;
			k = kend;
		}

	}

	if (class[1] == 'g')
	for (j = 0; j < n; ++j)
		if (iwork[j] != pp[j])
			goto done;

	ans = (exact) ? 1 : NA_LOGICAL;

done:
	if (class[1] == 'g')
	Matrix_Free(iwork, n);
	UNPROTECT(3); /* i_, p_, obj */
	return ans;
}

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
	UNPROTECT(1);
	return ans;
}


/* NB: requires symmetric nonzero pattern */
SEXP R_sparse_is_symmetric(SEXP s_obj,
                           SEXP s_trans, SEXP s_exact, SEXP s_checkDN)
{
	const char *class = Matrix_class(s_obj, valid_sparse, 6, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	int exact, checkDN;
	VALID_LOGIC2(s_exact  , exact  );
	VALID_LOGIC2(s_checkDN, checkDN);

	int ans_ = sparse_is_symmetric(s_obj, class, ct, exact, checkDN);
	SEXP ans = Rf_ScalarLogical(ans_);
	return ans;
}
