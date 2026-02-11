/* C implementation of methods for symmpart */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP sparse_as_kind(SEXP, const char *, char);

SEXP dense_symmpart(SEXP from, const char *class,
                    char op_ul, char op_ct)
{
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

SEXP sparse_symmpart(SEXP from, const char *class,
                     char op_ul, char op_ct)
{
	PROTECT(from = sparse_as_kind(from, class, ','));

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

	char cl[] = ".s.Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[2] = class[2];
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

	int up = (class[2] != 'R') == (op_ul == 'U' || ul == 'U');

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			j, k, kend, nnz0 = pp0[n], nnz1;
		pp0++;

		if (class[1] == 'g') {

			int *iwork = NULL;
			size_t liwork = (size_t) ((int_fast64_t) n + n + 1 + nnz0);
			Matrix_Calloc(iwork, liwork, int);

			int *pp0_ = iwork + n + 1, *pi0_ = iwork + n + n + 1;
			ncsptrans(pp0_ - 1, pi0_, NULL, pp0 - 1, pi0, NULL, n, n, 'T', iwork);
			memcpy(iwork, pp0 - 1, sizeof(int) * (size_t) n);

			SEXP p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0)));
			int *pp1 = INTEGER(p1), k_, kend_;
			*(pp1++) = 0;
			SET_SLOT(to, Matrix_pSym, p1);

			for (j = 0, k = 0, k_ = 0; j < n; ++j) {
				kend  = pp0 [j];
				kend_ = pp0_[j];
				pp1[j] = pp1[j - 1];
				while (k < kend) {
					if (pi0[k] > j)
						k = kend;
					else {
						while (k_ < kend_ && pi0_[k_] < pi0[k]) {
							++pp1[j];
							++k_;
						}
						++pp1[j];
						if (k_ < kend_ && pi0_[k_] == pi0[k])
							++k_;
						++k;
					}
				}
				while (k_ < kend_) {
					if (pi0_[k_] > j)
						k_ = kend_;
					else {
						++pp1[j];
						++k_;
					}
				}
			}
			nnz1 = pp1[n - 1];

			SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
				x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), l;
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				if (op_ct == 'C') \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (up && pi0[k] > j) \
							k = kend; \
						else if (!up && pi0[k] < j) \
							++k; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								*pi1 = pi0_[k_]; \
								c##ASSIGN_CONJ(*px1, px0[l]); \
								c##MULTIPLY(*px1, 0.5); \
								++k_; ++pi1; ++px1; \
							} \
							l = iwork[j]++; \
							*pi1 = pi0[k]; \
							if (pi0[k] == j) { \
								c##ASSIGN_PROJ_REAL(*px1, px0[k]); \
								++k_; \
							} else { \
								c##ASSIGN_IDEN(*px1, px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##INCREMENT_CONJ(*px1, px0[l]); \
									++k_; \
								} \
								c##MULTIPLY(*px1, 0.5); \
							} \
							++k; ++pi1; ++px1; \
						} \
					} \
					while (k_ < kend_) { \
						if (up && pi0_[k_] > j) \
							k_ = kend_; \
						else if (!up && pi0_[k_] < j) \
							++k_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							*pi1 = pi0_[k_]; \
							c##ASSIGN_CONJ(*px1, px0[l]); \
							c##MULTIPLY(*px1, 0.5); \
							++k_; ++pi1; ++px1; \
						} \
					} \
				} \
				else \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (up && pi0[k] > j) \
							k = kend; \
						else if (!up && pi0[k] < j) \
							++k; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								*pi1 = pi0_[k_]; \
								c##ASSIGN_IDEN(*px1, px0[l]); \
								c##MULTIPLY(*px1, 0.5); \
								++k_; ++pi1; ++px1; \
							} \
							l = iwork[j]++; \
							*pi1 = pi0[k]; \
							if (pi0[k] == j) { \
								c##ASSIGN_IDEN(*px1, px0[k]); \
								++k_; \
							} else { \
								c##ASSIGN_IDEN(*px1, px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##INCREMENT_IDEN(*px1, px0[l]); \
									++k_; \
								} \
								c##MULTIPLY(*px1, 0.5); \
							} \
							++k; ++pi1; ++px1; \
						} \
					} \
					while (k_ < kend_) { \
						if (up && pi0_[k_] > j) \
							k_ = kend_; \
						else if (!up && pi0_[k_] < j) \
							++k_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							*pi1 = pi0_[k_]; \
							c##ASSIGN_IDEN(*px1, px0[l]); \
							c##MULTIPLY(*px1, 0.5); \
							++k_; ++pi1; ++px1; \
						} \
					} \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			Matrix_Free(iwork, liwork);
			UNPROTECT(3); /* x1, i1, p1 */

		} else if (class[1] == 's') {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			/* Symmetric part of Hermitian matrix is real part */
			/* Hermitian part of symmetric matrix is real part */
			SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, nnz0));
			zvreal(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else if (nu == 'N') {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz0));
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					while (k < kend) { \
						if (pi0[k] != j) { \
							c##ASSIGN_IDEN(*px1, px0[k]); \
							c##MULTIPLY(*px1, 0.5); \
						} \
						else if (op_ct == 'C') \
							c##ASSIGN_PROJ_REAL(*px1, px0[k]); \
						else \
							c##ASSIGN_IDEN(*px1, px0[k]); \
						++k; ++px1; \
					} \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			UNPROTECT(1); /* x1 */

		} else {

			nnz1 = nnz0 + n;

			SEXP p1 = PROTECT(Rf_allocVector(INTSXP, XLENGTH(p0))),
				i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
				x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz1));
			int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1);
			*(pp1++) = 0;
			SET_SLOT(to, Matrix_pSym, p1);
			SET_SLOT(to,        iSym, i1);
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (j = 0, k = 0; j < n; ++j) { \
					kend = pp0[j]; \
					if (!up) { \
						*pi1 = j; \
						*px1 = c##UNIT; \
						++pi1; ++px1; \
					} \
					while (k < kend) { \
						*pi1 = pi0[k]; \
						c##ASSIGN_IDEN(*px1, px0[k]); \
						c##MULTIPLY(*px1, 0.5); \
						++k; ++pi1; ++px1; \
					} \
					if (up) { \
						*pi1 = j; \
						*px1 = c##UNIT; \
						++pi1; ++px1; \
					} \
					pp1[j] = kend + j + 1; \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			UNPROTECT(3); /* x1, i1, p1 */

		}

		UNPROTECT(3); /* x0, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1;

		if (class[1] == 'g') {

			SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz0)),
				j1 = PROTECT(Rf_allocVector(INTSXP, nnz0)),
				x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz0));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				if (op_ct == 'C') \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						if ((up) ? *pi0 < *pj0 : *pi0 > *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						c##ASSIGN_IDEN(*px1, *px0); \
						} else { \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
						c##ASSIGN_CONJ(*px1, *px0); \
						} \
						c##MULTIPLY(*px1, 0.5); \
					} else { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						c##ASSIGN_PROJ_REAL(*px1, *px0); \
					} \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
				else \
				for (k = 0; k < nnz0; ++k) { \
					if ((up) ? *pi0 <= *pj0 : *pi0 >= *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
					} else { \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
					} \
					c##ASSIGN_IDEN(*px1, *px0); \
					if (*pi0 != *pj0) \
					c##MULTIPLY(*px1, 0.5); \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			UNPROTECT(3); /* x1, j1, i1 */

		} else if (class[1] == 's') {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			/* Symmetric part of Hermitian matrix is real part */
			/* Hermitian part of symmetric matrix is real part */
			SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, nnz0));
			zvreal(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else if (nu == 'N') {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz0));
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						c##ASSIGN_IDEN(*px1, *px0); \
						c##MULTIPLY(*px1, 0.5); \
					} \
					else if (op_ct == 'C') \
						c##ASSIGN_PROJ_REAL(*px1, *px0); \
					else \
						c##ASSIGN_IDEN(*px1, *px0); \
					++pi0; ++pj0; ++px0; ++px1; \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			UNPROTECT(1); /* x1 */

		} else {

			nnz1 = nnz0 + n;

			SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
				j1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
				x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1), j;
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) { \
					*pi1 = *pi0; \
					*pj1 = *pj0; \
					c##ASSIGN_IDEN(*px1, *px0); \
					c##MULTIPLY(*px1, 0.5); \
					++pi0; ++pi1; ++pj0; ++pj1; ++px0; ++px1; \
				} \
				for (j = 0; j < n; ++j) { \
					*pi1 = *pj1 = j; \
					*px1 = c##UNIT; \
					++pi1; ++pj1; ++px1; \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			UNPROTECT(3); /* x1, j1, i1 */

		}

		UNPROTECT(3); /* x0, j0, i1 */

	}

	UNPROTECT(2); /* to, from */
	return to;
}

SEXP R_dense_symmpart(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ul, ct;
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return dense_symmpart(s_from, class, ul, ct);
}

SEXP R_sparse_symmpart(SEXP s_from, SEXP s_uplo, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	char ul, ct;
	VALID_UPLO (s_uplo , ul);
	VALID_TRANS(s_trans, ct);

	return sparse_symmpart(s_from, class, ul, ct);
}
