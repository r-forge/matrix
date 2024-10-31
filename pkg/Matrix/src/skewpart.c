/* C implementation of methods for skewpart */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP sparse_as_kind(SEXP, const char *, char);

SEXP dense_skewpart(SEXP from, const char *class, char op_ct)
{
	PROTECT(from = dense_as_kind(from, class, ',', 0));

	if (class[0] != 'z')
		op_ct = '\0';

	int packed = class[2] == 'p';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 's') ? 's' : 'g';
	cl[2] = (class[1] == 's') ? class[2] : 'e';
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error((op_ct == 'C')
		         ? _("attempt to get skew-Hermitian part of non-square matrix")
		         : _("attempt to get skew-symmetric part of non-square matrix"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));

	char ul = '\0', ct = '\0', nu = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U' && class[1] == 's')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && (nu = DIAG(from)) != 'N')
		;

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), x1;

	if (class[1] == 's') {

		/* Skew-symmetric part of Hermitian matrix is imaginary part */
		/* Skew-Hermitian part of symmetric matrix is imaginary part */
		if (op_ct == ct)
		PROTECT(x1 = allocZero(TYPEOF(x0), XLENGTH(x0)));
		else {
		PROTECT(x1 = Rf_allocVector(CPLXSXP, XLENGTH(x0)));
		zvimag(COMPLEX(x1), COMPLEX(x0), (size_t) XLENGTH(x0));
		}

	} else {

		if ((int_fast64_t) n * n > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		PROTECT(x1 = Rf_allocVector(TYPEOF(x0), (R_xlen_t) n * n));

		int i, j;

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = c##PTR(x0), *pu0 = px0, *pl0 = px0; \
			c##TYPE *px1 = c##PTR(x1), *pu1 = px1, *pl1 = px1; \
			if (class[1] == 'g') { \
				if (op_ct == 'C') \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##DECREMENT_CONJ(*pu1, *pl0); \
							c##ASSIGN_CONJ(*pl1, *pu1); \
							c##MULTIPLY(*pu1,  0.5); \
							c##MULTIPLY(*pl1, -0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl0 += n; \
							pl1 += n; \
						} \
						c##ASSIGN_PROJ_IMAG(*pu1, *pu0); \
						pu0 += n - j; \
						pu1 += n - j; \
						pl0 = px0 + j + 1; \
						pl1 = px1 + j + 1; \
					} \
				else \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##DECREMENT_IDEN(*pu1, *pl0); \
							c##ASSIGN_IDEN(*pl1, *pu1); \
							c##MULTIPLY(*pu1,  0.5); \
							c##MULTIPLY(*pl1, -0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl0 += n; \
							pl1 += n; \
						} \
						c##SET_ZERO(*pu1); \
						pu0 += n - j; \
						pu1 += n - j; \
						pl0 = px0 + j + 1; \
						pl1 = px1 + j + 1; \
					} \
			} else { \
				if (ul == 'U') \
					for (j = 0; j < n; ++j) { \
						for (i = 0; i < j; ++i) { \
							c##ASSIGN_IDEN(*pu1, *pu0); \
							c##ASSIGN_IDEN(*pl1, *pu0); \
							c##MULTIPLY(*pu1,  0.5); \
							c##MULTIPLY(*pl1, -0.5); \
							pu0 += 1; \
							pu1 += 1; \
							pl1 += n; \
						} \
						if (nu != 'N' || op_ct != 'C') \
							c##SET_ZERO(*pu1); \
						else \
							c##ASSIGN_PROJ_IMAG(*pu1, *pu0); \
						pu0 += 1; \
						pu1 += 1; \
						pl1 += n; \
						if (!packed) \
						pu0 += n - j - 1; \
						pu1 += n - j - 1; \
						pl1 = px1 + j + 1; \
					} \
				else \
					for (j = 0; j < n; ++j) { \
						if (!packed) \
						pl0 += j; \
						pl1 += j; \
						pu1 = pl1; \
						if (nu != 'N' || op_ct != 'C') \
							c##SET_ZERO(*pl1); \
						else \
							c##ASSIGN_PROJ_IMAG(*pl1, *pl0); \
						pl0 += 1; \
						pl1 += 1; \
						pu1 += n; \
						for (i = j + 1; i < n; ++i) { \
							c##ASSIGN_IDEN(*pl1, *pl0); \
							c##ASSIGN_IDEN(*pu1, *pl0); \
							c##MULTIPLY(*pl1,  0.5); \
							c##MULTIPLY(*pu1, -0.5); \
							pl0 += 1; \
							pl1 += 1; \
							pu1 += n; \
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

SEXP sparse_skewpart(SEXP from, const char *class, char op_ct)
{
	PROTECT(from = sparse_as_kind(from, class, ','));

	if (class[0] != 'z')
		op_ct = '\0';

	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 's') ? 's' : 'g';
	cl[2] = class[2];
	SEXP to = PROTECT(newObject(cl));

	int *pdim = DIM(from), n = pdim[1];
	if (pdim[0] != n)
		Rf_error((op_ct == 'C')
		         ? _("attempt to get skew-Hermitian part of non-square matrix")
		         : _("attempt to get skew-symmetric part of non-square matrix"));
	SET_DIM(to, n, n);
	SET_DIMNAMES(to, -(class[1] != 's'), DIMNAMES(from, 0));

	char ul = '\0', ct = '\0';
	if (class[1] != 'g' && (ul = UPLO(from)) != 'U' && class[1] == 's')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && (ct = TRANS(from)) != 'C')
		SET_TRANS(to);

	if (class[1] == 's' && op_ct == ct) {
		if (class[2] != 'T') {
			SEXP p1 = PROTECT(allocZero(INTSXP, (R_xlen_t) n + 1));
			SET_SLOT(to, Matrix_pSym, p1);
			UNPROTECT(1); /* p1 */
		}
		UNPROTECT(2); /* to, from */
		return to;
	}

	if (class[2] != 'T') {

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0),
			j, k, kend, nnz0 = pp0[n], nnz1;
		pp0++;

		if (class[1] == 's') {

			SET_SLOT(to, Matrix_pSym, p0);
			SET_SLOT(to,        iSym, i0);

			/* Skew-symmetric part of Hermitian matrix is imaginary part */
			/* Skew-Hermitian part of symmetric matrix is imaginary part */
			SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, nnz0));
			zvimag(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else {

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
				pp1[j] = 0;
				while (k < kend) {
					if (pi0[k] >= j)
						k = kend;
					else {
						while (k_ < kend_ && pi0_[k_] < pi0[k]) {
							++pp1[j];
							++pp1[pi0_[k_]];
							++k_;
						}
						++pp1[j];
						++pp1[pi0[k]];
						if (k_ < kend_ && pi0_[k_] == pi0[k])
							++k_;
						++k;
					}
				}
				while (k_ < kend_) {
					if (pi0_[k_] >= j)
						k_ = kend_;
					else {
						++pp1[j];
						++pp1[pi0_[k_]];
						++k_;
					}
				}
			}
			for (j = 0; j < n; ++j)
				pp1[j] += pp1[j - 1];
			nnz1 = pp1[n - 1];
			for (j = n - 1; j >= 0; --j)
				pp1[j] = pp1[j - 1];

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
						if (pi0[k] > j) \
							k = kend; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								c##ASSIGN_CONJ(px1[pp1[       j]], px0[l]); \
								c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
								c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
								c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
								pi1[pp1[       j]++] = pi0_[k_]; \
								pi1[pp1[pi0_[k_]]++] =        j; \
								++k_; \
							} \
							l = iwork[j]++; \
							if (pi0[k] == j) \
								++k_; \
							else { \
								c##ASSIGN_IDEN(px1[pp1[     j]], px0[k]); \
								c##ASSIGN_CONJ(px1[pp1[pi0[k]]], px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##DECREMENT_CONJ(px1[pp1[     j]], px0[l]); \
									c##DECREMENT_IDEN(px1[pp1[pi0[k]]], px0[l]); \
									++k_; \
								} \
								c##MULTIPLY   (px1[pp1[     j]],    0.5); \
								c##MULTIPLY   (px1[pp1[pi0[k]]],   -0.5); \
								pi1[pp1[     j]++] = pi0[k]; \
								pi1[pp1[pi0[k]]++] =      j; \
							} \
							++k; \
						} \
					} \
					while (k_ < kend_) { \
						if (pi0_[k_] >= j) \
							k_ = kend_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							c##ASSIGN_CONJ(px1[pp1[       j]], px0[l]); \
							c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
							c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
							c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
							pi1[pp1[       j]++] = pi0_[k_]; \
							pi1[pp1[pi0_[k_]]++] =        j; \
							++k_; \
						} \
					} \
				} \
				else \
				for (j = 0, k = 0, k_ = 0; j < n; ++j) { \
					kend  = pp0 [j]; \
					kend_ = pp0_[j]; \
					while (k < kend) { \
						if (pi0[k] > j) \
							k = kend; \
						else { \
							while (k_ < kend_ && pi0_[k_] < pi0[k]) { \
								l = iwork[pi0_[k_]]++; \
								c##ASSIGN_IDEN(px1[pp1[       j]], px0[l]); \
								c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
								c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
								c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
								pi1[pp1[       j]++] = pi0_[k_]; \
								pi1[pp1[pi0_[k_]]++] =        j; \
								++k_; \
							} \
							l = iwork[j]++; \
							if (pi0[k] == j) \
								++k_; \
							else { \
								c##ASSIGN_IDEN(px1[pp1[     j]], px0[k]); \
								c##ASSIGN_IDEN(px1[pp1[pi0[k]]], px0[k]); \
								if (k_ < kend_ && pi0_[k_] == pi0[k]) { \
									l = iwork[pi0[k]]++; \
									c##DECREMENT_IDEN(px1[pp1[     j]], px0[l]); \
									c##DECREMENT_IDEN(px1[pp1[pi0[k]]], px0[l]); \
									++k_; \
								} \
								c##MULTIPLY   (px1[pp1[     j]],    0.5); \
								c##MULTIPLY   (px1[pp1[pi0[k]]],   -0.5); \
								pi1[pp1[     j]++] = pi0[k]; \
								pi1[pp1[pi0[k]]++] =      j; \
							} \
							++k; \
						} \
					} \
					while (k_ < kend_) { \
						if (pi0_[k_] >= j) \
							k_ = kend_; \
						else { \
							l = iwork[pi0_[k_]]++; \
							c##ASSIGN_IDEN(px1[pp1[       j]], px0[l]); \
							c##ASSIGN_IDEN(px1[pp1[pi0_[k_]]], px0[l]); \
							c##MULTIPLY   (px1[pp1[       j]],   -0.5); \
							c##MULTIPLY   (px1[pp1[pi0_[k_]]],    0.5); \
							pi1[pp1[       j]++] = pi0_[k_]; \
							pi1[pp1[pi0_[k_]]++] =        j; \
							++k_; \
						} \
					} \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			Matrix_Free(iwork, liwork);
			UNPROTECT(3); /* x1, i1, p1 */

		}

		UNPROTECT(3); /* x0, i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym)),
			x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
		int *pi0 = INTEGER(i0), *pj0 = INTEGER(j0);
		R_xlen_t k, nnz0 = XLENGTH(i0), nnz1;

		if (class[1] == 's') {

			SET_SLOT(to, Matrix_iSym, i0);
			SET_SLOT(to, Matrix_jSym, j0);

			/* Skew-symmetric part of Hermitian matrix is imaginary part */
			/* Skew-Hermitian part of symmetric matrix is imaginary part */
			SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, nnz0));
			zvimag(COMPLEX(x1), COMPLEX(x0), (size_t) nnz0);
			SET_SLOT(to, Matrix_xSym, x1);
			UNPROTECT(1); /* x1 */

		} else {

			nnz1 = nnz0;
			for (k = 0; k < nnz0; ++k)
				if (pi0[k] == pj0[k])
					--nnz1;
			nnz1 *= 2;

			SEXP i1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
				j1 = PROTECT(Rf_allocVector(INTSXP, nnz1)),
				x1 = PROTECT(Rf_allocVector(TYPEOF(x0), nnz1));
			int *pi1 = INTEGER(i1), *pj1 = INTEGER(j1);
			SET_SLOT(to, Matrix_iSym, i1);
			SET_SLOT(to, Matrix_jSym, j1);
			SET_SLOT(to, Matrix_xSym, x1);

#define TEMPLATE(c) \
			do { \
				c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
				for (k = 0; k < nnz0; ++k) { \
					if (*pi0 != *pj0) { \
						*pi1 = *pi0; \
						*pj1 = *pj0; \
						c##ASSIGN_IDEN(*px1, *px0); \
						c##MULTIPLY(*px1,  0.5); \
						++pi1; ++pj1; ++px1; \
						*pi1 = *pj0; \
						*pj1 = *pi0; \
						if (op_ct == 'C') \
						c##ASSIGN_CONJ(*px1, *px0); \
						else \
						c##ASSIGN_IDEN(*px1, *px0); \
						c##MULTIPLY(*px1, -0.5); \
						++pi1; ++pj1; ++px1; \
					} \
					++pi0; ++pj0; ++px0; \
				} \
			} while (0)

			SWITCH2(cl[0], TEMPLATE);

#undef TEMPLATE

			UNPROTECT(3); /* x1, j1, i1 */

		}

		UNPROTECT(3); /* x0, j0, i0 */

	}

	UNPROTECT(2); /* to, from */
	return to;
}

SEXP R_dense_skewpart(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 6, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	return dense_skewpart(s_from, class, ct);
}

SEXP R_sparse_skewpart(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_sparse, 6, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	return sparse_skewpart(s_from, class, ct);
}
