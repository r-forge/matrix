/* C implementation of methods for expand */

#include "Mdefines.h"
#include "M5.h"

SEXP denseBunchKaufman_expand(SEXP s_trf)
{
	SEXP x = PROTECT(GET_SLOT(s_trf, Matrix_xSym));

	SEXP P_ = PROTECT(newObject("pMatrix"));
	char cl[] = "..CMatrix";
	cl[0] = (TYPEOF(x) == CPLXSXP) ? 'z' : 'd';
	cl[1] = 't';
	SEXP T_ = PROTECT(newObject(cl));
	cl[1] = 's';
	SEXP D_ = PROTECT(newObject(cl));

	int n = DIM(s_trf)[1], packed = XLENGTH(x) != (int_fast64_t) n * n;
	SET_DIM(P_, n, n);
	SET_DIM(T_, n, n);
	SET_DIM(D_, n, n);

	char ul = UPLO(s_trf),
		ct = (TYPEOF(x) == CPLXSXP) ? TRANS(s_trf) : 'C';
	if (ul != 'U') {
		SET_UPLO(T_);
		SET_UPLO(D_);
	}
	if (ct != 'C')
		SET_TRANS(D_);
	SET_DIAG(T_);

	int i, j, s;
	R_xlen_t n1a = (R_xlen_t) n + 1;

	SEXP pivot = PROTECT(GET_SLOT(s_trf, Matrix_permSym)),
		D_p = PROTECT(Rf_allocVector(INTSXP, n1a));
	int *ppivot = INTEGER(pivot), *D_pp = INTEGER(D_p),
		b = n, dp = (ul == 'U') ? 1 : 2;
	D_pp[0] = 0;
	j = 0;
	while (j < n) {
		if (ppivot[j] > 0) {
			D_pp[j+1] = D_pp[j] + 1;
			j += 1;
		} else {
			D_pp[j+1] = D_pp[j] + dp;
			D_pp[j+2] = D_pp[j] + 3;
			j += 2;
			b -= 1;
		}
	}
	SET_SLOT(D_, Matrix_pSym, D_p);
	UNPROTECT(1); /* D_p */

	SEXP P, P_perm, T, T_p, T_i, T_x,
		D_i = PROTECT(Rf_allocVector(INTSXP, D_pp[n])),
		D_x = PROTECT(Rf_allocVector(TYPEOF(x), D_pp[n]));
	int *P_pperm, *T_pp, *T_pi, *D_pi = INTEGER(D_i);

	R_xlen_t len = (R_xlen_t) 2 * b + 1, k = (ul == 'U') ? len - 2 : 0;
	SEXP ans = PROTECT(Rf_allocVector(VECSXP, len));

#define TEMPLATE(c) \
	do { \
		c##TYPE *T_px, *D_px = c##PTR(D_x), *px = c##PTR(x); \
		 \
		j = 0; \
		while (b--) { \
			s = (ppivot[j] > 0) ? 1 : 2; \
			dp = (ul == 'U') ? j : n - j - s; \
			 \
			PROTECT(P = Rf_duplicate(P_)); \
			PROTECT(P_perm = Rf_allocVector(INTSXP, n)); \
			PROTECT(T = Rf_duplicate(T_)); \
			PROTECT(T_p = Rf_allocVector(INTSXP, n1a)); \
			PROTECT(T_i = Rf_allocVector(INTSXP, (R_xlen_t) s * dp)); \
			PROTECT(T_x = Rf_allocVector(TYPEOF(x), (R_xlen_t) s * dp)); \
			 \
			P_pperm = INTEGER(P_perm); \
			T_pp = INTEGER(T_p); \
			T_pi = INTEGER(T_i); \
			T_px = c##PTR(T_x); \
			T_pp[0] = 0; \
			 \
			for (i = 0; i < j; ++i) { \
				T_pp[i+1] = 0; \
				P_pperm[i] = i + 1; \
			} \
			for (i = j; i < j+s; ++i) { \
				T_pp[i+1] = T_pp[i] + dp; \
				P_pperm[i] = i + 1; \
			} \
			for (i = j+s; i < n; ++i) { \
				T_pp[i+1] = T_pp[i]; \
				P_pperm[i] = i + 1; \
			} \
			 \
			if (s == 1) { \
				P_pperm[j] = ppivot[j]; \
				P_pperm[ppivot[j]-1] = j + 1; \
			} else if (ul == 'U') { \
				P_pperm[j] = -ppivot[j]; \
				P_pperm[-ppivot[j]-1] = j + 1; \
			} else { \
				P_pperm[j+1] = -ppivot[j]; \
				P_pperm[-ppivot[j]-1] = j + 2; \
			} \
			 \
			if (ul == 'U') { \
				for (i = 0; i < j; ++i) { \
					*(T_pi++) = i; \
					*(T_px++) = *(px++); \
				} \
				*(D_pi++) = j; \
				*(D_px++) = *(px++); \
				++j; \
				if (!packed) \
					px += n - j; \
				if (s == 2) { \
					for (i = 0; i < j-1; ++i) { \
						*(T_pi++) = i; \
						*(T_px++) = *(px++); \
					} \
					*(D_pi++) = j - 1; \
					*(D_pi++) = j; \
					*(D_px++) = *(px++); \
					*(D_px++) = *(px++); \
					++j; \
					if (!packed) \
						px += n - j; \
				} \
			} else { \
				if (s == 2) { \
					*(D_pi++) = j; \
					*(D_pi++) = j + 1; \
					*(D_px++) = *(px++); \
					*(D_px++) = *(px++); \
					for (i = j+2; i < n; ++i) { \
						*(T_pi++) = i; \
						*(T_px++) = *(px++); \
					} \
					++j; \
					if (!packed) \
						px += j; \
				} \
				*(D_pi++) = j; \
				*(D_px++) = *(px++); \
				for (i = j+1; i < n; ++i) { \
					*(T_pi++) = i; \
					*(T_px++) = *(px++); \
				} \
				++j; \
				if (!packed) \
					px += j; \
			} \
			 \
			SET_SLOT(P, Matrix_permSym, P_perm); \
			SET_SLOT(T, Matrix_pSym, T_p); \
			SET_SLOT(T, Matrix_iSym, T_i); \
			SET_SLOT(T, Matrix_xSym, T_x); \
			 \
			if (ul == 'U') { \
				SET_VECTOR_ELT(ans, k-1, P); \
				SET_VECTOR_ELT(ans, k  , T); \
				k -= 2; \
			} else { \
				SET_VECTOR_ELT(ans, k  , P); \
				SET_VECTOR_ELT(ans, k+1, T); \
				k += 2; \
			} \
			UNPROTECT(6); /* T_x, T_i, T_p, T, P_perm, P */ \
		} \
	} while (0)

	if (TYPEOF(x) == CPLXSXP)
		TEMPLATE(z);
	else
		TEMPLATE(d);

	SET_SLOT(D_, Matrix_iSym, D_i);
	SET_SLOT(D_, Matrix_xSym, D_x);
	SET_VECTOR_ELT(ans, len-1, D_);

	UNPROTECT(8); /* ans, D_x, D_i, pivot, D_, T_, P_, x */
	return ans;
}
