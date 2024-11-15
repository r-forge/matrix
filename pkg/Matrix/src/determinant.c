/* C implementation of methods for determinant */

#include "cholmod-etc.h"
#include "Mdefines.h"
#include <Rmath.h> /* logspace_add, logspace_sub */

static
SEXP det(double modulus, int logarithm, int sign)
{
	SEXP s_modulus = PROTECT(Rf_ScalarReal((logarithm) ? modulus : exp(modulus))),
		s_logarithm = PROTECT(Rf_ScalarLogical(logarithm)),
		s_sign = PROTECT(Rf_ScalarInteger(sign)),
		ans = PROTECT(Rf_allocVector(VECSXP, 2)),
		class = PROTECT(Rf_allocVector(STRSXP, 1)),
		names = PROTECT(Rf_allocVector(STRSXP, 2));
	static SEXP logarithmSym = NULL;
	if (!logarithmSym)
		logarithmSym = Rf_install("logarithm");
	SET_STRING_ELT(class, 0, Rf_mkChar("det"));
	SET_STRING_ELT(names, 0, Rf_mkChar("modulus"));
	SET_STRING_ELT(names, 1, Rf_mkChar("sign"));
	Rf_setAttrib(ans, R_ClassSymbol, class);
	Rf_setAttrib(ans, R_NamesSymbol, names);
	Rf_setAttrib(s_modulus, logarithmSym, s_logarithm);
	SET_VECTOR_ELT(ans, 0, s_modulus);
	SET_VECTOR_ELT(ans, 1, s_sign);
	UNPROTECT(6);
	return ans;
}

SEXP denseLU_determinant(SEXP s_trf, SEXP s_logarithm)
{

#define DETERMINANT_START(_F_) \
	int *pdim = DIM(_F_), n = pdim[1]; \
	if (pdim[0] != n) \
		Rf_error(_("matrix is not square")); \
	int givelog = Rf_asLogical(s_logarithm); \
	double modulus = 0.0; /* result for n == 0 */

	DETERMINANT_START(s_trf);

	SEXP x = PROTECT(GET_SLOT(s_trf, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	if (n > 0) {
	int j;
	R_xlen_t n1a = (R_xlen_t) n + 1;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			modulus += log(hypot((*px).r, (*px).i));
			px += n1a;
		}
	} else {
		SEXP pivot = GET_SLOT(s_trf, Matrix_permSym);
		int *ppivot = INTEGER(pivot);
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			if (ISNAN(*px) || *px >= 0.0) {
				modulus += log(*px);
				if (*ppivot != j + 1)
					sign = -sign;
			} else {
				modulus += log(-(*px));
				if (*ppivot == j + 1)
					sign = -sign;
			}
			px += n1a;
			ppivot += 1;
		}
	}
	}

	UNPROTECT(1); /* x */
	return det(modulus, givelog, sign);
}

SEXP denseBunchKaufman_determinant(SEXP s_trf, SEXP s_logarithm)
{
	DETERMINANT_START(s_trf);

	SEXP x = PROTECT(GET_SLOT(s_trf, Matrix_xSym));
	int sign = 1;

	char ct = (TYPEOF(x) == CPLXSXP) ? TRANS(s_trf) : 'C';
	if (ct != 'C')
		sign = NA_INTEGER;

	if (n > 0) {
	char ul = UPLO(s_trf);

	SEXP pivot = GET_SLOT(s_trf, Matrix_permSym);
	int *ppivot = INTEGER(pivot);

	int j = 0, packed = XLENGTH(x) != (int_fast64_t) n * n;
	R_xlen_t n1a = (R_xlen_t) n + 1;
	double logab, logcc;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x), a, b, c;
		while (j < n) {
			a = *px;
			if (ppivot[j] > 0) {
				if (ct != 'C')
					modulus += log(hypot(a.r, a.i));
				else if (ISNAN(a.r) || a.r >= 0.0)
					modulus += log(a.r);
				else {
					modulus += log(-a.r);
					sign = -sign;
				}
				px += (!packed) ? n1a : ((ul == 'U') ? j + 2 : n - j);
				j += 1;
			} else {
				if (ul == 'U') {
					px += (!packed) ? n1a : j + 2;
					b = *px;
					c = *(px - 1);
					px += (!packed) ? n1a : j + 3;
				} else {
					c = *(px + 1);
					px += (!packed) ? n1a : n - j;
					b = *px;
					px += (!packed) ? n1a : n - j - 1;
				}
				if (ct != 'C')
					modulus += log(hypot(a.r * b.r - a.i * b.i -
					                     c.r * c.r + c.i * c.i,
					                     a.r * b.i + a.i * b.r -
					                     2.0 * c.r * c.i));
				else {
				logab = log(fabs(a.r)) + log(fabs(b.r));
				logcc = 2.0 * log(hypot(c.r, c.i));
				if ((a.r < 0.0) != (b.r < 0.0)) {
					/* det = ab - cc = -(abs(ab) + cc) < 0 */
					modulus += Rf_logspace_add(logab, logcc);
					sign = -sign;
				} else if (logab < logcc) {
					/* det = ab - cc = -(cc - ab) < 0 */
					modulus += Rf_logspace_sub(logcc, logab);
					sign = -sign;
				} else {
					/* det = ab - cc > 0 */
					modulus += Rf_logspace_sub(logab, logcc);
				}
				}
				j += 2;
			}
		}
	} else {
		double *px = REAL(x), a, b, c;
		while (j < n) {
			a = *px;
			if (ppivot[j] > 0) {
				if (*px >= 0.0)
					modulus += log(a);
				else {
					modulus += log(-a);
					sign = -sign;
				}
				px += (!packed) ? n1a : ((ul == 'U') ? j + 2 : n - j);
				j += 1;
			} else {
				if (ul == 'U') {
					px += (!packed) ? n1a : j + 2;
					b = *px;
					c = *(px - 1);
					px += (!packed) ? n1a : j + 3;
				} else {
					c = *(px + 1);
					px += (!packed) ? n1a : n - j;
					b = *px;
					px += (!packed) ? n1a : n - j - 1;
				}
				logab = log(fabs(a)) + log(fabs(b));
				logcc = 2.0 * log(fabs(c));
				if ((a < 0.0) != (b < 0.0)) {
					/* det = ab - cc = -(abs(ab) + cc) < 0 */
					modulus += Rf_logspace_add(logab, logcc);
					sign = -sign;
				} else if (logab < logcc) {
					/* det = ab - cc = -(cc - ab) < 0 */
					modulus += Rf_logspace_sub(logcc, logab);
					sign = -sign;
				} else {
					/* det = ab - cc > 0 */
					modulus += Rf_logspace_sub(logab, logcc);
				}
				j += 2;
			}
		}
	}
	}

	UNPROTECT(1); /* x */
	return det(modulus, givelog, sign);
}

SEXP denseCholesky_determinant(SEXP s_trf, SEXP s_logarithm)
{
	DETERMINANT_START(s_trf);

	SEXP x = PROTECT(GET_SLOT(s_trf, Matrix_xSym));
	int sign = 1;

	if (n > 0) {
	char ul = UPLO(s_trf);

	int j, packed = XLENGTH(x) != (int_fast64_t) n * n;
	R_xlen_t n1a = (R_xlen_t) n + 1;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0; j < n; ++j) {
			modulus += log((*px).r);
			px += (!packed) ? n1a : ((ul == 'U') ? j + 2 : n - j);
		}
	} else {
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			modulus += log(*px);
			px += (!packed) ? n1a : ((ul == 'U') ? j + 2 : n - j);
		}
	}
	modulus *= 2.0;
	}

	UNPROTECT(1); /* x */
	return det(modulus, givelog, sign);
}

SEXP sparseQR_determinant(SEXP orf, SEXP s_logarithm)
{
	DETERMINANT_START(orf);

	SEXP R = PROTECT(GET_SLOT(orf, Matrix_RSym)),
		x = PROTECT(GET_SLOT(R, Matrix_xSym));
	int sign = 1;

	if (DIM(R)[0] > n)
		Rf_error(_("%s(<%s>) does not support structurally rank deficient case"),
		         "determinant", "sparseQR");

	if (n > 0) {
	SEXP p = PROTECT(GET_SLOT(R, Matrix_pSym)),
		i = PROTECT(GET_SLOT(R, Matrix_iSym));
	int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k, kend;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j) {
				if (ISNAN(px[kend - 1].r) || px[kend - 1].r >= 0.0)
					modulus += log(px[kend - 1].r);
				else {
					modulus += log(-px[kend - 1].r);
					sign = -sign;
				}
			} else {
				UNPROTECT(4); /* i, p, x, R */
				return det(R_NegInf, givelog, 1);
			}
			k = kend;
		}
	} else {
		double *px = REAL(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j) {
				if (ISNAN(px[kend - 1]) || px[kend - 1] >= 0.0)
					modulus += log(px[kend - 1]);
				else {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				}
			} else {
				UNPROTECT(4); /* i, p, x, R */
				return det(R_NegInf, givelog, 1);
			}
			k = kend;
		}
	}
	UNPROTECT(2); /* i, p */

	/* defined in ./perm.c : */
	int signPerm(const int *, int, int);

	p = GET_SLOT(orf, Matrix_pSym);
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
		sign = -sign;
	p = GET_SLOT(orf, Matrix_qSym);
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
		sign = -sign;
	if (n % 2)
		sign = -sign;
	}

	UNPROTECT(2); /* x, R */
	return det(modulus, givelog, sign);
}

SEXP sparseLU_determinant(SEXP s_trf, SEXP s_logarithm)
{
	DETERMINANT_START(s_trf);

	SEXP U = PROTECT(GET_SLOT(s_trf, Matrix_USym)),
		x = PROTECT(GET_SLOT(U, Matrix_xSym));
	int sign = (TYPEOF(x) == CPLXSXP) ? NA_INTEGER : 1;

	if (n > 0) {
	SEXP p = PROTECT(GET_SLOT(U, Matrix_pSym)),
		i = PROTECT(GET_SLOT(U, Matrix_iSym));
	int *pp = INTEGER(p) + 1, *pi = INTEGER(i), j, k, kend;
	if (TYPEOF(x) == CPLXSXP) {
		Rcomplex *px = COMPLEX(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j)
				modulus += log(hypot(px[kend - 1].r, px[kend - 1].i));
			else {
				UNPROTECT(4); /* i, p, x, U */
				return det(R_NegInf, givelog, 1);
			}
			k = kend;
		}
	} else {
		double *px = REAL(x);
		for (j = 0, k = 0; j < n; ++j) {
			kend = pp[j];
			if (k < kend && pi[kend - 1] == j) {
				if (ISNAN(px[kend - 1]) || px[kend - 1] >= 0.0)
					modulus += log(px[kend - 1]);
				else {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				}
			} else {
				UNPROTECT(4); /* i, p, x, U */
				return det(R_NegInf, givelog, 1);
			}
			k = kend;
		}
	}
	UNPROTECT(2); /* i, p */

	if (TYPEOF(x) != CPLXSXP) {
	/* defined in ./perm.c : */
	int signPerm(const int *, int, int);

	p = GET_SLOT(s_trf, Matrix_pSym);
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
		sign = -sign;
	p = GET_SLOT(s_trf, Matrix_qSym);
	if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
		sign = -sign;
	}
	}

	UNPROTECT(2); /* x, U */
	return det(modulus, givelog, sign);
}

SEXP sparseCholesky_determinant(SEXP s_trf, SEXP s_logarithm, SEXP s_root)
{
	DETERMINANT_START(s_trf);

	cholmod_factor *L = M2CHF(s_trf, 1);
	int sign = 1;

	if (n > 0) {
	int j, root = Rf_asLogical(s_root);
	if (L->is_super) {
		int k, nc,
			nsuper = (int) L->nsuper,
			*psuper = (int *) L->super,
			*ppi = (int *) L->pi,
			*ppx = (int *) L->px;
		R_xlen_t nr1a;
		if (L->xtype == CHOLMOD_COMPLEX) {
			Rcomplex *px = (Rcomplex *) L->x, *py;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k + 1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k + 1] - ppi[k]) + 1;
				py = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					modulus += log((*py).r);
					py += nr1a;
				}
			}
		} else {
			double *px = (double *) L->x, *py;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k + 1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k + 1] - ppi[k]) + 1;
				py = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					modulus += log(*py);
					py += nr1a;
				}
			}
		}
		modulus *= 2.0;
	} else {
		int *pp = (int *) L->p;
		if (L->xtype == CHOLMOD_COMPLEX) {
			Rcomplex *px = (Rcomplex *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j)
					modulus += log(px[pp[j]].r);
				modulus *= 2.0;
			} else {
				for (j = 0; j < n; ++j) {
					if (ISNAN(px[pp[j]].r) || px[pp[j]].r >= 0.0)
						modulus += log(px[pp[j]].r);
					else {
						if (root)
							return det(R_NaN, givelog, 1);
						modulus += log(-px[pp[j]].r);
						sign = -sign;
					}
				}
			}
		} else {
			double *px = (double *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j)
					modulus += log(px[pp[j]]);
				modulus *= 2.0;
			} else {
				for (j = 0; j < n; ++j) {
					if (ISNAN(px[pp[j]]) || px[pp[j]] >= 0.0)
						modulus += log(px[pp[j]]);
					else {
						if (root)
							return det(R_NaN, givelog, 1);
						modulus += log(-px[pp[j]]);
						sign = -sign;
					}
				}
			}
		}
	}
	if (root)
		modulus *= 0.5;
	}

	return det(modulus, givelog, sign);
}
