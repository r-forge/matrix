#include <Rmath.h> /* math.h, logspace_add, logspace_sub */
#include "Mdefines.h"
#include "cholmod-etc.h"
#include "determinant.h"

/* defined in ./perm.c : */
int signPerm(const int *, int, int);

static
SEXP mkDet(double modulus, int logarithm, int sign)
{
	SEXP nms = PROTECT(allocVector(STRSXP, 2)),
		cl = PROTECT(mkString("det")),
		det = PROTECT(allocVector(VECSXP, 2)),
		det0 = PROTECT(ScalarReal((logarithm) ? modulus : exp(modulus))),
		det1 = PROTECT(ScalarInteger(sign)),
		det0a = PROTECT(ScalarLogical(logarithm));
	SET_STRING_ELT(nms, 0, mkChar("modulus"));
	SET_STRING_ELT(nms, 1, mkChar("sign"));
	setAttrib(det, R_NamesSymbol, nms);
	setAttrib(det, R_ClassSymbol, cl);
	setAttrib(det0, install("logarithm"), det0a);
	SET_VECTOR_ELT(det, 0, det0);
	SET_VECTOR_ELT(det, 1, det1);
	UNPROTECT(6);
	return det;
}

SEXP denseLU_determinant(SEXP obj, SEXP logarithm)
{

#define DETERMINANT_START \
	SEXP dim = GET_SLOT(obj, Matrix_DimSym); \
	int *pdim = INTEGER(dim), m = pdim[0], n = pdim[1]; \
	if (m != n) \
		error(_("determinant of non-square matrix is undefined")); \
	int givelog = asLogical(logarithm) != 0, sign = 1; \
	double modulus = 0.0; /* result for n == 0 */

	DETERMINANT_START;
	if (n > 0) {
		SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		int j, *ppivot = INTEGER(pivot);
		R_xlen_t n1a = (R_xlen_t) n + 1;
		double *px = REAL(x);

		for (j = 0; j < n; ++j, px += n1a, ++ppivot) {
			if (*px < 0.0) {
				modulus += log(-(*px));
				if (*ppivot == j + 1)
					sign = -sign;
			} else {
				/* incl. 0, NaN cases */
				modulus += log(*px);
				if (*ppivot != j + 1)
					sign = -sign;
			}
		}
		UNPROTECT(2); /* x, pivot */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP BunchKaufman_determinant(SEXP obj, SEXP logarithm, SEXP packed)
{
	DETERMINANT_START;
	if (n > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
		UNPROTECT(1); /* uplo */

		SEXP pivot = PROTECT(GET_SLOT(obj, Matrix_permSym)),
			x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		int j = 0, unpacked = !asLogical(packed),
			*ppivot = INTEGER(pivot);
		R_xlen_t n1a = (R_xlen_t) n + 1;
		double *px = REAL(x), a, b, c, logab, logcc;
		while (j < n) {
			if (ppivot[j] > 0) {
				if (*px < 0.0) {
					modulus += log(-(*px));
					sign = -sign;
				} else {
					/* incl. 0, NaN cases */
					modulus += log(*px);
				}
				px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
				j += 1;
			} else {
				a = *px;
				if (upper) {
					px += (unpacked) ? n1a : j + 2;
					b = *px;
					c = *(px - 1);
					px += (unpacked) ? n1a : j + 3;
				} else {
					c = *(px + 1);
					px += (unpacked) ? n1a : n - j;
					b = *px;
					px += (unpacked) ? n1a : n - j - 1;
				}
				logab = log(fabs(a)) + log(fabs(b));
				logcc = 2.0 * log(fabs(c));
				if ((a < 0.0) != (b < 0.0)) {
					/* det = ab - cc = -(abs(ab) + cc) < 0 */
					modulus += logspace_add(logab, logcc);
					sign = -sign;
				} else if (logab < logcc) {
					/* det = ab - cc = -(cc - ab) < 0 */
					modulus += logspace_sub(logcc, logab);
					sign = -sign;
				} else {
					/* det = ab - cc > 0 */
					modulus += logspace_sub(logab, logcc);
				}
				j += 2;
			}
		}
		UNPROTECT(2); /* x, pivot */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP Cholesky_determinant(SEXP obj, SEXP logarithm, SEXP packed)
{
	DETERMINANT_START;
	if (n > 0) {
		SEXP uplo = PROTECT(GET_SLOT(obj, Matrix_uploSym));
		int upper = *CHAR(STRING_ELT(uplo, 0)) == 'U';
		UNPROTECT(1); /* uplo */

		SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
		int j, unpacked = !asLogical(packed);
		R_xlen_t n1a = (R_xlen_t) n + 1;
		double *px = REAL(x);
		for (j = 0; j < n; ++j) {
			if (*px < 0.0) {
				modulus += log(-(*px));
				sign = -sign;
			} else {
				/* incl. 0, NaN cases */
				modulus += log(*px);
			}
			px += (unpacked) ? n1a : ((upper) ? j + 2 : n - j);
		}
		modulus *= 2.0;
		UNPROTECT(1); /* x */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP sparseLU_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START;
	if (n > 0) {
		SEXP U = PROTECT(GET_SLOT(obj, Matrix_USym)),
			p = PROTECT(GET_SLOT(U, Matrix_pSym)),
			i = PROTECT(GET_SLOT(U, Matrix_iSym)),
			x = PROTECT(GET_SLOT(U, Matrix_xSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
		double *px = REAL(x);

		for (j = 0; j < n; ++j) {
			kend = *(++pp);
			if (kend > k && pi[kend - 1] == j) {
				if (px[kend - 1] < 0.0) {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				} else {
					/* incl. 0, NaN cases */
					modulus += log(px[kend - 1]);
				}
			} else {
				UNPROTECT(4); /* x, i, p, U */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}
		UNPROTECT(4); /* x, i, p, U */

		PROTECT(p = GET_SLOT(obj, Matrix_pSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
		PROTECT(p = GET_SLOT(obj, Matrix_qSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
	}
	return mkDet(modulus, givelog, sign);
}

SEXP sparseQR_determinant(SEXP obj, SEXP logarithm)
{
	DETERMINANT_START;
	if (n > 0) {
		SEXP R = PROTECT(GET_SLOT(obj, Matrix_RSym));
		PROTECT(dim = GET_SLOT(R, Matrix_DimSym));
		if (INTEGER(dim)[0] > n)
			error(_("%s(<%s>) does not support structurally rank deficient case"),
			      "determinant", "sparseQR");
		UNPROTECT(1); /* dim */

		SEXP p = PROTECT(GET_SLOT(R, Matrix_pSym)),
			i = PROTECT(GET_SLOT(R, Matrix_iSym)),
			x = PROTECT(GET_SLOT(R, Matrix_xSym));
		int *pp = INTEGER(p), *pi = INTEGER(i), j, k = 0, kend;
		double *px = REAL(x);

		for (j = 0; j < n; ++j) {
			kend = *(++pp);
			if (kend > k && pi[kend - 1] == j) {
				if (px[kend - 1] < 0.0) {
					modulus += log(-px[kend - 1]);
					sign = -sign;
				} else {
					/* incl. 0, NaN cases */
					modulus += log(px[kend - 1]);
				}
			} else {
				UNPROTECT(4); /* x, i, p, R */
				return mkDet(R_NegInf, givelog, 1);
			}
			k = kend;
		}
		UNPROTECT(4); /* x, i, p, U */

		PROTECT(p = GET_SLOT(obj, Matrix_pSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
		PROTECT(p = GET_SLOT(obj, Matrix_qSym));
		if (signPerm(INTEGER(p), LENGTH(p), 0) < 0)
			sign = -sign;
		UNPROTECT(1); /* p */
		if (n % 2)
			sign = -sign;
	}
	return mkDet(modulus, givelog, sign);
}

SEXP CHMfactor_determinant(SEXP obj, SEXP logarithm, SEXP sqrt)
{
	DETERMINANT_START;
	if (n > 0) {
		int sqrt_ = asLogical(sqrt);
		cholmod_factor *L = mf2cholmod(obj);
		if (L->is_super) {
			int k, j, nc,
				nsuper = (int) L->nsuper,
				*psuper = (int *) L->super,
				*ppi = (int *) L->pi,
				*ppx = (int *) L->px;
			double *px = (double *) L->x, *px_;
			R_xlen_t nr1a;
			for (k = 0; k < nsuper; ++k) {
				nc = psuper[k+1] - psuper[k];
				nr1a = (R_xlen_t) (ppi[k+1] - ppi[k]) + 1;
				px_ = px + ppx[k];
				for (j = 0; j < nc; ++j) {
					modulus += log(*px_);
					px_ += nr1a;
				}
			}
			modulus *= 2.0;
		} else {
			int j, *pp = (int *) L->p;
			double *px = (double *) L->x;
			if (L->is_ll) {
				for (j = 0; j < n; ++j)
					modulus += log(px[pp[j]]);
				modulus *= 2.0;
			} else {
				for (j = 0; j < n; ++j) {
					if (px[pp[j]] < 0.0) {
						if (sqrt_)
							return mkDet(R_NaN, givelog, 1);
						modulus += log(-px[pp[j]]);
						sign = -sign;
					} else {
						/* incl. 0, NaN cases */
						modulus += log(px[pp[j]]);
					}
				}
			}
		}
		if (sqrt_)
			modulus *= 0.5;
	}
	return mkDet(modulus, givelog, sign);
}
