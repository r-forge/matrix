#ifndef MATRIX_M5_H
#define MATRIX_M5_H

#include <R_ext/Error.h>

#define nTYPE int
#define lTYPE int
#define iTYPE int
#define dTYPE double
#define zTYPE Rcomplex

#define nTYPESXP LGLSXP
#define lTYPESXP LGLSXP
#define iTYPESXP INTSXP
#define dTYPESXP REALSXP
#define zTYPESXP CPLXSXP

#define nPTR LOGICAL
#define lPTR LOGICAL
#define iPTR INTEGER
#define dPTR REAL
#define zPTR COMPLEX

#define nZERO 0
#define lZERO 0
#define iZERO 0
#define dZERO 0.0
#define zZERO Matrix_zzero

#define nUNIT 1
#define lUNIT 1
#define iUNIT 1
#define dUNIT 1.0
#define zUNIT Matrix_zunit

#define nNA 1
#define lNA NA_LOGICAL
#define iNA NA_INTEGER
#define dNA NA_REAL
#define zNA Matrix_zna

#define nNAME(s) i ## s
#define lNAME(s) i ## s
#define iNAME(s) i ## s
#define dNAME(s) d ## s
#define zNAME(s) z ## s

#define nIF_NPATTERN(...)
#define lIF_NPATTERN(...) __VA_ARGS__
#define iIF_NPATTERN(...) __VA_ARGS__
#define dIF_NPATTERN(...) __VA_ARGS__
#define zIF_NPATTERN(...) __VA_ARGS__

#define nIFELSE_NPATTERN(x, y) y
#define lIFELSE_NPATTERN(x, y) x
#define iIFELSE_NPATTERN(x, y) x
#define dIFELSE_NPATTERN(x, y) x
#define zIFELSE_NPATTERN(x, y) x

#define nNOT_ZERO(x) ((x) != 0)
#define lNOT_ZERO(x) ((x) != 0)
#define iNOT_ZERO(x) ((x) != 0)
#define dNOT_ZERO(x) ((x) != 0.0)
#define zNOT_ZERO(x) ((x).r != 0.0 || (x).i != 0.0)

#define nNOT_ZERO_REAL(x) ((x) != 0)
#define lNOT_ZERO_REAL(x) ((x) != 0)
#define iNOT_ZERO_REAL(x) ((x) != 0)
#define dNOT_ZERO_REAL(x) ((x) != 0.0)
#define zNOT_ZERO_REAL(x) ((x).r != 0.0)

#define nNOT_ZERO_IMAG(x) (0)
#define lNOT_ZERO_IMAG(x) (0)
#define iNOT_ZERO_IMAG(x) (0)
#define dNOT_ZERO_IMAG(x) (0)
#define zNOT_ZERO_IMAG(x) ((x).i != 0.0)

#define nNOT_ZERO_TOL(x, tol) ((x) != 0)
#define lNOT_ZERO_TOL(x, tol) ((x) != 0)
#define iNOT_ZERO_TOL(x, tol) ((x) != 0)
#define dNOT_ZERO_TOL(x, tol) (ISNAN(x) || fabs(x) > tol)
#define zNOT_ZERO_TOL(x, tol) (ISNAN((x).r) || ISNAN((x).i) || hypot((x).r, (x).i) > tol)

#define nNOT_UNIT(x) ((x) == 0)
#define lNOT_UNIT(x) ((x) == 0 || (x) == NA_LOGICAL)
#define iNOT_UNIT(x) ((x) != 1)
#define dNOT_UNIT(x) ((x) != 1.0)
#define zNOT_UNIT(x) ((x).r != 1.0 || (x).i != 0.0)

#define nNOT_NA(x) (1)
#define lNOT_NA(x) ((x) != NA_LOGICAL)
#define iNOT_NA(x) ((x) != NA_INTEGER)
#define dNOT_NA(x) (!ISNAN(x))
#define zNOT_NA(x) (!ISNAN((x).r) && !ISNAN((x).i))

#define nNOT_IDEN(x, y) \
	(((x) != 0) != ((y) != 0))
#define lNOT_IDEN(x, y) \
	((x) != (y))
#define iNOT_IDEN(x, y) \
	((x) != (y))
#define dNOT_IDEN(x, y) \
	((ISNAN(x)) ? !ISNAN(y) : ISNAN(y) || (x) != (y))
#define zNOT_IDEN(x, y) \
	(((ISNAN((x).r)) ? !ISNAN((y).r) : ISNAN((y).r) || (x).r != (y).r) || \
	 ((ISNAN((x).i)) ? !ISNAN((y).i) : ISNAN((y).i) || (x).i != (y).i))

#define nNOT_CONJ(x, y) \
	((x != 0) != (y != 0))
#define lNOT_CONJ(x, y) \
	(x != y)
#define iNOT_CONJ(x, y) \
	(x != y)
#define dNOT_CONJ(x, y) \
	((ISNAN(x)) ? !ISNAN(y) : ISNAN(y) || x != y)
#define zNOT_CONJ(x, y) \
	(((ISNAN((x).r)) ? !ISNAN((y).r) : ISNAN((y).r) || (x).r !=  (y).r) || \
	 ((ISNAN((x).i)) ? !ISNAN((y).i) : ISNAN((y).i) || (x).i != -(y).i))

#define nSET_ZERO(x) ((x) = 0)
#define lSET_ZERO(x) ((x) = 0)
#define iSET_ZERO(x) ((x) = 0)
#define dSET_ZERO(x) ((x) = 0.0)
#define zSET_ZERO(x) ((x) = Matrix_zzero)

#define nSET_UNIT(x) ((x) = 1)
#define lSET_UNIT(x) ((x) = 1)
#define iSET_UNIT(x) ((x) = 1)
#define dSET_UNIT(x) ((x) = 1.0)
#define zSET_UNIT(x) ((x) = Matrix_zunit)

#define nSET_NA(x) ((x) = 1)
#define lSET_NA(x) ((x) = NA_LOGICAL)
#define iSET_NA(x) ((x) = NA_INTEGER)
#define dSET_NA(x) ((x) = NA_REAL)
#define zSET_NA(x) ((x) = Matrix_zna)

#define nSET_IDEN(x)
#define lSET_IDEN(x)
#define iSET_IDEN(x)
#define dSET_IDEN(x)
#define zSET_IDEN(x)

#define nSET_CONJ(x)
#define lSET_CONJ(x)
#define iSET_CONJ(x)
#define dSET_CONJ(x)
#define zSET_CONJ(x) ((x).i = -(x).i)

#define nSET_PROJ_REAL(x)
#define lSET_PROJ_REAL(x)
#define iSET_PROJ_REAL(x)
#define dSET_PROJ_REAL(x)
#define zSET_PROJ_REAL(x) ((x).i = 0.0)

#define nSET_PROJ_IMAG(x) ((x) = 0)
#define lSET_PROJ_IMAG(x) ((x) = 0)
#define iSET_PROJ_IMAG(x) ((x) = 0)
#define dSET_PROJ_IMAG(x) ((x) = 0.0)
#define zSET_PROJ_IMAG(x) ((x).r = 0.0)

#define nASSIGN_IDEN(x, y) ((x) = (y))
#define lASSIGN_IDEN(x, y) ((x) = (y))
#define iASSIGN_IDEN(x, y) ((x) = (y))
#define dASSIGN_IDEN(x, y) ((x) = (y))
#define zASSIGN_IDEN(x, y) ((x) = (y))

#define nASSIGN_CONJ(x, y) ((x) = (y))
#define lASSIGN_CONJ(x, y) ((x) = (y))
#define iASSIGN_CONJ(x, y) ((x) = (y))
#define dASSIGN_CONJ(x, y) ((x) = (y))
#define zASSIGN_CONJ(x, y) ((x).r = (y).r, (x).i = -(y).i)

#define nASSIGN_PROJ_REAL(x, y) ((x) = (y))
#define lASSIGN_PROJ_REAL(x, y) ((x) = (y))
#define iASSIGN_PROJ_REAL(x, y) ((x) = (y))
#define dASSIGN_PROJ_REAL(x, y) ((x) = (y))
#define zASSIGN_PROJ_REAL(x, y) ((x).r = (y).r, (x).i = 0.0)

#define nASSIGN_PROJ_IMAG(x, y) ((x) = 0)
#define lASSIGN_PROJ_IMAG(x, y) ((x) = 0)
#define iASSIGN_PROJ_IMAG(x, y) ((x) = 0)
#define dASSIGN_PROJ_IMAG(x, y) ((x) = 0.0)
#define zASSIGN_PROJ_IMAG(x, y) ((x).r = 0.0, (x).i = (y).i)

#define nINCREMENT_IDEN(x, y) \
	do { \
		(x) = 1; \
	} while (0)
#define lINCREMENT_IDEN(x, y) \
	do { \
		if ((y) == NA_LOGICAL) { \
			if ((x) == 0) \
				(x) = NA_LOGICAL; \
		} \
		else if ((y) != 0) \
			(x) = 1; \
	} while (0)
#define iINCREMENT_IDEN(x, y) \
	do { \
		if ((x) != NA_INTEGER) { \
			if ((y) == NA_INTEGER) \
				(x) = NA_INTEGER; \
			else if ((y) != 0) { \
				if (((y) > 0) \
				    ? ((x) >  INT_MAX - (y)) \
				    : ((x) <= INT_MIN - (y))) { \
					Rf_warning(_("NAs produced by integer overflow")); \
					(x) = NA_INTEGER; \
				} \
				else \
					(x) += (y); \
			} \
		} \
	} while (0)
#define dINCREMENT_IDEN(x, y) \
	do { \
		(x) += (y); \
	} while (0)
#define zINCREMENT_IDEN(x, y) \
	do { \
		(x).r += (y).r; \
		(x).i += (y).i; \
	} while (0)

#define nINCREMENT_CONJ(x, y) nINCREMENT_IDEN(x, y)
#define lINCREMENT_CONJ(x, y) lINCREMENT_IDEN(x, y)
#define iINCREMENT_CONJ(x, y) iINCREMENT_IDEN(x, y)
#define dINCREMENT_CONJ(x, y) dINCREMENT_IDEN(x, y)
#define zINCREMENT_CONJ(x, y) do { (x).r += (y).r; (x).i -= (y).i; } while (0)

#define nDECREMENT_IDEN(x, y)
#define lDECREMENT_IDEN(x, y) \
	do { \
		if ((y) == NA_LOGICAL) { \
			if ((x) == 0) \
				(x) = NA_LOGICAL; \
		} \
		else if ((y) == 0) \
			(x) = 1; \
	} while (0)
#define iDECREMENT_IDEN(x, y) \
	do { \
		if ((x) != NA_INTEGER) { \
			if ((y) == NA_INTEGER) \
				(x) = NA_INTEGER; \
			else if ((y) != 0) { \
				if (((y) < 0) \
				    ? ((x) >  INT_MAX + (y)) \
				    : ((x) <= INT_MIN + (y))) { \
					Rf_warning(_("NAs produced by integer overflow")); \
					(x) = NA_INTEGER; \
				} \
				else \
					(x) -= (y); \
			} \
		} \
	} while (0)
#define dDECREMENT_IDEN(x, y) \
	do { \
		(x) -= (y); \
	} while (0)
#define zDECREMENT_IDEN(x, y) \
	do { \
		(x).r -= (y).r; \
		(x).i -= (y).i; \
	} while (0)

#define nDECREMENT_CONJ(x, y) nDECREMENT_IDEN(x, y)
#define lDECREMENT_CONJ(x, y) lDECREMENT_IDEN(x, y)
#define iDECREMENT_CONJ(x, y) iDECREMENT_IDEN(x, y)
#define dDECREMENT_CONJ(x, y) dDECREMENT_IDEN(x, y)
#define zDECREMENT_CONJ(x, y) do { (x).r -= (y).r; (x).i += (y).i; } while (0)

#define nMULTIPLY(x, a) /* unused hence no-op for now */
#define lMULTIPLY(x, a) /* ditto */
#define iMULTIPLY(x, a) /* ditto */
#define dMULTIPLY(x, a) \
	do { (x)   *= a;             } while (0)
#define zMULTIPLY(x, a) \
	do { (x).r *= a; (x).i *= a; } while (0)

#define nDIVIDE(x, a) /* unused hence no-op for now */
#define lDIVIDE(x, a) /* ditto */
#define iDIVIDE(x, a) /* ditto */
#define dDIVIDE(x, a) \
	do { (x)   /= a;             } while (0)
#define zDIVIDE(x, a) \
	do { (x).r /= a; (x).i /= a; } while (0)

#define SWITCH2(c, template) \
do { \
	switch ((c)) { \
	case 'n': \
	case 'l': \
	case 'i': \
	case 'd': template(d); break; \
	case 'z': template(z); break; \
	default: break; \
	} \
} while (0)

#define SWITCH3(c, template) \
do { \
	switch ((c)) { \
	case 'n': \
	case 'l': \
	case 'i': template(i); break; \
	case 'd': template(d); break; \
	case 'z': template(z); break; \
	default: break; \
	} \
} while (0)

#define SWITCH4(c, template) \
do { \
	switch ((c)) { \
	case 'n': \
	case 'l': template(l); break; \
	case 'i': template(i); break; \
	case 'd': template(d); break; \
	case 'z': template(z); break; \
	default: break; \
	} \
} while (0)

#define SWITCH5(c, template) \
do { \
	switch ((c)) { \
	case 'n': template(n); break; \
	case 'l': template(l); break; \
	case 'i': template(i); break; \
	case 'd': template(d); break; \
	case 'z': template(z); break; \
	default: break; \
	} \
} while (0)

#endif /* MATRIX_M5_H */
