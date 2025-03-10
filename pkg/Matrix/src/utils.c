#include "Mdefines.h"

#if R_VERSION < R_Version(4, 5, 0)
int ANY_ATTRIB(SEXP x)
{
	return ATTRIB(x) != R_NilValue;
}

void CLEAR_ATTRIB(SEXP x)
{
	SET_ATTRIB(x, R_NilValue);
	SET_OBJECT(x, 0);
	UNSET_S4_OBJECT(x);
	return;
}
#endif

char *Matrix_sprintf(const char *format, ...)
{
	char *buf = R_alloc(Matrix_ErrorBufferSize, sizeof(char));
	va_list args;
	va_start(args, format);
	vsnprintf(buf, Matrix_ErrorBufferSize, format, args);
	va_end(args);
	return buf;
}

int equalString(SEXP s1, SEXP s2, R_xlen_t n)
{
	SEXP s1_, s2_;
	R_xlen_t j;
	for (j = 0; j < n; ++j) {
		s1_ = STRING_ELT(s1, j);
		s2_ = STRING_ELT(s2, j);
		if ((s1_ == NA_STRING) != (s2_ == NA_STRING) ||
		    strcmp(CHAR(s1_), CHAR(s2_)) != 0)
			return 0;
	}
	return 1;
}

SEXP duplicateVector(SEXP x)
{
	SEXPTYPE type = TYPEOF(x);
	R_xlen_t length = XLENGTH(x);
	SEXP ans = Rf_allocVector(type, length);
	switch (type) {
	case RAWSXP:
		memcpy(    RAW(ans),     RAW(x), sizeof(   Rbyte) * (size_t) length);
		break;
	case LGLSXP:
		memcpy(LOGICAL(ans), LOGICAL(x), sizeof(     int) * (size_t) length);
		break;
	case INTSXP:
		memcpy(INTEGER(ans), INTEGER(x), sizeof(     int) * (size_t) length);
		break;
	case REALSXP:
		memcpy(   REAL(ans),    REAL(x), sizeof(  double) * (size_t) length);
		break;
	case CPLXSXP:
		memcpy(COMPLEX(ans), COMPLEX(x), sizeof(Rcomplex) * (size_t) length);
		break;
	default:
		break;
	}
	return ans;
}

SEXP allocZero(SEXPTYPE type, R_xlen_t length)
{
	SEXP ans = Rf_allocVector(type, length);
	switch (type) {
	case RAWSXP:
		memset(    RAW(ans), 0, sizeof(   Rbyte) * (size_t) length);
		break;
	case LGLSXP:
		memset(LOGICAL(ans), 0, sizeof(     int) * (size_t) length);
		break;
	case INTSXP:
		memset(INTEGER(ans), 0, sizeof(     int) * (size_t) length);
		break;
	case REALSXP:
		memset(   REAL(ans), 0, sizeof(  double) * (size_t) length);
		break;
	case CPLXSXP:
		memset(COMPLEX(ans), 0, sizeof(Rcomplex) * (size_t) length);
		break;
	default:
		break;
	}
	return ans;
}

SEXP allocUnit(SEXPTYPE type, R_xlen_t length)
{
	SEXP ans = Rf_allocVector(type, length);
	R_xlen_t i;
	switch (type) {
	case RAWSXP:
	{
		Rbyte *pans = RAW(ans);
		for (i = 0; i < length; ++i)
			*(pans++) = 1;
		break;
	}
	case LGLSXP:
	{
		int *pans = LOGICAL(ans);
		for (i = 0; i < length; ++i)
			*(pans++) = 1;
		break;
	}
	case INTSXP:
	{
		int *pans = INTEGER(ans);
		for (i = 0; i < length; ++i)
			*(pans++) = 1;
		break;
	}
	case REALSXP:
	{
		double *pans = REAL(ans);
		for (i = 0; i < length; ++i)
			*(pans++) = 1.0;
		break;
	}
	case CPLXSXP:
	{
		Rcomplex *pans = COMPLEX(ans),
			u = { .r = 1., .i = 0. };
		for (i = 0; i < length; ++i)
			*(pans++) = u;
		break;
	}
	default:
		break;
	}
	return ans;
}

SEXP allocSeqInt(int from, R_xlen_t length)
{
	SEXP ans = Rf_allocVector(INTSXP, length);
	int *pans = INTEGER(ans);
	R_xlen_t i;
	for (i = 0; i < length; ++i)
		*(pans++) = from++;
	return ans;
}

void naToUnit(SEXP x)
{
	R_xlen_t i, length = XLENGTH(x);
	switch (TYPEOF(x)) {
	case LGLSXP:
	{
		int *px = LOGICAL(x);
		for (i = 0; i < length; ++i) {
			if (*px == NA_LOGICAL)
				*px = 1;
			px += 1;
		}
		break;
	}
	case INTSXP:
	{
		int *px = INTEGER(x);
		for (i = 0; i < length; ++i) {
			if (*px == NA_INTEGER)
				*px = 1;
			px += 1;
		}
		break;
	}
	case REALSXP:
	{
		double *px = REAL(x);
		for (i = 0; i < length; ++i) {
			if (ISNAN(*px))
				*px = 1.0;
			px += 1;
		}
		break;
	}
	case CPLXSXP:
	{
		Rcomplex *px = COMPLEX(x),
			u = { .r = 1., .i = 0. };
		for (i = 0; i < length; ++i) {
			if (ISNAN((*px).r) || ISNAN((*px).i))
				*px = u;
			px += 1;
		}
		break;
	}
	default:
		break;
	}
	return;
}
