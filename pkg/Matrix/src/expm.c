/* C implementation of methods for expm */

#include <complex.h> /* _Complex, cexp */
#include "Lapack-etc.h"
#include "Mdefines.h"
#include "idz.h"

#define SHOW(x, header) \
do { \
	Rprintf("%s ...\n", (header)); \
	for (k = 0; k < nn; ++k) \
		Rprintf("    % .8e\n", (x)[k]); \
} while (0)

/* defined in ./Schur.c : */
SEXP dense_schur(SEXP, const char *, int, int);

/* defined in ./coerce.c : */
SEXP dense_as_kind(SEXP, const char *, char, int);
SEXP dense_as_general(SEXP, const char *, int);
SEXP dense_as_packed(SEXP, const char *, char, char, char);

static double thetam[] = { 1.5e-2, 2.5e-1, 9.5e-1, 2.1e+0, 5.4e+0 };
static double padecm[][14] = {
	{ 120.0, 60.0, 12.0, 1.0 },
	{ 30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0 },
	{ 17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0,
	  1512.0, 56.0, 1.0 },
	{ 17643225600.0, 8821612800.0, 2075673600.0, 302702400.0,
	  30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0 },
	{ 64764752532480000.0, 32382376266240000.0,
	  7771770303897600.0, 1187353796428800.0, 129060195264000.0,
	  10559470521600.0, 670442572800.0, 33522128640.0,
	  1323241920.0, 40840800.0, 960960.0, 16380.0, 182.0, 1.0 } };

SEXP dense_expm(SEXP obj, const char *class)
{
	int *pdim = DIM(obj), n = pdim[1];
	if (pdim[0] != n)
		Rf_error(_("matrix is not square"));
	int cs = class[1] == 's' && class[0] == 'z' && TRANS(obj) != 'C';
	char cl[] = "...Matrix";
	cl[0] = (class[0] == 'z') ? 'z' : 'd';
	cl[1] = (class[1] == 'g' || cs) ? 'g' : (class[1] == 's') ? 'p' : 't';
	cl[2] = (class[1] == 'g' || cs) ? 'e' : (class[1] == 's') ? 'o' : 'r';

	SEXP ans = PROTECT(newObject(cl));
	SET_DIM(ans, n, n);
	SET_DIMNAMES(ans, -(class[1] == 's'), DIMNAMES(obj, 0));
	if (cl[1] != 'g' && UPLO(obj) != 'U')
		SET_UPLO(ans);
	if (n > 0) {

	PROTECT_INDEX pid;
	PROTECT_WITH_INDEX(obj, &pid);
	if (class[0] != 'z' && class[0] != 'd') {
		REPROTECT(obj = dense_as_kind(obj, class, ',', 1), pid);
		class = Matrix_class(obj, valid_dense, 6, __func__);
	}

	if (class[1] != 's' || cs) {

	/* Implementation of Algorithm 2.3 from Higham (2005). */

	REPROTECT(obj = dense_as_general(obj, class, 1), pid);

	SEXP x0 = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));
	size_t n1 = (size_t) n, nn = n1 * n1, dk = n1 + 1, k;
	int *pivot = (int *) R_alloc(n1, sizeof(int)),
		i, j, l, s, ilo, ihi, info;
	double *scale = (double *) R_alloc(n1, sizeof(double)),
		norm, *b;

	if (TYPEOF(x0) == CPLXSXP) {

	Rcomplex *a = COMPLEX(x0), *x = COMPLEX(x1), *work, *tmp,
		zero = Matrix_zzero, unit = Matrix_zunit, trace = zero;
	memcpy(x, a, sizeof(Rcomplex) * nn);
	k = 0;
	for (j = 0; j < n; ++j) {
		trace.r += x[k].r;
		trace.i += x[k].i;
		k += dk;
	}
	trace.r /= (double) n;
	trace.i /= (double) n;
	k = 0;
	for (j = 0; j < n; ++j) {
		x[k].r -= trace.r;
		x[k].i -= trace.i;
		k += dk;
	}
	F77_CALL(zgebal)("B", &n, x, &n, &ilo, &ihi, scale, &info FCONE);
	ilo -= 1; ihi -= 1;
	ERROR_LAPACK_1(zgebal, info);
	norm = F77_CALL(zlange)("O", &n, &n, x, &n, (double *) 0 FCONE);
	if (!R_FINITE(norm))
		Rf_error("matrix one-norm is not finite");
	s = (norm <= thetam[4]) ? 0 : (int) ceil(log2(norm / thetam[4]));
	for (l = 0; l < 4; ++l)
		if (norm <= thetam[l])
			break;
	b = &padecm[l][0];
	work = (Rcomplex *) R_alloc(nn * 2, sizeof(Rcomplex));
	Rcomplex *u = work, *v = u + nn;
	memset(u, 0, sizeof(Rcomplex) * nn * 2);
	k = 0;
	for (j = 0; j < n; ++j) {
		u[k].r = b[1];
		v[k].r = b[0];
		k += dk;
	}
	if (l < 4) {
		b += 2;
		work = (Rcomplex *) R_alloc((l > 0) ? nn * 3 : nn, sizeof(Rcomplex));
		Rcomplex *a2 = work, *a2m = work + nn, *a2m2 = a2m + nn;
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, x, &n,
		                x, &n, &zero, a2, &n FCONE FCONE);
		for (k = 0; k < nn; ++k) {
			u[k].r += b[1] * a2[k].r;
			u[k].i += b[1] * a2[k].i;
			v[k].r += b[0] * a2[k].r;
			v[k].i += b[0] * a2[k].i;
		}
		if (l > 0) {
			b += 2;
			memcpy(a2m, a2, sizeof(double) * nn);
			for (i = 0; i < l; ++i) {
				F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, a2, &n,
				                a2m, &n, &zero, a2m2, &n FCONE FCONE);
				for (k = 0; k < nn; ++k) {
					u[k].r += b[1] * a2m2[k].r;
					u[k].i += b[1] * a2m2[k].i;
					v[k].r += b[0] * a2m2[k].r;
					v[k].i += b[0] * a2m2[k].i;
				}
				b += 2;
				tmp = a2m2; a2m2 = a2m; a2m = tmp;
			}
		}
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, x, &n,
		                u, &n, &zero, a2, &n FCONE FCONE);
		memcpy(u, a2, sizeof(Rcomplex) * nn);
	} else {
		work = (Rcomplex *) R_alloc(nn * 5, sizeof(Rcomplex));
		Rcomplex *a2 = work, *a4 = a2 + nn, *a6 = a4 + nn,
			*upart = a6 + nn, *vpart = upart + nn;
		if (s > 0) {
			double e = ldexp(1.0, -s);
			for (k = 0; k < nn; ++k) {
				x[k].r *= e;
				x[k].i *= e;
			}
		}
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, x , &n,
		                x , &n, &zero, a2, &n FCONE FCONE);
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, a2, &n,
		                a2, &n, &zero, a4, &n FCONE FCONE);
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, a2, &n,
		                a4, &n, &zero, a6, &n FCONE FCONE);
		for (k = 0; k < nn; ++k) {
			u[k].r += b[3] * a2[k].r + b[5] * a4[k].r + b[7] * a6[k].r;
			u[k].i += b[3] * a2[k].i + b[5] * a4[k].i + b[7] * a6[k].i;
			v[k].r += b[2] * a2[k].r + b[4] * a4[k].r + b[6] * a6[k].r;
			v[k].i += b[2] * a2[k].i + b[4] * a4[k].i + b[6] * a6[k].i;
			upart[k].r = b[9] * a2[k].r + b[11] * a4[k].r + b[13] * a6[k].r;
			upart[k].i = b[9] * a2[k].i + b[11] * a4[k].i + b[13] * a6[k].i;
			vpart[k].r = b[8] * a2[k].r + b[10] * a4[k].r + b[12] * a6[k].r;
			vpart[k].i = b[8] * a2[k].i + b[10] * a4[k].i + b[12] * a6[k].i;
		}
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, a6, &n,
		                upart, &n, &unit, u, &n FCONE FCONE);
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, a6, &n,
		                vpart, &n, &unit, v, &n FCONE FCONE);
		F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, x, &n,
		                u, &n, &zero, a2, &n FCONE FCONE);
		memcpy(u, a2, sizeof(Rcomplex) * nn);
	}
	for (k = 0; k < nn; ++k) {
		x[k].r =  u[k].r + v[k].r;
		x[k].i =  u[k].i + v[k].i;
		u[k].r = -u[k].r + v[k].r;
		u[k].i = -u[k].i + v[k].i;
	}
	F77_CALL(zgetrf)(&n, &n, u, &n, pivot, &info);
	ERROR_LAPACK_2(zgetrf, info, 2, U);
	F77_CALL(zgetrs)("N", &n, &n, u, &n, pivot, x, &n, &info FCONE);
	ERROR_LAPACK_1(zgetrs, info);
	if (s > 0) {
		u = x;
		for (i = 0; i < s; ++i) {
			F77_CALL(zgemm)("N", "N", &n, &n, &n, &unit, u, &n,
			                u, &n, &zero, v, &n FCONE FCONE);
			tmp = v; v = u; u = tmp;
		}
		if (s % 2)
			memcpy(x, u, sizeof(Rcomplex) * nn);
	}
	tmp = x;
	for (j = 0; j < n; ++j) {
		for (i = ilo; i <= ihi; ++i) {
			tmp[i].r *= scale[i];
			tmp[i].i *= scale[i];
		}
		tmp += n;
	}
	tmp = x + n1 * (size_t) ilo;
	for (j = ilo; j <= ihi; ++j) {
		for (i = 0; i < n; ++i) {
			tmp[i].r /= scale[j];
			tmp[i].i /= scale[j];
		}
		tmp += n;
	}
	for (j = ilo - 1; j >= 0; ++j) {
		i = (int) scale[j] - 1;
		if (i != j) {
			zswap2(n1, x + n1 * (size_t) i, 1, x + n1 * (size_t) j, 1);
			zswap2(n1, x + i, n1, x + j, n1);
		}
	}
	for (j = ihi + 1; j < n; ++j) {
		i = (int) scale[j] - 1;
		if (i != j) {
			zswap2(n1, x + n1 * (size_t) i, 1, x + n1 * (size_t) j, 1);
			zswap2(n1, x + i, n1, x + j, n1);
		}
	}
#define asC(x) ((double _Complex *) &(x))[0]
#define asR(x) ((       Rcomplex *) &(x))[0]
	Rcomplex x___;
	double _Complex z = cexp(asC(trace));
	trace = asR(z);
	for (k = 0; k < nn; ++k) {
		x___ = x[k];
		x[k].r = x___.r * trace.r - x___.i * trace.i;
		x[k].i = x___.r * trace.i + x___.i * trace.r;
	}

	} else {

	double *a = REAL(x0), *x = REAL(x1), *work, *tmp,
		zero = 0.0, unit = 1.0, trace = zero;
	memcpy(x, a, sizeof(double) * nn);
	k = 0;
	for (j = 0; j < n; ++j) {
		trace += x[k];
		k += dk;
	}
	trace /= (double) n;
	k = 0;
	for (j = 0; j < n; ++j) {
		x[k] -= trace;
		k += dk;
	}
	F77_CALL(dgebal)("B", &n, x, &n, &ilo, &ihi, scale, &info FCONE);
	ERROR_LAPACK_1(dgebal, info);
	ilo -= 1; ihi -= 1;
	norm = F77_CALL(dlange)("O", &n, &n, x, &n, (double *) 0 FCONE);
	if (!R_FINITE(norm))
		Rf_error("matrix one-norm is not finite");
	s = (norm <= thetam[4]) ? 0 : (int) ceil(log2(norm / thetam[4]));
	for (l = 0; l < 4; ++l)
		if (norm <= thetam[l])
			break;
	b = &padecm[l][0];
	work = (double *) R_alloc(nn * 2, sizeof(double));
	double *u = work, *v = u + nn;
	memset(u, 0, sizeof(double) * nn * 2);
	k = 0;
	for (j = 0; j < n; ++j) {
		u[k] = b[1];
		v[k] = b[0];
		k += dk;
	}
	if (l < 4) {
		b += 2;
		work = (double *) R_alloc((l > 0) ? nn * 3 : nn, sizeof(double));
		double *a2 = work, *a2m = work + nn, *a2m2 = a2m + nn;
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, x, &n,
		                x, &n, &zero, a2, &n FCONE FCONE);
		for (k = 0; k < nn; ++k) {
			u[k] += b[1] * a2[k];
			v[k] += b[0] * a2[k];
		}
		if (l > 0) {
			b += 2;
			memcpy(a2m, a2, sizeof(double) * nn);
			for (i = 0; i < l; ++i) {
				F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, a2, &n,
				                a2m, &n, &zero, a2m2, &n FCONE FCONE);
				for (k = 0; k < nn; ++k) {
					u[k] += b[1] * a2m2[k];
					v[k] += b[0] * a2m2[k];
				}
				b += 2;
				tmp = a2m2; a2m2 = a2m; a2m = tmp;
			}
		}
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, x, &n,
		                u, &n, &zero, a2, &n FCONE FCONE);
		memcpy(u, a2, sizeof(double) * nn);
	} else {
		work = (double *) R_alloc(nn * 5, sizeof(double));
		double *a2 = work, *a4 = a2 + nn, *a6 = a4 + nn,
			*upart = a6 + nn, *vpart = upart + nn;
		if (s > 0) {
			double e = ldexp(1.0, -s);
			for (k = 0; k < nn; ++k)
				x[k] *= e;
		}
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, x , &n,
		                x , &n, &zero, a2, &n FCONE FCONE);
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, a2, &n,
		                a2, &n, &zero, a4, &n FCONE FCONE);
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, a2, &n,
		                a4, &n, &zero, a6, &n FCONE FCONE);
		for (k = 0; k < nn; ++k) {
			u[k] += b[3] * a2[k] + b[5] * a4[k] + b[7] * a6[k];
			v[k] += b[2] * a2[k] + b[4] * a4[k] + b[6] * a6[k];
			upart[k] = b[9] * a2[k] + b[11] * a4[k] + b[13] * a6[k];
			vpart[k] = b[8] * a2[k] + b[10] * a4[k] + b[12] * a6[k];
		}
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, a6, &n,
		                upart, &n, &unit, u, &n FCONE FCONE);
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, a6, &n,
		                vpart, &n, &unit, v, &n FCONE FCONE);
		F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, x, &n,
		                u, &n, &zero, a2, &n FCONE FCONE);
		memcpy(u, a2, sizeof(double) * nn);
	}
	for (k = 0; k < nn; ++k) {
		x[k] =  u[k] + v[k];
		u[k] = -u[k] + v[k];
	}
	F77_CALL(dgetrf)(&n, &n, u, &n, pivot, &info);
	ERROR_LAPACK_2(dgetrf, info, 2, U);
	F77_CALL(dgetrs)("N", &n, &n, u, &n, pivot, x, &n, &info FCONE);
	ERROR_LAPACK_1(dgetrs, info);
	if (l >= 4) {
		u = x;
		for (i = 0; i < s; ++i) {
			F77_CALL(dgemm)("N", "N", &n, &n, &n, &unit, u, &n,
			                u, &n, &zero, v, &n FCONE FCONE);
			tmp = v; v = u; u = tmp;
		}
		if (s % 2)
			memcpy(x, u, sizeof(double) * nn);
	}
	tmp = x;
	for (j = 0; j < n; ++j) {
		for (i = ilo; i <= ihi; ++i)
			tmp[i] *= scale[i];
		tmp += n;
	}
	tmp = x + n1 * (size_t) ilo;
	for (j = ilo; j <= ihi; ++j) {
		for (i = 0; i < n; ++i)
			tmp[i] /= scale[j];
		tmp += n;
	}
	for (j = ilo - 1; j >= 0; --j) {
		i = (int) scale[j] - 1;
		if (i != j) {
			dswap2(n1, x + n1 * (size_t) i, 1, x + n1 * (size_t) j, 1);
			dswap2(n1, x + i, n1, x + j, n1);
		}
	}
	for (j = ihi + 1; j < n; ++j) {
		i = (int) scale[j] - 1;
		if (i != j) {
			dswap2(n1, x + n1 * (size_t) i, 1, x + n1 * (size_t) j, 1);
			dswap2(n1, x + i, n1, x + j, n1);
		}
	}
	trace = exp(trace);
	for (k = 0; k < nn; ++k)
		x[k] *= trace;

	}

	SET_SLOT(ans, Matrix_xSym, x1);
	UNPROTECT(2); /* x1, x0 */

	} else {

	/* Use the Schur factorization as the argument is Hermitian. */

	REPROTECT(obj = dense_schur(obj, class, 2, 1), pid);
	SEXP vectors = PROTECT(GET_SLOT(obj, Matrix_vectorsSym)),
		values = PROTECT(GET_SLOT(obj, Matrix_valuesSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(vectors), XLENGTH(vectors)));
	double *w = REAL(values);
	size_t n1 = (size_t) n, nn = n1 * n1;
	int i, j;
	for (j = 0; j < n; ++j)
		w[j] = exp(w[j]);
	if (TYPEOF(vectors) == CPLXSXP) {
		Rcomplex *qw = (Rcomplex *) R_alloc(nn, sizeof(Rcomplex)),
			*q = COMPLEX(vectors), *x = COMPLEX(x1), *tmp = qw,
			zero = Matrix_zzero, unit = Matrix_zunit;
		for (j = 0; j < n; ++j) {
			for (i = 0; i < n; ++i) {
				qw[i].r = q[i].r * w[j];
				qw[i].i = q[i].i * w[j];
			}
			qw += n;
			q  += n;
		}
		qw = tmp; q = COMPLEX(vectors);
		F77_CALL(zgemm)("N", "C", &n, &n, &n, &unit, qw, &n,
		                q, &n, &zero, x, &n FCONE FCONE);
	} else {
		double *qw = (double *) R_alloc(nn, sizeof(double)),
			*q = REAL(vectors), *x = REAL(x1), *tmp = qw,
			zero = 0.0, unit = 1.0;
		for (j = 0; j < n; ++j) {
			for (i = 0; i < n; ++i)
				qw[i] = q[i] * w[j];
			qw += n;
			q  += n;
		}
		qw = tmp; q = REAL(vectors);
		F77_CALL(dgemm)("N", "C", &n, &n, &n, &unit, qw, &n,
		                q, &n, &zero, x, &n FCONE FCONE);
	}

	SET_SLOT(ans, Matrix_xSym, x1);
	UNPROTECT(3); /* x1, values, vectors */

	}

	UNPROTECT(1); /* obj */

	}

	if (cl[1] != 'g' && class[2] == 'p') {
		if (class[1] == 's')
		ans = dense_as_packed(ans, cl, '\0',  'C', '\0');
		else
		ans = dense_as_packed(ans, cl, '\0', '\0',  'N');
	}
	UNPROTECT(1); /* ans */
	return ans;
}

SEXP R_dense_expm(SEXP s_obj)
{
	const char *class = Matrix_class(s_obj, valid_dense, 6, __func__);
	return dense_expm(s_obj, class);
}
