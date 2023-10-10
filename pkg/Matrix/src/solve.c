#include "Lapack-etc.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "idz.h"
#include "solve.h"

static
void solveDN(SEXP rdn, SEXP adn, SEXP bdn)
{
	SEXP s;
	if (!isNull(s = VECTOR_ELT(adn, 1)))
		SET_VECTOR_ELT(rdn, 0, s);
	if (!isNull(s = VECTOR_ELT(bdn, 1)))
		SET_VECTOR_ELT(rdn, 1, s);
	PROTECT(adn = getAttrib(adn, R_NamesSymbol));
	PROTECT(bdn = getAttrib(bdn, R_NamesSymbol));
	if(!isNull(adn) || !isNull(bdn)) {
		PROTECT(s = allocVector(STRSXP, 2));
		if (!isNull(adn))
			SET_STRING_ELT(s, 0, STRING_ELT(adn, 1));
		if (!isNull(bdn))
			SET_STRING_ELT(s, 1, STRING_ELT(bdn, 1));
		setAttrib(rdn, R_NamesSymbol, s);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return;
}

SEXP denseLU_solve(SEXP a, SEXP b)
{

#define SOLVE_START \
	SEXP adim = GET_SLOT(a, Matrix_DimSym); \
	int *padim = INTEGER(adim), m = padim[0], n = padim[1]; \
	if (m != n) \
		error(_("'%s' is not square"), "a"); \
	if (!isNull(b)) { \
	SEXP bdim = GET_SLOT(b, Matrix_DimSym); \
	int *pbdim = INTEGER(bdim); \
	if (pbdim[0] != m) \
		error(_("dimensions of '%s' and '%s' are inconsistent"), \
		      "a", "b"); \
	n = pbdim[1]; \
	}

#define SOLVE_FINISH \
	SEXP rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym)), \
		adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)); \
	if (isNull(b)) \
		revDN(rdimnames, adimnames); \
	else { \
		SEXP bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)); \
		solveDN(rdimnames, adimnames, bdimnames); \
		UNPROTECT(1); /* bdimnames */ \
	} \
	UNPROTECT(2); /* adimnames, rdimnames */

	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;

	if (m > 0) {
		SEXP apivot = PROTECT(GET_SLOT(a, Matrix_permSym)), rx;
		int info;
		if (isNull(b)) {
			rx = duplicate(ax);
			PROTECT(rx);
			int lwork = -1;
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			Rcomplex work0, *work = &work0;
			F77_CALL(zgetri)(&m, COMPLEX(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_1(zgetri, info);
			lwork = (int) work0.r;
			work = (Rcomplex *) R_alloc((size_t) lwork, sizeof(Rcomplex));
			F77_CALL(zgetri)(&m, COMPLEX(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_2(zgetri, info, 2, U);
			} else {
#endif
			double   work0, *work = &work0;
			F77_CALL(dgetri)(&m,    REAL(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_1(dgetri, info);
			lwork = (int) work0;
			work = (double   *) R_alloc((size_t) lwork, sizeof(double  ));
			F77_CALL(dgetri)(&m,    REAL(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_2(dgetri, info, 2, U);
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			F77_CALL(zgetrs)("N", &m, &n, COMPLEX(ax), &m, INTEGER(apivot),
			                 COMPLEX(rx), &m, &info FCONE);
			ERROR_LAPACK_1(zgetrs, info);
			} else {
#endif
			F77_CALL(dgetrs)("N", &m, &n,    REAL(ax), &m, INTEGER(apivot),
			                    REAL(rx), &m, &info FCONE);
			ERROR_LAPACK_1(dgetrs, info);
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(2); /* rx, apivot */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP BunchKaufman_solve(SEXP a, SEXP b)
{
	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));
	int unpacked = (Matrix_int_fast64_t) m * m <= R_XLEN_T_MAX &&
		XLENGTH(ax) == (R_xlen_t) m * m;

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	if (!isNull(b)) {
		rcl[1] = 'g';
		rcl[2] = 'e';
	} else {
		rcl[1] = 's';
		rcl[2] = (unpacked) ? 'y' : 'p';
	}
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = CHAR(STRING_ELT(auplo, 0))[0];
	if (isNull(b) && aul != 'U') {
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	if (m > 0) {
		SEXP apivot = PROTECT(GET_SLOT(a, Matrix_permSym)), rx;
		int info;
		if (isNull(b)) {
			rx = duplicate(ax);
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			Rcomplex *work = (Rcomplex *) R_alloc((size_t) m, sizeof(Rcomplex));
			if (unpacked) {
				F77_CALL(zsytri)(&aul, &m, COMPLEX(rx), &m, INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(zsytri, info, 2, D);
			} else {
				F77_CALL(zsptri)(&aul, &m, COMPLEX(rx),     INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(zsptri, info, 2, D);
			}
			} else {
#endif
			double   *work = (double   *) R_alloc((size_t) m, sizeof(double  ));
			if (unpacked) {
				F77_CALL(dsytri)(&aul, &m, REAL(rx), &m, INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(dsytri, info, 2, D);
			} else {
				F77_CALL(dsptri)(&aul, &m, REAL(rx),     INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(dsptri, info, 2, D);
			}
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			if (unpacked) {
				F77_CALL(zsytrs)(&aul, &m, &n, COMPLEX(ax), &m, INTEGER(apivot),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zsytrs, info);
			} else {
				F77_CALL(zsptrs)(&aul, &m, &n, COMPLEX(ax),     INTEGER(apivot),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zsptrs, info);
			}
			} else {
#endif
			if (unpacked) {
				F77_CALL(dsytrs)(&aul, &m, &n,    REAL(ax), &m, INTEGER(apivot),
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dsytrs, info);
			} else {
				F77_CALL(dsptrs)(&aul, &m, &n,    REAL(ax),    INTEGER(apivot),
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dsptrs, info);
			}
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(2); /* rx, apivot */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP Cholesky_solve(SEXP a, SEXP b)
{
	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));
	int unpacked = (Matrix_int_fast64_t) m * m <= R_XLEN_T_MAX &&
		XLENGTH(ax) == (R_xlen_t) m * m;

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	if (!isNull(b)) {
		rcl[1] = 'g';
		rcl[2] = 'e';
	} else {
		rcl[1] = 'p';
		rcl[2] = (unpacked) ? 'o' : 'p';
	}
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = CHAR(STRING_ELT(auplo, 0))[0];
	if (isNull(b) && aul != 'U') {
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	if (m > 0) {
		SEXP rx;
		int info, nprotect = 1;

		SEXP aperm = NULL;
		if (HAS_SLOT(a, Matrix_permSym)) {
			SEXP tmp = GET_SLOT(a, Matrix_permSym);
			if (LENGTH(tmp) > 0) {
				PROTECT(aperm = tmp);
				++nprotect;
			}
		} /* else 'a' is a dtrMatrix or dtpMatrix, as in chol2inv */

		if (isNull(b)) {
			rx = duplicate(ax);
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			if (unpacked) {
				F77_CALL(zpotri)(&aul, &m, COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_2(zpotri, info, 2, L);
				if (aperm)
					zsymperm2(COMPLEX(rx), n, aul, INTEGER(aperm), 1, 1);
			} else {
				F77_CALL(zpptri)(&aul, &m, COMPLEX(rx),     &info FCONE);
				ERROR_LAPACK_2(zpptri, info, 2, L);
				if (aperm) {
					/* FIXME: zsymperm2 supporting packed matrices */
					double *work;
					size_t lwork = (size_t) n * n;
					Matrix_Calloc(work, lwork, Rcomplex);
					zunpack1 (work, COMPLEX(rx), n, aul, 'N');
					zsymperm2(work, n, aul, INTEGER(aperm), 1, 1);
					zpack2   (COMPLEX(rx), work, n, aul, 'N');
					Matrix_Free(work, lwork);
				}
			}
			} else {
#endif
			if (unpacked) {
				F77_CALL(dpotri)(&aul, &m,    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_2(dpotri, info, 2, L);
				if (aperm)
					dsymperm2(   REAL(rx), n, aul, INTEGER(aperm), 1, 1);
			} else {
				F77_CALL(dpptri)(&aul, &m,    REAL(rx),     &info FCONE);
				ERROR_LAPACK_2(dpptri, info, 2, L);
				if (aperm) {
					/* FIXME: dsymperm2 supporting packed matrices */
					double *work;
					size_t lwork = (size_t) n * n;
					Matrix_Calloc(work, lwork, double);
					dunpack1 (work,    REAL(rx), n, aul, 'N');
					dsymperm2(work, n, aul, INTEGER(aperm), 1, 1);
					dpack2   (   REAL(rx), work, n, aul, 'N');
					Matrix_Free(work, lwork);
				}
			}
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			if (aperm)
				zrowperm2(COMPLEX(rx), m, n, INTEGER(aperm), 1, 0);
			if (unpacked) {
				F77_CALL(zpotrs)(&aul, &m, &n, COMPLEX(ax), &m,
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zpotrs, info);
			} else {
				F77_CALL(zpptrs)(&aul, &m, &n, COMPLEX(ax),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zpptrs, info);
			}
			if (aperm)
				zrowperm2(COMPLEX(rx), m, n, INTEGER(aperm), 1, 1);
			} else {
#endif
			if (aperm)
				drowperm2(   REAL(rx), m, n, INTEGER(aperm), 1, 0);
			if (unpacked) {
				F77_CALL(dpotrs)(&aul, &m, &n,    REAL(ax), &m,
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dpotrs, info);
			} else {
				F77_CALL(dpptrs)(&aul, &m, &n,    REAL(ax),
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dpptrs, info);
			}
			if (aperm)
				drowperm2(   REAL(rx), m, n, INTEGER(aperm), 1, 1);
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(nprotect); /* rx, aperm */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP dtrMatrix_solve(SEXP a, SEXP b)
{
	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));
	int unpacked = (Matrix_int_fast64_t) m * m <= R_XLEN_T_MAX &&
		XLENGTH(ax) == (R_xlen_t) m * m;

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	if (!isNull(b)) {
		rcl[1] = 'g';
		rcl[2] = 'e';
	} else {
		rcl[1] = 't';
		rcl[2] = (unpacked) ? 'r' : 'p';
	}
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = m;
	prdim[1] = n;

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = CHAR(STRING_ELT(auplo, 0))[0];
	if (isNull(b) && aul != 'U') {
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	SEXP adiag = GET_SLOT(a, Matrix_diagSym);
	char adi = CHAR(STRING_ELT(adiag, 0))[0];
	if (isNull(b) && adi != 'N') {
		PROTECT(adiag);
		SET_SLOT(r, Matrix_diagSym, adiag);
		UNPROTECT(1); /* adiag */
	}

	if (m > 0) {
		SEXP rx;
		int info;
		if (isNull(b)) {
			rx = duplicate(ax);
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			if (unpacked) {
				F77_CALL(ztrtri)(&aul, &adi, &m, COMPLEX(rx), &m,
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(ztrtri, info, 2, A);
			} else {
				F77_CALL(ztptri)(&aul, &adi, &m, COMPLEX(rx),
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(ztptri, info, 2, A);
			}
			} else {
#endif
			if (unpacked) {
				F77_CALL(dtrtri)(&aul, &adi, &m,    REAL(rx), &m,
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(dtrtri, info, 2, A);
			} else {
				F77_CALL(dtptri)(&aul, &adi, &m,    REAL(rx),
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(dtptri, info, 2, A);
			}
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			rx = duplicate(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
#ifdef MATRIX_ENABLE_ZMATRIX
			if (TYPEOF(ax) == CPLXSXP) {
			if (unpacked) {
				F77_CALL(ztrtrs)(&aul, "N", &adi, &m, &n, COMPLEX(ax), &m,
				                 COMPLEX(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(ztrtrs, info);
			} else {
				F77_CALL(ztptrs)(&aul, "N", &adi, &m, &n, COMPLEX(ax),
				                 COMPLEX(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(ztptrs, info);
			}
			} else {
#endif
			if (unpacked) {
				F77_CALL(dtrtrs)(&aul, "N", &adi, &m, &n,    REAL(ax), &m,
				                    REAL(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(dtrtrs, info);
			} else {
				// https://bugs.r-project.org/show_bug.cgi?id=18534
				F77_CALL(dtptrs)(&aul, "N", &adi, &m, &n,    REAL(ax),
				                    REAL(rx), &m, &info
# ifdef usePR18534fix
				                 FCONE FCONE FCONE);
# else
				                 FCONE FCONE);
# endif
				ERROR_LAPACK_1(dtptrs, info);
			}
#ifdef MATRIX_ENABLE_ZMATRIX
			}
#endif
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP sparseLU_solve(SEXP a, SEXP b, SEXP sparse)
{

#define ERROR_SOLVE_OOM(_ACL_, _BCL_) \
	error(_("%s(<%s>, <%s>) failed: out of memory"), "solve", _ACL_, _BCL_)

	SOLVE_START;

	SEXP r,
		aL = PROTECT(GET_SLOT(a, Matrix_LSym)),
		aU = PROTECT(GET_SLOT(a, Matrix_USym)),
		ap = PROTECT(GET_SLOT(a, Matrix_pSym)),
		aq = PROTECT(GET_SLOT(a, Matrix_qSym));
	int j,
		*pap = INTEGER(ap),
		*paq = (LENGTH(aq)) ? INTEGER(aq) : NULL;
	double *work = (double *) R_alloc((size_t) m, sizeof(double));
	Matrix_cs *L = dgC2cs(aL, 1), *U = dgC2cs(aU, 1);
	if (L->xtype != U->xtype)
		error(_("'%s' and '%s' slots have different '%s'"),
		      "L", "U", "xtype");
	MCS_SET_XTYPE(L->xtype);
	if (!asLogical(sparse)) {
		PROTECT(r = newObject("dgeMatrix"));
		SEXP rdim = GET_SLOT(r, Matrix_DimSym);
		int *prdim = INTEGER(rdim);
		prdim[0] = m;
		prdim[1] = n;
		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP rx = PROTECT(allocVector(REALSXP, mn));
		double *prx = REAL(rx);
		if (isNull(b)) {
			Matrix_memset(prx, 0, mn, sizeof(double));
			for (j = 0; j < n; ++j) {
				prx[j] = 1.0;
				Matrix_cs_pvec(pap, prx, work, m);
				Matrix_cs_lsolve(L, work);
				Matrix_cs_usolve(U, work);
				Matrix_cs_ipvec(paq, work, prx, m);
				prx += m;
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			double *pbx = REAL(bx);
			for (j = 0; j < n; ++j) {
				Matrix_cs_pvec(pap, pbx, work, m);
				Matrix_cs_lsolve(L, work);
				Matrix_cs_usolve(U, work);
				Matrix_cs_ipvec(paq, work, prx, m);
				prx += m;
				pbx += m;
			}
			UNPROTECT(1); /* bx */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	} else {
		Matrix_cs *B, *X;
		int *papinv = Matrix_cs_pinv(pap, m);
		if (!papinv)
			ERROR_SOLVE_OOM("sparseLU", ".gCMatrix");
		if (isNull(b)) {
			B = Matrix_cs_spalloc(m, n, n, 1, 0);
			if (!B)
				ERROR_SOLVE_OOM("sparseLU", ".gCMatrix");
			double *B__x = (double *) B->x;
			for (j = 0; j < n; ++j) {
				B->p[j] = j;
				B->i[j] = j;
				B__x[j] = 1.0;
			}
			B->p[n] = n;
			X = Matrix_cs_permute(B, papinv, NULL, 1);
			B = Matrix_cs_spfree(B);
		} else
			X = Matrix_cs_permute(dgC2cs(b, 1), papinv, NULL, 1);
		papinv = Matrix_cs_free(papinv);
		if (!X)
			ERROR_SOLVE_OOM("sparseLU", ".gCMatrix");
		B = X;

		int i, k, top, nz, nzmax,
			*iwork = (int *) R_alloc((size_t) 2 * m, sizeof(int));

#define DO_TRIANGULAR_SOLVE(_A_, _ALO_, _BFR_, _ACL_, _BCL_) \
		do { \
			X = Matrix_cs_spalloc(m, n, B->nzmax, 1, 0); \
			double *X__x = (double *) X->x; \
			if (!X) { \
				if (_BFR_) \
					B = Matrix_cs_spfree(B); \
				ERROR_SOLVE_OOM(_ACL_, _BCL_); \
			} \
			X->p[0] = nz = 0; \
			nzmax = X->nzmax; \
			for (j = 0, k = 0; j < n; ++j) { \
				top = Matrix_cs_spsolve(_A_, B, j, iwork, work, NULL, _ALO_); \
				if (m - top > INT_MAX - nz) { \
					if (_BFR_) \
						B = Matrix_cs_spfree(B); \
					X = Matrix_cs_spfree(X); \
					error(_("attempt to construct sparse matrix with more than %s nonzero elements"), \
					      "2^31-1"); \
				} \
				nz += m - top; \
				if (nz > nzmax) { \
					nzmax = (nz <= INT_MAX / 2) ? 2 * nz : INT_MAX; \
					if (!Matrix_cs_sprealloc(X, nzmax)) { \
						if (_BFR_) \
							B = Matrix_cs_spfree(B); \
						X = Matrix_cs_spfree(X); \
						ERROR_SOLVE_OOM(_ACL_, _BCL_); \
					} \
					X__x = (double *) X->x; \
				} \
				X->p[j + 1] = nz; \
				if (_ALO_) { \
					for (i = top; i < m; ++i) { \
						X->i[k] =      iwork[i]; \
						X__x[k] = work[iwork[i]]; \
						++k; \
					} \
				} else { \
					for (i = m - 1; i >= top; --i) { \
						X->i[k] =      iwork[i]; \
						X__x[k] = work[iwork[i]]; \
						++k; \
					} \
				} \
			} \
			if (_BFR_) \
				B = Matrix_cs_spfree(B); \
			B = X; \
		} while (0)

		DO_TRIANGULAR_SOLVE(L, 1, 1, "sparseLU", ".gCMatrix");
		DO_TRIANGULAR_SOLVE(U, 0, 1, "sparseLU", ".gCMatrix");

		if (paq) {
			X = Matrix_cs_permute(B, paq, NULL, 1);
			B = Matrix_cs_spfree(B);
			if (!X)
				ERROR_SOLVE_OOM("sparseLU", ".gCMatrix");
			B = X;
		}

		/* Drop zeros from B and sort it : */
		Matrix_cs_dropzeros(B);
		X = Matrix_cs_transpose(B, 1);
		B = Matrix_cs_spfree(B);
		if (!X)
			ERROR_SOLVE_OOM("sparseLU", ".gCMatrix");
		B = Matrix_cs_transpose(X, 1);
		X = Matrix_cs_spfree(X);
		if (!B)
			ERROR_SOLVE_OOM("sparseLU", ".gCMatrix");

		PROTECT(r = cs2dgC(B, 1, 'g'));
		B = Matrix_cs_spfree(B);
	}

	SOLVE_FINISH;

	UNPROTECT(5); /* r, aq, ap, aU, aL */
	return r;
}

static
int strmatch(const char *x, const char **valid)
{
	int i = 0;
	while (valid[i][0] != '\0') {
		if (strcmp(x, valid[i]) == 0)
			return i;
		++i;
	}
	return -1;
}

SEXP CHMfactor_solve(SEXP a, SEXP b, SEXP sparse, SEXP system)
{
	static const char *valid[] = {
		"A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt", "" };
	int ivalid = -1;
	if (TYPEOF(system) != STRSXP || LENGTH(system) < 1 ||
	    (system = STRING_ELT(system, 0)) == NA_STRING ||
	    (ivalid = strmatch(CHAR(system), valid)) < 0)
		error(_("invalid '%s' to '%s'"), "system", __func__);

	SOLVE_START;

	SEXP r;
	int j;
	cholmod_factor *L = M2CF(a, 1);
	if (!asLogical(sparse)) {
		cholmod_dense *B = NULL, *X = NULL;
		if (isNull(b)) {
			B = cholmod_allocate_dense(m, n, m, L->xtype, &c);
			if (!B)
				ERROR_SOLVE_OOM("CHMfactor", ".geMatrix");

#define EYE(_CTYPE_, _ONE_) \
			do { \
				_CTYPE_ *B__x = (_CTYPE_ *) B->x; \
				Matrix_memset(B__x, 0, (R_xlen_t) m * n, sizeof(_CTYPE_)); \
				for (j = 0; j < n; ++j) { \
					*(B__x++) = _ONE_; \
					B__x += m; \
				} \
			} while (0)

#ifdef MATRIX_ENABLE_ZMATRIX
			if (L->xtype == CHOLMOD_COMPLEX)
			EYE(Rcomplex, Matrix_zone);
			else
#endif
			EYE(double, 1.0);

#undef EYE

			X = cholmod_solve(ivalid, L, B, &c);
			cholmod_free_dense(&B, &c);
			if (!X)
				ERROR_SOLVE_OOM("CHMfactor", ".geMatrix");
			PROTECT(r = CD2M(X, 0,
				(ivalid < 2) ? 'p' : ((ivalid < 7) ? 't' : 'g')));
		} else {
			B = M2CD(b, 0);
			X = cholmod_solve(ivalid, L, B, &c);
			if (!X)
				ERROR_SOLVE_OOM("CHMfactor", ".geMatrix");
			PROTECT(r = CD2M(X, 0, 'g'));
		}
		cholmod_free_dense(&X, &c);
	} else {
		cholmod_sparse *B = NULL, *X = NULL;
		if (isNull(b)) {
			B = cholmod_speye(m, n, L->xtype, &c);
			if (!B)
				ERROR_SOLVE_OOM("CHMfactor", ".gCMatrix");
			X = cholmod_spsolve(ivalid, L, B, &c);
			cholmod_free_sparse(&B, &c);
			if (X && ivalid < 7) {
				X->stype = (ivalid == 2 || ivalid == 4) ? -1 : 1;
				cholmod_sort(X, &c);
			}
			if (!X)
				ERROR_SOLVE_OOM("CHMfactor", ".gCMatrix");
			PROTECT(r = CS2M(X, 1,
				(ivalid < 2) ? 's' : ((ivalid < 7) ? 't' : 'g')));
		} else {
			B = M2CS(b, 1);
			X = cholmod_spsolve(ivalid, L, B, &c);
			if (!X)
				ERROR_SOLVE_OOM("CHMfactor", ".gCMatrix");
			PROTECT(r = CS2M(X, 1, 'g'));
		}
		cholmod_free_sparse(&X, &c);
	}
	if (isNull(b) && (ivalid == 2 || ivalid == 4)) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(r, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}

	SOLVE_FINISH;

	UNPROTECT(1); /* r */
	return r;
}

SEXP dtCMatrix_solve(SEXP a, SEXP b, SEXP sparse)
{
	SOLVE_START;

	SEXP r, auplo = PROTECT(GET_SLOT(a, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(auplo, 0));
	int j;
	Matrix_cs *A = dgC2cs(a, 1);
	MCS_SET_XTYPE(A->xtype);
	if (!asLogical(sparse)) {
		const char *cl = (isNull(b)) ? "dtrMatrix" : "dgeMatrix";
		PROTECT(r = newObject(cl));

		SEXP rdim = GET_SLOT(r, Matrix_DimSym);
		int *prdim = INTEGER(rdim);
		prdim[0] = m;
		prdim[1] = n;

		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP rx = PROTECT(allocVector(REALSXP, mn));
		double *prx = REAL(rx);
		if (isNull(b)) {
			Matrix_memset(prx, 0, mn, sizeof(double));
			for (j = 0; j < n; ++j) {
				prx[j] = 1.0;
				if (ul == 'U')
					Matrix_cs_usolve(A, prx);
				else
					Matrix_cs_lsolve(A, prx);
				prx += m;
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			double *pbx = REAL(bx);
			Matrix_memcpy(prx, pbx, mn, sizeof(double));
			UNPROTECT(1); /* bx */
			for (j = 0; j < n; ++j) {
				if (ul == 'U')
					Matrix_cs_usolve(A, prx);
				else
					Matrix_cs_lsolve(A, prx);
				prx += m;
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	} else {
		Matrix_cs *B, *X;

		if (isNull(b)) {
			B = Matrix_cs_spalloc(m, n, n, 1, 0);
			if (!B)
				ERROR_SOLVE_OOM(".tCMatrix", ".gCMatrix");
			double *B__x = (double *) B->x;
			for (j = 0; j < n; ++j) {
				B->p[j] = j;
				B->i[j] = j;
				B__x[j] = 1.0;
			}
			B->p[n] = n;
		} else
			B = dgC2cs(b, 1);

		int i, k, top, nz, nzmax,
			*iwork = (int *) R_alloc((size_t) 2 * m, sizeof(int));
		double *work = (double *) R_alloc((size_t) m, sizeof(double));

		DO_TRIANGULAR_SOLVE(A, ul != 'U', isNull(b), ".tCMatrix", ".gCMatrix");

		/* Drop zeros from B and sort it : */
		Matrix_cs_dropzeros(B);
		X = Matrix_cs_transpose(B, 1);
		B = Matrix_cs_spfree(B);
		if (!X)
			ERROR_SOLVE_OOM(".tCMatrix", ".gCMatrix");
		B = Matrix_cs_transpose(X, 1);
		X = Matrix_cs_spfree(X);
		if (!B)
			ERROR_SOLVE_OOM(".tCMatrix", ".gCMatrix");

		PROTECT(r = cs2dgC(B, 1, (isNull(b)) ? 't' : 'g'));
		B = Matrix_cs_spfree(B);
	}
	if (isNull(b))
		SET_SLOT(r, Matrix_uploSym, auplo);

	SOLVE_FINISH;

	UNPROTECT(2); /* r, auplo */
	return r;
}

SEXP sparseQR_matmult(SEXP qr, SEXP y, SEXP op, SEXP complete, SEXP yxjj)
{
	SEXP V = PROTECT(GET_SLOT(qr, Matrix_VSym)),
		beta = PROTECT(GET_SLOT(qr, Matrix_betaSym)),
		p = PROTECT(GET_SLOT(qr, Matrix_pSym));
	Matrix_cs *V_ = dgC2cs(V, 1);
	MCS_SET_XTYPE(V_->xtype);
	double *pbeta = REAL(beta);
	int m = V_->m, r = V_->n, n, i, j, op_ = asInteger(op),
		*pp = INTEGER(p), nprotect = 6;

	SEXP yx;
	double *pyx;
	if (isNull(y)) {
		n = (asLogical(complete)) ? m : r;

		R_xlen_t mn = (R_xlen_t) m * n, m1a = (R_xlen_t) m + 1;
		PROTECT(yx = allocVector(REALSXP, mn));
		pyx = REAL(yx);
		Matrix_memset(pyx, 0, mn, sizeof(double));

		if (isNull(yxjj)) {
			for (j = 0; j < n; ++j) {
				*pyx = 1.0;
				pyx += m1a;
			}
		} else if (TYPEOF(yxjj) == REALSXP && XLENGTH(yxjj) >= n) {
			double *pyxjj = REAL(yxjj);
			for (j = 0; j < n; ++j) {
				*pyx = *pyxjj;
				pyx += m1a;
				pyxjj += 1;
			}
		} else
			error(_("invalid '%s' to '%s'"), "yxjj", __func__);
	} else {
		SEXP ydim = PROTECT(GET_SLOT(y, Matrix_DimSym));
		int *pydim = INTEGER(ydim);
		if (pydim[0] != m)
			error(_("dimensions of '%s' and '%s' are inconsistent"),
			     "qr", "y");
		n = pydim[1];
		UNPROTECT(1); /* ydim */

		PROTECT(yx = GET_SLOT(y, Matrix_xSym));
	}
	pyx = REAL(yx);

	SEXP a = PROTECT(newObject("dgeMatrix")),
		adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
		ax = yx;
	int *padim = INTEGER(adim);
	padim[0] = (op_ != 0) ? m : r;
	padim[1] = n;
	if (!isNull(y) || padim[0] != m) {
		PROTECT(ax = allocVector(REALSXP, (R_xlen_t) padim[0] * padim[1]));
		++nprotect;
	}
	double *pax = REAL(ax), *work = NULL;
	if (op_ < 5)
		work = (double *) R_alloc((size_t) m, sizeof(double));

	switch (op_) {
	case 0: /* qr.coef : A = P2 R1^{-1} Q1' P1 y */
	{
		SEXP R = PROTECT(GET_SLOT(qr, Matrix_RSym)),
			q = PROTECT(GET_SLOT(qr, Matrix_qSym));
		Matrix_cs *R_ = dgC2cs(R, 1);
		if (V_->xtype != R_->xtype)
			error(_("'%s' and '%s' slots have different '%s'"),
			      "V", "R", "xtype");
		int *pq = (LENGTH(q)) ? INTEGER(q) : NULL;

		for (j = 0; j < n; ++j) {
			Matrix_cs_pvec(pp, pyx, work, m);
			for (i = 0; i < r; ++i)
				Matrix_cs_happly(V_, i, pbeta[i], work);
			Matrix_cs_usolve(R_, work);
			Matrix_cs_ipvec(pq, work, pax, r);
			pyx += m;
			pax += r;
		}

		UNPROTECT(2); /* q, R */
		break;
	}
	case 1: /* qr.fitted : A = P1' Q1 Q1' P1 y */
		for (j = 0; j < n; ++j) {
			Matrix_cs_pvec(pp, pyx, work, m);
			for (i = 0; i < r; ++i)
				Matrix_cs_happly(V_, i, pbeta[i], work);
			if (r < m)
				Matrix_memset(work + r, 0, m - r, sizeof(double));
			for (i = r - 1; i >= 0; --i)
				Matrix_cs_happly(V_, i, pbeta[i], work);
			Matrix_cs_ipvec(pp, work, pax, m);
			pyx += m;
			pax += m;
		}
		break;
	case 2: /* qr.resid : A = P1' Q2 Q2' P1 y */
		for (j = 0; j < n; ++j) {
			Matrix_cs_pvec(pp, pyx, work, m);
			for (i = 0; i < r; ++i)
				Matrix_cs_happly(V_, i, pbeta[i], work);
			if (r > 0)
				Matrix_memset(work, 0, r, sizeof(double));
			for (i = r - 1; i >= 0; --i)
				Matrix_cs_happly(V_, i, pbeta[i], work);
			Matrix_cs_ipvec(pp, work, pax, m);
			pyx += m;
			pax += m;
		}
		break;
	case 3: /* qr.qty {w/ perm.} : A = Q' P1 y */
		for (j = 0; j < n; ++j) {
			Matrix_cs_pvec(pp, pyx, work, m);
			Matrix_memcpy(pax, work, m, sizeof(double));
			for (i = 0; i < r; ++i)
				Matrix_cs_happly(V_, i, pbeta[i], pax);
			pyx += m;
			pax += m;
		}
		break;
	case 4: /* qr.qy {w/ perm.} : A = P1' Q y */
		for (j = 0; j < n; ++j) {
			Matrix_memcpy(work, pyx, m, sizeof(double));
			for (i = r - 1; i >= 0; --i)
				Matrix_cs_happly(V_, i, pbeta[i], work);
			Matrix_cs_ipvec(pp, work, pax, m);
			pyx += m;
			pax += m;
		}
	break;
	case 5: /* qr.qty {w/o perm.} : A = Q' y */
		if (ax != yx)
			Matrix_memcpy(pax, pyx, (R_xlen_t) m * n, sizeof(double));
		for (j = 0; j < n; ++j) {
			for (i = 0; i < r; ++i)
				Matrix_cs_happly(V_, i, pbeta[i], pax);
			pax += m;
		}
	break;
	case 6: /* qr.qy {w/o perm.} : A = Q y */
		if (ax != yx)
			Matrix_memcpy(pax, pyx, (R_xlen_t) m * n, sizeof(double));
		for (j = 0; j < n; ++j) {
			for (i = r - 1; i >= 0; --i)
				Matrix_cs_happly(V_, i, pbeta[i], pax);
			pax += m;
		}
	break;
	default:
		error(_("invalid '%s' to '%s'"), "op", __func__);
		break;
	}

	SET_SLOT(a, Matrix_xSym, ax);
	UNPROTECT(nprotect); /* ax, adim, a, yx, p, beta, V */
	return a;
}
