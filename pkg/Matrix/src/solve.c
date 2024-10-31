#include "Lapack-etc.h"
#include "cs-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "solve.h"

static
void solveDN(SEXP rdn, SEXP adn, SEXP bdn)
{
	SEXP s;
	if ((s = VECTOR_ELT(adn, 1)) != R_NilValue)
		SET_VECTOR_ELT(rdn, 0, s);
	if ((s = VECTOR_ELT(bdn, 1)) != R_NilValue)
		SET_VECTOR_ELT(rdn, 1, s);
	PROTECT(adn = Rf_getAttrib(adn, R_NamesSymbol));
	PROTECT(bdn = Rf_getAttrib(bdn, R_NamesSymbol));
	if (adn != R_NilValue || bdn != R_NilValue) {
		PROTECT(s = Rf_allocVector(STRSXP, 2));
		if (adn != R_NilValue)
			SET_STRING_ELT(s, 0, STRING_ELT(adn, 1));
		if (bdn != R_NilValue)
			SET_STRING_ELT(s, 1, STRING_ELT(bdn, 1));
		Rf_setAttrib(rdn, R_NamesSymbol, s);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return;
}

SEXP denseLU_solve(SEXP s_a, SEXP s_b)
{

#define SOLVE_START \
	int *padim = DIM(s_a), m = padim[0], n = padim[1]; \
	if (m != n) \
		Rf_error(_("'%s' is not square"), "a"); \
	if (s_b != R_NilValue) { \
	int *pbdim = DIM(s_b); \
	if (pbdim[0] != m) \
		Rf_error(_("dimensions of '%s' and '%s' are inconsistent"), \
		         "a", "b"); \
	n = pbdim[1]; \
	}

#define SOLVE_FINISH \
	SEXP rdimnames = PROTECT(DIMNAMES(r, 0)), \
		adimnames = PROTECT(DIMNAMES(s_a, 0)); \
	if (s_b == R_NilValue) \
		cpyDN(rdimnames, adimnames, 1); \
	else { \
		SEXP bdimnames = PROTECT(DIMNAMES(s_b, 0)); \
		solveDN(rdimnames, adimnames, bdimnames); \
		UNPROTECT(1); /* bdimnames */ \
	} \
	UNPROTECT(2); /* adimnames, rdimnames */

	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(s_a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, m, n);

	if (m > 0) {
		SEXP apivot = PROTECT(GET_SLOT(s_a, Matrix_permSym)), rx;
		int info;
		if (s_b == R_NilValue) {
			rx = duplicateVector(ax);
			PROTECT(rx);
			int lwork = -1;
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
			double   work0, *work = &work0;
			F77_CALL(dgetri)(&m,    REAL(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_1(dgetri, info);
			lwork = (int) work0;
			work = (double   *) R_alloc((size_t) lwork, sizeof(double  ));
			F77_CALL(dgetri)(&m,    REAL(rx), &m, INTEGER(apivot),
			                 work, &lwork, &info);
			ERROR_LAPACK_2(dgetri, info, 2, U);
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(s_b, Matrix_xSym));
			rx = duplicateVector(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			if (TYPEOF(ax) == CPLXSXP) {
			F77_CALL(zgetrs)("N", &m, &n, COMPLEX(ax), &m, INTEGER(apivot),
			                 COMPLEX(rx), &m, &info FCONE);
			ERROR_LAPACK_1(zgetrs, info);
			} else {
			F77_CALL(dgetrs)("N", &m, &n,    REAL(ax), &m, INTEGER(apivot),
			                    REAL(rx), &m, &info FCONE);
			ERROR_LAPACK_1(dgetrs, info);
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(2); /* rx, apivot */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP denseBunchKaufman_solve(SEXP s_a, SEXP s_b)
{
	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(s_a, Matrix_xSym));
	int packed = XLENGTH(ax) != (int_fast64_t) m * m;

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	if (s_b == R_NilValue) {
		rcl[1] = 's';
		rcl[2] = (packed) ? 'p' : 'y';
	} else {
		rcl[1] = 'g';
		rcl[2] = 'e';
	}
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, m, n);

	char aul = UPLO(s_a);
	if (s_b == R_NilValue && aul != 'U')
		SET_UPLO(r);

	char act = (TYPEOF(ax) == CPLXSXP) ? TRANS(s_a) : 'C';
	if (s_b == R_NilValue && act != 'C')
		SET_TRANS(r);

	if (m > 0) {
		SEXP apivot = PROTECT(GET_SLOT(s_a, Matrix_permSym)), rx;
		int info;
		if (s_b == R_NilValue) {
			rx = duplicateVector(ax);
			PROTECT(rx);
			if (TYPEOF(ax) == CPLXSXP) {
			Rcomplex *work = (Rcomplex *) R_alloc((size_t) m, sizeof(Rcomplex));
			if (act == 'C') {
			if (!packed) {
				F77_CALL(zhetri)(&aul, &m, COMPLEX(rx), &m, INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(zhetri, info, 2, D);
			} else {
				F77_CALL(zhptri)(&aul, &m, COMPLEX(rx),     INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(zhptri, info, 2, D);
			}
			} else {
			if (!packed) {
				F77_CALL(zsytri)(&aul, &m, COMPLEX(rx), &m, INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(zsytri, info, 2, D);
			} else {
				F77_CALL(zsptri)(&aul, &m, COMPLEX(rx),     INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(zsptri, info, 2, D);
			}
			}
			} else {
			double   *work = (double   *) R_alloc((size_t) m, sizeof(double  ));
			if (!packed) {
				F77_CALL(dsytri)(&aul, &m, REAL(rx), &m, INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(dsytri, info, 2, D);
			} else {
				F77_CALL(dsptri)(&aul, &m, REAL(rx),     INTEGER(apivot),
				                 work, &info FCONE);
				ERROR_LAPACK_2(dsptri, info, 2, D);
			}
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(s_b, Matrix_xSym));
			rx = duplicateVector(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			if (TYPEOF(ax) == CPLXSXP) {
			if (act == 'C') {
			if (!packed) {
				F77_CALL(zhetrs)(&aul, &m, &n, COMPLEX(ax), &m, INTEGER(apivot),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zhetrs, info);
			} else {
				F77_CALL(zhptrs)(&aul, &m, &n, COMPLEX(ax),     INTEGER(apivot),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zhptrs, info);
			}
			} else {
			if (!packed) {
				F77_CALL(zsytrs)(&aul, &m, &n, COMPLEX(ax), &m, INTEGER(apivot),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zsytrs, info);
			} else {
				F77_CALL(zsptrs)(&aul, &m, &n, COMPLEX(ax),     INTEGER(apivot),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zsptrs, info);
			}
			}
			} else {
			if (!packed) {
				F77_CALL(dsytrs)(&aul, &m, &n,    REAL(ax), &m, INTEGER(apivot),
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dsytrs, info);
			} else {
				F77_CALL(dsptrs)(&aul, &m, &n,    REAL(ax),    INTEGER(apivot),
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dsptrs, info);
			}
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(2); /* rx, apivot */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP denseCholesky_solve(SEXP s_a, SEXP s_b)
{
	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(s_a, Matrix_xSym));
	int packed = XLENGTH(ax) != (int_fast64_t) m * m;

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	if (s_b == R_NilValue) {
		rcl[1] = 'p';
		rcl[2] = (packed) ? 'p' : 'o';
	} else {
		rcl[1] = 'g';
		rcl[2] = 'e';
	}
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, m, n);

	char aul = UPLO(s_a);;
	if (s_b == R_NilValue && aul != 'U')
		SET_UPLO(r);

	if (m > 0) {
		SEXP rx, aperm = PROTECT(Rf_getAttrib(s_a, Matrix_permSym));
		int info, pivoted = TYPEOF(aperm) == INTSXP && LENGTH(aperm) > 0;
		if (s_b == R_NilValue) {
			PROTECT(rx = Rf_allocVector(TYPEOF(ax), XLENGTH(ax)));
			if (TYPEOF(ax) == CPLXSXP) {
			memcpy(COMPLEX(rx), COMPLEX(ax), sizeof(Rcomplex) * (size_t) XLENGTH(ax));
			if (!packed) {
				F77_CALL(zpotri)(&aul, &m, COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_2(zpotri, info, 2, L);
				zsymperm2(COMPLEX(rx), NULL,
				          m, aul, (pivoted) ? INTEGER(aperm) : NULL, 1, 1);
			} else {
				F77_CALL(zpptri)(&aul, &m, COMPLEX(rx),     &info FCONE);
				ERROR_LAPACK_2(zpptri, info, 2, L);
				zsymperm1(COMPLEX(rx), NULL,
				          m, aul, (pivoted) ? INTEGER(aperm) : NULL, 1, 1);
			}
			} else {
			memcpy(   REAL(rx),    REAL(ax), sizeof(  double) * (size_t) XLENGTH(ax));
			if (!packed) {
				F77_CALL(dpotri)(&aul, &m,    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_2(dpotri, info, 2, L);
				dsymperm2(   REAL(rx), NULL,
				          m, aul, (pivoted) ? INTEGER(aperm) : NULL, 1, 1);
			} else {
				F77_CALL(dpptri)(&aul, &m,    REAL(rx),     &info FCONE);
				ERROR_LAPACK_2(dpptri, info, 2, L);
				dsymperm1(   REAL(rx), NULL,
				          m, aul, (pivoted) ? INTEGER(aperm) : NULL, 1, 1);
			}
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(s_b, Matrix_xSym));
			PROTECT(rx = Rf_allocVector(TYPEOF(ax), XLENGTH(bx)));
			if (TYPEOF(ax) == CPLXSXP) {
			zrowperm2(COMPLEX(rx), COMPLEX(bx),
			          m, n, (pivoted) ? INTEGER(aperm) : NULL, 1, 0);
			if (!packed) {
				F77_CALL(zpotrs)(&aul, &m, &n, COMPLEX(ax), &m,
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zpotrs, info);
			} else {
				F77_CALL(zpptrs)(&aul, &m, &n, COMPLEX(ax),
				                 COMPLEX(rx), &m, &info FCONE);
				ERROR_LAPACK_1(zpptrs, info);
			}
			zrowperm2(COMPLEX(rx), NULL,
			          m, n, (pivoted) ? INTEGER(aperm) : NULL, 1, 1);
			} else {
			drowperm2(   REAL(rx),    REAL(bx),
			          m, n, (pivoted) ? INTEGER(aperm) : NULL, 1, 0);
			if (!packed) {
				F77_CALL(dpotrs)(&aul, &m, &n,    REAL(ax), &m,
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dpotrs, info);
			} else {
				F77_CALL(dpptrs)(&aul, &m, &n,    REAL(ax),
				                    REAL(rx), &m, &info FCONE);
				ERROR_LAPACK_1(dpptrs, info);
			}
			drowperm2(   REAL(rx), NULL,
			          m, n, (pivoted) ? INTEGER(aperm) : NULL, 1, 1);
			}
			UNPROTECT(1); /* bx */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(2); /* rx, aperm */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP trMatrix_solve(SEXP s_a, SEXP s_b)
{
	SOLVE_START;

	SEXP ax = PROTECT(GET_SLOT(s_a, Matrix_xSym));
	int packed = XLENGTH(ax) != (int_fast64_t) m * m;

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	if (s_b == R_NilValue) {
		rcl[1] = 't';
		rcl[2] = (packed) ? 'p' : 'r';
	} else {
		rcl[1] = 'g';
		rcl[2] = 'e';
	}
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, m, n);

	char aul = UPLO(s_a);
	if (s_b == R_NilValue && aul != 'U')
		SET_UPLO(r);

	char anu = DIAG(s_a);
	if (s_b == R_NilValue && anu != 'N')
		SET_DIAG(r);

	if (m > 0) {
		SEXP rx;
		int info;
		if (s_b == R_NilValue) {
			rx = duplicateVector(ax);
			PROTECT(rx);
			if (TYPEOF(ax) == CPLXSXP) {
			if (!packed) {
				F77_CALL(ztrtri)(&aul, &anu, &m, COMPLEX(rx), &m,
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(ztrtri, info, 2, A);
			} else {
				F77_CALL(ztptri)(&aul, &anu, &m, COMPLEX(rx),
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(ztptri, info, 2, A);
			}
			} else {
			if (!packed) {
				F77_CALL(dtrtri)(&aul, &anu, &m,    REAL(rx), &m,
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(dtrtri, info, 2, A);
			} else {
				F77_CALL(dtptri)(&aul, &anu, &m,    REAL(rx),
				                 &info FCONE FCONE);
				ERROR_LAPACK_2(dtptri, info, 2, A);
			}
			}
		} else {
			SEXP bx = PROTECT(GET_SLOT(s_b, Matrix_xSym));
			rx = duplicateVector(bx);
			UNPROTECT(1); /* bx */
			PROTECT(rx);
			if (TYPEOF(ax) == CPLXSXP) {
			if (!packed) {
				F77_CALL(ztrtrs)(&aul, "N", &anu, &m, &n, COMPLEX(ax), &m,
				                 COMPLEX(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(ztrtrs, info);
			} else {
				F77_CALL(ztptrs)(&aul, "N", &anu, &m, &n, COMPLEX(ax),
				                 COMPLEX(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(ztptrs, info);
			}
			} else {
			if (!packed) {
				F77_CALL(dtrtrs)(&aul, "N", &anu, &m, &n,    REAL(ax), &m,
				                    REAL(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(dtrtrs, info);
			} else {
				F77_CALL(dtptrs)(&aul, "N", &anu, &m, &n,    REAL(ax),
				                    REAL(rx), &m, &info FCONE FCONE FCONE);
				ERROR_LAPACK_1(dtptrs, info);
			}
			}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	}

	SOLVE_FINISH;

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP sparseLU_solve(SEXP s_a, SEXP s_b, SEXP s_sparse)
{
	SOLVE_START;

	SEXP r,
		aL = PROTECT(GET_SLOT(s_a, Matrix_LSym)),
		aU = PROTECT(GET_SLOT(s_a, Matrix_USym)),
		ap = PROTECT(GET_SLOT(s_a, Matrix_pSym)),
		aq = PROTECT(GET_SLOT(s_a, Matrix_qSym));
	int i, j,
		*pap = (LENGTH(ap)) ? INTEGER(ap) : NULL,
		*paq = (LENGTH(aq)) ? INTEGER(aq) : NULL;
	Matrix_cs *L = M2CXS(aL, 1), *U = M2CXS(aU, 1);
	CXSPARSE_XTYPE_SET(L->xtype);
	if (!Rf_asLogical(s_sparse)) {
		if ((int_fast64_t) m * n > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		char rcl[] = ".geMatrix";
		rcl[0] = (L->xtype == CXSPARSE_COMPLEX) ? 'z' : 'd';
		PROTECT(r = newObject(rcl));
		SET_DIM(r, m, n);
		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP rx = PROTECT(Rf_allocVector((L->xtype == CXSPARSE_COMPLEX) ? CPLXSXP : REALSXP, mn));
		if (s_b == R_NilValue) {

#define SOLVE_DENSE(c) \
			do { \
				c##TYPE *prx = c##PTR(rx), \
					*work = (c##TYPE *) R_alloc((size_t) m, sizeof(c##TYPE)); \
				memset(prx, 0, sizeof(c##TYPE) * (size_t) mn); \
				for (j = 0; j < n; ++j) { \
					prx[j] = c##UNIT; \
					Matrix_cs_pvec(pap, prx, work, m); \
					Matrix_cs_lsolve(L, work); \
					Matrix_cs_usolve(U, work); \
					Matrix_cs_ipvec(paq, work, prx, m); \
					prx += m; \
				} \
			} while (0)

			if (L->xtype == CXSPARSE_COMPLEX)
			SOLVE_DENSE(z);
			else
			SOLVE_DENSE(d);

#undef SOLVE_DENSE

		} else {
			SEXP bx = PROTECT(GET_SLOT(s_b, Matrix_xSym));

#define SOLVE_DENSE(c) \
			do { \
				c##TYPE *prx = c##PTR(rx), *pbx = c##PTR(bx), \
					*work = (c##TYPE *) R_alloc((size_t) m, sizeof(c##TYPE)); \
				for (j = 0; j < n; ++j) { \
					Matrix_cs_pvec(pap, pbx, work, m); \
					Matrix_cs_lsolve(L, work); \
					Matrix_cs_usolve(U, work); \
					Matrix_cs_ipvec(paq, work, prx, m); \
					prx += m; \
					pbx += m; \
				} \
			} while (0)

			if (L->xtype == CXSPARSE_COMPLEX)
			SOLVE_DENSE(z);
			else
			SOLVE_DENSE(d);

#undef SOLVE_DENSE

			UNPROTECT(1); /* bx */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	} else {
		Matrix_cs *B = NULL, *X = NULL;
		if (s_b == R_NilValue) {
			B = Matrix_cs_speye(m, m, 1, 0);
			if (B && pap)
				for (i = 0; i < m; ++i)
					B->i[pap[i]] = i;
		} else {
			B = M2CXS(s_b, 1);
			if (B && pap) {
				int *papinv = Matrix_cs_pinv(pap, m);
				if (!papinv)
					ERROR_OOM(__func__);
				B = Matrix_cs_permute(B, papinv, NULL, 1);
				papinv = Matrix_cs_free(papinv);
			}
		}
		if (!B)
			ERROR_OOM(__func__);
		int Bfr = s_b == R_NilValue || pap;

		int k, top, nz, nzmax,
			*iwork = (int *) R_alloc((size_t) m * 2, sizeof(int));

#define SOLVE_SPARSE_TRIANGULAR(c, _A_, _ALO_, _BFR_) \
		do { \
			X = Matrix_cs_spalloc(m, n, B->nzmax, 1, 0); \
			if (!X) { \
				if (_BFR_) \
					B = Matrix_cs_spfree(B); \
				ERROR_OOM(__func__); \
			} \
			c##TYPE *X__x = (c##TYPE *) X->x; \
			X->p[0] = nz = 0; \
			nzmax = X->nzmax; \
			for (j = 0, k = 0; j < n; ++j) { \
				top = Matrix_cs_spsolve(_A_, B, j, iwork, work, NULL, _ALO_); \
				if (m - top > INT_MAX - nz) { \
					if (_BFR_) \
						B = Matrix_cs_spfree(B); \
					X = Matrix_cs_spfree(X); \
					Rf_error(_("attempt to construct %s with more than %s nonzero elements"), \
					         "sparseMatrix", "2^31-1"); \
				} \
				nz += m - top; \
				if (nz > nzmax) { \
					nzmax = (nz <= INT_MAX / 2) ? 2 * nz : INT_MAX; \
					if (!Matrix_cs_sprealloc(X, nzmax)) { \
						if (_BFR_) \
							B = Matrix_cs_spfree(B); \
						X = Matrix_cs_spfree(X); \
						ERROR_OOM(__func__); \
					} \
					X__x = (c##TYPE *) X->x; \
				} \
				X->p[j + 1] = nz; \
				if (_ALO_) { \
					for (i = top; i < m; ++i) { \
						X->i[k] =      iwork[i] ; \
						X__x[k] = work[iwork[i]]; \
						++k; \
					} \
				} else { \
					for (i = m - 1; i >= top; --i) { \
						X->i[k] =      iwork[i] ; \
						X__x[k] = work[iwork[i]]; \
						++k; \
					} \
				} \
			} \
			if (_BFR_) \
				B = Matrix_cs_spfree(B); \
			B = X; \
		} while (0)

#define SOLVE_SPARSE(c) \
		do { \
			c##TYPE *work = (c##TYPE *) R_alloc((size_t) m, sizeof(c##TYPE)); \
			SOLVE_SPARSE_TRIANGULAR(c, L, 1, Bfr); \
			SOLVE_SPARSE_TRIANGULAR(c, U, 0,   1); \
		} while (0)

		if (L->xtype == CXSPARSE_COMPLEX)
		SOLVE_SPARSE(z);
		else
		SOLVE_SPARSE(d);

#undef SOLVE_SPARSE

		if (paq) {
			X = Matrix_cs_permute(B, paq, NULL, 1);
			B = Matrix_cs_spfree(B);
			if (!X)
				ERROR_OOM(__func__);
			B = X;
		}

		/* Drop zeros from B and sort it : */
		Matrix_cs_dropzeros(B);
		X = Matrix_cs_transpose(B, 1);
		B = Matrix_cs_spfree(B);
		if (!X)
			ERROR_OOM(__func__);
		B = Matrix_cs_transpose(X, 1);
		X = Matrix_cs_spfree(X);
		if (!B)
			ERROR_OOM(__func__);

		PROTECT(r = CXS2M(B, 1, 'g'));
		B = Matrix_cs_spfree(B);
	}

	SOLVE_FINISH;

	UNPROTECT(5); /* r, aq, ap, aU, aL */
	return r;
}

static
int strmatch(const char *s, const char **nms)
{
	int i = 0;
	while (nms[i][0] != '\0') {
		if (strcmp(s, nms[i]) == 0)
			return i;
		++i;
	}
	return -1;
}

SEXP sparseCholesky_solve(SEXP s_a, SEXP s_b, SEXP s_sparse, SEXP s_system)
{
	static const char *valid[] = {
		"A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt", "" };
	int ivalid = -1;
	if (TYPEOF(s_system) != STRSXP || LENGTH(s_system) < 1 ||
	    (s_system = STRING_ELT(s_system, 0)) == NA_STRING ||
	    (ivalid = strmatch(CHAR(s_system), valid)) < 0)
		Rf_error(_("invalid '%s' to '%s'"), "system", __func__);

	SOLVE_START;

	SEXP r;
	size_t m_ = (size_t) m, n_ = (size_t) n;
	cholmod_factor *L = M2CHF(s_a, 1);
	if (!Rf_asLogical(s_sparse)) {
		if ((int_fast64_t) m * n > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		cholmod_dense *B = NULL, *X = NULL;
		if (s_b == R_NilValue) {
			B = cholmod_eye(m_, n_, L->xtype + L->dtype, &c);
			if (!B)
				ERROR_OOM(__func__);
			X = cholmod_solve(ivalid, L, B, &c);
			cholmod_free_dense(&B, &c);
			if (!X)
				ERROR_OOM(__func__);
			PROTECT(r = CHD2M(X, 'N',
				(ivalid < 2) ? 'p' : ((ivalid < 7) ? 't' : 'g')));
		} else {
			B = M2CHD(s_b, 'N');
			X = cholmod_solve(ivalid, L, B, &c);
			if (!X)
				ERROR_OOM(__func__);
			PROTECT(r = CHD2M(X, 'N', 'g'));
		}
		cholmod_free_dense(&X, &c);
	} else {
		cholmod_sparse *B = NULL, *X = NULL;
		if (s_b == R_NilValue) {
			B = cholmod_speye(m_, n_, L->xtype + L->dtype, &c);
			if (!B)
				ERROR_OOM(__func__);
			X = cholmod_spsolve(ivalid, L, B, &c);
			cholmod_free_sparse(&B, &c);
			if (X && ivalid < 7) {
				if (!X->sorted)
					cholmod_sort(X, &c);
				B = cholmod_copy(X, (ivalid == 2 || ivalid == 4) ? -1 : 1, 1, &c);
				cholmod_free_sparse(&X, &c);
				X = B;
			}
			if (!X)
				ERROR_OOM(__func__);
			PROTECT(r = CHS2M(X, 1,
				(ivalid < 2) ? 'p' : ((ivalid < 7) ? 't' : 'g')));
		} else {
			B = M2CHS(s_b, 1);
			X = cholmod_spsolve(ivalid, L, B, &c);
			if (!X)
				ERROR_OOM(__func__);
			PROTECT(r = CHS2M(X, 1, 'g'));
		}
		cholmod_free_sparse(&X, &c);
	}
	if (s_b == R_NilValue && (ivalid == 2 || ivalid == 4))
		SET_UPLO(r);

	SOLVE_FINISH;

	UNPROTECT(1); /* r */
	return r;
}

SEXP tCMatrix_solve(SEXP s_a, SEXP s_b, SEXP s_sparse)
{
	SOLVE_START;

	SEXP r;
	char aul = UPLO(s_a);
	int i, j;
	Matrix_cs *A = M2CXS(s_a, 1);
	CXSPARSE_XTYPE_SET(A->xtype);
	if (!Rf_asLogical(s_sparse)) {
		if ((int_fast64_t) m * n > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		char rcl[] = "...Matrix";
		rcl[0] = (A->xtype == CXSPARSE_COMPLEX) ? 'z' : 'd';
		rcl[1] = (s_b == R_NilValue) ? 't' : 'g';
		rcl[2] = (s_b == R_NilValue) ? 'r' : 'e';
		PROTECT(r = newObject(rcl));
		SET_DIM(r, m, n);
		R_xlen_t mn = (R_xlen_t) m * n;
		SEXP rx = PROTECT(Rf_allocVector((A->xtype == CXSPARSE_COMPLEX) ? CPLXSXP : REALSXP, mn));
		if (s_b == R_NilValue) {

#define SOLVE_DENSE(c) \
			do { \
				c##TYPE *prx = c##PTR(rx); \
				memset(prx, 0, sizeof(c##TYPE) * (size_t) mn); \
				for (j = 0; j < n; ++j) { \
					prx[j] = c##UNIT; \
					if (aul == 'U') \
						Matrix_cs_usolve(A, prx); \
					else \
						Matrix_cs_lsolve(A, prx); \
					prx += m; \
				} \
			} while (0)

			if (A->xtype == CXSPARSE_COMPLEX)
			SOLVE_DENSE(z);
			else
			SOLVE_DENSE(d);

#undef SOLVE_DENSE

		} else {
			SEXP bx = PROTECT(GET_SLOT(s_b, Matrix_xSym));

#define SOLVE_DENSE(c) \
			do { \
				c##TYPE *prx = c##PTR(rx), *pbx = c##PTR(bx); \
				memcpy(prx, pbx, sizeof(c##TYPE) * (size_t) mn); \
				for (j = 0; j < n; ++j) { \
					if (aul == 'U') \
						Matrix_cs_usolve(A, prx); \
					else \
						Matrix_cs_lsolve(A, prx); \
					prx += m; \
				} \
			} while (0)

			if (A->xtype == CXSPARSE_COMPLEX)
			SOLVE_DENSE(z);
			else
			SOLVE_DENSE(d);

#undef SOLVE_DENSE

			UNPROTECT(1); /* bx */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
	} else {
		Matrix_cs *B = NULL, *X = NULL;
		if (s_b == R_NilValue)
			B = Matrix_cs_speye(m, m, 1, 0);
		else
			B = M2CXS(s_b, 1);
		if (!B)
			ERROR_OOM(__func__);

		int k, top, nz, nzmax,
			*iwork = (int *) R_alloc((size_t) m * 2, sizeof(int));

#define SOLVE_SPARSE(c) \
		do { \
			c##TYPE *work = (c##TYPE *) R_alloc((size_t) m, sizeof(c##TYPE)); \
			SOLVE_SPARSE_TRIANGULAR(c, A, aul != 'U', s_b == R_NilValue); \
		} while (0)

		if (A->xtype == CXSPARSE_COMPLEX)
		SOLVE_SPARSE(z);
		else
		SOLVE_SPARSE(d);

#undef SOLVE_SPARSE

		/* Drop zeros from B and sort it : */
		Matrix_cs_dropzeros(B);
		X = Matrix_cs_transpose(B, 1);
		B = Matrix_cs_spfree(B);
		if (!X)
			ERROR_OOM(__func__);
		B = Matrix_cs_transpose(X, 1);
		X = Matrix_cs_spfree(X);
		if (!B)
			ERROR_OOM(__func__);

		PROTECT(r = CXS2M(B, 1, (s_b == R_NilValue) ? 't' : 'g'));
		B = Matrix_cs_spfree(B);
	}
	if (s_b == R_NilValue && aul != 'U')
		SET_UPLO(r);

	SOLVE_FINISH;

	UNPROTECT(1); /* r */
	return r;
}

SEXP sparseQR_matmult(SEXP s_qr, SEXP s_y, SEXP s_op,
                      SEXP s_complete, SEXP s_yxjj)
{
	SEXP V = PROTECT(GET_SLOT(s_qr, Matrix_VSym));
	Matrix_cs *V_ = M2CXS(V, 1);
	CXSPARSE_XTYPE_SET(V_->xtype);

	SEXP beta = PROTECT(GET_SLOT(s_qr, Matrix_betaSym));
	double *pbeta = REAL(beta);

	SEXP p = PROTECT(GET_SLOT(s_qr, Matrix_pSym));
	int *pp = (LENGTH(p) > 0) ? INTEGER(p) : NULL;

	int m = V_->m, r = V_->n, n, i, j, op = Rf_asInteger(s_op), nprotect = 5;

	SEXP yx;
	if (s_y == R_NilValue) {
		n = (Rf_asLogical(s_complete)) ? m : r;
		if ((int_fast64_t) m * n > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");
		R_xlen_t mn = (R_xlen_t) m * n, m1a = (R_xlen_t) m + 1;
		PROTECT(yx = Rf_allocVector((V_->xtype == CXSPARSE_COMPLEX) ? CPLXSXP : REALSXP, mn));

#define EYE(c) \
		do { \
			c##TYPE *pyx = c##PTR(yx); \
			memset(pyx, 0, sizeof(c##TYPE) * (size_t) mn); \
			if (s_yxjj == R_NilValue) { \
				for (j = 0; j < n; ++j) { \
					*pyx = c##UNIT; \
					pyx += m1a; \
				} \
			} else if (TYPEOF(s_yxjj) == TYPEOF(yx) && XLENGTH(s_yxjj) >= n) { \
				c##TYPE *pyxjj = c##PTR(s_yxjj); \
				for (j = 0; j < n; ++j) { \
					*pyx = *pyxjj; \
					pyx += m1a; \
					pyxjj += 1; \
				} \
			} else \
				Rf_error(_("invalid '%s' to '%s'"), "yxjj", __func__); \
		} while (0)

		if (V_->xtype == CXSPARSE_COMPLEX)
		EYE(z);
		else
		EYE(d);

#undef EYE

	} else {
		int *pydim = DIM(s_y);
		if (pydim[0] != m)
			Rf_error(_("dimensions of '%s' and '%s' are inconsistent"),
			         "qr", "y");
		n = pydim[1];

		PROTECT(yx = GET_SLOT(s_y, Matrix_xSym));
	}

	char acl[] = ".geMatrix";
	acl[0] = (V_->xtype == CXSPARSE_COMPLEX) ? 'z' : 'd';
	SEXP a = PROTECT(newObject(acl));

	int *padim = DIM(a);
	padim[0] = (op != 0) ? m : r;
	padim[1] = n;

	SEXP ax;
	if (s_y == R_NilValue && padim[0] == m)
		ax = yx;
	else {
		R_xlen_t mn = (R_xlen_t) padim[0] * padim[1];
		PROTECT(ax = Rf_allocVector((V_->xtype == CXSPARSE_COMPLEX) ? CPLXSXP : REALSXP, mn));
		++nprotect;
	}
	SET_SLOT(a, Matrix_xSym, ax);

#define MATMULT(c) \
	do { \
		c##TYPE *pyx = c##PTR(yx), *pax = c##PTR(ax), *work = NULL; \
		if (op < 5) \
			work = (c##TYPE *) R_alloc((size_t) m, sizeof(c##TYPE)); \
		switch (op) { \
		case 0: /* qr.coef : A = P2 R1^{-1} Q1' P1 y */ \
		{ \
			SEXP R = PROTECT(GET_SLOT(s_qr, Matrix_RSym)), \
				q = PROTECT(GET_SLOT(s_qr, Matrix_qSym)); \
			Matrix_cs *R_ = M2CXS(R, 1); \
			int *pq = (LENGTH(q) > 0) ? INTEGER(q) : NULL; \
			for (j = 0; j < n; ++j) { \
				Matrix_cs_pvec(pp, pyx, work, m); \
				for (i = 0; i < r; ++i) \
					Matrix_cs_happly(V_, i, pbeta[i], work); \
				Matrix_cs_usolve(R_, work); \
				Matrix_cs_ipvec(pq, work, pax, r); \
				pyx += m; \
				pax += r; \
			} \
			UNPROTECT(2); /* q, R */ \
			break; \
		} \
		case 1: /* qr.fitted : A = P1' Q1 Q1' P1 y */ \
			for (j = 0; j < n; ++j) { \
				Matrix_cs_pvec(pp, pyx, work, m); \
				for (i = 0; i < r; ++i) \
					Matrix_cs_happly(V_, i, pbeta[i], work); \
				if (r < m) \
					memset(work + r, 0, sizeof(c##TYPE) * (size_t) (m - r)); \
				for (i = r - 1; i >= 0; --i) \
					Matrix_cs_happly(V_, i, pbeta[i], work); \
				Matrix_cs_ipvec(pp, work, pax, m); \
				pyx += m; \
				pax += m; \
			} \
			break; \
		case 2: /* qr.resid : A = P1' Q2 Q2' P1 y */ \
			for (j = 0; j < n; ++j) { \
				Matrix_cs_pvec(pp, pyx, work, m); \
				for (i = 0; i < r; ++i) \
					Matrix_cs_happly(V_, i, pbeta[i], work); \
				if (r > 0) \
					memset(work, 0, sizeof(c##TYPE) * (size_t) r); \
				for (i = r - 1; i >= 0; --i) \
					Matrix_cs_happly(V_, i, pbeta[i], work); \
				Matrix_cs_ipvec(pp, work, pax, m); \
				pyx += m; \
				pax += m; \
			} \
			break; \
		case 3: /* qr.qty {w/ perm.} : A = Q' P1 y */ \
			for (j = 0; j < n; ++j) { \
				Matrix_cs_pvec(pp, pyx, work, m); \
				memcpy(pax, work, sizeof(c##TYPE) * (size_t) m); \
				for (i = 0; i < r; ++i) \
					Matrix_cs_happly(V_, i, pbeta[i], pax); \
				pyx += m; \
				pax += m; \
			} \
			break; \
		case 4: /* qr.qy {w/ perm.} : A = P1' Q y */ \
			for (j = 0; j < n; ++j) { \
				memcpy(work, pyx, sizeof(c##TYPE) * (size_t) m); \
				for (i = r - 1; i >= 0; --i) \
					Matrix_cs_happly(V_, i, pbeta[i], work); \
				Matrix_cs_ipvec(pp, work, pax, m); \
				pyx += m; \
				pax += m; \
			} \
		break; \
		case 5: /* qr.qty {w/o perm.} : A = Q' y */ \
			if (ax != yx) \
				memcpy(pax, pyx, sizeof(c##TYPE) * (size_t) m * (size_t) n); \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < r; ++i) \
					Matrix_cs_happly(V_, i, pbeta[i], pax); \
				pax += m; \
			} \
		break; \
		case 6: /* qr.qy {w/o perm.} : A = Q y */ \
			if (ax != yx) \
				memcpy(pax, pyx, sizeof(c##TYPE) * (size_t) m * (size_t) n); \
			for (j = 0; j < n; ++j) { \
				for (i = r - 1; i >= 0; --i) \
					Matrix_cs_happly(V_, i, pbeta[i], pax); \
				pax += m; \
			} \
		break; \
		default: \
			Rf_error(_("invalid '%s' to '%s'"), "op", __func__); \
			break; \
		} \
	} while (0)

	if (V_->xtype == CXSPARSE_COMPLEX)
	MATMULT(z);
	else
	MATMULT(d);

#undef MATMULT

	UNPROTECT(nprotect); /* ax, a, yx, p, beta, V */
	return a;
}
