/* C implementation of methods for %&%, %*%, crossprod, tcrossprod */

#include "Lapack-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "M5.h"
#include "idz.h"
#include "coerce.h"

SEXP dense_transpose(SEXP, const char *, char);
SEXP sparse_transpose(SEXP, const char *, char, int);
SEXP sparse_dropzero(SEXP, const char *, double);
SEXP sparse_diag_U2N(SEXP, const char *);

static const char *valid_matmult[] = {
VALID_DIAGONAL, VALID_SPARSE, VALID_DENSE, "" };
/* ptr:      0,            5,          56, 87 */

static
void matmultDim(SEXP x, SEXP y, char *xtrans, char *ytrans, char *ztrans,
                int *m, int *n, int *v)
{
	int xt = *xtrans == 'C' || *xtrans == 'T'; if (!xt) *xtrans = 'N';
	int yt = *ytrans == 'C' || *ytrans == 'T'; if (!yt) *ytrans = 'N';
	int zt = *ztrans == 'C' || *ztrans == 'T'; if (!zt) *ztrans = 'N';
	if (y == R_NilValue) {
		if (xt == yt)
			Rf_error(_("should never happen ..."));
		SEXP
			xdim = (TYPEOF(x) == OBJSXP)
			? GET_SLOT(x, Matrix_DimSym) : Rf_getAttrib(x, R_DimSymbol);
		if (TYPEOF(xdim) == INTSXP && LENGTH(xdim) == 2) {
			*v = 0;
			*m = *n = INTEGER(xdim)[(xt) ? 1 : 0];
		} else if (XLENGTH(x) <= INT_MAX) {
			*v = 1;
			*m = *n = (xt) ? 1 : LENGTH(x);
		} else
			Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
	} else {
		/* MJ: So that I don't lose my mind ... : */
		if (zt) {
			SWAP(x, y, SEXP, );
			SWAP(xt, yt, int, !);
		}
		SEXP
			xdim = (TYPEOF(x) == OBJSXP)
			? GET_SLOT(x, Matrix_DimSym) : Rf_getAttrib(x, R_DimSymbol),
			ydim = (TYPEOF(y) == OBJSXP)
			? GET_SLOT(y, Matrix_DimSym) : Rf_getAttrib(y, R_DimSymbol);
		int xm, xn, ym, yn, x2, y2;
		xm = xn = ym = yn = -1;
		x2 = TYPEOF(xdim) == INTSXP && LENGTH(xdim) == 2;
		y2 = TYPEOF(ydim) == INTSXP && LENGTH(ydim) == 2;
		if (x2) {
			int *pxdim = INTEGER(xdim);
			xm = pxdim[0];
			xn = pxdim[1];
		} else if (XLENGTH(x) > INT_MAX)
			Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		if (y2) {
			int *pydim = INTEGER(ydim);
			ym = pydim[0];
			yn = pydim[1];
		} else if (XLENGTH(y) > INT_MAX)
			Rf_error(_("dimensions cannot exceed %s"), "2^31-1");
		/* MJ: R's do_matprod behaves quite asymmetrically ... what a pain */
		if (x2 && y2)
			*v = 0;
		else if (y2) {
			*v = (zt) ? 2 : 1;
			int k = (yt) ? yn : ym, xl = LENGTH(x);
			if (k == xl || (k == 1 && !(xt))) {
				xm = (int) xl;
				xn = 1;
				xt = (k == xl) ? 1 : 0;
			}
		} else if (x2) {
			*v = (zt) ? 1 : 2;
			int k = (xt) ? xm : xn, yl = LENGTH(y);
			if (yt) {
				if (xm == 1 || xn == 1) {
					ym = (int) yl;
					yn = 1;
					yt = (((xt) ? xn : xm) == 1) ? 0 : 1;
				}
			} else {
				if (k == yl || k == 1) {
					ym = (int) yl;
					yn = 1;
					yt = (k == yl) ? 0 : 1;
				}
			}
		} else {
			*v = 3;
			int xl = LENGTH(x), yl = LENGTH(y);
			if (xt) {
				xm = xl;
				xn = 1;
				ym = yl;
				yn = 1;
				yt = xl == 1;
			} else if (yt) {
				xm = xl;
				xn = 1;
				ym = yl;
				yn = 1;
				/* xt = 0; */
			} else {
				xm = 1;
				xn = xl;
				ym = (xl == 1) ? 1 : yl;
				yn = (xl == 1) ? yl : 1;
			}
		}
		if (((xt) ? xm : xn) != ((yt) ? yn : ym))
			Rf_error(_("non-conformable arguments"));
		*m = (xt) ? xn : xm;
		*n = (yt) ? ym : yn;
		if (zt) {
			SWAP(*m, *n, int, );
			SWAP(xt, yt, int, !);
		}
	}
	if (*v % 2) *xtrans = (xt) ? 'T' : 'N';
	if (*v > 1) *ytrans = (yt) ? 'T' : 'N';
	return;
}

static
void matmultDN(SEXP dest, SEXP asrc, int ai, SEXP bsrc, int bi) {
	SEXP s;
	if ((s = VECTOR_ELT(asrc, ai)) != R_NilValue)
		SET_VECTOR_ELT(dest, 0, s);
	if ((s = VECTOR_ELT(bsrc, bi)) != R_NilValue)
		SET_VECTOR_ELT(dest, 1, s);
	PROTECT(asrc = Rf_getAttrib(asrc, R_NamesSymbol));
	PROTECT(bsrc = Rf_getAttrib(bsrc, R_NamesSymbol));
	if (asrc != R_NilValue || bsrc != R_NilValue) {
		SEXP destnms = PROTECT(Rf_allocVector(STRSXP, 2));
		if (asrc != R_NilValue && CHAR(s = STRING_ELT(asrc, ai))[0] != '\0')
			SET_STRING_ELT(destnms, 0, s);
		if (bsrc != R_NilValue && CHAR(s = STRING_ELT(bsrc, bi))[0] != '\0')
			SET_STRING_ELT(destnms, 1, s);
		Rf_setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return;
}

#define CONJ2(_X_, _M_, _N_) \
do { \
	size_t m = (size_t) _M_, n = (size_t) _N_, xlen = m * n; \
	Rcomplex *x = (Rcomplex *) _X_; \
	Rcomplex *y = (Rcomplex *) R_alloc(xlen, sizeof(Rcomplex)); \
	zvconj(x, y, xlen); \
	_X_ = y; \
} while (0)

#define CONJ1(_X_, _N_) \
do { \
	size_t n = (size_t) _N_, xlen = PACKED_LENGTH(n); \
	Rcomplex *x = (Rcomplex *) _X_; \
	Rcomplex *y = (Rcomplex *) R_alloc(xlen, sizeof(Rcomplex)); \
	zvconj(x, y, xlen); \
	_X_ = y; \
} while (0)

/* op(<,ge>) * op(<,ge>) */
static
SEXP geMatrix_matmult(SEXP a, SEXP b, char atrans, char btrans)
{
	int *padim = DIM(a), am = padim[0], an = padim[1],
		rm = (atrans != 'N') ? an : am, rk = (atrans != 'N') ? am : an;

	if (b == R_NilValue) {

		if ((int_fast64_t) rm * rm > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");

		SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

		char rct = (TYPEOF(ax) != CPLXSXP || ((atrans != 'N') ? atrans : btrans) == 'C') ? 'C' : 'T';

		char rcl[] = "...Matrix";
		rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
		rcl[1] = (rct == 'C') ? 'p' : 's';
		rcl[2] = (rct == 'C') ? 'o' : 'y';

		SEXP r = PROTECT(newObject(rcl));

		SET_DIM(r, rm, rm);

		SEXP adimnames = PROTECT(DIMNAMES(a, 0)),
			rdimnames = PROTECT(DIMNAMES(r, 0));
		symDN(rdimnames, adimnames, (atrans != 'N') ? 1 : 0);
		UNPROTECT(2); /* rdimnames, adimnames */

		if (rct != 'C')
		SET_TRANS(r);

		if (rm > 0) {
		SEXP rx = PROTECT(Rf_allocVector(TYPEOF(ax), (R_xlen_t) rm * rm));
		if (TYPEOF(ax) == CPLXSXP) {
		Rcomplex *prx = COMPLEX(rx);
		memset(prx, 0, sizeof(Rcomplex) * (size_t) rm * (size_t) rm);
		if (rk > 0) {
		Rcomplex *pax = COMPLEX(ax),
			zero = Matrix_zzero, one = Matrix_zunit;
		if (rct == 'C')
		F77_CALL(zherk)("U", &atrans, &rm, &rk,
		                & one.r, pax, &am,
		                &zero.r, prx, &rm FCONE FCONE);
		else
		F77_CALL(zsyrk)("U", &atrans, &rm, &rk,
		                & one  , pax, &am,
		                &zero  , prx, &rm FCONE FCONE);
		}
		} else {
		double *prx = REAL(rx);
		memset(prx, 0, sizeof(double) * (size_t) rm * (size_t) rm);
		if (rk > 0) {
		double *pax = REAL(ax),
			zero = 0.0, one = 1.0;
		F77_CALL(dsyrk)("U", &atrans, &rm, &rk,
		                & one  , pax, &am,
		                &zero  , prx, &rm FCONE FCONE);
		}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
		}

		UNPROTECT(2); /* r, ax */
		return r;

	} else {

		int *pbdim = DIM(b), bm = pbdim[0], bn = pbdim[1],
			rn = (btrans != 'N') ? bm : bn;

		if (rk != ((btrans != 'N') ? bn : bm))
			Rf_error(_("non-conformable arguments"));
		if ((int_fast64_t) rm * rn > R_XLEN_T_MAX)
			Rf_error(_("attempt to allocate vector of length exceeding %s"),
			         "R_XLEN_T_MAX");

		SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

		char rcl[] = ".geMatrix";
		rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
		SEXP r = PROTECT(newObject(rcl));

		SET_DIM(r, rm, rn);

		SEXP adimnames = PROTECT(DIMNAMES(a, 0)),
			bdimnames = PROTECT(DIMNAMES(b, 0)),
			rdimnames = PROTECT(DIMNAMES(r, 0));
		matmultDN(rdimnames,
		          adimnames, (atrans != 'N') ? 1 : 0,
		          bdimnames, (btrans != 'N') ? 0 : 1);
		UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

		if (rm > 0 && rn > 0) {
		SEXP rx = PROTECT(Rf_allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
		if (TYPEOF(ax) == CPLXSXP) {
		Rcomplex *prx = COMPLEX(rx);
		if (rk == 0)
		memset(prx, 0, sizeof(Rcomplex) * (size_t) rm * (size_t) rn);
		else {
		SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
		Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx),
			zero = Matrix_zzero, one = Matrix_zunit;
		F77_CALL(zgemm)(&atrans, &btrans, &rm, &rn, &rk,
		                & one, pax, &am, pbx, &bm,
		                &zero, prx, &rm FCONE FCONE);
		UNPROTECT(1); /* bx */
		}
		} else {
		double *prx = REAL(rx);
		if (rk == 0)
		memset(prx, 0, sizeof(double) * (size_t) rm * (size_t) rn);
		else {
		SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
		double *pax = REAL(ax), *pbx = REAL(bx),
			zero = 0.0, one = 1.0;
		F77_CALL(dgemm)(&atrans, &btrans, &rm, &rn, &rk,
		                & one, pax, &am, pbx, &bm,
		                &zero, prx, &rm FCONE FCONE);
		UNPROTECT(1); /* bx */
		}
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
		}

		UNPROTECT(2); /* r, ax */
		return r;

	}
}

/* op(<,sy>) * op(<,ge>)  or  op(<,ge>) * op(<,sy>) */
static
SEXP syMatrix_matmult(SEXP a, SEXP b, char atrans, char btrans, char aside)
{
	int *padim = DIM(a),
		rk = padim[0];

	int *pbdim = DIM(b), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		Rf_error(_("non-conformable arguments"));
	if ((int_fast64_t) rm * rn > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, rm, rn);

	SEXP adimnames = PROTECT(DIMNAMES(a, -1)),
		bdimnames = PROTECT(DIMNAMES(b, 0)),
		rdimnames = PROTECT(DIMNAMES(r, 0));
	if (aside == 'L')
	matmultDN(rdimnames, adimnames,             0, bdimnames, btrans == 'N');
	else
	matmultDN(rdimnames, bdimnames, btrans != 'N', adimnames,             1);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	/* use *mm                     */
	/*   R := A  B                 */
	/*   R := A' B                 */
	/*   R := B  A                 */
	/*   R := B  A'                */
	/* use *mv and access B by row */
	/*   R := A  B.                */
	/*   R := A' B.                */
	/* use *mv and access R by row */
	/*   R := B. A  = (A.  B).     */
	/*   R := B. A' = (A'. B).     */

	if (rm > 0 && rn > 0) {
	SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(Rf_allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
	char aul = UPLO(a);
	int i,
		d     = (aside == 'L') ? rn : rm,
		binc  = (aside == 'L') ? bm :  1,
		bincp = (aside == 'L') ?  1 : bm,
		rinc  = (aside == 'L') ?  1 : rm,
		rincp = (aside == 'L') ? rm :  1;

	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		zero = Matrix_zzero, one = Matrix_zunit;
	char act = TRANS(a);
	if (btrans == 'N') {
	if (atrans != 'N' && atrans != act)
		CONJ2(pax, rk, rk);
	if (act == 'C')
	F77_CALL(zhemm)(&aside, &aul, &rm, &rn,
	                & one, pax, &rk, pbx, &bm,
	                &zero, prx, &rm FCONE FCONE);
	else
	F77_CALL(zsymm)(&aside, &aul, &rm, &rn,
	                & one, pax, &rk, pbx, &bm,
	                &zero, prx, &rm FCONE FCONE);
	} else {
	if (aside == 'L') {
	if (atrans != 'N' && atrans != act)
		CONJ2(pax, rk, rk);
	if (btrans == 'C')
		CONJ2(pbx, bm, bn);
	} else {
	if (((atrans != 'N') ? atrans : act) != btrans)
		CONJ2(pax, rk, rk);
	}
	for (i = 0; i < d; ++i) {
	if (act == 'C')
	F77_CALL(zhemv)(       &aul, &rk,
	                & one, pax, &rk, pbx, &binc,
	                &zero, prx, &rinc FCONE);
	else
	F77_CALL(zsymv)(       &aul, &rk,
	                & one, pax, &rk, pbx, &binc,
	                &zero, prx, &rinc FCONE);
	pbx += bincp;
	prx += rincp;
	}
	if (aside != 'L' && btrans == 'C')
		zvconj(pax, NULL, (size_t) rk * (size_t) rk);
	}
	} else {
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		zero = 0.0, one = 1.0;
	if (btrans == 'N') {
	F77_CALL(dsymm)(&aside, &aul, &rm, &rn,
	                & one, pax, &rk, pbx, &bm,
	                &zero, prx, &rm FCONE FCONE);
	} else {
	for (i = 0; i < d; ++i) {
	F77_CALL(dsymv)(       &aul, &rk,
	                &one, pax, &rk, pbx, &binc,
	                &zero, prx, &rinc FCONE);
	pbx += bincp;
	prx += rincp;
	}
	}
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(2); /* rx, bx */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

/* op(<,sp>) * op(<,ge>)  or  op(<,ge>) * op(<,sp>) */
static
SEXP spMatrix_matmult(SEXP a, SEXP b, char atrans, char btrans, char aside)
{
	int *padim = DIM(a),
		rk = padim[0];

	int *pbdim = DIM(b), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		Rf_error(_("non-conformable arguments"));
	if ((int_fast64_t) rm * rn > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, rm, rn);

	SEXP adimnames = PROTECT(DIMNAMES(a, -1)),
		bdimnames = PROTECT(DIMNAMES(b, 0)),
		rdimnames = PROTECT(DIMNAMES(r, 0));
	if (aside == 'L')
	matmultDN(rdimnames, adimnames,             0, bdimnames, btrans == 'N');
	else
	matmultDN(rdimnames, bdimnames, btrans != 'N', adimnames,             1);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	/* use *mv                           */
	/*   R := A  B                       */
	/*   R := A' B                       */
	/* use *mv and access B       by row */
	/*   R := A  B.                      */
	/*   R := A' B.                      */
	/* use *mv and access       R by row */
	/*   R := B. A  = (A.  B ).          */
	/*   R := B. A' = (A'. B ).          */
	/* use *mv and access B and R by row */
	/*   R := B  A  = (A*  B*)*          */
	/*   R := B  A' = (A   B')'          */

	if (rm > 0 && rn > 0) {
	SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) rm * rn));
	char aul = UPLO(a);
	int i,
		d     = ( aside == 'L'                    ) ? rn : rm,
		binc  = ((aside == 'L') == (btrans != 'N')) ? bm :  1,
		bincp = ((aside == 'L') == (btrans != 'N')) ?  1 : bm,
		rinc  = ( aside == 'L'                    ) ?  1 : rm,
		rincp = ( aside == 'L'                    ) ? rm :  1;

	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		zero = Matrix_zzero, one = Matrix_zunit;
	char act = TRANS(a);
	if (aside == 'L') {
	if (atrans != 'N' && atrans != act)
		CONJ1(pax, rk);
	if (btrans == 'C')
		CONJ2(pbx, bm, bn);
	} else {
	if (btrans != 'N' && ((atrans != 'N') ? atrans : act) != btrans)
		CONJ1(pax, rk);
	if (btrans == 'N' && ((atrans != 'N') ? atrans : act) == 'C')
		CONJ2(pbx, bm, bn);
	}
	for (i = 0; i < d; ++i) {
	if (act == 'C')
	F77_CALL(zhpmv)(&aul, &rk,
	                & one, pax, pbx, &binc,
	                &zero, prx, &rinc FCONE);
	else
	F77_CALL(zspmv)(&aul, &rk,
	                & one, pax, pbx, &binc,
	                &zero, prx, &rinc FCONE);
	pbx += bincp;
	prx += rincp;
	}
	if (aside != 'L' &&
	    ((btrans != 'N') ? btrans : ((atrans != 'N') ? atrans : act)) == 'C')
		zvconj(prx, NULL, (size_t) rm * (size_t) rn);
	} else {
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		zero = 0.0, one = 1.0;
	for (i = 0; i < d; ++i) {
	F77_CALL(dspmv)(&aul, &rk,
	                & one, pax, pbx, &binc,
	                &zero, prx, &rinc FCONE);
	pbx += bincp;
	prx += rincp;
	}
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(2); /* rx, bx */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

/* op(<,tr>) * op(<,ge>)  or  op(<,ge>) * op(<,tr>) */
static
SEXP trMatrix_matmult(SEXP a, SEXP b, char atrans, char btrans, char aside,
                      int triangular)
{
	int *padim = DIM(a),
		rk = padim[0];

	int *pbdim = DIM(b), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		Rf_error(_("non-conformable arguments"));
	if ((int_fast64_t) rm * rn > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	rcl[1] = (triangular != 0) ? 't' : 'g';
	rcl[2] = (triangular != 0) ? 'r' : 'e';
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, rm, rn);

	SEXP adimnames = PROTECT(DIMNAMES(a, 0)),
		bdimnames = PROTECT(DIMNAMES(b, 0)),
		rdimnames = PROTECT(DIMNAMES(r, 0));
	if (aside == 'L')
	matmultDN(rdimnames, adimnames, atrans != 'N', bdimnames, btrans == 'N');
	else
	matmultDN(rdimnames, bdimnames, btrans != 'N', adimnames, atrans == 'N');
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	char aul = UPLO(a);
	if (triangular < 0)
		SET_UPLO(r);

	char anu = DIAG(a);
	if (triangular < -1 || triangular > 1)
		SET_DIAG(r);

	/* use *mm after     copying B into R */
	/*   R := A  B                        */
	/*   R := A' B                        */
	/*   R := B  A                        */
	/*   R := B  A'                       */
	/* use *mm after transposing B into R */
	/*   R := A  B.                       */
	/*   R := A' B.                       */
	/*   R := B. A                        */
	/*   R := B. A'                       */

	if (rm > 0 && rn > 0) {
	SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(Rf_allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		one = Matrix_zunit;
	ztrans2(prx, pbx, (size_t) bm, (size_t) bn, btrans);
	F77_CALL(ztrmm)(&aside, &aul, &atrans, &anu, &rm, &rn,
	                &one, pax, &rk, prx, &rm FCONE FCONE FCONE FCONE);
	} else {
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		one = 1.0;
	dtrans2(prx, pbx, (size_t) bm, (size_t) bn, btrans);
	F77_CALL(dtrmm)(&aside, &aul, &atrans, &anu, &rm, &rn,
	                &one, pax, &rk, prx, &rm FCONE FCONE FCONE FCONE);
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(2); /* rx, bx */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

/* op(<,tp>) * op(<,ge>)  or  op(<,ge>) * op(<,tp>) */
static
SEXP tpMatrix_matmult(SEXP a, SEXP b, char atrans, char btrans, char aside,
                      int triangular)
{
	int *padim = DIM(a),
		rk = padim[0];

	int *pbdim = DIM(b), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		Rf_error(_("non-conformable arguments"));
	if ((int_fast64_t) rm * rn > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	rcl[1] = (triangular != 0) ? 't' : 'g';
	rcl[2] = (triangular != 0) ? 'r' : 'e';
	SEXP r = PROTECT(newObject(rcl));

	SET_DIM(r, rm, rn);

	SEXP adimnames = PROTECT(DIMNAMES(a, 0)),
		bdimnames = PROTECT(DIMNAMES(b, 0)),
		rdimnames = PROTECT(DIMNAMES(r, 0));
	if (aside == 'L')
	matmultDN(rdimnames, adimnames, atrans != 'N', bdimnames, btrans == 'N');
	else
	matmultDN(rdimnames, bdimnames, btrans != 'N', adimnames, atrans == 'N');
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	char aul = UPLO(a);
	if (triangular < 0)
		SET_UPLO(r);

	char anu = DIAG(a);
	if (triangular < -1 || triangular > 1)
		SET_DIAG(r);

	/* use *mv after     copying B into R                     */
	/*   R := A  B                                            */
	/*   R := A' B                                            */
	/* use *mv after transposing B into R                     */
	/*   R := A  B.                                           */
	/*   R := A' B.                                           */
	/* use *mv after transposing B into R and access R by row */
	/*   R := B. A  = (A.  B ).                               */
	/*   R := B. A' = (A'. B ).                               */
	/* use *mv after     copying B into R and access R by row */
	/*   R := B  A  = (A*  B*)*                               */
	/*   R := B  A' = (A   B')'                               */

	if (rm > 0 && rn > 0) {
	SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) rm * rn));
	int i, rinc = (aside == 'L') ? 1 : rm, rincp = (aside == 'L') ? rm : 1;

	char
		atransp = (aside == 'L')
		? atrans : ((atrans != 'N') ? 'N' : ((btrans != 'N') ? btrans : 'T')),
		btransp = (aside == 'L')
		? btrans : ((btrans != 'N') ? 'T' : 'N');

	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx);
	if (aside != 'L' && atrans != 'N' && btrans != 'N' && atrans != btrans)
		CONJ1(pax, rk);
	ztrans2(prx, pbx, (size_t) bm, (size_t) bn, btransp);
	if (aside != 'L' && atrans == 'C' && btrans == 'N')
		zvconj(prx, NULL, (size_t) rm * (size_t) rn);
	for (i = 0; i < rn; ++i) {
	F77_CALL(ztpmv)(&aul, &atransp, &anu, &rk,
	                pax, prx, &rinc FCONE FCONE FCONE);
	prx += rincp;
	}
	if (aside != 'L' &&
	    ((btrans != 'N') ? btrans : ((atrans != 'N') ? atrans : 'T')) == 'C')
		zvconj(prx, NULL, (size_t) rm * (size_t) rn);
	} else {
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx);
	dtrans2(prx, pbx, (size_t) bm, (size_t) bn, btransp);
	for (i = 0; i < rn; ++i) {
	F77_CALL(dtpmv)(&aul, &atransp, &anu, &rk,
	                pax, prx, &rinc FCONE FCONE FCONE);
	prx += rincp;
	}
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(2); /* rx, bx */
	}

	UNPROTECT(2); /* r, ax */
	return r;
}

SEXP R_dense_matmult(SEXP s_x, SEXP s_y,
                     SEXP s_xtrans, SEXP s_ytrans)
{
	char
		xtrans = CHAR(STRING_ELT(s_xtrans, 0))[0],
		ytrans = CHAR(STRING_ELT(s_ytrans, 0))[0],
		ztrans = 'N';
	int m, n, v;
	matmultDim(s_x, s_y, &xtrans, &ytrans, &ztrans, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(s_x, &xpid);
	PROTECT_WITH_INDEX(s_y, &ypid);

#define DO_S3(s_a, atrans, apid, isv) \
	do { \
	if (TYPEOF(s_a) != OBJSXP) { \
	REPROTECT(s_a = matrix_as_dense(s_a, ",ge", '\0', '\0', '\0', (atrans != 'N') ? 0 : 1, 0), apid); \
	if (isv) { \
		/* Vector: discard names and don't transpose again */ \
		SET_VECTOR_ELT(DIMNAMES(s_a, 0), \
		               (atrans != 'N') ? 1 : 0, R_NilValue); \
		atrans = 'N'; \
	} \
	} \
	} while (0)

	DO_S3(s_x, xtrans, xpid, v % 2);
	if (s_y != R_NilValue)
	DO_S3(s_y, ytrans, ypid, v > 1);

#undef DO_S3

	const char **valid = valid_matmult + 56;
	const char *xclass = NULL, *yclass = NULL;
	xclass = Matrix_class(s_x, valid, 6, __func__);
	if (s_y != R_NilValue)
	yclass = Matrix_class(s_y, valid, 6, __func__);

	char kind = (xclass[0] == 'z' || (s_y != R_NilValue && yclass[0] == 'z'))
		? 'z' : 'd';

#define DO_AS(s_a, aclass, atrans, apid) \
	do { \
	if (aclass[0] != kind) { \
		REPROTECT(s_a = dense_as_kind(s_a, aclass, kind, 0), apid); \
		aclass = Matrix_class(s_a, valid, 6, __func__); \
	} \
	} while (0)

	DO_AS(s_x, xclass, xtrans, xpid);
	if (s_y != R_NilValue)
	DO_AS(s_y, yclass, ytrans, ypid);

#undef DO_AS

	if (s_y == R_NilValue) {
		REPROTECT(s_x = dense_as_general(s_x, xclass, 1), xpid);
		s_x = geMatrix_matmult(s_x, s_y, xtrans, ytrans);
	} else if (xclass[1] == 'g' && yclass[1] == 'g') {
		s_x = geMatrix_matmult(s_x, s_y, xtrans, ytrans);
	} else if (xclass[1] == 'g' || yclass[1] == 'g') {
		s_x = (xclass[1] == 'g')
			? ((yclass[1] == 's')
			   ? ((yclass[2] != 'p')
			      ? syMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R')
			      : spMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R'))
			   : ((yclass[2] != 'p')
			      ? trMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R', 0)
			      : tpMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R', 0)))
			: ((xclass[1] == 's')
			   ? ((xclass[2] != 'p')
			      ? syMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L')
			      : spMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L'))
			   : ((xclass[2] != 'p')
			      ? trMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L', 0)
			      : tpMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L', 0)));
	} else if (xclass[1] == 's' && yclass[1] == 's') {
		if (xclass[2] == 'p' && yclass[2] == 'p') {
			REPROTECT(s_y = dense_as_general(s_y, yclass, 1), ypid);
			s_x = spMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L');
		} else if (xclass[2] == 'p') {
			REPROTECT(s_x = dense_as_general(s_x, xclass, 1), xpid);
			s_x = syMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R');
		} else {
			REPROTECT(s_y = dense_as_general(s_y, yclass, 1), ypid);
			s_x = syMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L');
		}
	} else if (xclass[1] == 's' || yclass[1] == 's') {
		if (xclass[1] == 's') {
			REPROTECT(s_x = dense_as_general(s_x, xclass, 1), xpid);
			s_x = (yclass[2] != 'p')
				? trMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R', 0)
				: tpMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R', 0);
		} else {
			REPROTECT(s_y = dense_as_general(s_y, yclass, 1), ypid);
			s_x = (xclass[2] != 'p')
				? trMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L', 0)
				: tpMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L', 0);
		}
	} else {
		int triangular = 0;

#define DO_TR \
		do { \
		char xul = UPLO(s_x), xnu = DIAG(s_x), \
			yul = UPLO(s_y), ynu = DIAG(s_y); \
		if (xtrans != 'N') \
			xul = (xul == 'U') ? 'L' : 'U'; \
		if (ytrans != 'N') \
			yul = (yul == 'U') ? 'L' : 'U'; \
		triangular = (xul != yul) ? 0 : ((xnu != ynu || xnu == 'N') ? 1 : 2); \
		if (xul != 'U') \
			triangular = -triangular; \
		} while (0)

		DO_TR;

		if (xclass[2] == 'p' && yclass[2] == 'p') {
			REPROTECT(s_y = dense_as_general(s_y, yclass, 1), ypid);
			s_x = tpMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L', triangular);
		} else if (xclass[2] == 'p') {
			REPROTECT(s_x = dense_as_general(s_x, xclass, 1), xpid);
			s_x = trMatrix_matmult(s_y, s_x, ytrans, xtrans, 'R', triangular);
		} else {
			REPROTECT(s_y = dense_as_general(s_y, yclass, 1), ypid);
			s_x = trMatrix_matmult(s_x, s_y, xtrans, ytrans, 'L', triangular);
		}
	}

	UNPROTECT(2); /* s_y, s_x */
	return s_x;
}

/* boolean: op(op(<.gC>) & op(<.gC>)) */
/* numeric: op(op(<,gC>) * op(<,gC>)) */
static
SEXP gCgCMatrix_matmult(SEXP x, SEXP y, char xtrans, char ytrans, char ztrans,
                        int triangular, int boolean)
{

#define ASMODE(_TRANS_) ((boolean) ? 0 : (((_TRANS_) != 'C') ? 1 : 2))

	cholmod_sparse *X = M2CHS(x, !boolean), *Y = NULL, *Z = NULL;

	SEXP z;
	char zclass[] = "..CMatrix";

	if (y == R_NilValue) {

		zclass[0] = (boolean) ? 'n' : ((X->xtype != CHOLMOD_COMPLEX) ? 'd' : 'z');
		zclass[1] = (boolean) ? 's' : ((X->xtype != CHOLMOD_COMPLEX) ? 'p' : ((((xtrans != 'N') ? xtrans : ytrans) == 'C') ? 'p' : 's'));

		if (xtrans != 'N')
			X = cholmod_transpose(X, ASMODE(xtrans), &c);
		Z = cholmod_aat(X, NULL, 0, ASMODE((xtrans != 'N') ? xtrans : ytrans), &c);
		if (xtrans != 'N')
			cholmod_free_sparse(&X, &c);
		if (!Z->sorted)
			cholmod_sort(Z, &c);
		X = cholmod_copy(Z, (ztrans != 'N') ? -1 : 1, !boolean, &c);
		cholmod_free_sparse(&Z, &c);
		Z = X;
		PROTECT(z = CHS2M(Z, !boolean, zclass[1]));
		cholmod_free_sparse(&Z, &c);

		SEXP xdimnames = PROTECT(DIMNAMES(x, 0)),
			zdimnames = PROTECT(DIMNAMES(z, 0));
		symDN(zdimnames, xdimnames, (xtrans != 'N') ? 1 : 0);
		UNPROTECT(2); /* zdimnames, xdimnames */

		if (ztrans != 'N')
			SET_UPLO(z);
		if (zclass[1] == 's' && zclass[0] == 'z')
			SET_TRANS(z);

	} else {

		Y = M2CHS(y, !boolean);

		if (((xtrans != 'N') ? X->nrow : X->ncol) !=
		    ((ytrans != 'N') ? Y->ncol : Y->nrow))
			Rf_error(_("non-conformable arguments"));

		zclass[0] = (boolean) ? 'n' : ((X->xtype != CHOLMOD_COMPLEX && Y->xtype != CHOLMOD_COMPLEX) ? 'd' : 'z');
		zclass[1] = (triangular != 0) ? 't' : 'g';

		if (xtrans != 'N')
			X = cholmod_transpose(X, ASMODE(xtrans), &c);
		if (ytrans != 'N')
			Y = cholmod_transpose(Y, ASMODE(ytrans), &c);
		Z = cholmod_ssmult(X, Y, 0, !boolean, 1, &c);
		if (xtrans != 'N')
			cholmod_free_sparse(&X, &c);
		if (ytrans != 'N')
			cholmod_free_sparse(&Y, &c);
		if (triangular < -1)
			cholmod_band_inplace(-((int) Z->nrow), -1, !boolean, Z, &c);
		else if (triangular > 1)
			cholmod_band_inplace(1, (int) Z->ncol, !boolean, Z, &c);
		PROTECT(z = CHS2M(Z, !boolean, zclass[1]));
		cholmod_free_sparse(&Z, &c);

		SEXP xdimnames = PROTECT(DIMNAMES(x, 0)),
			ydimnames = PROTECT(DIMNAMES(y, 0)),
			zdimnames = PROTECT(DIMNAMES(z, 0));
		matmultDN(zdimnames,
		          xdimnames, (xtrans != 'N') ? 1 : 0,
		          ydimnames, (ytrans != 'N') ? 0 : 1);
		UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

		if (triangular < 0)
			SET_UPLO(z);
		if (triangular < -1 || triangular > 1)
			SET_DIAG(z);

	}

	if (ztrans != 'N')
		z = sparse_transpose(z, zclass, ztrans, 1);

	UNPROTECT(1); /* z */
	return z;
}

/* op(op(<,[gs]C>) * op(<,ge>)) */
/* NB: currently, caller must coerce complex symmetric zsC to zgC */
static
SEXP gCgeMatrix_matmult(SEXP x, SEXP y, int xtrans, char ytrans, char ztrans,
                        int triangular, int symmetric)
{
	cholmod_sparse *X = M2CHS(x, 1);
	X->stype = symmetric;

	cholmod_dense *Y = M2CHD(y, ytrans);

	if (((xtrans != 'N') ? X->nrow : X->ncol) != Y->nrow)
		Rf_error(_("non-conformable arguments"));
	size_t m = (xtrans != 'N') ? X->ncol : X->nrow, n = Y->ncol;
	if (m * n > R_XLEN_T_MAX)
		Rf_error(_("attempt to allocate vector of length exceeding %s"),
		         "R_XLEN_T_MAX");

	if (X->xtype == CHOLMOD_COMPLEX && xtrans == 'T')
		CONJ2(X->x, X->nrow, X->ncol);

	char zclass[] = "...Matrix";
	zclass[0] = (X->xtype != CHOLMOD_COMPLEX && Y->xtype != CHOLMOD_COMPLEX) ? 'd' : 'z';
	zclass[1] = (triangular != 0) ? 't' : 'g';
	zclass[2] = (triangular != 0) ? 'r' : 'e';
	SEXP z = PROTECT(newObject(zclass));

	SET_DIM(z,
	        (int) ((ztrans != 'N') ? n : m),
	        (int) ((ztrans != 'N') ? m : n));

	SEXP xdimnames = PROTECT(DIMNAMES(x, -(symmetric != 0))),
		ydimnames = PROTECT(DIMNAMES(y, 0)),
		zdimnames = PROTECT(DIMNAMES(z, 0));
	if (ztrans != 'N')
	matmultDN(zdimnames,
	          ydimnames, (ytrans != 'N') ? 0 : 1,
	          xdimnames, (xtrans != 'N') ? 1 : 0);
	else
	matmultDN(zdimnames,
	          xdimnames, (xtrans != 'N') ? 1 : 0,
	          ydimnames, (ytrans != 'N') ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

	if (triangular != 0 && (ztrans != 'N') == (triangular > 0))
		SET_UPLO(z);
	if (triangular < -1 || triangular > 1)
		SET_DIAG(z);

	double alpha[2] = { 1.0, 0.0 }, beta[2] = { 0.0, 0.0 };
	cholmod_dense *Z = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(Z, 0, sizeof(cholmod_dense));
	Z->nrow = m;
	Z->ncol = n;
	Z->d = Z->nrow;
	Z->nzmax = Z->nrow * Z->ncol;
	Z->xtype = X->xtype;
	Z->dtype = X->dtype;

	SEXP zx;
	if (zclass[0] == 'z') {
	PROTECT(zx = Rf_allocVector(CPLXSXP, (R_xlen_t) (m * n)));
	Z->x = (ztrans != 'N')
		? (Rcomplex *) R_alloc(m * n, sizeof(Rcomplex))
		: COMPLEX(zx);
	} else {
	PROTECT(zx = Rf_allocVector(REALSXP, (R_xlen_t) (m * n)));
	Z->x = (ztrans != 'N')
		? (double *) R_alloc(m * n, sizeof(double))
		: REAL(zx);
	}
	cholmod_sdmult(X, xtrans != 'N', alpha, beta, Y, Z, &c);
	if (ztrans != 'N') {
	if (zclass[0] == 'z')
	ztrans2(COMPLEX(zx), (Rcomplex *) Z->x, m, n, ztrans);
	else
	dtrans2(   REAL(zx), (  double *) Z->x, m, n, ztrans);
	}
	SET_SLOT(z, Matrix_xSym, zx);
	UNPROTECT(1); /* zx */

	UNPROTECT(1); /* z */
	return z;
}

SEXP R_sparse_matmult(SEXP s_x, SEXP s_y,
                      SEXP s_xtrans, SEXP s_ytrans, SEXP s_ztrans,
                      SEXP s_boolean)
{
	if (TYPEOF(s_boolean) != LGLSXP || LENGTH(s_boolean) < 1)
		Rf_error(_("invalid '%s' to '%s'"), "boolean", __func__);
	int boolean = LOGICAL(s_boolean)[0];

	char
		xtrans = CHAR(STRING_ELT(s_xtrans, 0))[0],
		ytrans = CHAR(STRING_ELT(s_ytrans, 0))[0],
		ztrans = CHAR(STRING_ELT(s_ztrans, 0))[0];
	int m, n, v;
	matmultDim(s_x, s_y, &xtrans, &ytrans, &ztrans, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(s_x, &xpid);
	PROTECT_WITH_INDEX(s_y, &ypid);

#define DO_S3(s_a, atrans, apid, isv) \
	do { \
	if (TYPEOF(s_a) != OBJSXP) { \
	if (boolean == NA_LOGICAL || !boolean) \
	REPROTECT(s_a = matrix_as_dense (s_a, ",ge", '\0', '\0', '\0', (atrans != 'N') ? 0 : 1, 0), apid); \
	else if (atrans != 'N') \
	REPROTECT(s_a = matrix_as_sparse(s_a, "ngR", '\0', '\0', '\0', 0), apid); \
	else \
	REPROTECT(s_a = matrix_as_sparse(s_a, "ngC", '\0', '\0', '\0', 0), apid); \
	if (isv) { \
		/* Discard names and don't transpose again */ \
		SET_VECTOR_ELT(DIMNAMES(s_a, 0), \
		               (atrans != 'N') ? 1 : 0, R_NilValue); \
		atrans = 'N'; \
	} \
	} \
	} while (0)

	DO_S3(s_x, xtrans, xpid, v % 2);
	if (s_y != R_NilValue)
	DO_S3(s_y, ytrans, ypid, v > 1);

#undef DO_S3

	const char **valid = valid_matmult + 5;
	const char *xclass = NULL, *yclass = NULL;
	xclass = Matrix_class(s_x, valid, 6, __func__);
	if (s_y != R_NilValue)
	yclass = Matrix_class(s_y, valid, 6, __func__);

	if (boolean == NA_LOGICAL)
		boolean = xclass[0] == 'n' && (s_y == R_NilValue || yclass[0] == 'n');
	char kind = (boolean) ? 'n' :
		((xclass[0] != 'z' && (s_y == R_NilValue || yclass[0] != 'z')) ? 'd' : 'z');

#define DO_AS(s_a, s_b, aclass, atrans, btrans, apid) \
	do { \
	if (aclass[2] != 'C' && atrans != 'N') { \
		if (aclass[2] != 'R' && aclass[2] != 'T') { \
			REPROTECT(s_a = dense_as_sparse(s_a, aclass, 'R'), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
		REPROTECT(s_a = sparse_transpose(s_a, aclass, atrans, 1), apid); \
		aclass = Matrix_class(s_a, valid, 6, __func__); \
		if (s_b == R_NilValue) \
			SWAP(atrans, btrans, char, ); \
		else \
			atrans = 'N'; \
	} \
	if (aclass[2] != 'C') { \
		if (aclass[2] != 'R' && aclass[2] != 'T') \
			REPROTECT(s_a = dense_as_sparse(s_a, aclass, 'C'), apid); \
		else \
			REPROTECT(s_a = sparse_as_Csparse(s_a, aclass), apid); \
		aclass = Matrix_class(s_a, valid, 6, __func__); \
	} \
	if (atrans != 'N' && aclass[1] == 's' && \
	    (aclass[0] != 'z' || TRANS(s_a) == atrans)) \
		atrans = 'N'; \
	if (aclass[0] != kind) { \
		if (boolean) \
			REPROTECT(s_a = sparse_dropzero(s_a, aclass, 0.0), apid); \
		else { \
			REPROTECT(s_a = sparse_as_kind(s_a, aclass, kind), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
	} \
	} while (0)

	DO_AS(s_x, s_y, xclass, xtrans, ytrans, xpid);

	if (s_y == R_NilValue) {
		REPROTECT(s_x = sparse_as_general(s_x, xclass), xpid);
		s_x = gCgCMatrix_matmult(
			s_x, s_y, xtrans, ytrans, ztrans, 0, boolean);
		UNPROTECT(2); /* s_y, s_x */
		return s_x;
	}

	int triangular = 0;
	if (xclass[1] == 't' && yclass[1] == 't')
		DO_TR;

	if (!boolean && yclass[2] != 'C' && yclass[2] != 'R' && yclass[2] != 'T') {
		if (xclass[1] == 's' && xclass[0] == 'z' && TRANS(s_x) != 'C') {
			REPROTECT(s_x = sparse_as_general(s_x, xclass), xpid);
			xclass = Matrix_class(s_x, valid, 6, __func__);
		}
		int symmetric = xclass[1] == 's';
		if (symmetric && UPLO(s_x) != 'U')
			symmetric = -symmetric;
		if (yclass[0] != kind) {
			REPROTECT(s_y = dense_as_kind(s_y, yclass, kind, 0), ypid);
			yclass = Matrix_class(s_y, valid, 6, __func__);
		}
		REPROTECT(s_x = sparse_diag_U2N(s_x, xclass), xpid);
		REPROTECT(s_y = dense_as_general(s_y, yclass, 1), ypid);
		s_x = gCgeMatrix_matmult(
			s_x, s_y, xtrans, ytrans, ztrans, triangular, symmetric);
		UNPROTECT(2); /* s_y, s_x */
		return s_x;
	}

	DO_AS(s_y, s_x, yclass, ytrans, xtrans, ypid);

#undef DO_AS

	REPROTECT(s_x = sparse_as_general(s_x, xclass), xpid);
	REPROTECT(s_y = sparse_as_general(s_y, yclass), ypid);
	s_x = gCgCMatrix_matmult(
		s_x, s_y, xtrans, ytrans, ztrans, triangular, boolean);
	UNPROTECT(2); /* s_y, s_x */
	return s_x;
}

#define zSCALE(x, y) \
	do { \
		Rcomplex tmp = (x); \
		(x).r = tmp.r * (y).r - tmp.i * (y).i; \
		(x).i = tmp.r * (y).i + tmp.i * (y).r; \
	} while (0)
#define dSCALE(x, y) \
	(x) = (x)  * (y)
#define lSCALE(x, y) \
	(x) = (x) && (y)

#define SCALE3(t, index) \
	do { \
		switch (t) { \
		case CPLXSXP: SCALE(z, index); break; \
		case REALSXP: SCALE(d, index); break; \
		case  LGLSXP: SCALE(l, index); break; \
		default: break; \
		} \
	} while (0)

static
void dense_colscale(SEXP obj, SEXP d, int m, int n, char ul, char nu)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = XLENGTH(x) != (int_fast64_t) m * n;

#define SCALE(c, index) \
	do { \
		c##TYPE *px = c##PTR(x), *pd = c##PTR(d); \
		if (ul == '\0') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < m; ++i) { \
					c##SCALE(*px, pd[index]); \
					++px; \
				} \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i <= j; ++i) { \
					c##SCALE(*px, pd[index]); \
					++px; \
				} \
				if (!packed) \
					px += n - j - 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				for (i = j; i < n; ++i) { \
					c##SCALE(*px, pd[index]); \
					++px; \
				} \
			} \
		} \
		if (nu != '\0' && nu != 'N') { \
			size_t n_ = (size_t) n; \
			if (!packed) \
				c##NAME(copy2)(n_, c##PTR(x), n_ + 1, pd, 1); \
			else if (ul == 'U') \
				c##NAME(copy1)(n_, c##PTR(x), 2 , 1, 0, pd, 1, 0, 0); \
			else \
				c##NAME(copy1)(n_, c##PTR(x), n_, 1, 1, pd, 1, 0, 0); \
		} \
	} while (0)

	SCALE3(TYPEOF(d), j);
	return;
}

static
void dense_rowscale(SEXP obj, SEXP d, int m, int n, char ul, char nu)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = XLENGTH(x) != (int_fast64_t) m * n;
	SCALE3(TYPEOF(d), i);

#undef SCALE

	return;
}

/* boolean: <lgC> & <ldi>  or  <ldi> & <lgR> */
/* numeric: <,gC> * <,di>  or  <,di> * <,gR> */
static
void Csparse_colscale(SEXP obj, SEXP d)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym));
	int *pp = INTEGER(p) + 1, n = (int) (XLENGTH(p) - 1), j, k, kend;
	UNPROTECT(2); /* p, x */

#define SCALE(c, index) \
	do { \
		c##TYPE *px = c##PTR(x), *pd = c##PTR(d); \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				c##SCALE(*px, *pd); \
				++px; \
				++k; \
			} \
			++pd; \
		} \
	} while (0)

	SCALE3(TYPEOF(d), );

#undef SCALE

	return;
}

/* boolean: <ldi> & <lgC>  or  <lgR> & <ldi> */
/* numeric: <,di> * <,gC>  or  <,gR> * <,di> */
static
void Csparse_rowscale(SEXP obj, SEXP d, SEXP iSym)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		p = PROTECT(GET_SLOT(obj, Matrix_pSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	int *pi = INTEGER(i), k, kend = INTEGER(p)[XLENGTH(p) - 1];
	UNPROTECT(3); /* i, p, x */

#define SCALE(c, index) \
	do { \
		c##TYPE *px = c##PTR(x), *pd = c##PTR(d); \
		for (k = 0; k < kend; ++k) { \
			c##SCALE(*px, pd[*pi]); \
			++px; \
			++pi; \
		} \
	} while (0)

	SCALE3(TYPEOF(d), );
	return;
}

/* boolean: <ldi> & <lgT>  or  <lgT> & <ldi> */
/* numeric: <,di> * <,gT>  or  <,gT> * <,di> */
static
void Tsparse_rowscale(SEXP obj, SEXP d, SEXP iSym)
{
	SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym)),
		i = PROTECT(GET_SLOT(obj, iSym));
	int *pi = INTEGER(i);
	R_xlen_t k, kend = XLENGTH(i);
	UNPROTECT(2); /* i, x */
	SCALE3(TYPEOF(d), );

#undef SCALE

	return;
}

#undef SCALE3

SEXP R_diagonal_matmult(SEXP s_x, SEXP s_y,
                        SEXP s_xtrans, SEXP s_ytrans,
                        SEXP s_boolean)
{
	SEXP _s_x = s_x, _s_y = s_y; /* for later pointer comparison */

	if (TYPEOF(s_boolean) != LGLSXP || LENGTH(s_boolean) < 1)
		Rf_error(_("invalid '%s' to '%s'"), "boolean", __func__);
	int boolean = LOGICAL(s_boolean)[0];

	char
		xtrans = CHAR(STRING_ELT(s_xtrans, 0))[0],
		ytrans = CHAR(STRING_ELT(s_ytrans, 0))[0],
		ztrans = 'N';
	int m, n, v;
	matmultDim(s_x, s_y, &xtrans, &ytrans, &ztrans, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(s_x, &xpid);
	PROTECT_WITH_INDEX(s_y, &ypid);

#define DO_S3(s_a, atrans, apid, isv) \
	do { \
	if (TYPEOF(s_a) != OBJSXP) { \
	if (boolean == NA_LOGICAL || !boolean) \
	REPROTECT(s_a = matrix_as_dense(s_a, ",ge", '\0', '\0', '\0', (atrans != 'N') ? 0 : 1, 2), apid); \
	else \
	REPROTECT(s_a = matrix_as_dense(s_a, "nge", '\0', '\0', '\0', (atrans != 'N') ? 0 : 1, 2), apid); \
	if (isv) { \
		/* Vector: discard names and don't transpose again */ \
		SET_VECTOR_ELT(DIMNAMES(s_a, 0), \
		               (atrans != 'N') ? 1 : 0, R_NilValue); \
		atrans = 'N'; \
	} \
	} \
	} while (0)

	DO_S3(s_x, xtrans, xpid, v % 2);
	if (s_y != R_NilValue)
	DO_S3(s_y, ytrans, ypid, v > 1);

#undef DO_S3

	const char **valid = valid_matmult;
	const char *xclass = NULL, *yclass = NULL;
	xclass = Matrix_class(s_x, valid, 6, __func__);
	if (s_y != R_NilValue)
	yclass = Matrix_class(s_y, valid, 6, __func__);

	if (boolean == NA_LOGICAL)
		boolean = xclass[0] == 'n' && yclass[0] == 'n';
	char kind = (boolean) ? 'n' :
		((xclass[0] != 'z' && yclass[0] != 'z') ? 'd' : 'z');
	char ks = (boolean) ? 'l' : kind, kd = kind;

	int mg = -1, id = -1;
	if (xclass[2] == 'i') {
		mg = 0;
		id = DIAG(s_x) != 'N';
	} else if (yclass[2] == 'i') {
		mg = 1;
		id = DIAG(s_y) != 'N';
	} else
		Rf_error(_("should never happen ..."));

#define DO_AS(s_a, aclass, atrans, apid) \
	do { \
	if (atrans != 'N' && aclass[1] == 's' && \
	    (aclass[0] != 'z' || TRANS(s_a) == atrans)) \
		atrans = 'N'; \
	switch (aclass[2]) { \
	case 'i': \
		if (!id && aclass[0] != ks) { \
			REPROTECT(s_a = diagonal_as_kind(s_a, aclass, ks), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
		break; \
	case 'C': \
	case 'R': \
	case 'T': \
		if (aclass[0] != ks) { \
			REPROTECT(s_a = sparse_as_kind(s_a, aclass, ks), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
		if (!id && aclass[1] == 's') { \
			REPROTECT(s_a = sparse_as_general(s_a, aclass), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
		if (!id && aclass[1] == 't') \
			REPROTECT(s_a = sparse_diag_U2N(s_a, aclass), apid); \
		if (atrans != 'N') { \
			REPROTECT(s_a = sparse_transpose(s_a, aclass, atrans, 0), apid); \
			atrans = 'N'; \
		} \
		break; \
	default: \
		if (aclass[0] != kd) { \
			REPROTECT(s_a = dense_as_kind(s_a, aclass, kd, 1), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
		if (!id && aclass[1] == 's') { \
			REPROTECT(s_a = dense_as_general(s_a, aclass, s_a == _ ## s_a), apid); \
			aclass = Matrix_class(s_a, valid, 6, __func__); \
		} \
		if (atrans != 'N') { \
			REPROTECT(s_a = dense_transpose(s_a, aclass, atrans), apid); \
			atrans = 'N'; \
		} \
		break; \
	} \
	} while (0)

	DO_AS(s_x, xclass, xtrans, xpid);
	DO_AS(s_y, yclass, ytrans, ypid);

	PROTECT_INDEX zpid;
	const char *zclass = (mg == 0) ? yclass : xclass;
	SEXP z = newObject(zclass);
	PROTECT_WITH_INDEX(z, &zpid);

	SET_DIM(z, m, n);

	SEXP xdimnames = PROTECT(DIMNAMES(s_x, 0)),
		ydimnames = PROTECT(DIMNAMES(s_y, 0)),
		zdimnames = PROTECT(DIMNAMES(z, 0));
	matmultDN(zdimnames,
	          xdimnames, (xtrans != 'N') ? 1 : 0,
	          ydimnames, (ytrans != 'N') ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

	char ul = '\0', nu = '\0';
	if (zclass[1] != 'g' && (ul = UPLO((mg == 0) ? s_y : s_x)) != 'U')
		SET_UPLO(z);
	if (zclass[1] == 't' && (nu = DIAG((mg == 0) ? s_y : s_x)) != 'N' && id)
		SET_DIAG(z);

	if (zclass[2] == 'C' || zclass[2] == 'R' || zclass[2] == 'T') {
	if (zclass[2] != 'T') {
		SEXP p = PROTECT(GET_SLOT((mg == 0) ? s_y : s_x, Matrix_pSym));
		SET_SLOT(z, Matrix_pSym, p);
		UNPROTECT(1); /* p */
	}
	if (zclass[2] != 'R') {
		SEXP i = PROTECT(GET_SLOT((mg == 0) ? s_y : s_x, Matrix_iSym));
		SET_SLOT(z, Matrix_iSym, i);
		UNPROTECT(1); /* i */
	}
	if (zclass[2] != 'C') {
		SEXP j = PROTECT(GET_SLOT((mg == 0) ? s_y : s_x, Matrix_jSym));
		SET_SLOT(z, Matrix_jSym, j);
		UNPROTECT(1); /* j */
	}
	}

	SEXP x0 = PROTECT(GET_SLOT((mg == 0) ? s_y : s_x, Matrix_xSym));
	if (id || ((mg == 0) ? s_y != _s_y : s_x != _s_x))
	SET_SLOT(z, Matrix_xSym, x0);
	else {
	SEXP x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));
	switch (kind) {
	case 'z':
		memcpy(COMPLEX(x1), COMPLEX(x0), sizeof(Rcomplex) * (size_t) XLENGTH(x0));
		break;
	case 'd':
		memcpy(   REAL(x1),    REAL(x0), sizeof(  double) * (size_t) XLENGTH(x0));
		break;
	default:
		memcpy(LOGICAL(x1), LOGICAL(x0), sizeof(     int) * (size_t) XLENGTH(x0));
		break;
	}
	SET_SLOT(z, Matrix_xSym, x1);
	UNPROTECT(1); /* x1 */
	}
	UNPROTECT(1); /* x0 */

	if (!id) {
	SEXP d = PROTECT(GET_SLOT((mg == 0) ? s_x : s_y, Matrix_xSym));
	switch (zclass[2]) {
	case 'C':
		if (mg == 0)
			Csparse_rowscale(z, d, Matrix_iSym);
		else
			Csparse_colscale(z, d);
		break;
	case 'R':
		if (mg == 0)
			Csparse_colscale(z, d);
		else
			Csparse_rowscale(z, d, Matrix_jSym);
		break;
	case 'T':
		if (mg == 0)
			Tsparse_rowscale(z, d, Matrix_iSym);
		else
			Tsparse_rowscale(z, d, Matrix_jSym);
		break;
	default:
		if (mg == 0)
			dense_rowscale(z, d, m, n, ul, nu);
		else
			dense_colscale(z, d, m, n, ul, nu);
		break;
	}
	UNPROTECT(1); /* d */
	}

	if (boolean && (zclass[2] == 'C' || zclass[2] == 'R' || zclass[2] == 'T')) {
		REPROTECT(z = sparse_dropzero(z, zclass, 0.0), zpid);
		REPROTECT(z = sparse_as_kind(z, zclass, 'n'), zpid);
	}

	UNPROTECT(3); /* z, s_y, s_x */
	return z;
}
