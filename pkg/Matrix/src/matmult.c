#include "Lapack-etc.h"
#include "cholmod-etc.h"
#include "Mdefines.h"
#include "idz.h"
#include "coerce.h"
#include "dense.h"
#include "sparse.h"
#include "matmult.h"

static
void matmultDim(SEXP x, SEXP y, char *xtrans, char *ytrans, char *ztrans,
                int *m, int *n, int *v)
{
	int xt = *xtrans == 'C' || *xtrans == 'T'; if (!xt) *xtrans = 'N';
	int yt = *ytrans == 'C' || *ytrans == 'T'; if (!yt) *ytrans = 'N';
	int zt = *ztrans == 'C' || *ztrans == 'T'; if (!zt) *ztrans = 'N';
	if (y == R_NilValue) {
		if (xt == yt)
			error(_("should never happen"));
		SEXP
			xdim = (TYPEOF(x) == S4SXP)
			? GET_SLOT(x, Matrix_DimSym) : getAttrib(x, R_DimSymbol);
		if (TYPEOF(xdim) == INTSXP && LENGTH(xdim) == 2) {
			*v = 0;
			*m = *n = INTEGER(xdim)[(xt) ? 1 : 0];
		} else if (XLENGTH(x) <= INT_MAX) {
			*v = 1;
			*m = *n = (xt) ? 1 : LENGTH(x);
		} else
			error(_("dimensions cannot exceed %s"), "2^31-1");
	} else {
		/* MJ: So that I don't lose my mind ... : */
		if (zt) {
			int tmp = !xt; xt = !yt; yt = tmp;
			SEXP s = x; x = y; y = s;
		}
		SEXP
			xdim = (TYPEOF(x) == S4SXP)
			? GET_SLOT(x, Matrix_DimSym) : getAttrib(x, R_DimSymbol),
			ydim = (TYPEOF(y) == S4SXP)
			? GET_SLOT(y, Matrix_DimSym) : getAttrib(y, R_DimSymbol);
		int xm, xn, ym, yn, x2, y2;
		xm = xn = ym = yn = -1;
		x2 = TYPEOF(xdim) == INTSXP && LENGTH(xdim) == 2;
		y2 = TYPEOF(ydim) == INTSXP && LENGTH(ydim) == 2;
		if (x2) {
			int *pxdim = INTEGER(xdim);
			xm = pxdim[0];
			xn = pxdim[1];
		} else if (XLENGTH(x) > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
		if (y2) {
			int *pydim = INTEGER(ydim);
			ym = pydim[0];
			yn = pydim[1];
		} else if (XLENGTH(y) > INT_MAX)
			error(_("dimensions cannot exceed %s"), "2^31-1");
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
			error(_("non-conformable arguments"));
		*m = (xt) ? xn : xm;
		*n = (yt) ? ym : yn;
		if (zt) {
			int tmp = !xt; xt = !yt; yt = tmp;
			tmp = *m; *m = *n; *n = tmp;
		}
	}
	if (*v % 2) *xtrans = (xt) ? 'T' : 'N';
	if (*v > 1) *ytrans = (yt) ? 'T' : 'N';
	return;
}

static
void matmultDN(SEXP dest, SEXP asrc, int ai, SEXP bsrc, int bi) {
	SEXP s;
	if (!isNull(s = VECTOR_ELT(asrc, ai)))
		SET_VECTOR_ELT(dest, 0, s);
	if (!isNull(s = VECTOR_ELT(bsrc, bi)))
		SET_VECTOR_ELT(dest, 1, s);
	PROTECT(asrc = getAttrib(asrc, R_NamesSymbol));
	PROTECT(bsrc = getAttrib(bsrc, R_NamesSymbol));
	if (!isNull(asrc) || !isNull(bsrc)) {
		SEXP destnms = PROTECT(allocVector(STRSXP, 2));
		if (!isNull(asrc) && *CHAR(s = STRING_ELT(asrc, ai)) != '\0')
			SET_STRING_ELT(destnms, 0, s);
		if (!isNull(bsrc) && *CHAR(s = STRING_ELT(bsrc, bi)) != '\0')
			SET_STRING_ELT(destnms, 1, s);
		setAttrib(dest, R_NamesSymbol, destnms);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return;
}

#define CONJ2(_X_, _M_, _N_) \
do { \
	int m = (int) _M_, n = (int) _N_; \
	size_t mn = (size_t) m * n; \
	Rcomplex *dest = (Rcomplex *) R_alloc(mn, sizeof(Rcomplex)); \
	Rcomplex *src = (Rcomplex *) _X_; \
	for (size_t k = 0; k < mn; ++k) \
		ASSIGN2_COMPLEX_CJ(dest[k], src[k]); \
	_X_ = dest; \
} while (0)

#define CONJ1(_X_, _N_) \
do { \
	int n = (int) _N_; \
	size_t mn = (size_t) PACKED_LENGTH(n); \
	Rcomplex *dest = (Rcomplex *) R_alloc(mn, sizeof(Rcomplex)); \
	Rcomplex *src = (Rcomplex *) _X_; \
	for (size_t k = 0; k < mn; ++k) \
		ASSIGN2_COMPLEX_CJ(dest[k], src[k]); \
	_X_ = dest; \
} while (0)

/* op(<,ge>) * op(<,ge>) */
static
SEXP geMatrix_matmult(SEXP a, SEXP b, char atrans, char btrans)
{
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int *padim = INTEGER(adim), am = padim[0], an = padim[1],
		rm = (atrans != 'N') ? an : am, rk = (atrans != 'N') ? am : an;

	if (b == R_NilValue) {

		if ((Matrix_int_fast64_t) rm * rm > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");

		SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

		char rct = (TYPEOF(ax) != CPLXSXP || ((atrans != 'N') ? atrans : btrans) == 'C') ? 'C' : 'T';

		char rcl[] = "...Matrix";
		rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
		rcl[1] = (rct == 'C') ? 'p' : 's';
		rcl[2] = (rct == 'C') ? 'o' : 'y';

		SEXP r = PROTECT(newObject(rcl));

		SEXP rdim = GET_SLOT(r, Matrix_DimSym);
		int *prdim = INTEGER(rdim);
		prdim[0] = prdim[1] = rm;

		SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
			rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
		symDN(rdimnames, adimnames, (atrans != 'N') ? 1 : 0);
		UNPROTECT(2); /* rdimnames, adimnames */

		if (rct != 'C') {
		SEXP rtrans = PROTECT(mkString("T"));
		SET_SLOT(r, Matrix_transSym, rtrans);
		UNPROTECT(1); /* rtrans */
		}

		if (rm > 0) {
		SEXP rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rm));
		if (TYPEOF(ax) == CPLXSXP) {
		Rcomplex *prx = COMPLEX(rx);
		Matrix_memset(prx, 0, (R_xlen_t) rm * rm, sizeof(Rcomplex));
		if (rk > 0) {
		Rcomplex *pax = COMPLEX(ax),
			zero = Matrix_zzero, one = Matrix_zone;
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
		Matrix_memset(prx, 0, (R_xlen_t) rm * rm, sizeof(double));
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

		SEXP bdim = GET_SLOT(b, Matrix_DimSym);
		int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
			rn = (btrans != 'N') ? bm : bn;

		if (rk != ((btrans != 'N') ? bn : bm))
			error(_("non-conformable arguments"));
		if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");

		SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

		char rcl[] = ".geMatrix";
		rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
		SEXP r = PROTECT(newObject(rcl));

		SEXP rdim = GET_SLOT(r, Matrix_DimSym);
		int *prdim = INTEGER(rdim);
		prdim[0] = rm;
		prdim[1] = rn;

		SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
			bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
			rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
		matmultDN(rdimnames,
		          adimnames, (atrans != 'N') ? 1 : 0,
		          bdimnames, (btrans != 'N') ? 0 : 1);
		UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

		if (rm > 0 && rn > 0) {
		SEXP rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
		if (TYPEOF(ax) == CPLXSXP) {
		Rcomplex *prx = COMPLEX(rx);
		if (rk == 0)
		Matrix_memset(prx, 0, (R_xlen_t) rm * rn, sizeof(Rcomplex));
		else {
		SEXP bx = PROTECT(GET_SLOT(b, Matrix_xSym));
		Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx),
			zero = Matrix_zzero, one = Matrix_zone;
		F77_CALL(zgemm)(&atrans, &btrans, &rm, &rn, &rk,
		                & one, pax, &am, pbx, &bm,
		                &zero, prx, &rm FCONE FCONE);
		UNPROTECT(1); /* bx */
		}
		} else {
		double *prx = REAL(rx);
		if (rk == 0)
		Matrix_memset(prx, 0, (R_xlen_t) rm * rn, sizeof(double));
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
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(get_symmetrized_DimNames(a, -1)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
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
		rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
	char aul = *CHAR(STRING_ELT(GET_SLOT(a, Matrix_uploSym), 0));
	int i,
		d     = (aside == 'L') ? rn : rm,
		binc  = (aside == 'L') ? bm :  1,
		bincp = (aside == 'L') ?  1 : bm,
		rinc  = (aside == 'L') ?  1 : rm,
		rincp = (aside == 'L') ? rm :  1;

	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		zero = Matrix_zzero, one = Matrix_zone;
	char act = *CHAR(STRING_ELT(GET_SLOT(a, Matrix_transSym), 0));
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
		zvconj(pax, (R_xlen_t) rk * rk); /* in place */
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
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = ".geMatrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(get_symmetrized_DimNames(a, -1)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
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
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	char aul = *CHAR(STRING_ELT(GET_SLOT(a, Matrix_uploSym), 0));
	int i,
		d     = ( aside == 'L'                    ) ? rn : rm,
		binc  = ((aside == 'L') == (btrans != 'N')) ? bm :  1,
		bincp = ((aside == 'L') == (btrans != 'N')) ?  1 : bm,
		rinc  = ( aside == 'L'                    ) ?  1 : rm,
		rincp = ( aside == 'L'                    ) ? rm :  1;

	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		zero = Matrix_zzero, one = Matrix_zone;
	char act = *CHAR(STRING_ELT(GET_SLOT(a, Matrix_transSym), 0));
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
		zvconj(prx, (R_xlen_t) rm * rn); /* in place */
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
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	rcl[1] = (triangular > 0) ? 't' : 'g';
	rcl[2] = (triangular > 0) ? 'r' : 'e';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aside == 'L')
	matmultDN(rdimnames, adimnames, atrans != 'N', bdimnames, btrans == 'N');
	else
	matmultDN(rdimnames, bdimnames, btrans != 'N', adimnames, atrans == 'N');
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = *CHAR(STRING_ELT(auplo, 0));
	if (triangular > 0 && ((atrans != 'N') ? aul == 'U' : aul != 'U')) {
		if (atrans != 'N')
			auplo = mkString("L");
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	SEXP adiag = GET_SLOT(a, Matrix_diagSym);
	char adi = *CHAR(STRING_ELT(adiag, 0));
	if (triangular > 1 && adi != 'N') {
		PROTECT(adiag);
		SET_SLOT(r, Matrix_diagSym, adiag);
		UNPROTECT(1); /* adiag */
	}

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
		rx = PROTECT(allocVector(TYPEOF(ax), (R_xlen_t) rm * rn));
	if (TYPEOF(ax) == CPLXSXP) {
	Rcomplex *pax = COMPLEX(ax), *pbx = COMPLEX(bx), *prx = COMPLEX(rx),
		one = Matrix_zone;
	ztrans2(prx, pbx, bm, bn, btrans);
	F77_CALL(ztrmm)(&aside, &aul, &atrans, &adi, &rm, &rn,
	                &one, pax, &rk, prx, &rm FCONE FCONE FCONE FCONE);
	} else {
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		one = 1.0;
	dtrans2(prx, pbx, bm, bn, btrans);
	F77_CALL(dtrmm)(&aside, &aul, &atrans, &adi, &rm, &rn,
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
	SEXP adim = GET_SLOT(a, Matrix_DimSym);
	int rk = INTEGER(adim)[0];

	SEXP bdim = GET_SLOT(b, Matrix_DimSym);
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans != 'N') ? bn : bm, rn = (btrans != 'N') ? bm : bn;

	if (rk != (((aside == 'L') == (btrans != 'N')) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));

	char rcl[] = "...Matrix";
	rcl[0] = (TYPEOF(ax) == CPLXSXP) ? 'z' : 'd';
	rcl[1] = (triangular > 0) ? 't' : 'g';
	rcl[2] = (triangular > 0) ? 'r' : 'e';
	SEXP r = PROTECT(newObject(rcl));

	SEXP rdim = GET_SLOT(r, Matrix_DimSym);
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;

	SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aside == 'L')
	matmultDN(rdimnames, adimnames, atrans != 'N', bdimnames, btrans == 'N');
	else
	matmultDN(rdimnames, bdimnames, btrans != 'N', adimnames, atrans == 'N');
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	SEXP auplo = GET_SLOT(a, Matrix_uploSym);
	char aul = *CHAR(STRING_ELT(auplo, 0));
	if (triangular > 0 && ((atrans != 'N') ? aul == 'U' : aul != 'U')) {
		if (atrans != 'N')
			auplo = mkString("L");
		PROTECT(auplo);
		SET_SLOT(r, Matrix_uploSym, auplo);
		UNPROTECT(1); /* auplo */
	}

	SEXP adiag = GET_SLOT(a, Matrix_diagSym);
	char adi = *CHAR(STRING_ELT(adiag, 0));
	if (triangular > 1 && adi != 'N') {
		PROTECT(adiag);
		SET_SLOT(r, Matrix_diagSym, adiag);
		UNPROTECT(1); /* adiag */
	}

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
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
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
	ztrans2(prx, pbx, bm, bn, btransp);
	if (aside != 'L' && atrans == 'C' && btrans == 'N')
		zvconj(prx, (R_xlen_t) rm * rn); /* in place */
	for (i = 0; i < rn; ++i) {
	F77_CALL(ztpmv)(&aul, &atransp, &adi, &rk,
	                pax, prx, &rinc FCONE FCONE FCONE);
	prx += rincp;
	}
	if (aside != 'L' &&
	    ((btrans != 'N') ? btrans : ((atrans != 'N') ? atrans : 'T')) == 'C')
		zvconj(prx, (R_xlen_t) rm * rn); /* in place */
	} else {
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx);
	dtrans2(prx, pbx, bm, bn, btransp);
	for (i = 0; i < rn; ++i) {
	F77_CALL(dtpmv)(&aul, &atransp, &adi, &rk,
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

SEXP R_dense_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans)
{
	char
		xtrans_ = *CHAR(STRING_ELT(xtrans, 0)),
		ytrans_ = *CHAR(STRING_ELT(ytrans, 0)),
		ztrans_ = 'N';
	int m, n, v;
	matmultDim(x, y, &xtrans_, &ytrans_, &ztrans_, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (TYPEOF(x) != S4SXP) {
		REPROTECT(x = matrix_as_dense(x, ",ge", '\0', '\0', '\0', (xtrans_ != 'N') ? 0 : 1, 0), xpid);
		if (v % 2) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (xtrans_ != 'N') ? 1 : 0, R_NilValue);
			xtrans_ = 'N';
		}
	}
	if (TYPEOF(y) != S4SXP && y != R_NilValue) {
		REPROTECT(y = matrix_as_dense(y, ",ge", '\0', '\0', '\0', (ytrans_ != 'N') ? 0 : 1, 0), ypid);
		if (v > 1) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (ytrans_ != 'N') ? 1 : 0, R_NilValue);
			ytrans_ = 'N';
		}
	}

	static const char *valid[] = { VALID_DENSE, "" };
	const char *xcl = NULL, *ycl = NULL;
	int ivalid;
	ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	xcl = valid[ivalid];
	if (y != R_NilValue) {
	ivalid = R_check_class_etc(y, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(y, __func__);
	ycl = valid[ivalid];
	}

	char kind = (xcl[0] == 'z' || (y != R_NilValue && ycl[0] == 'z'))
		? 'z' : 'd';
	if (xcl[0] != kind) {
		REPROTECT(x = dense_as_kind(x, xcl, kind, 0), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
	}
	if (y != R_NilValue) {
	if (ycl[0] != kind) {
		REPROTECT(y = dense_as_kind(y, ycl, kind, 0), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
	}
	}

	if (y == R_NilValue) {
		REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
		x = geMatrix_matmult(x, y, xtrans_, ytrans_);
	} else if (xcl[1] == 'g' && ycl[1] == 'g') {
		x = geMatrix_matmult(x, y, xtrans_, ytrans_);
	} else if (xcl[1] == 'g' || ycl[1] == 'g') {
		x = (xcl[1] == 'g')
			? ((ycl[1] == 's')
			   ? ((ycl[2] != 'p')
			      ? syMatrix_matmult(y, x, ytrans_, xtrans_, 'R')
			      : spMatrix_matmult(y, x, ytrans_, xtrans_, 'R'))
			   : ((ycl[2] != 'p')
			      ? trMatrix_matmult(y, x, ytrans_, xtrans_, 'R', 0)
			      : tpMatrix_matmult(y, x, ytrans_, xtrans_, 'R', 0)))
			: ((xcl[1] == 's')
			   ? ((xcl[2] != 'p')
			      ? syMatrix_matmult(x, y, xtrans_, ytrans_, 'L')
			      : spMatrix_matmult(x, y, xtrans_, ytrans_, 'L'))
			   : ((xcl[2] != 'p')
			      ? trMatrix_matmult(x, y, xtrans_, ytrans_, 'L', 0)
			      : tpMatrix_matmult(x, y, xtrans_, ytrans_, 'L', 0)));
	} else if (xcl[1] == 's' && ycl[1] == 's') {
		if (xcl[2] == 'p' && ycl[2] == 'p') {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = spMatrix_matmult(x, y, xtrans_, ytrans_, 'L');
		} else if (xcl[2] == 'p') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = syMatrix_matmult(y, x, ytrans_, xtrans_, 'R');
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = syMatrix_matmult(x, y, xtrans_, ytrans_, 'L');
		}
	} else if (xcl[1] == 's' || ycl[1] == 's') {
		if (xcl[1] == 's') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = (ycl[2] != 'p')
				? trMatrix_matmult(y, x, ytrans_, xtrans_, 'R', 0)
				: tpMatrix_matmult(y, x, ytrans_, xtrans_, 'R', 0);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = (xcl[2] != 'p')
				? trMatrix_matmult(x, y, xtrans_, ytrans_, 'L', 0)
				: tpMatrix_matmult(x, y, xtrans_, ytrans_, 'L', 0);
		}
	} else {
		SEXP
			xuplo = PROTECT(GET_SLOT(x, Matrix_uploSym)),
			yuplo = PROTECT(GET_SLOT(y, Matrix_uploSym)),
			xdiag = PROTECT(GET_SLOT(x, Matrix_diagSym)),
			ydiag = PROTECT(GET_SLOT(y, Matrix_diagSym));
		char
			xul = *CHAR(STRING_ELT(xuplo, 0)),
			yul = *CHAR(STRING_ELT(yuplo, 0)),
			xdi = *CHAR(STRING_ELT(xdiag, 0)),
			ydi = *CHAR(STRING_ELT(ydiag, 0));
		if (xtrans_ != 'N')
			xul = (xul == 'U') ? 'L' : 'U';
		if (ytrans_ != 'N')
			yul = (yul == 'U') ? 'L' : 'U';
		int triangular = (xul != yul) ? 0 : ((xdi != ydi || xdi == 'N') ? 1 : 2);
		UNPROTECT(4); /* ydiag, xdiag, yuplo, xuplo */

		if (xcl[2] == 'p' && ycl[2] == 'p') {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = tpMatrix_matmult(x, y, xtrans_, ytrans_, 'L', triangular);
		} else if (xcl[2] == 'p') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = trMatrix_matmult(y, x, ytrans_, xtrans_, 'R', triangular);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = trMatrix_matmult(x, y, xtrans_, ytrans_, 'L', triangular);
		}
	}

	UNPROTECT(2); /* y, x */
	return x;
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
	char zcl[] = "..CMatrix";

	if (y == R_NilValue) {

		zcl[0] = (boolean) ? 'n' : ((X->xtype != CHOLMOD_COMPLEX) ? 'd' : 'z');
		zcl[1] = (boolean) ? 's' : ((X->xtype != CHOLMOD_COMPLEX) ? 'p' : ((((xtrans != 'N') ? xtrans : ytrans) == 'C') ? 'p' : 's'));
#ifndef MATRIX_ENABLE_POSDEF
		if (zcl[1] == 'p')
			zcl[1] = 's';
#endif

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
		PROTECT(z = CHS2M(Z, !boolean, zcl[1]));
		cholmod_free_sparse(&Z, &c);

		SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
			zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
		symDN(zdimnames, xdimnames, (xtrans != 'N') ? 1 : 0);
		UNPROTECT(2); /* zdimnames, xdimnames */

		if (ztrans != 'N') {
			SEXP zuplo = PROTECT(mkString("L"));
			SET_SLOT(z, Matrix_uploSym, zuplo);
			UNPROTECT(1); /* zuplo */
		}

		if (zcl[1] == 's' && zcl[0] == 'z') {
			SEXP ztrans = PROTECT(mkString("T"));
			SET_SLOT(z, Matrix_transSym, ztrans);
			UNPROTECT(1); /* ztrans */
		}

	} else {

		Y = M2CHS(y, !boolean);

		if (((xtrans != 'N') ? X->nrow : X->ncol) !=
		    ((ytrans != 'N') ? Y->ncol : Y->nrow))
			error(_("non-conformable arguments"));

		zcl[0] = (boolean) ? 'n' : ((X->xtype != CHOLMOD_COMPLEX && Y->xtype != CHOLMOD_COMPLEX) ? 'd' : 'z');
		zcl[1] = (triangular != 0) ? 't' : 'g';

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
			cholmod_band_inplace(-Z->nrow, -1, !boolean, Z, &c);
		else if (triangular > 1)
			cholmod_band_inplace( 1,  Z->ncol, !boolean, Z, &c);
		PROTECT(z = CHS2M(Z, !boolean, zcl[1]));
		cholmod_free_sparse(&Z, &c);

		SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
			ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
			zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
		matmultDN(zdimnames,
		          xdimnames, (xtrans != 'N') ? 1 : 0,
		          ydimnames, (ytrans != 'N') ? 0 : 1);
		UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

		if (triangular < 0) {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(z, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}

		if (triangular < -1 || triangular > 1) {
			SEXP diag = PROTECT(mkString("U"));
			SET_SLOT(z, Matrix_diagSym, diag);
			UNPROTECT(1); /* diag */
		}

	}

	if (ztrans != 'N')
		z = sparse_transpose(z, zcl, 1);

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
		error(_("non-conformable arguments"));
	int m = (int) ((xtrans != 'N') ? X->ncol : X->nrow), n = (int) Y->ncol;
	if ((Matrix_int_fast64_t) m * n > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	if (X->xtype == CHOLMOD_COMPLEX && xtrans == 'T')
		CONJ2(X->x, X->nrow, X->ncol);

	char zcl[] = "...Matrix";
	zcl[0] = (X->xtype != CHOLMOD_COMPLEX && Y->xtype != CHOLMOD_COMPLEX) ? 'd' : 'z';
	zcl[1] = (triangular != 0) ? 't' : 'g';
	zcl[2] = (triangular != 0) ? 'r' : 'e';
	SEXP z = PROTECT(newObject(zcl));

	SEXP zdim = GET_SLOT(z, Matrix_DimSym);
	INTEGER(zdim)[0] = (ztrans != 'N') ? n : m;
	INTEGER(zdim)[1] = (ztrans != 'N') ? m : n;

	SEXP xdimnames = (symmetric != 0)
		? PROTECT(get_symmetrized_DimNames(x, -1))
		: PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
		ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
		zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
	if (ztrans != 'N')
	matmultDN(zdimnames,
	          ydimnames, (ytrans != 'N') ? 0 : 1,
	          xdimnames, (xtrans != 'N') ? 1 : 0);
	else
	matmultDN(zdimnames,
	          xdimnames, (xtrans != 'N') ? 1 : 0,
	          ydimnames, (ytrans != 'N') ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

	if (triangular != 0 && (ztrans != 'N') == (triangular > 0)) {
		SEXP zuplo = PROTECT(mkString("L"));
		SET_SLOT(z, Matrix_uploSym, zuplo);
		UNPROTECT(1); /* zuplo */
	}

	if (triangular < -1 || triangular > 1) {
		SEXP zdiag = PROTECT(mkString("U"));
		SET_SLOT(z, Matrix_diagSym, zdiag);
		UNPROTECT(1); /* zdiag */
	}

	double alpha[2] = { 1.0, 0.0 }, beta[2] = { 0.0, 0.0 };
	cholmod_dense *Z = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(Z, 0, sizeof(cholmod_dense));
	Z->nrow = (size_t) m;
	Z->ncol = (size_t) n;
	Z->d = Z->nrow;
	Z->nzmax = Z->nrow * Z->ncol;
	Z->xtype = X->xtype;
	Z->dtype = X->dtype;

	SEXP zx;
	if (zcl[0] == 'z') {
	PROTECT(zx = allocVector(CPLXSXP, (R_xlen_t) m * n));
	Z->x = (ztrans != 'N')
		? (Rcomplex *) R_alloc((size_t) m * n, sizeof(Rcomplex))
		: COMPLEX(zx);
	} else {
	PROTECT(zx = allocVector(REALSXP, (R_xlen_t) m * n));
	Z->x = (ztrans != 'N')
		? (double *) R_alloc((size_t) m * n, sizeof(double))
		: REAL(zx);
	}
	cholmod_sdmult(X, xtrans != 'N', alpha, beta, Y, Z, &c);
	if (ztrans != 'N') {
	if (zcl[0] == 'z')
	ztrans2(COMPLEX(zx), (Rcomplex *) Z->x, m, n, ztrans);
	else
	dtrans2(   REAL(zx), (  double *) Z->x, m, n, ztrans);
	}
	SET_SLOT(z, Matrix_xSym, zx);
	UNPROTECT(1); /* zx */

	UNPROTECT(1); /* z */
	return z;
}

SEXP R_sparse_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans, SEXP ztrans,
                      SEXP boolean)
{
	if (TYPEOF(boolean) != LGLSXP || LENGTH(boolean) < 1)
		error(_("invalid '%s' to '%s'"), "boolean", __func__);
	int boolean_ = LOGICAL(boolean)[0];

	char
		xtrans_ = *CHAR(STRING_ELT(xtrans, 0)),
		ytrans_ = *CHAR(STRING_ELT(ytrans, 0)),
		ztrans_ = *CHAR(STRING_ELT(ztrans, 0));
	int m, n, v;
	matmultDim(x, y, &xtrans_, &ytrans_, &ztrans_, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (TYPEOF(x) != S4SXP) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(x = matrix_as_dense (x, ",ge", '\0', '\0', '\0', (xtrans_ != 'N') ? 0 : 1, 0), xpid);
		else if (xtrans_ != 'N')
		REPROTECT(x = matrix_as_sparse(x, "ngR", '\0', '\0', '\0', 0), xpid);
		else
		REPROTECT(x = matrix_as_sparse(x, "ngC", '\0', '\0', '\0', 0), xpid);
		if (v % 2) {
			/* Discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (xtrans_ != 'N') ? 1 : 0, R_NilValue);
			xtrans_ = 'N';
		}
	}
	if (TYPEOF(y) != S4SXP && y != R_NilValue) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(y = matrix_as_dense (y, ",ge", '\0', '\0', '\0', (ytrans_ != 'N') ? 0 : 1, 0), ypid);
		else if (ytrans_ != 'N')
		REPROTECT(y = matrix_as_sparse(y, "ngR", '\0', '\0', '\0', 0), ypid);
		else
		REPROTECT(y = matrix_as_sparse(y, "ngC", '\0', '\0', '\0', 0), ypid);
		if (v > 1) {
			/* Discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (ytrans_ != 'N') ? 1 : 0, R_NilValue);
			ytrans_ = 'N';
		}
	}

	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, VALID_DENSE, "" };
	const char *xcl = NULL, *ycl = NULL;
	int ivalid;
	ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	xcl = valid[ivalid];
	if (y != R_NilValue) {
	ivalid = R_check_class_etc(y, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(y, __func__);
	ycl = valid[ivalid];
	}
	if (boolean_ == NA_LOGICAL)
		boolean_ = xcl[0] == 'n' && (y == R_NilValue || ycl[0] == 'n');
	char kind = (boolean_) ? 'n' :
		((xcl[0] != 'z' && (y == R_NilValue || ycl[0] != 'z')) ? 'd' : 'z');

	if (xcl[2] != 'C' && xtrans_ != 'N') {
		if (xcl[2] != 'R' && xcl[2] != 'T') {
			REPROTECT(x = dense_as_sparse(x, xcl, 'R'), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		REPROTECT(x = sparse_transpose(x, xcl, 1), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
		xtrans_ = 'N';
	}
	if (xcl[2] != 'C') {
		if (xcl[2] != 'R' && xcl[2] != 'T')
			REPROTECT(x = dense_as_sparse(x, xcl, 'C'), xpid);
		else
			REPROTECT(x = sparse_as_Csparse(x, xcl), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
	}
	if (xtrans_ != 'N' && xcl[1] == 's' &&
	    (xcl[0] != 'z' || *CHAR(STRING_ELT(GET_SLOT(x, Matrix_transSym), 0)) == xtrans_))
		xtrans_ = 'N';
	if (xcl[0] != kind) {
		if (boolean_)
			REPROTECT(x = sparse_drop0(x, xcl, 0.0), xpid);
		else {
			REPROTECT(x = sparse_as_kind(x, xcl, kind), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
	}

	if (y == R_NilValue) {
		REPROTECT(x = sparse_as_general(x, xcl), xpid);
		x = gCgCMatrix_matmult(
			x, y, xtrans_, ytrans_, ztrans_, 0, boolean_);
		UNPROTECT(2); /* y, x */
		return x;
	}

	int triangular = 0;
	if (xcl[1] == 't' && ycl[1] == 't') {
		SEXP
			xuplo = PROTECT(GET_SLOT(x, Matrix_uploSym)),
			yuplo = PROTECT(GET_SLOT(y, Matrix_uploSym)),
			xdiag = PROTECT(GET_SLOT(x, Matrix_diagSym)),
			ydiag = PROTECT(GET_SLOT(y, Matrix_diagSym));
		char
			xul = *CHAR(STRING_ELT(xuplo, 0)),
			yul = *CHAR(STRING_ELT(yuplo, 0)),
			xdi = *CHAR(STRING_ELT(xdiag, 0)),
			ydi = *CHAR(STRING_ELT(ydiag, 0));
		if (xtrans_ != 'N')
			xul = (xul == 'U') ? 'L' : 'U';
		if (ytrans_ != 'N')
			yul = (yul == 'U') ? 'L' : 'U';
		triangular = (xul != yul) ? 0 : ((xdi != ydi || xdi == 'N') ? 1 : 2);
		if (xul != 'U')
			triangular = -triangular;
		UNPROTECT(4); /* ydiag, xdiag, yuplo, xuplo */
	}

	if (!boolean_ && ycl[2] != 'C' && ycl[2] != 'R' && ycl[2] != 'T') {
		if (xcl[1] == 's' && xcl[0] == 'z') {
			SEXP xtrans = GET_SLOT(x, Matrix_transSym);
			char xct = *CHAR(STRING_ELT(xtrans, 0));
			if (xct != 'C') {
				REPROTECT(x = sparse_as_general(x, xcl), xpid);
				xcl = valid[R_check_class_etc(x, valid)];
			}
		}
		int symmetric = xcl[1] == 's';
		if (symmetric) {
			SEXP xuplo = GET_SLOT(x, Matrix_uploSym);
			char xul = *CHAR(STRING_ELT(xuplo, 0));
			if (xul != 'U')
				symmetric = -symmetric;
		}
		if (ycl[0] != kind) {
			REPROTECT(y = dense_as_kind(y, ycl, kind, 0), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		REPROTECT(x = sparse_diag_U2N(x, xcl), xpid);
		REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
		x = gCgeMatrix_matmult(
			x, y, xtrans_, ytrans_, ztrans_, triangular, symmetric);
		UNPROTECT(2); /* y, x */
		return x;
	}

	if (ycl[2] != 'C' && ytrans_ != 'N') {
		if (ycl[2] != 'R' && ycl[2] != 'T') {
			REPROTECT(y = dense_as_sparse(y, ycl, 'R'), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		REPROTECT(y = sparse_transpose(y, ycl, 1), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
		ytrans_ = 'N';
	}
	if (ycl[2] != 'C') {
		if (ycl[2] != 'R' && ycl[2] != 'T')
			REPROTECT(y = dense_as_sparse(y, ycl, 'C'), ypid);
		else
			REPROTECT(y = sparse_as_Csparse(y, ycl), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
	}
	if (ytrans_ != 'N' && ycl[1] == 's' &&
	    (ycl[0] != 'z' || *CHAR(STRING_ELT(GET_SLOT(y, Matrix_transSym), 0)) == ytrans_))
		ytrans_ = 'N';
	if (ycl[0] != kind) {
		if (boolean_)
			REPROTECT(y = sparse_drop0(y, ycl, 0.0), ypid);
		else {
			REPROTECT(y = sparse_as_kind(y, ycl, kind), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
	}

	REPROTECT(x = sparse_as_general(x, xcl), xpid);
	REPROTECT(y = sparse_as_general(y, ycl), ypid);
	x = gCgCMatrix_matmult(
		x, y, xtrans_, ytrans_, ztrans_, triangular, boolean_);
	UNPROTECT(2); /* y, x */
	return x;
}

#define MULTIPLY_COMPLEX(_X_, _D_) \
	do { \
		tmp = (_X_); \
		(_X_).r = tmp.r * (_D_).r - tmp.i * (_D_).i; \
		(_X_).i = tmp.r * (_D_).i + tmp.i * (_D_).r; \
	} while (0)
#define MULTIPLY_REAL(_X_, _D_) \
	(_X_) = (_X_)  * (_D_)
#define MULTIPLY_LOGICAL(_X_, _D_) \
	(_X_) = (_X_) && (_D_)

#define SCALE_CASES(_J_) \
	do { \
		switch (TYPEOF(d)) { \
		case CPLXSXP: \
		{ \
			Rcomplex tmp; \
			SCALE(z, Rcomplex, COMPLEX, MULTIPLY_COMPLEX, _J_); \
			break; \
		} \
		case REALSXP: \
			SCALE(d, double, REAL, MULTIPLY_REAL, _J_); \
			break; \
		case LGLSXP: \
			SCALE(i, int, LOGICAL, MULTIPLY_LOGICAL, _J_); \
			break; \
		default: \
			break; \
		} \
	} while (0)

static
void dense_colscale(SEXP obj, SEXP d, int m, int n, char ul, char di)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = XLENGTH(x) != (Matrix_int_fast64_t) m * n;

#define SCALE(_PREFIX_, _CTYPE_, _PTR_, _OP_, _J_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pd = _PTR_(d); \
		if (ul == '\0') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < m; ++i) { \
					_OP_(*px, pd[_J_]); \
					++px; \
				} \
			} \
		} else if (ul == 'U') { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i <= j; ++i) { \
					_OP_(*px, pd[_J_]); \
					++px; \
				} \
				if (!packed) \
					px += m - j - 1; \
			} \
		} else { \
			for (j = 0; j < n; ++j) { \
				if (!packed) \
					px += j; \
				for (i = j; i < m; ++i) { \
					_OP_(*px, pd[_J_]); \
					++px; \
				} \
			} \
		} \
		if (di != '\0' && di != 'N') { \
			if (!packed) \
				_PREFIX_ ## dcopy2(_PTR_(x), pd, n, n,     'U', 'N'); \
			else \
				_PREFIX_ ## dcopy1(_PTR_(x), pd, n, n, ul, 'U', 'N'); \
		} \
	} while (0)

	SCALE_CASES(j);
	return;
}

static
void dense_rowscale(SEXP obj, SEXP d, int m, int n, char ul, char di)
{
	SEXP x = GET_SLOT(obj, Matrix_xSym);
	int i, j, packed = XLENGTH(x) != (Matrix_int_fast64_t) m * n;
	SCALE_CASES(i);

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
	int *pp = INTEGER(p) + 1, n = (int) (XLENGTH(p) - 1), j, k = 0, kend;
	UNPROTECT(2); /* p, x */

#define SCALE(_PREFIX_, _CTYPE_, _PTR_, _OP_, _J_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pd = _PTR_(d); \
		for (j = 0; j < n; ++j) { \
			kend = pp[j]; \
			while (k < kend) { \
				_OP_(*px, *pd); \
				++px; \
				++k; \
			} \
			++pd; \
		} \
	} while (0)

	SCALE_CASES();

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
	int *pi = INTEGER(i), k, nnz = INTEGER(p)[XLENGTH(p) - 1];
	UNPROTECT(3); /* i, p, x */

#define SCALE(_PREFIX_, _CTYPE_, _PTR_, _OP_, _J_) \
	do { \
		_CTYPE_ *px = _PTR_(x), *pd = _PTR_(d); \
		for (k = 0; k < nnz; ++k) { \
			_OP_(*px, pd[*pi]); \
			++px; \
			++pi; \
		} \
	} while (0)

	SCALE_CASES();
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
	R_xlen_t k, nnz = XLENGTH(i);
	UNPROTECT(2); /* i, x */
	SCALE_CASES();

#undef SCALE

	return;
}

SEXP R_diagonal_matmult(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans,
                        SEXP boolean)
{
	SEXP x_ = x, y_ = y; /* for later pointer comparison */

	if (TYPEOF(boolean) != LGLSXP || LENGTH(boolean) < 1)
		error(_("invalid '%s' to '%s'"), "boolean", __func__);
	int boolean_ = LOGICAL(boolean)[0];

	char
		xtrans_ = *CHAR(STRING_ELT(xtrans, 0)),
		ytrans_ = *CHAR(STRING_ELT(ytrans, 0)),
		ztrans_ = 'N';
	int m, n, v;
	matmultDim(x, y, &xtrans_, &ytrans_, &ztrans_, &m, &n, &v);

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (TYPEOF(x) != S4SXP) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(x = matrix_as_dense(x, ",ge", '\0', '\0', '\0', (xtrans_ != 'N') ? 0 : 1, 2), xpid);
		else
		REPROTECT(x = matrix_as_dense(x, "nge", '\0', '\0', '\0', (xtrans_ != 'N') ? 0 : 1, 2), xpid);
		if (v % 2) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (xtrans_ != 'N') ? 1 : 0, R_NilValue);
			xtrans_ = 'N';
		}
	}
	if (TYPEOF(y) != S4SXP) {
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(y = matrix_as_dense(y, ",ge", '\0', '\0', '\0', (ytrans_ != 'N') ? 0 : 1, 2), ypid);
		else
		REPROTECT(y = matrix_as_dense(y, "nge", '\0', '\0', '\0', (ytrans_ != 'N') ? 0 : 1, 2), ypid);
		if (v > 1) {
			/* Vector: discard names and don't transpose again */
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (ytrans_ != 'N') ? 1 : 0, R_NilValue);
			ytrans_ = 'N';
		}
	}

	static const char *valid[] = {
		VALID_DIAGONAL,
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, VALID_DENSE, "" };
	const char *xcl = NULL, *ycl = NULL;
	int ivalid;
	ivalid = R_check_class_etc(x, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(x, __func__);
	xcl = valid[ivalid];
	ivalid = R_check_class_etc(y, valid);
	if (ivalid < 0)
		ERROR_INVALID_CLASS(y, __func__);
	ycl = valid[ivalid];
	if (boolean_ == NA_LOGICAL)
		boolean_ = xcl[0] == 'n' && ycl[0] == 'n';
	char kind = (boolean_) ? 'n' :
		((xcl[0] != 'z' && ycl[0] != 'z') ? 'd' : 'z');

	int mg = -1, id = -1;
	if (xcl[2] == 'i') {
		mg = 0;
		id = *CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0)) != 'N';
	} else if (ycl[2] == 'i') {
		mg = 1;
		id = *CHAR(STRING_ELT(GET_SLOT(y, Matrix_diagSym), 0)) != 'N';
	} else
		error(_("should never happen"));

	if (xtrans_ != 'N' && xcl[1] == 's' &&
	    (xcl[0] != 'z' || *CHAR(STRING_ELT(GET_SLOT(x, Matrix_transSym), 0)) == xtrans_))
		xtrans_ = 'N';
	if (ytrans_ != 'N' && ycl[1] == 's' &&
	    (ycl[0] != 'z' || *CHAR(STRING_ELT(GET_SLOT(y, Matrix_transSym), 0)) == ytrans_))
		ytrans_ = 'N';

	char ks = (boolean_) ? 'l' : kind, kd = kind;
	switch (xcl[2]) {
	case 'i':
		if (!id && xcl[0] != ks) {
			REPROTECT(x = diagonal_as_kind(x, xcl, ks), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		break;
	case 'C':
	case 'R':
	case 'T':
		if (xcl[0] != ks) {
			REPROTECT(x = sparse_as_kind(x, xcl, ks), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (!id && xcl[1] == 's') {
			REPROTECT(x = sparse_as_general(x, xcl), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (!id && xcl[1] == 't')
			REPROTECT(x = sparse_diag_U2N(x, xcl), xpid);
		if (xtrans_ != 'N') {
			REPROTECT(x = sparse_transpose(x, xcl, 0), xpid);
			xtrans_ = 'N';
		}
		break;
	default:
		if (xcl[0] != kd) {
			REPROTECT(x = dense_as_kind(x, xcl, kd, 1), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (!id && xcl[1] == 's') {
			REPROTECT(x = dense_as_general(x, xcl, x == x_), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (xtrans_ != 'N') {
			REPROTECT(x = dense_transpose(x, xcl, xtrans_), xpid);
			xtrans_ = 'N';
		}
		break;
	}
	switch (ycl[2]) {
	case 'i':
		if (!id && ycl[0] != ks) {
			REPROTECT(y = diagonal_as_kind(y, ycl, ks), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		break;
	case 'C':
	case 'R':
	case 'T':
		if (ycl[0] != ks) {
			REPROTECT(y = sparse_as_kind(y, ycl, ks), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (ytrans_ != 'N') {
			REPROTECT(y = sparse_transpose(y, ycl, 0), ypid);
			ytrans_ = 'N';
		}
		if (!id && ycl[1] == 's') {
			REPROTECT(y = sparse_as_general(y, ycl), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (!id && ycl[1] == 't')
			REPROTECT(y = sparse_diag_U2N(y, ycl), ypid);
		break;
	default:
		if (ycl[0] != kd) {
			REPROTECT(y = dense_as_kind(y, ycl, kd, 1), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (ytrans_ != 'N') {
			REPROTECT(y = dense_transpose(y, ycl, ytrans_), ypid);
			ytrans_ = 'N';
		}
		if (!id && ycl[1] == 's') {
			REPROTECT(y = dense_as_general(y, ycl, y == y_), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		break;
	}

	const char *zcl = (mg == 0) ? ycl : xcl;
	PROTECT_INDEX zpid;
	SEXP z = newObject(zcl);
	PROTECT_WITH_INDEX(z, &zpid);

	SEXP zdim = PROTECT(GET_SLOT(z, Matrix_DimSym));
	int *pzdim = INTEGER(zdim);
	pzdim[0] = m;
	pzdim[1] = n;
	UNPROTECT(1); /* zdim */

	SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
		ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
		zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
	matmultDN(zdimnames,
	          xdimnames, (xtrans_ != 'N') ? 1 : 0,
	          ydimnames, (ytrans_ != 'N') ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */

	char ul = '\0', di = '\0';
	if (zcl[1] != 'g') {
		SEXP uplo = PROTECT(GET_SLOT((mg == 0) ? y : x, Matrix_uploSym));
		ul = *CHAR(STRING_ELT(uplo, 0));
		if (ul != 'U')
			SET_SLOT(z, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (zcl[1] == 't') {
		SEXP diag = PROTECT(GET_SLOT((mg == 0) ? y : x, Matrix_diagSym));
		di = *CHAR(STRING_ELT(diag, 0));
		if (di != 'N' && id)
			SET_SLOT(z, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}

	if (zcl[2] == 'C' || zcl[2] == 'R' || zcl[2] == 'T') {
		if (zcl[2] != 'T') {
			SEXP p = PROTECT(GET_SLOT((mg == 0) ? y : x, Matrix_pSym));
			SET_SLOT(z, Matrix_pSym, p);
			UNPROTECT(1); /* p */
		}
		if (zcl[2] != 'R') {
			SEXP i = PROTECT(GET_SLOT((mg == 0) ? y : x, Matrix_iSym));
			SET_SLOT(z, Matrix_iSym, i);
			UNPROTECT(1); /* i */
		}
		if (zcl[2] != 'C') {
			SEXP j = PROTECT(GET_SLOT((mg == 0) ? y : x, Matrix_jSym));
			SET_SLOT(z, Matrix_jSym, j);
			UNPROTECT(1); /* j */
		}
	}

	SEXP x0 = PROTECT(GET_SLOT((mg == 0) ? y : x, Matrix_xSym));
	if (id || ((mg == 0) ? y != y_ : x != x_))
		SET_SLOT(z, Matrix_xSym, x0);
	else {
		SEXP x1 = PROTECT(allocVector(TYPEOF(x0), XLENGTH(x0)));
		switch (kind) {
		case 'z':
			Matrix_memcpy(COMPLEX(x1), COMPLEX(x0), XLENGTH(x0), sizeof(Rcomplex));
			break;
		case 'd':
			Matrix_memcpy(   REAL(x1),    REAL(x0), XLENGTH(x0), sizeof(  double));
			break;
		default:
			Matrix_memcpy(LOGICAL(x1), LOGICAL(x0), XLENGTH(x0), sizeof(     int));
			break;
		}
		SET_SLOT(z, Matrix_xSym, x1);
		UNPROTECT(1); /* x1 */
	}
	UNPROTECT(1); /* x0 */

	if (!id) {
		SEXP d = PROTECT(GET_SLOT((mg == 0) ? x : y, Matrix_xSym));
		switch (zcl[2]) {
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
				  dense_rowscale(z, d, m, n, ul, di);
			else
				  dense_colscale(z, d, m, n, ul, di);
			break;
		}
		UNPROTECT(1); /* d */
	}

	if (boolean_ && (zcl[2] == 'C' || zcl[2] == 'R' || zcl[2] == 'T')) {
		REPROTECT(z = sparse_drop0(z, zcl, 0.0), zpid);
		REPROTECT(z = sparse_as_kind(z, zcl, 'n'), zpid);
	}

	UNPROTECT(3); /* z, y, x */
	return z;
}
