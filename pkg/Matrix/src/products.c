#include "products.h"
#include "coerce.h"
#include "sparse.h"

/* defined in ./factorizations.c : */
cholmod_sparse *dgC2cholmod(SEXP, int);
SEXP cholmod2dgC(      cholmod_sparse *, const char *, int);
cholmod_dense  *dge2cholmod(SEXP, int);
SEXP cholmod2dge(const cholmod_dense  *, const char *, int);

/* op(<dge>) * op(<dge>) */
static
SEXP dgeMatrix_prod(SEXP a, SEXP b, int atrans, int btrans)
{
	SEXP adim = PROTECT(GET_SLOT(a, Matrix_DimSym));
	int *padim = INTEGER(adim), am = padim[0], an = padim[1],
		rm = (atrans) ? an : am, rk = (atrans) ? am : an;
	UNPROTECT(1); /* adim */

	if (b == R_NilValue) {

		if ((Matrix_int_fast64_t) rm * rm > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");

		SEXP r = PROTECT(NEW_OBJECT_OF_CLASS("dpoMatrix"));

		SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
		int *prdim = INTEGER(rdim);
		prdim[0] = prdim[1] = rm;
		UNPROTECT(1); /* rdim */

		SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
			rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
		symmDN(rdimnames, adimnames, (atrans) ? 1 : 0);
		UNPROTECT(2); /* rdimnames, adimnames */

		if (rm > 0) {
		SEXP rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rm));
		double *prx = REAL(rx);
		Matrix_memset(prx, 0, (R_xlen_t) rm * rm, sizeof(double));
		if (rk > 0) {
			SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym));
			double *pax = REAL(ax), zero = 0.0, one = 1.0;
			F77_CALL(dsyrk)(
				"U", (atrans) ? "T" : "N", &rm, &rk,
				&one, pax, &am, &zero, prx, &rm FCONE FCONE);
			UNPROTECT(1); /* ax */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
		}

		UNPROTECT(1); /* r */
		return r;

	} else {

		SEXP bdim = PROTECT(GET_SLOT(b, Matrix_DimSym));
		int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
			rn = (btrans) ? bm : bn;
		UNPROTECT(1); /* bdim */

		if (rk != ((btrans) ? bn : bm))
			error(_("non-conformable arguments"));
		if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding %s"),
			      "R_XLEN_T_MAX");

		SEXP r = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix"));

		SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
		int *prdim = INTEGER(rdim);
		prdim[0] = rm;
		prdim[1] = rn;
		UNPROTECT(1); /* rdim */

		SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
			bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
			rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
		mmultDN(rdimnames,
		        adimnames, (atrans) ? 1 : 0,
		        bdimnames, (btrans) ? 0 : 1);
		UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

		if (rm > 0 && rn > 0) {
		SEXP rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
		double *prx = REAL(rx);
		if (rk == 0)
			Matrix_memset(prx, 0, (R_xlen_t) rm * rn, sizeof(double));
		else {
			SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
				bx = PROTECT(GET_SLOT(b, Matrix_xSym));
			double *pax = REAL(ax), *pbx = REAL(bx),
				zero = 0.0, one = 1.0;
			F77_CALL(dgemm)(
				(atrans) ? "T" : "N", (btrans) ? "T" : "N", &rm, &rn, &rk,
				&one, pax, &am, pbx, &bm, &zero, prx, &rm FCONE FCONE);
			UNPROTECT(2); /* bx, ax */
		}
		SET_SLOT(r, Matrix_xSym, rx);
		UNPROTECT(1); /* rx */
		}

		UNPROTECT(1); /* r */
		return r;

	}
}

/* <dsy> * op(<dge>)  or  op(<dge>) * <dsy> */
static
SEXP dsyMatrix_prod(SEXP a, SEXP b, int aleft, int btrans)
{
	SEXP adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
		bdim = PROTECT(GET_SLOT(b, Matrix_DimSym));
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn,
		rk = INTEGER(adim)[0];
	UNPROTECT(2); /* bdim, adim */

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix"));

	SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;
	UNPROTECT(1); /* rdim */

	SEXP adimnames = PROTECT(get_symmetrized_DimNames(a, -1)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	mmultDN(rdimnames, adimnames,      0, bdimnames, !btrans);
	else
	mmultDN(rdimnames, bdimnames, btrans, adimnames,       1);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	if (rm > 0 && rn > 0) {
	SEXP uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
		ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
		bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		zero = 0.0, one = 1.0;
	char ul = *CHAR(STRING_ELT(uplo, 0));
	if (!btrans)
		F77_CALL(dsymm)(
			(aleft) ? "L" : "R", &ul, &rm, &rn,
			&one, pax, &rk, pbx, &bm, &zero, prx, &rm FCONE FCONE);
	else {
		int i, d = (aleft) ? rn : rm,
			binc = (aleft) ? bm : 1, bincp = (aleft) ? 1 : bm,
			rinc = (aleft) ? 1 : rm, rincp = (aleft) ? rm : 1;
		for (i = 0; i < d; ++i) {
			F77_CALL(dsymv)(
				&ul, &rk, &one, pax, &rk, pbx, &binc, &zero, prx, &rinc FCONE);
			pbx += bincp;
			prx += rincp;
		}
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(4); /* rx, bx, ax, uplo */
	}

	UNPROTECT(1); /* r */
	return r;
}

/* op(<dtr>) * op(<dge>)  or  op(<dge>) * op(<dtr>) */
static
SEXP dtrMatrix_prod(SEXP a, SEXP b, int aleft, int atrans, int btrans,
                    int triangular)
{
	SEXP adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
		bdim = PROTECT(GET_SLOT(b, Matrix_DimSym));
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn,
		rk = INTEGER(adim)[0];
	UNPROTECT(2); /* bdim, adim */

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS(
		(triangular > 0) ? "dtrMatrix" : "dgeMatrix"));

	SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;
	UNPROTECT(1); /* rdim */

	SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	mmultDN(rdimnames, adimnames, atrans, bdimnames, !btrans);
	else
	mmultDN(rdimnames, bdimnames, btrans, adimnames, !atrans);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	SEXP uplo = PROTECT(GET_SLOT(a, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	if (triangular > 0 && ((atrans) ? ul == 'U' : ul != 'U')) {
		if (atrans) {
			UNPROTECT(1); /* uplo */
			PROTECT(uplo = mkString("L"));
		}
		SET_SLOT(r, Matrix_uploSym, uplo);
	}
	UNPROTECT(1); /* uplo */

	SEXP diag = PROTECT(GET_SLOT(a, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (triangular > 1 && di != 'N')
		SET_SLOT(r, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	if (rm > 0 && rn > 0) {
	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
		bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		one = 1.0;
	if (!btrans)
		Matrix_memcpy(prx, pbx, (R_xlen_t) rm * rn, sizeof(double));
	else {
		int i, j;
		R_xlen_t mn1s = (R_xlen_t) bm * bn - 1;
		for (j = 0; j < bm; ++j, pbx -= mn1s)
			for (i = 0; i < bn; ++i, pbx += bm)
				*(prx++) = *pbx;
		prx -= mn1s + 1;
	}
	F77_CALL(dtrmm)(
		(aleft) ? "L" : "R", &ul, (atrans) ? "T" : "N", &di, &rm, &rn,
		&one, pax, &rk, prx, &rm FCONE FCONE FCONE FCONE);
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(3); /* rx, bx, ax */
	}

	UNPROTECT(1); /* r */
	return r;
}

/* <dsp> * op(<dge>)  or  op(<dge>) * <dsp> */
static
SEXP dspMatrix_prod(SEXP a, SEXP b, int aleft, int btrans)
{
	SEXP adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
		bdim = PROTECT(GET_SLOT(b, Matrix_DimSym));
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn,
		rk = INTEGER(adim)[0];
	UNPROTECT(2); /* bdim, adim */

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS("dgeMatrix"));

	SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;
	UNPROTECT(1); /* rdim */

	SEXP adimnames = PROTECT(get_symmetrized_DimNames(a, -1)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	mmultDN(rdimnames, adimnames,      0, bdimnames, !btrans);
	else
	mmultDN(rdimnames, bdimnames, btrans, adimnames,       1);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	if (rm > 0 && rn > 0) {
	SEXP uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
		ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
		bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx),
		zero = 0.0, one = 1.0;
	char ul = *CHAR(STRING_ELT(uplo, 0));
	int i, d = (aleft) ? rn : rm,
		binc = (aleft == btrans) ? bm : 1, bincp = (aleft == btrans) ? 1 : bm,
		rinc = (aleft          ) ? 1 : rm, rincp = (aleft          ) ? rm : 1;
	for (i = 0; i < d; ++i) {
		F77_CALL(dspmv)(
			&ul, &rk, &one, pax, pbx, &binc, &zero, prx, &rinc FCONE);
		pbx += bincp;
		prx += rincp;
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(4); /* rx, bx, ax, uplo */
	}

	UNPROTECT(1); /* r */
	return r;
}

/* op(<dtp>) * op(<dge>)  or  op(<dge>) * op(<dtp>) */
static
SEXP dtpMatrix_prod(SEXP a, SEXP b, int aleft, int atrans, int btrans,
                    int triangular)
{
	SEXP adim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
		bdim = PROTECT(GET_SLOT(b, Matrix_DimSym));
	int *pbdim = INTEGER(bdim), bm = pbdim[0], bn = pbdim[1],
		rm = (btrans) ? bn : bm, rn = (btrans) ? bm : bn,
		rk = INTEGER(adim)[0];
	UNPROTECT(2); /* bdim, adim */

	if (rk != ((aleft == btrans) ? bn : bm))
		error(_("non-conformable arguments"));
	if ((Matrix_int_fast64_t) rm * rn > R_XLEN_T_MAX)
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");

	SEXP r = PROTECT(NEW_OBJECT_OF_CLASS(
		(triangular > 0) ? "dtrMatrix" : "dgeMatrix"));

	SEXP rdim = PROTECT(GET_SLOT(r, Matrix_DimSym));
	int *prdim = INTEGER(rdim);
	prdim[0] = rm;
	prdim[1] = rn;
	UNPROTECT(1); /* rdim */

	SEXP adimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
		bdimnames = PROTECT(GET_SLOT(b, Matrix_DimNamesSym)),
		rdimnames = PROTECT(GET_SLOT(r, Matrix_DimNamesSym));
	if (aleft)
	mmultDN(rdimnames, adimnames, atrans, bdimnames, !btrans);
	else
	mmultDN(rdimnames, bdimnames, btrans, adimnames, !atrans);
	UNPROTECT(3); /* rdimnames, bdimnames, adimnames */

	SEXP uplo = PROTECT(GET_SLOT(a, Matrix_uploSym));
	char ul = *CHAR(STRING_ELT(uplo, 0));
	if (triangular > 0 && ((atrans) ? ul == 'U' : ul != 'U')) {
		if (atrans) {
			UNPROTECT(1); /* uplo */
			PROTECT(uplo = mkString("L"));
		}
		SET_SLOT(r, Matrix_uploSym, uplo);
	}
	UNPROTECT(1); /* uplo */

	SEXP diag = PROTECT(GET_SLOT(a, Matrix_diagSym));
	char di = *CHAR(STRING_ELT(diag, 0));
	if (triangular > 1 && di != 'N')
		SET_SLOT(r, Matrix_diagSym, diag);
	UNPROTECT(1); /* diag */

	if (rm > 0 && rn > 0) {
	SEXP ax = PROTECT(GET_SLOT(a, Matrix_xSym)),
		bx = PROTECT(GET_SLOT(b, Matrix_xSym)),
		rx = PROTECT(allocVector(REALSXP, (R_xlen_t) rm * rn));
	double *pax = REAL(ax), *pbx = REAL(bx), *prx = REAL(rx);
	if (!btrans)
		Matrix_memcpy(prx, pbx, (R_xlen_t) rm * rn, sizeof(double));
	else {
		int i, j;
		R_xlen_t mn1s = (R_xlen_t) bm * bn - 1;
		for (j = 0; j < bm; ++j, pbx -= mn1s)
			for (i = 0; i < bn; ++i, pbx += bm)
				*(prx++) = *pbx;
		prx -= mn1s + 1;
	}
	int i, rinc = (aleft) ? 1 : rm, rincp = (aleft) ? rm : 1;
	for (i = 0; i < rn; ++i) {
		F77_CALL(dtpmv)(
			&ul, (aleft == atrans) ? "T" : "N", &di, &rk,
			pax, prx, &rinc FCONE FCONE FCONE);
		prx += rincp;
	}
	SET_SLOT(r, Matrix_xSym, rx);
	UNPROTECT(3); /* rx, bx, ax */
	}

	UNPROTECT(1); /* r */
	return r;
}

SEXP R_dense_prod(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans)
{
	static const char *valid[] = { VALID_DENSE, "" };
	int xtrans_ = LOGICAL(xtrans)[0], ytrans_ = LOGICAL(ytrans)[0],
		ivalid;

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (!IS_S4_OBJECT(x)) {
		int isM = isMatrix(x), trans = 0;
		if (!isM) {
			int k = INTEGER(GET_SLOT(y, Matrix_DimSym))[(ytrans_) ? 1 : 0];
			if (k != 1) {
				trans = 1;
				if (XLENGTH(x) != k)
					error(_("non-conformable arguments"));
			}
			xtrans_ = 0;
		}
		REPROTECT(x = matrix_as_dense(x, ".ge", '\0', '\0', trans, 0), xpid);
		if (!isM)
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (trans) ? 1 : 0, R_NilValue);
	}
	if (!IS_S4_OBJECT(y) && y != R_NilValue) {
		int isM = isMatrix(y), trans = 0;
		if (!isM) {
			int k = INTEGER(GET_SLOT(x, Matrix_DimSym))[(xtrans_) ? 0 : 1];
			if (k == 1)
				trans = 1;
			else if (XLENGTH(y) != k)
				error(_("non-conformable arguments"));
			ytrans_ = 0;
		}
		REPROTECT(y = matrix_as_dense(y, ".ge", '\0', '\0', trans, 0), ypid);
		if (!isM)
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (trans) ? 1 : 0, R_NilValue);
	}

	const char *xcl = NULL, *ycl = NULL;
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

	if (xcl[0] != 'd') {
		REPROTECT(x = dense_as_kind(x, xcl, 'd'), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
	}
	if (y != R_NilValue) {
	if (ycl[0] != 'd') {
		REPROTECT(y = dense_as_kind(y, ycl, 'd'), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
	}
	}

	if (y == R_NilValue) {
		REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
		x = dgeMatrix_prod(x, y, xtrans_, !xtrans_);
	} else if (xcl[1] == 'g' && ycl[1] == 'g') {
		x = dgeMatrix_prod(x, y, xtrans_, ytrans_);
	} else if (xcl[1] == 'g' || ycl[1] == 'g') {
		x = (xcl[1] == 'g')
			? ((ycl[1] == 's')
			   ? ((ycl[2] != 'p')
			      ? dsyMatrix_prod(y, x, 0, xtrans_)
			      : dspMatrix_prod(y, x, 0, xtrans_))
			   : ((ycl[2] != 'p')
			      ? dtrMatrix_prod(y, x, 0, ytrans_, xtrans_, 0)
			      : dtpMatrix_prod(y, x, 0, ytrans_, xtrans_, 0)))
			: ((xcl[1] == 's')
			   ? ((xcl[2] != 'p')
			      ? dsyMatrix_prod(x, y, 1, ytrans_)
			      : dspMatrix_prod(x, y, 1, ytrans_))
			   : ((xcl[2] != 'p')
			      ? dtrMatrix_prod(x, y, 1, xtrans_, ytrans_, 0)
			      : dtpMatrix_prod(x, y, 1, xtrans_, ytrans_, 0)));
	} else if (xcl[1] == 's' && ycl[1] == 's') {
		if (xcl[2] == 'p' && ycl[2] == 'p') {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dspMatrix_prod(x, y, 1, ytrans_);
		} else if (xcl[2] == 'p') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = dsyMatrix_prod(y, x, 0, xtrans_);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dsyMatrix_prod(x, y, 1, ytrans_);
		}
	} else if (xcl[1] == 's' || ycl[1] == 's') {
		if (xcl[1] == 's') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = (ycl[2] != 'p')
				? dtrMatrix_prod(y, x, 0, ytrans_, 0, 0)
				: dtpMatrix_prod(y, x, 0, ytrans_, 0, 0);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = (xcl[2] != 'p')
				? dtrMatrix_prod(x, y, 1, xtrans_, 0, 0)
				: dtpMatrix_prod(x, y, 1, xtrans_, 0, 0);
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
		if (xtrans_)
			xul = (xul == 'U') ? 'L' : 'U';
		if (ytrans_)
			yul = (yul == 'U') ? 'L' : 'U';
		int triangular = (xul != yul) ? 0 : ((xdi != ydi || xdi == 'N') ? 1 : 2);
		UNPROTECT(4); /* ydiag, xdiag, yuplo, xuplo */

		if (xcl[2] == 'p' && ycl[2] == 'p') {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dtpMatrix_prod(x, y, 1, xtrans_, ytrans_, triangular);
		} else if (xcl[2] == 'p') {
			REPROTECT(x = dense_as_general(x, xcl, 1), xpid);
			x = dtrMatrix_prod(y, x, 0, ytrans_, xtrans_, triangular);
		} else {
			REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
			x = dtrMatrix_prod(x, y, 1, xtrans_, ytrans_, triangular);
		}
	}

	UNPROTECT(2); /* y, x */
	return x;
}

/* boolean: op(op(<.gC>) * op(<.gC>)) */
/* numeric: op(op(<dgC>) * op(<dgC>)) */
static
SEXP dgCMatrix_dgCMatrix_prod(SEXP x, SEXP y, int xtrans, int ytrans,
                              int ztrans, int triangular, int boolean)
{
	PROTECT_INDEX zpid;
	SEXP z;
	char zcl[] = "..CMatrix";
	zcl[0] = (boolean) ? 'n' : 'd';
	if (y == R_NilValue) {
		zcl[1] = 's';
		cholmod_sparse *X = dgC2cholmod(x, !boolean);
		if (xtrans)
			X = cholmod_transpose(X, !boolean, &c);
		cholmod_sparse *Z = cholmod_aat(X, (int *) NULL, 0, !boolean, &c);
		if (xtrans)
			cholmod_free_sparse(&X, &c);
		Z->stype = (ztrans) ? -1 : 1;
		cholmod_sort(Z, &c);
		PROTECT_WITH_INDEX(z = cholmod2dgC(Z, zcl, !boolean), &zpid);
		cholmod_free_sparse(&Z, &c);
		SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
			zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
		symmDN(zdimnames, xdimnames, (xtrans) ? 1 : 0);
		UNPROTECT(2); /* zdimnames, xdimnames */
		if (ztrans) {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(z, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
	} else {
		zcl[1] = (triangular != 0) ? 't' : 'g';
		cholmod_sparse
			*X = dgC2cholmod(x, !boolean),
			*Y = dgC2cholmod(y, !boolean);
		if (((xtrans) ? X->nrow : X->ncol) != ((ytrans) ? Y->ncol : Y->nrow))
			error(_("non-conformable arguments"));
		if (xtrans)
			X = cholmod_transpose(X, !boolean, &c);
		if (ytrans)
			Y = cholmod_transpose(Y, !boolean, &c);
		cholmod_sparse *Z = cholmod_ssmult(X, Y, 0, !boolean, 1, &c);
		if (xtrans)
			cholmod_free_sparse(&X, &c);
		if (ytrans)
			cholmod_free_sparse(&Y, &c);
		PROTECT_WITH_INDEX(z = cholmod2dgC(Z, zcl, !boolean), &zpid);
		cholmod_free_sparse(&Z, &c);
		SEXP xdimnames = PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
			ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
			zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
		mmultDN(zdimnames,
		        xdimnames, (xtrans) ? 1 : 0,
		        ydimnames, (ytrans) ? 0 : 1);
		UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */
		if (triangular < 0) {
			SEXP uplo = PROTECT(mkString("L"));
			SET_SLOT(z, Matrix_uploSym, uplo);
			UNPROTECT(1); /* uplo */
		}
		if (triangular < -1 || triangular > 1)
			REPROTECT(z = sparse_diag_N2U(z, zcl), zpid);
	}
	if (ztrans)
		REPROTECT(z = sparse_transpose(z, zcl, 1), zpid);
	UNPROTECT(1); /* z */
	return z;
}

/* op(op(<d[gs]C>) * op(<dge>)) */
static
SEXP dgCMatrix_dgeMatrix_prod(SEXP x, SEXP y, int xtrans, int ytrans,
                              int ztrans, int triangular, int symmetric)
{
	SEXP z;
	char zcl[] = "d..Matrix";
	zcl[1] = (triangular) ? 't' : 'g';
	zcl[2] = (triangular) ? 'r' : 'e';
	double alpha[2] = { 1.0, 0.0 }, beta[2] = { 0.0, 0.0 };
	cholmod_sparse *X = dgC2cholmod(x, 1);
	cholmod_dense *Y = dge2cholmod(y, ytrans);
	X->stype = (symmetric) ? 1 : 0;
	int m = (int) ((xtrans) ? X->ncol : X->nrow), n = (int) Y->ncol;
	if ((Matrix_int_fast64_t) m * n > R_XLEN_T_MAX) {
		if (ytrans)
			R_Free(Y->x);
		error(_("attempt to allocate vector of length exceeding %s"),
		      "R_XLEN_T_MAX");
	}
	cholmod_dense *Z = (cholmod_dense *) R_alloc(1, sizeof(cholmod_dense));
	memset(Z, 0, sizeof(cholmod_dense));
	Z->nrow = (size_t) m;
	Z->ncol = (size_t) n;
	Z->d = Z->nrow;
	Z->nzmax = Z->nrow * Z->ncol;
	Z->xtype = CHOLMOD_REAL;
	Z->dtype = CHOLMOD_DOUBLE;
	if (ztrans) {
		Z->x = R_Calloc(Z->nzmax, double);
		cholmod_sdmult(X, xtrans, alpha, beta, Y, Z, &c);
		PROTECT(z = cholmod2dge(Z, zcl, ztrans));
		R_Free(Z->x);
	} else {
		PROTECT(z = NEW_OBJECT_OF_CLASS(zcl));
		SEXP zdim = PROTECT(GET_SLOT(z, Matrix_DimSym));
		INTEGER(zdim)[0] = m;
		INTEGER(zdim)[1] = n;
		UNPROTECT(1); /* zdim */
		SEXP zx = PROTECT(allocVector(REALSXP, (R_xlen_t) m * n));
		Z->x = REAL(zx);
		cholmod_sdmult(X, xtrans, alpha, beta, Y, Z, &c);
		SET_SLOT(z, Matrix_xSym, zx);
		UNPROTECT(1); /* zx */
	}
	if (ytrans)
		R_Free(Y->x);
	SEXP xdimnames = (symmetric)
		? PROTECT(get_symmetrized_DimNames(x, -1))
		: PROTECT(GET_SLOT(x, Matrix_DimNamesSym)),
		ydimnames = PROTECT(GET_SLOT(y, Matrix_DimNamesSym)),
		zdimnames = PROTECT(GET_SLOT(z, Matrix_DimNamesSym));
	if (ztrans)
	mmultDN(zdimnames,
			ydimnames, (ytrans) ? 0 : 1,
			xdimnames, (xtrans) ? 1 : 0);
	else
	mmultDN(zdimnames,
			xdimnames, (xtrans) ? 1 : 0,
			ydimnames, (ytrans) ? 0 : 1);
	UNPROTECT(3); /* zdimnames, ydimnames, xdimnames */
	if (triangular != 0 && ztrans == (triangular > 0)) {
		SEXP uplo = PROTECT(mkString("L"));
		SET_SLOT(z, Matrix_uploSym, uplo);
		UNPROTECT(1); /* uplo */
	}
	if (triangular < -1 || triangular > 1) {
		SEXP diag = PROTECT(mkString("U"));
		SET_SLOT(z, Matrix_diagSym, diag);
		UNPROTECT(1); /* diag */
	}
	UNPROTECT(1); /* z */
	return z;
}

SEXP R_sparse_prod(SEXP x, SEXP y, SEXP xtrans, SEXP ytrans, SEXP ztrans,
                   SEXP boolean)
{
	if (TYPEOF(boolean) != LGLSXP || LENGTH(boolean) < 1)
		error(_("invalid '%s' to %s()"), "boolean", __func__);

	static const char *valid[] = {
		VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE, VALID_DENSE, "" };
	int xtrans_ = LOGICAL(xtrans)[0], ytrans_ = LOGICAL(ytrans)[0],
		ztrans_ = LOGICAL(ztrans)[0], boolean_ = LOGICAL(boolean)[0],
		ivalid, triangular = 0;

	PROTECT_INDEX xpid, ypid;
	PROTECT_WITH_INDEX(x, &xpid);
	PROTECT_WITH_INDEX(y, &ypid);

	if (!IS_S4_OBJECT(x)) {
		int isM = isMatrix(x), trans = 0;
		if (!isM) {
			int k = INTEGER(GET_SLOT(y, Matrix_DimSym))[(ytrans_) ? 1 : 0];
			if (k != 1) {
				trans = 1;
				if (XLENGTH(x) != k)
					error(_("non-conformable arguments"));
			}
			xtrans_ = 0;
		}
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(x = matrix_as_dense( x, ".ge", '\0', '\0', trans, 0), xpid);
		else if (!xtrans_)
		REPROTECT(x = matrix_as_sparse(x, "ngC", '\0', '\0', trans   ), xpid);
		else
		REPROTECT(x = matrix_as_sparse(x, "ngR", '\0', '\0', trans   ), xpid);
		if (!isM)
			SET_VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym),
			               (trans) ? 1 : 0, R_NilValue);
	}
	if (!IS_S4_OBJECT(y) && y != R_NilValue) {
		int isM = isMatrix(y), trans = 0;
		if (!isM) {
			int k = INTEGER(GET_SLOT(x, Matrix_DimSym))[(xtrans_) ? 0 : 1];
			if (k == 1)
				trans = 1;
			else if (XLENGTH(y) != k)
				error(_("non-conformable arguments"));
			ytrans_ = 0;
		}
		if (boolean_ == NA_LOGICAL || !boolean_)
		REPROTECT(y = matrix_as_dense( y, ".ge", '\0', '\0', trans, 0), ypid);
		else if (!ytrans_)
		REPROTECT(y = matrix_as_sparse(y, "ngC", '\0', '\0', trans   ), ypid);
		else
		REPROTECT(y = matrix_as_sparse(y, "ngR", '\0', '\0', trans   ), ypid);
		if (!isM)
			SET_VECTOR_ELT(GET_SLOT(y, Matrix_DimNamesSym),
			               (trans) ? 1 : 0, R_NilValue);
	}

	const char *xcl = NULL, *ycl = NULL;
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

	if (xcl[2] != 'C' && xtrans_) {
		if (xcl[2] != 'R' && xcl[2] != 'T') {
			REPROTECT(x = dense_as_sparse(x, xcl, 'R'), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		if (xcl[1] != 's' || xcl[1] != 'T') {
			REPROTECT(x = sparse_transpose(x, xcl, 1), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
		xtrans_ = 0;
	}
	if (xcl[2] != 'C') {
		if (xcl[2] != 'R' && xcl[2] != 'T')
			REPROTECT(x = dense_as_sparse(x, xcl, 'C'), xpid);
		else
			REPROTECT(x = sparse_as_Csparse(x, xcl), xpid);
		xcl = valid[R_check_class_etc(x, valid)];
	}
	if (xcl[1] == 's')
		xtrans_ = 0;
	if (xcl[0] != (boolean_) ? 'n' : 'd') {
		if (boolean_)
			REPROTECT(x = sparse_drop0(x, xcl, 0.0), xpid);
		else {
			REPROTECT(x = sparse_as_kind(x, xcl, 'd'), xpid);
			xcl = valid[R_check_class_etc(x, valid)];
		}
	}

	if (y == R_NilValue) {
		REPROTECT(x = sparse_as_general(x, xcl), xpid);
		x = dgCMatrix_dgCMatrix_prod(
			x, y, xtrans_, !xtrans_, ztrans_, triangular, boolean_);
		UNPROTECT(2); /* y, x */
		return x;
	}

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
		if (xtrans_)
			xul = (xul == 'U') ? 'L' : 'U';
		if (ytrans_)
			yul = (yul == 'U') ? 'L' : 'U';
		triangular = (xul != yul) ? 0 : ((xdi != ydi || xdi == 'N') ? 1 : 2);
		if (xul != 'U')
			triangular = -triangular;
		UNPROTECT(4); /* ydiag, xdiag, yuplo, xuplo */
	}

	if (!boolean_ && ycl[2] != 'C' && ycl[2] != 'R' && ycl[2] != 'T') {
		if (ycl[0] != 'd') {
			REPROTECT(y = dense_as_kind(y, ycl, 'd'), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (xcl[1] == 't')
			REPROTECT(x = sparse_diag_U2N(x, xcl), xpid);
		REPROTECT(y = dense_as_general(y, ycl, 1), ypid);
		x = dgCMatrix_dgeMatrix_prod(
			x, y, xtrans_, ytrans_, ztrans_, triangular, xcl[1] == 's');
		UNPROTECT(2); /* y, x */
		return x;
	}

	if (ycl[2] != 'C' && ytrans_) {
		if (ycl[2] != 'R' && ycl[2] != 'T') {
			REPROTECT(y = dense_as_sparse(y, ycl, 'R'), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		if (ycl[1] != 's' || ycl[1] != 'T') {
			REPROTECT(y = sparse_transpose(y, ycl, 1), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
		ytrans_ = 0;
	}
	if (ycl[2] != 'C') {
		if (ycl[2] != 'R' && ycl[2] != 'T')
			REPROTECT(y = dense_as_sparse(y, ycl, 'C'), ypid);
		else
			REPROTECT(y = sparse_as_Csparse(y, ycl), ypid);
		ycl = valid[R_check_class_etc(y, valid)];
	}
	if (ycl[1] == 's')
		ytrans_ = 0;
	if (ycl[0] != (boolean_) ? 'n' : 'd') {
		if (boolean_)
			REPROTECT(y = sparse_drop0(y, ycl, 0.0), ypid);
		else {
			REPROTECT(y = sparse_as_kind(y, ycl, 'd'), ypid);
			ycl = valid[R_check_class_etc(y, valid)];
		}
	}

	REPROTECT(x = sparse_as_general(x, xcl), xpid);
	REPROTECT(y = sparse_as_general(y, ycl), ypid);
	x = dgCMatrix_dgCMatrix_prod(
		x, y, xtrans_, ytrans_, ztrans_, triangular, boolean_);
	UNPROTECT(2); /* y, x */
	return x;
}
