/* double (precision) TRiangular Matrices */

/* MJ: no longer needed ... nothing below */
#if 0
#include "dtrMatrix.h"
#endif /* MJ */

/* MJ: no longer needed ... prefer Cholesky_solve() */
#if 0

SEXP dtrMatrix_chol2inv(SEXP a)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dpoMatrix")),
	dim = PROTECT(GET_SLOT(a, Matrix_DimSym)),
	dimnames = PROTECT(GET_SLOT(a, Matrix_DimNamesSym)),
	uplo = PROTECT(GET_SLOT(a, Matrix_uploSym)),
	x;
    PROTECT_INDEX pid;
    PROTECT_WITH_INDEX(x = GET_SLOT(a, Matrix_xSym), &pid);
    REPROTECT(x = duplicate(x), pid);

    SET_SLOT(val, Matrix_DimSym, dim);
    SET_SLOT(val, Matrix_DimNamesSym, dimnames);
    SET_SLOT(val, Matrix_uploSym, uplo);
    SET_SLOT(val, Matrix_xSym, x);

    int *pdim = INTEGER(dim), info;
    double *px = REAL(x);
    const char *ul = CHAR(STRING_ELT(uplo, 0));
    
    F77_CALL(dpotri)(ul, pdim, px, pdim, &info FCONE);

    UNPROTECT(5);
    return val;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_diag_[gs]et() */
#if 0

#define GET_trMatrix_Diag(_C_TYPE_, _SEXPTYPE_, _SEXP_, _ONE_)		\
    int i, n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];			\
    SEXP x_x = GET_SLOT(x, Matrix_xSym);				\
									\
    SEXP ret = PROTECT(allocVector(_SEXPTYPE_, n));			\
    _C_TYPE_ *rv = _SEXP_(ret),						\
	     *xv = _SEXP_(x_x);						\
									\
    if ('U' == diag_P(x)[0]) {						\
	for (i = 0; i < n; i++) rv[i] = _ONE_;				\
    } else {								\
	for (i = 0; i < n; i++) rv[i] = xv[i * (n + 1)];		\
    }									\
    UNPROTECT(1);							\
    return ret


SEXP dtrMatrix_getDiag(SEXP x) {
    GET_trMatrix_Diag(double, REALSXP, REAL, 1.);
}

SEXP ltrMatrix_getDiag(SEXP x) {
    GET_trMatrix_Diag(  int, LGLSXP, LOGICAL, 1);
}

#define SET_trMatrix_Diag(_C_TYPE_, _SEXP_)				\
    if ('U' == diag_P(x)[0])						\
	error(_("cannot set diag() as long as 'diag = \"U\"'"));	\
			    /* careful to recycle RHS value: */		\
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];			\
    int l_d = LENGTH(d); Rboolean d_full = (l_d == n);			\
    if (!d_full && l_d != 1)						\
	error(_("replacement diagonal has wrong length"));		\
    SEXP ret = PROTECT(duplicate(x)),					\
	r_x = GET_SLOT(ret, Matrix_xSym);				\
    _C_TYPE_ *dv = _SEXP_(d),						\
	     *rv = _SEXP_(r_x);						\
									\
    if(d_full) for (int i = 0; i < n; i++)				\
	rv[i * (n + 1)] = dv[i];					\
    else for (int i = 0; i < n; i++)					\
	rv[i * (n + 1)] = *dv;						\
									\
    UNPROTECT(1);							\
    return ret

SEXP dtrMatrix_setDiag(SEXP x, SEXP d) {
    SET_trMatrix_Diag(double, REAL);
}

SEXP ltrMatrix_setDiag(SEXP x, SEXP d) {
    SET_trMatrix_Diag(  int, LOGICAL);
}

SEXP dtrMatrix_addDiag(SEXP x, SEXP d) {
    int n = INTEGER(GET_SLOT(x, Matrix_DimSym))[0];
    SEXP ret = PROTECT(duplicate(x)),
	r_x = GET_SLOT(ret, Matrix_xSym);
    double *dv = REAL(d), *rv = REAL(r_x);

    if ('U' == diag_P(x)[0])
	error(_("cannot add diag() as long as 'diag = \"U\"'"));
    for (int i = 0; i < n; i++) rv[i * (n + 1)] += dv[i];

    UNPROTECT(1);
    return ret;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general unpackedMatrix_pack() */
#if 0

SEXP dtrMatrix_as_dtpMatrix(SEXP from)
{
    SEXP val = PROTECT(NEW_OBJECT_OF_CLASS("dtpMatrix")),
	uplo = GET_SLOT(from, Matrix_uploSym),
	diag = GET_SLOT(from, Matrix_diagSym),
	dimP = GET_SLOT(from, Matrix_DimSym);
    int n = *INTEGER(dimP);

    SET_SLOT(val, Matrix_DimSym, duplicate(dimP));
    SET_SLOT(val, Matrix_diagSym, duplicate(diag));
    SET_SLOT(val, Matrix_uploSym, duplicate(uplo));
    ddense_pack(
	REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, (n*(n+1))/2)),
	REAL(GET_SLOT(from, Matrix_xSym)),
	n,
	*CHAR(STRING_ELT(uplo, 0)) == 'U' ? UPP : LOW,
	*CHAR(STRING_ELT(diag, 0)) == 'N' ? NUN : UNT);
    SET_SLOT(val, Matrix_DimNamesSym,
	     duplicate(GET_SLOT(from, Matrix_DimNamesSym)));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */

/* MJ: no longer needed ... prefer more general R_dense_as_matrix() */
#if 0

SEXP dtrMatrix_as_matrix(SEXP from, SEXP keep_dimnames)
{
    int *Dim = INTEGER(GET_SLOT(from, Matrix_DimSym));
    int m = Dim[0], n = Dim[1];
    SEXP val = PROTECT(allocMatrix(REALSXP, m, n));
    ddense_unpacked_make_triangular(Memcpy(REAL(val),
					   REAL(GET_SLOT(from, Matrix_xSym)),
					   m * n),
				    from);
    if(asLogical(keep_dimnames))
	setAttrib(val, R_DimNamesSymbol, GET_SLOT(from, Matrix_DimNamesSym));
    UNPROTECT(1);
    return val;
}

#endif /* MJ */
