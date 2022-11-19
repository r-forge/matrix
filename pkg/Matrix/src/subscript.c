/* unfinished and not-yet-used, along with ../R/subscript.R */
#include "subscript.h"

static SEXP unpackedMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

static SEXP   packedMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

static SEXP  CsparseMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

static SEXP  RsparseMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

static SEXP  TsparseMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

static SEXP diagonalMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

static SEXP      indMatrix_subscript_1ary(SEXP x, SEXP i, const char *cl)
{

}

SEXP R_subscript_1ary(SEXP x, SEXP i)
{
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(obj, valid);
    const char *cl = valid[ivalid];
    if (ivalid < 0)
	ERROR_INVALID_CLASS(X, "R_subscript_1ary");

    //SEXPTYPE type = (cl[2] == 'd') ? LGLSXP : kind2type(cl[0]);
    //R_xlen_t len = XLENGTH(i);
    //SEXP res = PROTECT(allocVector(type, len));
    
    switch (cl[2]) {
    case 'e':
    case 'y':
    case 'r':
	return unpackedMatrix_subscript_1ary(x, i, cl);
    case 'p':
	return   packedMatrix_subscript_1ary(x, i, cl);
    case 'C':
	return  CsparseMatrix_subscript_1ary(x, i, cl);
    case 'R':
	return  RsparseMatrix_subscript_1ary(x, i, cl);
    case 'T':
        return  TsparseMatrix_subscript_1ary(x, i, cl);
    case 'i':
	return diagonalMatrix_subscript_1ary(x, i, cl);
    default:
	return      indMatrix_subscript_1ary(x, i, cl);
    }
}

static SEXP unpackedMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

static SEXP   packedMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

static SEXP  CsparseMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

static SEXP  RsparseMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

static SEXP  TsparseMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

static SEXP diagonalMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

static SEXP      indMatrix_subscript_1ary_mat(SEXP x, SEXP i, const char *cl)
{

}

SEXP R_subscript_1ary_mat(SEXP x, SEXP i)
{
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(obj, valid);
    const char *cl = valid[ivalid];
    if (ivalid < 0)
	ERROR_INVALID_CLASS(x, "R_subscript_1ary_mat");

    //SEXPTYPE type = (cl[2] == 'd') ? LGLSXP : kind2type(cl[0]);
    //R_xlen_t len = XLENGTH(i) / 2;
    //SEXP res = PROTECT(allocVector(type, len));
    
    switch (cl[2]) {
    case 'e':
    case 'y':
    case 'r':
	return unpackedMatrix_subscript_1ary_mat(x, i, cl);
    case 'p':
	return   packedMatrix_subscript_1ary_mat(x, i, cl);
    case 'C':
	return  CsparseMatrix_subscript_1ary_mat(x, i, cl);
    case 'R':
	return  RsparseMatrix_subscript_1ary_mat(x, i, cl);
    case 'T':
        return  TsparseMatrix_subscript_1ary_mat(x, i, cl);
    case 'i':
	return diagonalMatrix_subscript_1ary_mat(x, i, cl);
    default:
	return      indMatrix_subscript_1ary_mat(x, i, cl);
    }
}

static SEXP unpackedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP   packedMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP  CsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP  RsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP  TsparseMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP diagonalMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

static SEXP      indMatrix_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop,
					  const char *cl)
{

}

SEXP R_subscript_2ary(SEXP x, SEXP i, SEXP j, SEXP drop)
{
    static const char *valid[] = { VALID_NONVIRTUAL, "" };
    int ivalid = R_check_class_etc(obj, valid);
    const char *cl = valid[ivalid];
    if (ivalid < 0)
	ERROR_INVALID_CLASS(x, "R_subscript_2ary");

    switch (cl[2]) {
    case 'e':
    case 'y':
    case 'r':
	return unpackedMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'p':
	return   packedMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'C':
	return  CsparseMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'R':
	return  RsparseMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'T':
	return  TsparseMatrix_subscript_2ary(x, i, j, drop, cl);
    case 'i':
	return diagonalMatrix_subscript_2ary(x, i, j, drop, cl);
    default:
	return      indMatrix_subscript_2ary(x, i, j, drop, cl);
    }
}
