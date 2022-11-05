#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h> /* C99 for int64_t */
#include <ctype.h>
#include <R.h> /* includes Rconfig.h, R_ext/RS.h */
#include <Rinternals.h>
#include <Rversion.h>
#include "Mdefines.h"

SEXP NEW_OBJECT_OF_CLASS(const char* what);

Rboolean DimNames_is_trivial(SEXP dn);
Rboolean DimNames_is_symmetric(SEXP dn);
SEXP R_DimNames_is_symmetric(SEXP dn);
    
void symmDN(SEXP dest, SEXP src, int J);
SEXP R_symmDN(SEXP dn);
SEXP get_symmetrized_DimNames(SEXP obj, int J);
void set_symmetrized_DimNames(SEXP obj, SEXP dn, int J);

void revDN(SEXP dest, SEXP src);
SEXP R_revDN(SEXP dn);
SEXP get_reversed_DimNames(SEXP obj);
void set_reversed_DimNames(SEXP obj, SEXP dn);

void set_DimNames(SEXP obj, SEXP dn);

SEXP get_factor(SEXP obj, const char *nm);
void set_factor(SEXP obj, const char *nm, SEXP val);
SEXP R_set_factor(SEXP obj, SEXP val, SEXP nm, SEXP warn);
SEXP R_empty_factors(SEXP obj, SEXP warn);

#define PACK(_PREFIX_, _CTYPE_)						\
void _PREFIX_ ## dense_pack(_CTYPE_ *dest, const _CTYPE_ *src, int n,	\
			    char uplo, char diag)
PACK(d, double);
PACK(i, int);
PACK(z, Rcomplex);
#undef PACK

#define UNPACK(_PREFIX_, _CTYPE_)					\
void _PREFIX_ ## dense_unpack(_CTYPE_ *dest, const _CTYPE_ *src, int n, \
			      char uplo, char diag)
UNPACK(d, double);
UNPACK(i, int);
UNPACK(z, Rcomplex);
#undef UNPACK

#define UNPACKED_MAKE_SYMMETRIC(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_make_symmetric(_CTYPE_ *x, int n, char uplo)
UNPACKED_MAKE_SYMMETRIC(d, double);
UNPACKED_MAKE_SYMMETRIC(i, int);
UNPACKED_MAKE_SYMMETRIC(z, Rcomplex);
#undef UNPACKED_MAKE_SYMMETRIC

#define UNPACKED_MAKE_TRIANGULAR(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_make_triangular(_CTYPE_ *x, int m, int n, \
						char uplo, char diag)
UNPACKED_MAKE_TRIANGULAR(d, double);
UNPACKED_MAKE_TRIANGULAR(i, int);
UNPACKED_MAKE_TRIANGULAR(z, Rcomplex);
#undef UNPACKED_MAKE_TRIANGULAR

#define UNPACKED_MAKE_BANDED(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_unpacked_make_banded(_CTYPE_ *x,			\
					    int m, int n, int a, int b,	\
					    char diag)
UNPACKED_MAKE_BANDED(d, double);
UNPACKED_MAKE_BANDED(i, int);
UNPACKED_MAKE_BANDED(z, Rcomplex);
#undef UNPACKED_MAKE_BANDED

#define PACKED_MAKE_BANDED(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_make_banded(_CTYPE_ *x,			\
					  int n, int a, int b,		\
					  char uplo, char diag)
PACKED_MAKE_BANDED(d, double);
PACKED_MAKE_BANDED(i, int);
PACKED_MAKE_BANDED(z, Rcomplex);
#undef PACKED_MAKE_BANDED

#define UNPACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_)			\
void _PREFIX_ ## dense_unpacked_copy_diagonal(_CTYPE_ *dest,	        \
					      const _CTYPE_ *src,	\
					      int n, R_xlen_t len,	\
					      char uplo, char diag)
UNPACKED_COPY_DIAGONAL(d, double);
UNPACKED_COPY_DIAGONAL(i, int);
UNPACKED_COPY_DIAGONAL(z, Rcomplex);
#undef UNPACKED_COPY_DIAGONAL

#define PACKED_COPY_DIAGONAL(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_copy_diagonal(_CTYPE_ *dest,		\
					    const _CTYPE_ *src,		\
					    int n, R_xlen_t len,	\
					    char uplo_dest,		\
					    char uplo_src,		\
					    char diag)
PACKED_COPY_DIAGONAL(d, double);
PACKED_COPY_DIAGONAL(i, int);
PACKED_COPY_DIAGONAL(z, Rcomplex);
#undef PACKED_COPY_DIAGONAL

#define UNPACKED_IS_SYMMETRIC(_PREFIX_, _CTYPE_)			\
Rboolean _PREFIX_ ## dense_unpacked_is_symmetric(const _CTYPE_ *x, int n)
UNPACKED_IS_SYMMETRIC(d, double);
UNPACKED_IS_SYMMETRIC(l, int);
UNPACKED_IS_SYMMETRIC(n, int);
UNPACKED_IS_SYMMETRIC(i, int);
UNPACKED_IS_SYMMETRIC(z, Rcomplex);
#undef UNPACKED_IS_SYMMETRIC
    
#define UNPACKED_IS_TRIANGULAR(_PREFIX_, _CTYPE_)			\
Rboolean _PREFIX_ ## dense_unpacked_is_triangular(const _CTYPE_ *x, int n, \
						  char uplo)
UNPACKED_IS_TRIANGULAR(d, double);
UNPACKED_IS_TRIANGULAR(i, int);
UNPACKED_IS_TRIANGULAR(z, Rcomplex);
#undef UNPACKED_IS_TRIANGULAR

#define UNPACKED_IS_DIAGONAL(_PREFIX_, _CTYPE_)				\
Rboolean _PREFIX_ ## dense_unpacked_is_diagonal(const _CTYPE_ *x, int n)
UNPACKED_IS_DIAGONAL(d, double);
UNPACKED_IS_DIAGONAL(i, int);
UNPACKED_IS_DIAGONAL(z, Rcomplex);
#undef UNPACKED_IS_DIAGONAL
    
#define PACKED_IS_DIAGONAL(_PREFIX_, _CTYPE_)			       \
Rboolean _PREFIX_ ## dense_packed_is_diagonal(const _CTYPE_ *x, int n, \
					      char uplo)
PACKED_IS_DIAGONAL(d, double);
PACKED_IS_DIAGONAL(i, int);
PACKED_IS_DIAGONAL(z, Rcomplex);
#undef PACKED_IS_DIAGONAL

#define PACKED_TRANSPOSE(_PREFIX_, _CTYPE_)				\
void _PREFIX_ ## dense_packed_transpose(_CTYPE_ *dest, const _CTYPE_ *src, \
					int n, char uplo)
PACKED_TRANSPOSE(d, double);
PACKED_TRANSPOSE(i, int);
PACKED_TRANSPOSE(z, Rcomplex);
#undef PACKED_TRANSPOSE

SEXP packed_transpose(SEXP x, int n, char uplo);
SEXP unpacked_force(SEXP x, int n, char uplo, char diag);
    
char type2kind(SEXPTYPE type);
SEXPTYPE kind2type(char kind);
size_t kind2size(char kind);

char Matrix_kind(SEXP obj, int i2d);
SEXP R_Matrix_kind(SEXP obj, SEXP i2d);
char Matrix_shape(SEXP obj);
SEXP R_Matrix_shape(SEXP obj);
char Matrix_repr(SEXP obj);
SEXP R_Matrix_repr(SEXP obj);

SEXP matrix_as_dense(SEXP from, const char *code, char uplo, char diag,
		     int new, int transpose_if_vector);
SEXP R_matrix_as_dense(SEXP from, SEXP code, SEXP uplo, SEXP diag);

SEXP dense_as_general(SEXP from, char kind, int new, int transpose_if_vector);
SEXP R_dense_as_general(SEXP from, SEXP kind);

SEXP R_index_triangle(SEXP n_, SEXP upper_, SEXP diag_, SEXP packed_);
SEXP R_index_diagonal(SEXP n_, SEXP upper_,             SEXP packed_);

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax);

void conjugate(SEXP x);
void zeroRe(SEXP x);
void zeroIm(SEXP x);
void na2one(SEXP x);

SEXP v2spV(SEXP from);

Rboolean equal_string_vectors(SEXP s1, SEXP s2, int n);
R_xlen_t strmatch(const char *nm, SEXP s);
SEXP append_to_named_list(SEXP x, const char *nm, SEXP val);

char La_norm_type(const char *typstr);
char La_rcond_type(const char *typstr);
SEXP as_det_obj(double mod, int log, int sign);

/* MJ: no longer needed ... prefer more general (un)?packedMatrix_diag_[gs]et() */
#if 0
void d_packed_getDiag(double *dest, SEXP x, int n);
void l_packed_getDiag(   int *dest, SEXP x, int n);
SEXP d_packed_setDiag(double *diag, int l_d, SEXP x, int n);
SEXP l_packed_setDiag(   int *diag, int l_d, SEXP x, int n);
void tr_d_packed_getDiag(double *dest, SEXP x, int n);
void tr_l_packed_getDiag(   int *dest, SEXP x, int n);
SEXP tr_d_packed_setDiag(double *diag, int l_d, SEXP x, int n);
SEXP tr_l_packed_setDiag(   int *diag, int l_d, SEXP x, int n);
/* were unused, not replaced: */
SEXP d_packed_addDiag(double *diag, int l_d, SEXP x, int n); 
SEXP tr_d_packed_addDiag(double *diag, int l_d, SEXP x, int n);
#endif /* MJ */

#if 0 /* unused */
double get_double_by_name(SEXP obj, char *nm);
SEXP set_double_by_name(SEXP obj, double val, char *nm);
SEXP dgCMatrix_set_Dim(SEXP x, int nrow);
SEXP new_dgeMatrix(int nrow, int ncol);
#endif /* unused */

SEXP Matrix_expand_pointers(SEXP pP);
SEXP m_encodeInd (SEXP ij,        SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP Mmatrix(SEXP args);

SEXP R_rbind2_vector(SEXP a, SEXP b);
SEXP R_all0(SEXP x);
SEXP R_any0(SEXP x);

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 * NOTE:  GET_SLOT(x, what)        :== R_do_slot       (x, what)
 * ----   SET_SLOT(x, what, value) :== R_do_slot_assign(x, what, value)
 * and the R_do_slot* are in src/main/attrib.c
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, R_xlen_t length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}

/**
 * Expand compressed pointers in the array mp into a full set of indices
 * in the array mj.
 *
 * @param ncol number of columns (or rows)
 * @param mp column pointer vector of length ncol + 1
 * @param mj vector of length mp[ncol] to hold the result
 *
 * @return mj
 */
static R_INLINE
int* expand_cmprPt(int ncol, const int mp[], int mj[])
{
    int j;
    for (j = 0; j < ncol; j++) {
	int j2 = mp[j+1], jj;
	for (jj = mp[j]; jj < j2; jj++) mj[jj] = j;
    }
    return mj;
}

/**
 * Check if slot(obj, "x") contains any NA (or NaN).
 *
 * @param obj   a 'Matrix' object with a (double precision) 'x' slot.
 *
 * @return Rboolean :== any(is.na(slot(obj, "x") )
 */
static R_INLINE
Rboolean any_NA_in_x(SEXP obj)
{
    SEXP x = PROTECT(GET_SLOT(obj, Matrix_xSym));
    R_xlen_t i, n = XLENGTH(x);
    double *px = REAL(x);
    Rboolean res = FALSE;
    for (i = 0; i < n; ++i) {
	if (ISNAN(px[i])) {
	    res = TRUE;
	    break;
	}
    }
    UNPROTECT(1);
    return res;
}

/** Inverse Permutation
 * C version of   .inv.perm.R <- function(p) { p[p] <- seq_along(p) ; p }
 */
static R_INLINE
SEXP inv_permutation(SEXP p_, SEXP zero_p, SEXP zero_res)
{
    int np = 1;
    if(!isInteger(p_)) {p_ = PROTECT(coerceVector(p_, INTSXP)); np++; }
    int *p = INTEGER(p_), n = LENGTH(p_);
    SEXP val = PROTECT(allocVector(INTSXP, n));
    int *v = INTEGER(val), p_0 = asLogical(zero_p), r_0 = asLogical(zero_res);
    if(!p_0) v--; // ==> use 1-based indices
    // shorter (but not 100% sure if ok: is LHS always eval'ed *before* RHS ?) :
    // for(int i=0; i < n; ) v[p[i]] = ++i;
    for(int i=0; i < n; ) {
	int j = p[i]; v[j] = (r_0) ? i++ : ++i;
    }
    UNPROTECT(np);
    return val;
}

/**
 * Return the 0-based index of a string match in a vector of strings
 * terminated by an empty string.  Returns -1 for no match.
 * Is  __cheap__ :  __not__ looking at superclasses --> better use  R_check_class_etc(obj, *)
 *
 * @param x string to match
 * @param valid vector of possible matches terminated by an empty string
 *
 * @return index of match or -1 for no match
 */
static R_INLINE
int Matrix_check_class_(char *x, const char **valid)
{
    int ans = 0;
    while (strlen(valid[ans]) > 0)
	if (strcmp(x, valid[ans]) == 0)
	    return ans;
	else
	    ++ans;
    return -1;
}

static R_INLINE
int Matrix_check_class(SEXP x, const char **valid)
{
    return Matrix_check_class_((char *) class_P(x), valid);
}

/** Accessing  *sparseVectors :  fast (and recycling)  v[i] for v = ?sparseVector:
 * -> ./sparseVector.c  -> ./t_sparseVector.c :
 */
// Type_ans sparseVector_sub(int64_t i, int nnz_v, int* v_i, Type_ans* v_x, int len_v):

/* Define all of
 *  dsparseVector_sub(....)
 *  isparseVector_sub(....)
 *  lsparseVector_sub(....)
 *  nsparseVector_sub(....)
 *  zsparseVector_sub(....)
 */
#define _dspV_
#include "t_sparseVector.c"

#define _ispV_
#include "t_sparseVector.c"

#define _lspV_
#include "t_sparseVector.c"

#define _nspV_
#include "t_sparseVector.c"

#define _zspV_
#include "t_sparseVector.c"


#ifdef __cplusplus
}
#endif

#endif /* MATRIX_UTILS_H */
