#ifndef LME4_UTILS_H
#define LME4_UTILS_H
#include <R_ext/Constants.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rversion.h>
#include <Rconfig.h>
#include "Matrix.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
#else
# define attr_hidden
#endif

extern
#include "Syms.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

extern cholmod_common c;

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}

/**
 * Create a named vector of type TYP
 *
 * @param TYP a vector SEXP type (e.g. REALSXP)
 * @param names names of list elements with null string appended
 *
 * @return pointer to a named vector of type TYP
 */
static R_INLINE SEXP
internal_make_named(int TYP, char *names[])
{
    SEXP ans, nms;
    int n, i;

    for (n = 0; strlen(names[n]) > 0; n++) {}
    ans = PROTECT(allocVector(TYP, n));
    nms = PROTECT(allocVector(STRSXP, n));
    for (i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static R_INLINE double*
internal_symmetrize(double *a, int nc)
{
    int i,j;
    for (i = 1; i < nc; i++)
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

#endif
