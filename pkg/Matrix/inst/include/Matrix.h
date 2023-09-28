#ifndef MATRIX_H
#define MATRIX_H

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#else
# define bool Rboolean
#endif

/* Copy and paste from Defn.h : */
/* 'alloca' is neither C99 nor POSIX */
#ifdef __GNUC__
/* This covers GNU, Clang and Intel compilers */
/* #undef needed in case some other header, e.g. malloc.h, already did this */
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
# ifdef HAVE_ALLOCA_H
/* This covers native compilers on Solaris and AIX */
#  include <alloca.h>
# endif
/* It might have been defined via some other standard header, e.g. stdlib.h */
# if !HAVE_DECL_ALLOCA
extern void *alloca(size_t);
# endif
#endif

/* Declarations of M_cholmod_*(), which are implemented in ./Matrix_stubs.c : */
#include "cholmod.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

#define MATRIX_VALID_ge_dense \
"dmatrix", "dgeMatrix", \
"lmatrix", "lgeMatrix", \
"nmatrix", "ngeMatrix", \
"zmatrix", "zgeMatrix"

#define MATRIX_VALID_Csparse \
"dgCMatrix", "dsCMatrix", "dtCMatrix", \
"lgCMatrix", "lsCMatrix", "ltCMatrix", \
"ngCMatrix", "nsCMatrix", "ntCMatrix", \
"zgCMatrix", "zsCMatrix", "ztCMatrix"

bool Matrix_isclass_ge_dense(SEXP);
bool Matrix_isclass_Csparse (SEXP);

CHM_FR M_as_cholmod_factor(CHM_FR, SEXP);
CHM_SP M_as_cholmod_sparse(CHM_SP, SEXP, Rboolean, Rboolean);
CHM_DN M_as_cholmod_dense (CHM_DN, SEXP);

SEXP M_chm_factor_to_SEXP(const_CHM_FR, int);
SEXP M_chm_sparse_to_SEXP(const_CHM_SP, int, int, int, const char *, SEXP);

double M_chm_factor_ldetL2(const_CHM_FR);
CHM_FR M_chm_factor_update(CHM_FR, const_CHM_SP, double);
CHM_DN M_numeric_as_chm_dense(CHM_DN, double *, int, int);

#define AS_CHM_FR(x) \
	M_as_cholmod_factor((CHM_FR) alloca(sizeof(cholmod_factor)), x)
#define AS_CHM_DN(x) \
	M_as_cholmod_dense ((CHM_DN) alloca(sizeof(cholmod_dense )), x)
#define AS_CHM_SP(x) \
	M_as_cholmod_sparse((CHM_SP) alloca(sizeof(cholmod_sparse)), x, \
	                    (Rboolean)  TRUE, (Rboolean) FALSE)
/* Non-(diag == "U")-checking version : */
#define AS_CHM_SP__(x) \
	M_as_cholmod_sparse((CHM_SP) alloca(sizeof(cholmod_sparse)), x, \
	                    (Rboolean) FALSE, (Rboolean) FALSE)

#define N_AS_CHM_DN(x, nr, nc) \
	M_numeric_as_chm_dense((CHM_DN) alloca(sizeof(cholmod_dense)), x , nr, nc)

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H */
