#ifndef MATRIX_H
#define MATRIX_H
#include <Rdefines.h>
#include <Rconfig.h>
#include "cholmod.h"

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

// Copied from ../../src/Mutils.h ----------------------------------------
#define MATRIX_valid_dense			\
        "dmatrix", "dgeMatrix",			\
	"lmatrix", "lgeMatrix",			\
	"nmatrix", "ngeMatrix",			\
	"zmatrix", "zgeMatrix"

#define MATRIX_valid_Csparse			\
 "dgCMatrix", "dsCMatrix", "dtCMatrix",		\
 "lgCMatrix", "lsCMatrix", "ltCMatrix",		\
 "ngCMatrix", "nsCMatrix", "ntCMatrix",		\
 "zgCMatrix", "zsCMatrix", "ztCMatrix"

#define MATRIX_valid_Tsparse			\
 "dgTMatrix", "dsTMatrix", "dtTMatrix",		\
 "lgTMatrix", "lsTMatrix", "ltTMatrix",		\
 "ngTMatrix", "nsTMatrix", "ntTMatrix",		\
 "zgTMatrix", "zsTMatrix", "ztTMatrix"

#define MATRIX_valid_Rsparse			\
 "dgRMatrix", "dsRMatrix", "dtRMatrix",		\
 "lgRMatrix", "lsRMatrix", "ltRMatrix",		\
 "ngRMatrix", "nsRMatrix", "ntRMatrix",		\
 "zgRMatrix", "zsRMatrix", "ztRMatrix"

#define MATRIX_valid_CHMfactor "dCHMsuper", "dCHMsimpl", "nCHMsuper", "nCHMsimpl"

CHM_SP M_as_cholmod_sparse (CHM_SP ans, SEXP x, Rboolean check_Udiag, Rboolean sort_in_place);
CHM_TR M_as_cholmod_triplet(CHM_TR ans, SEXP x, Rboolean check_Udiag);
CHM_DN M_as_cholmod_dense(CHM_DN ans, SEXP x);
CHM_DN M_numeric_as_chm_dense(CHM_DN ans, double *v, int nr, int nc);
CHM_FR M_as_cholmod_factor(CHM_FR ans, SEXP x);
double M_chm_factor_ldetL2(const_CHM_FR f);
CHM_FR M_chm_factor_update(CHM_FR f, CHM_SP A, double mult);

#define AS_CHM_DN(x) M_as_cholmod_dense((CHM_DN)alloca(sizeof(cholmod_dense)), x )
#define AS_CHM_FR(x) M_as_cholmod_factor((CHM_FR)alloca(sizeof(cholmod_factor)), x )

#define AS_CHM_SP(x) M_as_cholmod_sparse ((CHM_SP)alloca(sizeof(cholmod_sparse)), x, (Rboolean)TRUE, (Rboolean)FALSE)
#define AS_CHM_TR(x) M_as_cholmod_triplet((CHM_TR)alloca(sizeof(cholmod_triplet)),x, (Rboolean)TRUE)
/* the non-diagU2N-checking versions : */
#define AS_CHM_SP__(x) M_as_cholmod_sparse ((CHM_SP)alloca(sizeof(cholmod_sparse)), x, (Rboolean)FALSE, (Rboolean)FALSE)
#define AS_CHM_TR__(x) M_as_cholmod_triplet((CHM_TR)alloca(sizeof(cholmod_triplet)), x, (Rboolean)FALSE)

#define N_AS_CHM_DN(x,nr,nc) M_numeric_as_chm_dense((CHM_DN)alloca(sizeof(cholmod_dense)), x , nr, nc )

SEXP M_Csparse_diagU2N(SEXP x);
SEXP M_chm_factor_to_SEXP(CHM_FR f, int dofree);
SEXP M_chm_sparse_to_SEXP(CHM_SP a, int dofree, int uploT, int Rkind,
			  char *diag, SEXP dn);
SEXP M_chm_triplet_to_SEXP(CHM_TR a, int dofree, int uploT, int Rkind,
			   const char* diag, SEXP dn);

SEXP M_dpoMatrix_chol(SEXP x);

/* TODO:  Utilities for C level of model_matrix(*, sparse) */

#ifdef	__cplusplus
}
#endif

#endif /* MATRIX_H */
