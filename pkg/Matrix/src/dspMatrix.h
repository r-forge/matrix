#ifndef MATRIX_SPMATRIX_H
#define MATRIX_SPMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

/* defined in factorizations.c : */
SEXP dspMatrix_trf_(SEXP, int);

double get_norm_dsp(SEXP obj, const char *typstr);
SEXP dspMatrix_norm(SEXP obj, SEXP type);
SEXP dspMatrix_rcond(SEXP obj);

/* MJ: no longer needed ... prefer more general packedMatrix_diag_[gs]et() */
#if 0
SEXP dspMatrix_getDiag(SEXP x);
SEXP lspMatrix_getDiag(SEXP x);
SEXP dspMatrix_setDiag(SEXP x, SEXP d);
SEXP lspMatrix_setDiag(SEXP x, SEXP d);
#endif

/* MJ: no longer needed ... prefer more general packedMatrix_unpack() */
#if 0
SEXP dspMatrix_as_dsyMatrix(SEXP from);
#endif
/* MJ */

#endif
