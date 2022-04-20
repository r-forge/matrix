#ifndef MATRIX_PAMATRIX_H
#define MATRIX_PAMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP packedMatrix_unpack(SEXP from);

SEXP packedMatrix_is_symmetric(SEXP obj, SEXP checkDN);
SEXP packedMatrix_is_triangular(SEXP obj, SEXP upper);
SEXP packedMatrix_is_diagonal(SEXP obj);

SEXP packedMatrix_t(SEXP obj);

SEXP packedMatrix_diag_get(SEXP obj, SEXP nms);
SEXP packedMatrix_diag_set(SEXP obj, SEXP val);

SEXP packedMatrix_sub1(SEXP obj, SEXP index);
SEXP packedMatrix_sub1_mat(SEXP obj, SEXP index);
SEXP packedMatrix_sub2(SEXP obj, SEXP index1, SEXP index2, SEXP drop);

/* int i, j, n; R_xlen_t n2; */
#define PM_AR21_UP(i, j) (i) + ((j) * (((R_xlen_t) (j)) + 1)) / 2
#define PM_AR21_LO(i, j, n2) (i) + ((j) * ((n2) - (j) - 1)) / 2

#define PM_ERROR_INVALID_CLASS(_CLASS_, _METHOD_)			\
    error(_("invalid class \"%s\" to 'packedMatrix_%s()'"),		\
	  _CLASS_, _METHOD_)
#define PM_ERROR_INVALID_SLOT_TYPE(_SLOT_, _SEXPTYPE_, _METHOD_)	\
    error(_("'%s' slot of invalid type \"%s\" in 'packedMatrix_%s()'"), \
	  _SLOT_, type2char(_SEXPTYPE_), _METHOD_)

#endif
