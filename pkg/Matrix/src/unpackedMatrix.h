#ifndef MATRIX_UNPAMATRIX_H
#define MATRIX_UNPAMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP unpackedMatrix_pack(SEXP from, SEXP tr_if_ge, SEXP up_if_ge);
SEXP matrix_pack(SEXP from, SEXP tr, SEXP up);

SEXP unpackedMatrix_is_symmetric(SEXP obj);
SEXP matrix_is_symmetric(SEXP obj);

SEXP unpackedMatrix_is_triangular(SEXP obj, SEXP upper);
SEXP matrix_is_triangular(SEXP obj, SEXP upper);

SEXP unpackedMatrix_is_diagonal(SEXP obj);
SEXP matrix_is_diagonal(SEXP obj);

SEXP unpackedMatrix_t(SEXP obj);

SEXP unpackedMatrix_diag_get(SEXP obj, SEXP nms);
SEXP unpackedMatrix_diag_set(SEXP obj, SEXP val);

#define UPM_ERROR_INVALID_CLASS(_CLASS_, _METHOD_)			\
    error(_("invalid class \"%s\" to 'unpackedMatrix_%s()'"),		\
	  _CLASS_, _METHOD_)
#define UPM_ERROR_INVALID_SLOT_TYPE(_SLOT_, _SEXPTYPE_, _METHOD_)	\
    error(_("'%s' slot of invalid type \"%s\" in 'unpackedMatrix_%s()'"), \
	  _SLOT_, type2char(_SEXPTYPE_), _METHOD_)

#endif
