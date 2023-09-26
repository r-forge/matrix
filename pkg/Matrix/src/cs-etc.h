#ifndef MATRIX_CS_ETC_H
#define MATRIX_CS_ETC_H

/* GOAL: move from CSparse to CXSparse to support complex LU and QR */

#include <Rinternals.h>
#include "cs.h"

cs *dgC2cs(SEXP, int);
SEXP cs2dgC(const cs *, const char *, int);

#endif /* MATRIX_CS_ETC_H */
