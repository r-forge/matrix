#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stddef.h> /* size_t */
#include <Rinternals.h>
#include <Rversion.h>

#if R_VERSION < R_Version(4, 5, 0)
int ANY_ATTRIB(SEXP);
void CLEAR_ATTRIB(SEXP);
#endif

char *Matrix_sprintf(const char *, ...);

int equalString(SEXP, SEXP, R_xlen_t);
SEXP duplicateVector(SEXP);
SEXP allocZero(SEXPTYPE, R_xlen_t);
SEXP allocUnit(SEXPTYPE, R_xlen_t);
SEXP allocSeqInt(int, R_xlen_t);
void naToUnit(SEXP);

#endif /* MATRIX_UTILS_H */
