#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stddef.h> /* size_t */
#include <Rinternals.h>

#if 0 /* MJ: used for a while but not clearly needed */
void *Matrix_memset(void *, int, R_xlen_t, size_t);
void *Matrix_memcpy(void *, const void *, R_xlen_t, size_t);
#endif
char *Matrix_sprintf(const char *, ...);

int equalString(SEXP, SEXP, R_xlen_t);
void naToOne(SEXP);

#endif /* MATRIX_UTILS_H */
