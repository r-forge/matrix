#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stddef.h> /* size_t */
#include <Rinternals.h>

void *Matrix_memset(void *, int, R_xlen_t, size_t);
void *Matrix_memcpy(void *, const void *, R_xlen_t, size_t);
char *Matrix_sprintf(const char *, ...);

int equal_character_vectors(SEXP, SEXP, int);

void naToOne(SEXP);

#endif /* MATRIX_UTILS_H */
