#ifndef MATRIX_MUTILS_H
#define MATRIX_MUTILS_H

#include "Mdefines.h"
#include "Minlines.h"

#ifdef __cplusplus
extern "C" {
/* NB: this block must not include system or R headers */
#endif

SEXP newObject(const char *what);

void *Matrix_memset(void *dest,        int   ch, R_xlen_t length, size_t size);
void *Matrix_memcpy(void *dest, const void *src, R_xlen_t length, size_t size);

Rboolean DimNames_is_trivial(SEXP dn);
Rboolean DimNames_is_symmetric(SEXP dn);
SEXP R_DimNames_is_symmetric(SEXP dn);

void symDN(SEXP dest, SEXP src, int J);
SEXP R_symDN(SEXP dn);
SEXP get_symmetrized_DimNames(SEXP obj, int J);
void set_symmetrized_DimNames(SEXP obj, SEXP dn, int J);

void revDN(SEXP dest, SEXP src);
SEXP R_revDN(SEXP dn);
SEXP get_reversed_DimNames(SEXP obj);
void set_reversed_DimNames(SEXP obj, SEXP dn);

void set_DimNames(SEXP obj, SEXP dn);

SEXP get_factor(SEXP obj, const char *nm);
void set_factor(SEXP obj, const char *nm, SEXP val);
SEXP R_set_factor(SEXP obj, SEXP nm, SEXP val, SEXP warn);

int isPerm(const int *p, int n, int off);
int signPerm(const int *p, int n, int off);
void invertPerm(const int *p, int *ip, int n, int off, int ioff);
void asPerm(const int *p, int *ip, int m, int n, int off, int ioff);

SEXP R_isPerm(SEXP p, SEXP off);
SEXP R_signPerm(SEXP p, SEXP off);
SEXP R_invertPerm(SEXP p, SEXP off, SEXP ioff);
SEXP R_asPerm(SEXP p, SEXP off, SEXP ioff, SEXP n);

char typeToKind(SEXPTYPE type);
SEXPTYPE kindToType(char kind);
size_t kindToSize(char kind);

const char *Matrix_nonvirtual(SEXP obj, int strict);
SEXP R_Matrix_nonvirtual(SEXP obj, SEXP strict);
char Matrix_kind(SEXP obj, int i2d);
SEXP R_Matrix_kind(SEXP obj, SEXP i2d);
char Matrix_shape(SEXP obj);
SEXP R_Matrix_shape(SEXP obj);
char Matrix_repr(SEXP obj);
SEXP R_Matrix_repr(SEXP obj);

SEXP R_index_triangle(SEXP n, SEXP packed, SEXP upper, SEXP diag);
SEXP R_index_diagonal(SEXP n, SEXP packed, SEXP upper);

SEXP R_nnz(SEXP x, SEXP countNA, SEXP nnzmax);

void conjugate(SEXP x);
void zeroRe(SEXP x);
void zeroIm(SEXP x);
void naToOne(SEXP x);

Rboolean equal_character_vectors(SEXP s1, SEXP s2, int n);
SEXP append_to_named_list(SEXP x, const char *nm, SEXP val);

SEXP Matrix_expand_pointers(SEXP pP);
SEXP compressed_non_0_ij(SEXP x, SEXP colP);
SEXP m_encodeInd (SEXP ij,        SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP m_encodeInd2(SEXP i, SEXP j, SEXP di, SEXP orig_1, SEXP chk_bnds);
SEXP Mmatrix(SEXP args);

SEXP R_all0(SEXP x);
SEXP R_any0(SEXP x);


/* ================================================================== */
/* Defined elsewhere but used in a few places, hence "exported" here: */
/* ================================================================== */

#include "idz.h"

SEXP unpacked_force(SEXP, int, char, char);
SEXP packed_transpose(SEXP, int, char);

void validObject(SEXP, const char *);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_MUTILS_H */
