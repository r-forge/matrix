#ifndef MATRIX_COERCE_H
#define MATRIX_COERCE_H

#include "Mutils.h"

SEXP MJ_matrix_as_dense(SEXP from, const char *zzz, char ul, char di,
                     int transpose_if_vector, int new);

SEXP MJ_R_matrix_as_dense(SEXP from, SEXP class, SEXP uplo, SEXP diag);

SEXP MJ_sparse_as_dense(SEXP from, const char *class, int packed);

SEXP MJ_R_sparse_as_dense(SEXP from, SEXP packed);

SEXP MJ_diagonal_as_dense(SEXP from, const char *class,
                       char shape, int packed, char ul);

SEXP MJ_R_diagonal_as_dense(SEXP from, SEXP shape, SEXP packed, SEXP uplo);

SEXP MJ_index_as_dense(SEXP from, const char *class, char kind);

SEXP MJ_R_index_as_dense(SEXP from, SEXP kind);

SEXP MJ_matrix_as_sparse(SEXP from, const char *zzz, char ul, char di,
                      int transpose_if_vector);

SEXP MJ_R_matrix_as_sparse(SEXP from, SEXP class, SEXP uplo, SEXP diag);

SEXP MJ_dense_as_sparse(SEXP from, const char *class, char repr);

SEXP MJ_R_dense_as_sparse(SEXP from, SEXP repr);

SEXP MJ_diagonal_as_sparse(SEXP from, const char *class,
                        char shape, char repr, char ul);

SEXP MJ_R_diagonal_as_sparse(SEXP from, SEXP shape, SEXP repr, SEXP uplo);

SEXP MJ_index_as_sparse(SEXP from, const char *class, char kind, char repr);

SEXP MJ_R_index_as_sparse(SEXP from, SEXP kind, SEXP repr);

SEXP MJ_dense_as_kind(SEXP from, const char *class, char kind);

SEXP MJ_R_dense_as_kind(SEXP from, SEXP kind);

SEXP MJ_sparse_as_kind(SEXP from, const char *class, char kind);

SEXP MJ_R_sparse_as_kind(SEXP from, SEXP kind);

SEXP MJ_diagonal_as_kind(SEXP from, const char *class, char kind);

SEXP MJ_R_diagonal_as_kind(SEXP from, SEXP kind);

SEXP MJ_index_as_kind(SEXP from, const char *class, char kind);

SEXP MJ_R_index_as_kind(SEXP from, SEXP kind);

SEXP MJ_dense_as_general(SEXP from, const char *class, int new);

SEXP MJ_R_dense_as_general(SEXP from);

SEXP MJ_sparse_as_general(SEXP from, const char *class);

SEXP MJ_R_sparse_as_general(SEXP from);

SEXP MJ_dense_as_unpacked(SEXP from, const char *class);

SEXP MJ_R_dense_as_unpacked(SEXP from);

SEXP MJ_dense_as_packed(SEXP from, const char *class, char ul, char di);

SEXP MJ_R_dense_as_packed(SEXP from, SEXP uplo, SEXP diag);

SEXP MJ_sparse_as_Csparse(SEXP from, const char *class);

SEXP MJ_R_sparse_as_Csparse(SEXP from);

SEXP MJ_sparse_as_Rsparse(SEXP from, const char *class);

SEXP MJ_R_sparse_as_Rsparse(SEXP from);

SEXP MJ_sparse_as_Tsparse(SEXP from, const char *class);

SEXP MJ_R_sparse_as_Tsparse(SEXP from);

SEXP MJ_R_Matrix_as_vector(SEXP from);

SEXP MJ_R_Matrix_as_matrix(SEXP from);

SEXP MJ_R_Matrix_as_unpacked(SEXP from);

SEXP MJ_R_Matrix_as_packed(SEXP from);

SEXP MJ_R_Matrix_as_Csparse(SEXP from);

SEXP MJ_R_Matrix_as_Rsparse(SEXP from);

SEXP MJ_R_Matrix_as_Tsparse(SEXP from);

SEXP MJ_R_Matrix_as_kind(SEXP from, SEXP kind, SEXP sparse);

SEXP MJ_R_Matrix_as_general(SEXP from, SEXP kind);

#endif
