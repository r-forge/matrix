#ifndef MATRIX_CSCBLOCKED_H
#define MATRIX_CSCBLOCKED_H

#include "Mutils.h"
#include "R_ldl.h"
#include <R_ext/Lapack.h>

SEXP cscBlocked_validate(SEXP x);
void cscBlocked_mm(char side, char transa, int m, int n, int k,
		   double alpha, int nr, int nc,
		   const int ap[], const int ai[],
		   const double ax[],
		   const double b[], int ldb,
		   double beta, double c[], int ldc);
void cscb_tri(char upper, char unit, SEXP A, const int Parent[], SEXP AI);
void cscb_syrk(char uplo, char trans,
	       double alpha, SEXP A,
	       double beta, SEXP C);
int cscb_ldl(SEXP A, const int Parent[], SEXP L, SEXP D);
void cscb_trmm(char side, char uplo, char transa, char diag,
	       double alpha, SEXP A, double B[], int m, int n, int ldb);
void cscb_trcbm(char side, char uplo, char transa, char diag,
		double alpha, SEXP A, SEXP B);
void cscb_cscbm(char transa, char transb, double alpha, SEXP A, SEXP B,
	     double beta, SEXP C);
void cscb_mm(char side, char transa, int m, int n, int k,
	     double alpha, SEXP A,
	     const double B[], int ldb,
	     double beta, double C[], int ldc);

#endif
