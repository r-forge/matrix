#ifndef MATRIX_CSCBLOCKED_H
#define MATRIX_CSCBLOCKED_H

#include "Mutils.h"

SEXP cscBlocked_validate(SEXP x);
void cscBlocked_mm(char side, char transa, int m, int n, int k,
		   double alpha, int nr, int nc,
		   const int ap[], const int ai[],
		   const double ax[],
		   double b[], int ldb,
		   double beta, double c[], int ldc);

#endif
