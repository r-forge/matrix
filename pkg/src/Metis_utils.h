#ifndef MATRIX_METIS_UTILS_H
#define MATRIX_METIS_UTILS_H

#include <Rdefines.h>
#include "metis.h"

void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     idxtype* perm, idxtype* iperm);

void col_metis_order(int j0, int j1, int i2,
		     const int Tp[], const int Ti[], int ans[]);

#endif

