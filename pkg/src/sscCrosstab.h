#ifndef MATRIX_SSCCROSSTAB_H
#define MATRIX_SSCCROSSTAB_H

#include "Mutils.h"
#include "ldl.h"

SEXP sscCrosstab(SEXP flist, SEXP upper);
extern void ssc_metis_order(int n, const int Tp [], const int Ti [],
			    int Perm[], int iPerm[]);
SEXP sscCrosstab_groupedPerm(SEXP ctab);
/* Only used for testing - will be removed */
extern void ssclme_fill_LIp(int n, const int Parent[], int LIp[]);
SEXP sscCrosstab_L_LI_sizes(SEXP ctab, SEXP permexp);

#endif
