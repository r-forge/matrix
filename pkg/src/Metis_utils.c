#include "Metis_utils.h"

void ssc_metis_order(int n, const int Tp [], const int Ti [],
		     idxtype* perm, idxtype* iperm)
{
    int  j, num_flag = 0, options_flag = 0;
    idxtype
	*xadj = Calloc(n+1, idxtype),
	*adj = Calloc(2 * (Tp[n] - n), idxtype);

				/* temporarily use perm to store lengths */
    memset(perm, 0, sizeof(idxtype) * n);
    for (j = 0; j < n; j++) {
	int ip, p2 = Tp[j+1];
	for (ip = Tp[j]; ip < p2; ip++) {
	    int i = Ti[ip];
	    if (i != j) {
		perm[i]++;
		perm[j]++;
	    }
	}
    }
    xadj[0] = 0;
    for (j = 0; j < n; j++) xadj[j+1] = xadj[j] + perm[j];
				/* temporarily use perm to store pointers */
    Memcpy(perm, xadj, n);
    for (j = 0; j < n; j++) {
	int ip, p2 = Tp[j+1];
	for (ip = Tp[j]; ip < p2; ip++) {
	    int i = Ti[ip];
	    if (i != j) {
		adj[perm[i]] = j;
		adj[perm[j]] = i;
		perm[i]++;
		perm[j]++;
	    }
	}
    }
    METIS_NodeND(&n, xadj, adj, &num_flag, &options_flag, perm, iperm);
    Free(xadj); Free(adj);
}

void col_metis_order(int j0, int j1, int i2, const int Tp[], const int Ti[],
		     int ans[])
{
    int j, nz = 0;		/* count off-diagonal pairs */
    for (j = j0; j < j1; j++) {	/* columns of interest */
	int ii, nr = 0, p2 = Tp[j + 1];
	for (ii = Tp[j]; ii < p2; ii++) {
	    int i = Ti[ii];
	    if (j1 <= i && i < i2) nr++; /* verify row index */
	}
	nz += (nr * (nr - 1))/2; /* add number of pairs of rows */
    }
    if (nz > 0) {		/* Form an ssc Matrix */
	int j, n = i2 - j1,	/* number of rows */
	    nnz = n + nz, pos;
	int *Ap = Calloc(n + 1, int),
	    *Ai = Calloc(nnz, int),
	    *Tj = Calloc(nnz, int),
	    *TTi = Calloc(nnz, int);
	double			/* needed for triplet_to_col */
	    *Ax = Calloc(nnz, double), /* FIXME: change triplet_to_col */
	    *Tx = Calloc(nnz, double); /* to check for null pointers */
	idxtype *perm = Calloc(n, idxtype),
	    *iperm = Calloc(n, idxtype);

	for (j = 0; j < n; j++) { /* diagonals */
	    TTi[j] = Tj[j] = j;
	    Tx[j] = 1.;
	}
	pos = n;
	for (j = j0; j < j1; j++) { /* create the pairs */
	    int ii, nr = 0, p2 = Tp[j + 1];
	    for (ii = Tp[j]; ii < p2; ii++) {
		int r1 = Ti[ii], i1;
		if (j1 <= r1 && r1 < i2) {
		    for (i1 = ii + 1; i1 < p2; i1++) {
			int r2 = Ti[i1];
			if (r2 < i2) {
			    TTi[pos] = r2 - j1;
			    Tj[pos] = r1 - j1;
			    Tx[pos] = 1.;
			    pos++;
			}
		    }
		}
	    }
	}
	triplet_to_col(n, n, nnz, TTi, Tj, Tx, Ap, Ai, Ax);
	ssc_metis_order(n, Ap, Ai, perm, iperm);
	for (j = j1; j < i2; j++) ans[j] = j1 + iperm[j - j1];
	Free(Tx); Free(Ax); Free(TTi); Free(Tj); Free(Ai); Free(Ap);
	Free(perm); Free(iperm);
    }
}
