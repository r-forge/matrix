// CXSparse/Source/cs_randperm: generate a random permutation
// CXSparse, Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: LGPL-2.1+
#include "cs.h"
#include <R_ext/Random.h>

/* return a random permutation vector, the identity perm, or p = n-1:-1:0.
 * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
 * p = random permutation.  */
CS_INT *cs_randperm (CS_INT n, CS_INT seed)
{
    CS_INT *p, k, j, t ;
    if (seed == 0) return (NULL) ;      /* return p = NULL (identity) */
    p = cs_malloc (n, sizeof (CS_INT)) ;   /* allocate result */
    if (!p) return (NULL) ;             /* out of memory */
    for (k = 0 ; k < n ; k++) p [k] = n-k-1 ;
    if (seed == -1) return (p) ;        /* return reverse permutation */
    GetRNGstate();
    for (k = 0 ; k < n ; k++)
    {
        j = k + (CS_INT) (unif_rand() * (RAND_MAX + 1.0)) % (n-k) ;    /* j = rand integer in range k to n-1 */
        t = p [j] ;                     /* swap p[k] and p[j] */
        p [j] = p [k] ;
        p [k] = t ;
    }
    PutRNGstate();
    return (p) ;
}
