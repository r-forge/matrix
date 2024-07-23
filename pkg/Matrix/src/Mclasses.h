#ifndef MATRIX_MCLASSES_H
#define MATRIX_MCLASSES_H

#define VALID_UNPACKED \
"ngeMatrix", "lgeMatrix", "igeMatrix", "dgeMatrix", "zgeMatrix", \
"nsyMatrix", "lsyMatrix", "isyMatrix", "dsyMatrix", "zsyMatrix", \
                                       "dpoMatrix", "zpoMatrix", \
                                       "corMatrix",              \
"ntrMatrix", "ltrMatrix", "itrMatrix", "dtrMatrix", "ztrMatrix"

#define VALID_PACKED \
"nspMatrix", "lspMatrix", "ispMatrix", "dspMatrix", "zspMatrix", \
                                       "dppMatrix", "zppMatrix", \
                                       "copMatrix",              \
"ntpMatrix", "ltpMatrix", "itpMatrix", "dtpMatrix", "ztpMatrix"

#define VALID_CSPARSE \
"ngCMatrix", "lgCMatrix", "igCMatrix", "dgCMatrix", "zgCMatrix", \
"nsCMatrix", "lsCMatrix", "isCMatrix", "dsCMatrix", "zsCMatrix", \
                                       "dpCMatrix", "zpCMatrix", \
"ntCMatrix", "ltCMatrix", "itCMatrix", "dtCMatrix", "ztCMatrix"

#define VALID_RSPARSE \
"ngRMatrix", "lgRMatrix", "igRMatrix", "dgRMatrix", "zgRMatrix", \
"nsRMatrix", "lsRMatrix", "isRMatrix", "dsRMatrix", "zsRMatrix", \
                                       "dpRMatrix", "zpRMatrix", \
"ntRMatrix", "ltRMatrix", "itRMatrix", "dtRMatrix", "ztRMatrix"

#define VALID_TSPARSE \
"ngTMatrix", "lgTMatrix", "igTMatrix", "dgTMatrix", "zgTMatrix", \
"nsTMatrix", "lsTMatrix", "isTMatrix", "dsTMatrix", "zsTMatrix", \
                                       "dpTMatrix", "zpTMatrix", \
"ntTMatrix", "ltTMatrix", "itTMatrix", "dtTMatrix", "ztTMatrix"

#define VALID_DENSE \
VALID_UNPACKED, VALID_PACKED

#define VALID_SPARSE \
VALID_CSPARSE, VALID_RSPARSE, VALID_TSPARSE

#define VALID_SPARSE_COMPRESSED \
VALID_CSPARSE, VALID_RSPARSE

#define VALID_SPARSE_TRIPLET \
VALID_TSPARSE

#define VALID_DIAGONAL \
"ndiMatrix", "ldiMatrix", "idiMatrix", "ddiMatrix", "zdiMatrix"

#define VALID_INDEX \
"indMatrix",   "pMatrix"

#define VALID_MATRIX \
VALID_DENSE, VALID_SPARSE, VALID_DIAGONAL, VALID_INDEX

#define VALID_VECTOR \
"nsparseVector", "lsparseVector", "isparseVector", "dsparseVector", "zsparseVector"

#define VALID_MATRIX_OR_VECTOR \
VALID_MATRIX, VALID_VECTOR

#endif /* MATRIX_MCLASSES_H */
