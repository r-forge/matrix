## new() does not work at build time because native symbols such as
## 'Matrix_validate' needed for validity checking are not available ...
.new <- function(cl, ...) {
    def <- getClassDef(cl)
    structure(def@prototype, class = def@className, ...)
}


########################################################################
##  1. Matrix
########################################################################

## ====== Virtual Subclasses ===========================================

## ------ The Mother Class 'Matrix' ------------------------------------

## Virtual class of all Matrix objects
setClass("Matrix",
         contains = "VIRTUAL",
         slots = c(Dim = "integer", Dimnames = "list"),
         prototype = list(Dim = integer(2L), Dimnames = list(NULL, NULL)),
         validity = function(object) .Call(Matrix_validate, object))


## ------ Virtual by data type -----------------------------------------

## Virtual class of nonzero pattern matrices
setClass("nMatrix",
         contains = c("VIRTUAL", "Matrix"))

## Virtual class of logical matrices,
setClass("lMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "logical"),
         validity = function(object) .Call(lMatrix_validate, object))

## Virtual class of integer matrices
setClass("iMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "integer"),
         validity = function(object) .Call(iMatrix_validate, object))

## Virtual class of double matrices
setClass("dMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "numeric"),
         validity = function(object) .Call(dMatrix_validate, object))

## Virtual class of complex matrices
setClass("zMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "complex"),
         validity = function(object) .Call(zMatrix_validate, object))


## ------ Virtual by structure -----------------------------------------

## Virtual class of general matrices
setClass("generalMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(factors = "list"),
         validity = function(object) .Call(generalMatrix_validate, object))

## Virtual class of Hermitian or symmetric matrices
setClass("symmetricMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(uplo = "character", factors = "list"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(symmetricMatrix_validate, object))

## Virtual class of positive semidefinite matrices
setClass("posdefMatrix",
         contains = c("VIRTUAL", "symmetricMatrix"))

## Virtual class of triangular matrices
setClass("triangularMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(uplo = "character", diag = "character"),
         prototype = list(uplo = "U", diag = "N"),
         validity = function(object) .Call(triangularMatrix_validate, object))


## ------ Virtual by storage format ------------------------------------

## Virtual class of dense format matrices
setClass("denseMatrix",
         contains = c("VIRTUAL", "Matrix"))

## Virtual class of conventional, dense format matrices
setClass("unpackedMatrix",
         contains = c("VIRTUAL", "denseMatrix"),
         validity = function(object) .Call(unpackedMatrix_validate, object))

## Virtual class of packed, dense format matrices
setClass("packedMatrix",
         contains = c("VIRTUAL", "denseMatrix"),
         slots = c(uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(packedMatrix_validate, object))

## Virtual class of sparse (:= not dense) format matrices
setClass("sparseMatrix",
         contains = c("VIRTUAL", "Matrix"))

## Virtual class of compressed sparse column (CSC) format matrices
setClass("CsparseMatrix",
         contains = c("VIRTUAL", "sparseMatrix"),
         slots = c(i = "integer", p = "integer"),
         prototype = list(p = 0L), # to be valid
         validity = function(object) .Call(CsparseMatrix_validate, object))

## Virtual class of compressed sparse row (CSR) format matrices
setClass("RsparseMatrix",
         contains = c("VIRTUAL", "sparseMatrix"),
         slots = c(p = "integer", j = "integer"),
         prototype = list(p = 0L), # to be valid
         validity = function(object) .Call(RsparseMatrix_validate, object))

## Virtual class of triplet format matrices
setClass("TsparseMatrix",
         contains = c("VIRTUAL", "sparseMatrix"),
         slots = c(i = "integer", j = "integer"),
         validity = function(object) .Call(TsparseMatrix_validate, object))

## Virtual class of diagonal matrices
setClass("diagonalMatrix",
         contains = c("VIRTUAL", "sparseMatrix"),
         slots = c(diag = "character"),
         prototype = list(diag = "N"),
         validity = function(object) .Call(diagonalMatrix_validate, object))


## ------ Virtual intersections ----------------------------------------

setClass("ndenseMatrix",
         contains = c("VIRTUAL", "nMatrix", "denseMatrix"),
         slots = c(x = "logical"),
         validity = function(object) .Call(nMatrix_validate, object))

setClass("ldenseMatrix",
         contains = c("VIRTUAL", "lMatrix", "denseMatrix"))

setClass("idenseMatrix",
         contains = c("VIRTUAL", "iMatrix", "denseMatrix"))

setClass("ddenseMatrix",
         contains = c("VIRTUAL", "dMatrix", "denseMatrix"))

setClass("zdenseMatrix",
         contains = c("VIRTUAL", "zMatrix", "denseMatrix"))

setClass("nsparseMatrix",
         contains = c("VIRTUAL", "nMatrix", "sparseMatrix"))

setClass("lsparseMatrix",
         contains = c("VIRTUAL", "lMatrix", "sparseMatrix"))

setClass("isparseMatrix",
         contains = c("VIRTUAL", "iMatrix", "sparseMatrix"))

setClass("dsparseMatrix",
         contains = c("VIRTUAL", "dMatrix", "sparseMatrix"))

setClass("zsparseMatrix",
         contains = c("VIRTUAL", "zMatrix", "sparseMatrix"))


## ====== Non-Virtual Subclasses =======================================

## ------ Non-Virtual Dense --------------------------------------------

## ...... Dense, nonzero pattern .......................................

## Unpacked, general
setClass("ngeMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("nsyMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "symmetricMatrix"))

## Unpacked, triangular
setClass("ntrMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("nspMatrix",
         contains = c("packedMatrix", "ndenseMatrix", "symmetricMatrix"))

## Packed, triangular
setClass("ntpMatrix",
         contains = c("packedMatrix", "ndenseMatrix", "triangularMatrix"))


## ...... Dense, logical ...............................................

## Unpacked, general
setClass("lgeMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("lsyMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "symmetricMatrix"))

## Unpacked, triangular
setClass("ltrMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("lspMatrix",
         contains = c("packedMatrix", "ldenseMatrix", "symmetricMatrix"))

## Packed, triangular
setClass("ltpMatrix",
         contains = c("packedMatrix", "ldenseMatrix", "triangularMatrix"))


## ...... Dense, integer ...............................................

## Unpacked, general
setClass("igeMatrix",
         contains = c("unpackedMatrix", "idenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("isyMatrix",
         contains = c("unpackedMatrix", "idenseMatrix", "symmetricMatrix"))

## Unpacked, triangular
setClass("itrMatrix",
         contains = c("unpackedMatrix", "idenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("ispMatrix",
         contains = c("packedMatrix", "idenseMatrix", "symmetricMatrix"))

## Packed, triangular
setClass("itpMatrix",
         contains = c("packedMatrix", "idenseMatrix", "triangularMatrix"))


## ...... Dense, double ................................................

## Unpacked, general
setClass("dgeMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("dsyMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "symmetricMatrix"))

## Unpacked, symmetric, positive semidefinite
setClass("dpoMatrix",
         contains = c("dsyMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpoMatrix_validate, object))

## Unpacked, symmetric, positive semidefinite, correlation
setClass("corMatrix",
         contains = "dpoMatrix",
         slots = c(sd = "numeric"),
         validity = function(object) .Call(corMatrix_validate, object))

## Unpacked, triangular
setClass("dtrMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("dspMatrix",
         contains = c("packedMatrix", "ddenseMatrix", "symmetricMatrix"))

## Packed, symmetric, positive semidefinite
setClass("dppMatrix",
         contains = c("dspMatrix", "posdefMatrix"),
         validity = function(object) .Call(xppMatrix_validate, object))

## Packed, symmetric, positive semidefinite, correlation
setClass("copMatrix",
         contains = "dppMatrix",
         slots = c(sd = "numeric"),
         validity = function(object) .Call(copMatrix_validate, object))

## Packed, triangular
setClass("dtpMatrix",
         contains = c("packedMatrix", "ddenseMatrix", "triangularMatrix"))


## ...... Dense, complex ...............................................

## Unpacked, general
setClass("zgeMatrix",
         contains = c("unpackedMatrix", "zdenseMatrix", "generalMatrix"))

## Unpacked, Hermitian or symmetric
setClass("zsyMatrix",
         slots = c(trans = "character"),
         prototype = list(trans = "C"),
         contains = c("unpackedMatrix", "zdenseMatrix", "symmetricMatrix"))

## Unpacked, Hermitian, positive semidefinite
setClass("zpoMatrix",
         contains = c("zsyMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpoMatrix_validate, object))

## Unpacked, triangular
setClass("ztrMatrix",
         contains = c("unpackedMatrix", "zdenseMatrix", "triangularMatrix"))

## Packed, Hermitian or symmetric
setClass("zspMatrix",
         slots = c(trans = "character"),
         prototype = list(trans = "C"),
         contains = c("packedMatrix", "zdenseMatrix", "symmetricMatrix"))

## Packed, Hermitian, positive semidefinite
setClass("zppMatrix",
         contains = c("zspMatrix", "posdefMatrix"),
         validity = function(object) .Call(xppMatrix_validate, object))

## Packed, triangular
setClass("ztpMatrix",
         contains = c("packedMatrix", "zdenseMatrix", "triangularMatrix"))


## ------ Non-Virtual Sparse -------------------------------------------

## ...... Sparse, nonzero pattern ......................................

## CSC, general
setClass("ngCMatrix",
         contains = c("CsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSC, symmetric
setClass("nsCMatrix",
         contains = c("CsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(sCMatrix_validate, object))

## CSC, triangular
setClass("ntCMatrix",
         contains = c("CsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(tCMatrix_validate, object))

## CSR, general
setClass("ngRMatrix",
         contains = c("RsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSR, symmetric
setClass("nsRMatrix",
         contains = c("RsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(sRMatrix_validate, object))

## CSR, triangular
setClass("ntRMatrix",
         contains = c("RsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(tRMatrix_validate, object))

## Triplet general
setClass("ngTMatrix",
         contains = c("TsparseMatrix", "nsparseMatrix", "generalMatrix"))

## Triplet, symmetric
setClass("nsTMatrix",
         contains = c("TsparseMatrix", "nsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(sTMatrix_validate, object))

## Triplet, triangular
setClass("ntTMatrix",
         contains = c("TsparseMatrix", "nsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(tTMatrix_validate, object))

## Diagonal
setClass("ndiMatrix",
         contains = c("diagonalMatrix", "nMatrix"),
         slots = c(x = "logical"),
         validity = function(object) .Call(nMatrix_validate, object))


## ...... Sparse, logical ..............................................

## CSC, general
setClass("lgCMatrix",
         contains = c("CsparseMatrix", "lsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, symmetric
setClass("lsCMatrix",
         contains = c("CsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSC, triangular
setClass("ltCMatrix",
         contains = c("CsparseMatrix", "lsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtCMatrix_validate, object))

## CSR, general
setClass("lgRMatrix",
         contains = c("RsparseMatrix", "lsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, symmetric
setClass("lsRMatrix",
         contains = c("RsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## CSR, triangular
setClass("ltRMatrix",
         contains = c("RsparseMatrix", "lsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtRMatrix_validate, object))

## Triplet, general
setClass("lgTMatrix",
         contains = c("TsparseMatrix", "lsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, symmetric
setClass("lsTMatrix",
         contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, triangular
setClass("ltTMatrix",
         contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("ldiMatrix",
         contains = c("diagonalMatrix", "lMatrix"))


## ...... Sparse, integer ..............................................

## CSC, general
setClass("igCMatrix",
         contains = c("CsparseMatrix", "isparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, symmetric
setClass("isCMatrix",
         contains = c("CsparseMatrix", "isparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSC, triangular
setClass("itCMatrix",
         contains = c("CsparseMatrix", "isparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtCMatrix_validate, object))

## CSR, general
setClass("igRMatrix",
         contains = c("RsparseMatrix", "isparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, symmetric
setClass("isRMatrix",
         contains = c("RsparseMatrix", "isparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## CSR, triangular
setClass("itRMatrix",
         contains = c("RsparseMatrix", "isparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtRMatrix_validate, object))

## Triplet, general
setClass("igTMatrix",
         contains = c("TsparseMatrix", "isparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, symmetric
setClass("isTMatrix",
         contains = c("TsparseMatrix", "isparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, triangular
setClass("itTMatrix",
         contains = c("TsparseMatrix", "isparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("idiMatrix",
         contains = c("diagonalMatrix", "iMatrix"))


## ...... Sparse, double ...............................................

## CSC, general
setClass("dgCMatrix",
         contains = c("CsparseMatrix", "dsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, symmetric
setClass("dsCMatrix",
         contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSC, symmetric, positive semidefinite
setClass("dpCMatrix",
         contains = c("dsCMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpCMatrix_validate, object))

## CSC, triangular
setClass("dtCMatrix",
         contains = c("CsparseMatrix", "dsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtCMatrix_validate, object))

## CSR, general
setClass("dgRMatrix",
         contains = c("RsparseMatrix", "dsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, symmetric
setClass("dsRMatrix",
         contains = c("RsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## CSR, symmetric, positive semidefinite
setClass("dpRMatrix",
         contains = c("dsRMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpRMatrix_validate, object))

## CSR, triangular
setClass("dtRMatrix",
         contains = c("RsparseMatrix", "dsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtRMatrix_validate, object))

## Triplet, general
setClass("dgTMatrix",
         contains = c("TsparseMatrix", "dsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, symmetric
setClass("dsTMatrix",
         contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, symmetric, positive semidefinite
setClass("dpTMatrix",
         contains = c("dsTMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpTMatrix_validate, object))

## Triplet, triangular
setClass("dtTMatrix",
         contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("ddiMatrix",
         contains = c("diagonalMatrix", "dMatrix"))


## ...... Sparse, complex ..............................................

## CSC, general
setClass("zgCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, Hermitian or symmetric
setClass("zsCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "symmetricMatrix"),
         slots = c(trans = "character"),
         prototype = list(trans = "C"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSC, Hermitian, positive semidefinite
setClass("zpCMatrix",
         contains = c("zsCMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpCMatrix_validate, object))

## CSC, triangular
setClass("ztCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtCMatrix_validate, object))

## CSR, general
setClass("zgRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, Hermitian or symmetric
setClass("zsRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "symmetricMatrix"),
         slots = c(trans = "character"),
         prototype = list(trans = "C"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## CSR, Hermitian, positive semidefinite
setClass("zpRMatrix",
         contains = c("zsRMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpRMatrix_validate, object))

## CSR, triangular
setClass("ztRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtRMatrix_validate, object))

## Triplet, general
setClass("zgTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, Hermitian or symmetric
setClass("zsTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "symmetricMatrix"),
         slots = c(trans = "character"),
         prototype = list(trans = "C"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, Hermitian, positive semidefinite
setClass("zpTMatrix",
         contains = c("zsTMatrix", "posdefMatrix"),
         validity = function(object) .Call(xpTMatrix_validate, object))

## Triplet, triangular
setClass("ztTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("zdiMatrix",
         contains = c("diagonalMatrix", "zMatrix"))


## ...... Sparse, index ................................................

## Row or column index
setClass("indMatrix",
         contains = "sparseMatrix",
         slots = c(perm = "integer", margin = "integer"),
         prototype = list(margin = 1L), # to be valid
         validity = function(object) .Call(indMatrix_validate, object))

## Row or column permutation
setClass("pMatrix",
         contains = "indMatrix",
         validity = function(object) .Call(pMatrix_validate, object))


########################################################################
##  2. MatrixFactorization
########################################################################

## ------ The Mother Class "MatrixFactorization" -----------------------

setClass("MatrixFactorization",
         contains = "VIRTUAL",
         slots = c(Dim = "integer", Dimnames = "list"),
         prototype = list(Dim = integer(2L), Dimnames = list(NULL, NULL)),
         validity = function(object) .Call(MatrixFactorization_validate, object))


## ------ Schur --------------------------------------------------------

setClass("Schur",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseSchur",
         contains = c("VIRTUAL", "Schur"),
         slots = c(x = "vector", values = "vector", vectors = "vector"),
         prototype = list(values = double(0L)),
         validity = function(object) .Call(denseSchur_validate, object))

setClass("ddenseSchur",
         contains = "denseSchur",
         slots = c(x = "numeric", vectors = "numeric"),
         prototype = list(x = double(0L), vectors = double(0L)))

setClass("zdenseSchur",
         contains = "denseSchur",
         slots = c(x = "complex", vectors = "complex"),
         prototype = list(x = complex(0L), vectors = complex(0L)))


## ------ QR -----------------------------------------------------------

setClass("QR",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseQR",
         contains = c("VIRTUAL", "QR"),
         slots = c(x = "vector", beta = "vector", perm = "integer"),
         validity = function(object) .Call(denseQR_validate, object))

setClass("ddenseQR",
         contains = "denseQR",
         slots = c(x = "numeric", beta = "numeric"),
         prototype = list(x = double(0L), beta = double(0L)))

setClass("zdenseQR",
         contains = "denseQR",
         slots = c(x = "complex", beta = "complex"),
         prototype = list(x = complex(0L), beta = complex(0L)))

setClass("sparseQR",
         contains = c("VIRTUAL", "QR"),
         slots = c(V = "Matrix", R = "Matrix", beta = "vector",
                   p = "integer", q = "integer"),
         validity = function(object) .Call(sparseQR_validate, object))

setClass("dsparseQR",
         contains = "sparseQR",
         slots = c(V = "dgCMatrix", R = "dgCMatrix", beta = "numeric"),
         prototype = list(V = .new("dgCMatrix"),
                          R = .new("dgCMatrix"),
                          beta = double(0L)))

setClass("zsparseQR",
         contains = "sparseQR",
         slots = c(V = "zgCMatrix", R = "zgCMatrix", beta = "numeric"),
         prototype = list(V = .new("zgCMatrix"),
                          R = .new("zgCMatrix"),
                          beta = double(0L)))


## ------ LU -----------------------------------------------------------

setClass("LU",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseLU",
         contains = c("VIRTUAL", "LU"),
         slots = c(x = "vector", perm = "integer"),
         validity = function(object) .Call(denseLU_validate, object))

setClass("ddenseLU",
         contains = "denseLU",
         slots = c(x = "numeric"),
         prototype = list(x = double(0L)))

setClass("zdenseLU",
         contains = "denseLU",
         slots = c(x = "complex"),
         prototype = list(x = complex(0L)))

setClass("sparseLU",
         contains = c("VIRTUAL", "LU"),
         slots = c(L = "Matrix", U = "Matrix",
                   p = "integer", q = "integer"),
         validity = function(object) .Call(sparseLU_validate, object))

setClass("dsparseLU",
         contains = "sparseLU",
         slots = c(L = "dtCMatrix", U = "dtCMatrix"),
         prototype = list(L = .new("dtCMatrix", uplo = "L"),
                          U = .new("dtCMatrix")))

setClass("zsparseLU",
         contains = "sparseLU",
         slots = c(L = "ztCMatrix", U = "ztCMatrix"),
         prototype = list(L = .new("ztCMatrix", uplo = "L"),
                          U = .new("ztCMatrix")))


## ------ Bunch-Kaufman ------------------------------------------------

setClass("BunchKaufman",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseBunchKaufman",
         contains = c("VIRTUAL", "BunchKaufman"),
         slots = c(x = "vector", perm = "integer", uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(denseBunchKaufman_validate, object))

setClass("ddenseBunchKaufman",
         contains = "denseBunchKaufman",
         slots = c(x = "numeric"),
         prototype = list(x = double(0L)))

setClass("zdenseBunchKaufman",
         contains = "denseBunchKaufman",
         slots = c(x = "complex", trans = "character"),
         prototype = list(x = complex(0L), trans = "C"))


## ------ Cholesky -----------------------------------------------------

setClass("Cholesky",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseCholesky",
         contains = c("VIRTUAL", "Cholesky"),
         slots = c(x = "vector", perm = "integer", uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(denseCholesky_validate, object))

setClass("ddenseCholesky",
         contains = "denseCholesky",
         slots = c(x = "numeric"),
         prototype = list(x = double(0L)))

setClass("zdenseCholesky",
         contains = "denseCholesky",
         slots = c(x = "complex"),
         prototype = list(x = complex(0L)))

setClass("sparseCholesky",
         contains = c("VIRTUAL", "Cholesky"),
         slots = c(perm = "integer", colcount = "integer",
                   ordering = "integer"),
         prototype = list(ordering = 0L),
         validity = function(object) .Call(sparseCholesky_validate, object))

setClass("nsparseCholesky",
         contains = c("VIRTUAL", "sparseCholesky"))

setClass("dsparseCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(minor = "integer", x = "numeric"),
         prototype = list(minor = 0L, x = double(0L)))

setClass("zsparseCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(minor = "integer", x = "complex"),
         prototype = list(minor = 0L, x = complex(0L)))

setClass("simplicialCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         validity = function(object) .Call(simplicialCholesky_validate, object))

setClass("nsimplicialCholesky",
         contains = c("nsparseCholesky", "simplicialCholesky"))

setClass("dsimplicialCholesky",
         contains = c("dsparseCholesky", "simplicialCholesky"),
         slots = c(p = "integer", i = "integer", nz = "integer",
                   `next` = "integer", prev = "integer",
                   is_ll = "logical", is_monotonic = "logical"),
         prototype = list(p = 0L, `next` = c(-1L, 0L), prev = c(1L, -1L),
                          is_ll = FALSE, is_monotonic = TRUE))

setClass("zsimplicialCholesky",
         contains = c("zsparseCholesky", "simplicialCholesky"),
         slots = c(p = "integer", i = "integer", nz = "integer",
                   `next` = "integer", prev = "integer",
                   is_ll = "logical", is_monotonic = "logical"),
         prototype = list(p = 0L, `next` = c(-1L, 0L), prev = c(1L, -1L),
                          is_ll = FALSE, is_monotonic = TRUE))

setClass("supernodalCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(maxcsize = "integer", maxesize = "integer",
                   super = "integer", pi = "integer", px = "integer",
                   s = "integer"),
         prototype = list(maxcsize = 0L, maxesize = 0L,
                          super = 0L, pi = 0L, px = 0L),
         validity = function(object) .Call(supernodalCholesky_validate, object))

setClass("nsupernodalCholesky",
         contains = c("nsparseCholesky", "supernodalCholesky"))

setClass("dsupernodalCholesky",
         contains = c("dsparseCholesky", "supernodalCholesky"))

setClass("zsupernodalCholesky",
         contains = c("zsparseCholesky", "supernodalCholesky"))


########################################################################
##  3. sparseVector
########################################################################

## ------ The Mother Class 'sparseVector' ------------------------------

setClass("sparseVector",
         contains = "VIRTUAL",
         slots = c(length = "numeric", i = "numeric"),
         prototype = list(length = 0L, i = integer(0L)),
         validity = function(object) .Call(sparseVector_validate, object))


## ------ Non-Virtual Subclasses ---------------------------------------

setClass("nsparseVector",
         contains = "sparseVector")

setClass("lsparseVector",
         contains = "sparseVector",
         slots = c(x = "logical"),
         validity = function(object) .Call(lsparseVector_validate, object))

setClass("isparseVector",
         contains = "sparseVector",
         slots = c(x = "integer"),
         validity = function(object) .Call(isparseVector_validate, object))

setClass("dsparseVector",
         contains = "sparseVector",
         slots = c(x = "numeric"),
         validity = function(object) .Call(dsparseVector_validate, object))

setClass("zsparseVector",
         contains = "sparseVector",
         slots = c(x = "complex"),
         validity = function(object) .Call(zsparseVector_validate, object))


########################################################################
##  4. Index and more "miscellaneous" classes
########################################################################

## Idea: represent x = c(seq(from1, to1, by1), seq(from2, to2, by2), ...)
##       as list(first = x[1L], rle = rle(diff(x)))
setClass("rleDiff",
         ## MJ: simpler would be slots = c(first=, lengths=, values=) ...
         slots = c(first = "numeric", rle = "rle"),
         prototype = list(first = integer(0L), rle = rle(integer(0L))),
         validity = function(object) {
             if(length(object@first) != 1L)
                 "'first' slot does not have length 1"
             else if(!is.list(rle <- object@rle))
                 "'rle' slot is not a list"
             else if(length(rle) != 2L)
                 "'rle' slot does not have length 2"
             else if(is.null(nms <- names(rle)) ||
                     anyNA(match(nms, c("lengths", "values"))))
                 "'rle' slot does not have names \"lengths\", \"values\""
             else if(!is.numeric(lens <- rle$lengths))
                 "'lengths' is not numeric"
             else if(!is.numeric(vals <- rle$values))
                 "'values' is not numeric"
             else if(length(lens) != length(vals))
                 "'lengths' and 'values' do not have equal length"
             else if(length(lens) == 0L)
                 TRUE
             else if(anyNA(lens))
                 "'lengths' contains NA"
             else if(is.double(lens)) {
                 if(!(all(is.finite(r <- range(lens))) &&
                      all(lens == trunc(lens))))
                     "'lengths' is not integer-valued"
                 else if(r[1L] < 1)
                     "'lengths' is not positive"
                 else TRUE
             } else {
                 if(min(lens) < 1L)
                     "'lengths' is not positive"
                 else TRUE
             }
         })

## Idea: represent x = c(seq(from1, to1, by1), seq(from2, to2, by2), ...)
##       as rbind(c(from1, from2, ...), c(to1, to2, ...), c(by1, by2, ...))
## MM: (2010-03-04) more efficient than "rleDiff" [TODO: write rleDiff<->seqMat]
## MJ: (2022-09-06) data.frame(from, to, by) could be _handled_ more efficiently
setClass("seqMat",
         contains = "matrix",
         prototype = matrix(integer(0L), nrow = 3L, ncol = 0L),
         validity = function(object) {
             if(!is.numeric(object))
                 "matrix is not numeric"
             else if(nrow(object) != 3L)
                 "matrix does not have 3 rows"
             else if(anyNA(object))
                 "matrix contains NA"
             else if(is.double(object) && !(all(is.finite(range(object))) &&
                                            all(object == trunc(object))))
                 "matrix is not integer-valued"
             else {
                 from <- object[1L, ]
                 to   <- object[2L, ]
                 by   <- object[3L, ]
                 if(any((from < to & by <= 0) | (from > to & by >= 0)))
                     "degenerate sequence(s): sign(to - from) != sign(by)"
                 else TRUE
             }
         })

## Idea: _ab_stract index
## MJ: (2022-09-06) why not just
##     setClassUnion("abIndex", members = c("numeric", "rleDiff", "seqMat")) ?
setClass("abIndex",
         slots = c(kind = "character", x = "numeric", rleD = "rleDiff"),
         prototype = list(kind = "int32", x = integer(0L)),
         validity = function(object) {
             ## MJ: should 'rleD' be "empty" if kind != "rleDiff" ?
             if(length(kind <- object@kind) != 1L)
                 "'kind' slot does not have length 1"
             else switch(kind,
                         "int32" =
                             if(is.integer(object@x))
                                 TRUE
                             else "kind=\"int32\" but 'x' slot is not of type \"integer\"",
                         "double" =
                             if(is.double(object@x))
                                 TRUE
                             else "kind=\"double\" but 'x' slot is not of type \"double\"",
                         "rleDiff" =
                             if(length(object@x) == 0L)
                                 TRUE
                             else "kind=\"rleDiff\" but 'x' slot is nonempty",
                         ## otherwise:
                         "'kind' is not \"int32\", \"double\", or \"rleDiff\"")
         })


########################################################################
##  5. Class unions
########################################################################

## MJ: aim to deprecate and eventually remove these
##     (except perhaps 'index')
## MM: agreed that "raw" is non-sense; couldn't we have method-signatures (& dispatch),
##     e.g., for "Ops" / "Arith" ... using
## --  genNum3 := {logical,numeric,complex}  and/or
## --  genNum4 := {logical,numeric,complex , integer} ?

setClassUnion("atomicVector",
              members = c("raw", "logical", "numeric", "complex", "character"))
setClassUnion("index",
              members = c(       "logical", "numeric",            "character"))
setClassUnion("numLike",
              members = c(       "logical", "numeric"                        ))
setClassUnion("number",
              members = c(                  "numeric", "complex"             ))
setClassUnion("replValue",
              members = c("raw", "logical", "numeric", "complex"             ))

## Removing these entirely in Matrix 1.7-0 invalidates class definitions
## serialized in existing installations of the following packages:
##
##  [1] CVXR                 FinNet               GENLIB
##  [4] MSnbase              MachineShop          MatrixModels
##  [7] SeuratObject         SingleCellExperiment WoodburyMatrix
## [10] apcluster            arules               chromVAR
## [13] conText              copula               destiny
## [16] distrom              genomation           hypr
## [19] iGraphMatch          kebabs               mcompanion
## [22] pcts                 podkat               qpgraph
## [25] quadrupen            quanteda             quanteda.textstats
## [28] saeRobust            scuttle              snpStats
## [31] softImpute           spflow               xcms
##
## Define stubs so that the serialized class definitions do not cause
## S4 machinery to throw warnings or errors.  Remove the stubs once
## binaries in most repositories seem to have been rebuilt under 1.7-0.
setClass("compMatrix")
setClass("pcorMatrix")
setClassUnion("replValueSp")

## Matrix 1.8-0
setClass("BunchKaufmanFactorization")
setClass("CHMfactor")
setClass("CHMsimpl")
setClass("CHMsuper")
setClass("CholeskyFactorization")
setClass("SchurFactorization")
setClass("dCHMsimpl")
setClass("dCHMsuper")
setClass("nCHMsimpl")
setClass("nCHMsuper")
setClass("pBunchKaufman")
setClass("pCholesky")

rm(.new)
