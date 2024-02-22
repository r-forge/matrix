## README:
##
## Validity methods should assume that methods for superclasses have passed,
## following validObject(). We should _not_ be testing, e.g., length(Dim),
## typeof(Dimnames), etc. repeatedly ...
##
## When checking whether a class is validated correctly, take care to follow
## the 'contains' recursively!!


## initialize() method for Matrix and MatrixFactorization, which both
## allow Dimnames[[i]] to be a vector of type other than "character"
## and furthermore to be a vector of length zero rather than NULL ...
.initialize <- function(.Object, ...) {
    .Object <- callNextMethod()
    ## Suboptimal if ...names() is NULL but that will "never"
    ## happen if ...length() is nonzero:
    if(...length() && any(...names() == "Dimnames"))
        .Object@Dimnames <- fixupDN(.Object@Dimnames)
    .Object
}

## new() does not work at build time because native symbols such as
## 'Matrix_validate' needed for validity checking are not available ...
.new <- function(cl, ...) {
    def <- getClassDef(cl)
    structure(def@prototype, class = def@className, ...)
}

.M2 <- nzchar(Sys.getenv("R_MATRIX_MATRIX2"))


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

setMethod("initialize", c(.Object = "Matrix"),
          .initialize)


## ------ Virtual by structure -----------------------------------------

## Virtual class of composite matrices, having "nontrivial" factorizations
setClass("compMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(factors = "list"),
         validity = function(object) .Call(compMatrix_validate, object))

## Virtual class of general matrices
setClass("generalMatrix",
         contains = c("VIRTUAL", "compMatrix"))

## Virtual class of symmetric matrices
setClass("symmetricMatrix",
         contains = c("VIRTUAL", "compMatrix"),
         slots = c(uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(symmetricMatrix_validate, object))

if(.M2) {
## Virtual class of Hermitian matrices
setClass("HermitianMatrix",
         contains = c("VIRTUAL", "compMatrix"),
         slots = c(uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(HermitianMatrix_validate, object))
}

## Virtual class of triangular matrices
setClass("triangularMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(uplo = "character", diag = "character"),
         prototype = list(uplo = "U", diag = "N"),
         validity = function(object) .Call(triangularMatrix_validate, object))


## ------ Virtual by kind ----------------------------------------------

## Virtual class of nonzero pattern matrices
setClass("nMatrix",
         contains = c("VIRTUAL", "Matrix"))

## Virtual class of logical matrices,
setClass("lMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "logical"),
         validity = function(object) .Call(lMatrix_validate, object))

if(.M2 || TRUE) {
## Virtual class of integer matrices
setClass("iMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "integer"),
         validity = function(object) .Call(iMatrix_validate, object))
}

## Virtual class of double matrices
setClass("dMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "numeric"),
         validity = function(object) .Call(dMatrix_validate, object))

if(.M2 || TRUE) {
## Virtual class of complex matrices
setClass("zMatrix",
         contains = c("VIRTUAL", "Matrix"),
         slots = c(x = "complex"),
         validity = function(object) .Call(zMatrix_validate, object))
}


## ------ Virtual Dense ------------------------------------------------

## Virtual class of dense matrices
setClass("denseMatrix",
         contains = c("VIRTUAL", "Matrix"))


## ...... Virtual Dense ... by storage .................................

## Virtual class of dense, unpacked matrices
setClass("unpackedMatrix",
         contains = c("VIRTUAL", "denseMatrix"),
         validity = function(object) .Call(unpackedMatrix_validate, object))

## Virtual class of dense, packed matrices
setClass("packedMatrix",
         contains = c("VIRTUAL", "denseMatrix"),
         slots = c(uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(packedMatrix_validate, object))


## ...... Virtual Dense ... by kind ....................................

## Virtual class of dense, nonzero pattern matrices
setClass("ndenseMatrix",
         contains = c("VIRTUAL", "nMatrix", "denseMatrix"),
         slots = c(x = "logical"),
         validity = function(object) .Call(nMatrix_validate, object))

## Virtual class of dense, logical matrices
setClass("ldenseMatrix",
         contains = c("VIRTUAL", "lMatrix", "denseMatrix"))

if(.M2) {
## Virtual class of dense, integer matrices
setClass("idenseMatrix",
         contains = c("VIRTUAL", "iMatrix", "denseMatrix"))
}

## Virtual class of dense, double matrices
setClass("ddenseMatrix",
         contains = c("VIRTUAL", "dMatrix", "denseMatrix"))

if(.M2) {
## Virtual class of dense, complex matrices
setClass("zdenseMatrix",
         contains = c("VIRTUAL", "zMatrix", "denseMatrix"))
}


## ------ Virtual Sparse -----------------------------------------------

## Virtual class of sparse matrices
setClass("sparseMatrix",
         contains = c("VIRTUAL", "Matrix"))


## ...... Virtual Sparse ... by storage ................................

## Virtual class of sparse matrices in compressed sparse column (CSC) format
setClass("CsparseMatrix",
         contains = c("VIRTUAL", "sparseMatrix"),
         slots = c(i = "integer", p = "integer"),
         prototype = list(p = 0L), # to be valid
         validity = function(object) .Call(CsparseMatrix_validate, object))

## Virtual class of sparse matrices in compressed sparse row (CSR) format
setClass("RsparseMatrix",
         contains = c("VIRTUAL", "sparseMatrix"),
         slots = c(p = "integer", j = "integer"),
         prototype = list(p = 0L), # to be valid
         validity = function(object) .Call(RsparseMatrix_validate, object))

## Virtual class of sparse matrices in triplet format
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

if(FALSE) { # --NOT YET--
## These methods would allow initialization of zero matrices _without_ 'p',
## as in the call new("dgCMatrix", Dim = c(6L, 6L)).  However, they would
## also incur a small performance penalty on all other new("..[CR]Matrix")
## calls.
setMethod("initialize", c(.Object = "CsparseMatrix"),
          function(.Object, ...) {
              ## Suboptimal if ...names() is NULL or if 'Dim' is missing
              ## but that will "never" happen if ...length() is nonzero:
              if(...length() &&
                 all((nms <- ...names()) != "p") &&
                 length(w <- which(nms == "Dim")) &&
                 !is.character(validDim(d <- ...elt(w[1L]))))
                  callNextMethod(.Object, ..., p = integer(d[2L] + 1))
              else callNextMethod()
          })

setMethod("initialize", c(.Object = "RsparseMatrix"),
          function(.Object, ...) {
              ## Suboptimal if ...names() is NULL or if 'Dim' is missing
              ## but that will "never" happen if ...length() is nonzero:
              if(...length() &&
                 all((nms <- ...names()) != "p") &&
                 length(w <- which(nms == "Dim")) &&
                 !is.character(validDim(d <- ...elt(w[1L]))))
                  callNextMethod(.Object, ..., p = integer(d[1L] + 1))
              else callNextMethod()
          })
} # --NOT YET--


## ...... Virtual Sparse ... by kind ...................................

## Virtual class of sparse, nonzero pattern matrices
setClass("nsparseMatrix",
         contains = c("VIRTUAL", "nMatrix", "sparseMatrix"))

## Virtual class of sparse, logical matrices
setClass("lsparseMatrix",
         contains = c("VIRTUAL", "lMatrix", "sparseMatrix"))

if(.M2) {
## Virtual class of sparse, integer matrices
setClass("isparseMatrix",
         contains = c("VIRTUAL", "iMatrix", "sparseMatrix"))
}

## Virtual class of sparse, double matrices
setClass("dsparseMatrix",
         contains = c("VIRTUAL", "dMatrix", "sparseMatrix"))

if(.M2) {
## Virtual class of sparse, complex matrices
setClass("zsparseMatrix",
         contains = c("VIRTUAL", "zMatrix", "sparseMatrix"))
}


## ====== Non-Virtual Subclasses =======================================

## ------ Non-Virtual Dense --------------------------------------------

## ...... Dense, nonzero pattern .......................................

## Unpacked, general
setClass("ngeMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("nsyMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Unpacked, triangular
setClass("ntrMatrix",
         contains = c("unpackedMatrix", "ndenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("nspMatrix",
         contains = c("packedMatrix", "ndenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Packed, triangular
setClass("ntpMatrix",
         contains = c("packedMatrix", "ndenseMatrix", "triangularMatrix"))


## ...... Dense, logical ...............................................

## Unpacked, general
setClass("lgeMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("lsyMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Unpacked, triangular
setClass("ltrMatrix",
         contains = c("unpackedMatrix", "ldenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("lspMatrix",
         contains = c("packedMatrix", "ldenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Packed, triangular
setClass("ltpMatrix",
         contains = c("packedMatrix", "ldenseMatrix", "triangularMatrix"))


## ...... Dense, integer ...............................................

if(.M2) {
## Unpacked, general
setClass("igeMatrix",
         contains = c("unpackedMatrix", "idenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("isyMatrix",
         contains = c("unpackedMatrix", "idenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Unpacked, triangular
setClass("itrMatrix",
         contains = c("unpackedMatrix", "idenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("ispMatrix",
         contains = c("packedMatrix", "idenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Packed, triangular
setClass("itpMatrix",
         contains = c("packedMatrix", "idenseMatrix", "triangularMatrix"))
}


## ...... Dense, double ................................................

## Unpacked, general
setClass("dgeMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("dsyMatrix",
         contains = c("unpackedMatrix", "ddenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Unpacked, symmetric, positive semidefinite
setClass("dpoMatrix",
         contains = "dsyMatrix",
         validity = function(object) .Call(dpoMatrix_validate, object))

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
         contains = c("packedMatrix", "ddenseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"))

## Packed, symmetric, positive semidefinite
setClass("dppMatrix",
         contains = "dspMatrix",
         validity = function(object) .Call(dppMatrix_validate, object))

## Packed, symmetric, positive semidefinite, correlation
setClass("pcorMatrix",
         contains = "dppMatrix",
         slots = c(sd = "numeric"),
         validity = function(object) .Call(pcorMatrix_validate, object))

## Packed, triangular
setClass("dtpMatrix",
         contains = c("packedMatrix", "ddenseMatrix", "triangularMatrix"))


## ...... Dense, complex ...............................................

if(.M2) {
## Unpacked, general
setClass("zgeMatrix",
         contains = c("unpackedMatrix", "zdenseMatrix", "generalMatrix"))

## Unpacked, symmetric
setClass("zsyMatrix",
         contains = c("unpackedMatrix", "zdenseMatrix", "symmetricMatrix"))

## Unpacked, Hermitian
setClass("zheMatrix",
         contains = c("unpackedMatrix", "zdenseMatrix", "HermitianMatrix"))

## Unpacked, Hermitian, positive semidefinite
setClass("zpoMatrix",
         contains = "zheMatrix",
         validity = function(object) .Call(zpoMatrix_validate, object))

## Unpacked, triangular
setClass("ztrMatrix",
         contains = c("unpackedMatrix", "zdenseMatrix", "triangularMatrix"))

## Packed, symmetric
setClass("zspMatrix",
         contains = c("packedMatrix", "zdenseMatrix", "symmetricMatrix"))

## Packed, Hermitian
setClass("zhpMatrix",
         contains = c("packedMatrix", "zdenseMatrix", "HermitianMatrix"))

## Packed, Hermitian, positive semidefinite
setClass("zppMatrix",
         contains = "zhpMatrix",
         validity = function(object) .Call(zppMatrix_validate, object))

## Packed, triangular
setClass("ztpMatrix",
         contains = c("packedMatrix", "zdenseMatrix", "triangularMatrix"))
}


## ------ Non-Virtual Sparse -------------------------------------------

## ...... Sparse, nonzero pattern ......................................

## CSC, general
setClass("ngCMatrix",
         contains = c("CsparseMatrix", "nsparseMatrix", "generalMatrix"))

## CSC, symmetric
setClass("nsCMatrix",
         contains = c("CsparseMatrix", "nsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("RsparseMatrix", "nsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("TsparseMatrix", "nsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("CsparseMatrix", "lsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("RsparseMatrix", "lsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, triangular
setClass("ltTMatrix",
         contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("ldiMatrix",
         contains = c("diagonalMatrix", "lMatrix"))


## ...... Sparse, integer ..............................................

if(.M2) {
## CSC, general
setClass("igCMatrix",
         contains = c("CsparseMatrix", "isparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, symmetric
setClass("isCMatrix",
         contains = c("CsparseMatrix", "isparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("RsparseMatrix", "isparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
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
         contains = c("TsparseMatrix", "isparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, triangular
setClass("itTMatrix",
         contains = c("TsparseMatrix", "isparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("idiMatrix",
         contains = c("diagonalMatrix", "iMatrix"))
}


## ...... Sparse, double ...............................................

## CSC, general
setClass("dgCMatrix",
         contains = c("CsparseMatrix", "dsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, symmetric
setClass("dsCMatrix",
         contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

if(.M2) {
## CSC, symmetric, positive semidefinite
setClass("dpCMatrix",
         contains = "dsCMatrix",
         validity = function(object) .Call(xpCMatrix_validate, object))
}

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
         contains = c("RsparseMatrix", "dsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

if(.M2) {
## CSR, symmetric, positive semidefinite
setClass("dpRMatrix",
         contains = "dsRMatrix",
         validity = function(object) .Call(xpRMatrix_validate, object))
}

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
         contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix",
                      if(.M2) "HermitianMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

if(.M2) {
## Triplet, symmetric, positive semidefinite
setClass("dpTMatrix",
         contains = "dsTMatrix",
         validity = function(object) .Call(xpTMatrix_validate, object))
}

## Triplet, triangular
setClass("dtTMatrix",
         contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("ddiMatrix",
         contains = c("diagonalMatrix", "dMatrix"))


## ...... Sparse, complex ..............................................

if(.M2) {
## CSC, general
setClass("zgCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgCMatrix_validate, object))

## CSC, symmetric
setClass("zsCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSC, Hermitian
setClass("zhCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "HermitianMatrix"),
         validity = function(object) .Call(xsCMatrix_validate, object))

## CSC, Hermitian, positive semidefinite
setClass("zpCMatrix",
         contains = "zhCMatrix",
         validity = function(object) .Call(xpCMatrix_validate, object))

## CSC, triangular
setClass("ztCMatrix",
         contains = c("CsparseMatrix", "zsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtCMatrix_validate, object))

## CSR, general
setClass("zgRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgRMatrix_validate, object))

## CSR, symmetric
setClass("zsRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## CSR, Hermitian
setClass("zhRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "HermitianMatrix"),
         validity = function(object) .Call(xsRMatrix_validate, object))

## CSR, Hermitian, positive semidefinite
setClass("zpRMatrix",
         contains = "zhRMatrix",
         validity = function(object) .Call(xpRMatrix_validate, object))

## CSR, triangular
setClass("ztRMatrix",
         contains = c("RsparseMatrix", "zsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtRMatrix_validate, object))

## Triplet, general
setClass("zgTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "generalMatrix"),
         validity = function(object) .Call(xgTMatrix_validate, object))

## Triplet, symmetric
setClass("zsTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "symmetricMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, Hermitian
setClass("zhTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "HermitianMatrix"),
         validity = function(object) .Call(xsTMatrix_validate, object))

## Triplet, Hermitian, positive semidefinite
setClass("zpTMatrix",
         contains = "zhTMatrix",
         validity = function(object) .Call(xpTMatrix_validate, object))

## Triplet, triangular
setClass("ztTMatrix",
         contains = c("TsparseMatrix", "zsparseMatrix", "triangularMatrix"),
         validity = function(object) .Call(xtTMatrix_validate, object))

## Diagonal
setClass("zdiMatrix",
         contains = c("diagonalMatrix", "zMatrix"))
}


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

setMethod("initialize", c(.Object = "MatrixFactorization"),
          .initialize)


## ------ LU -----------------------------------------------------------

setClass("LU",
         contains = c("VIRTUAL", "MatrixFactorization"))

if(.M2) {
setClass("denseLU",
         contains = c("VIRTUAL", "LU"),
         slots = c(x = "vector", perm = "integer"),
         validity = function(object) .Call(denseLU_validate, object))

setClass("ddenseLU",
         contains = "denseLU",
         slots = c(x = "numeric"))

setClass("zdenseLU",
         contains = "denseLU",
         slots = c(x = "complex"))

setClass("sparseLU",
         contains = c("VIRTUAL", "LU"),
         slots = c(L = "Matrix", U = "Matrix",
                   p = "integer", q = "integer"),
         validity = function(object) .Call(sparseLU_validate, object))

setClass("dsparseLU",
         contains = "sparseLU",
         slots = c(L = "dtCMatrix", U = "dtCMatrix"),
         prototype = list(L = .new("dtCMatrix", uplo = "L")))

setClass("zsparseLU",
         contains = "sparseLU",
         slots = c(L = "ztCMatrix", U = "ztCMatrix"),
         prototype = list(L = .new("ztCMatrix", uplo = "L")))
} else {
setClass("denseLU",
         contains = "LU",
         slots = c(x = "numeric", perm = "integer"),
         validity = function(object) {
             object. <- new("dgeMatrix")
             object.@Dim <- object@Dim
             object.@Dimnames <- object@Dimnames
             object.@x <- object@x
             if(is.character(valid <- validObject(object., test = TRUE)))
                 valid
             else .Call(denseLU_validate, object)
         })

setClass("sparseLU",
         contains = "LU",
         slots = c(L = "dtCMatrix", U = "dtCMatrix",
                   p = "integer", q = "integer"),
         prototype = list(L = .new("dtCMatrix", uplo = "L")),
         validity = function(object) .Call(sparseLU_validate, object))
}


## ------ QR -----------------------------------------------------------

setClass("QR",
         contains = c("VIRTUAL", "MatrixFactorization"))

if(.M2) {
setClass("denseQR",
         contains = c("VIRTUAL", "QR"),
         slots = c(x = "vector", beta = "vector", perm = "integer"),
         validity = function(object) .Call(denseQR_validate, object))

setClass("ddenseQR",
         contains = "denseQR",
         slots = c(x = "numeric", beta = "numeric"))

setClass("zdenseQR",
         contains = "denseQR",
         slots = c(x = "complex", beta = "complex"))

setClass("sparseQR",
         contains = c("VIRTUAL", "QR"),
         slots = c(V = "Matrix", R = "Matrix", beta = "vector",
                   p = "integer", q = "integer"),
         validity = function(object) .Call(sparseQR_validate, object))

setClass("dsparseQR",
         contains = "sparseQR",
         slots = c(V = "dgCMatrix", R = "dgCMatrix", beta = "numeric"))

setClass("zsparseQR",
         contains = "sparseQR",
         slots = c(V = "zgCMatrix", R = "zgCMatrix", beta = "complex"))
} else {
setClass("sparseQR",
         contains = "QR",
         slots = c(V = "dgCMatrix", R = "dgCMatrix", beta = "numeric",
                   p = "integer", q = "integer"),
         validity = function(object) .Call(sparseQR_validate, object))
}


## ------ Bunch-Kaufman ------------------------------------------------

if(.M2) {
setClass("BunchKaufman",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseBunchKaufman",
         contains = c("VIRTUAL", "BunchKaufman"),
         slots = c(x = "vector", perm = "integer", uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(denseBunchKaufman_validate, object))

setClass("ddenseBunchKaufman",
         contains = "denseBunchKaufman",
         slots = c(x = "numeric"))

setClass("zdenseBunchKaufman",
         contains = "denseBunchKaufman",
         slots = c(x = "complex"))
} else {
setClass("BunchKaufmanFactorization",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("BunchKaufman",
         contains = "BunchKaufmanFactorization",
         slots = c(uplo = "character", x = "numeric", perm = "integer"),
         prototype = list(uplo = "U"),
         validity = function(object) {
             object. <- new("dtrMatrix")
             object.@Dim <- object@Dim
             object.@Dimnames <- object@Dimnames
             object.@uplo <- object@uplo
             object.@x <- object@x
             if(is.character(valid <- validObject(object., test = TRUE)))
                 valid
             else .Call(BunchKaufman_validate, object)
         })

setClass("pBunchKaufman",
         contains = "BunchKaufmanFactorization",
         slots = c(uplo = "character", x = "numeric", perm = "integer"),
         prototype = list(uplo = "U"),
         validity = function(object) {
             object. <- new("dtpMatrix")
             object.@Dim <- object@Dim
             object.@Dimnames <- object@Dimnames
             object.@uplo <- object@uplo
             object.@x <- object@x
             if(is.character(valid <- validObject(object., test = TRUE)))
                 valid
             else .Call(pBunchKaufman_validate, object)
         })
}


## ------ Cholesky -----------------------------------------------------

if(.M2) {
setClass("Cholesky",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseCholesky",
         contains = c("VIRTUAL", "Cholesky"),
         slots = c(x = "vector", perm = "integer", uplo = "character"),
         prototype = list(uplo = "U"),
         validity = function(object) .Call(denseCholesky_validate, object))

setClass("ddenseCholesky",
         contains = "denseCholesky",
         slots = c(x = "numeric"))

setClass("zdenseCholesky",
         contains = "denseCholesky",
         slots = c(x = "complex"))

setClass("sparseCholesky",
         contains = c("VIRTUAL", "Cholesky"),
         slots = c(perm = "integer", colcount = "integer",
                   i = "integer", x = "vector", ordering = "integer"),
         prototype = list(ordering = 0L),
         validity = function(object) .Call(sparseCholesky_validate, object))

setClass("dsparseCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(x = "numeric"))

setClass("zsparseCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(x = "complex"))

setClass("simplicialCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(p = "integer", nz = "integer",
                   prev = "integer", `next` = "integer",
                   is.ll = "logical", is.monotonic = "logical"),
         prototype = list(p = 0L, prev = c(1L, -1L), `next` = c(-1L, 0L),
                          is.ll = TRUE, is.monotonic = TRUE),
         validity = function(object) .Call(simplicialCholesky_validate, object))

setClass("supernodalCholesky",
         contains = c("VIRTUAL", "sparseCholesky"),
         slots = c(super = "integer", pi = "integer", px = "integer",
                   maxcsize = "integer", maxesize = "integer"),
         prototype = list(super = 0L, pi = 0L, px = 0L,
                          maxcsize = 0L, maxesize = 0L),
         validity = function(object) .Call(supernodalCholesky_validate, object))

setClass("dsimplicialCholesky",
         contains = c("dsparseCholesky", "simplicialCholesky"))

setClass("zsimplicialCholesky",
         contains = c("zsparseCholesky", "simplicialCholesky"))

setClass("dsupernodalCholesky",
         contains = c("dsparseCholesky", "supernodalCholesky"))

setClass("zsupernodalCholesky",
         contains = c("zsparseCholesky", "supernodalCholesky"))
} else {
setClass("CholeskyFactorization",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("Cholesky",
         contains = "CholeskyFactorization",
         slots = c(uplo = "character", x = "numeric", perm = "integer"),
         prototype = list(uplo = "U"),
         validity = function(object) {
             object. <- new("dtrMatrix")
             object.@Dim <- object@Dim
             object.@Dimnames <- object@Dimnames
             object.@uplo <- object@uplo
             object.@x <- object@x
             if(is.character(valid <- validObject(object., test = TRUE)))
                 valid
             else .Call(Cholesky_validate, object)
         })

setClass("pCholesky",
         contains = "CholeskyFactorization",
         slots = c(uplo = "character", x = "numeric", perm = "integer"),
         prototype = list(uplo = "U"),
         validity = function(object) {
             object. <- new("dtpMatrix")
             object.@Dim <- object@Dim
             object.@Dimnames <- object@Dimnames
             object.@uplo <- object@uplo
             object.@x <- object@x
             if(is.character(valid <- validObject(object., test = TRUE)))
                 valid
             else .Call(pCholesky_validate, object)
         })

setClass("CHMfactor",
         contains = c("VIRTUAL", "CholeskyFactorization"),
         slots = c(type = "integer", colcount = "integer", perm = "integer"),
         validity = function(object) .Call(CHMfactor_validate, object))

setClass("CHMsimpl",
         contains = c("VIRTUAL", "CHMfactor"),
         slots = c(p = "integer", i = "integer", nz = "integer",
                   nxt = "integer", prv = "integer"),
         prototype = list(type = c(0L, 1L, 0L, 1L, 0L, 0L),
                          p = 0L, nxt = c(-1L, 0L), prv = c(1L, -1L)),
         validity = function(object) .Call(CHMsimpl_validate, object))

setClass("nCHMsimpl",
         contains = "CHMsimpl")
setClass("dCHMsimpl",
         contains = "CHMsimpl",
         slots = c(x = "numeric"),
         validity = function(object) .Call(dCHMsimpl_validate, object))

setClass("CHMsuper",
         contains = c("VIRTUAL", "CHMfactor"),
         slots = c(super = "integer", pi = "integer", px = "integer",
                   s = "integer"),
         prototype = list(type = c(0L, 1L, 1L, 1L, 0L, 0L),
                          super = 0L, pi = 0L, px = 0L),
         validity = function(object) .Call(CHMsuper_validate, object))

setClass("nCHMsuper",
         contains = "CHMsuper")
setClass("dCHMsuper",
         contains = "CHMsuper",
         slots = c(x = "numeric"),
         validity = function(object) .Call(dCHMsuper_validate, object))
}


## ------ Schur --------------------------------------------------------

if(.M2) {
setClass("Schur",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("denseSchur",
         contains = c("VIRTUAL", "Schur"),
         slots = c(x = "vector", vectors = "vector", values = "vector"),
         prototype = list(values = double(0L)),
         validity = function(object) .Call(denseSchur_validate, object))

setClass("ddenseSchur",
         contains = "SchurFactorization",
         slots = c(x = "numeric", vectors = "numeric"))

setClass("zdenseSchur",
         contains = "SchurFactorization",
         slots = c(x = "complex", vectors = "complex"))
} else {
setClass("SchurFactorization",
         contains = c("VIRTUAL", "MatrixFactorization"))

setClass("Schur",
         contains = "SchurFactorization",
         slots = c(Q = "Matrix", T = "Matrix", EValues = "vector"),
         prototype = list(Q = .new("dgeMatrix"), T = .new("dgeMatrix")),
         validity = function(object) .Call(Schur_validate, object))
}


########################################################################
##  3. sparseVector
########################################################################

## ------ The Mother Class 'sparseVector' ------------------------------

setClass("sparseVector",
         contains = "VIRTUAL",
         slots = c(length = "numeric", i = "numeric"),
         prototype = list(length = 0),
         validity = function(object) .Call(sparseVector_validate, object))

## Allow users to do new("[nlidz]sparseVector", i=, x=) with unsorted 'i'
setMethod("initialize", c(.Object = "sparseVector"),
          function(.Object, i, x, ...) {
              if(has.x <- !missing(x))
                  x <- x # MJ: why is this necessary?
              if(!missing(i)) {
                  i.uns <- is.unsorted(i, strictly = TRUE)
                  i <-
                      if(is.na(i.uns) || !i.uns)
                          i
                      else {
                          ## we know that there are no NA, and the order of
                          ## ties does not matter (since ties are an error),
                          ## hence it is safe to use "quick" here
                          m <- if(is.integer(length(i))) "radix" else "quick"
                          if(.hasSlot(.Object, "x") && has.x) {
                              s <- sort.int(i, method = m, index.return = TRUE)
                              x <- x[s$ix]
                              s$x
                          } else sort.int(i, method = m)
                      }
              }
              callNextMethod()
          })


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
##  4. Index and more "miscellaneous" classes, but _not_ class unions
########################################################################

## Idea: represent x = c(seq(from1, to1, by1), seq(from2, to2, by2), ...)
##       as list(first = x[1L], rle = rle(diff(x)))
setOldClass("rle")
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

## MJ: aim to deprecate and eventually remove these ...

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
setClassUnion("replValueSp",
              members = c("replValue", "sparseVector", "matrix", "Matrix"))

rm(.new, .initialize)
