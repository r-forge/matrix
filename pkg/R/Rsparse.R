#### Sparse Matrices in Compressed row-oriented format
####                               --- "R"

### ``mainly for completeness'' --- we *do* favour Csparse
##    - - - - - - - - - - - -   hence only "minimal" methods here !
##  see also ./SparseM-conv.R

### contains = "dMatrix"

## compressed_to_TMatrix -- fails on 32bit--enable-R-shlib with segfault {Kurt}
## ------------ --> ../src/dgCMatrix.c
##_SF_ .R.2.T <- function(from) .Call(compressed_to_TMatrix, from, FALSE)
## slow R-level workaround
## this is cheap; alternative: going there directly, using
##	i <- .Call(Matrix_expand_pointers, from@p),
.R.2.T <- function(from) as(.R.2.C(from), "TsparseMatrix")

## R_to_CMatrix -- fails on 32bit--enable-R-shlib with segfault {Kurt}
## ------------ --> ../src/dgCMatrix.c
##_SF_ .R.2.C <- function(from) .Call(R_to_CMatrix, from)
## "slow" R-level workaround
.R.2.C <- function(from)
{
    cl <- class(from)
    valid <- c("dgRMatrix", "dsRMatrix", "dtRMatrix",
               "lgRMatrix", "lsRMatrix", "ltRMatrix",
               "ngRMatrix", "nsRMatrix", "ntRMatrix",
               "zgRMatrix", "zsRMatrix", "ztRMatrix")
    icl <- match(cl, valid) - 1:1
    if(is.na(icl)) stop("invalid class:", cl)
    Ccl <- sub("^(..)R","\\1C", cl)  # corresponding Csparse class name
    r <- new(Ccl)
    r@Dim <- rev(from@Dim)
    if(icl %/% 3 != 2) ## not "n..Matrix" --> has 'x' slot
        r@x <- from@x
    if(icl %% 3 != 0) {                 # symmetric or triangular
        r@uplo <- from@uplo
        if(icl %% 3 == 2)               # triangular
            r@diag <- from@diag
    }
    r@i <- from@j
    r@p <- from@p
    r <- t(r)
    r@Dimnames <- from@Dimnames
    r
}

## coercion to other virtual classes --- the functionality we want to encourage

setAs("RsparseMatrix", "TsparseMatrix", .R.2.T)
setAs("RsparseMatrix", "CsparseMatrix", .R.2.C)

setAs("RsparseMatrix", "denseMatrix",
      function(from) as(.R.2.C(from), "denseMatrix"))

setAs("RsparseMatrix", "dsparseMatrix",
      function(from) as(.R.2.C(from), "dsparseMatrix"))
##_SF_       function(from) as(.Call(R_to_CMatrix, from), "dsparseMatrix"))
setAs("RsparseMatrix", "lsparseMatrix",
      function(from) as(.R.2.C(from), "lsparseMatrix"))
##_SF_      function(from) as(.Call(R_to_CMatrix, from), "lsparseMatrix"))
setAs("RsparseMatrix", "nsparseMatrix",
      function(from) as(.R.2.C(from), "nsparseMatrix"))
##_SF_      function(from) as(.Call(R_to_CMatrix, from), "nsparseMatrix"))

setAs("RsparseMatrix", "dMatrix",
      function(from) as(.R.2.C(from), "dMatrix"))
##_SF_      function(from) as(.Call(R_to_CMatrix, from), "dMatrix"))
setAs("RsparseMatrix", "lMatrix",
      function(from) as(.R.2.C(from), "lMatrix"))
##_SF_      function(from) as(.Call(R_to_CMatrix, from), "lMatrix"))
setAs("RsparseMatrix", "nMatrix",
      function(from) as(.R.2.C(from), "nMatrix"))
##_SF_      function(from) as(.Call(R_to_CMatrix, from), "nMatrix"))


## for printing etc:
setAs("RsparseMatrix", "dgeMatrix",
      function(from) as(.R.2.C(from), "dgeMatrix"))
setAs("RsparseMatrix", "matrix",
      function(from) as(.R.2.C(from), "matrix"))


##--- and all these are just "the essential low-level coercions" : ----------

## setAs("dgRMatrix", "matrix",
##       function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))
## setAs("lgRMatrix", "matrix",
##       function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))
## setAs("ngRMatrix", "matrix",
##       function(from) as(.Call(compressed_to_TMatrix, from, FALSE), "matrix"))

## setAs("dgRMatrix", "dgeMatrix",
##       function(from) as(.R.2.C(from), "dgeMatrix"))
## setAs("lgRMatrix", "lgeMatrix",
##       function(from) as(.R.2.C(from), "lgeMatrix"))
## setAs("ngRMatrix", "ngeMatrix",
##       function(from) as(.R.2.C(from), "ngeMatrix"))

## setAs("dgRMatrix", "dgCMatrix", .R.2.C)
## setAs("lgRMatrix", "lgCMatrix", .R.2.C)
## setAs("ngRMatrix", "ngCMatrix", .R.2.C)
## ## really needed? :
## setAs("dgRMatrix", "CsparseMatrix", .R.2.C)


## setAs("dgRMatrix", "dgTMatrix", .R.2.T)
## setAs("lgRMatrix", "lgTMatrix", .R.2.T)
## setAs("ngRMatrix", "ngTMatrix", .R.2.T)

##=== Now the same stories for the "s" (symmetric) and "t" (triangular) ones ===

## setAs("dsRMatrix", "dsCMatrix", .R.2.C)
## setAs("lsRMatrix", "lsCMatrix", .R.2.C)
## setAs("nsRMatrix", "nsCMatrix", .R.2.C)

## setAs("dsRMatrix", "dsTMatrix", .R.2.T)
## setAs("lsRMatrix", "lsTMatrix", .R.2.T)
## setAs("nsRMatrix", "nsTMatrix", .R.2.T)

## setAs("dsRMatrix", "dsyMatrix",
##       function(from) as(.R.2.C(from), "dsyMatrix"))
## setAs("lsRMatrix", "lsyMatrix",
##       function(from) as(.R.2.C(from), "lsyMatrix"))
## setAs("nsRMatrix", "nsyMatrix",
##       function(from) as(.R.2.C(from), "nsyMatrix"))

## setAs("dtRMatrix", "dtCMatrix", .R.2.C)
## setAs("ltRMatrix", "ltCMatrix", .R.2.C)
## setAs("ntRMatrix", "ntCMatrix", .R.2.C)

## setAs("dtRMatrix", "dtTMatrix", .R.2.T)
## setAs("ltRMatrix", "ltTMatrix", .R.2.T)
## setAs("ntRMatrix", "ntTMatrix", .R.2.T)

## setAs("dtRMatrix", "dtrMatrix",
##       function(from) as(.R.2.C(from), "dtrMatrix"))
## setAs("ltRMatrix", "ltrMatrix",
##       function(from) as(.R.2.C(from), "ltrMatrix"))
## setAs("ntRMatrix", "ntrMatrix",
##       function(from) as(.R.2.C(from), "ntrMatrix"))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })

## **VERY** cheap substitutes:  work via dgC and t(.)
.viaC.to.dgR <- function(from) {
    m <- as(t(from), "dgCMatrix")
    new("dgRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
	p = m@p, j = m@i, x = m@x)
}

setAs("matrix",    "dgRMatrix", .viaC.to.dgR)## one of the few coercions "to specific"
setAs("matrix",    "RsparseMatrix", .viaC.to.dgR)
setAs("ddenseMatrix", "RsparseMatrix", .viaC.to.dgR)
setAs("dsparseMatrix", "RsparseMatrix", .viaC.to.dgR)

## symmetric: can use same 'p' slot
setAs("dsCMatrix", "dsRMatrix",
      function(from) new("dsRMatrix", Dim = dim(from), Dimnames = .M.DN(from),
	      p = from@p, j = from@i, x = from@x,
	      uplo = if (from@uplo == "U") "L" else "U"))
## FIXME: if this makes sense, do it for "l" and "n" as well as "d"

## setAs("dtCMatrix", "dtRMatrix", .viaC.to.dgR) # should work; can NOT use 'p'


##setAs("dgRMatrix", "dgeMatrix",
##      function(from) .Call(csc_to_dgeMatrix, from))

##setAs("matrix", "dgRMatrix",
##      function(from) {
##          storage.mode(from) <- "double"
##          .Call(matrix_to_csc, from)
##      })


##setMethod("diag", signature(x = "dgRMatrix"),
##          function(x = 1, nrow, ncol = n) .Call(csc_getDiag, x))

## try to define for "Matrix" -- once and for all -- but that fails -- why? __ FIXME __
## setMethod("dim", signature(x = "dgRMatrix"),
##           function(x) x@Dim, valueClass = "integer")

##setMethod("t", signature(x = "dgRMatrix"),
##          function(x) .Call(csc_transpose, x),
##          valueClass = "dgRMatrix")

setMethod("image", "dgRMatrix",
          function(x, ...) {
              x <- as(x, "TsparseMatrix")
              callGeneric()
          })

setMethod("t", "RsparseMatrix", function(x) as(t(.R.2.T(x)), "RsparseMatrix"))


## Want tril(), triu(), band() --- just as "indexing" ---
## return a "close" class:
setMethod("tril", "RsparseMatrix",
	  function(x, k = 0, ...)
	  as(tril(.R.2.C(x), k = k, ...), "RsparseMatrix"))
setMethod("triu", "RsparseMatrix",
	  function(x, k = 0, ...)
	  as(triu(.R.2.C(x), k = k, ...), "RsparseMatrix"))
setMethod("band", "RsparseMatrix",
	  function(x, k1, k2, ...)
	  as(band(.R.2.C(x), k1 = k1, k2 = k2, ...), "RsparseMatrix"))
