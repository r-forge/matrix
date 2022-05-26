#### Positive-definite Symmetric Packed Matrices -- Coercion and Methods

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsp2dpp <- function(from) {
    if(is.null(tryCatch(.Call(dppMatrix_chol, from), error = function(e) NULL)))
        stop("not a positive definite matrix")
    ## FIXME: check=FALSE
    copyClass(from, "dppMatrix",
              sNames = c("Dim", "Dimnames", "uplo", "x", "factors"))
}
.M2dpp.via.dsp <- function(from) {
    .dsp2dpp(as(from, "dspMatrix"))
}
.M2dpp.via.virtual <- function(from) {
    ## still needs as(<ds[py]Matrix>, "dppMatrix") to work
    as(as(as(as(from,"symmetricMatrix"),"dMatrix"),"denseMatrix"),"dppMatrix")
}

setAs("dspMatrix", "dppMatrix", .dsp2dpp)
setAs("dsyMatrix", "dppMatrix", .M2dpp.via.dsp)
setAs(   "matrix", "dppMatrix", .M2dpp.via.dsp)
setAs(   "Matrix", "dppMatrix", .M2dpp.via.virtual)

rm(.M2dpp.via.dsp, .M2dpp.via.virtual)


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

..unpack <- function(from) unpack(from)
setAs("dppMatrix", "dpoMatrix", ..unpack)
setAs("dppMatrix", "dsyMatrix", ..unpack)
rm(..unpack)

## MJ: no longer needed ... prefer variant above going via unpack()
if(FALSE) {
setAs("dppMatrix", "dpoMatrix",
      function(from) {
          ## FIXME: check=FALSE
          copyClass(.Call(dspMatrix_as_dsyMatrix, from),
                    "dpoMatrix",
                    sNames = c("x", "Dim", "Dimnames", "uplo", "factors"))
      })
} ## MJ

## MJ: redundant, as coercions are inherited from superclass dspMatrix
if(FALSE) {
dpp2sC <- function(from) as(.Call(dspMatrix_as_dsyMatrix, from), "dsCMatrix")
## setAs("dppMatrix", "dsCMatrix", dpp2sC)
setAs("dppMatrix", "CsparseMatrix", dpp2sC)
setAs("dppMatrix", "sparseMatrix", dpp2sC)
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
## (was infelicitous anyway because result did not have packed storage)
if(FALSE) {
setAs("dppMatrix", "lMatrix",
      function(from) as(as(from, "dsyMatrix"), "lMatrix"))
setAs("dppMatrix", "nMatrix",
      function(from) as(as(from, "dsyMatrix"), "nMatrix"))
} ## MJ


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", signature(x = "dppMatrix"),
	  function(x, pivot, LINPACK) .Call(dppMatrix_chol, x))

setMethod("determinant", signature(x = "dppMatrix", logarithm = "logical"), mkDet.via.chol)
setMethod("determinant", signature(x = "dppMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) mkDet.via.chol(x, logarithm=TRUE))

setMethod("rcond", signature(x = "dppMatrix", norm = "character"),
	  function(x, norm, ...)
	  .Call(dppMatrix_rcond, x, norm),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dppMatrix", norm = "missing"),
          function(x, norm, ...)
          .Call(dppMatrix_rcond, x, "O"),
          valueClass = "numeric")

setMethod("solve", signature(a = "dppMatrix", b = "missing"),
          function(a, b, ...)
          .Call(dppMatrix_solve, a),
          valueClass = "dppMatrix")

setMethod("solve", signature(a = "dppMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call(dppMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dppMatrix", b = "matrix"),
          function(a, b, ...)
          .Call(dppMatrix_matrix_solve, a, b),
          valueClass = "dgeMatrix")

##setMethod("solve", signature(a = "dppMatrix", b = "numeric"),
##          function(a, b, ...)
##          .Call(dppMatrix_matrix_solve, a, as.matrix(b)),
##          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dppMatrix", b = "integer"),
          function(a, b, ...) {
              storage.mode(b) <- "double"
              .Call(dppMatrix_matrix_solve, a, cbind(b, deparse.level=0L))
          }, valueClass = "dgeMatrix")

## MJ: no longer needed ... replacement in ./packedMatrix.R
if (FALSE) {
setMethod("t", signature(x = "dppMatrix"),
          function(x) as(t(as(x, "dspMatrix")), "dppMatrix"),
          valueClass = "dppMatrix")
}
