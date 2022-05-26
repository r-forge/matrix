#### Positive-definite Symmetric Matrices -- Coercion and Methods

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsy2dpo <- function(from) {
    if(is.null(tryCatch(.Call(dpoMatrix_chol, from), error = function(e) NULL)))
        stop("not a positive definite matrix")
    ## FIXME: check=FALSE
    copyClass(from, "dpoMatrix",
              sNames = c("Dim", "Dimnames", "uplo", "x", "factors"))
}
.M2dpo.via.dsy <- function(from) {
    .dsy2dpo(as(from, "dsyMatrix"))
}
.M2dpo.via.virtual <- function(from) {
    ## still needs as(<ds[py]Matrix>, "dpoMatrix") to work
    as(as(as(as(from,"symmetricMatrix"),"dMatrix"),"denseMatrix"),"dpoMatrix")
}

setAs("dsyMatrix", "dpoMatrix", .dsy2dpo)
setAs("dspMatrix", "dpoMatrix", .M2dpo.via.dsy)
setAs(   "matrix", "dpoMatrix", .M2dpo.via.dsy)
setAs(   "Matrix", "dpoMatrix", .M2dpo.via.virtual)

rm(.M2dpo.via.dsy, .M2dpo.via.virtual)


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

..pack <- function(from) pack(from)
setAs("dpoMatrix", "dppMatrix", ..pack)
setAs("dpoMatrix", "dspMatrix", ..pack)
rm(..pack)

## MJ: no longer needed ... prefer variant above going via pack()
if(FALSE) {
setAs("dpoMatrix", "dppMatrix",
      function(from) {
          ## FIXME: check=FALSE
          copyClass(.Call(dsyMatrix_as_dspMatrix, from), "dppMatrix",
                    sNames = c("x", "Dim", "Dimnames", "uplo", "factors"))
      })
} ## MJ

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("dpoMatrix", "lMatrix",
      function(from) as(as(from, "dsyMatrix"), "lMatrix"))
setAs("dpoMatrix", "nMatrix",
      function(from) as(as(from, "dsyMatrix"), "nMatrix"))
} ## MJ


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", signature(x = "dpoMatrix"),
	  function(x, pivot, ...) .Call(dpoMatrix_chol, x))

setMethod("rcond", signature(x = "dpoMatrix", norm = "character"),
          function(x, norm, ...) .Call(dpoMatrix_rcond, x, norm))

setMethod("rcond", signature(x = "dpoMatrix", norm = "missing"),
          function(x, norm, ...) .Call(dpoMatrix_rcond, x, "O"))

setMethod("solve", signature(a = "dpoMatrix", b = "missing"),
          function(a, b, ...) .Call(dpoMatrix_solve, a),
          valueClass = "dpoMatrix")

setMethod("solve", signature(a = "dpoMatrix", b = "dgeMatrix"),
          function(a, b, ...) .Call(dpoMatrix_dgeMatrix_solve, a, b),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dpoMatrix", b = "matrix"),
          function(a, b, ...) .Call(dpoMatrix_matrix_solve, a, b),
          valueClass = "matrix")

mkDet.via.chol <- function(x, logarithm, ...)
    mkDet(logarithm, ldet = 2*sum(log(abs(diag(chol(x))))), sig = 1L)

setMethod("determinant", signature(x = "dpoMatrix", logarithm = "logical"), mkDet.via.chol)
setMethod("determinant", signature(x = "dpoMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) mkDet.via.chol(x, logarithm=TRUE))

## Is this usable / necessary?  -- FIXME!
## setMethod("solve", signature(a = "dpoMatrix", b = "numeric"),
##          function(a, b, ...)
##          as.numeric(.Call(dpoMatrix_matrix_solve,
##                           a, as.matrix(b))),
##          valueClass = "numeric")
