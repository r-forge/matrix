## METHODS FOR CLASS: dpoMatrix
## dense (unpacked) symmetric positive definite matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dsy2dpo <- function(from) {
    if(is.null(tryCatch(.Call(dpoMatrix_chol, from), error = function(e) NULL)))
        stop("not a positive definite matrix")
    ## FIXME: check=FALSE
    copyClass(from, "dpoMatrix",
              sNames = c("Dim", "Dimnames", "uplo", "x", "factors"))
}

setAs("dsyMatrix", "dpoMatrix", .dsy2dpo)

setAs("dspMatrix", "dpoMatrix",
      function(from) unpack(.dsp2dpp(from)))

setAs("matrix", "dpoMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .dsy2dpo(.M2symm(from))
      })

setAs("Matrix", "dpoMatrix",
      function(from) {
          ## still needs as(<ds[yp]Matrix>, "dpoMatrix") to work
          as(as(as(as(from,"symmetricMatrix"),"dMatrix"),"denseMatrix"),
             "dpoMatrix")
      })


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("dpoMatrix", "dppMatrix", function(from) pack(from))

## MJ: no longer needed ... prefer above
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

## MJ: now inherited from ANY
if(FALSE) {
setMethod("rcond", signature(x = "dpoMatrix", norm = "missing"),
          function(x, norm, ...) .Call(dpoMatrix_rcond, x, "O"))
} ## MJ

setMethod("chol", signature(x = "dpoMatrix"),
	  function(x, ...) .Call(dpoMatrix_chol, x))

setMethod("rcond", signature(x = "dpoMatrix", norm = "character"),
          function(x, norm, ...) .Call(dpoMatrix_rcond, x, norm))

setMethod("solve", signature(a = "dpoMatrix", b = "missing"),
          function(a, b, ...) .Call(dpoMatrix_solve, a))

setMethod("solve", signature(a = "dpoMatrix", b = "Matrix"),
          function(a, b, ...)
              .Call(dpoMatrix_matrix_solve, a, as(b, "denseMatrix")))

setMethod("solve", signature(a = "dpoMatrix", b = "matrix"),
          function(a, b, ...) .Call(dpoMatrix_matrix_solve, a, b))

setMethod("solve", signature(a = "dpoMatrix", b = "numLike"),
          function(a, b, ...) .Call(dpoMatrix_matrix_solve, a, b))
