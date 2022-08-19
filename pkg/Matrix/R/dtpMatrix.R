## METHODS FOR CLASS: dtpMatrix
## dense (packed) triangular matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("dtpMatrix", "dtrMatrix",
      dtp2dtr <- function(from) .Call(dtpMatrix_as_dtrMatrix, from))
setAs("matrix", "dtpMatrix",
      function(from) as(as(from, "dtrMatrix"), "dtpMatrix"))
setAs("dtpMatrix", "matrix",
      function(from) as(dtp2dtr(from), "matrix"))
setAs("pCholesky", "lMatrix",
      function(from) as(as(from, "dtpMatrix"), "lMatrix"))
setAs("pBunchKaufman", "lMatrix",
      function(from) as(as(from, "dtpMatrix"), "lMatrix"))
} ## MJ

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
## Is this needed?  already have coercion to "TsparseMatrix" {FIXME}
setAs("dtpMatrix", "dtTMatrix",
      function(from) {
	  x <- as(from, "TsparseMatrix")
          cld <- getClassDef(class(x))
	  if(extends(cld, "dtTMatrix"))
	      x
	  else { ## triangularity lost: should not have happened
	      warning("inefficient coercion (lost triangularity); please report")
	      gT2tT(as(x, "dgTMatrix"), uplo = from@uplo, diag = from@diag,
		    toClass = "dtTMatrix", do.n = FALSE)
	  }
      })
} ## MJ

## MJ: no longer needed ... replacement in ./packedMatrix.R
if (FALSE) {
setMethod("t", "dtpMatrix",
	  function(x) dtr2dtp(t(dtp2dtr(x))), valueClass = "dtpMatrix")
setMethod("diag", signature(x = "dtpMatrix"),
	  function(x, nrow, ncol) .Call(dtpMatrix_getDiag, x),
	  valueClass = "numeric")
setMethod("diag<-", signature(x = "dtpMatrix"),
	  function(x, value) {
	      .Call(dtpMatrix_setDiag,
		    if(x@diag == "U") .dense.diagU2N(x, "d", isPacked=TRUE) else x,
		    value)
	  })
} ## MJ

setMethod("solve", signature(a = "dtpMatrix", b = "missing"),
	  function(a, b, ...) .Call(dtpMatrix_solve, a))

setMethod("solve", signature(a = "dtpMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dtpMatrix_matrix_solve, a, as(b, "denseMatrix")))

setMethod("solve", signature(a = "dtpMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b))

setMethod("solve", signature(a = "dtpMatrix", b = "numLike"),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b))
