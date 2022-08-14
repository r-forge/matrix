## METHODS FOR CLASS: dtrMatrix
## dense (unpacked) triangular matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
setAs("dgeMatrix", "dtrMatrix", function(from) asTri(from, "dtrMatrix"))
setAs("dtrMatrix", "dtpMatrix",
      dtr2dtp <- function(from) .Call(dtrMatrix_as_dtpMatrix, from))
setAs("matrix", "dtrMatrix",
      function(from) as(..2dge(from), "dtrMatrix"))

.dtr2mat <- function(from, keep.dimnames=TRUE)
    .Call(dtrMatrix_as_matrix, from, keep.dimnames)
## needed for t() method
setAs("dtrMatrix", "matrix",
      function(from) .Call(dtrMatrix_as_matrix, from, TRUE))
setAs("Cholesky", "lMatrix",
      function(from) as(as(from, "dtrMatrix"), "lMatrix"))
setAs("BunchKaufman", "lMatrix",
      function(from) as(as(from, "dtrMatrix"), "lMatrix"))
setAs("dtrMatrix", "sparseMatrix",
      function(from) .dense2C(from, kind = "tri", uplo = from@uplo))
setAs("dtrMatrix", "CsparseMatrix",
      function(from) .dense2C(from, kind = "tri", uplo = from@uplo))
} ## MJ

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if (FALSE) {
setMethod("t", signature(x = "dtrMatrix"), t_trMatrix)
setMethod("diag", signature(x = "dtrMatrix"),
          .mkSpec.diag(quote(dtrMatrix_getDiag)),
          valueClass = "numeric")
setMethod("diag<-", signature(x = "dtrMatrix"),
	  function(x, value) {
	      .Call(dtrMatrix_setDiag,
		    if(x@diag == "U") .dense.diagU2N(x, "d", isPacked=FALSE) else x,
		    value)
	  })
} ## MJ

## MJ: now inherited from ANY
if(FALSE) {
setMethod("norm", signature(x = "dtrMatrix", type = "missing"),
	  function(x, type, ...) .Call(dtrMatrix_norm, x, "O"))

setMethod("rcond", signature(x = "dtrMatrix", norm = "missing"),
	  function(x, norm, ...) .Call(dtrMatrix_rcond, x, "O"))
} ## MJ

setMethod("norm", signature(x = "dtrMatrix", type = "character"),
	  function(x, type, ...) {
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dtrMatrix_norm, x, type)
          })

setMethod("rcond", signature(x = "dtrMatrix", norm = "character"),
	  function(x, norm, ...) .Call(dtrMatrix_rcond, x, norm))

setMethod("chol2inv", signature(x = "dtrMatrix"),
	  function (x, ...) {
	      if(x@diag != "N")
                  x <- .dense.diagU2N(x)
	      .Call(dtrMatrix_chol2inv, x)
	  })

setMethod("solve", signature(a = "dtrMatrix", b = "missing"),
	  function(a, b, ...) .Call(dtrMatrix_solve, a))

setMethod("solve", signature(a = "dtrMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dtrMatrix_matrix_solve, a, as(b, "denseMatrix")))

setMethod("solve", signature(a = "dtrMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dtrMatrix_matrix_solve, a, b))

setMethod("solve", signature(a = "dtrMatrix", b = "numLike"),
	  function(a, b, ...) .Call(dtrMatrix_matrix_solve, a, b))

.is.na <- .is.infinite <- function(x) {
    n <- (d <- x@Dim)[1L]
    i <- is.na(x@x)
    i[indTri(n, (uplo <- x@uplo) != "U", x@diag != "N")] <- FALSE
    if(any(i))
        new("ntrMatrix", Dim = d, Dimnames = x@Dimnames,
            uplo = uplo, diag = "N", x = i)
    else is.na_nsp(x)
}
body(.is.infinite) <-
    do.call(substitute, list(body(.is.infinite),
                             list(is.na = quote(is.infinite))))

.is.finite <- function(x) {
    n <- (d <- x@Dim)[1L]
    i <- is.finite(x@x)
    i[indTri(n, x@uplo != "U", x@diag != "N")] <- TRUE
    new("ngeMatrix", Dim = d, Dimnames = x@Dimnames, x = i)
}

for(.cl in paste0(c("d", "l"), "trMatrix"))
    setMethod("is.na", signature(x = .cl), .is.na)
setMethod("is.infinite", signature(x = "dtrMatrix"), .is.infinite)
setMethod("is.finite", signature(x = "dtrMatrix"), .is.finite)

rm(.cl, .is.na, .is.infinite, .is.finite)
