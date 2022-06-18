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

## MJ: now inherited from ANY
if(FALSE) {
setMethod("determinant", signature(x = "dtpMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) mkDet(diag(x), TRUE))

setMethod("norm", signature(x = "dtpMatrix", type = "missing"),
	  function(x, type, ...) .Call(dtpMatrix_norm, x, "O"))

setMethod("rcond", signature(x = "dtpMatrix", norm = "missing"),
	  function(x, norm, ...) .Call(dtpMatrix_rcond, x, "O"))
} ## MJ

setMethod("determinant", signature(x = "dtpMatrix", logarithm = "logical"),
	  function(x, logarithm, ...) mkDet(diag(x), logarithm))

setMethod("norm", signature(x = "dtpMatrix", type = "character"),
	  function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dtpMatrix_norm, x, type))

setMethod("rcond", signature(x = "dtpMatrix", norm = "character"),
	  function(x, norm, ...) .Call(dtpMatrix_rcond, x, norm))

setMethod("solve", signature(a = "dtpMatrix", b = "missing"),
	  function(a, b, ...) .Call(dtpMatrix_solve, a))

setMethod("solve", signature(a = "dtpMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dtpMatrix_matrix_solve, a, as(b, "denseMatrix")))

setMethod("solve", signature(a = "dtpMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b))

setMethod("solve", signature(a = "dtpMatrix", b = "numLike"),
	  function(a, b, ...) .Call(dtpMatrix_matrix_solve, a, b))

.is.na <- .is.infinite <- function(x) {
    n <- (d <- x@Dim)[1L]
    i <- logical(n^2)
    if(n > 0L) {
        uploU <- (uplo <- x@uplo) == "U"
        diagN <- x@diag == "N"
        i[indTri(n, uploU, TRUE)] <- i.t <- is.na(x@x)
        if(!diagN) {
            i[indDiag(n)] <- FALSE
            i.t <- i[indTri(n, uploU, FALSE)]
        }
        if(any(i.t))
            return(new("ntpMatrix", Dim = d, Dimnames = x@Dimnames,
                       uplo = uplo, diag = "N", x = i))
    }
    is.na_nsp(x)
}
body(.is.infinite) <-
    do.call(substitute, list(body(.is.infinite),
                             list(is.na = quote(is.infinite))))

.is.finite <- function(x) {
    n <- (d <- x@Dim)[1L]
    i <- rep.int(TRUE, n^2)
    if(n > 0L) {
        uploU <- x@uplo == "U"
        diagN <- x@diag == "N"
        i[indTri(n, uploU, TRUE)] <- is.finite(x@x)
        if(!diagN)
            i[indDiag(n)] <- TRUE
    }
    new("ngeMatrix", Dim = d, Dimnames = x@Dimnames, x = i)
}

for(.cl in paste0(c("d", "l"), "tpMatrix"))
    setMethod("is.na", signature(x = .cl), .is.na)
setMethod("is.infinite", signature(x = "dtpMatrix"), .is.infinite)
setMethod("is.finite", signature(x = "dtpMatrix"), .is.finite)

rm(.cl, .is.na, .is.infinite, .is.finite)
