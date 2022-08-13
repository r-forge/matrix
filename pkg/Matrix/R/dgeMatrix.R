## METHODS FOR CLASS: dgeMatrix
## dense general matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./denseMatrix.R
if(FALSE) {
## ..2dge() -> ./Auxiliaries.R
setAs("matrix",  "dgeMatrix", function(from) ..2dge(from))
setAs("numLike", "dgeMatrix", function(from) ..2dge(from))

ge2mat <- function(from) array(from@x, dim = from@Dim, dimnames = from@Dimnames)
setAs("dgeMatrix", "matrix", ge2mat)

setMethod("as.vector", "dgeMatrix",
          function(x, mode) as.vector(x@x, mode))
} ## MJ

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if (FALSE) {
..get.diag <- function(x, nrow, ncol, names=TRUE) {
    ##         vvvvvvvvvvvvvvvvv here just a place holder, replaced in .mkSpec.diag()
    y <- .Call(dgeMatrix_getDiag, x) # double or logical
    if(names) {
        nms <- dimnames(x)
        if(is.list(nms) && !any(vapply(nms, is.null, NA)) &&
           identical((nm <- nms[[1L]][im <- seq_len(min(dim(x)))]), nms[[2L]][im]))
            names(y) <- nm
    }
    y
}
.mkSpec.diag <- function(symb) {
    rr <- ..get.diag
    body(rr)[[2]][[3]][[2]] <- symb
    rr
}
.dge.diag <- .mkSpec.diag(quote(dgeMatrix_getDiag))

setMethod("t", signature(x = "dgeMatrix"), t_geMatrix)
setMethod("diag", signature(x = "dgeMatrix"), .dge.diag)
setMethod("diag<-", signature(x = "dgeMatrix"),
	  function(x, value) .Call(dgeMatrix_setDiag, x, value))
} ## MJ

## MJ: now inherited from ANY
if(FALSE) {
setMethod("norm", signature(x = "dgeMatrix", type = "missing"),
	  function(x, type, ...) .Call(dgeMatrix_norm, x, "O"))

setMethod("rcond", signature(x = "dgeMatrix", norm = "missing"),
	  function(x, norm, ...) rcond(x, norm = "O", ...))
} ## MJ

## MJ: redundant as ddenseMatrix has the same
if(FALSE) {
setMethod("chol", signature(x = "dgeMatrix"), cholMat)
} ## MJ

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, warnSing = TRUE, ...) .Call(dgeMatrix_LU, x, warnSing))

setMethod("norm", signature(x = "dgeMatrix", type = "character"),
	  function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dgeMatrix_norm, x, type))

setMethod("rcond", signature(x = "dgeMatrix", norm = "character"),
	  function(x, norm, ...) {
              d <- x@Dim
	      if(d[1L] != d[2L])
		  rcond(qr.R(qr(if(d[1L] < d[2L]) t(x) else x)),
                        norm = norm, ...)
              else .Call(dgeMatrix_rcond, x, norm)
	  })

setMethod("solve", signature(a = "dgeMatrix", b = "missing"),
	  function(a, b, ...) .Call(dgeMatrix_solve, a))

setMethod("solve", signature(a = "dgeMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dgeMatrix_matrix_solve, a, as(b, "denseMatrix")))

setMethod("solve", signature(a = "dgeMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b))

setMethod("solve", signature(a = "dgeMatrix", b = "numLike"),
	  function(a, b, ...) .Call(dgeMatrix_matrix_solve, a, b))

.is.na <- .is.infinite <- .is.finite <- function(x) {
    if(any(i <- is.na(x@x)))
        new("ngeMatrix", Dim = x@Dim, Dimnames = x@Dimnames, x = i)
    else is.na_nsp(x)
}
body(.is.infinite) <-
    do.call(substitute, list(body(.is.infinite),
                             list(is.na = quote(is.infinite))))
body(.is.finite) <-
    do.call(substitute, list(body(.is.finite),
                             list(is.na = quote(is.finite))))

for(.cl in paste0(c("d", "l"), "geMatrix"))
    setMethod("is.na", signature(x = .cl), .is.na)
setMethod("is.infinite", signature(x = "dgeMatrix"), .is.infinite)
setMethod("is.finite", signature(x = "dgeMatrix"), .is.finite)

rm(.cl, .is.na, .is.infinite, .is.finite)
