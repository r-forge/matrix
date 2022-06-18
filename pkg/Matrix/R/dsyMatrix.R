## METHODS FOR CLASS: dsyMatrix
## dense (unpacked) symmetric matrices with 'x' slot of type "double"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: no longer needed ... replacement in ./symmetricMatrix.R, ./denseMatrix.R
##     { where .dense2sy(), .dsy2mat() and .dsy2dsp(), which
##     have all been exported, are defined as simple aliases ... }
if(FALSE) {
.dense2sy <- function(from, ...) {
    if(isSymmetric(from, ...)) # < with tolerance!
	.Call(dense_to_symmetric, from, "U", FALSE)
    else
	stop("not a symmetric matrix; consider forceSymmetric() or symmpart()")
}
..dense2sy <- function(from) {
    ## NB: The alternative, 'zero tolerance' { <=> isSymmetric(*, tol=0) }
    ##     breaks too much previous code -- though it would be much faster --
    if(isSymmetric(from)) # < with tolerance!
	.Call(dense_to_symmetric, from, "U", FALSE)
    else
	stop("not a symmetric matrix; consider forceSymmetric() or symmpart()")
}
.dsy2mat <- function(from, keep.dimnames = TRUE) {
    .Call(dsyMatrix_as_matrix, from, keep.dimnames)
}
..dsy2mat <- function(from) {
    .Call(dsyMatrix_as_matrix, from, TRUE)
}
.dsy2dsp <- function(from) {
    .Call(dsyMatrix_as_dspMatrix, from)
}

setAs("dgeMatrix", "dsyMatrix", ..dense2sy)
setAs(   "matrix", "dsyMatrix", function(from) .dense2sy(..2dge(from)))
setAs("dsyMatrix",    "matrix", ..dsy2mat)
setAs("dsyMatrix", "dspMatrix", .dsy2dsp)

dsy2T <- function(from) { # 'dsT': only store upper *or* lower
    uplo <- from@uplo
    if(any0(dim(from))) {
	ij <- matrix(0L, 0,2) ; m <- from@x
    } else {
	## FIXME!	 working via "matrix" is *not* efficient:
	## the "other triangle" is filled, compared with 0, and then trashed:
	m <- .dense2m(from)
	ij <- which(m != 0, arr.ind = TRUE, useNames = FALSE)
	ij <- ij[if(uplo == "U") ij[,1] <= ij[,2] else ij[,1] >= ij[,2], , drop = FALSE]
    }
    new("dsTMatrix", i = ij[,1] - 1L, j = ij[,2] - 1L,
	x = as.vector(m[ij]), uplo = uplo,
	Dim = from@Dim, Dimnames = from@Dimnames)
}
dsy2C <- function(from) .T2Cmat(dsy2T(from), isTri=FALSE)

setAs("dsyMatrix", "dsTMatrix", dsy2T)
setAs("dsyMatrix", "dsCMatrix", dsy2C)
} ## MJ

## MJ: no longer needed ... replacement in ./unpackedMatrix.R
if(FALSE) {
.dsy.diag <- function(x, nrow, ncol, names=TRUE) {
    if(min(dim(x)) == 0L) return(numeric(0L))
    y <- .Call(dgeMatrix_getDiag, x)
    if(names) {
        nms <- dimnames(x) # rely on method for "symmetricMatrix" to symmetrize
        if(is.list(nms) && length(nms) == 2L)
            names(y) <- nms[[1L]]
    }
    y
}
## *Should* create the opposite storage format:  "U" -> "L"  and vice-versa:
setMethod("t", signature(x = "dsyMatrix"), t_trMatrix,
          valueClass = "dsyMatrix")
setMethod("diag", signature(x = "dsyMatrix"), .dsy.diag)
setMethod("diag<-", signature(x = "dsyMatrix"),
	  function(x, value) .Call(dgeMatrix_setDiag, x, value))
} ## MJ

## MJ: now inherited from ANY
if(FALSE) {
setMethod("norm", signature(x = "dsyMatrix", type = "missing"),
          function(x, type, ...) .Call(dsyMatrix_norm, x, "O"))

setMethod("rcond", signature(x = "dsyMatrix", norm = "missing"),
          function(x, norm, ...) .Call(dsyMatrix_rcond, x, "O"))
} ## MJ

## Now that we have "chol", we can define  "determinant" methods,
## exactly like in ./dsCMatrix.R  [[ also in ./dpoMatrix.R ]]
## DB - Probably figure out how to use the BunchKaufman decomposition instead
## {{FIXME: Shouldn't it be possible to have "determinant" work by
## default automatically for "Matrix"es  when there's a "chol" method available?
## ..> work with ss <- selectMethod("chol", signature("dgCMatrix"))
## -- not have to define showMethod("determinant", ...) for all classes

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
	  function(x, ...) .Call(dsyMatrix_trf, x))

setMethod("BunchKaufman", signature(x = "matrix"),
	  function(x, uplo = NULL, ...) .Call(matrix_trf, x, uplo))

## TODO: currently inherited from ddenseMatrix which goes via dgeMatrix
if(FALSE) {
setMethod("determinant", signature(x = "dsyMatrix", logarithm = "logical"),
          function(x, logarithm, ...) .)
}

setMethod("norm", signature(x = "dsyMatrix", type = "character"),
          function(x, type, ...)
              if(identical(type, "2"))
                  norm2(x)
              else .Call(dsyMatrix_norm, x, type))

setMethod("rcond", signature(x = "dsyMatrix", norm = "character"),
          function(x, norm, ...) .Call(dsyMatrix_rcond, x, norm))

setMethod("solve", signature(a = "dsyMatrix", b = "missing"),
          function(a, b, ...) .Call(dsyMatrix_solve, a))

setMethod("solve", signature(a = "dsyMatrix", b = "Matrix"),
	  function(a, b, ...)
              .Call(dsyMatrix_matrix_solve, a, as(b, "denseMatrix")))

setMethod("solve", signature(a = "dsyMatrix", b = "matrix"),
          function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b))

setMethod("solve", signature(a = "dsyMatrix", b = "numLike"),
	  function(a, b, ...) .Call(dsyMatrix_matrix_solve, a, b))

.is.na <- .is.infinite <- .is.finite <- function(x) {
    n <- (d <- x@Dim)[1L]
    i <- is.na(x@x)
    uplo <- x@uplo
    if(any(i[indTri(n, uplo == "U", TRUE)]))
        new("nsyMatrix", Dim = d, Dimnames = x@Dimnames,
            uplo = uplo, x = i)
    else is.na_nsp(x)
}
body(.is.infinite) <-
    do.call(substitute, list(body(.is.infinite),
                             list(is.na = quote(is.infinite))))
body(.is.finite) <-
    do.call(substitute, list(body(.is.finite),
                             list(is.na = quote(is.finite))))

for(.cl in paste0(c("d", "l"), "syMatrix"))
    setMethod("is.na", signature(x = .cl), .is.na)
setMethod("is.infinite", signature(x = "dsyMatrix"), .is.infinite)
setMethod("is.finite", signature(x = "dsyMatrix"), .is.finite)

rm(.cl, .is.na, .is.infinite, .is.finite)
