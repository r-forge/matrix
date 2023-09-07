## These are used in ./diagMatrix.R, ./indMatrix.R, ./pMatrix.R :

mmultDim <- function(d.a, d.b, type = 1L) {
    ## Return the 'dim' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b   { crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    ## after asserting that ncol(<left operand>) == nrow(<right operand>)
    i.a <- 1L + (type != 2L)
    i.b <- 1L + (type == 3L)
    if(d.a[i.a] != d.b[i.b])
        stop(gettextf("non-conformable arguments in %s",
                      deparse(sys.call(sys.parent()))),
             call. = FALSE, domain = NA)
    c(d.a[-i.a], d.b[-i.b])
}

mmultDimnames <- function(dn.a, dn.b, type = 1L)
    ## Return the 'dimnames' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b   { crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    c(if(is.null(dn.a)) list(NULL) else dn.a[2L - (type != 2L)],
      if(is.null(dn.b)) list(NULL) else dn.b[2L - (type == 3L)])


## METHODS FOR GENERIC: %*%
## NB: x %*% y == t(t(y) %*% t(x))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(.cl in c("Matrix", "sparseVector")) {
setMethod("%*%", signature(x = .cl, y = "ANY"),
          function(x, y)
              x %*% (if(length(dim(y)) == 2L) as.matrix else as.vector)(y))

setMethod("%*%", signature(x = "ANY", y = .cl),
          function(x, y)
              (if(length(dim(x)) == 2L) as.matrix else as.vector)(x) %*% y)
}


setMethod("%*%", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y)
              .Call( R_dense_prod, x, y, FALSE, FALSE))

for(.cl in c("matrix", "vector")) {
setMethod("%*%", signature(x = "denseMatrix", y = .cl),
          function(x, y)
              .Call( R_dense_prod, x, y, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "denseMatrix"),
          function(x, y)
              .Call( R_dense_prod, x, y, FALSE, FALSE))
}


setMethod("%*%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "CsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))
}


setMethod("%*%", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))

setMethod("%*%", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))

setMethod("%*%", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "RsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))
}


setMethod("%*%", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))

setMethod("%*%", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%*%", signature(x = "TsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE, FALSE))

setMethod("%*%", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE, FALSE))
}


setMethod("%*%", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y)
              if((nx <- length(x)) == 1L)
                  .tCRT(.V2C(x * y))
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L)
           ##           .V2C(x * y)
           ## else if(nx == ny)
                  .m2sparse(sum(x * y), ".gR")
              else stop("non-conformable arguments"))

for(.cl in c("Matrix", "matrix")) {
setMethod("%*%", signature(x = "sparseVector", y = .cl),
          function(x, y)
              if((k <- dim(y)[1L]) == 1L)
                        .V2C(x)  %*% y
              else if(length(x) == k)
                  .tCRT(.V2C(x)) %*% y
              else stop("non-conformable arguments"))

setMethod("%*%", signature(x = .cl, y = "sparseVector"),
          function(x, y)
              if((k <- dim(x)[2L]) == 1L)
                  x %*% .tCRT(.V2C(y))
              else if(length(y) == k)
                  x %*%       .V2C(y)
              else stop("non-conformable arguments"))
}

setMethod("%*%", signature(x = "sparseVector", y = "vector"),
          function(x, y)
              if((nx <- length(x)) == 1L)
                  .m2dense(.V2v(x * y), ".ge", trans = TRUE)
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L)
           ##     .m2dense(.V2v(x * y), ".ge")
           ## else if(nx == ny)
                  .m2dense(sum(x * y), ".ge")
              else stop("non-conformable arguments"))

setMethod("%*%", signature(x = "vector", y = "sparseVector"),
          function(x, y)
              if((nx <- length(x)) == 1L)
                  .m2dense(.V2v(x * y), ".ge", trans = TRUE)
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L)
           ##     .m2dense(.V2v(x * y), ".ge")
           ## else if(nx == ny)
                  .m2dense(sum(x * y), ".ge")
              else stop("non-conformable arguments"))


## METHODS FOR GENERIC: %&%
## NB: x %*% y == t(t(y) %*% t(x))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("%&%", signature(x = "ANY", y = "ANY"),
          function(x, y)
              (if(length(dim(x)) == 2L) as.matrix else as.vector)(x) %&%
              (if(length(dim(y)) == 2L) as.matrix else as.vector)(y))

for(.cl in c("Matrix", "sparseVector", "matrix", "vector")) {
setMethod("%&%", signature(x = .cl, y = "ANY"),
          function(x, y)
              x %&% (if(length(dim(y)) == 2L) as.matrix else as.vector)(y))

setMethod("%&%", signature(x = "ANY", y = .cl),
          function(x, y)
              (if(length(dim(x)) == 2L) as.matrix else as.vector)(x) %&% y)
}


setMethod("%&%", signature(x = "matrix", y = "matrix"),
          function(x, y)
              .m2sparse(x, "ngC") %&% .m2sparse(y, "ngC"))

setMethod("%&%", signature(x = "matrix", y = "vector"),
          function(x, y)
              .m2sparse(x, "ngC") %&%           y        )

setMethod("%&%", signature(x = "vector", y = "matrix"),
          function(x, y)
                        x         %&% .m2sparse(y, "ngC"))

setMethod("%&%", signature(x = "vector", y = "vector"),
          function(x, y) {
              r <-
              if((nx <- length(x)) == 1L)
                  .m2sparse(x, "ngC") %&% .m2sparse(y, "ngR", trans = TRUE)
              else if(nx == length(y))
           ## else if((ny <- length(y)) == 1L || nx == ny)
                  .m2sparse(x, "ngR", trans = TRUE) %&% .m2sparse(y, "ngC")
              else stop("non-conformable arguments")
              r@Dimnames <- list(NULL, NULL)
              r
          })


setMethod("%&%", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

for(.cl in c("matrix", "vector")) {
setMethod("%&%", signature(x = "denseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "denseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))
}


setMethod("%&%", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "CsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))
}


setMethod("%&%", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "RsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE,  TRUE))
}


setMethod("%&%", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, y, x,  TRUE,  TRUE,  TRUE,  TRUE))

setMethod("%&%", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("%&%", signature(x = "TsparseMatrix", y = .cl),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))

setMethod("%&%", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y)
              .Call(R_sparse_prod, x, y, FALSE, FALSE, FALSE,  TRUE))
}


setMethod("%&%", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y) {
              x <- .V2kind(.drop0(x, isM = FALSE), "n")
              y <- .V2kind(.drop0(y, isM = FALSE), "n")
              if((nx <- length(x)) == 1L) {
                  if(length(x@i) == 0L)
                      y@i <- integer(0L)
                  .tCRT(.V2C(y))
              } else if(nx == length(y))
           ## } else if((ny <- length(y)) == 1L) {
           ##     if(length(y@i))
           ##         x@i <- integer(0L)
           ##           .V2C(x)
           ## } else if(nx == ny)
                  .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
              else stop("non-conformable arguments")
          })

for(.cl in c("Matrix", "matrix")) {
setMethod("%&%", signature(x = "sparseVector", y = .cl),
          function(x, y) {
              x <- .V2kind(.drop0(x, isM = FALSE), "n")
              if((k <- dim(y)[1L]) == 1L)
                        .V2C(x)  %&% y
              else if(length(x) == k)
                  .tCRT(.V2C(x)) %&% y
              else stop("non-conformable arguments")
          })

setMethod("%&%", signature(x = .cl, y = "sparseVector"),
          function(x, y) {
              y <- .V2kind(.drop0(y, isM = FALSE), "n")
              if((k <- dim(x)[2L]) == 1L)
                  x %&% .tCRT(.V2C(y))
              else if(length(y) == k)
                  x %&%       .V2C(y)
              else stop("non-conformable arguments")
          })
}

setMethod("%&%", signature(x = "sparseVector", y = "vector"),
          function(x, y)
              .V2kind(.drop0(x, isM = FALSE), "n") %&% .m2V(y, "n"))

setMethod("%&%", signature(x = "vector", y = "sparseVector"),
          function(x, y)
              .m2V(x, "n") %&% .V2kind(.drop0(y, isM = FALSE), "n"))


## METHODS FOR GENERIC: crossprod
## NB: t(x) %*% y == t(t(y) %*% x)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(.cl in c("Matrix", "sparseVector")) {
setMethod("crossprod", signature(x = .cl, y = "ANY"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(t(x), y))

setMethod("crossprod", signature(x = "ANY", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(t(x), y))
}


setMethod("crossprod", signature(x = "denseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" else boolArith)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y,  TRUE, FALSE))

setMethod("crossprod", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" && .M.kind(y) == "n" else boolArith)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y,  TRUE, FALSE))

for(.cl in c("matrix", "vector")) {
setMethod("crossprod", signature(x = "denseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y,  TRUE, FALSE))

setMethod("crossprod", signature(x = .cl, y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y,  TRUE, FALSE))
}


setMethod("crossprod", signature(x = "CsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "CsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  TRUE, FALSE,  TRUE, boolArith))
}


setMethod("crossprod", signature(x = "RsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE,  TRUE, boolArith))

setMethod("crossprod", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  TRUE, FALSE,  TRUE, boolArith))

setMethod("crossprod", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  TRUE, FALSE,  TRUE, boolArith))

setMethod("crossprod", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  TRUE, FALSE,  TRUE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "RsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  TRUE, FALSE,  TRUE, boolArith))
}


setMethod("crossprod", signature(x = "TsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("crossprod", signature(x = "TsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  TRUE, FALSE, FALSE, boolArith))

setMethod("crossprod", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  TRUE, FALSE,  TRUE, boolArith))
}


setMethod("crossprod", signature(x = "sparseVector", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .V.kind(x) == "n" else boolArith) {
                  if(!is.na(boolArith))
                      x <- .V2kind(.drop0(x, isM = FALSE), "n")
                  .m2sparse(length(x@i) > 0L, "nsR")
              } else .m2sparse(sum(x * x), ".sR"))

setMethod("crossprod", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .V.kind(x) == "n" && .V.kind(y) == "n" else boolArith) {
                  if(!is.na(boolArith)) {
                      x <- .V2kind(.drop0(x, isM = FALSE), "n")
                      y <- .V2kind(.drop0(y, isM = FALSE), "n")
                  }
                  if((nx <- length(x)) == 1L) {
                      if(length(x@i) == 0L)
                          y@i <- integer(0L)
                      .tCRT(.V2C(y))
                  } else if(nx == length(y))
                      .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
                  else stop("non-conformable arguments")
              } else {
                  if((nx <- length(x)) == 1L)
                      .tCRT(.V2C(x * y))
                  else if(nx == length(y))
                      .m2sparse(sum(x * y), ".gR")
                  else stop("non-conformable arguments")
              })

for(.cl in c("Matrix", "matrix")) {
setMethod("crossprod", signature(x = "sparseVector", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              crossprod(.tCRT(.V2C(x)), y, boolArith = boolArith, ...))

setMethod("crossprod", signature(x = .cl, y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(dim(x)[1L] == 1L)
              crossprod(x, .tCRT(.V2C(y)), boolArith = boolArith, ...)
              else
              crossprod(x,       .V2C(y) , boolArith = boolArith, ...))
}

setMethod("crossprod", signature(x = "sparseVector", y = "vector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith) {
                  x <- .V2kind(.drop0(x, isM = FALSE), "n")
                  y <- .m2V(y, "n")
                  if((nx <- length(x)) == 1L) {
                      if(length(x@i) == 0L)
                          y@i <- integer(0L)
                      .tCRT(.V2C(y))
                  } else if(nx == length(y))
                      .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
                  else stop("non-conformable arguments")
              }
              else if((nx <- length(x)) == 1L)
                  .m2dense(.V2v(x * y), ".ge", trans = TRUE)
              else if(nx == length(y))
                  .m2dense(sum(x * y), ".ge")
              else stop("non-conformable arguments"))

setMethod("crossprod", signature(x = "vector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith) {
                  x <- .m2V(x, "n")
                  y <- .V2kind(.drop0(y, isM = FALSE), "n")
                  if((nx <- length(x)) == 1L) {
                      if(length(x@i) == 0L)
                          y@i <- integer(0L)
                      .tCRT(.V2C(y))
                  } else if(nx == length(y))
                      .m2sparse(any(match(x@i, y@i, 0L)), "ngR")
                  else stop("non-conformable arguments")
              }
              else if((nx <- length(x)) == 1L)
                  .m2dense(.V2v(x * y), ".ge", trans = TRUE)
              else if(nx == length(y))
                  .m2dense(sum(x * y), ".ge")
              else stop("non-conformable arguments"))


## METHODS FOR GENERIC: tcrossprod
## NB: x %*% t(y) == t(y %*% t(x))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(.cl in c("Matrix", "sparseVector")) {
setMethod("tcrossprod", signature(x = .cl, y = "ANY"),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(x, t(y)))

setMethod("tcrossprod", signature(x = "ANY", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              (if(!is.na(boolArith) && boolArith) `%&%` else `%*%`)(x, t(y)))
}


setMethod("tcrossprod", signature(x = "denseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" else boolArith)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y, FALSE,  TRUE))

setMethod("tcrossprod", signature(x = "denseMatrix", y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .M.kind(x) == "n" && .M.kind(y) == "n" else boolArith)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y, FALSE,  TRUE))

for(.cl in c("matrix", "vector")) {
setMethod("tcrossprod", signature(x = "denseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y, FALSE,  TRUE))

setMethod("tcrossprod", signature(x = .cl, y = "denseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE,  TRUE)
              else
              .Call( R_dense_prod, x, y, FALSE,  TRUE))
}


setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "CsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "CsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x, FALSE,  TRUE,  TRUE, boolArith))
}


setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y,  FALSE,  TRUE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  FALSE,  TRUE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  FALSE,  TRUE,  TRUE, boolArith))

setMethod("tcrossprod", signature(x = "RsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x,  FALSE,  TRUE,  TRUE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "RsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x, FALSE,  TRUE,  TRUE, boolArith))
}


setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "CsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "RsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = "TsparseMatrix", y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

for(.cl in c("denseMatrix", "matrix", "vector")) {
setMethod("tcrossprod", signature(x = "TsparseMatrix", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, x, y, FALSE,  TRUE, FALSE, boolArith))

setMethod("tcrossprod", signature(x = .cl, y = "TsparseMatrix"),
          function(x, y = NULL, boolArith = NA, ...)
              .Call(R_sparse_prod, y, x, FALSE,  TRUE,  TRUE, boolArith))
}


setMethod("tcrossprod", signature(x = "sparseVector", y = "missing"),
          function(x, y = NULL, boolArith = NA, ...)
              tcrossprod(.V2C(x), boolArith = boolArith, ...))

setMethod("tcrossprod", signature(x = "sparseVector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(if(is.na(boolArith)) .V.kind(x) == "n" && .V.kind(y) == "n" else boolArith)
                  x %&% .tCRT(.V2C(y))
              else
                  x %*% .tCRT(.V2C(y)))

for(.cl in c("Matrix", "matrix")) {
setMethod("tcrossprod", signature(x = "sparseVector", y = .cl),
          function(x, y = NULL, boolArith = NA, ...)
              if(dim(y)[2L] == 1L)
              tcrossprod(      .V2C(x) , y, boolArith = boolArith, ...)
              else
              tcrossprod(.tCRT(.V2C(x)), y, boolArith = boolArith, ...))

setMethod("tcrossprod", signature(x = .cl, y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              tcrossprod(x, .tCRT(.V2C(y)), boolArith = boolArith, ...))
}

setMethod("tcrossprod", signature(x = "sparseVector", y = "vector"),
          function(x, y = NULL, boolArith = NA, ...) {
              r <-
              if(!is.na(boolArith) && boolArith)
                  x %&% .m2sparse(y, "ngR", trans = TRUE)
              else
                  x %*% .m2dense (y, ".ge", trans = TRUE)
              r@Dimnames[[2L]] <- NULL
              r
          })

setMethod("tcrossprod", signature(x = "vector", y = "sparseVector"),
          function(x, y = NULL, boolArith = NA, ...)
              if(!is.na(boolArith) && boolArith)
                  x %&% .tCRT(.V2C(y))
              else
                  x %*% .tCRT(.V2C(y)))
