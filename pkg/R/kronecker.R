#### Collect methods for  kronecker() here.
####			  ===========

### ... all but the ``fall back methods'' which are in ./Matrix.R ...

### Request: Should be *fast* particularly when used with Diagonal() !

setMethod("kronecker", signature(X="diagonalMatrix", Y="ANY"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
	      X <- as(X, "TsparseMatrix")
	      callGeneric()
	  })
setMethod("kronecker", signature(X="ANY", Y="diagonalMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
	      Y <- as(Y, "TsparseMatrix")
	      callGeneric()
	  })
## also needs
setMethod("kronecker", signature(X="sparseMatrix", Y="ANY"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
	      Y <- as(Y, "TsparseMatrix")
	      callGeneric()
	  })
setMethod("kronecker", signature(X="ANY", Y="sparseMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
	      X <- as(X, "TsparseMatrix")
	      callGeneric()
	  })

## from ./dgTMatrix.R :
setMethod("kronecker", signature(X = "dgTMatrix", Y = "dgTMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
      {
	  if (FUN != "*") stop("kronecker method must use default 'FUN'")
	  ## otherwise we don't know that many results will be zero
	  ydim <- Y@Dim
	  xi <- X@i
	  xnnz <- length(xi)
	  yi <- Y@i
	  ynnz <- length(yi)
	  new("dgTMatrix", Dim = X@Dim * ydim,
	      i = rep.int(yi, xnnz) + ydim[1] * rep.int(xi, rep.int(ynnz, xnnz)),
	      j = rep.int(Y@j, xnnz) + ydim[2] * rep.int(X@j, rep.int(ynnz, xnnz)),
	      ## x = as.vector(outer(Y@x, X@x, FUN = FUN)
	      x = as.vector(Y@x %*% t(X@x)))
      }, valueClass = "dgTMatrix")

## from ./dsparseMatrix.R :
setMethod("kronecker", signature(X = "dsparseMatrix", Y = "dsparseMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
	  callGeneric(as(X, "dgTMatrix"), as(Y, "dgTMatrix")))
