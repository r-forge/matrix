## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a dgeMatrix.

setMethod("%*%", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
          function(x, y) callGeneric(x, as(y, "dgeMatrix")))

setMethod("%*%", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
          function(x, y) callGeneric(as(x, "dgeMatrix"), y))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "ddenseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "dsparseMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix"), y))

## and coerce dsparse* to dgC*
setMethod("%*%", signature(x = "dsparseMatrix", y = "dgeMatrix"),
          function(x, y) callGeneric(as(x, "dgCMatrix"), y))

setMethod("%*%", signature(x = "dgeMatrix", y = "dsparseMatrix"),
          function(x, y) callGeneric(x, as(y, "dgCMatrix")))

setMethod("crossprod", signature(x = "dsparseMatrix", y = "dgeMatrix"),
          function(x, y = NULL) callGeneric(as(x, "dgCMatrix"), y))

## NB: there's already
##     ("CsparseMatrix", "missing") and ("TsparseMatrix", "missing") methods

setMethod("crossprod", signature(x = "dgeMatrix", y = "dsparseMatrix"),
          function(x, y = NULL) callGeneric(x, as(y, "dgCMatrix")))

setMethod("image", "dsparseMatrix",
          function(x, ...) image(as(x, "dgTMatrix"), ...))

setMethod("kronecker", signature(X = "dsparseMatrix", Y = "dsparseMatrix"),
          function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
          callGeneric(as(X, "dgTMatrix"),as(Y, "dgTMatrix")))


## Group Methods, see ?Arith (e.g.)
## -----

## NOT YET (first need 'dgCMatrix' method):
## setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
##           signature(e1 = "dsparseMatrix", e2 = "dsparseMatrix"),
##           function(e1, e2) callGeneric(as(e1, "dgCMatrix"),
##                                        as(e2, "dgCMatrix")))
## setMethod("Arith",
##           signature(e1 = "dsparseMatrix", e2 = "numeric"),
##           function(e1, e2) callGeneric(as(e1, "dgCMatrix"), e2))
## setMethod("Arith",
##           signature(e1 = "numeric", e2 = "dsparseMatrix"),
##           function(e1, e2) callGeneric(e1, as(e2, "dgCMatrix")))

setMethod("Math",
	  signature(x = "dsparseMatrix"),
	  function(x) {
	      r <- callGeneric(as(x, "dgCMatrix"))
	      if(is(r, "dsparseMatrix")) as(r, class(x))
	  })

if(FALSE) ## unneeded with "Math2" in ./dMatrix.R
setMethod("Math2",
	  signature(x = "dsparseMatrix", digits = "numeric"),
	  function(x, digits) {
	      r <- callGeneric(as(x, "dgCMatrix"), digits = digits)
	      if(is(r, "dsparseMatrix")) as(r, class(x))
	  })


### cbind2 / rbind2
if(paste(R.version$major, R.version$minor, sep=".") >= "2.2") {
    ## for R 2.2.x (and later):

### cbind2
    setMethod("cbind2", signature(x = "dsparseMatrix", y = "numeric"),
	      function(x, y) {
		  d <- dim(x); nr <- d[1]; nc <- d[2]; cl <- class(x)
                  x <- as(x, "dgCMatrix")
                  if(nr > 0) {
                      y <- rep(y, length.out = nr) # 'silent procrustes'
                      n0y <- y != 0
                      n.e <- length(x@i)
                      x@i <- c(x@i, (0:(nr-1))[n0y])
                      x@p <- c(x@p, n.e + sum(n0y))
                      x@x <- c(x@x, y[n0y])
                  } else { ## nr == 0

                  }
		  x@Dim[2] <- nc + 1:1
		  if(is.character(dn <- x@Dimnames[[2]]))
		      x@Dimnames[[2]] <- c(dn, "")
		  x
	      })
    ## the same, (x,y) <-> (y,x):
    setMethod("cbind2", signature(x = "numeric", y = "dsparseMatrix"),
	      function(x, y) {
		  d <- dim(y); nr <- d[1]; nc <- d[2]; cl <- class(y)
                  y <- as(y, "dgCMatrix")
                  if(nr > 0) {
                      x <- rep(x, length.out = nr) # 'silent procrustes'
                      n0x <- x != 0
                      y@i <- c((0:(nr-1))[n0x], y@i)
                      y@p <- c(0:0, sum(n0x) + y@p)
                      y@x <- c(x[n0x], y@x)
                  } else { ## nr == 0

                  }
		  y@Dim[2] <- nc + 1:1
		  if(is.character(dn <- y@Dimnames[[2]]))
		      y@Dimnames[[2]] <- c(dn, "")
		  y
	      })


    setMethod("cbind2", signature(x = "dsparseMatrix", y = "matrix"),
	      function(x, y) callGeneric(x, as(y, "dgCMatrix")))
    setMethod("cbind2", signature(x = "matrix", y = "dsparseMatrix"),
	      function(x, y) callGeneric(as(x, "dgCMatrix"), y))

    setMethod("cbind2", signature(x = "dsparseMatrix", y = "dsparseMatrix"),
	      function(x, y) {
		  nr <- rowCheck(x,y)
		  ncx <- x@Dim[2]
		  ncy <- y@Dim[2]
		  ## beware of (packed) triangular, symmetric, ...
		  hasDN <- !is.null(dnx <- dimnames(x)) |
			   !is.null(dny <- dimnames(y))
                  x <- as(x, "dgCMatrix")
                  y <- as(y, "dgCMatrix")
                  ne.x <- length(x@i)
                  x@i <- c(x@i, y@i)
                  x@p <- c(x@p, ne.x + y@p[-1])
		  x@x <- c(x@x, y@x)
		  x@Dim[2] <- ncx + ncy
		  if(hasDN) {
		      ## R and S+ are different in which names they take
		      ## if they differ -- but there's no warning in any case
		      rn <- if(!is.null(dnx[[1]])) dnx[[1]] else dny[[1]]
		      cx <- dnx[[2]] ; cy <- dny[[2]]
		      cn <- if(is.null(cx) && is.null(cy)) NULL
		      else c(if(!is.null(cx)) cx else rep.int("", ncx),
			     if(!is.null(cy)) cy else rep.int("", ncy))
		      x@Dimnames <- list(rn, cn)
		  }
		  x
	      })

### rbind2 -- analogous to cbind2 --- more to do for @x though:


}## R-2.2.x ff
