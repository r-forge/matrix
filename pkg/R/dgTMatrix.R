setAs("dgTMatrix", "dgCMatrix",
      function(from) .Call("dgTMatrix_to_dgCMatrix", from, PACKAGE = "Matrix")
      )

setAs("dgTMatrix", "dgeMatrix",
      function(from) .Call("dgTMatrix_to_dgeMatrix", from, PACKAGE = "Matrix")
      )

setAs("dgTMatrix", "matrix",
      function(from) .Call("dgTMatrix_to_matrix", from, PACKAGE = "Matrix")
      )

setAs("dgeMatrix", "dgTMatrix",
      function(from) as(as(from, "dgCMatrix"), "dgTMatrix"))

setAs("dgTMatrix", "dsCMatrix",
      function(from) {
          if (!isSymmetric(from))
              stop("cannot coerce non-symmetric dgTMatrix to dsCMatrix class")
          upper <- from@i <= from@j
          uC <- as(new("dgTMatrix", Dim = from@Dim, i = from@i[upper],
                       j = from@j[upper], x = from@x[upper]), "dgCMatrix")
          new("dsCMatrix", Dim = uC@Dim, p = uC@p, i = uC@i, x = uC@x, uplo = "U")
      })

setAs("dgTMatrix", "dtTMatrix",
      function(from) gt2tT(from, uplo = from@uplo, diag = from@diag))

setAs("matrix", "dgTMatrix",
      function(from) {
	  x <- as.double(from)
	  nz <- as.logical(x)
	  new("dgTMatrix", Dim = dim(from),
	      i = row(from)[nz] - 1:1,
	      j = col(from)[nz] - 1:1,
	      x = x[nz])
      })



## "[" methods are now in ./Tsparse.R

## FIXME? -- should these be moved to /Tsparse.R -- for *all* Tsparse ones?

## setReplaceMethod("[", signature(x = "dgTMatrix", i = "index", j = "missing",
##                                 value = "numeric"),
##                  function (x, i, value) {

##                      stop("NOT YET")

##                      as(r, class(x))
##                  })

## setReplaceMethod("[", signature(x = "dgTMatrix", i = "missing", j = "index",
##                                 value = "numeric"),
##                  function (x, j, value) {

##                      stop("NOT YET")

##                      as(r, class(x))
##                  })

setReplaceMethod("[", signature(x = "dgTMatrix", i = "index", j = "index",
                                value = "numeric"),
                 function (x, i, j, value)
             {
                 di <- dim(x)
                 dn <- dimnames(x)
                 ## FIXME: .ind.prep() is not sufficient;
                 ## ----- since we need to "know" which (i,j) have a match
                 ##       and which have *not*
                 ip1 <- .ind.prep(x@i, i, 1, di, dn)
                 ip2 <- .ind.prep(x@j, j, 2, di, dn)
                 sel <- ip1$m > 0:0  &  ip2$m > 0:0
                 ## the simplest case
                 if(all(value == 0)) { ## just drop the non-zero entries
                     if(any(sel)) { ## non-zero there
                         x@i <- x@i[!sel]
                         x@j <- x@j[!sel]
                         x@x <- x@x[!sel]
                     }
                     return(x)
                 }

                 ## else --  some( value != 0 ) --

                 ## another simple, typical case:
                 if(is.numeric(i) && length(i) == 1 &&
                    is.numeric(j) && length(j) == 1) {
                     if(length(value) != 1) {
                         if(length(value == 0))
                             stop("nothing to replace")
                         else
                             stop("too many replacement values")
                     }
                     if(any(sel)) { ## non-zero there
                         x@x[sel] <- value
                     } else {
                         x@i <- c(x@i, as.integer(i - 1))
                         x@j <- c(x@j, as.integer(j - 1))
                         x@x <- c(x@x, value)
                     }
                     return(x)
                 }

                 ## 1) replace those that are already non-zero

                 ## 2) add     those that were structural 0

                 stop("not-yet-implemented 'dgTMatrix[<-' method")

                 if(any(sel)) { ## non-zero there
                     x@x[sel] <- value
                 }

             })




setMethod("crossprod", signature(x = "dgTMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_crossprod", as(x, "dgCMatrix"), PACKAGE = "Matrix"))

setMethod("crossprod", signature(x = "dgTMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", as(x, "dgCMatrix"), y, PACKAGE = "Matrix"))

##setMethod("crossprod", signature(x = "dgTMatrix", y = "numeric"),
##          function(x, y = NULL)
##          .Call("csc_matrix_crossprod", as(x, "dgCMatrix"), as.matrix(y), PACKAGE = "Matrix"))

setMethod("tcrossprod", signature(x = "dgTMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_tcrossprod", as(x, "dgCMatrix"), PACKAGE = "Matrix"))

setMethod("image", "dgTMatrix",
          function(x,
                   xlim = c(-0.5, matdim[2]-0.5),
                   ylim = c(matdim[1]-0.5, -0.5),
                   sub = sprintf("Dimensions: %d x %d", matdim[1], matdim[2]),
                   xlab = "Column", ylab = "Row",
                   cuts = 20,
                   col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
                   ...)
      {
          matdim <- x@Dim
          levelplot(abs(x@x) ~ x@j * x@i,
                    sub = sub,
                    xlab = xlab, ylab = ylab,
                    xlim = xlim, ylim = ylim,
                    col.regions = col.regions,
                    par.settings = list(background = list(col = "transparent")),
                    panel = function(x, y, z, subscripts, at, ..., col.regions)
                {
                    x <- as.numeric(x[subscripts])
                    y <- as.numeric(y[subscripts])

                    numcol <- length(at) - 1
                    numcol.r <- length(col.regions)
                    col.regions <-
                        if (numcol.r <= numcol)
                            rep(col.regions, length = numcol)
                        else col.regions[floor(1+(1:numcol-1)*(numcol.r-1)/
                                               (numcol-1))]
                    zcol <- rep(NA, length(z)) #numeric(length(z))
                    for (i in seq(along = col.regions))
                        zcol[!is.na(x) & !is.na(y) & !is.na(z) &
                             z>=at[i] & z<at[i+1]] <- i

                    zcol <- as.numeric(zcol[subscripts])
                    if (any(subscripts))
                        grid.rect(x = x, y = y, width = 1, height = 1,
                                  default.units = "native",
                                  gp = gpar(fill = col.regions[zcol],
                                  col = NULL))
                }, ...)
      })

## Uses the triplet convention of *adding* entries with same (i,j):
setMethod("+", signature(e1 = "dgTMatrix", e2 = "dgTMatrix"),
          function(e1, e2) {
              dimCheck(e1, e2)
              new("dgTMatrix", i = c(e1@i, e2@i), j = c(e1@j, e2@j),
                  x = c(e1@x, e2@x), Dim = e1@Dim)
          })

setMethod("t", signature(x = "dgTMatrix"),
          function(x)
          new("dgTMatrix", i = x@j, j = x@i, x = x@x, Dim = rev(x@Dim)))

setMethod("kronecker", signature(X = "dgTMatrix", Y = "dgTMatrix"),
          function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
      {
          if (FUN != "*") stop("kronecker method must use default 'FUN'")
          ydim <- Y@Dim
          xi <- X@i
          xnnz <- length(xi)
          yi <- Y@i
          ynnz <- length(yi)
          new("dgTMatrix", Dim = X@Dim * ydim,
              i = rep.int(yi, xnnz) + ydim[1] * rep.int(xi, rep.int(ynnz, xnnz)),
              j = rep.int(Y@j, xnnz) + ydim[2] * rep.int(X@j, rep.int(ynnz, xnnz)),
              x = as.vector(outer(Y@x, X@x)))
      }, valueClass = "dgTMatrix")

setMethod("writeHB", signature(obj = "dgTMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeHarwellBoeing", obj, as.character(file), "DGT", PACKAGE = "Matrix"))

setMethod("writeMM", signature(obj = "dgTMatrix"),
          function(obj, file, ...)
          .Call("Matrix_writeMatrixMarket", obj, as.character(file), "DGT", PACKAGE = "Matrix"))
