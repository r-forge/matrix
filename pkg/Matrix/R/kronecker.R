## METHODS FOR GENERIC: kronecker
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Fallback methods, typically going via base::kronecker
setMethod("kronecker", signature(X = "Matrix", Y = "ANY"),
	  function(X, Y, FUN, make.dimnames, ...) {
	      if(is(X, "sparseMatrix"))
		  warning("kronecker() coercing from sparse to dense")
	      Matrix(kronecker(as(X, "matrix"), Y, FUN, make.dimnames, ...))
          })

setMethod("kronecker", signature(X = "ANY", Y = "Matrix"),
	  function(X, Y, FUN, make.dimnames, ...) {
	      if(is(Y, "sparseMatrix"))
		  warning("kronecker() coercing from sparse to dense")
	      Matrix(kronecker(X, as(Y, "matrix"), FUN, make.dimnames, ...))
          })

## FIXME: use this everywhere below
.kroneckerDimnames <- function(dnx, dx = lengths(dnx, FALSE),
                               dny, dy = lengths(dny, FALSE),
                               sep = ":") {
    dnr <- list(NULL, NULL)
    if(identical(dnx, dnr) && identical(dny, dnr))
        return(NULL)
    for(i in 1:2) {
        if(have.sx <- (ny <- dy[i]) && !is.null(sx <- dnx[[i]]))
            sx <- rep    (sx,  each = ny)
        if(have.sy <- (nx <- dx[i]) && !is.null(sy <- dny[[i]]))
            sy <- rep.int(sy, times = nx)
        dnr[[i]] <-
            if(have.sx && have.sy)
                paste0(sx, sep, sy)
            else if(have.sx)
                paste0(sx, sep    )
            else if(have.sy)
                paste0(    sep, sy)
            else if(nx && ny)
                rep.int(sep, nx * ny)
    }
    dnr
}

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "diagonalMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              r <- new("ddiMatrix")
              r@Dim <- dr <- dX * dY
              if((uX <- X@diag != "N") & (uY <- Y@diag != "N"))
                  r@diag <- "U"
              else if(uX)
                  r@x <- rep.int(as.double(Y@x), dX[1L])
              else if(uY)
                  r@x <- rep(as.double(X@x), each = dY[1L])
              else r@x <- rep(X@x, each = dY[1L]) * Y@x
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(X@Dimnames, dX,
                                                    Y@Dimnames, dY)))
                  r@Dimnames <- dnr
              r
	  })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "matrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, .m2ge(Y), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "matrix", Y = "diagonalMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(.m2ge(X), Y, FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "denseMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              uX <- x@diag != "N"
              uY <- FALSE
              shape <- .M.shape(Y)
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- uplo <- Y@uplo
                  if(shape == "t")
                      uY <- Y@diag != "N"
              }
              if(!all(dr)) {
                  r@p <- integer(dr[2L] + 1)
                  if(uX && uY)
                      r@diag <- "U"
              } else {
                  m <- dY[1L]
                  nx <- dX[1L]
                  ny <- length(y <- Y@x)
                  if(shape != "g" && ny > 1L && ny == prod(dy)) {
                      Y <- pack(Y)
                      ny <- length(y <- Y@x)
                  }
                  if(as.double(nx) * ny > .Machine$integer.max)
                      stop("number of nonzero entries would exceed 2^31-1")
                  if(!uX && uY) {
                      diag(Y) <- TRUE
                      ny <- length(y <- Y@x)
                  }
                  if(is.logical(y) && .M.kind(Y) == "n")
                      y <- y | is.na(y)
                  if(shape == "g") {
                      r@p <- seq.int(0L, by = m, length.out = dr[2L] + 1)
                      r@i <-
                          rep(seq.int(0L, by = m, length.out = nx),
                              each = ny) +
                          0:(m-1L)
                  } else if(uplo == "U") {
                      r@p <- c(0L, cumsum(rep.int(1:m, nx)))
                      r@i <-
                          rep(seq.int(0L, by = m, length.out = nx),
                              each = ny) +
                          sequence.default(nvec = 1:m, from = 0L)
                  } else {
                      r@p <- c(0L, cumsum(rep.int(m:1, nx)))
                      r@i <-
                          rep(seq.int(0L, by = m, length.out = nx),
                              each = ny) +
                          sequence.default(nvec = m:1, from = 0:(m-1L))
                  }
                  r@x <-
                      if(uX)
                          rep.int(as.double(y), nx)
                      else as.double(y) * rep(X@x, each = ny)
                  if(uX && uY)
                      r <- ..diagN2U(r)
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(X@Dimnames, dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
	  })

setMethod("kronecker", signature(X = "denseMatrix", Y = "diagonalMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              shape <- .M.shape(X)
              uX <- FALSE
              uY <- Y@diag != "N"
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- uplo <- X@uplo
                  if(shape == "t")
                      uX <- X@diag != "N"
              }
              if(!all(dr)) {
                  r@p <- integer(dr[2L] + 1)
                  if(uX && uY)
                      r@diag <- "U"
              } else {
                  m <- dX[1L]
                  n <- dX[2L]
                  nx <- length(x <- X@x)
                  ny <- dY[1L]
                  if(shape != "g" && nx > 1L && nx == prod(dx)) {
                      X <- pack(X)
                      nx <- length(x <- X@x)
                  }
                  if(as.double(nx) * ny > .Machine$integer.max)
                      stop("number of nonzero entries would exceed 2^31-1")
                  if(uX && !uY) {
                      diag(X) <- TRUE
                      nx <- length(x <- X@x)
                  }
                  if(is.logical(x) && .M.kind(X) == "n")
                      x <- x | is.na(x)
                  if(shape == "g") {
                      x. <- function() {
                          j <- rep(1:n, each = ny)
                          as.double(`dim<-`(as.double(x), dX)[, j])
                      }
                      r@p <- seq.int(0L, by = m, length.out = dr[2L] + 1)
                      r@i <- rep.int(sequence.default(nvec = rep.int(m, ny),
                                                      from = 0:(k-1L),
                                                      by = ny),
                                     n)
                      r@x <-
                          if(uY)
                              x.()
                          else x.() * rep(Y@x, each = m)
                  } else if(uplo == "U") {
                      rep.1.n <- rep(1:n, each = ny)
                      x. <- function() {
                          k <- sequence.default(
                              nvec = rep.1.n,
                              from = rep(1L + cumsum(0:(n-1L)), each = ny),
                              by = 1L)
                          as.double(x)[k]
                      }
                      r@p <- c(0L, cumsum(rep.1.n))
                      r@i <- sequence.default(nvec = rep.1.n,
                                              from = rep.int(0:(ny-1L), n),
                                              by = ny)
                      r@x <-
                          if(uY)
                              x.()
                          else x.() * rep.int(rep.int(Y@x, n), rep.1.n)
                  } else {
                      rep.n.1 <- rep(n:1, each = ny)
                      x. <- function() {
                          k <- sequence.default(
                              nvec = rep.n.1,
                              from = rep(1L + cumsum(c(0L, if(n > 1L) n:2)),
                                         each = ny),
                              by = 1L)
                          as.double(x)[k]
                      }
                      r@p <- c(0L, cumsum(rep.n.1))
                      r@i <- sequence.default(nvec = rep.n.1,
                                              from = 0:(n*ny-1L),
                                              by = ny)
                      r@x <-
                          if(uY)
                              x.()
                          else x.() * rep.int(rep.int(Y@x, n), rep.n.1)
                  }
                  if(uX && uY)
                      r <- ..diagN2U(r)
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    Y@Dimnames, dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
	  })

setMethod("kronecker", signature(X = "denseMatrix", Y = "denseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {

          })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "CsparseMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              shape <- .M.shape(Y)
              need.U2N <- FALSE
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- Y@uplo
                  if(shape == "t" && Y@diag != "N") {
                      if(X@diag != "N")
                          r@diag <- "U"
                      else need.U2N <- TRUE
                  }
              }
              if(!all(dr))
                  r@p <- integer(dr[2L] + 1)
              else {
                  if(need.U2N)
                      Y <- ..diagU2N(Y, "C")
                  if((nnz <- (p <- Y@p)[length(p)]) == 0L)
                      r@p <- integer(dr[2L] + 1)
                  else if(as.double(k <- dX[1L]) * nnz > .Machine$integer.max)
                      stop("number of nonzero entries would exceed 2^31-1")
                  else {
                      head. <-
                          if(length(Y@i) > nnz)
                              function(x) x[seq_len(nnz)]
                          else identity
                      r@p <- c(0L, cumsum(rep.int(p[-1L] - p[-length(p)], k)))
                      r@i <-
                          rep(seq.int(0L, by = dY[1L], length.out = k),
                              each = nnz) +
                          head.(Y@i)
                      r@x <-
                          if(X@diag == "N") {
                              if(.M.kind(Y) != "n")
                                  rep(as.double(X@x), each = nnz) * head.(Y@x)
                              else rep(as.double(X@x), each = nnz)
                          } else {
                              if(.M.kind(Y) != "n")
                                  rep.int(as.double(head.(Y@x)), k)
                              else rep.int(1, k * nnz)
                          }
                  }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(X@Dimnames, dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
	  })

setMethod("kronecker", signature(X = "CsparseMatrix", Y = "diagonalMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              shape <- .M.shape(X)
              need.U2N <- FALSE
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- X@uplo
                  if(shape == "t" && X@diag != "N") {
                      if(Y@diag != "N")
                          r@diag <- "U"
                      else need.U2N <- TRUE
                  }
              }
              if(!all(dr))
                  r@p <- integer(dr[2L] + 1)
              else {
                  if(need.U2N)
                      X <- ..diagU2N(X, "C")
                  if((nnz <- (p <- X@p)[length(p)]) == 0L)
                      r@p <- integer(dr[2L] + 1)
                  else if(as.double(k <- dY[1L]) * nnz > .Machine$integer.max)
                      stop("number of nonzero entries would exceed 2^31-1")
                  else {
                      dp <- p[-1L] - p[-length(p)]
                      j. <- which(dp > 0L)
                      F. <- function(f) unlist(lapply(j., f), FALSE, FALSE)
                      r@p <- c(0L, cumsum(rep(dp, each = k)))
                      r@i <- F.({
                          i <- X@i
                          .seq.0.km1 <- 0:(k-1L)
                          function(j) k * i[(p[j]+1L):p[j+1L]] +
                                          rep(.seq.0.km1, each = dp[j])
                      })
                      r@x <-
                          if(Y@diag == "N") {
                              if(.M.kind(X) != "n")
                                  as.double(F.({
                                      x <- X@x
                                      y <- Y@x
                                      function(j) x[(p[j]+1L):p[j+1L]] *
                                                      rep(y, each = dp[j])
                                  }))
                              else
                                  as.double(F.({
                                      y <- Y@x
                                      function(j) rep(y, each = dp[j])
                                  }))
                          } else {
                              if(.M.kind(X) != "n")
                                  as.double(F.({
                                      x <- X@x
                                      function(j) rep.int(x[(p[j]+1L):p[j+1L]], k)
                                  }))
                              else
                                  rep.int(1, k * nnz)
                          }
                  }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    Y@Dimnames, dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
	  })

setMethod("kronecker", signature(X = "CsparseMatrix", Y = "CsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              shape.X <- .M.shape(X)
              shape.Y <- .M.shape(Y)
              shape <- switch(shape.X,
                              g = "g",
                              t = if(shape.Y == "t" && X@uplo == Y@uplo)
                                      "t"
                                  else "g",
                              s = if(shape.Y == "s")
                                      "s"
                                  else "g")
              uX <- uY <- FALSE
              r <- new(`substr<-`("d.CMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- uplo <- X@uplo
                  if(shape == "t" &&
                     (uX <- X@diag != "N") && (uY <- Y@diag != "N"))
                      r@diag <- "U"
              }
              if(!all(dr))
                  r@p <- integer(dr[2L] + 1)
              else {
                  if(uX != uY) {
                      if(uX)
                          X <- ..diagU2N(X, "C")
                      else Y <- ..diagU2N(Y, "C")
                  }
                  if((nnzX <- (pX <- X@p)[length(pX)]) == 0L ||
                     (nnzY <- (pY <- Y@p)[length(pY)]) == 0L)
                      r@p <- integer(dr[2L] + 1)
                  else if(as.double(nnzX) * nnzY > .Machine$integer.max)
                      stop("number of nonzero entries would exceed 2^31-1")
                  else {
                      if(shape == "g") {
                          if(shape.Y != "g")
                              Y <- ..sparse2g(Y)
                      } else if(shape == "s") {
                          if(Y@uplo != uplo)
                              Y <- t(Y)
                      }

                      dpX <- pX[-1L] -         pX[-length(pX)]
                      dpY <- pY[-1L] - (pY. <- pY[-length(pY)])

                      pY. <- pY. + 1L

                      m <- dY[1L]
                      n <- dY[2L]

                      X.i <- X@i
                      X.x <- X@x
                      Y.i <- Y@i
                      Y.x <- Y@x

                      j2ix <- function(j) {
                          k <- (pX[j]+1L):pX[j+1L]
                          r <- rep(dpY, each = dpX[j])
                          s <- sequence.default(nvec = r,
                                                from = rep(pY., each = dpX[j]),
                                                by = 1L)
                          list(rep.int(rep.int(m * X.i[k], n), r) + Y.i[s],
                               rep.int(rep.int(    X.x[k], n), r) * Y.x[s])
                      }
                      L <- unlist(lapply(which(dpX > 0L), j2ix), FALSE, FALSE)

                      r@p <- c(0L, cumsum(rep(dpX, each = n) * dpY))
                      r@i <- unlist(L[c(TRUE, FALSE)], FALSE, FALSE)
                      r@x <- unlist(L[c(FALSE, TRUE)], FALSE, FALSE)
                  }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "RsparseMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .tCR2RC(kronecker(t(X), .tCR2RC(Y), FUN, make.dimnames, ...)))

setMethod("kronecker", signature(X = "RsparseMatrix", Y = "diagonalMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .tCR2RC(kronecker(.tCR2RC(X), t(Y), FUN, make.dimnames, ...)))

setMethod("kronecker", signature(X = "RsparseMatrix", Y = "RsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .tCR2RC(kronecker(.tCR2RC(X), .tCR2RC(Y), FUN, make.dimnames, ...)))

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "TsparseMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              shape <- .M.shape(Y)
              need.U2N <- FALSE
              r <- new(`substr<-`("d.TMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- Y@uplo
                  if(shape == "t" && Y@diag != "N") {
                      if(X@diag != "N")
                          r@diag <- "U"
                      else need.U2N <- TRUE
                  }
              }
              if(all(dr)) {
                  if(any((kind <- .M.kind(Y)) == c("n", "l")))
                      Y <- .Call(Tsparse_aggregate, Y)
                  if(need.U2N)
                      Y <- ..diagU2N(Y, "T")
                  k <- dX[1L]
                  nnz <- length(Y@i)
                  r@i <-
                      rep(seq.int(0L, by = dY[1L], length.out = k),
                          each = nnz) +
                      Y@i
                  r@j <-
                      rep(seq.int(0L, by = dY[2L], length.out = k),
                          each = nnz) +
                      Y@j
                  r@x <-
                      if(X@diag == "N") {
                          if(kind != "n")
                              rep(as.double(X@x), each = nnz) * Y@x
                          else rep(as.double(X@x), each = nnz)
                      } else {
                          if(kind != "n")
                              rep.int(as.double(Y@x), k)
                          else rep.int(1, k * nnz)
                      }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(X@Dimnames, dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
	  })

setMethod("kronecker", signature(X = "TsparseMatrix", Y = "diagonalMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              shape <- .M.shape(X)
              need.U2N <- FALSE
              r <- new(`substr<-`("d.TMatrix", 2L, 2L, shape))
              r@Dim <- dr <- dX * dY
              if(shape != "g") {
                  r@uplo <- X@uplo
                  if(shape == "t" && X@diag != "N") {
                      if(Y@diag != "N")
                          r@diag <- "U"
                      else need.U2N <- TRUE
                  }
              }
              if(all(dr)) {
                  if(any((kind <- .M.kind(X)) == c("n", "l")))
                      X <- .Call(Tsparse_aggregate, X)
                  if(need.U2N)
                      X <- ..diagU2N(X, "T")
                  k <- dY[1L]
                  nnz <- length(X@i)
                  r@i <- rep(k * X@i, each = k) + 0:(k-1L)
                  r@j <- rep(k * X@j, each = k) + 0:(k-1L)
                  r@x <-
                      if(Y@diag == "N") {
                          if(kind != "n")
                              rep(as.double(X@x), each = k) * Y@x
                          else rep(as.double(X@x), each = k)
                      } else {
                          if(kind != "n")
                              rep.int(as.double(Y@x), k)
                          else rep.int(1, k * nnz)
                      }
              }
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    Y@Dimnames, dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
	  })

setMethod("kronecker", signature(X = "TsparseMatrix", Y = "TsparseMatrix"),
          function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")

              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(dimnames(X), dX,
                                                    dimnames(Y), dY))) {
                  if(shape == "s" && !isSymmetricDN(dnr))
                      r <- .sparse2g(r)
                  r@Dimnames <- dnr
              }
              r
          })

setMethod("kronecker", signature(X = "diagonalMatrix", Y = "indMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(X, as(Y, "nsparseMatrix"), FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "indMatrix", Y = "diagonalMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              kronecker(as(X, "nsparseMatrix"), Y, FUN, make.dimnames, ...))

setMethod("kronecker", signature(X = "indMatrix", Y = "indMatrix"),
	  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(!(missing(FUN) || identical(FUN, "*")))
                  stop("method for kronecker() must use default FUN=\"*\"")
              if(any(as.double(dX <- X@Dim) * (dY <- Y@Dim) >
                     .Machine$integer.max))
		  stop("dimensions cannot exceed 2^31-1")
              r <- new("indMatrix")
              r@Dim <- dX * dY
              r@perm <- dY[2L] * rep(X@perm - 1L, each = dY[1L]) +
                  rep.int(Y@perm, dX[1L])
              if(make.dimnames &&
                 !is.null(dnr <- .kroneckerDimnames(X@Dimnames, dX,
                                                    Y@Dimnames, dY)))
                  r@Dimnames <- dnr
              r
	  })

tmp <- function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
    kronecker(as(X, "TsparseMatrix"), Y,
	      FUN = FUN, make.dimnames = make.dimnames, ...)
}
if(FALSE) {
setMethod("kronecker", signature(X="diagonalMatrix", Y="ANY"          ), tmp)
setMethod("kronecker", signature(X="diagonalMatrix", Y="Matrix"       ), tmp)
}
setMethod("kronecker", signature(X="ANY",            Y="sparseMatrix" ), tmp)
## the above could recurse infinitely :
setMethod("kronecker", signature(X="sparseMatrix",   Y="TsparseMatrix"), tmp)

tmp <- function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
    kronecker(X, as(Y, "TsparseMatrix"),
	      FUN = FUN, make.dimnames = make.dimnames, ...)
}
if(FALSE) {
setMethod("kronecker", signature(X="ANY",           Y="diagonalMatrix"), tmp)
setMethod("kronecker", signature(X="Matrix",        Y="diagonalMatrix"), tmp)
}
setMethod("kronecker", signature(X="sparseMatrix",  Y="ANY"           ), tmp)
setMethod("kronecker", signature(X="TsparseMatrix", Y="sparseMatrix"  ), tmp)
rm(tmp)

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
	      ## faster than x = as.vector(outer(Y@x, X@x, FUN = FUN)
	      x = as.vector(Y@x %*% t(X@x)))
      })

## triangularity -- should be preserved "when obvious":
setMethod("kronecker", signature(X = "dtTMatrix", Y = "dtTMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
      {
	  if (FUN != "*") stop("kronecker method must use default 'FUN'")
	  ## otherwise we don't know that many results will be zero
	  if(X@uplo != Y@uplo) { ## result not triangular
	      X <- .sparse2g(X)
	      Y <- .sparse2g(Y)
	      return(callGeneric())
	  }
	  ## else: both 'uplo' are the same -- result *is* triangular
	  ## d.U <- (dX <- X@diag == "U") && (dY <- Y@diag == "U")
	  if(Y@diag == "U")
	      Y <- .diagU2N(Y, "dtTMatrix")
	  ydim <- Y@Dim
	  if(X@diag != "U") {
	      xi <- X@i
	      xj <- X@j
	      xx <- X@x
	  } else { ## X@diag == "U"
	      nx <- X@Dim[1] # triangular matrices are square
	      ii <- seq_len(nx) - 1L
	      xi <- c(X@i, ii)
	      xj <- c(X@j, ii)
	      xx <- c(X@x, rep.int(1, nx))
	  }
	  xnnz <- length(xi)
	  yi <- Y@i
	  ynnz <- length(yi)
	  new("dtTMatrix", Dim = X@Dim * ydim,
	      i = rep.int(yi,  xnnz) + ydim[1] * rep.int(xi, rep.int(ynnz, xnnz)),
	      j = rep.int(Y@j, xnnz) + ydim[2] * rep.int(xj, rep.int(ynnz, xnnz)),
	      ## faster than x = as.vector(outer(Y@x, X@x, FUN = FUN)
	      x = as.vector(Y@x %*% t(xx)),
	      uplo = X@uplo,
	      diag = "N" # if(d.U) { "U" , but drop the entries}  else "N"
	      )
      })

setMethod("kronecker", signature(X = "dtTMatrix", Y = "dgTMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(it <- isTriangular(Y))
                  ## improve: also test for unit diagonal
                  Y <- if(attr(it, "kind") == "U") triu(Y) else tril(Y)
	      else
                  X <- .sparse2g(X)
	      callGeneric() #-> dtT o dtT   or	 dgT o dgT
	  })

setMethod("kronecker", signature(X = "dgTMatrix", Y = "dtTMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(it <- isTriangular(X))
                  ## improve: also test for unit diagonal
                  X <- if(attr(it, "kind") == "U") triu(X) else tril(X)
	      else
		  Y <- .sparse2g(Y)
	      callGeneric() #-> dtT o dtT   or	 dgT o dgT
	  })

setMethod("kronecker", signature(X = "TsparseMatrix", Y = "TsparseMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(.hasSlot(X, "uplo") && !.hasSlot(X, "diag"))
                  X <- .sparse2g(X)
              if(.hasSlot(Y, "uplo") && !.hasSlot(Y, "diag"))
                  Y <- .sparse2g(Y)
              X <- ..sparse2d(X)
              Y <- ..sparse2d(Y)
              callGeneric()
	  })

setMethod("kronecker", signature(X = "dsparseMatrix", Y = "dsparseMatrix"),
	  function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if(.hasSlot(X, "uplo") && !.hasSlot(X, "diag"))
                  X <- .sparse2g(X)
              if(.hasSlot(Y, "uplo") && !.hasSlot(Y, "diag"))
                  Y <- .sparse2g(Y)
              if(.hasSlot(X, "p"))
                  X <- .CR2T(X)
              if(.hasSlot(Y, "p"))
                  Y <- .CR2T(Y)
              callGeneric()
	  })
