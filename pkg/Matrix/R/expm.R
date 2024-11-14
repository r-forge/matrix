## METHODS FOR GENERIC: expm
## the matrix exponential
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expm", c(x = "denseMatrix"),
          function(x)
              .Call(R_dense_expm, x))

setMethod("expm", c(x = "sparseMatrix"),
          function(x)
              .Call(R_dense_expm, .M2unpacked(x)))

setMethod("expm", c(x = "diagonalMatrix"),
          function(x) {
              x <- .M2kind(x, ",")
              z <- is.complex(y <- x@x)
              r <- new(if (z) "zdiMatrix" else "ddiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <-
                  if (x@diag == "N")
                      exp(y)
                  else rep.int(exp(if (z) 1+0i else 1), d[1L])
              r
          })

setMethod("expm", c(x = "matrix"),
          function(x) {
              d <- dim(x)
              if (d[1L] != d[2L])
                  stop("matrix is not square")
              .Call(R_dense_expm, .m2dense(x, ",ge"))
          })
