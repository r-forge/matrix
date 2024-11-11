## METHODS FOR GENERIC: chol2inv
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol2inv", c(x = "denseMatrix"),
          function(x, uplo = "U", ...) {
              d <- x@Dim
              switch(.M.shape(x),
                     "g" =
                         if (d[1L] != d[2L])
                             stop("matrix is not square"),
                     "s" =
                         uplo <- x@uplo,
                     "t" =
                         {
                             uplo <- x@uplo
                             if (x@diag != "N")
                                 diag(x) <- TRUE
                         })
              i <- if (uplo == "U") 2L else 1L
              x <- .M2kind(x, ",")
              y <- x@x
              z <- is.complex(y)
              trf <- new(if (z) "zdenseCholesky" else "ddenseCholesky")
              trf@Dim <- d
              trf@Dimnames <- dimnames(x)[c(i, i)]
              trf@uplo <- uplo
              trf@x <- y
              .Call(denseCholesky_solve, trf, NULL)
          })

setMethod("chol2inv", c(x = "sparseMatrix"),
          function(x, uplo = "U", ...) {
              d <- x@Dim
              switch(.M.shape(x),
                     "g" =, "i" =
                         if (d[1L] != d[2L])
                             stop("matrix is not square"),
                     "s" =, "t" =
                         uplo <- x@uplo)
              if (uplo == "U")
                  tcrossprod(solve(triu(x)))
              else crossprod(solve(tril(x)))
          })
