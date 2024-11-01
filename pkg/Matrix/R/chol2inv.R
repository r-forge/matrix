## METHODS FOR GENERIC: chol2inv
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol2inv", c(x = "generalMatrix"),
          function(x, uplo = "U", ...) {
              d <- x@Dim
              if (d[1L] != d[2L])
                  stop("matrix is not square")
              chol2inv((if (  uplo == "U") triu else tril)(x), ...)
          })

setMethod("chol2inv", c(x = "symmetricMatrix"),
          function(x, ...)
              chol2inv((if (x@uplo == "U") triu else tril)(x), ...))

setMethod("chol2inv", c(x = "triangularMatrix"),
          function(x, ...)
              chol2inv(.M2kind(x, ","), ...))

setMethod("chol2inv", c(x = "diagonalMatrix"),
          function(x, ...)
              chol2inv(.M2kind(x, ","), ...))

for (.cl in paste0("dt", c("r", "p"), "Matrix"))
setMethod("chol2inv", c(x = .cl),
          function(x, ...) {
              if (x@diag != "N")
                  x <- ..diagU2N(x)
              r <- .Call(denseCholesky_solve, x, NULL)
              i <- if (x@uplo == "U") 2L else 1L
              r@Dimnames <- x@Dimnames[c(i, i)]
              r
          })

for (.cl in paste0("dt", c("C", "R", "T"), "Matrix"))
setMethod("chol2inv", c(x = .cl),
          function(x, ...)
              (if (x@uplo == "U") tcrossprod else crossprod)(solve(x)))

## Argument 'uplo' can affect the 'Dimnames' of the result here :

setMethod("chol2inv", c(x = "ddiMatrix"),
          function(x, uplo = "U", ...)
              (if (  uplo == "U") tcrossprod else crossprod)(solve(x)))

rm(.cl)
