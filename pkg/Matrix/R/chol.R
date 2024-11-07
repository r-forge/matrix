## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", c(x = "generalMatrix"),
          function(x, uplo = "U", ...) {
              ch <- chol(forceSymmetric(x, uplo), ...)
              ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("chol", c(x = "symmetricMatrix"),
          function(x, ...)
              chol(.M2kind(x, ","), ...))

setMethod("chol", c(x = "triangularMatrix"),
          function(x, uplo = "U", ...) {
              if (identical(uplo, x@uplo)) {
                  ch <- chol(forceSymmetric(x, uplo), ...)
                  ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
                  ch
              } else chol(forceDiagonal(x, diag = x@diag), ...)
          })

setMethod("chol", c(x = "diagonalMatrix"),
          function(x, ...)
              chol(.M2kind(x, ","), ...))

setMethod("chol", c(x = "dsyMatrix"),
          function(x, pivot = FALSE, tol = -1, ...) {
              ch <- as(Cholesky(x, perm = pivot, tol = tol), "dtrMatrix")
              ch@Dimnames <- dimnames(x)
              if (ch@uplo != "U") t(ch) else ch
          })

setMethod("chol", c(x = "dspMatrix"),
          function(x, ...) {
              ch <- as(Cholesky(x), "dtpMatrix")
              ch@Dimnames <- dimnames(x)
              if (ch@uplo != "U") t(ch) else ch
          })

for(.cl in paste0("ds", c("C", "R", "T"), "Matrix"))
setMethod("chol", c(x = .cl),
          function(x, pivot = FALSE, ...) {
              ch <- t(as(Cholesky(x, perm = pivot, LDL = FALSE, super = FALSE),
                         "dtCMatrix")) # FIXME? give dtRMatrix, dtTMatrix?
              ch@Dimnames <- dimnames(x)
              ch
          })
rm(.cl)

setMethod("chol", c(x = "ddiMatrix"),
          function(x, ...) {
              if (length(y <- x@x)) {
                  if (is.na(min.y <- min(y)) || min.y < 0)
                      stop(gettextf("%1$s(%2$s) is undefined: '%2$s' is not positive semidefinite",
                                    "chol", "x"),
                           domain = NA)
                  x@x <- sqrt(y)
              }
              x
          })
