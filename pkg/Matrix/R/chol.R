## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", c(x = "denseMatrix"),
          function(x, pivot = FALSE, tol = -1, ...) {
              trf <- Cholesky(x, perm = pivot, tol = tol)
              r <- expand1(trf, "L.")
              dimnames(r) <- dimnames(x)
              r
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("chol", c(x = .cl),
          function(x, pivot = FALSE, ...) {
              trf <- Cholesky(x, perm = pivot, LDL = FALSE, super = FALSE)
              r <- expand1(trf, "L.")

              dimnames(r) <- dimnames(x)
              r
          })

setMethod("chol", c(x = "diagonalMatrix"),
          function(x, direct = TRUE, ...) {
              if (direct) {
                  z <- is.complex(x@x)
                  r <- new(if (z) "zdiMatrix" else "ddiMatrix")
                  r@Dim <- d <- x@Dim
                  r@Dimnames <- x@Dimnames
                  if (x@diag != "N")
                      r@diag <- "U"
                  else if (d[1L] > 0L) {
                      x <- .M2kind(x, ",")
                      y <- x@x
                      if (z && any(is.na(tmp <- range(Im(y))) | tmp != 0))
                          stop("matrix is not Hermitian")
                      if (is.na(tmp <- min(if (z) Re(y) else y)) || tmp < 0)
                          stop("matrix is not positive semidefinite")
                      r@x <- sqrt(y)
                  }
                  r
              } else chol(.diag2sparse(x, ",", "s", "C"), ...)
          })

setMethod("chol", c(x = "indMatrix"),
          function(x, ...)
              chol(.ind2sparse(x, ",", "C"), ...))
