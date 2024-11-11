## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", c(x = "denseMatrix"),
          function(x, pivot = FALSE, tol = -1, ...) {
              trf <- .Call(R_dense_cholesky, x, if (pivot) 1L else 2L, pivot, tol, NULL)
              r <- expand1(trf, "L.")
              dimnames(r) <- dimnames(x)
              r
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("chol", c(x = .cl),
          function(x, pivot = FALSE, ...) {
              trf <- .Call(R_sparse_cholesky, x, 2L, pivot, FALSE, FALSE, 0, TRUE, NULL)
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
