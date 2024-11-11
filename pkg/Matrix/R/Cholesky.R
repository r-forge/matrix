## METHODS FOR GENERIC: Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Cholesky", c(A = "denseMatrix"),
          function(A, perm = TRUE, tol = -1, uplo = "U", ...)
              .Call(R_dense_cholesky, A, if (perm) 1L else 2L, perm, tol, uplo))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("Cholesky", c(A = .cl),
          function(A, perm = TRUE, LDL = !super, super = FALSE,
                   Imult = 0, force = TRUE, uplo = "U", ...)
              .Call(R_sparse_cholesky, A, 2L, perm, !LDL, super, Imult, force, uplo))

setMethod("Cholesky", c(A = "diagonalMatrix"),
          function(A, direct = TRUE, ...) {
              if (direct) {
                  A <- .M2kind(A, ",")
                  y <- A@x
                  z <- is.complex(y)
                  if (z && length(y) &&
                      any(is.na(tmp <- range(Im(y))) | tmp != 0))
                      stop("matrix is not Hermitian")
                  n <- (d <- A@Dim)[1L]
                  r <- new(if (z) "zsimplicialCholesky" else "dsimplicialCholesky")
                  r@Dim <- d
                  r@Dimnames <- A@Dimnames
                  r@minor <- n
                  r@colcount <- r@nz <- rep.int(1L, n)
                  r@`next` <- c(seq_len(n), -1L, 0L)
                  r@ prev  <- c(n + 1L, s, -1L)
                  r@p <- 0L:n
                  r@i <- s <- seq.int(0L, length.out = n)
                  r@x <- if (length(y)) y else rep.int(if (z) 1+0i else 1, n)
                  r@ordering <- 0L
                  r@is_ll <- FALSE
                  r@is_monotonic <- TRUE
              } else Cholesky(.diag2sparse(A, ",", "s", "C"), ...)
          })

setMethod("Cholesky", c(A = "indMatrix"),
          function(A, ...)
              Cholesky(.ind2sparse(A, ",", "C"), ...))

setMethod("Cholesky", c(A = "matrix"),
          function(A, ...)
              Cholesky(.m2dense(A, ",ge"), ...))

rm(.cl)
