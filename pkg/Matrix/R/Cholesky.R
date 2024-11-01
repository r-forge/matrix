## METHODS FOR GENERIC: Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Cholesky", c(A = "generalMatrix"),
          function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              ch@Dimnames <- A@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("Cholesky", c(A = "symmetricMatrix"),
          function(A, ...)
              Cholesky(.M2kind(A, ","), ...))

setMethod("Cholesky", c(A = "triangularMatrix"),
          function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              ch@Dimnames <- A@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("Cholesky", c(A = "diagonalMatrix"),
          function(A, ...)
              Cholesky(.M2kind(A, ","), ...))

setMethod("Cholesky", c(A = "dsyMatrix"),
          function(A, perm = TRUE, tol = -1, ...)
              .Call(poMatrix_trf, A, if(perm) 1L else 2L, perm, tol))

setMethod("Cholesky", c(A = "dspMatrix"),
          function(A, ...)
              .Call(ppMatrix_trf, A, 2L))

setMethod("Cholesky", c(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = !super, super = FALSE,
                   Imult = 0, ...)
              .Call(pCMatrix_trf, A, 2L, perm, !LDL, super, Imult))

setMethod("Cholesky", c(A = "dsRMatrix"),
          function(A, ...)
              Cholesky(.tCRT(A), ...))

setMethod("Cholesky", c(A = "dsTMatrix"),
          function(A, ...)
              Cholesky(.M2C(A), ...))

setMethod("Cholesky", c(A = "ddiMatrix"),
          function(A, ...) {
              if(length(y <- A@x) && (is.na(min.y <- min(y)) || min.y < 0))
                  stop(gettextf("%1$s(%2$s) is undefined: '%2$s' is not positive semidefinite",
                                "Cholesky", "x"),
                       domain = NA)
              n <- (d <- A@Dim)[1L]
              r <- new("dsimplicialCholesky")
              r@Dim <- d
              r@Dimnames <- A@Dimnames
              r@minor <- n
              r@colcount <- r@nz <- rep.int(1L, n)
              r@`next` <- c(seq_len(n), -1L, 0L)
              r@ prev  <- c(n + 1L, s, -1L)
              r@p <- 0:n
              r@i <- s <- seq.int(0L, length.out = n)
              r@x <- if(length(y)) y else rep.int(1, n)
              r@ordering <- 0L
              r@is_ll <- FALSE
              r@is_monotonic <- TRUE
              r
          })

setMethod("Cholesky", c(A = "matrix"),
          function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              if(!is.null(dn <- dimnames(A)))
                  ch@Dimnames <- dn # restore asymmetric 'Dimnames'
              ch
          })
