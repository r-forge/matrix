## METHODS FOR GENERIC: toeplitz
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("toeplitz", c(x = "sparseVector"),
          function(x, r = NULL, symmetric = is.null(r),
                   repr = c("C", "R", "T"), giveCsparse, ...) {
              m <- n <- length(x)
              if (!symmetric && !is.null(r)) {
                  stopifnot(is(r, "sparseVector"))
                  n <- length(r)
              }
              if (is.double(m) || is.double(n))
                  stop(gettextf("dimensions cannot exceed %s",
                                "2^31-1"),
                       domain = NA)
              d <- c(m, n)
              y <-
              if (symmetric || is.null(r))
                  x[as.integer(abs(.row(d) - .col(d))) + 1L]
              else c(r[if (n >= 2L) n:2L else 0L], x)[.row(d) - .col(d) + n]
              
              repr <-
              if (!missing(repr) || missing(giveCsparse))
                  match.arg(repr)
              else if (giveCsparse)
                  "C"
              else "T"
              if (!missing(giveCsparse))
                  warning(gettextf("'%s' is deprecated; using %s=\"%s\"",
                                   "giveCsparse", "repr", repr),
                          domain = NA)
              .V2sparse(y, paste0(".", if (symmetric) "s" else "g", repr),
                        uplo = "U", trans = "T", nrow = m, ncol = n)
          })
