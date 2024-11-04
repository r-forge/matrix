## METHODS FOR GENERIC: toeplitz
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("toeplitz", c(x = "sparseVector"),
          function(x, symmetric = TRUE, repr = c("C", "R", "T"),
                   giveCsparse, ...) {
              n <- length(x)
              if (n > .Machine[["integer.max"]])
                  stop(gettextf("dimensions cannot exceed %s",
                                "2^31-1"),
                       domain = NA)
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
              d <- c(n, n)
              .V2sparse(x[as.integer(abs(.col(d) - .row(d))) + 1L],
                        paste0(".", if (symmetric) "s" else "g", repr),
                        uplo = "U", trans = "T", nrow = n, ncol = n)
          })
