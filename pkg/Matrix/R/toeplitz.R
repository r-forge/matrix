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
              nn <- c(n, n)
              r <- spV2M(x[as.integer(abs(.col(nn) - .row(nn))) + 1L],
                         nrow = n, ncol = n, symmetric = symmetric,
                         check = FALSE)
              repr <- # keep in sync with sparseMatrix
              if (missing(giveCsparse))
                  match.arg(repr)
              else if (!missing(repr)) {
                  warning(gettextf("'%s' is deprecated; using '%s' instead",
                                   "giveCsparse", "repr"),
                          domain = NA)
                  match.arg(repr)
              }
              else if (giveCsparse)
                  "C"
              else {
                  warning(gettextf("'%s' is deprecated; using %s=\"%s\"",
                                   "giveCsparse", "repr", "T"),
                          domain = NA)
                  "T"
              }
              switch(repr, "C" = .M2C(r), "R" = .M2R(r), "T" = r)
          })
