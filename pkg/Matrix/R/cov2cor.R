## METHODS FOR GENERIC: cov2cor
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("cov2cor", c(V = "denseMatrix"),
          function(V) {
              d <- V@Dim
              if (d[1L] != d[2L] || .M.kind(V) == "z")
                  stop(gettextf("'%s' is not a square numeric matrix",
                                "V"),
                       domain = NA)
              as(forceSymmetric(V),
                 if (!.isPacked(V)) "corMatrix" else "copMatrix")
          })

setMethod("cov2cor", c(V = "sparseMatrix"),
          function(V) {
              d <- V@Dim
              if (d[1L] != d[2L] || .M.kind(V) == "z")
                  stop(gettextf("'%s' is not a square numeric matrix",
                                "V"),
                       domain = NA)
              dn <- symDN(V@Dimnames)
              V <- .M2kind(V, "d")
              Vjj <- diag(V, names = FALSE)
              if (length(Vjj) > 0L && is.na(m <- min(Vjj)) || m <= 0) {
                  warning(gettextf("non-positive or non-finite entries in diag(%s) replaced with NaN",
                                   "V"),
                          domain = NA)
                  Vjj[is.na(Vjj) | Vjj < 0] <- NaN
              }
              D <- Diagonal(x = sqrt(1/Vjj))
              r <- forceSymmetric(D %*% V %*% D)
              diag(r) <- 1
              r@Dimnames <- dn
              r
          })
