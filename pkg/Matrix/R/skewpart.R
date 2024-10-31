## METHODS FOR GENERIC: skewpart
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("skewpart", c(x = "matrix"),
          function(x, trans = "C", ...) {
              op <- if (is.complex(x) && identical(trans, "C"))
                        Conj
                    else identity
              symmetrizeDN(0.5 * (x - op(t(x))))
          })

setMethod("skewpart", c(x = "denseMatrix"),
          function(x, trans = "C", ...)
              .Call(R_dense_skewpart, x, trans))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("skewpart", c(x = .cl),
          function(x, trans = "C", ...)
              .Call(R_sparse_skewpart, x, trans))

setMethod("skewpart", c(x = "diagonalMatrix"),
          function(x, trans = "C", ...) {
              kind <- .M.kind(x)
              r <- new(if (kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- symDN(x@Dimnames)
              r@x <-
              if (kind != "z")
                  double(d[1L])
              else if (x@diag != "N" || !identical(trans, "C"))
                  complex(d[1L])
              else complex(real = 0, imaginary = Im(x@x))
              r
          })

setMethod("skewpart", c(x = "indMatrix"),
          function(x, ...)
              skewpart(.M2kind(x, "d"), ...))

rm(.cl)
