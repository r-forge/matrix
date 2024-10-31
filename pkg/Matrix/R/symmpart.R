## METHODS FOR GENERIC: symmpart
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("symmpart", c(x = "matrix"),
          function(x, trans = "C", ...) {
              op <- if (is.complex(x) && identical(trans, "C"))
                        Conj
                    else identity
              symmetrizeDN(0.5 * (x + op(t(x))))
          })

setMethod("symmpart", c(x = "denseMatrix"),
          function(x, uplo = "U", trans = "C", ...)
              .Call(R_dense_symmpart, x, uplo, trans))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("symmpart", c(x = .cl),
          function(x, uplo = "U", trans = "C", ...)
              .Call(R_sparse_symmpart, x, uplo, trans))

setMethod("symmpart", c(x = "diagonalMatrix"),
          function(x, trans = "C", ...) {
              kind <- .M.kind(x)
              r <- new(if (kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- symDN(x@Dimnames)
              if (x@diag != "N")
                  r@diag <- "U"
              else {
                  y <- x@x
                  r@x <-
                  switch(kind,
                         "n" = as.double(y | is.na(y)),
                         "l" =,
                         "i" =,
                         "d" = as.double(y),
                         "z" =
                             if (identical(trans, "C"))
                                 complex(real = Re(y), imaginary = 0)
                             else as.complex(y))
              }
              r
          })

setMethod("symmpart", c(x = "indMatrix"),
          function(x, ...)
              symmpart(.M2kind(x, "d"), ...))

rm(.cl)
