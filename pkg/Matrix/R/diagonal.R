## METHODS FOR CLASS: diagonalMatrix (virtual)
## diagonal matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceSymmetric", c(x = "diagonalMatrix"),
          function(x, uplo = "U", trans = "C", ...)
              .diag2sparse(x, ".", "s", "C", uplo, trans))

setMethod("symmpart", c(x = "diagonalMatrix"),
          function(x, trans = "C", ...) {
              kind <- .M.kind(x)
              r <- new(if(kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- symDN(x@Dimnames)
              if(x@diag != "N")
                  r@diag <- "U"
              else {
                  y <- x@x
                  r@x <- switch(kind,
                                "n" = as.double(y | is.na(y)),
                                "l" = ,
                                "i" = ,
                                "d" = as.double(y),
                                "z" = if(identical(trans, "C")) complex(real = Re(y), imaginary = 0) else as.complex(y))
              }
              r
          })

setMethod("skewpart", c(x = "diagonalMatrix"),
          function(x, trans = "C", ...) {
              kind <- .M.kind(x)
              r <- new(if(kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- symDN(x@Dimnames)
              r@x <-
                  if(kind != "z")
                      double(d[1L])
                  else if(x@diag != "N" || !identical(trans, "C"))
                      complex(d[1L])
                  else complex(real = 0, imaginary = Im(x@x))
              r
          })
