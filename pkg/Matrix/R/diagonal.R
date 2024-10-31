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

setMethod("isDiagonal", c(object = "diagonalMatrix"),
          function(object) TRUE)

setMethod("isTriangular", c(object = "diagonalMatrix"),
          function(object, upper = NA)
              if(is.na(upper)) `attr<-`(TRUE, "kind", "U") else TRUE)

setMethod("isSymmetric", c(object = "diagonalMatrix"),
          function(object,
                   tol = 100 * .Machine$double.eps,
                   trans = "C", checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              ae <- function(target, current, tolerance, scale = NULL, ...)
                  all.equal.numeric(target = target, current = current,
                                    tolerance = tolerance, scale = scale,
                                    check.attributes = FALSE, check.class = FALSE)
              !is.complex(x <- object@x) || !identical(trans, "C") ||
                  object@diag != "N" || isTRUE(ae(x, Conj(x), tolerance = tol, ...))
          })
