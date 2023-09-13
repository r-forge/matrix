## METHODS FOR CLASS: diagonalMatrix (virtual)
## diagonal matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("band", signature(x = "diagonalMatrix"),
          function(x, k1, k2, ...) {
              if(k1 <= 0L && k2 >= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("triu", signature(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if(k <= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("tril", signature(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if(k >= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("diag", signature(x = "diagonalMatrix"),
          function(x, nrow, ncol, names = TRUE) {
              r <- .diag.x(x)
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 {
                     i <- seq_len(min(x@Dim))
                     identical(nms <- dn[[1L]][i], dn[[2L]][i])
                 })
                  names(r) <- nms
              r
          })

setMethod("diag<-", signature(x = "diagonalMatrix"),
          function(x, value) {
              n <- x@Dim[2L]
              nv <- length(value)
              if(nv != 1L && nv != n)
                  stop("replacement diagonal has wrong length")
              x@x <-
                  if(is.logical(x@x))
                      switch(typeof(value),
                             logical = rep_len(value, n),
                             integer =,
                             double =
                                 {
                                     x <- .M2kind(x, "d")
                                     rep_len(as.double(x), n)
                                 },
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"", typeof(value)),
                                  domain = NA))
                  else
                      switch(typeof(value),
                             logical =,
                             integer =,
                             double = rep_len(as.double(value), n),
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"", typeof(value)),
                                  domain = NA))
              x@diag <- "N"
              x
          })

setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1]; x })

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "missing"),
          function(x, uplo) .diag2sparse(x, "s", "C",  "U"))

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "character"),
          function(x, uplo) .diag2sparse(x, "s", "C", uplo))

setMethod("symmpart", signature(x = "diagonalMatrix"),
          function(x) {
              kind <- .M.kind(x)
              r <- new(if(kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- symmDN(x@Dimnames)
              if(x@diag != "N")
                  r@diag <- "U"
              else
                  r@x <- switch(kind,
                                "n" = as.double(x@x | is.na(x@x)),
                                "l" = ,
                                "i" = ,
                                "d" = as.double(x@x),
                                "z" = complex(real = Re(x@x), imaginary = 0))
              r
          })

setMethod("skewpart", signature(x = "diagonalMatrix"),
          function(x) {
              kind <- .M.kind(x)
              r <- new(if(kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- symmDN(x@Dimnames)
              r@x <- if(kind != "z")
                         double(d[1L])
                     else if(x@diag != "N")
                         complex(d[1L])
                     else complex(real = 0, imaginary = Im(x@x))
              r
          })

setMethod("isDiagonal", signature(object = "diagonalMatrix"),
          function(object) TRUE)

setMethod("isTriangular", signature(object = "diagonalMatrix"),
          function(object, upper = NA, ...)
              if(is.na(upper)) `attr<-`(TRUE, "kind", "U") else TRUE)

setMethod("isSymmetric", signature(object = "diagonalMatrix"),
          function(object, checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              TRUE
          })
