## METHODS FOR CLASS: diagonalMatrix (virtual)
## diagonal matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("band", c(x = "diagonalMatrix"),
          function(x, k1, k2, ...) {
              if(k1 <= 0L && k2 >= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("triu", c(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if(k <= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("tril", c(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if(k >= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("diag", c(x = "diagonalMatrix"),
          function(x = 1, nrow, ncol, names = TRUE) {
              kind <- .M.kind(x)
              r <-
                  if(x@diag != "N") {
                      one <- switch(kind, "n" = , "l" = TRUE, "i" = 1L, "d" = 1, "z" = 1+0i)
                      rep.int(one, x@Dim[1L])
                  } else {
                      y <- x@x
                      if(kind == "n" && anyNA(y)) y | is.na(y) else y
                  }
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 {
                     i <- seq_len(min(x@Dim))
                     identical(nms <- dn[[1L]][i], dn[[2L]][i])
                 })
                  names(r) <- nms
              r
          })

setMethod("diag<-", c(x = "diagonalMatrix"),
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
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"",
                                           typeof(value)),
                                  domain = NA))
                  else
                      switch(typeof(value),
                             logical =,
                             integer =,
                             double = rep_len(as.double(value), n),
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"",
                                           typeof(value)),
                                  domain = NA))
              x@diag <- "N"
              x
          })

setMethod( "t", c(x = "diagonalMatrix"),
          function(x) {
              x@Dimnames <- x@Dimnames[2:1]
              x
          })

setMethod("ct", c(x = "diagonalMatrix"),
          function(x) {
              x@Dimnames <- x@Dimnames[2:1]
              if(is.complex(y <- x@x))
                  x@x <- Conj(y)
              x
          })

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
          function(object, trans = "C", checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              .M.kind(object) != "z" || object@diag != "N" || !identical(trans, "C") ||
                  { x <- object@x; isTRUE(all.equal.numeric(x, Conj(x), ...)) }
          })
