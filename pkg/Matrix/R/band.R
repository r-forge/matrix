## METHODS FOR GENERIC: band, triu, tril
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("matrix", "denseMatrix")) {
setMethod("band", c(x = .cl),
          function(x, k1, k2, ...)
              .Call(R_dense_band, x, k1, k2))

setMethod("triu", c(x = .cl),
          function(x, k = 0L, ...)
              .Call(R_dense_band, x, k, NULL))

setMethod("tril", c(x = .cl),
          function(x, k = 0L, ...)
              .Call(R_dense_band, x, NULL, k))
}

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix")) {
setMethod("band", c(x = .cl),
          function(x, k1, k2, ...)
              .Call(R_sparse_band, x, k1, k2))

setMethod("triu", c(x = .cl),
          function(x, k = 0L, ...)
              .Call(R_sparse_band, x, k, NULL))

setMethod("tril", c(x = .cl),
          function(x, k = 0L, ...)
              .Call(R_sparse_band, x, NULL, k))
}

setMethod("band", c(x = "diagonalMatrix"),
          function(x, k1, k2, ...) {
              if (k1 <= 0L && k2 >= 0L)
                  return(x)
              r <- new(.M.class(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("triu", c(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if (k <= 0L)
                  return(x)
              r <- new(.M.class(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("tril", c(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if (k >= 0L)
                  return(x)
              r <- new(.M.class(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("band", c(x = "indMatrix"),
          function(x, k1, k2, ...)
              band(.ind2sparse(x), k1, k2, ...))

setMethod("triu", c(x = "indMatrix"),
          function(x, k = 0L, ...)
              triu(.ind2sparse(x), k, ...))

setMethod("tril", c(x = "indMatrix"),
          function(x, k = 0L, ...)
              tril(.ind2sparse(x), k, ...))

rm(.cl)
