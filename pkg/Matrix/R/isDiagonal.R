## METHODS FOR GENERIC: isDiagonal
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("matrix", "denseMatrix"))
setMethod("isDiagonal", c(object = .cl),
          function(object, ...)
              .Call(R_dense_is_diagonal, object))

setMethod("isDiagonal", c(object = "table"),
          function(object, ...) {
              if (length(dim(object)) != 2L)
                  return(FALSE)
              .Call(R_dense_is_diagonal, object)
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("isDiagonal", c(object = .cl),
          function(object, ...)
              .Call(R_sparse_is_diagonal, object))

setMethod("isDiagonal", c(object = "diagonalMatrix"),
          function(object, ...)
              TRUE)

setMethod("isDiagonal", c(object = "indMatrix"),
          function(object, ...) {
              d <- object@Dim
              if ((n <- d[2L]) != d[1L])
                  return(FALSE)
              all(object@perm == seq_len(n))
          })

rm(.cl)
