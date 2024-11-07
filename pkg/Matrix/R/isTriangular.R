## METHODS FOR GENERIC: isTriangular
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("matrix", "denseMatrix"))
setMethod("isTriangular", c(object = .cl),
          function(object, upper = NA, ...)
              .Call(R_dense_is_triangular, object, upper))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("isTriangular", c(object = .cl),
          function(object, upper = NA, ...)
              .Call(R_sparse_is_triangular, object, upper))

setMethod("isTriangular", c(object = "diagonalMatrix"),
          function(object, upper = NA, ...)
              if (is.na(upper)) `attr<-`(TRUE, "kind", "U") else TRUE)

setMethod("isTriangular", c(object = "indMatrix"),
          function(object, upper = NA, ...) {
              d <- object@Dim
              if ((n <- d[2L]) != d[1L])
                  return(FALSE)
              if (object@margin == 1L) {
                  i <- seq_len(n)
                  j <- object@perm
              } else {
                  i <- object@perm
                  j <- seq_len(n)
              }
              if (is.na(upper)) {
                  if (all(j >= i))
                      return(`attr<-`(TRUE, "kind", "U"))
                  if (all(i <= j))
                      return(`attr<-`(TRUE, "kind", "L"))
                  FALSE
              }
              else if (upper)
                  all(j >= i)
              else
                  all(i <= j)
          })

rm(.cl)
