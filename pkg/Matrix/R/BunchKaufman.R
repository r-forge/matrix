## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
	  function(x, ...) .Call(dsyMatrix_trf, x, 2L))

setMethod("BunchKaufman", signature(x = "dspMatrix"),
	  function(x, ...) .Call(dspMatrix_trf, x, 2L))

setMethod("BunchKaufman", signature(x = "matrix"),
	  function(x, uplo = "U", ...) {
              if(is.logical(x) || is.integer(x))
                  storage.mode(x) <- "double"
              .Call(matrix_trf, x, 2L, uplo)
          })


## METHODS FOR CLASS: p?BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("BunchKaufman", "dtrMatrix",
      function(from) {
          to <- new("dtrMatrix")
          to@Dim <- from@Dim
          to@Dimnames <- from@Dimnames
          to@uplo <- from@uplo
          to@x <- from@x
          to
      })

setAs("pBunchKaufman", "dtpMatrix",
      function(from) {
          to <- new("dtpMatrix")
          to@Dim <- from@Dim
          to@Dimnames <- from@Dimnames
          to@uplo <- from@uplo
          to@x <- from@x
          to
      })

## FIXME: need all 2*b+1 factors, not only b+1; and need Dimnames
##
## returning:
##
## list(D, P[1], U[1], ..., P[b], U[b])
##     where A = U' D U and U = P[1] U[1] ... P[b] U[b]
##
## OR
##
## list(P[1], L[1], ..., P[b], L[b], D)
##     where A = L D L' and L = P[1] L[1] ... P[b] L[b]
##
## as described in the documentation for LAPACK 'ds[yp]trf'
setMethod("expand2", signature(x = "BunchKaufman"),
          function(x, ...) .Call(BunchKaufman_expand, x))

setMethod("expand2", signature(x = "pBunchKaufman"),
          function(x, ...) .Call(BunchKaufman_expand, x))

if(FALSE) {
library(Matrix)
set.seed(1)

n <- 1000L
X <- new("dsyMatrix", Dim = c(n, n), x = rnorm(n * n))
Y <- t(X)

as(bkX <- BunchKaufman(X), "dtrMatrix")
as(bkY <- BunchKaufman(Y), "dtrMatrix")

DU <- .Call("BunchKaufman_expand", bkX)
D <- DU[[1L]]
U <- Reduce(`%*%`, DU[-1L])
## FIXME: 'DU' looks correct ... but is actually wrong {second test fails}??
stopifnot(all.equal(as(t(U) %*% D %*% U, "matrix"), as(X, "matrix")))

LD <- .Call("BunchKaufman_expand", bkY)
D <- LD[[length(LD)]]
L <- Reduce(`%*%`, LD[-length(LD)])
stopifnot(all.equal(as(L %*% D %*% t(L), "matrix"), as(Y, "matrix")))
}
