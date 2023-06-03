## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", signature(x = "dsyMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(dsyMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", signature(x = "dspMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(dspMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", signature(x = "matrix"),
          function(x, uplo = "U", ...)
              BunchKaufman(.m2dense(x, "dsy", uplo), ...))


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

.def.unpacked <- .def.packed <- function(x, ...) {
    r <- .Call(BunchKaufman_expand, x, .PACKED)
    dn <- x@Dimnames
    if(b <- length(r) - 1L)
        r <- if(x@uplo == "U")
                 c(lapply(r[(b + 1L):2], t), r)
             else c(r, lapply(r[b:1L], t))
    r[[1L]]@Dimnames <- c(dn[1L], list(NULL))
    r[[length(r)]]@Dimnames <- c(list(NULL), dn[2L])
    r
}
body(.def.unpacked) <-
    do.call(substitute, list(body(.def.unpacked), list(.PACKED = FALSE)))
body(.def.packed) <-
    do.call(substitute, list(body(.def.packed  ), list(.PACKED =  TRUE)))

## returning:
##
## list(U[1]', P[1]', ..., U[b]', P[b]', D, P[b], U[b], ..., P[1], U[1])
##     where A = U' D U and U = P[b] U[b] ... P[1] U[1]
##
## OR
##
## list(P[1], L[1], ..., P[b], L[b], D, L[b]', P[b]', ..., L[1]', P[1]')
##     where A = L D L' and L = P[1] L[1] ... P[b] L[b]
setMethod("expand2", signature(x =  "BunchKaufman"), .def.unpacked)
setMethod("expand2", signature(x = "pBunchKaufman"), .def.packed)

rm(.def.unpacked, .def.packed)


if(FALSE) {
## FIXME: but expansion is wrong for uplo = "U" ??
library(Matrix)
set.seed(1)
n <- 6L

X <- new("dsyMatrix", Dim = c(n, n), x = rnorm(n * n))
Y <- t(X) # same but with uplo = "L"

bkX <- BunchKaufman(X)
bkY <- BunchKaufman(Y)

bkX. <- as(bkX, "dtrMatrix")
bkY. <- as(bkY, "dtrMatrix")

eX <- expand2(bkX)
eY <- expand2(bkY)

## 'eX' seems correct ... but is actually wrong {first test fails}
stopifnot(exprs = {
    all.equal(Reduce(`%*%`, lapply(eX, as, "matrix")), as(X, "matrix"))
    all.equal(Reduce(`%*%`, lapply(eY, as, "matrix")), as(Y, "matrix"))
})
}
