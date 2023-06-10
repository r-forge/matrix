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
          to@uplo <- from@uplo
          to@x <- from@x
          to
      })

setAs("pBunchKaufman", "dtpMatrix",
      function(from) {
          to <- new("dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          to
      })

.def.unpacked <- .def.packed <- function(x, ...) {
    r <- .Call(BunchKaufman_expand, x, .PACKED)
    dn <- x@Dimnames
    if(b <- length(r) - 1L) {
        r <- c(r, lapply(r[b:1L], t))
        r[[1L]]@Dimnames <- c(dn[1L], list(NULL))
        r[[length(r)]]@Dimnames <- c(list(NULL), dn[2L])
    } else r@Dimnames <- dn
    r
}
body(.def.unpacked) <-
    do.call(substitute, list(body(.def.unpacked), list(.PACKED = FALSE)))
body(.def.packed) <-
    do.call(substitute, list(body(.def.packed  ), list(.PACKED =  TRUE)))

## returning:
##
## list(P[b], U[b], ..., P[1], U[1], D, U[1]', P[1]', ..., U[b]', P[b]')
##     where A = U D U' and U = P[b] U[b] ... P[1] U[1]
##
## OR
##
## list(P[1], L[1], ..., P[b], L[b], D, L[b]', P[b]', ..., L[1]', P[1]')
##     where A = L D L' and L = P[1] L[1] ... P[b] L[b]
setMethod("expand2", signature(x =  "BunchKaufman"), .def.unpacked)
setMethod("expand2", signature(x = "pBunchKaufman"), .def.packed)

rm(.def.unpacked, .def.packed)
