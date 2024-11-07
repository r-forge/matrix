## METHODS FOR GENERIC: isCanonical
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.isCanonical <-
function(x)
    switch(.M.repr(x),
           "n" =, "p" =
               .Call(R_dense_is_canonical, x),
           "C" =, "R" =, "T" =
               .Call(R_sparse_is_canonical, x),
           stop("should never happen ..."))

setMethod("isCanonical", c(x = "denseMatrix"),
          function(x, ...)
              .Call(R_dense_is_canonical, x))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("isCanonical", c(x = .cl),
          function(x, ...)
              .Call(R_sparse_is_canonical, x))

setMethod("isCanonical", c(x = "diagonalMatrix"),
          function(x, ...)
              x@diag == "N" && (.M.kind(x) != "n" || !anyNA(x@x)))

setMethod("isCanonical", c(x = "indMatrix"),
          function(x, ...)
              TRUE)

setMethod("isCanonical", c(x = "sparseVector"),
          function(x, ...)
              (is.integer(length. <- x@length) ||
               (length. - 1 >= .Machine[["integer.max"]] &&
                length. == trunc(length.))) &&
              (is.integer(i. <- x@i) ||
               (length(i.) > 0L &&
                i.[length(i.)] - 1 >= .Machine[["integer.max"]] &&
                all(i. == trunc(i.)))))

rm(.cl)
