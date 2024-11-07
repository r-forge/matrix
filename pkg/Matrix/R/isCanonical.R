## METHODS FOR GENERIC: isCanonical
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.isCanonical <-
function(object)
    switch(.M.repr(object),
           "n" =, "p" =
               .Call(R_dense_is_canonical, object),
           "C" =, "R" =, "T" =
               .Call(R_sparse_is_canonical, object),
           stop("should never happen ..."))

setMethod("isCanonical", c(object = "denseMatrix"),
          function(object, ...)
              .Call(R_dense_is_canonical, object))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("isCanonical", c(object = .cl),
          function(object, ...)
              .Call(R_sparse_is_canonical, object))

setMethod("isCanonical", c(object = "diagonalMatrix"),
          function(object, ...)
              object@diag == "N" &&
              (.M.kind(object) != "n" || !anyNA(object@x)))

setMethod("isCanonical", c(object = "indMatrix"),
          function(object, ...)
              TRUE)

setMethod("isCanonical", c(object = "sparseVector"),
          function(object, ...)
              (is.integer(length. <- object@length) ||
               (length. - 1 >= .Machine[["integer.max"]] &&
                length. == trunc(length.))) &&
              (is.integer(i. <- object@i) ||
               (length(i.) > 0L &&
                i.[length(i.)] - 1 >= .Machine[["integer.max"]] &&
                all(i. == trunc(i.)))))

rm(.cl)
