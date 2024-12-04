## METHODS FOR GENERIC: isCanonical
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("isCanonical", c(object = "denseMatrix"),
          function(object, diag = NULL, ...)
              .Call(R_dense_is_canonical, object) &&
              (is.null(diag) || diag != "N" ||
               .M.shape(object) != "t" || object@diag == "N"))

setMethod("isCanonical", c(object = "CsparseMatrix"),
          function(object, diag = NULL, ...)
              length(object@i) == object@p[length(object@p)] &&
              .Call(R_sparse_is_canonical, object) &&
              (is.null(diag) || diag != "N" ||
               .M.shape(object) != "t" || object@diag == "N"))

setMethod("isCanonical", c(object = "RsparseMatrix"),
          function(object, diag = NULL, ...)
              length(object@j) == object@p[length(object@p)] &&
              .Call(R_sparse_is_canonical, object) &&
              (is.null(diag) || diag != "N" ||
               .M.shape(object) != "t" || object@diag == "N"))

setMethod("isCanonical", c(object = "TsparseMatrix"),
          function(object, diag = NULL, ...)
              !anyDuplicatedT(object) &&
              .Call(R_sparse_is_canonical, object) &&
              (is.null(diag) || diag != "N" ||
               .M.shape(object) != "t" || object@diag == "N"))

setMethod("isCanonical", c(object = "diagonalMatrix"),
          function(object, diag = NULL, ...)
              (.M.kind(object) != "n" || !anyNA(object@x)) &&
              (is.null(diag) || diag != "N" || object@diag == "N"))

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
