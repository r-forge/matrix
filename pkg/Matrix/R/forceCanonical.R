## METHODS FOR GENERIC: forceCanonical
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceCanonical", c(x = "denseMatrix"),
          function(x, check = TRUE, diag = NULL, ...) {
              r <- .Call(R_dense_force_canonical, x, check)
              if (!is.null(diag) && diag == "N" &&
                  .M.shape(x) == "t" && x@diag != "N")
                  r@diag <- "N"
              r
          })

setMethod("forceCanonical", c(x = "CsparseMatrix"),
          function(x, check = TRUE, diag = NULL, ...) {
              if (length(x@i) > (nnz <- x@p[length(x@p)])) {
                  x@i <- x@i[seq_len(nnz)]
                  if (.M.kind(x) != "n")
                  x@x <- x@x[seq_len(nnz)]
              }
              r <- .Call(R_sparse_force_canonical, x, check)
              if (!is.null(diag) && diag == "N" &&
                  .M.shape(x) == "t" && x@diag != "N")
                  diag(r) <- TRUE
              r
          })

setMethod("forceCanonical", c(x = "RsparseMatrix"),
          function(x, check = TRUE, diag = NULL, ...) {
              if (length(x@j) > (nnz <- x@p[length(x@p)])) {
                  x@j <- x@j[seq_len(nnz)]
                  if (.M.kind(x) != "n")
                  x@x <- x@x[seq_len(nnz)]
              }
              r <- .Call(R_sparse_force_canonical, x, check)
              if (!is.null(diag) && diag == "N" &&
                  .M.shape(x) == "t" && x@diag != "N")
                  diag(r) <- TRUE
              r
          })

setMethod("forceCanonical", c(x = "TsparseMatrix"),
          function(x, check = TRUE, diag = NULL, ...) {
              x <- aggregateT(x)
              r <- .Call(R_sparse_force_canonical, x, check)
              if (!is.null(diag) && diag == "N" &&
                  .M.shape(x) == "t" && x@diag != "N")
                  diag(r) <- TRUE
              r
          })

setMethod("forceCanonical", c(x = "diagonalMatrix"),
          function(x, diag = NULL, ...) {
              kind <- .M.kind(x)
              if (if (x@diag == "N")
                      kind != "n" || !anyNA(y <- x@x)
                  else is.null(diag) || diag != "N")
                  return(x)
              r <- new(paste0(kind, "diMatrix"))
              r@Dim <- d <- x@Dim
              r@x <-
                  if (x@diag == "N")
                      y | is.na(y)
                  else rep.int(switch(kind, "z" = 1+0i, "d" = 1, "i" = 1L, TRUE), d[1L])
              r
          })

setMethod("forceCanonical", c(x = "indMatrix"),
          function(x, ...)
              x)

setMethod("forceCanonical", c(x = "sparseVector"),
          function(x, ...) {
              length. <- x@length
              i. <- x@i
              if (is.integer(length.) && is.integer(i.))
                  return(x)
              kind <- .M.kind(x)
              r <- new(paste0(kind, "sparseVector"))
              r@length <-
                  if (is.integer(length.))
                      length.
                  else if (length. - 1 < .Machine[["integer.max"]])
                      as.integer(length.)
                  else trunc(length.)
              r@i <-
                  if (is.integer(i.))
                      i.
                  else if (length(i.) == 0L || i.[length(i.)] - 1 < .Machine[["integer.max"]])
                      as.integer(i.)
                  else trunc(i.)
              if (kind != "n")
                  r@x <- x@x
              r
          })
