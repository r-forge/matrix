## METHODS FOR GENERIC: mean
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("mean", c(x = "denseMatrix"),
          function(x, ...) mean(.M2v(x), ...))

setMethod("mean", c(x = "sparseMatrix"),
          function(x, ...) mean(.M2V(x), ...))

setMethod("mean", c(x = "sparseVector"),
          function(x, trim = 0, na.rm = FALSE, ...) {
              kind <- .M.kind(x)
              if (kind == "z" && trim > 0)
                  stop("trimmed means are not defined for complex data")
              n <- length(x)
              if (kind != "n" && n > 0L && anyNA(x@x)) {
                  if (!na.rm)
                      return(NA_real_)
                  n <- n - sum(is.na(x@x))
              }
              if (n == 0L)
                  return(if (kind == "z") NaN * 0i else NaN)
              if (kind == "n") {
                  nnz <- length(x@i)
                  if (trim <= 0)
                      return(nnz / n)
                  ntrim <- trunc(n * min(trim, 0.5))
                  if (nnz < ntrim)
                      0
                  else if (nnz == ntrim) {
                      if (n - 2 * ntrim > 0)
                           0
                      else 0.5
                  } else {
                      if (n - 2 * ntrim > 0)
                          (nnz - ntrim - max(ntrim - (n - nnz), 0)) /
                              (n - 2 * ntrim)
                      else 1
                  }
              } else {
                  if (trim <= 0)
                      return(sum(x@x, na.rm = na.rm) / n)
                  ntrim <- trunc(n * min(trim, 0.5))
                  x <- .V.sort(x, na.last = NA)[(ntrim + 1):(n - ntrim)]
                  sum(x@x) / length(x)
              }
          })
