## METHODS FOR GENERIC: colSums, colMeans, rowSums, rowMeans
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("colSums",  c(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_marginsum, x, 2L, na.rm, FALSE))
setMethod("colMeans", c(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_marginsum, x, 2L, na.rm,  TRUE))
setMethod("rowSums",  c(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_marginsum, x, 1L, na.rm, FALSE))
setMethod("rowMeans", c(x = "denseMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...)
              .Call(R_dense_marginsum, x, 1L, na.rm,  TRUE))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix")) {
setMethod("colSums",  c(x = .cl),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(R_sparse_marginsum, x, 2L, na.rm, FALSE, sparseResult))
setMethod("colMeans", c(x = .cl),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(R_sparse_marginsum, x, 2L, na.rm,  TRUE, sparseResult))
setMethod("rowSums",  c(x = .cl),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(R_sparse_marginsum, x, 1L, na.rm, FALSE, sparseResult))
setMethod("rowMeans", c(x = .cl),
          function(x, na.rm = FALSE, dims = 1L, sparseResult = FALSE, ...)
              .Call(R_sparse_marginsum, x, 1L, na.rm,  TRUE, sparseResult))
}

.marSums.diagonal <-
function(x, na.rm = FALSE, dims = 1L, ...) {
    kind <- .M.kind(x)
    n <- x@Dim[1L]
    if (n == 0L)
        return(vector(switch(kind, "z" = "complex", "d" = , "i" = "double", "integer"), 0L))
    if (x@diag != "N")
        r <- rep.int(switch(kind, "z" = 1+0i, "d" = , "i" = 1, 1L), n)
    else {
        r <- switch(kind, "z" = , "d" = x@x, "i" = as.double(x@x), as.integer(x@x))
        if ((na.rm || kind == "n") && anyNA(r))
            r[is.na(r)] <- switch(kind, "z" = 0+0i, "d" = , "i" = 0, "n" = 1L, 0L)
    }
    if (!is.null(nms <- x@Dimnames[[.MARGIN]]))
        names(r) <- nms
    r
}
.marMeans.diagonal <-
function(x, na.rm = FALSE, dims = 1L, ...) {
    kind <- .M.kind(x)
    n <- x@Dim[1L]
    if (n == 0L)
        return(vector(switch(kind, "z" = "complex", "double"), 0L))
    else if (x@diag != "N")
        r <- rep.int(switch(kind, "z" = 1+0i, 1) / n, n)
    else {
        r <- x@x / n
        if ((na.rm || kind == "n") && anyNA(r))
            r[is.na(r)] <-
                switch(kind,
                       "z" = if (n == 1L) NaN * (0+0i) else 0+0i,
                       "n" = 1 / n,
                       if (n == 1L) NaN else 0)
    }
    if (!is.null(nms <- x@Dimnames[[.MARGIN]]))
        names(r) <- nms
    r
}

setMethod("colSums",  c(x = "diagonalMatrix"),
          updateBody(.marSums.diagonal, .MARGIN = 2L))
setMethod("colMeans", c(x = "diagonalMatrix"),
          updateBody(.marMeans.diagonal, .MARGIN = 2L))
setMethod("rowSums",  c(x = "diagonalMatrix"),
          updateBody(.marSums.diagonal, .MARGIN = 1L))
setMethod("rowMeans", c(x = "diagonalMatrix"),
          updateBody(.marMeans.diagonal, .MARGIN = 1L))

rm(.marSums.diagonal, .marMeans.diagonal)

setMethod("colSums",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              n <- x@Dim[2L]
              r <- if (x@margin == 1L)
                       tabulate(x@perm, n)
                   else rep.int(1L, n)
              if (!is.null(nms <- x@Dimnames[[2L]]))
                  names(r) <- nms
              r
          })
setMethod("colMeans",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              n <- (d <- x@Dim)[2L]
              r <- if (x@margin == 1L)
                       tabulate(x@perm, n) / d[1L]
                   else rep.int(1 / d[1L], n)
              if (!is.null(nms <- x@Dimnames[[2L]]))
                  names(r) <- nms
              r
          })
setMethod("rowSums",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              m <- x@Dim[1L]
              r <- if (x@margin == 1L)
                       rep.int(1L, m)
                   else tabulate(x@perm, m)
              if (!is.null(nms <- x@Dimnames[[1L]]))
                  names(r) <- nms
              r
          })
setMethod("rowMeans",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              m <- (d <- x@Dim)[1L]
              r <- if (x@margin == 1L)
                       rep.int(1 / d[2L], m)
                   else tabulate(x@perm, m) / d[2L]
              if (!is.null(nms <- x@Dimnames[[1L]]))
                  names(r) <- nms
              r
          })

rm(.cl)
