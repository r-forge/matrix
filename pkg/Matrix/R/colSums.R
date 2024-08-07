## METHODS FOR GENERIC: colSums, rowSums, colMeans, rowMeans
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ==== denseMatrix ====================================================

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


## ==== sparseMatrix ===================================================

## ---- [CRT]sparseMatrix ----------------------------------------------

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
rm(.cl)


## ---- diagonalMatrix -------------------------------------------------

.diag.cS <- .diag.rS <- function(x, na.rm = FALSE, dims = 1L, ...) {
    kind <- .M.kind(x)
    if((n <- x@Dim[1L]) == 0L)
        return(vector(switch(kind, "z" = "complex", "d" = , "i" = "double", "integer"), 0L))
    else if(x@diag != "N")
        r <- rep.int(switch(kind, "z" = 1+0i, "d" = , "i" = 1, 1L), n)
    else {
        r <- switch(kind, "z" = , "d" = x@x, "i" = as.double(x@x), as.integer(x@x))
        if((na.rm || kind == "n") && anyNA(r))
            r[is.na(r)] <- switch(kind, "z" = 0+0i, "d" = , "i" = 0, "n" = 1L, 0L)
    }
    if(!is.null(nms <- x@Dimnames[[.MARGIN]]))
        names(r) <- nms
    r
}
body(.diag.cS) <- do.call(substitute, list(body(.diag.cS), list(.MARGIN = 2L)))
body(.diag.rS) <- do.call(substitute, list(body(.diag.rS), list(.MARGIN = 1L)))

.diag.cM <- .diag.rM <- function(x, na.rm = FALSE, dims = 1L, ...) {
    kind <- .M.kind(x)
    if((n <- x@Dim[1L]) == 0L)
        return(vector(switch(kind, "z" = "complex", "double"), 0L))
    else if(x@diag != "N")
        r <- rep.int(switch(kind, "z" = 1+0i, 1) / n, n)
    else {
        r <- x@x / n
        if((na.rm || kind == "n") && anyNA(r))
            r[is.na(r)] <- switch(kind,
                                  "z" = if(n == 1L) NaN * (0+0i) else 0+0i,
                                  "n" = 1 / n,
                                  if(n == 1L) NaN else 0)
    }
    if(!is.null(nms <- x@Dimnames[[.MARGIN]]))
        names(r) <- nms
    r
}
body(.diag.cM) <- do.call(substitute, list(body(.diag.cM), list(.MARGIN = 2L)))
body(.diag.rM) <- do.call(substitute, list(body(.diag.rM), list(.MARGIN = 1L)))

setMethod("colSums",  c(x = "diagonalMatrix"), .diag.cS)
setMethod("colMeans", c(x = "diagonalMatrix"), .diag.cM)
setMethod("rowSums",  c(x = "diagonalMatrix"), .diag.rS)
setMethod("rowMeans", c(x = "diagonalMatrix"), .diag.rM)

rm(.diag.cS, .diag.cM, .diag.rS, .diag.rM)


## ---- indMatrix (incl. pMatrix) --------------------------------------

setMethod("colSums",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              n <- x@Dim[2L]
              r <- if(x@margin == 1L)
                       tabulate(x@perm, n)
                   else rep.int(1L, n)
              if(!is.null(nms <- x@Dimnames[[2L]]))
                  names(r) <- nms
              r
          })
setMethod("colMeans",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              n <- (d <- x@Dim)[2L]
              r <- if(x@margin == 1L)
                       tabulate(x@perm, n) / d[1L]
                   else rep.int(1 / d[1L], n)
              if(!is.null(nms <- x@Dimnames[[2L]]))
                  names(r) <- nms
              r
          })
setMethod("rowSums",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              m <- x@Dim[1L]
              r <- if(x@margin == 1L)
                       rep.int(1L, m)
                   else tabulate(x@perm, m)
              if(!is.null(nms <- x@Dimnames[[1L]]))
                  names(r) <- nms
              r
          })
setMethod("rowMeans",  c(x = "indMatrix"),
          function(x, na.rm = FALSE, dims = 1L, ...) {
              m <- (d <- x@Dim)[1L]
              r <- if(x@margin == 1L)
                       rep.int(1 / d[2L], m)
                   else tabulate(x@perm, m) / d[2L]
              if(!is.null(nms <- x@Dimnames[[1L]]))
                  names(r) <- nms
              r
          })
