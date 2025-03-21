## METHODS FOR GENERIC: Summary (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Summary")
## [1] "max"   "min"   "range" "prod"  "sum"   "any"   "all"

## NB: Summary depends on the existence, _not_ count, of zeros and ones.
##     The only exception is 'sum' which ignores zeros and counts ones.

setMethod("Summary", c(x = "denseMatrix"),
          function(x, ..., na.rm = FALSE) {
              ## Avoid wrong overflow :
              if (.Generic == "sum")
                  return(sum (.Call(R_dense_sum , x, na.rm),
                              ..., na.rm = na.rm))
              if (.Generic == "prod")
                  return(prod(.Call(R_dense_prod, x, na.rm),
                              ..., na.rm = na.rm))
              g <- get(.Generic, mode = "function")
              x <- forceCanonical(x)
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              zero <- switch(kind, "z" = 0+0i, "d" = 0, "i" = 0L, FALSE)
              n <- x@Dim[2L]
              y1 <- x@x
              y2 <- if (shape == "t" && n > 1L)
                        zero
              g(y1, y2, ..., na.rm = na.rm)
          })

setMethod("Summary", c(x = "sparseMatrix"), # (*not* diagonalM*, nor indMatrix (done below))
          function(x, ..., na.rm = FALSE) {
              ## Avoid wrong overflow :
              if (.Generic == "sum")
                  return(sum (.Call(R_sparse_sum , x, na.rm),
                              ..., na.rm = na.rm))
              if (.Generic == "prod")
                  return(prod(.Call(R_sparse_prod, x, na.rm),
                              ..., na.rm = na.rm))
              g <- get(.Generic, mode = "function")
              x <- forceCanonical(x)
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              zero <- switch(kind, "z" = 0+0i, "d" = 0, "i" = 0L, FALSE)
              unit <- switch(kind, "z" = 1+0i, "d" = 1, "i" = 1L,  TRUE)
              n <- (d <- x@Dim)[2L]
              nnz <- length(switch(repr, "C" = x@i, "R"=, "T" = x@j, stop("invalid 'repr': ", repr)))
              if(shape == "t" && x@diag == "U") nnz <- nnz + n
              nnz.max <- if (shape == "s") 0.5 * (prod(d) + n) else prod(d)
              y1 <- if (kind != "n")
                        x@x
                    else if (nnz > 0L)
                        TRUE
                    else logical(0L)
              y2 <- if (nnz < nnz.max)
                        zero
              y3 <- if (shape == "t" && n > 0L && x@diag != "N")
                        unit
              g(y1, y2, y3, ..., na.rm = na.rm)
          })

setMethod("Summary", c(x = "diagonalMatrix"),
          function(x, ..., na.rm = FALSE) {
              g <- get(.Generic, mode = "function")
              kind <- .M.kind(x)
              zero <- switch(kind, "z" = 0+0i, "d" = 0, "i" = 0L, FALSE)
              unit <- switch(kind, "z" = 1+0i, "d" = 1, "i" = 1L,  TRUE)
              n <- x@Dim[2L]
              y1 <- if (x@diag == "N") {
                        y <- x@x
                        if (kind != "n") {
                            if (.Generic == "prod" && n > 1L)
                                ## Avoid wrong overflow :
                                c(y[1L], zero, y[-1L])
                            else y
                        }
                        else y | is.na(y)
                    }
              y2 <- if (n > 1L)
                        zero
              y3 <- if (x@diag != "N") {
                        if (.Generic == "sum")
                            unit * n
                        else if (n > 0L)
                            unit
                        else unit[0L]
                    }
              g(y1, y2, y3, ..., na.rm = na.rm)
          })

setMethod("Summary", c(x = "indMatrix"),
          function(x, ..., na.rm = FALSE) {
              g <- get(.Generic, mode = "function")
              nnz <- length(x@perm)
              nnz.max <- prod(x@Dim)
              y1 <- if (.Generic == "sum")
                        nnz
                    else if (nnz > 0L)
                        TRUE
                    else logical(0L)
              y2 <- if (nnz < nnz.max)
                        FALSE
              g(y1, y2, ..., na.rm = na.rm)
          })

setMethod("Summary", c(x = "sparseVector"),
          function(x, ..., na.rm = FALSE) {
              g <- get(.Generic, mode = "function")
              kind <- .M.kind(x)
              zero <- switch(kind, "z" = 0+0i, "d" = 0, "i" = 0L, FALSE)
              nnz <- length(i <- x@i)
              nnz.max <- length(x)
              y1 <- if (kind != "n") {
                        y <- x@x
                        if (.Generic == "prod" && nnz > 0L && nnz < nnz.max) {
                            ## Avoid wrong overflow :
                            if (i[1L] > 1L)
                                c(zero, y)
                            else if (nnz >= (q <- which.min(i == seq_along(i))))
                                c(y[1L:(q - 1L)], zero, y[q:nnz])
                            else y
                        }
                        else y
                    }
                    else if (.Generic == "sum")
                        nnz
                    else if (nnz > 0L)
                        TRUE
                    else logical(0L)
              y2 <- if (nnz < nnz.max)
                        zero
              g(y1, y2, ..., na.rm = na.rm)
          })
