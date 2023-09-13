## METHODS FOR GENERIC: Summary (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Summary")
## [1] "max"   "min"   "range" "prod"  "sum"   "any"   "all"

## NB: Summary depends on the existence, _not_ count, of zeros and ones.
##     The only exception is 'sum' which ignores zeros and counts ones.

setMethod("Summary", signature(x = "denseMatrix"),
          function(x, ..., na.rm = FALSE) {
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              switch(kind,
                     "n" = ,
                     "l" = { zero <- FALSE; one <- TRUE },
                     "i" = { zero <- 0L   ; one <- 1L   },
                     "d" = { zero <- 0    ; one <- 1    },
                     "z" = { zero <- 0+0i ; one <- 1+0i })
              if(shape == "s" && any(.Generic == c("sum", "prod"))) {
                  x <- .M2gen(x)
                  shape <- "g"
              }
              if(shape != "g") {
                  if(repr != "p")
                      x <- .M2packed(x)
                  if(shape == "t" && x@diag != "N")
                      diag(x) <- TRUE
              }
              n <- (d <- x@Dim)[2L]
              y <- x@x
              y1 <- if(kind != "n" || !anyNA(y))
                        y
                    else y | is.na(y)
              y2 <- if(shape == "t" && n > 1L)
                        zero
              get(.Generic, mode = "function")(y1, y2, ..., na.rm = na.rm)
          })

setMethod("Summary", signature(x = "sparseMatrix"),
          function(x, ..., na.rm = FALSE) {
              cl <- .M.nonvirtual(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              switch(kind,
                     "n" = ,
                     "l" = { zero <- FALSE; one <- TRUE },
                     "i" = { zero <- 0L   ; one <- 1L   },
                     "d" = { zero <- 0    ; one <- 1    },
                     "z" = { zero <- 0+0i ; one <- 1+0i })
              if(shape == "s" && any(.Generic == c("sum", "prod"))) {
                  x <- .M2gen(x)
                  shape <- "g"
              }
              if(repr == "T") {
                  x <- aggregateT(x)
                  nnz <- length(x@i)
              } else {
                  nnz <- { p <- x@p; p[length(p)] }
                  if(length(if(repr == "C") x@i else x@j) > nnz) {
                      h <- seq_len(nnz)
                      if(repr == "C")
                          x@i <- x@i[h]
                      else
                          x@j <- x@j[h]
                      if(kind != "n")
                          x@x <- x@x[h]
                  }
              }
              n <- (d <- x@Dim)[2L]
              nnz.max <- if(shape == "s") 0.5 * (prod(d) + n) else prod(d)
              y1 <- if(kind != "n")
                        x@x
                    else if(.Generic == "sum")
                        nnz
                    else if(nnz > 0L)
                        TRUE
                    else logical(0L)
              y2 <- if(nnz < nnz.max)
                        zero
              y3 <- if(shape == "t" && x@diag != "N") {
                        if(.Generic == "sum")
                            one * n
                        else if(n > 0L)
                            one
                        else one[0L]
                    }
              get(.Generic, mode = "function")(y1, y2, y3, ..., na.rm = na.rm)
          })

setMethod("Summary", signature(x = "diagonalMatrix"),
          function(x, ..., na.rm = FALSE) {
              kind <- .M.kind(x)
              switch(kind,
                     "n" = ,
                     "l" = { zero <- FALSE; one <- TRUE },
                     "i" = { zero <- 0L   ; one <- 1L   },
                     "d" = { zero <- 0    ; one <- 1    },
                     "z" = { zero <- 0+0i ; one <- 1+0i })
              n <- x@Dim[2L]
              y <- x@x
              y1 <- if(kind != "n" || !anyNA(y))
                        y
                    else y | is.na(y)
              y2 <- if(n > 1L)
                        zero
              y3 <- if(x@diag != "N") {
                        if(.Generic == "sum")
                            one * n
                        else if(n > 0L)
                            one
                        else one[0L]
                    }
              get(.Generic, mode = "function")(y1, y2, y3, ..., na.rm = na.rm)
          })

setMethod("Summary", signature(x = "indMatrix"),
          function(x, ..., na.rm = FALSE) {
              nnz <- length(x@perm)
              y1 <- if(.Generic == "sum")
                        nnz
                    else if(nnz > 0L)
                        TRUE
                    else logical(0L)
              y2 <- if(nnz < prod(x@Dim))
                        FALSE
              get(.Generic, mode = "function")(y1, y2, ..., na.rm = na.rm)
          })

setMethod("Summary", signature(x = "sparseVector"),
          function(x, ..., na.rm = FALSE) {
              kind <- .M.kind(x)
              zero <- switch(kind, "n" = , "l" = FALSE, "i" = 0L, "d" = 0, "z" = 0+0i)
              nnz <- length(x@i)
              y1 <- if(kind != "n")
                        x@x
                    else if(.Generic == "sum")
                        nnz
                    else if(nnz > 0L)
                        TRUE
                    else logical(0L)
              y2 <- if(nnz < length(x))
                        zero
              get(.Generic, mode = "function")(y1, y2, ..., na.rm = na.rm)
          })
