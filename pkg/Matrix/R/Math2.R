## METHODS FOR GENERIC: Math2 (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Math2")
## [1] "round"  "signif"

setMethod("Math2", c(x = "denseMatrix"),
          function(x, digits) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              if (shape == "t" && x@diag != "N" &&
                  .Generic == "round" && digits <= -1L)
                  diag(x) <- TRUE
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g")
                  r@uplo <- x@uplo
              if (shape == "s" && kind == "z")
                  r@trans <- x@trans
              if (shape == "t")
                  r@diag <- x@diag
              r@x <- g({ y <- x@x; if (kind == "n") y | is.na(y) else y }, digits = digits)
              r
          })

setMethod("Math2", c(x = "sparseMatrix"),
          function(x, digits) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              if (kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              if (shape == "t" && x@diag != "N" &&
                  .Generic == "round" && digits <= -1L)
                  diag(x) <- TRUE
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g")
                  r@uplo <- x@uplo
              if (shape == "s" && kind == "z")
                  r@trans <- x@trans
              if (shape == "t")
                  r@diag <- x@diag
              nnz <- length(
                  switch(repr,
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { x <- aggregateT(x)
                                 r@i <- x@i; r@j <- x@j }))
              r@x <- if (kind == "n") rep.int(g(TRUE, digits = digits), nnz) else g(x@x, digits = digits)
              r
          })

setMethod("Math2", c(x = "diagonalMatrix"),
          function(x, digits) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              if (kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (x@diag == "N")
                  r@x <- g({ y <- x@x; if (kind == "n") y | is.na(y) else y }, digits = digits)
              else if (.Generic != "round" || digits > -1L)
                  r@diag <- "U"
              else r@x <- rep.int(g(switch(kind, "z" = 1+0i, "d" = 1, 1L), digits = digits), d[1L])
              r
          })

setMethod("Math2", c(x = "indMatrix"),
          function(x, digits) {
              g <- get(.Generic, mode = "function")
              g(.ind2sparse(x), digits = digits)
          })

setMethod("Math2", c(x = "sparseVector"),
          function(x, digits) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              r <- new(if (kind == "z") "zsparseVector" else "dsparseVector")
              r@length <- x@length
              r@i <- i <- x@i
              r@x <- if (kind == "n") rep.int(g(TRUE, digits = digits), length(i)) else g(x@x, digits = digits)
              r
          })
