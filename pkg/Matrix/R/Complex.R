## METHODS FOR GENERIC: Complex (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Complex")
## [1] "Arg"  "Conj" "Im"   "Mod"  "Re"

setMethod("Complex", c(z = "denseMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(z, 6L)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (.Generic != "Conj" || kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              direct <- FALSE
              if (any(.Generic == c("Arg", "Im"))) {
                  if (shape == "s" && (kind != "z" || z@trans == "C")) {
                      if (kind != "z" && .Generic == "Im")
                          direct <- TRUE
                      else z <- .M2gen(z)
                      shape <- "g"
                      substr(cl, 2L, 3L) <- "ge"
                  }
                  if (shape == "t" && z@diag != "N")
                      diag(z) <- TRUE
              }
              r <- new(cl)
              r@Dim <- z@Dim
              r@Dimnames <- dimnames(z)
              if (shape != "g")
                  r@uplo <- z@uplo
              if (shape == "s" && kind == "z" && .Generic == "Conj")
                  r@trans <- z@trans
              if (shape == "t")
                  r@diag <- z@diag
              r@x <-
                  if (direct)
                      double(prod(z@Dim))
                  else g({ y <- z@x; if (kind == "n") y | is.na(y) else y })
              r
          })

setMethod("Complex", c(z = "sparseMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(z, 6L)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              if (repr == "T")
                  z <- agregateT(z)
              if (.Generic != "Conj" || kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              if (any(.Generic == c("Arg", "Im"))) {
                  if (shape == "s" && (kind != "z" || z@trans == "C")) {
                      z <- .M2gen(z)
                      shape <- "g"
                      substr(cl, 2L, 2L) <- "g"
                  }
                  if (shape == "t" && z@diag != "N")
                      diag(z) <- TRUE
              }
              r <- new(cl)
              r@Dim <- z@Dim
              r@Dimnames <- dimnames(z)
              if (shape != "g")
                  r@uplo <- z@uplo
              if (shape == "s" && kind == "z" && .Generic == "Conj")
                  r@trans <- z@trans
              if (shape == "t")
                  r@diag <- z@diag
              switch(repr,
                     "C" = { r@p <- z@p; i <- r@i <- z@i },
                     "R" = { r@p <- z@p; i <- r@j <- z@j },
                     "T" = { r@i <- z@i; i <- r@j <- z@j })
              r@x <-
                  if (kind == "n")
                      rep.int(g(TRUE), length(i))
                  else g(z@x)
              r
          })

setMethod("Complex", c(z = "diagonalMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              kind <- .M.kind(z)
              r <- new(if (.Generic != "Conj" || kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- d <- z@Dim
              r@Dimnames <- z@Dimnames
              if (z@diag == "N")
                  r@x <- g({ y <- z@x; if (kind == "n") y | is.na(y) else y })
              else if (any(.Generic == c("Arg", "Im")))
                  r@x <- double(d[1L])
              else r@diag <- "U"
              r
          })

setMethod("Complex", c(z = "indMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              g(.ind2sparse(z))
          })

setMethod("Complex", c(z = "sparseVector"),
          function(z) {
              g <- get(.Generic, mode = "function")
              kind <- .M.kind(z)
              r <- new(if (.Generic != "Conj" || kind != "z") "dsparseVector" else "zsparseVector")
              r@length <- z@length
              r@i <- i <- z@i
              r@x <-
                  if (kind == "n")
                      rep.int(g(TRUE), length(i))
                  else g(z@x)
              r
          })
