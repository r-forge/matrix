## METHODS FOR GENERIC: Complex (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Complex")
## [1] "Arg"  "Conj" "Im"   "Mod"  "Re"

.nm <- getGroupMembers("Complex")
.fn <- mget(.nm, mode = "function", inherits = TRUE)
.cc <- c("Conj", "Mod", "Re")
suppressWarnings({
## is.integer(g(integer())
.Complex.integer <-
vapply(.fn, function(g) is.integer(g(integer())), FALSE)
## is.complex(g(complex())
.Complex.complex <-
vapply(.fn, function(g) tryCatch(is.complex(g(complex())), error = function(e) TRUE), FALSE)
## g(0) == 0
.Complex.identity0 <-
vapply(.fn, function(g) identical(g(0), 0), FALSE)
## g(1) == 1
.Complex.identity1 <-
vapply(.fn, function(g) identical(g(1), 1), FALSE)
## g(Conj(z)) == Conj(g(z))
.Complex.identityH <-
`names<-`(match(.nm, .cc, 0L) > 0L, .nm)
})
rm(.fn, .nm, .cc)

# Assumed below:
stopifnot(!any(.Complex.integer), all(.Complex.identity0))

setMethod("Complex", c(z = "denseMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(z)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (!(kind == "z" && .Complex.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- "d"
              if (shape == "s" && kind == "z" && z@trans == "C" &&
                  !.Complex.identityH[[.Generic]]) {
                  z <- .M2gen(z)
                  shape <- "g"
                  substr(cl, 2L, 3L) <- "ge"
              }
              if (shape == "t" && z@diag != "N" &&
                  !.Complex.identity1[[.Generic]])
                  diag(z) <- TRUE
              r <- new(cl)
              r@Dim <- z@Dim
              r@Dimnames <- z@Dimnames
              if (shape != "g")
                  r@uplo <- z@uplo
              if (shape == "s" && kind == "z" &&
                  .Complex.complex[[.Generic]])
                  r@trans <- z@trans
              if (shape == "t")
                  r@diag <- z@diag
              r@x <- g({ y <- z@x; if (kind == "n") y | is.na(y) else y })
              r
          })

setMethod("Complex", c(z = "sparseMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(z)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              if (!(kind == "z" && .Complex.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- "d"
              if (shape == "s" && kind == "z" && z@trans == "C" &&
                  !.Complex.identityH[[.Generic]]) {
                  z <- .M2gen(z)
                  shape <- "g"
                  substr(cl, 2L, 2L) <- "g"
              }
              if (shape == "t" && z@diag != "N" &&
                  !.Complex.identity1[[.Generic]])
                  diag(z) <- TRUE
              r <- new(cl)
              r@Dim <- z@Dim
              r@Dimnames <- z@Dimnames
              if (shape != "g")
                  r@uplo <- z@uplo
              if (shape == "s" && kind == "z" &&
                  .Complex.complex[[.Generic]])
                  r@trans <- z@trans
              if (shape == "t")
                  r@diag <- z@diag
              nnz <- length(
                  switch(repr,
                         "C" = { r@p <- z@p; r@i <- z@i },
                         "R" = { r@p <- z@p; r@j <- z@j },
                         "T" = { z <- aggregateT(z)
                                 r@i <- z@i; r@j <- z@j }))
              r@x <- if (kind == "n") rep.int(g(TRUE), nnz) else g(z@x)
              r
          })

setMethod("Complex", c(z = "diagonalMatrix"),
          function(z) {
              g <- get(.Generic, mode = "function")
              cl <- .M.class(z)
              kind <- substr(cl, 1L, 1L)
              if (!(kind == "z" && .Complex.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- "d"
              r <- new(cl)
              r@Dim <- d <- z@Dim
              r@Dimnames <- z@Dimnames
              if (z@diag == "N")
                  r@x <- g({ y <- z@x; if (kind == "n") y | is.na(y) else y })
              else if (.Complex.identity1[[.Generic]])
                  r@diag <- "U"
              else r@x <- rep.int(g(switch(kind, "z" = 1+0i, "d" = 1, 1L)), d[1L])
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
              cl <- .M.class(z)
              kind <- substr(cl, 1L, 1L)
              if (!(kind == "z" && .Complex.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- "d"
              r <- new(cl)
              r@length <- z@length
              r@i <- i <- z@i
              r@x <- if (kind == "n") rep.int(g(TRUE), length(i)) else g(z@x)
              r
          })
