## METHODS FOR GENERIC: Math (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Math")
##  [1] "abs"      "sign"     "sqrt"     "ceiling"  "floor"    "trunc"
##  [7] "cummax"   "cummin"   "cumprod"  "cumsum"   "exp"      "expm1"
## [13] "log"      "log10"    "log2"     "log1p"    "cos"      "cosh"
## [19] "sin"      "sinh"     "tan"      "tanh"     "acos"     "acosh"
## [25] "asin"     "asinh"    "atan"     "atanh"    "cospi"    "sinpi"
## [31] "tanpi"    "gamma"    "lgamma"   "digamma"  "trigamma"

.nm <- getGroupMembers("Math")
.fn <- mget(.nm, mode = "function", inherits = TRUE)
.cc <- c("abs", "sqrt", "exp", "expm1", "log", "log10", "log2", "log1p",
         "cos", "cosh", "sin", "sinh", "tan", "tanh",
         "acos", "acosh", "asin", "asinh", "atan", "atanh",
         "cospi", "sinpi", "tanpi",
         "gamma", "lgamma", "digamma", "trigamma")
suppressWarnings({
## is.integer(g(integer())
.Math.integer <-
vapply(.fn, function(g) is.integer(g(integer())), FALSE)
## is.complex(g(complex())
.Math.complex <-
vapply(.fn, function(g) tryCatch(is.complex(g(complex())), error = function(e) TRUE), FALSE)
## g(0) == 0
.Math.identity0 <-
vapply(.fn, function(g) identical(g(0), 0), FALSE)
## g(1) == 1
.Math.identity1 <-
vapply(.fn, function(g) identical(g(1), 1), FALSE)
## g(Conj(x)) == Conj(g(x))
.Math.identityH <-
`names<-`(match(.nm, .cc, 0L) > 0L, .nm)
})
rm(.fn, .nm, .cc)

setMethod("Math", c(x = "denseMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if (startsWith(.Generic, "cum"))
                  return(g(.M2v(x)))
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (!(kind == "z" && .Math.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- if (kind == "z" || kind == "d" || !.Math.integer[[.Generic]]) "d" else "i"
              if ((shape == "s" && kind == "z" && x@trans == "C" &&
                   !.Math.identityH[[.Generic]]) ||
                  (shape == "t" && !.Math.identity0[[.Generic]])) {
                  x <- .M2gen(x)
                  shape <- "g"
                  substr(cl, 2L, 3L) <- "ge"
              }
              if (shape == "t" && x@diag != "N" &&
                  !.Math.identity1[[.Generic]])
                  diag(x) <- TRUE
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g")
                  r@uplo <- x@uplo
              if (shape == "s" && kind == "z" &&
                  .Math.complex[[.Generic]])
                  r@trans <- x@trans
              if (shape == "t")
                  r@diag <- x@diag
              r@x <- g({ y <- x@x; if (kind == "n") y | is.na(y) else y })
              r
          })

setMethod("Math", c(x = "sparseMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if (startsWith(.Generic, "cum"))
                  return(g(.M2v(x)))
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              repr <- substr(cl, 3L, 3L)
              if (!(kind == "z" && .Math.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- if (kind == "z" || kind == "d" || !.Math.integer[[.Generic]]) "d" else "i"
              if ((shape == "s" && kind == "z" && x@trans == "C" &&
                   !.Math.identityH[[.Generic]]) ||
                  (shape == "t" && !.Math.identity0[[.Generic]])) {
                  x <- .M2gen(x)
                  shape <- "g"
                  substr(cl, 2L, 2L) <- "g"
              }
              if (!.Math.identity0[[.Generic]])
                  substr(cl, 3L, 3L) <- switch(shape, "g" = "e", "s" = "y", "t" = "r")
              if (shape == "t" && x@diag != "N" &&
                  !.Math.identity1[[.Generic]])
                  diag(x) <- TRUE
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g")
                  r@uplo <- x@uplo
              if (shape == "s" && kind == "z" &&
                  .Math.complex[[.Generic]])
                  r@trans <- x@trans
              if (shape == "t")
                  r@diag <- x@diag
              if (!.Math.identity0[[.Generic]]) {
              x <- .M2V(x)
              r@x <- replace(rep.int(g(switch(kind, "z" = 0+0i, "d" = 0, 0L)), x@length),
                             x@i,
                             g(if (kind == "n") TRUE else x@x))
              } else {
              nnz <- length(
                  switch(repr,
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { x <- aggregateT(x)
                                 r@i <- x@i; r@j <- x@j }))
              r@x <- if (kind == "n") rep.int(g(TRUE), nnz) else g(x@x)
              }
              r
          })

setMethod("Math", c(x = "diagonalMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if (startsWith(.Generic, "cum"))
                  return(g(.M2v(x)))
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              if (!(kind == "z" && .Math.complex[[.Generic]]))
                  substr(cl, 1L, 1L) <- if (kind == "z" || kind == "d" || !.Math.integer[[.Generic]]) "d" else "i"
              if (!.Math.identity0[[.Generic]])
                  substr(cl, 2L, 3L) <- "ge"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (!.Math.identity0[[.Generic]])
                  r@x <- replace(rep.int(g(switch(kind, "z" = 0+0i, "d" = 0, 0L)), prod(d)),
                                 seq.int(from = 1L, by = d[1L] + 1, length.out = d[1L]),
                                 g(if (x@diag == "N") { y <- x@x; if (kind == "n") y | is.na(y) else y } else switch(kind, "z" = 1+0i, "d" = 1, 1L)))
              else if (x@diag == "N")
                  r@x <- g({ y <- x@x; if (kind == "n") y | is.na(y) else y })
              else if (.Math.identity1[[.Generic]])
                  r@diag <- "U"
              else r@x <- rep.int(g(switch(kind, "z" = 1+0i, "d" = 1, 1L)), d[1L])
              r
          })

setMethod("Math", c(x = "indMatrix"),
          function(x) {
              g <- get(.Generic, mode = "function")
              g(.ind2sparse(x))
          })

setMethod("Math", c(x = "sparseVector"),
          function(x) {
              g <- get(.Generic, mode = "function")
              if (startsWith(.Generic, "cum"))
                  return(g(.V2v(x)))
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              if (!.Math.identity0[[.Generic]])
                  replace(rep.int(g(switch(kind, "z" = 0+0i, "d" = 0, 0L)), x@length),
                          x@i,
                          g(if (kind == "n") TRUE else x@x))
              else {
                  if (!(kind == "z" && .Math.complex[[.Generic]]))
                      substr(cl, 1L, 1L) <- if (kind == "z" || kind == "d" || !.Math.integer[[.Generic]]) "d" else "i"
                  r <- new(cl)
                  r@length <- x@length
                  r@i <- i <- x@i
                  r@x <- if (kind == "n") rep.int(g(TRUE), length(i)) else g(x@x)
                  r
              }
          })


## METHODS FOR GENERIC: log
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("log", c(x = "denseMatrix"),
          function(x, ...) {
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              if (shape == "t") {
                  x <- .M2gen(x)
                  shape <- "g"
                  substr(cl, 2L, 3L) <- "ge"
              }
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g")
                  r@uplo <- x@uplo
              if (shape == "s" && kind == "z")
                  r@trans <- x@trans
              r@x <- log({ y <- x@x; if (kind == "n") y | is.na(y) else y }, ...)
              r
          })

setMethod("log", c(x = "sparseMatrix"),
          function(x, ...) {
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              if (kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              if (shape == "t") {
                  shape <- "g"
                  substr(cl, 2L, 2L) <- "g"
              }
              substr(cl, 3L, 3L) <- if (shape == "g") "e" else "y"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g")
                  r@uplo <- x@uplo
              if (shape == "s" && kind == "z")
                  r@trans <- x@trans
              x <- .M2V(x)
              r@x <- replace(rep.int(log(switch(kind, "z" = 0+0i, "d" = 0, 0L), ...), x@length),
                             x@i,
                             log(if (kind == "n") TRUE else x@x, ...))
              r
          })

setMethod("log", c(x = "diagonalMatrix"),
          function(x, ...) {
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              if (kind != "z")
                  substr(cl, 1L, 1L) <- "d"
              substr(cl, 2L, 3L) <- "ge"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- replace(rep.int(log(switch(kind, "z" = 0+0i, "d" = 0, 0L), ...), prod(d)),
                             seq.int(from = 1L, by = d[1L] + 1, length.out = d[1L]),
                             log(if (x@diag == "N") { y <- x@x; if (kind == "n") y | is.na(y) else y } else switch(kind, "z" = 1+0i, "d" = 1, 1L), ...))
              r
          })

setMethod("log", c(x = "indMatrix"),
          function(x, ...)
              log(.ind2sparse(x), ...))

setMethod("log", c(x = "sparseVector"),
          function(x, ...) {
              cl <- .M.class(x)
              kind <- substr(cl, 1L, 1L)
              replace(rep.int(log(switch(kind, "z" = 0+0i, "d" = 0, 0L), ...), x@length),
                      x@i,
                      log(if (kind == "n") TRUE else x@x, ...))
          })
