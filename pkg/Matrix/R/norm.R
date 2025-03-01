## METHODS FOR GENERIC: norm
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("norm", c(x = "ANY", type = "missing"),
          function(x, type, ...)
              norm(x, type = "O", ...))

setMethod("norm", c(x = "denseMatrix", type = "character"),
          function(x, type, ...) {
              if (identical(type, "2"))
                  return(base::norm(.M2m(x), type = "2"))
              x <- .M2kind(x, ",")
              switch(substr(.M.class(x), 2L, 3L),
                     "ge" = .Call(geMatrix_norm, x, type),
                     "sy" = .Call(syMatrix_norm, x, type),
                     "sp" = .Call(spMatrix_norm, x, type),
                     "tr" = .Call(trMatrix_norm, x, type),
                     "tp" = .Call(tpMatrix_norm, x, type))
          })

setMethod("norm", c(x = "sparseMatrix", type = "character"),
          function(x, type, ...) {
              if (any(x@Dim == 0L))
                  return(0)
              switch(type,
                     "2" =
                         base::norm(.M2m(x), type = "2"),
                     "M" =, "m" =
                         max(range(abs(x))),
                     "O" =, "o" =, "1" =
                         max(colSums(abs(x))),
                     "I" =, "i" =
                         max(rowSums(abs(x))),
                     "F" =, "f" =, "E" =, "e" =
                         {
                             if (.M.kind(x) == "z")
                                 x <- abs(x)
                             sqrt(sum(x * x))
                         },
                     stop(gettextf("invalid %s=\"%s\"",
                                   "type", type),
                          domain = NA))
          })

setMethod("norm", c(x = "diagonalMatrix", type = "character"),
          function(x, type, ...) {
              if ((n <- x@Dim[1L]) == 0L)
                  return(0)
              if (nonunit <- x@diag == "N") {
                  y <- x@x
                  if (.M.kind(x) == "n" && anyNA(y))
                      y <- y | is.na(y)
              }
              switch(type,
                     "2" =,
                     "M" =, "m" =,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =
                         if (nonunit) max(abs(y)) else 1,
                     "F" =, "f" =, "E" =, "e" =
                         if (nonunit) {
                             if (is.complex(y))
                                 y <- abs(y)
                             sqrt(sum(y * y))
                         } else sqrt(n),
                     stop(gettextf("invalid %s=\"%s\"",
                                   "type", type),
                          domain = NA))
          })

setMethod("norm", c(x = "indMatrix", type = "character"),
          function(x, type, ...) {
              d <- x@Dim
              if ((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(0)
              switch(type,
                     "2" =
                         sqrt(max(tabulate(x@perm, if (x@margin == 1L) n else m))),
                     "M" =, "m" =
                         1,
                     "O" =, "o" =, "1" =
                         if (x@margin == 1L) max(tabulate(x@perm, n)) else 1,
                     "I" =, "i" =
                         if (x@margin == 1L) 1 else max(tabulate(x@perm, m)),
                     "F" =, "f" =, "E" =, "e" =
                         if (x@margin == 1L) sqrt(m) else sqrt(n),
                     stop(gettextf("invalid %s=\"%s\"",
                                   "type", type),
                          domain = NA))
          })

setMethod("norm", c(x = "pMatrix", type = "character"),
          function(x, type, ...) {
              if ((n <- x@Dim[1L]) == 0L)
                  return(0)
              switch(type,
                     "2" =,
                     "M" =, "m" =,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =
                         1,
                     "F" =, "f" =, "E" =, "e" =
                         sqrt(n),
                     stop(gettextf("invalid %s=\"%s\"",
                                   "type", type),
                          domain = NA))
          })
