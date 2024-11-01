## METHODS FOR GENERIC: rcond
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("rcond", c(x = "ANY", norm = "missing"),
          function(x, norm, ...)
              rcond(x, norm = "O", ...))

setMethod("rcond", c(x = "denseMatrix", norm = "character"),
          function(x, norm, ...) {
              x <- .M2kind(x, ",")
              switch(substr(.M.class(x, 2L), 2L, 3L),
                     "ge" =
                         {
                             d <- x@Dim
                             m <- d[1L]
                             n <- d[2L]
                             if (m == n) {
                                 trf <- lu(x, warnSing = FALSE)
                                 .Call(geMatrix_rcond, x, trf, norm)
                             } else {
                                 ## MJ: norm(A = P1' Q R P2') = norm(R) holds
                                 ##     in general only for norm == "2", but
                                 ##     La_rcond_type() disallows norm == "2"
                                 ##     ... FIXME ??
                                 if (m < n) {
                                     x <- t(x)
                                     n <- m
                                 }
                                 R <- triu(qr(x)[["qr"]][seq_len(n), , drop = FALSE])
                                 rcond(R, norm = norm, ...)
                             }
                         },
                     "sy" =
                         {
                             trf <- BunchKaufman(x, warnSing = FALSE)
                             .Call(syMatrix_rcond, x, trf, norm)
                         },
                     "sp" =
                         {
                             trf <- BunchKaufman(x, warnSing = FALSE)
                             .Call(spMatrix_rcond, x, trf, norm)
                         },
                     "po" =
                         {
                             ok <- TRUE
                             trf <- tryCatch(
                                 Cholesky(x, perm = FALSE),
                                 error = function(e) {
                                     ok <<- FALSE
                                     BunchKaufman(x, warnSing = FALSE)
                                 })
                             if (ok)
                                 .Call(poMatrix_rcond, x, trf, norm)
                             else .Call(syMatrix_rcond, x, trf, norm)
                         },
                     "pp" =
                         {
                             ok <- TRUE
                             trf <- tryCatch(
                                 Cholesky(x, perm = FALSE),
                                 error = function(e) {
                                     ok <<- FALSE
                                     BunchKaufman(x, warnSing = FALSE)
                                 })
                             if (ok)
                                 .Call(ppMatrix_rcond, x, trf, norm)
                             else .Call(spMatrix_rcond, x, trf, norm)
                         },
                     "tr" = .Call(trMatrix_rcond, x, norm),
                     "tp" = .Call(tpMatrix_rcond, x, norm))
          })

setMethod("rcond", c(x = "sparseMatrix", norm = "character"),
          function(x, norm, useInv = FALSE, ...) {
              d <- x@Dim
              if ((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(Inf)
              if (m == n) {
                  if (isS4(useInv) || useInv) {
                      if (!isS4(useInv))
                          useInv <- solve(x)
                      1/(norm(x, type = norm) * norm(useInv, type = norm))
                  } else {
                      warning(gettextf("'%s' via sparse -> dense coercion",
                                       "rcond"),
                              domain = NA)
                      rcond(.M2unpacked(x), norm = norm, ...)
                  }
              } else {
                  ## MJ: norm(A = P1' Q R P2') = norm(R) holds in general
                  ##     only for norm == "2", but La_rcond_type() disallows
                  ##     norm == "2" ... FIXME ??
                  if (m < n) {
                      x <- t(x)
                      n <- m
                  }
                  R <- triu(qr(x)@R[seq_len(n), , drop = FALSE])
                  rcond(R, norm = norm, ...)
              }
           })

setMethod("rcond", c(x = "diagonalMatrix", norm = "character"),
          function(x, norm, ...) {
              if ((n <- x@Dim[1L]) == 0L)
                  return(Inf)
              if (nonunit <- x@diag == "N") {
                  y <- x@x
                  if (.M.kind(x) == "n" && anyNA(y))
                      y <- y | is.na(y)
              }
              switch(EXPR = norm[1L],
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =,
                     "2" =,
                     "M" =, "m" =
                         if (nonunit) {
                             ry <- range(abs(y))
                             ry[1L]/ry[2L]
                         } else 1,
                     "F" =, "f" =, "E" =, "e" =
                         if (nonunit) {
                             if (is.complex(y))
                                 y <- abs(y)
                             yy <- y * y
                             1/sqrt(sum(yy) * sum(1/yy))
                         } else 1/n,
                     stop(gettext("invalid %s=\"%s\"",
                                  "norm", norm[1L]),
                          domain = NA))
          })

setMethod("rcond", c(x = "indMatrix", norm = "character"),
          function(x, norm, ...) {
              d <- x@Dim
              if ((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(Inf)
              if (m == n) {
                  if (anyDuplicated.default(x@perm))
                      return(0)
                  switch(EXPR = norm[1L],
                         "O" =, "o" =, "1" =,
                         "I" =, "i" =,
                         "2" =,
                         "M" =, "m" =
                             1,
                         "F" =, "f" =, "E" =, "e" =
                             1/n,
                         stop(gettext("invalid %s=\"%s\"",
                                      "norm", norm[1L]),
                              domain = NA))
              } else {
                  if (m < n) {
                      x <- t(x)
                      n <- m
                  }
                  R <- triu(qr(x)@R[seq_len(n), , drop = FALSE])
                  rcond(R, norm = norm, ...)
              }
          })

setMethod("rcond", c(x = "pMatrix", norm = "character"),
          function(x, norm, ...) {
              if ((n <- x@Dim[1L]) == 0L)
                  return(Inf)
              switch(EXPR = norm[1L],
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =,
                     "2" =,
                     "M" =, "m" =
                         1,
                     "F" =, "f" =, "E" =, "e" =
                         1/n,
                     stop(gettext("invalid %s=\"%s\"",
                                  "norm", norm[1L]),
                          domain = NA))
          })
