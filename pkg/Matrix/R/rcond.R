## METHODS FOR GENERIC: kappa
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("kappa", c(z = "Matrix"),
          function(z, ...) 1/rcond(z, ...))


## METHODS FOR GENERIC: rcond
## NB: approximation by QR factorization is exact for
##     rcond(<square>, norm = <"2" or "F">)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("rcond", c(x = "ANY", norm = "missing"),
          function(x, norm, ...)
              rcond(x, norm = "O", ...))

setMethod("rcond", c(x = "denseMatrix", norm = "character"),
          function(x, norm, exact = FALSE, inverse = solve(x),
                   warn = TRUE, ...) {
              switch(norm,
                     "2" =,
                     "M" =, "m" =,
                     "F" =, "f" =, "E" =, "e" =
                         exact <- TRUE,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =
                         NULL,
                     stop(gettext("invalid %s=\"%s\"",
                                  "norm", norm),
                          domain = NA))
              d <- x@Dim
              if ((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(Inf)
              x <- .M2kind(x, ",")
              if (exact)
              switch(norm,
                     "2" =,
                         {
                             s <- svd(x, nu = 0L, nv = 0L)[["d"]]
                             if (s[1L] > 0) s[length(s)]/s[1L] else 0
                         },
                     "M" =, "m" =,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =,
                     "F" =, "f" =, "E" =, "e" =
                         1/(norm(x, type = norm) * norm(inverse, type = norm)))
              else
              switch(substr(.M.class(x, 2L), 2L, 3L),
                     "ge" =
                         {
                             if (m == n) {
                             trf <- lu(x, warnSing = FALSE)
                             .Call(geMatrix_rcond, x, trf, norm)
                             } else {
                             if (warn)
                             warning(gettextf("rcond(%s, %s) via QR factorization; dubious for non-square '%s', %s=\"%s\"",
                                              "x", "norm", "x", "norm", norm),
                                     domain = NA)
                             R <- triu(qr(if (m < n) ct(x) else x)[["qr"]][seq_len(if (m < n) m else n), , drop = FALSE])
                             .Call(trMatrix_rcond, R, norm)
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
                             trf <- tryCatch(Cholesky(x, perm = FALSE),
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
                             trf <- tryCatch(Cholesky(x, perm = FALSE),
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
          function(x, norm, exact = FALSE, inverse = solve(x),
                   warn = TRUE, ...) {
              switch(norm,
                     "2" =,
                     "M" =, "m" =,
                     "F" =, "f" =, "E" =, "e" =
                         exact <- TRUE,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =
                         NULL,
                     stop(gettext("invalid %s=\"%s\"",
                                  "norm", norm),
                          domain = NA))
              d <- x@Dim
              if ((m <- d[1L]) == 0L || (n <- d[2L]) == 0L)
                  return(Inf)
              x <- .M2kind(x, ",")
              if (exact)
              switch(norm,
                     "2" =,
                         {
                             s <- svd(x, nu = 0L, nv = 0L)[["d"]]
                             if (s[1L] > 0) s[length(s)]/s[1L] else 0
                         },
                     "M" =, "m" =,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =,
                     "F" =, "f" =, "E" =, "e" =
                         1/(norm(x, type = norm) * norm(inverse, type = norm)))
              else if (m == n)
              rcond(.sparse2dense(x), norm = norm, exact = FALSE, ...)
              else {
              if (warn)
              warning(gettextf("rcond(%s, %s) via QR factorization; dubious for non-square '%s', %s=\"%s\"",
                               "x", "norm", "x", "norm", norm),
                      domain = NA)
              R <- triu(qr(if (m < n) ct(x) else x)@R[seq_len(if (m < n) m else n), , drop = FALSE])
              rcond(R, norm = norm, exact = FALSE, ...)
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
              switch(norm[1L],
                     "2" =,
                     "M" =, "m" =,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =
                         if (nonunit) {
                             r <- range(abs(y))
                             if (r[2L] > 0) r[1L]/r[2L] else 0
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
                  switch(norm,
                         "2" =,
                         "M" =, "m" =,
                         "O" =, "o" =, "1" =,
                         "I" =, "i" =
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
              switch(norm,
                     "2" =,
                     "M" =, "m" =,
                     "O" =, "o" =, "1" =,
                     "I" =, "i" =
                         1,
                     "F" =, "f" =, "E" =, "e" =
                         1/n,
                     stop(gettext("invalid %s=\"%s\"",
                                  "norm", norm[1L]),
                          domain = NA))
          })
