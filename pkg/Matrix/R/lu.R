## METHODS FOR GENERIC: lu
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lu", c(x = "denseMatrix"),
          function(x, warnSing = TRUE, direct = FALSE, ...) {
              if (direct && .M.shape(x) == "t" &&
                  (x@uplo == "U" || x@diag == "U")) {
                  y <- .M2gen(x, ",")
                  z <- is.complex(y@x)
                  n <- (d <- y@Dim)[1L]
                  r <- new(if (z) "zdenseLU" else "ddenseLU")
                  r@Dim <- d
                  r@Dimnames <- y@Dimnames
                  r@perm <- seq_len(n)
                  r@x <- y@x
                  r
              } else .Call(R_dense_lu, x, if (warnSing) 1L else 0L)
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("lu", c(x = .cl),
          function(x, errSing = TRUE, order = NA_integer_, tol = 1,
                   direct = FALSE, ...) {
              if (direct && .M.shape(x) == "t" &&
                  ((up <- x@uplo == "U") || x@diag == "U")) {
                  y <- .M2kind(x, ",")
                  z <- is.complex(y@x)
                  n <- (d <- y@Dim)[1L]
                  r <- new(if (z) "zsparseLU" else "dsparseLU")
                  r@Dim <- d
                  r@Dimnames <- y@Dimnames
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
                  tmp <- new(if (z) "ztCMatrix" else "dtCMatrix")
                  tmp@Dim <- d
                  if (up)
                  tmp@uplo <- "L"
                  tmp@diag <- "U"
                  tmp@p <- integer(n + 1)
                  r@L <- if (up) tmp else .M2C(y)
                  r@U <- if (up) .M2C(y) else tmp
                  r
              } else .Call(R_sparse_lu, x, if (errSing) 2L else 0L, order, tol)
          })

setMethod("lu", c(x = "diagonalMatrix"),
          function(x, direct = TRUE, ...) {
              if (direct) {
                  y <- .M2kind(x, ",")
                  z <- is.complex(y@x)
                  n <- (d <- y@Dim)[1L]
                  r <- new(if (z) "zsparseLU" else "dsparseLU")
                  r@Dim <- d
                  r@Dimnames <- y@Dimnames
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
                  tmp <- new(if (z) "ztCMatrix" else "dtCMatrix")
                  tmp@Dim <- d
                  tmp@uplo <- "L"
                  tmp@diag <- "U"
                  tmp@p <- integer(n + 1)
                  r@L <- tmp
                  tmp@uplo <- "U"
                  if (x@diag == "N") {
                      tmp@diag <- "N"
                      tmp@p <- seq.int(from = 0L, length.out = n + 1)
                      tmp@x <- y@x
                  }
                  r@U <- tmp
                  r
              } else lu(.diag2sparse(x, ",", "g", "C"), ...)
          })

setMethod("lu", c(x = "indMatrix"),
          function(x, ...)
              lu(.ind2sparse(x, ",", "C"), ...))

setMethod("lu", c(x = "matrix"),
          function(x, ...)
              lu(.m2dense(x, ",ge"), ...))

rm(.cl)
