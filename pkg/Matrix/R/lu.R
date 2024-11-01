## METHODS FOR GENERIC: lu
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lu", c(x = "matrix"),
          function(x, ...) lu(.m2dense(x, ",ge"), ...))

setMethod("lu", c(x = "denseMatrix"),
          function(x, ...) lu(.M2kind(x, ","), ...))

setMethod("lu", c(x = "dgeMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(geMatrix_trf, x, as.logical(warnSing)))

setMethod("lu", c(x = "dsyMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["denseLU"]]))
                  return(ch)
              r <- lu(.M2gen(x), ...)
              if(cache) .set.factor(x, "denseLU", r) else r
          })

setMethod("lu", c(x = "dspMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["denseLU"]]))
                  return(ch)
              r <- lu(.M2gen(x), ...)
              if(cache) .set.factor(x, "denseLU", r) else r
          })

for(.cl in c("dtrMatrix", "dtpMatrix"))
setMethod("lu", c(x = .cl),
          function(x, ...) {
              if(x@uplo == "U" || x@diag == "U") {
                  r <- new("ddenseLU")
                  r@Dim <- d <- x@Dim
                  r@Dimnames <- x@Dimnames
                  r@perm <- seq_len(d[1L])
                  r@x <- .M2gen(x)@x
                  r
              } else lu(.M2gen(x), ...)
          })
rm(.cl)

setMethod("lu", c(x = "sparseMatrix"),
          function(x, ...)
              lu(.M2kind(.M2C(x), ","), ...))

setMethod("lu", c(x = "dgCMatrix"),
          function(x, errSing = TRUE, order = NA_integer_, tol = 1, ...)
              .Call(gCMatrix_trf, x, if(errSing) 2L else 0L, order, tol))

setMethod("lu", c(x = "dsCMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2gen(x), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", c(x = "dtCMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("dsparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@Dim <- d
                  r@Dimnames <- x@Dimnames
                  r@L <- if(upper) y else x
                  r@U <- if(upper) x else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.M2gen(x), ...)
          })

setMethod("lu", c(x = "dgRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2C(x), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", c(x = "dsRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2gen(.tCRT(x)), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", c(x = "dtRMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("dsparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@Dim <- d
                  r@Dimnames <- x@Dimnames
                  r@L <- if(upper) y else .M2C(x)
                  r@U <- if(upper) .M2C(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.M2gen(.M2C(x)), ...)
          })

setMethod("lu", c(x = "dgTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2C(x), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", c(x = "dsTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["sparseLU"]]))
                  return(ch)
              r <- lu(.M2gen(.M2C(x)), ...)
              if(cache) .set.factor(x, "sparseLU", r) else r
          })

setMethod("lu", c(x = "dtTMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("dsparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@Dim <- d
                  r@Dimnames <- x@Dimnames
                  r@L <- if(upper) y else .M2C(x)
                  r@U <- if(upper) .M2C(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
                  r
              } else lu(.M2gen(.M2C(x)), ...)
          })

setMethod("lu", c(x = "diagonalMatrix"),
          function(x, ...) {
              x <- .M2kind(x, ",")
              n <- (d <- x@Dim)[1L]
              L <- new("dtCMatrix")
              r <- new("dsparseLU")
              L@Dim <- d
              L@uplo <- "L"
              L@diag <- "U"
              L@p <- integer(n + 1L)
              r@L <- L
              if(x@diag == "N") {
                  L@diag <- "N"
                  L@p <- seq.int(from = 0L, length.out = n + 1L)
                  L@x <- x@x
              }
              r@U <- L
              r@Dim <- d
              r@Dimnames <- x@Dimnames
              r@p <- r@q <- seq.int(from = 0L, length.out = n)
              r
          })
