## METHODS FOR GENERIC: lu
## pivoted LU factorization, returning denseLU or sparseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("lu", signature(x = "matrix"),
          function(x, ...) lu(.m2ge(x, "d"), ...))

setMethod("lu", signature(x = "denseMatrix"),
	  function(x, ...) lu(..dense2d(x), ...))

setMethod("lu", signature(x = "dgeMatrix"),
	  function(x, warnSing = TRUE, ...)
              .Call(dgeMatrix_trf, x, warnSing))

setMethod("lu", signature(x = "dsyMatrix"),
	  function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.dense2g(x), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dspMatrix"),
	  function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.dense2g(x), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

for(.cl in c("dtrMatrix", "dtpMatrix"))
setMethod("lu", signature(x = .cl),
	  function(x, ...) {
              if(x@uplo == "U" || x@diag == "U") {
                  r <- new("denseLU")
                  r@Dim <- d <- x@Dim
                  r@perm <- seq_len(d[1L])
                  r@x <- .dense2g(x)@x
                  r
              } else lu(.dense2g(x), ...)
          })
rm(.cl)

setMethod("lu", signature(x = "sparseMatrix"),
	  function(x, ...)
              lu(..sparse2d(as(x, "CsparseMatrix")), ...))

setMethod("lu", signature(x = "dgCMatrix"),
          function(x, errSing = TRUE, order = NA_integer_, tol = 1,
                   keep.dimnames = TRUE, ...)
              .Call(dgCMatrix_trf, x, errSing, keep.dimnames, order, tol))

setMethod("lu", signature(x = "dsCMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.sparse2g(x), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

setMethod("lu", "dtCMatrix",
	  function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- r@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@L <- if(upper) y else x
                  r@U <- if(upper) x else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.sparse2g(x), ...)
	  })

setMethod("lu", signature(x = "dgRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.CR2RC(x), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dsRMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.sparse2g(.tCR2RC(x)), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dtRMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- r@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@L <- if(upper) y else .CR2RC(x)
                  r@U <- if(upper) .CR2RC(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
              } else lu(.sparse2g(.CR2RC(x)), ...)
          })

setMethod("lu", signature(x = "dgTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.T2C(x), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dsTMatrix"),
          function(x, cache = TRUE, ...) {
              if(!is.null(ch <- x@factors[["LU"]]))
                  return(ch)
              r <- lu(.sparse2g(.T2C(x)), ...)
              if(cache) .set.factor(x, "LU", r) else r
          })

setMethod("lu", signature(x = "dtTMatrix"),
          function(x, ...) {
              if((upper <- x@uplo == "U") || x@diag == "U") {
                  n <- (d <- x@Dim)[1L]
                  r <- new("sparseLU")
                  y <- new("dtCMatrix")
                  y@Dim <- r@Dim <- d
                  y@uplo <- if(upper) "L" else "U"
                  y@diag <- "U"
                  y@p <- integer(n + 1L)
                  r@L <- if(upper) y else .T2C(x)
                  r@U <- if(upper) .T2C(x) else y
                  r@p <- r@q <- seq.int(from = 0L, length.out = n)
                  r
              } else lu(.sparse2g(.T2C(x)), ...)
          })

setMethod("lu", "diagonalMatrix",
	  function(x, ...) {
              n <- (d <- x@Dim)[1L]
              L <- new("dtCMatrix")
              r <- new("sparseLU")
              L@Dim <- r@Dim <- d
              L@uplo <- "L"
              L@diag <- "U"
              L@p <- integer(n + 1L)
              r@L <- L
              if(x@diag == "N") {
                  L@diag <- "N"
                  L@p <- seq.int(from = 0L, length.out = n + 1L)
                  L@x <- as.double(x@x)
              }
              r@U <- L
              r@p <- r@q <- seq.int(from = 0L, length.out = n)
              r
	  })


## METHODS FOR CLASS: denseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("denseLU", "dgeMatrix",
      function(from) {
          to <- new("dgeMatrix")
          to@Dim <- from@Dim
          to@Dimnames <- from@Dimnames
          to@x <- from@x
          to
      })

## returning list(P1', L, U), where A = P1' L U
setMethod("expand2", signature(x = "denseLU"),
          function(x, ...) .Call(denseLU_expand, x))

## returning list(L, U, P), where A = P L U
## MJ: for backwards compatibility
setMethod("expand", signature(x = "denseLU"),
          function(x, ...) {
              r <- .Call(denseLU_expand, x)[c(2L, 3L, 1L)]
              names(r) <- c("L", "U", "P")
              NN <- list(NULL, NULL)
              for(i in 1:3)
                  r[[i]]@Dimnames <- NN
              r
          })


## METHODS FOR CLASS: sparseLU
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(P1', L, U, P2'), where A = P1' L U P2'
setMethod("expand2", signature(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              p1 <- x@p
              p2 <- x@q
              P1. <- new("pMatrix",
                         Dim = d,
                         Dimnames = c(dn[1L], list(NULL)),
                         margin = 1L,
                         perm = invPerm(p1, zero.p = TRUE, zero.res = FALSE))
              P2. <- new("pMatrix",
                         Dim = d,
                         Dimnames = c(list(NULL), dn[2L]),
                         margin = 2L,
                         perm = if(length(p2)) invPerm(p2, zero.p = TRUE, zero.res = FALSE) else seq_len(d[1L]))
              L <- x@L
              U <- x@U
              if(L@diag == "N")
                  L <- ..diagN2U(L, sparse = TRUE)
              list(P1. = P1., L = L, U = U, P2. = P2.)
          })

## returning list(P, L, U, Q), where A = P' L U Q
## MJ: for backwards compatibility
setMethod("expand", signature(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              p <- x@p + 1L
              q <- x@q + 1L
              if(length(p) && !is.null(rn <- dn[[1L]]))
                  dn[[1L]] <- rn[p]
              if(length(q) && !is.null(cn <- dn[[2L]]))
                  dn[[2L]] <- cn[q]
              P <- new("pMatrix",
                       Dim = d,
                       perm = p)
              Q <- new("pMatrix",
                       Dim = d,
                       perm = if(length(q)) q else seq_len(d[1L]))
              L <- x@L
              U <- x@U
              L@Dimnames <- c(dn[1L], list(NULL))
              U@Dimnames <- c(list(NULL), dn[2L])
              list(P = P, L = L, U = U, Q = Q)
          })
