## METHODS FOR GENERIC: diag
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diag", c(x = "denseMatrix"),
          function(x = 1, nrow, ncol, names = TRUE)
              .Call(R_dense_diag_get, x, names))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("diag", c(x = .cl),
          function(x = 1, nrow, ncol, names = TRUE)
              .Call(R_sparse_diag_get, x, names))

setMethod("diag", c(x = "diagonalMatrix"),
          function(x = 1, nrow, ncol, names = TRUE) {
              kind <- .M.kind(x)
              ans <-
              if (x@diag != "N") {
                  one <- switch(kind, "n" = , "l" = TRUE, "i" = 1L, "d" = 1, "z" = 1+0i)
                  rep.int(one, x@Dim[1L])
              } else {
                  y <- x@x
                  if (kind == "n" && anyNA(y)) y | is.na(y) else y
              }
              if (names &&
                  !any(vapply(dn <- x@Dimnames, is.null, FALSE)) &&
                  {
                      i <- seq_len(min(x@Dim))
                      identical(nms <- dn[[1L]][i], dn[[2L]][i])
                  })
                  names(ans) <- nms
              ans
          })

setMethod("diag", c(x = "indMatrix"),
          function(x = 1, nrow, ncol, names = TRUE) {
              if ((r <- min(x@Dim)) == 0L)
                  return(logical(0L))
              i <- seq_len(r)
              ans <- x@perm[i] == i
              if (names &&
                  !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                  identical(nms <- dn[[1L]][i], dn[[2L]][i]))
                  names(ans) <- nms
              ans
          })

setMethod("diag", c(x = "denseCholesky"),
          function(x = 1, nrow, ncol, names = TRUE)
              .Call(denseCholesky_diag_get, x, FALSE))

setMethod("diag", c(x = "sparseCholesky"),
          function(x = 1, nrow, ncol, names = TRUE)
              .Call(sparseCholesky_diag_get, x, FALSE))

rm(.cl)


## METHODS FOR GENERIC: diag<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diag<-", c(x = "denseMatrix"),
          function(x, value)
              .Call(R_dense_diag_set, x, value))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("diag<-", c(x = .cl),
          function(x, value)
              .Call(R_sparse_diag_set, x, value))

setMethod("diag<-", c(x = "diagonalMatrix"),
          function(x, value) {
              tt <- c(typeof(x@x), typeof(value))
              tt. <- match(tt, c("logical", "integer", "double", "complex"), 0L)
              if (!tt.[2L])
                  stop(gettextf("replacement diagonal has incompatible type \"%s\"",
                                tt[2L]),
                       domain = NA)
              n <- (d <- x@Dim)[2L]
              if (all(length(value) != c(1L, n)))
                  stop("replacement diagonal has wrong length")
 ### FIXME: keep  x@diag == "U"  if(all(value == as1(x@x)))
             if (tt.[1L] < tt.[2L])
                  kind <- c("l", "i", "d", "z")[tt.[2L]]
              else {
                  kind <- .M.kind(x)
                  value <- as.vector(value, tt[1L])
              }
              ans <- new(paste0(kind, "diMatrix"))
              ans@Dim <- d
              ans@Dimnames <- x@Dimnames
              ans@x <- rep_len(value, n)
              ans
          })

setMethod("diag<-", c(x = "indMatrix"),
          function(x, value)
              `diag<-`(.M2kind(x, "n"), value))

rm(.cl)
