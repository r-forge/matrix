## METHODS FOR GENERIC: forceDiagonal
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceDiagonal", c(x = "matrix"),
          function(x, diag = "N", ...) {
              kind <- .M.kind(x)
              if (kind == "i")
                  kind <- "d"
              n <- min(d <- dim(x))
              dn <- dimnames(x)
              r <- new(paste0(kind, "diMatrix"))
              r@Dim <- c(n, n)
              if (!is.null(dn)) {
              if (max(d) != n &&
                  !is.null(tmp <- dn[[w <- which.max(d)]]))
                  dn[[w]] <- tmp[seq_len(n)]
              r@Dimnames <- dn
              }
              if (diag == "N") {
                  y <- diag(x, names = FALSE)
                  r@x <- if (is.integer(y)) as.double(y) else y
              } else r@diag <- "U"
              r
          })

for (.cl in c("denseMatrix", paste0(c("C", "R", "T"), "sparseMatrix")))
setMethod("forceDiagonal", c(x = .cl),
          function(x, diag = NULL, ...) {
              diag <-
                  if (is.null(diag))
                      (if (.M.shape(x) != "t") "N" else x@diag)
                  else if (diag == "N") "N" else "U"
              kind <- .M.kind(x)
              n <- min(d <- x@Dim)
              dn <- dimnames(x)
              r <- new(paste0(kind, "diMatrix"))
              r@Dim <- c(n, n)
              if (max(d) != n &&
                  !is.null(tmp <- dn[[w <- which.max(d)]]))
                  dn[[w]] <- tmp[seq_len(n)]
              r@Dimnames <- dn
              if (diag == "N")
                  r@x <- diag(x, names = FALSE)
              else r@diag <- "U"
              r
          })

setMethod("forceDiagonal", c(x = "diagonalMatrix"),
          function(x, diag = NULL, ...) {
              if (is.null(diag) || diag == x@diag)
                  return(x)
              r <- new(.M.class(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (diag == "N")
                  r@x <- rep.int(as.vector(TRUE, typeof(r@x)), d[1L])
              else r@diag <- "U"
              r
          })

setMethod("forceDiagonal", c(x = "indMatrix"),
          function(x, ...)
              forceDiagonal(.M2kind(x, "n"), ...))

rm(.cl)
