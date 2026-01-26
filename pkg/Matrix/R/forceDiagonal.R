## METHODS FOR GENERIC: forceDiagonal
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceDiagonal", c(x = "matrix"),
          function(x, diag = NULL, ...) {
              d <- dim(x)
              if (d[1L] != d[2L])
                  stop("matrix is not square")
              if (is.null(diag) || diag == "N")
                  y <- diag(x, names = FALSE)
              if (is.null(diag))
                  diag <- if (anyNA(y) || any(y != 1L)) "N" else "U"
              r <- new(paste0(.M.kind(x), "diMatrix"))
              r@Dim <- d
              if (!is.null(dn <- dimnames(x)))
                  r@Dimnames <- dn
              if (diag == "N")
                  r@x <- y
              else r@diag <- "U"
              r
          })

for (.cl in c("denseMatrix", paste0(c("C", "R", "T"), "sparseMatrix")))
setMethod("forceDiagonal", c(x = .cl),
          function(x, diag = NULL, ...) {
              d <- x@Dim
              if (d[1L] != d[2L])
                  stop("matrix is not square")
              if (is.null(diag))
                  diag <- if (.M.shape(x) != "t") "N" else x@diag
              r <- new(paste0(.M.kind(x), "diMatrix"))
              r@Dim <- d
              r@Dimnames <- dimnames(x)
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
              forceDiagonal(.ind2sparse(x), ...))

rm(.cl)
