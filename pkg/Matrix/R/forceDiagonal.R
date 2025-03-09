## METHODS FOR GENERIC: forceDiagonal
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceDiagonal", c(x = "matrix"),
          function(x, diag = NULL, ...) {
              d <- dim(x)
              if (d[1L] != d[2L])
                  stop("matrix is not square")
              diax <- diag(x, names = FALSE)
              if(is.null(diag))
                  diag <- if(!anyNA(diax) && all(diax == as1(x[0L]))) "U" else "N"
              r <- new(paste0(.M.kind(x), "diMatrix"), Dim = d,
                       x = if (diag == "N") diax else x[0L], diag = diag)
              if(!is.null(dn <- dimnames(x)))
                  r@Dimnames <- dn
              r
          })

for (.cl in c("denseMatrix", paste0(c("C", "R", "T"), "sparseMatrix")))
setMethod("forceDiagonal", c(x = .cl),
          function(x, diag = NULL, ...) {
              d <- x@Dim
              if (d[1L] != d[2L])
                  stop("matrix is not square")
              dn <- dimnames(x)
              diag <-
                  if (is.null(diag))
                      (if (.M.shape(x) != "t") "N" else x@diag)
                  else if (diag == "N") "N" else "U"
              r <- new(paste0(.M.kind(x), "diMatrix"))
              r@Dim <- d
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
              forceDiagonal(.ind2sparse(x), ...))

rm(.cl)
