## METHODS FOR GENERIC: forceTriangular
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceTriangular", c(x = "matrix"),
          function(x, uplo = "U", diag = "N", ...)
              .m2dense(x, ".tr", uplo = uplo, diag = diag))

setMethod("forceTriangular", c(x = "denseMatrix"),
          function(x, uplo = NULL, diag = NULL, ...) {
              shape <- .M.shape(x)
              uplo <-
                  if (is.null(uplo))
                      (if (shape == "g") "U" else x@uplo)
                  else if (uplo == "U") "U" else "L"
              diag <-
                  if (is.null(diag))
                      (if (shape != "t") "N" else x@diag)
                  else if (diag == "N") "N" else "U"
              x <-
                  if (uplo == "U")
                      .Call(R_dense_band, x, 0L, NULL)
                  else .Call(R_dense_band, x, NULL, 0L)
              if (diag != x@diag) {
                  if (diag == "N")
                      diag(x) <- TRUE
                  else x@diag <- diag
              }
              x
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("forceTriangular", c(x = .cl),
          function(x, uplo = NULL, diag = NULL, ...) {
              shape <- .M.shape(x)
              uplo <-
                  if (is.null(uplo))
                      (if (shape == "g") "U" else x@uplo)
                  else if (uplo == "U") "U" else "L"
              diag <-
                  if (is.null(diag))
                      (if (shape != "t") "N" else x@diag)
                  else if (diag == "N") "N" else "U"
              k <-
                  if (diag == "N")
                      0L
                  else if (uplo == "U")
                      (if (x@Dim[2L] > 0L) 1L else 0L)
                  else (if (x@Dim[1L] > 0L) -1L else 0L)
              x <-
                  if (uplo == "U")
                      .Call(R_sparse_band, x, k, NULL)
                  else .Call(R_sparse_band, x, NULL, k)
              if (diag != x@diag) {
                  if (diag == "N")
                      diag(x) <- TRUE
                  else x@diag <- diag
              }
              x
          })

setMethod("forceTriangular", c(x = "diagonalMatrix"),
          function(x, uplo = "U", diag = NULL, ...) {
              x <- forceDiagonal(x, diag = diag)
              .diag2sparse(x, ".", "t", "C", uplo)
          })

setMethod("forceTriangular", c(x = "indMatrix"),
          function(x, ...)
              forceTriangular(.M2kind(x, "n"), ...))

rm(.cl)
