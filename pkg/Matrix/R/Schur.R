## METHODS FOR GENERIC: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Schur", c(x = "denseMatrix"),
          function(x, vectors = TRUE, direct = FALSE, ...) {
              if (direct && .M.shape(x) == "t") {
                  x <- .M2kind(x, ",")
                  z <- is.complex(x@x)
                  r <- new(if (z) "zdenseSchur" else "ddenseSchur")
                  r@Dim <- d <- x@Dim
                  r@Dimnames <- x@Dimnames
                  if ((n <- d[1L]) > 0L) {
                  if (vectors)
                  r@vectors <- `[<-`(vector(typeof(x@x), prod(d)), if (x@uplo == "U") seq.int(from = 1L, by = n + 1, length.out = n) else seq.int(from = n, by = n - 1, length.out = n), 1)
                  r@values <- diag(x, names = FALSE)
                  r@x <- if (x@uplo == "U") .M2v(x) else { perm <- n:1L; as.vector(.M2m(x)[perm, perm]) }
                  }
                  r
              } else .Call(R_dense_schur, x, 1L, vectors)
          })

setMethod("Schur", c(x = "sparseMatrix"),
          function(x, vectors = TRUE, ...)
              Schur(.M2unpacked(x), vectors = vectors, ...))

setMethod("Schur", c(x = "matrix"),
          function(x, vectors = TRUE, ...)
              Schur(.m2dense(x, ",ge"), vectors = vectors, ...))
