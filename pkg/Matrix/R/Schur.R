## METHODS FOR GENERIC: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Schur", c(x = "denseMatrix"),
          function(x, vectors = TRUE, ...) {
              x <- .M2kind(x, ",")
              switch(substr(.M.class(x), 2L, 3L),
                     "ge" = .Call(geMatrix_scf, x, TRUE, vectors),
                     "sy" = .Call(syMatrix_scf, x, TRUE, vectors),
                     "sp" = .Call(spMatrix_scf, x, TRUE, vectors),
                     "tr" = , "tp" =
                         {
                             r <- new("ddenseSchur")
                             r@Dim <- d <- x@Dim
                             r@Dimnames <- x@Dimnames
                             if((n <- d[1L]) > 0L) {
                             if(vectors)
                             r@vectors <- as.vector(new("pMatrix", Dim = d, perm = if(x@uplo == "U") 1L:n else n:1L), typeof(x@x))
                             r@values <- diag(x, names = FALSE)
                             r@x <- .M2v(if(x@uplo == "U") x else { perm <- n:1L; x[perm, perm] })
                             }
                             r
                         })
          })

setMethod("Schur", c(x = "sparseMatrix"),
          function(x, vectors = TRUE, ...)
              Schur(.M2unpacked(x), vectors = vectors, ...))

setMethod("Schur", c(x = "matrix"),
          function(x, vectors = TRUE, ...)
              Schur(.m2dense(x, ",ge"), vectors = vectors, ...))
