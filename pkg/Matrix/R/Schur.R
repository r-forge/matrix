## METHODS FOR GENERIC: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Schur", c(x = "denseMatrix"),
          function(x, vectors = TRUE, ...) {
              x <- .M2kind(x, ",")
              switch(substr(.M.nonvirtual(x), 2L, 3L),
                     "ge" = .Call(geMatrix_scf, x, TRUE, vectors),
                     "sy" = .Call(syMatrix_scf, x, TRUE, vectors),
                     "sp" = .Call(spMatrix_scf, x, TRUE, vectors),
                     "tr" = , "tp" =
                         {
                             r <- new("Schur")
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


## METHODS FOR CLASS: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expand1", c(x = "Schur"),
          function(x, which, ...) {
              d <- x@Dim
              switch(which,
                     "T" =
                         {
                             values <- x@values
                             x <- x@x
                             if(d[1L] > 0L && length(x) == 0L)
                                 new("ddiMatrix", Dim = d, x = values)
                             else if(is.double(values))
                                 new("dtrMatrix", Dim = d, x = x)
                             else new("dgeMatrix", Dim = d, x = x)
                         },
                     "Q" = , "Q." =
                         {
                             vectors <- x@vectors
                             if(d[1L] > 0L && length(vectors) == 0L)
                                 stop("missing requisite Schur vectors")
                             Q <- new("dgeMatrix", Dim = d, x = vectors)
                             switch(which, "Q" = Q, "Q." = t(Q))
                         },
                     stop(gettextf("'%1$s' is not \"%2$s\", \"%3$s\", or \"%2$s.\"",
                                   "which", "Q", "T"),
                          domain = NA))
          })

setMethod("expand2", c(x = "Schur"),
          function(x, ...) {
              T <- expand1(x, "T")
              Q <- expand1(x, "Q")
              Q. <- t(Q)
              dn <- x@Dimnames
              Q @Dimnames <- c(dn[1L], list(NULL))
              Q.@Dimnames <- c(list(NULL), dn[2L])
              list(Q = Q, T = T, Q. = Q.)
          })
