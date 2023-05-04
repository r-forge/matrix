## METHODS FOR GENERIC: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FIXME: C-level code should scan first for non-finite values
##        since R-level test needs an allocation

setMethod("Schur", signature(x = "Matrix", vectors = "missing"),
          function(x, vectors = TRUE, ...) Schur(x, TRUE, ...))

setMethod("Schur", signature(x = "matrix", vectors = "missing"),
          function(x, vectors = TRUE, ...) Schur(x, TRUE, ...))

setMethod("Schur", signature(x = "dgeMatrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...) {
              if(!all(is.finite(x@x)))
                  stop("'x' has non-finite values")
              cl <- .Call(dgeMatrix_Schur, x, vectors, TRUE)
              if(all(cl$WI == 0)) {
                  vals <- cl$WR
                  T <- triu(cl$T)
              } else {
                  vals <- complex(real = cl$WR, imaginary = cl$WI)
                  T <- .m2ge(cl$T)
              }
              if(vectors)
                  new("Schur", Dim = x@Dim, Dimnames = x@Dimnames,
                      Q = .m2ge(cl$Z), T = T, EValues = vals)
              else list(T = T, EValues = vals)
          })

setMethod("Schur", signature(x = "dsyMatrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...) {
              e <- eigen(x, symmetric = TRUE, only.values = !vectors)
              vals <- as.double(e$values)
              T <- new("ddiMatrix", Dim = x@Dim, x = vals)
              if(vectors)
                  new("Schur", Dim = x@Dim, Dimnames = symmDN(x@Dimnames),
                      Q = .m2ge(e$vectors), T = T, EValues = vals)
              else list(T = T, EValues = vals)
          })

setMethod("Schur", signature(x = "matrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...) {
              storage.mode(x) <- "double"
              if(!all(is.finite(x)))
                  stop("'x' has non-finite values")
              cl <- .Call(dgeMatrix_Schur, x, vectors, FALSE)
              vals <-
                  if(all(cl$WI == 0))
                      cl$WR
                  else complex(real = cl$WR, imaginary = cl$WI)
              if(vectors)
                  list(Q = cl$Z, T = cl$T, EValues = vals)
              else list(T = cl$T, EValues = vals)
          })

## FIXME: don't coerce from sparse to dense
setMethod("Schur", signature(x = "generalMatrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...)
              Schur(as(as(x, "dMatrix"), "unpackedMatrix"), vectors, ...))

## FIXME: don't coerce from sparse to dense
setMethod("Schur", signature(x = "symmetricMatrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...)
              Schur(as(as(x, "dMatrix"), "unpackedMatrix"), vectors, ...))

## Giving the _unsorted_ Schur factorization
setMethod("Schur", signature(x = "diagonalMatrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...) {
              d <- x@Dim
              if(x@diag != "N") {
                  vals <- rep.int(1, d[1L])
                  T <- new("ddiMatrix", Dim = d, diag = "U")
              } else {
                  vals <- x@x
                  if(!all(is.finite(vals)))
                      stop("'x' has non-finite values")
                  T <- new("ddiMatrix", Dim = d, x = vals)
              }
              if(vectors) {
                  Q <- new("ddiMatrix", Dim = d, diag = "U")
                  new("Schur", Dim = d, Dimnames = x@Dimnames,
                      Q = Q, T = T, EValues = vals)
              } else list(T = T, EValues = vals)
          })

setMethod("Schur", signature(x = "triangularMatrix", vectors = "logical"),
	  function(x, vectors = TRUE, ...) {
              cld <- getClassDef(class(x))
              if(!extends(cld, "nMatrix") &&
                 (anyNA(x) || (extends(cld, "dMatrix") && any(is.infinite(x)))))
                  ## any(is.finite(<sparseMatrix>)) would allocate too much
                  stop("'x' has non-finite values")
              n <- (d <- x@Dim)[1L]
              vals <- diag(x, names = FALSE)
              if(x@uplo == "U" || n == 0L) {
                  if(vectors) {
                      Q <- new("ddiMatrix", Dim = d, diag = "U")
                      new("Schur", Dim = d, Dimnames = x@Dimnames,
                          Q = Q, T = x, EValues = vals)
                  } else list(T = x, EValues = vals)
              } else {
                  perm <- n:1L
                  vals <- vals[perm]
                  T <- triu(x[perm, perm, drop = FALSE])
                  if(vectors) {
                      Q <- new("pMatrix", Dim = d, perm = perm)
                      new("Schur", Dim = d, Dimnames = x@Dimnames,
                          Q = Q, T = T, EValues = vals)
                  } else list(T = x, EValues = vals)
              }
          })


## METHODS FOR CLASS: Schur
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expand2", signature(x = "Schur"),
          function(x, ...) {
              Q  <- x@Q
              Q. <- t(Q)
              dn <- x@Dimnames
              Q @Dimnames <- c(dn[1L], list(NULL))
              Q.@Dimnames <- c(list(NULL), dn[2L])
              list(Q = Q, T = x@T, Q. = Q.)
          })
