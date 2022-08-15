## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", signature(x = "generalMatrix"), # -> symmetricMatrix
	  function(x, pivot, ...)
              chol(.M2symm(x, checkDN = FALSE), pivot, ...))

setMethod("chol", signature(x = "symmetricMatrix"), # -> ds[ypCRT]Matrix
	  function(x, pivot, ...) {
              if(is(x, "nMatrix"))
                  stop("symbolic factorization of nMatrix via chol() is not yet implemented") # TODO
              chol(as(x, "dMatrix"), pivot, ...)
          })

setMethod("chol", signature(x = "triangularMatrix"), # ->diagonalMatrix
	  function(x, pivot, ...) {
              if(isDiagonal(x))
                  chol(.M2diag(x, check = FALSE), pivot, ...)
              else stop("chol(x) is undefined: 'x' is not symmetric")
          })

setMethod("chol", signature(x = "diagonalMatrix"),
	  function(x, pivot, ...) {
              x <- .diag2kind(x, "d")
              if(x@diag == "N") {
                  if(any(x@x < 0))
                      stop("chol(x) is undefined: 'x' is not positive definite")
                  x@x <- sqrt(x@x)
              }
              x
          })

setMethod("chol", signature(x = "dsyMatrix"),
          function(x, pivot, ...) {
              if(!is.null(ch <- x@factors[["Cholesky"]]))
                  return(ch) # use the cache
              tryCatch(.Call(dpoMatrix_chol, x),
                       error = function(e) stop("chol(x) is undefined: 'x' is not positive definite"))
          })

setMethod("chol", signature(x = "dspMatrix"),
          function(x, pivot, ...) {
              if(!is.null(ch <- x@factors[["pCholesky"]]))
                  return(ch) # use the cache
              tryCatch(.Call(dppMatrix_chol, x),
                       error = function(e) stop("chol(x) is undefined: 'x' is not positive definite"))
          })

setMethod("chol", signature(x = "dsCMatrix"),
	  function(x, pivot = FALSE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
		   return(ch) # use the cache
              tryCatch(.Call(dsCMatrix_chol, x, pivot),
                       error = function(e) stop("chol(x) is undefined: 'x' is not positive definite"))
          })

setMethod("chol", signature(x = "dsRMatrix"),
	  function(x, pivot = FALSE, cache = TRUE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
                  return(ch) # use the cache
              ch <- chol(.tCR2RC(x), pivot, ...)
              if(cache) .set.factors(x, nm, ch) else ch
          })

setMethod("chol", signature(x = "dsTMatrix"),
	  function(x, pivot = FALSE, cache = TRUE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
                  return(ch) # use the cache
              ch <- chol(.T2C(x), pivot, ...)
              if(cache) .set.factors(x, nm, ch) else ch
          })


## METHODS FOR GENERIC: Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Cholesky", signature(A = "denseMatrix"),
	  function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              stop("Cholesky(A) is implemented for sparseMatrix 'A' only; consider chol(A) instead"))

setMethod("Cholesky", signature(A = "sparseMatrix"), # ->dsCMatrix
	  function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              Cholesky(..sparse2d(.M2symm(as(A, "CsparseMatrix"))),
                       perm = perm, LDL = LDL, super = super, Imult = Imult,
                       ...))

setMethod("Cholesky", signature(A = "nsparseMatrix"),
	  function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              stop("symbolic factorization of nsparseMatrix via Cholesky() is not yet implemented")) # TODO

setMethod("Cholesky", signature(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              .Call(dsCMatrix_Cholesky, A, perm, LDL, super, Imult))


## METHODS FOR GENERIC: chol2inv
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol2inv", signature(x = "denseMatrix"), # ->dtrMatrix
	  function(x, ...)
              chol2inv(unpack(.M2tri(.dense2kind(x, "d"))), ...))

setMethod("chol2inv", signature(x = "dtrMatrix"),
	  function (x, ...) {
	      if(x@diag == "U")
                  x <- .dense.diagU2N(x)
	      .Call(dtrMatrix_chol2inv, x)
	  })

setMethod("chol2inv", signature(x = "sparseMatrix"),
	  function (x, ...) {
	      chkDots(..., which.call = -2L)
	      tcrossprod(solve(.M2tri(x)))
	  })

setMethod("chol2inv", signature(x = "diagonalMatrix"),
	  function (x, ...) {
	      chkDots(..., which.call = -2L)
	      tcrossprod(solve(x))
	  })

setMethod("chol2inv", signature(x = "CHMfactor"),
	  function (x, ...) {
	      chkDots(..., which.call = -2L)
	      solve(x, system = "A")
	  })
