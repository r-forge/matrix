## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: all of the C-level *_chol() functions need a second look ...
##     these methods feel much more complicated than they need to be ...

setMethod("chol", signature(x = "generalMatrix"),
	  function(x, ...) {
              ch <- chol(.M2sym(x, checkDN = FALSE), ...)
              ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("chol", signature(x = "symmetricMatrix"),
	  function(x, ...)
              chol(as(x, "dMatrix"), ...))

setMethod("chol", signature(x = "triangularMatrix"),
	  function(x, ...) {
              if(isDiagonal(x))
                  chol(forceDiagonal(x), ...)
              else stop("chol(x) is undefined: 'x' is not symmetric")
          })

setMethod("chol", signature(x = "diagonalMatrix"),
	  function(x, ...) {
              x <- ..diag2d(x)
              if(x@diag == "N") {
                  if(any(x@x < 0))
                      stop("chol(x) is undefined: 'x' is not positive definite")
                  x@x <- sqrt(x@x)
              }
              x
          })

setMethod("chol", signature(x = "dsyMatrix"),
          function(x, ...) {
              if(!is.null(ch <- x@factors[["Cholesky"]]))
                  return(ch) # use the cache
              tryCatch(.Call(dpoMatrix_trf, if(x@uplo == "U") x else t(x), 2L),
                       error = function(e) stop("chol(x) is undefined: 'x' is not positive definite"))
          })

setMethod("chol", signature(x = "dspMatrix"),
          function(x, ...) {
              if(!is.null(ch <- x@factors[["pCholesky"]]))
                  return(ch) # use the cache
              tryCatch(.Call(dppMatrix_trf, if(x@uplo == "U") x else t(x), 2L),
                       error = function(e) stop("chol(x) is undefined: 'x' is not positive definite"))
          })

## FIXME: Is there is a simple way at C-level to cache L' rather than L,
##        given that chol() is documented to return L'?  Then we wouldn't
##        need t() here.  Of course, we could construct L' as a dtRMatrix
##        "for free" with .tCR2RC() ...
for(.cl in paste0("dg", c("C", "R", "T"), "Matrix"))
setMethod("chol", signature(x = .cl),
	  function(x, pivot = FALSE, cache = TRUE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
                  return(t(as(ch, "CsparseMatrix")))
              ch <- chol(y <- .M2sym(x, checkDN = FALSE), pivot = pivot, ...)
              ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
              if(cache)
                  ## dsCMatrix_chol() caches CHMfactor and returns dtCMatrix
                  .set.factors(x, nm, y@factors[[nm]])
              ch
          })
rm(.cl)

setMethod("chol", signature(x = "dsCMatrix"),
	  function(x, pivot = FALSE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
                  return(t(as(ch, "CsparseMatrix")))
              tryCatch(.Call(dsCMatrix_chol, x, pivot),
                       error = function(e) stop("chol(x) is undefined: 'x' is not positive definite"))
          })

setMethod("chol", signature(x = "dsRMatrix"),
	  function(x, pivot = FALSE, cache = TRUE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
                  return(t(as(ch, "CsparseMatrix")))
              ch <- chol(y <- .tCR2RC(x), pivot = pivot, ...)
              if(cache)
                  ## dsCMatrix_chol() caches CHMfactor and returns dtCMatrix
                  .set.factors(x, nm, y@factors[[nm]])
              ch
          })

setMethod("chol", signature(x = "dsTMatrix"),
	  function(x, pivot = FALSE, cache = TRUE, ...) {
              nm <- if(pivot) "sPdCholesky" else "spdCholesky"
              if(!is.null(ch <- x@factors[[nm]]))
                  return(t(as(ch, "CsparseMatrix")))
              ch <- chol(y <- .T2C(x), pivot = pivot, ...)
              if(cache)
                  ## dsCMatrix_chol() caches CHMfactor and returns dtCMatrix
                  .set.factors(x, nm, y@factors[[nm]])
              ch
          })


## METHODS FOR GENERIC: Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Cholesky", signature(A = "denseMatrix"),
	  function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              stop("Cholesky(A) is implemented for sparseMatrix 'A' only; consider chol(A) instead"))

setMethod("Cholesky", signature(A = "sparseMatrix"), # ->dsCMatrix
	  function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              Cholesky(..sparse2d(.M2sym(as(A, "CsparseMatrix"))),
                       perm = perm, LDL = LDL, super = super, Imult = Imult,
                       ...))

setMethod("Cholesky", signature(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...)
              .Call(dsCMatrix_Cholesky, A, perm, LDL, super, Imult))


## METHODS FOR GENERIC: chol2inv
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol2inv", signature(x = "denseMatrix"), # ->dtrMatrix
	  function(x, ...)
              chol2inv(unpack(.M2tri(..dense2d(x))), ...))

setMethod("chol2inv", signature(x = "dtrMatrix"),
	  function (x, ...) {
	      if(x@diag != "N")
                  x <- ..diagU2N(x)
	      .Call(dtrMatrix_chol2inv, x)
	  })

setMethod("chol2inv", signature(x = "sparseMatrix"),
	  function (x, ...)
	      tcrossprod(solve(.M2tri(x))))

setMethod("chol2inv", signature(x = "diagonalMatrix"),
	  function (x, ...)
	      tcrossprod(solve(x)))

setMethod("chol2inv", signature(x = "CHMfactor"),
	  function (x, ...)
	      solve(x, system = "A"))


## METHODS FOR CLASS: p?Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(L, L') or list(L1, D, L1'),
## where  A = L L' = L1 D L1'  and  L = L1 sqrt(D)
.def.unpacked <- .def.packed <- function(x, LDL = TRUE, ...) {
    x <- as(x, .CL)
    dn <- x@Dimnames
    up <- x@uplo == "U"
    nu <- x@diag == "N"
    if(LDL && nu) {
        L.ii <- diag(x, names = FALSE)
        x@x <- x@x / if(up) .UP else .LO
        x@diag <- "U"
    }
    L  <- if(up) t(x) else   x
    L. <- if(up)   x  else t(x)
    L @Dimnames <- c(dn[1L], list(NULL))
    L.@Dimnames <- c(list(NULL), dn[2L])
    if(LDL) {
        D <- new("ddiMatrix")
        D@Dim <- x@Dim
        if(nu)
            D@x <- L.ii * L.ii
        else D@diag <- "U"
        list(L1 = L, D = D, L1. = L.)
    } else list(L = L, L. = L.)
}

body(.def.unpacked) <-
    do.call(substitute,
            list(body(.def.unpacked),
                 list(.CL = "dtrMatrix",
                      .UP = quote(L.ii),
                      .LO = quote(rep(L.ii, each = x@Dim[1L])))))
body(.def.packed) <-
    do.call(substitute,
            list(body(.def.packed),
                 list(.CL = "dtpMatrix",
                      .UP = quote(L.ii[sequence.default(seq_len(x@Dim[1L]))]),
                      .LO = quote(rep.int(L.ii, seq.int(to = 1L, by = -1L, length.out = x@Dim[1L]))))))

## returning list(L, L') or list(L1, D, L1'), where A = L L' = L1 D L1'
setMethod("expand2", signature(x =  "Cholesky"), .def.unpacked)
setMethod("expand2", signature(x = "pCholesky"), .def.packed)

rm(.def.unpacked, .def.packed)
