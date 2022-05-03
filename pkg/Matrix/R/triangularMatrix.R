#### Methods for the virtual class 'triangularMatrix' of triangular matrices
#### Note that specific methods are in (8 different) ./?t?Matrix.R

setAs("triangularMatrix", "symmetricMatrix",
      function(from) as(as(from, "generalMatrix"), "symmetricMatrix"))

setAs("dgeMatrix", "triangularMatrix", function(from) asTri(from, "dtrMatrix"))
setAs("lgeMatrix", "triangularMatrix", function(from) asTri(from, "ltrMatrix"))
setAs("ngeMatrix", "triangularMatrix", function(from) asTri(from, "ntrMatrix"))

setAs("matrix", "triangularMatrix", function(from) mat2tri(from))

.tril.tr <- function(x, k = 0, ...) {  # are always square
    k <- as.integer(k[1])
    dd <- dim(x)
    stopifnot(-dd[1] <= k, k <= dd[1])  # had k <= 0
    if(k == 0 && x@uplo == "L") x
    else { ## more to do
        if(x@diag == "U") x <- .diagU2N(x, class(x), checkDense = TRUE)
        callNextMethod()
    }
}

.triu.tr <- function(x, k = 0, ...) {  # are always square
    k <- as.integer(k[1])
    dd <- dim(x)
    stopifnot(-dd[1] <= k, k <= dd[1])  # had k >= 0
    if(k == 0 && x@uplo == "U") x
    else { ## more to do
        if(x@diag == "U") x <- .diagU2N(x, class(x), checkDense = TRUE)
        callNextMethod()
    }
}

## In order to evade method dispatch ambiguity (with [CTR]sparse* and ddense*),
## but still remain "general"
## we use this hack instead of signature  x = "triangularMatrix" :

trCls <- names(getClass("triangularMatrix")@subclasses)
trCls. <- grep("^.t.Matrix$", trCls, value = TRUE) # not "p?Cholesky", etc.
for(cls in trCls.) {
    setMethod("tril", cls, .tril.tr)
    setMethod("triu", cls, .triu.tr)
}
rm(trCls, trCls., cls)

setMethod("isTriangular", signature(object = "triangularMatrix"),
          function(object, upper = NA, ...) {
              if(is.na(upper))
                  `attr<-`(TRUE, "kind", object@uplo)
              else
                  object@uplo == (if(upper) "U" else "L") || isDiagonal(object)
          })

## NB: [dz]t.Matrix should _not_ use this method as it does not
## tolerate numerical fuzz
setMethod("isSymmetric", signature(object = "triangularMatrix"),
	  function(object, checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...) check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              isDiagonal(object)
          })

cholTrimat <- function(x, ...) {
    if(isDiagonal(x))
        cholDiag(as(x, "diagonalMatrix"), ...)
    else stop("'x' is not symmetric -- chol() undefined.")
}
setMethod("chol", signature(x = "dtCMatrix"), cholTrimat)
setMethod("chol", signature(x = "dtTMatrix"), cholTrimat)
setMethod("chol", signature(x = "dtRMatrix"), cholTrimat)
## setMethod("chol", signature(x = "triangularMatrix"), cholTrimat)

