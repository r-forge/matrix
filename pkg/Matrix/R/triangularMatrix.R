## Methods for virtual class "triangularMatrix" of triangular matrices
.tM.subclasses <- names(getClassDef("triangularMatrix")@subclasses)

## ~~~~ COERCIONS TO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("Matrix", "triangularMatrix", ..M2tri)
setAs("matrix", "triangularMatrix", ..M2tri)

## MJ: prefer more general methods above going via triu() and tril(),
##     for which specialized methods should be defined
if(FALSE) {
setAs(   "matrix", "triangularMatrix", function(from) mat2tri(from))
setAs("dgeMatrix", "triangularMatrix", function(from) asTri(from, "dtrMatrix"))
setAs("lgeMatrix", "triangularMatrix", function(from) asTri(from, "ltrMatrix"))
setAs("ngeMatrix", "triangularMatrix", function(from) asTri(from, "ntrMatrix"))
}


## ~~~~ COERCIONS FROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: prefer more general method in ./symmetricMatrix going via
##     forceSymmetric(), for which specialized methods should be defined
if(FALSE) {
setAs("triangularMatrix", "symmetricMatrix",
      function(from) as(as(from, "generalMatrix"), "symmetricMatrix"))
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## In order to evade method dispatch ambiguity (with [CTR]sparse* and dense*),
## but still remain "general"
## we use this hack instead of signature  x = "triangularMatrix" :

## NB: denseMatrix goes via C utility R_dense_band() from ../src/dense.c
##     which obviates the need for .tri[ul].tr() above ...
.dM.subclasses <- names(getClassDef("denseMatrix")@subclasses)
for (.cl in grep("^.t.Matrix$",
                 setdiff(.tM.subclasses, .dM.subclasses), value = TRUE)) {
    setMethod("tril", signature(x = .cl), .tril.tr)
    setMethod("triu", signature(x = .cl), .triu.tr)
}
rm(.cl, .dM.subclasses)

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

rm(.tM.subclasses)
