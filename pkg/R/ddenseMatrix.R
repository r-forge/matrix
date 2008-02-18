### Define Methods that can be inherited for all subclasses

## This replaces many "d..Matrix" -> "dgeMatrix" ones
## >> but << needs all sub(sub(sub)) classes of "ddenseMatrix" listed
##   -----  in  ../src/Mutils.c

setAs("ddenseMatrix", "dgeMatrix",
      function(from) {
	  if (class(from) != "dgeMatrix")
	      .Call(dup_mMatrix_as_dgeMatrix, from)
	  else from
      })


## d(ouble) to l(ogical):
setAs("dgeMatrix", "lgeMatrix", function(from) d2l_Matrix(from, "dgeMatrix"))
setAs("dsyMatrix", "lsyMatrix", function(from) d2l_Matrix(from, "dsyMatrix"))
setAs("dspMatrix", "lspMatrix", function(from) d2l_Matrix(from, "dspMatrix"))
setAs("dtrMatrix", "ltrMatrix", function(from) d2l_Matrix(from, "dtrMatrix"))
setAs("dtpMatrix", "ltpMatrix", function(from) d2l_Matrix(from, "dtpMatrix"))

setAs("ddenseMatrix", "CsparseMatrix",
      function(from) {
	  if (class(from) != "dgeMatrix") # don't lose symmetry/triangularity/...
	      as_Csparse(from)
	  else .Call(dense_to_Csparse, from)
      })

## special case
setAs("dgeMatrix", "dgCMatrix",
      function(from) .Call(dense_to_Csparse, from))

setAs("matrix", "CsparseMatrix",
      function(from) {
	    if(is.numeric(from))
		.Call(dense_to_Csparse, .Call(dup_mMatrix_as_dgeMatrix, from))
	    else if(is.logical(from)) ## FIXME: this works, but maybe wastefully
                as(Matrix(from, sparse=TRUE), "CsparseMatrix")
	    else stop('not-yet-implemented coercion to "CsparseMatrix"')
      })


## special case needed in the Matrix function
setAs("matrix", "dgCMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .Call(dense_to_Csparse, from)
      })

setAs("numeric", "CsparseMatrix",
      function(from)
      .Call(dense_to_Csparse, .Call(dup_mMatrix_as_dgeMatrix, from)))

setMethod("as.numeric", signature(x = "ddenseMatrix"),
	  function(x, ...) as(x, "dgeMatrix")@x)

## -- see also ./Matrix.R  e.g., for a show() method

## These methods are the 'fallback' methods for all dense numeric
## matrices in that they simply coerce the ddenseMatrix to a
## dgeMatrix. Methods for special forms override these.

setMethod("norm", signature(x = "ddenseMatrix", type = "missing"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("norm", signature(x = "ddenseMatrix", type = "character"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix"), type))

setMethod("rcond", signature(x = "ddenseMatrix", type = "missing"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("rcond", signature(x = "ddenseMatrix", type = "character"),
          function(x, type, ...) callGeneric(as(x, "dgeMatrix"), type))

## Not really useful; now require *identical* class for result:
## setMethod("t", signature(x = "ddenseMatrix"),
## 	  function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("tcrossprod", signature(x = "ddenseMatrix", y = "missing"),
	  function(x, y = NULL) callGeneric(as(x, "dgeMatrix")))

setMethod("crossprod", signature(x = "ddenseMatrix", y = "missing"),
          function(x, y = NULL) callGeneric(as(x, "dgeMatrix")))

setMethod("diag", signature(x = "ddenseMatrix"),
          function(x, nrow, ncol) callGeneric(as(x, "dgeMatrix")))

setMethod("solve", signature(a = "ddenseMatrix", b = "missing"),
          function(a, b, ...) callGeneric(as(a, "dgeMatrix")))

setMethod("solve", signature(a = "ddenseMatrix", b = "ANY"),
          function(a, b, ...) callGeneric(as(a, "dgeMatrix"), b))

## General method for dense matrix multiplication in case specific methods
## have not been defined.
setMethod("%*%", signature(x = "ddenseMatrix", y = "ddenseMatrix"),
          function(x, y) .Call(dgeMatrix_matrix_mm,
                               .Call(dup_mMatrix_as_dgeMatrix, x), y, FALSE),
          valueClass = "dgeMatrix")

setMethod("lu", signature(x = "ddenseMatrix"),
          function(x, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("chol", signature(x = "ddenseMatrix", pivot = "ANY"), cholMat)

setMethod("determinant", signature(x = "ddenseMatrix", logarithm = "missing"),
          function(x, logarithm, ...) callGeneric(as(x, "dgeMatrix")))

setMethod("determinant", signature(x = "ddenseMatrix", logarithm = "logical"),
          function(x, logarithm, ...)
          callGeneric(as(x, "dgeMatrix"), logarithm))

## now done for "dMatrix":
## setMethod("expm", signature(x = "ddenseMatrix"),
##           function(x) callGeneric(as(x, "dgeMatrix")))

setMethod("Math",
          signature(x = "ddenseMatrix"),
          function(x) callGeneric(as(x, "dgeMatrix")))


### FIXME: band() et al should be extended from "ddense" to "dense" !
###        However, needs much work to generalize dup_mMatrix_as_dgeMatrix()
### --> use workaround below: go via "d"(ouble) and back

.trilDense <- function(x, k = 0, ...) {
    k <- as.integer(k[1])
    dd <- dim(x); sqr <- dd[1] == dd[2]
    stopifnot(-dd[1] <= k, k <= dd[1]) # had k <= 0
    ## returns "lower triangular" if k <= 0 && sqr
    .Call(ddense_band, x, -dd[1], k)
}
## NB: have extra tril(), triu() methods for symmetric ["dsy" and "dsp"] and
##     for triangular ["dtr" and "dtp"]
setMethod("tril", "ddenseMatrix", .trilDense)
setMethod("tril",	"matrix",
	  function(x, k = 0, ...) {
	      if(is.double(x)) .trilDense(x, k)
	      else {
		  r <- .trilDense(x, k)
		  storage.mode(r) <- storage.mode(x)
		  r
	      }})
setMethod("tril", "denseMatrix",# all but ddense*
	  function(x, k = 0, ...)
	      as(.trilDense(as(x, "dMatrix"), k), class(x)))

.triuDense <- function(x, k = 0, ...) {
    k <- as.integer(k[1])
    dd <- dim(x); sqr <- dd[1] == dd[2]
    stopifnot(-dd[1] <= k, k <= dd[1]) # had k >= 0
    ## returns "upper triangular" if k >= 0
    .Call(ddense_band, x, k, dd[2])
}
setMethod("triu", "ddenseMatrix", .triuDense)
setMethod("triu",	"matrix",
	  function(x, k = 0, ...) {
	      if(is.double(x)) .triuDense(x, k)
	      else {
		  r <- .triuDense(x, k)
		  storage.mode(r) <- storage.mode(x)
		  r
	      }})
setMethod("triu", "denseMatrix",# all but ddense*
	  function(x, k = 0, ...)
	      as(.triuDense(as(x, "dMatrix"), k), class(x)))

.bandDense <- function(x, k1, k2, ...) {
    k1 <- as.integer(k1[1])
    k2 <- as.integer(k2[1])
    dd <- dim(x); sqr <- dd[1] == dd[2]
    stopifnot(-dd[1] <= k1, k1 <= k2, k2 <= dd[1])
    r <- .Call(ddense_band, x, k1, k2)
    if (k1 < 0	&&  k1 == -k2  && isSymmetric(x)) ## symmetric
	as(r, paste(.M.kind(x), "syMatrix", sep=''))
    else
	r
}

setMethod("band", "ddenseMatrix", .bandDense)
setMethod("band",	"matrix",
	  function(x, k1, k2, ...) {
	      if(is.double(x)) .bandDense(x, k1, k2)
	      else {
		  r <- .bandDense(x, k1, k2)
		  storage.mode(r) <- storage.mode(x)
		  r
	      }})
setMethod("band", "denseMatrix",# all but ddense*
	  function(x, k1, k2, ...)
	      as(.bandDense(as(x, "dMatrix"), k1, k2), class(x)))

setMethod("symmpart", signature(x = "ddenseMatrix"),
	  function(x) .Call(ddense_symmpart, x))
setMethod("skewpart", signature(x = "ddenseMatrix"),
	  function(x) .Call(ddense_skewpart, x))

