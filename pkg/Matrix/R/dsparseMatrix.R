## METHODS FOR CLASS: dsparseMatrix (virtual)
## sparse matrices with 'x' slot of type "double"
## ... but _excluding_ ddiMatrix (FIXME?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.finite", signature(x = "dsparseMatrix"),
	  function(x) {
	      if(any(!is.finite(x@x))) {
                  r <- allTrueMat(x, packed = FALSE)
                  x <- as(as(as(x,"CsparseMatrix"), "dgCMatrix"),"dgTMatrix")
		  notF <- which(!is.finite(x@x))
		  r[cbind(x@i[notF], x@j[notF]) + 1L] <- FALSE
                  r
	      }
              else allTrueMat(x)
          })

setMethod("is.infinite", signature(x = "dsparseMatrix"),
	  function(x) {
	      if(any((isInf <- is.infinite(x@x)))) {
		  cld <- getClassDef(class(x))
		  if(extends(cld, "triangularMatrix") && x@diag == "U")
		      isInf <- is.infinite((x <- .diagU2N(x, cld))@x)
		  r <- as(x, "lMatrix") # will be "lsparseMatrix" - *has* x slot
		  r@x <- if(length(isInf) == length(r@x)) isInf else is.infinite(r@x)
		  if(!extends(cld, "CsparseMatrix"))
		      r <- as(r, "CsparseMatrix")
		  as(.Call(Csparse_drop, r, 0), "nMatrix") # a 'pattern matrix
	      }
	      else is.na_nsp(x)
	  })

## MJ: no longer needed ... now inherited from [CRT]sparseMatrix
if(FALSE) {
setMethod("symmpart", signature(x = "dsparseMatrix"),
          function(x) forceSymmetric(x + t(x)) / 2)
setMethod("skewpart", signature(x = "dsparseMatrix"),
          function(x) symmetrizeDimnames(x - t(x)) / 2)
} ## MJ

## MJ: no longer needed ... now inherited from Matrix
if(FALSE) {
setMethod("image", "dsparseMatrix",
	  function(x, ...) image(as(x, "dgTMatrix"), ...))
} ## MJ
