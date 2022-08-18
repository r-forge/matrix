## METHODS FOR CLASS: dsparseMatrix (virtual)
## sparse matrices with 'x' slot of type "double"
## ... but _excluding_ ddiMatrix (FIXME?)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NB: classes without 'x' slot must have their own methods;
## see, e.g., ./nsparseMatrix.R and ./sparseVector.R
setMethod("anyNA", signature(x = "xMatrix"), function(x) anyNA(x@x))




setMethod("is.finite", signature(x = "dsparseMatrix"),
	  function(x) {
              cld <- getClassDef(class(x))
              symmetric <- extends(cld, "symmetricMatrix")
              if(!all(is.finite(x@x))) {
                  r <- allTrueMat(x, symmetric = symmetric, packed = FALSE)
                  x <- .sparse2g(if(extends(cld, "TsparseMatrix"))
                                     x
                                 else .CR2T(x))
                  notF <- which(!is.finite(x@x))
		  r[cbind(x@i[notF], x@j[notF]) + 1L] <- FALSE
                  r
	      } else allTrueMat(x, symmetric = symmetric, packed = TRUE)
          })

setMethod("is.infinite", signature(x = "dsparseMatrix"),
	  function(x) {
	      if(any(isInf <- is.infinite(x@x))) {
		  cld <- getClassDef(class(x))
		  if(extends(cld, "triangularMatrix") && x@diag == "U")
		      isInf <- is.infinite((x <- .diagU2N(x, cld))@x)
		  r <- .sparse2kind(x, "l", drop0 = FALSE)
		  r@x <- isInf
		  .sparse2kind(r, "n", drop0 = TRUE)
	      } else is.na_nsp(x)
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
