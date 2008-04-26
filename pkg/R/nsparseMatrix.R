#### Superclass Methods for all sparse nonzero-pattern matrices

setAs("CsparseMatrix", "nsparseMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from,
			   is(from, "triangularMatrix")))
setAs("CsparseMatrix", "nMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from,
			   is(from, "triangularMatrix")))

setAs("nsparseMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))

###------- Work via  as(*, lgC) : ------------

## For multiplication operations, sparseMatrix overrides other method
## selections.  Coerce a ddensematrix argument to a lsparseMatrix.

setMethod("%*%", signature(x = "nsparseMatrix", y = "ndenseMatrix"),
	  function(x, y) x %*% as(y, "nsparseMatrix"))

setMethod("%*%", signature(x = "ndenseMatrix", y = "nsparseMatrix"),
	  function(x, y) as(x, "nsparseMatrix") %*% y)

setMethod("crossprod", signature(x = "nsparseMatrix", y = "ndenseMatrix"),
	  function(x, y = NULL) crossprod(x, as(y, "nsparseMatrix")))

setMethod("crossprod", signature(x = "ndenseMatrix", y = "nsparseMatrix"),
	  function(x, y = NULL) crossprod(as(x, "nsparseMatrix"), y))

## and coerce lsparse* to lgC*
setMethod("%*%", signature(x = "nsparseMatrix", y = "nsparseMatrix"),
	  function(x, y) as(x, "ngCMatrix") %*% as(y, "ngCMatrix"))

setMethod("crossprod", signature(x = "nsparseMatrix", y = "nsparseMatrix"),
	  function(x, y = NULL)
	  crossprod(as(x, "ngCMatrix"), as(y, "ngCMatrix")))

setMethod("is.na", signature(x = "nsparseMatrix"), is.na_nsp)

setMethod("image", "lsparseMatrix", function(x, ...) image(as(x,"dMatrix")))
