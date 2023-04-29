## METHODS FOR CLASS: MatrixFactorization (virtual)
## mother class containing all matrix factorizations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(FALSE) {
## MJ: These would make a lot of sense, but we already have
##     * as(      <denseLU>, "Matrix")
##     * as( <BunchKaufman>, "Matrix")
##     * as(<pBunchKaufman>, "Matrix")
##     * as(    <CHMfactor>, "Matrix")
##     _not_ returning the factorized matrix ...
##
##     Remove those and replace them with more specific coercions ??
##     * as(      <denseLU>, "dgeMatrix")
##     * as( <BunchKaufman>, "dtrMatrix")
##     * as(<pBunchKaufman>, "dtpMatrix")
##     * as(    <CHMfactor>, "dtCMatrix")
setAs("MatrixFactorization", "Matrix",
      function(from) {
          n <- length(x <- expand2(from))
          to <- x[[1L]]
          if(n >= 2L) for(i in 2:n) to <- to %*% x[[i]]
          to
      })

setAs("MatrixFactorization", "matrix",
      function(from) as(as(from, "Matrix"), "matrix"))
}

setMethod("dim", signature(x = "MatrixFactorization"),
          function(x) x@Dim)

setMethod("length", "MatrixFactorization",
          function(x) prod(x@Dim))

setMethod("dimnames", signature(x = "MatrixFactorization"),
          function(x) x@Dimnames)

setMethod("dimnames<-", signature(x = "MatrixFactorization", value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "MatrixFactorization", value = "NULL"),
	  function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("unname", signature(obj = "MatrixFactorization"),
	  function(obj, force = FALSE) {
              obj@Dimnames <- list(NULL, NULL)
              obj
          })

setMethod("show", "MatrixFactorization",
	  function(object) {
	      cat("matrix factorization of ")
	      str(object)
	  })

setMethod("show", "CholeskyFactorization",
          function(object) {
              cat("Cholesky factorization of ")
              str(object)
          })

setMethod("show", "BunchKaufmanFactorization",
          function(object) {
              cat("Bunch-Kaufman factorization of ")
              str(object)
          })

setMethod("show", "SchurFactorization",
	  function(object) {
	      cat("Schur factorization of ")
	      str(object)
	  })

setMethod("show", "LU",
          function(object) {
              cat("LU factorization of ")
              str(object)
          })

setMethod("show", "QR",
          function(object) {
              cat("QR factorization of ")
              str(object)
          })

## 'Cholesky' and 'pCholesky' represent matrices and we "show" them accordingly

setMethod("show", "Cholesky",
          function(object) prMatrix(object))

setMethod("show", "pCholesky",
          function(object) prMatrix(object))
