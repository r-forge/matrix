## METHODS FOR CLASS: MatrixFactorization (virtual)
## mother class containing all matrix factorizations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim", "MatrixFactorization", function(x) x@Dim)

setMethod("show", "MatrixFactorization",
	  function(object) {
	      cat("matrix factorization of ")
	      str(object)
	  })

setMethod("show", "LU",
          function(object) {
              cat("LU factorization of ")
              str(object)
          })

setMethod("show", "sparseQR", # no virtual class "QR" yet
          function(object) {
              cat("QR factorization of ")
              str(object)
          })

setMethod("show", "CholeskyFactorization",
          function(object) {
              cat("Cholesky factorization of ")
              str(object)
          })

## 'Cholesky' and 'pCholesky' represent matrices and we "show" them accordingly

setMethod("show", "Cholesky",
          function(object) prMatrix(object))

setMethod("show", "pCholesky",
          function(object) prMatrix(object))

setMethod("show", "BunchKaufman",
	  function(object) {
	      cat("Bunch-Kaufman factorization of ")
	      str(object)
	  })

setMethod("show", "pBunchKaufman",
	  function(object) {
	      cat("Bunch-Kaufman factorization of ")
	      str(object)
	  })

setMethod("show", "Schur",
	  function(object) {
	      cat("Schur factorization of ")
	      str(object)
	  })
