## METHODS FOR GENERIC: qr
## pivoted QR factorization of dense and sparse matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## FIXME? We could have methods for generalMatrix, symmetricMatrix,
##        triangularMatrix, and diagonalMatrix instead?  We could
##        construct the result "directly" in the upper triangular
##        (incl. diagonal) cases, and cache it in the general and
##        symmetric cases.  See the methods for chol() in ./chol.R ...

setMethod("qr", signature(x = "sparseMatrix"),
	  function(x, ...)
              qr(..sparse2d(.sparse2g(as(x, "CsparseMatrix"))), ...))

setMethod("qr", signature(x = "dgCMatrix"),
          function(x, order = 3L, ...) {
              r <- .Call(dgCMatrix_orf, x, order, TRUE)
              if(n <- r@V@Dim[1L] - r@Dim[1L])
                  Matrix.msg(gettextf("matrix is structurally rank deficient; returning QR factorization of matrix augmented with %d rows of zeros", n),
                             domain = NA)
              r
          })
