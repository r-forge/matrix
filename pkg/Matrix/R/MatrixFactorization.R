## METHODS FOR CLASS: MatrixFactorization (virtual)
## mother class containing all matrix factorizations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(FALSE) {
## MJ: not yet ... existing as(<CHMfactor>, "Matrix") must become defunct first?
setAs("MatrixFactorization", "Matrix",
      function(from) {
          n <- length(x <- expand2(from))
          to <- x[[1L]]
          if(n >= 2L) for(i in 2:n) to <- to %*% x[[i]]
          to
      })
}
