## METHODS FOR GENERIC: forceSymmetric
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("forceSymmetric", c(x = "matrix"),
          function(x, uplo = "U", trans = "C", ...)
              .m2dense(x, ".sy", uplo = uplo, trans = trans))

setMethod("forceSymmetric", c(x = "denseMatrix"),
          function(x, uplo = NULL, trans = "C", ...)
              .Call(R_dense_force_symmetric, x, uplo, trans))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("forceSymmetric", c(x = .cl),
          function(x, uplo = NULL, trans = "C", ...)
              .Call(R_sparse_force_symmetric, x, uplo, trans))

setMethod("forceSymmetric", c(x = "diagonalMatrix"),
          function(x, uplo = "U", trans = "C", ...)
              .diag2sparse(x, ".", "s", "C", uplo, trans))

setMethod("forceSymmetric", c(x = "indMatrix"),
          function(x, ...)
              forceSymmetric(.ind2sparse(x), ...))

rm(.cl)
