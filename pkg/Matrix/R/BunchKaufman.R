## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", c(x = "denseMatrix"),
          function(x, warnSing = TRUE, uplo = "U", trans = "C", ...)
              .Call(R_dense_bunchkaufman, x, if (warnSing) 1L else 0L, uplo, trans))

setMethod("BunchKaufman", c(x = "sparseMatrix"),
          function(x, ...)
              BunchKaufman(.M2unpacked(x), ...))

setMethod("BunchKaufman", c(x = "matrix"),
          function(x, ...)
              BunchKaufman(.m2dense(x, ",ge"), ...))
