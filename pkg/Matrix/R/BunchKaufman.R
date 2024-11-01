## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", c(x = "dsyMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(syMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", c(x = "dspMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(spMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", c(x = "matrix"),
          function(x, uplo = "U", ...)
              BunchKaufman(.m2dense(x, ",sy", uplo), ...))
