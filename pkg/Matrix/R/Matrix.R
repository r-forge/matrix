## METHODS FOR CLASS: Matrix (virtual)
## mother class containing all matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diff", signature(x = "Matrix"),
          ## Mostly cut and paste of 'base::diff.default' :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) > 1L ||
                  lag < 1L || differences < 1L)
                  stop("'lag' and 'differences' must be integers >= 1")
              if(lag * differences >= x@Dim[1L])
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences)) {
                  m <- x@Dim[1L]
                  x <- x[i1, , drop = FALSE] -
                      x[-m:-(m - lag + 1L), , drop = FALSE]
              }
              x
          })

if(FALSE) { ## still does not work for c(1, Matrix(2))
## For the same reason (and just in case) also do both S3 and S4 here:
c.Matrix <- function(...) unlist(lapply(list(...), as.vector))
## NB: Must use   signature  '(x, ..., recursive = FALSE)' :
setMethod("c", "Matrix", function(x, ..., recursive) c.Matrix(x, ...))
## The above is not sufficient for  c(NA, 3:2, <Matrix>, <matrix>)
setMethod("c", "numMatrixLike", function(x, ..., recursive) c.Matrix(x, ...))
}# not yet
