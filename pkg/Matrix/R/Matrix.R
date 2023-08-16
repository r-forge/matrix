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
