## METHODS FOR GENERIC: diff
## NB: mostly cut and paste of base::diff.default
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("denseMatrix", "sparseMatrix"))
setMethod("diff", c(x = .cl),
          function(x, lag = 1L, differences = 1L, ...) {
              if (length(lag) != 1L || length(differences) != 1L ||
                  lag < 1L || differences < 1L)
                  stop(gettextf("'%s' and '%s' must be positive integers",
                                "lag", "differences"),
                       domain = NA)
              if (lag * differences >= x@Dim[1L])
                  return(x[0L])
              i1 <- -seq_len(lag)
              for (i in seq_len(differences)) {
                  m <- x@Dim[1L]
                  x <- x[i1, , drop = FALSE] -
                      x[-m:-(m - lag + 1L), , drop = FALSE]
              }
              x
          })

setMethod("diff", c(x = "sparseVector"),
          function(x, lag = 1L, differences = 1L, ...) {
              if (length(lag) != 1L || length(differences) != 1L ||
                  lag < 1L || differences < 1L)
                  stop(gettextf("'%s' and '%s' must be positive integers",
                                "lag", "differences"),
                       domain = NA)
              if (lag * differences >= length(x))
                  return(x[0L])
              i1 <- -seq_len(lag)
              for (i in seq_len(differences))
                  x <- x[i1] - x[-length(x):-(length(x) - lag + 1L)]
              x
          })

rm(.cl)
