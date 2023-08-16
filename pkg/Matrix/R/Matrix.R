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

## We want to use all.equal.numeric() *and* make sure that uses
## not just base::as.vector but the generic with our methods:
all.equal_num <- base::all.equal.numeric
##               ^^^^^^<R>/src/library/base/R/all.equal.R
environment(all.equal_num) <- environment() # our namespace

all.equal_Mat <- function(target, current, check.attributes = TRUE,
                          factorsCheck = FALSE, ...)
{
    msg <- attr.all_Mat(target, current, check.attributes=check.attributes,
                        factorsCheck=factorsCheck, ...)
    if(is.list(msg)) msg[[1]]
    else .a.e.comb(msg,
                   all.equal_num(as.vector(target), as.vector(current),
                                 check.attributes=check.attributes, ...))
}

## The all.equal() methods for dense matrices (and fallback):
setMethod("all.equal", c(target = "Matrix", current = "Matrix"),
          all.equal_Mat)
setMethod("all.equal", c(target = "Matrix", current = "ANY"),
          all.equal_Mat)
setMethod("all.equal", c(target = "ANY", current = "Matrix"),
          all.equal_Mat)
rm(all.equal_Mat)
## -> ./sparseMatrix.R, ./sparseVector.R  have specific methods
