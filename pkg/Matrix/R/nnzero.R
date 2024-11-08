## METHODS FOR GENERIC: nnzero
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## na.counted:
## FALSE ... NA is treated as    zero and so excluded from count
##  TRUE ... NA is treated as nonzero and so included   in count
##    NA ... NA is indeterminate (zero or nonzero) hence count is NA

.nnzero <-
function(x, na.counted = TRUE, nnzmax = length(x))
    .Call(R_nnz, x, na.counted, nnzmax)

.sparseDefault <-
function(x, na.counted = TRUE, tol = 0.5)
    .nnzero(x, na.counted = na.counted) < tol * length(x)

 sparseDefault <-
function(x, na.counted = TRUE, tol = 0.5)
     nnzero(x, na.counted = na.counted) < tol * length(x)

setMethod("nnzero", c(x = "ANY"),
          function(x, na.counted = NA, ...)
              switch(typeof(x),
                     logical =, integer =, double =, complex =
                         .nnzero(x, na.counted),
                     ## Else hope that methods exist for 'is.na', '!=' :
                     sum(if(is.na(na.counted))
                             x != 0
                         else if(na.counted)
                             is.na(x) | x != 0
                         else !is.na(x) & x != 0)))

setMethod("nnzero", c(x = "denseMatrix"),
          function(x, na.counted = NA, ...) {
              d <- x@Dim
              if(any(d == 0L))
                  return(0L)
              if(.M.kind(x) == "n")
                  na.counted <- TRUE
              if((shape <- .M.shape(x)) != "g")
                  x <- .M2packed(x)
              N <- .nnzero(x@x, na.counted)
              switch(shape,
                     "g" = N,
                     "s" = N + N - .nnzero(diag(x, names = FALSE), na.counted),
                     "t" = if(x@diag == "N") N else N + d[1L] - .nnzero(x@x[indDiag(d[1L], upper = x@uplo == "U", packed = TRUE)], na.counted))
          })

setMethod("nnzero", c(x = "sparseMatrix"),
          function(x, na.counted = NA, ...) {
              d <- x@Dim
              if(any(d == 0L))
                  return(0L)
              N <- switch(.M.repr(x),
                          "C" = x@p[d[2L]+1L],
                          "R" = x@p[d[1L]+1L],
                          "T" = length((x <- aggregateT(x))@i))
              if(.M.kind(x) != "n")
                  N <- .nnzero(x@x, na.counted, N)
              switch(.M.shape(x),
                     "g" = N,
                     "s" = N + N - .nnzero(diag(x, names = FALSE), na.counted),
                     "t" = if(x@diag == "N") N else N + d[1L])
          })

setMethod("nnzero", c(x = "diagonalMatrix"),
          function(x, na.counted = NA, ...) {
              if(x@diag != "N")
                  x@Dim[1L]
              else {
                  if(.M.kind(x) == "n")
                      na.counted <- TRUE
                  .nnzero(x@x, na.counted)
              }
          })

setMethod("nnzero", c(x = "indMatrix"),
          function(x, ...)
              length(x@perm))

setMethod("nnzero", c(x = "sparseVector"),
          function(x, na.counted = NA, ...)
              if(.M.kind(x) == "n")
                  length(x@i)
              else .nnzero(x@x, na.counted))

setMethod("nnzero", c(x = "sparseCholesky"),
          function(x, ...)
              nnzero(as(x, "CsparseMatrix"), ...))
