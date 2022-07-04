## METHODS FOR GENERIC: nnzero
## * used to retrieve number of nonzero elements,
##   i.e., number of elements excl. both structural and non-structural zeros
## * like MATLAB's nnz() but more sophisticated due to handling of NA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.nnz <- function(x, countNA = NA, nnzmax = length(x))
    .Call(R_nnz, x, countNA, nnzmax)

## na.counted:
## FALSE ... NA is treated as    zero and so excluded from count
##  TRUE ... NA is treated as nonzero and so included   in count
##    NA ... NA is indeterminate (could be zero or nonzero) hence count is NA

## A fallback for any class with methods for 'is.na' and '!=' ...
setMethod("nnzero", "ANY",
	  function(x, na.counted = NA) {
              sum(if(is.na(na.counted))
                      x != 0
                  else if(na.counted)
                      is.na(x) | x != 0
                  else !is.na(x) & x != 0)
          })

setMethod("nnzero", "vector",
	  function(x, na.counted = NA) .nnz(x, na.counted))

setMethod("nnzero", "array",
	  function(x, na.counted = NA) .nnz(x, na.counted))

setMethod("nnzero", "CHMfactor",
	  function(x, na.counted = NA)
              callGeneric(as(x, "sparseMatrix"), na.counted))

setMethod("nnzero", "diagonalMatrix",
	  function(x, na.counted = NA)
              if(x@diag == "N") .nnz(x@x, na.counted) else x@Dim[1L])

setMethod("nnzero", "indMatrix",
          function(x, na.counted = NA) x@Dim[1L])

setMethod("nnzero", "sparseMatrix",
	  function(x, na.counted = NA) {
	      d <- x@Dim
              if(any(d == 0L))
                  return(0L)
              cld <- getClassDef(class(x))
              N <- if(extends(cld, "CsparseMatrix"))
                       x@p[d[2L]+1L]
                   else if(extends(cld, "RsparseMatrix"))
                       x@p[d[1L]+1L]
                   else length((x <- .Call(Tsparse_aggregate, x))@i)
              if(!extends(cld, "nsparseMatrix"))
                  N <- .nnz(x@x, na.counted, N)
              if(extends(cld, "generalMatrix"))
                  N
              else if(extends(cld, "symmetricMatrix"))
                  N + N - .nnz(diag(x), na.counted)
              else if(x@diag != "N")
                  N + d[1L]
              else N
          })

setMethod("nnzero", "denseMatrix",
	  function(x, na.counted = NA) {
              d <- x@Dim
              if(any(d == 0L))
                  return(0L)
              xx <- x@x
              cld <- getClassDef(class(x))
              if(extends(cld, "ndenseMatrix"))
                  na.counted <- TRUE
              if(extends(cld, "generalMatrix"))
                  return(.nnz(xx, na.counted))
              n <- d[1L]
              upper <- x@uplo == "U"
              if(extends(cld, "unpackedMatrix"))
                  xx <- xx[indTri(n, upper, diag = TRUE, packed = FALSE)]
              N <- .nnz(xx, na.counted)
              if(extends(cld, "symmetricMatrix"))
                  N + N - .nnz(diag(x), na.counted)
              else if(x@diag != "N")
                  N + n - .nnz(xx[indDiag(n, upper, packed = TRUE)], na.counted)
              else N
          })
