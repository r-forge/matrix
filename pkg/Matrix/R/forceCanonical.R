## METHODS FOR GENERIC: forceCanonical
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.forceCanonical <-
function(x, check = TRUE)
    switch(.M.repr(x),
           "n" =, "p" =
               .Call(R_dense_force_canonical, x, check),
           "C" =, "R" =, "T" =
               .Call(R_sparse_force_canonical, x, check),
           stop("should never happen ..."))

setMethod("forceCanonical", c(x = "denseMatrix"),
          function(x, check = TRUE, ...)
              .Call(R_dense_force_canonical, x, check))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("forceCanonical", c(x = .cl),
          function(x, check = TRUE, ...)
              .Call(R_sparse_force_canonical, x, check))

setMethod("forceCanonical", c(x = "diagonalMatrix"),
          function(x, ...) {
              diag. <- x@diag
              x. <- x@x
              kind. <- .M.kind(x)
              if (diag. == "N" && (kind. != "n" || !anyNA(x.)))
                  return(x)
              r <- new(paste0(kind., "diMatrix"))
              r@Dim <- d <- x@Dim
              r@x <-
                  if (diag. == "N")
                      x. | is.na(x.)
                  else rep.int(switch(typeof(x.), logical = TRUE, integer = 1L, double = 1, complex = 1+0i), d[1L])
              r
          })

setMethod("forceCanonical", c(x = "indMatrix"),
          function(x, ...)
              x)

setMethod("forceCanonical", c(x = "sparseVector"),
          function(x, ...) {
              length. <- x@length
              i. <- x@i
              if (is.integer(length.) && is.integer(i.))
                  return(x)
              kind. <- .M.kind(x)
              r <- new(paste0(kind., "sparseVector"))
              r@length <-
                  if (is.integer(length.))
                      length.
                  else if (length. - 1 < .Machine[["integer.max"]])
                      as.integer(length.)
                  else trunc(length.)
              r@i <-
                  if (is.integer(i.))
                      i.
                  else if (length(i.) == 0L || i.[length(i.)] - 1 < .Machine[["integer.max"]])
                      as.integer(i.)
                  else trunc(i.)
              if (kind. != "n")
                  r@x <- x@x
              r
          })

rm(.cl)
