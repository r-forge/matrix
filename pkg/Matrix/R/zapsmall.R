## METHODS FOR GENERIC: zapsmall
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: not quite optimal as we lose diag="U" even if 1 is not zapped

setMethod("zapsmall", signature(x = "Matrix"),
          function(x, ...) {
              r <- .M2kind(forceCanonical(.promote(x)), ",")
              r@x <- zapsmall(x@x, ...)
              if (.M.shape(r) != "t")
                  r@factors <- list()
              r
          })

setMethod("zapsmall", signature(x = "sparseVector"),
          function(x, ...) {
              r <- .V2kind(.promote(x), ",")
              r@x <- zapsmall(x@x, ...)
              r
          })
