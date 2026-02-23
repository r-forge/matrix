## METHODS FOR GENERIC: zapsmall
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: not quite optimal as we lose diag="U" even if 1 is not zapped

setMethod("zapsmall", signature(x = "Matrix"),
          function(x, digits = getOption("digits"),
                   mFUN = function(x, ina) max(abs(x[!ina])),
                   min.d = 0L, ...) {
              r <- .M2kind(forceCanonical(.promote(x), diag = "N"), ",")
              r@x <- zapsmall(x@x, digits = digits, mFUN = mFUN, min.d = min.d)
              if (.M.shape(r) != "t")
                  r@factors <- list()
              r
          })

setMethod("zapsmall", signature(x = "sparseVector"),
          function(x, digits = getOption("digits"),
                   mFUN = function(x, ina) max(abs(x[!ina])),
                   min.d = 0L, ...) {
              r <- .V2kind(.promote(x), ",")
              r@x <- zapsmall(x@x, digits = digits, mFUN = mFUN, min.d = min.d)
              r
          })
