## METHODS FOR GENERIC: updown
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("updown",
          c(update = "character", C = "ANY", L = "ANY"),
          function(update, C, L)
              updown(identical(update, "+"), C, L))

setMethod("updown",
          c(update = "logical", C = "Matrix", L = "sparseCholesky"),
          function(update, C, L) {
              C <- .M2gen(.M2C(C), ",")
              if (length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .Call(sparseCholesky_updown, L, C, update)
          })

setMethod("updown",
          c(update = "logical", C = "matrix", L = "sparseCholesky"),
          function(update, C, L)
              updown(update, .m2sparse(C, ",gC"), L))


## METHODS FOR GENERIC: update
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Exported:
.updateCHMfactor <-
function(object, parent, mult = 0)
    .Call(sparseCholesky_update, object, parent, mult)

setMethod("update", c(object = "sparseCholesky"),
          function(object, parent, mult = 0, ...) {
              parent <- .M2kind(.M2C(parent), ",")
              if (.M.shape(parent) == "t" && parent@diag != "N")
                  parent <- ..diagU2N(parent)
              .Call(sparseCholesky_update, object, parent, mult)
          })
