## METHODS FOR GENERIC: updown
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("updown",
          c(update = "character", C = "ANY", L = "ANY"),
          function(update, C, L)
              updown(identical(update, "+"), C, L))

setMethod("updown",
          c(update = "logical", C = "Matrix", L = "sparseCholesky"),
          function(update, C, L)
              updown(update, .M2kind(.M2C(C), ","), L))

setMethod("updown",
          c(update = "logical", C = "matrix", L = "sparseCholesky"),
          function(update, C, L)
              updown(update, .m2sparse(C, ",gC"), L))

for (.cl in c("dgCMatrix", "dsCMatrix"))
setMethod("updown",
          c(update = "logical", C = .cl, L = "sparseCholesky"),
          function(update, C, L) {
              if (length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .Call(sparseCholesky_updown, L, C, update)
          })

setMethod("updown",
          c(update = "logical", C = "dtCMatrix", L = "sparseCholesky"),
          function(update, C, L) {
              if (C@diag != "N")
                  C <- ..diagU2N(C)
              if (length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .Call(sparseCholesky_updown, L, C, update)
          })

rm(.cl)


## METHODS FOR GENERIC: update
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Exported:
.updateCHMfactor <-
function(object, parent, mult = 0)
    .Call(sparseCholesky_update, object, parent, mult)

setMethod("update", c(object = "sparseCholesky"),
          function(object, parent, mult = 0, ...) {
              parent <- .M2kind(.M2C(parent), ",")
              if ((shape <- .M.shape(parent)) != "s") {
                  Matrix.message(gettextf("'%1$s' is not formally symmetric; factorizing tcrossprod(%1$s)",
                                          "parent"),
                                 domain = NA)
                  if (shape == "t" && parent@diag != "N")
                      parent <- ..diagU2N(parent)
              }
              .Call(sparseCholesky_update, object, parent, mult)
          })
