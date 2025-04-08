## METHODS FOR GENERIC: Arith (group)                   WORK IN PROGRESS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Arith")
## [1] "+"   "-"   "*"   "^"   "%%"  "%/%" "/"


## .... denseMatrix ....................................................

setMethod("Arith", c(e1 = "denseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "denseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "denseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "denseMatrix", e2 = "vector"),
          function(e1, e2) {
              if(all(typeof(e2) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e2), domain = NA)

          })


## .... sparseMatrix ...................................................

setMethod("Arith", c(e1 = "sparseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "sparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "sparseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "sparseMatrix", e2 = "vector"),
          function(e1, e2) {
              if(all(typeof(e2) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e2), domain = NA)

          })


## .... sparseVector ...................................................

setMethod("Arith", c(e1 = "sparseVector", e2 = "Matrix"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Arith", c(e1 = "sparseVector", e2 = "vector"),
          function(e1, e2) {
              if(all(typeof(e2) != c("logical", "integer", "double", "complex")))
                  stop(.Ops.invalid(e2), domain = NA)

          })


## .... vector .........................................................

setMethod("Arith", c(e1 = "vector", e2 = "denseMatrix"),
          function(e1, e2) {
              if(all(typeof(e1) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e1), domain = NA)

          })

setMethod("Arith", c(e1 = "vector", e2 = "sparseMatrix"),
          function(e1, e2) {
              if(all(typeof(e1) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e1), domain = NA)

          })

setMethod("Arith", c(e1 = "vector", e2 = "sparseVector"),
          function(e1, e2) {
              if(all(typeof(e1) != c("logical", "integer", "double", "complex")))
                  stop(.Ops.invalid(e1), domain = NA)

          })

rm(.cl)
