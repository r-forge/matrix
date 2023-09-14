## METHODS FOR GENERIC: Arith (group)                   WORK IN PROGRESS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Arith")
## [1] "+"   "-"   "*"   "^"   "%%"  "%/%" "/"


## .... unary ..........................................................

setMethod("+", signature(e1 = "Matrix", e2 = "missing"),
          function(e1, e2)
              if(any(.M.kind(e1) == c("n", "l"))) .M2kind(e1, "d") else e1)

setMethod("+", signature(e1 = "sparseVector", e2 = "missing"),
          function(e1, e2)
              if(any(.V.kind(e1) == c("n", "l"))) .V2kind(e1, "i") else e1)

for(.cl in c("generalMatrix", "symmetricMatrix"))
setMethod("-", signature(e1 = .cl, e2 = "missing"),
          function(e1, e2) {
              if(any(.M.kind(e1) == c("n", "l")))
                  e1 <- .M2kind(e1, "d")
              if(length(e1@factors) > 0L)
                  e1@factors <- list()
              e1@x <- -e1@x
              e1
          })

for(.cl in c("triangularMatrix", "diagonalMatrix"))
setMethod("-", signature(e1 = .cl, e2 = "missing"),
          function(e1, e2) {
              if(any(.M.kind(e1) == c("n", "l")))
                  e1 <- .M2kind(e1, "d")
              if(e1@diag != "N")
                  diag(e1) <- TRUE
              e1@x <- -e1@x
              e1
          })

setMethod("-", signature(e1 = "sparseVector", e2 = "missing"),
          function(e1, e2) {
              if(any(.V.kind(e1) == c("n", "l")))
                  e1 <- .V2kind(e1, "i")
              e1@x <- -e1@x
              e1
          })


## .... denseMatrix ....................................................

setMethod("Arith", signature(e1 = "denseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "denseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "denseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "denseMatrix", e2 = "vector"),
          function(e1, e2) {
              if(all(typeof(e2) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e2), domain = NA)

          })


## .... sparseMatrix ...................................................

setMethod("Arith", signature(e1 = "sparseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "sparseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "sparseMatrix", e2 = "vector"),
          function(e1, e2) {
              if(all(typeof(e2) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e2), domain = NA)

          })


## .... sparseVector ...................................................

setMethod("Arith", signature(e1 = "sparseVector", e2 = "Matrix"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Arith", signature(e1 = "sparseVector", e2 = "vector"),
          function(e1, e2) {
              if(all(typeof(e2) != c("logical", "integer", "double", "complex")))
                  stop(.Ops.invalid(e2), domain = NA)

          })


## .... vector .........................................................

setMethod("Arith", signature(e1 = "vector", e2 = "denseMatrix"),
          function(e1, e2) {
              if(all(typeof(e1) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e1), domain = NA)

          })

setMethod("Arith", signature(e1 = "vector", e2 = "sparseMatrix"),
          function(e1, e2) {
              if(all(typeof(e1) != c("logical", "integer", "double")))
                  stop(.Ops.invalid(e1), domain = NA)

          })

setMethod("Arith", signature(e1 = "vector", e2 = "sparseVector"),
          function(e1, e2) {
              if(all(typeof(e1) != c("logical", "integer", "double", "complex")))
                  stop(.Ops.invalid(e1), domain = NA)

          })

rm(.cl)
