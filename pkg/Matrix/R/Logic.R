## METHODS FOR GENERIC: Logic (group)                   WORK IN PROGRESS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Logic") # excluding unary "!" -> ./not.R
## [1] "&" "|"

.Logic.swap <- function(e1, e2)
    switch(.Generic, "&" = e2 & e1, "|" = e2 | e1)


## .... denseMatrix ....................................................

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "vector"),
          function(e1, e2) {

          })


## .... sparseMatrix ...................................................

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "vector"),
          function(e1, e2) {

          })


## .... sparseVector ...................................................

setMethod("Logic", signature(e1 = "sparseVector", e2 = "Matrix"),
          .Logic.swap)

setMethod("Logic", signature(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseVector", e2 = "vector"),
          function(e1, e2) {

          })


## .... vector .........................................................

setMethod("Logic", signature(e1 = "vector", e2 = "Matrix"),
          .Logic.swap)

setMethod("Logic", signature(e1 = "vector", e2 = "sparseVector"),
          .Logic.swap)
