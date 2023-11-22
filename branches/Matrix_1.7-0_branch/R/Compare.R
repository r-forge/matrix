## METHODS FOR GENERIC: Compare (group)                 WORK IN PROGRESS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Compare")
## [1] "==" ">"  "<"  "!=" "<=" ">="

.Compare.swap <- function(e1, e2)
    switch(.Generic,
           "==" = e2 == e1, "<"  = e2 > e1, "<=" = e2 >= e1,
           "!=" = e2 != e1, ">"  = e2 < e1, ">=" = e2 <= e1,
           stop(gettextf("unexpected %s=\"%s\" in '%s' method",
                         ".Generic", .Generic, "Compare"),
                domain = NA))


## .... denseMatrix ....................................................

setMethod("Compare", signature(e1 = "denseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "denseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "denseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "denseMatrix", e2 = "vector"),
          function(e1, e2) {

          })


## .... sparseMatrix ...................................................

setMethod("Compare", signature(e1 = "sparseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "sparseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "sparseMatrix", e2 = "vector"),
          function(e1, e2) {

          })


## .... sparseVector ...................................................

setMethod("Compare", signature(e1 = "sparseVector", e2 = "Matrix"),
          .Compare.swap)

setMethod("Compare", signature(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Compare", signature(e1 = "sparseVector", e2 = "vector"),
          function(e1, e2) {

          })


## .... vector .........................................................

setMethod("Compare", signature(e1 = "vector", e2 = "Matrix"),
          .Compare.swap)

setMethod("Compare", signature(e1 = "vector", e2 = "sparseVector"),
          .Compare.swap)
