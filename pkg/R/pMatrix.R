#### Permutation Matrices -- Coercion and Methods

setAs("integer", "pMatrix",
      function(from) {
          n <- length(from)
          new("pMatrix", Dim = rep.int(n, 2), perm = from)
      })

setAs("pMatrix", "matrix",
      function(from) {
          fp <- from@perm
          diag(nrow = length(fp))[fp , ]
      })

setMethod("solve", signature(a = "pMatrix", b = "missing"),
          function(a, b) {
              bp <- ap <- a@perm
              bp[ap] <- seq(along = ap)
              new("pMatrix", Dim = a@Dim, perm = bp)
          }, valueClass = "pMatrix")

setMethod("t", signature(x = "pMatrix"), function(x) solve(x))

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
	  function(x, y) x[ , y@perm], valueClass = "matrix")

setMethod("%*%", signature(x = "pMatrix", y = "matrix"),
	  function(x, y) y[x@perm ,], valueClass = "matrix")

## the following methods can be rewritten when "[" methods for
## dgeMatrix are available  

setMethod("%*%", signature(x = "dgeMatrix", y = "pMatrix"),
	  function(x, y) as(callGeneric(x, as(y, "matrix")), "dgeMatrix"),
          valueClass = "dgeMatrix")

setMethod("%*%", signature(x = "pMatrix", y = "dgeMatrix"),
          function(x, y) as(callGeneric(as(x, "matrix"), y), "dgeMatrix"),
          valueClass = "dgeMatrix")
