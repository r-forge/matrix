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

setMethod("%*%", signature(x = "matrix", y = "pMatrix"),
	  function(x, y) x[ , y@perm], valueClass = "matrix")

setMethod("%*%", signature(x = "pMatrix", y = "matrix"),
	  function(x, y) y[x@perm ,], valueClass = "matrix")

setMethod("t", signature(x = "pMatrix"), function(x) solve(x))
