setMethod("eigen", signature(x = "dgeMatrix", only.values = "missing"),
          function(x, symmetric, only.values, EISPACK) {
              nCall = match.call()
              nCall$only.values = FALSE
              eval(nCall, parent.frame())
          })

setMethod("eigen", signature(x = "dgeMatrix", only.values = "logical"),
          function(x, symmetric, only.values, EISPACK)
          .Call("dgeMatrix_eigen", x, only.values)
          )

setMethod("Schur", signature(x = "dgeMatrix", vectors = "missing"),
          function(x, vectors, ...) {
              nCall = match.call()
              nCall$vectors = FALSE
              eval(nCall, parent.frame())
          })
         
setMethod("Schur", signature(x = "dgeMatrix", vectors = "logical"),
          function(x, vectors, ...)
          .Call("dgeMatrix_Schur", x, vectors)
          )

