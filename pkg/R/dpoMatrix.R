setMethod("chol", signature(x = "poMatrix"),
          function(x, pivot, LINPACK)
          .Call("poMatrix_chol", x))

setMethod("rcond", signature(x = "poMatrix", type = "character"),
          function(x, type, ...)
          .Call("poMatrix_rcond", x, type),
          valueClass = "numeric")

setMethod("rcond", signature(x = "poMatrix", type = "missing"),
          function(x, type, ...)
          .Call("poMatrix_rcond", x, "O"),
          valueClass = "numeric")

setMethod("solve", signature(a = "poMatrix", b = "missing"),
          function(a, b, ...)
          .Call("poMatrix_solve", a),
          valueClass = "poMatrix")

setMethod("solve", signature(a = "poMatrix", b = "dgeMatrix"),
          function(a, b, ...)
          .Call("poMatrix_dgeMatrix_solve", a, b),
          valueClass = "dgeMatrix")

setMethod("solve", signature(a = "poMatrix", b = "matrix"),
          function(a, b, ...)
          .Call("poMatrix_matrix_solve", a, b),
          valueClass = "matrix")

setMethod("solve", signature(a = "poMatrix", b = "numeric"),
          function(a, b, ...)
          as.numeric(.Call("poMatrix_matrix_solve",
                           a, as.matrix(b))),
          valueClass = "numeric")
