#### Triangular Packed Matrices -- Coercion and Methods

setAs("dtpMatrix", "dtrMatrix",
      function(from) .Call("dtpMatrix_as_dtrMatrix", from) )

setAs("dtpMatrix", "dgeMatrix",
      function(from) as(as(from, "dtrMatrix"), "dgeMatrix"))

setAs("dtpMatrix", "matrix",
      function(from) as(as(from, "dtrMatrix"), "matrix"))

setMethod("%*%", signature(x = "dtpMatrix", y = "dgeMatrix"),
	  function(x, y) .Call("dtpMatrix_dgeMatrix_mm", x, y))

setMethod("determinant", signature(x = "dtpMatrix", logarithm = "missing"),
	  function(x, logarithm, ...) determinant(x, TRUE))

setMethod("diag", signature(x = "dtpMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("dtpMatrix_getDiag", x),
          valueClass = "numeric")
          
setMethod("determinant", signature(x = "dtpMatrix", logarithm = "logical"),
	  function(x, logarithm, ...) {
	      dg <- diag(x)
	      if (logarithm) {
		  modulus <- sum(log(abs(dg)))
		  sgn <- prod(sign(dg))
	      } else {
		  modulus <- prod(dg)
		  sgn <- sign(modulus)
		  modulus <- abs(modulus)
	      }
	      attr(modulus, "logarithm") <- logarithm
	      val <- list(modulus = modulus, sign = sgn)
	      class(val) <- "det"
	      val
	  })

setMethod("norm", signature(x = "dtpMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dtpMatrix_norm", x, type),
	  valueClass = "numeric")

setMethod("norm", signature(x = "dtpMatrix", type = "missing"),
	  function(x, type, ...)
	  .Call("dtpMatrix_norm", x, "O"),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtpMatrix", type = "character"),
	  function(x, type, ...)
	  .Call("dtpMatrix_rcond", x, type),
	  valueClass = "numeric")

setMethod("rcond", signature(x = "dtpMatrix", type = "missing"),
	  function(x, type, ...)
	  .Call("dtpMatrix_rcond", x, "O"),
	  valueClass = "numeric")

setMethod("solve", signature(a = "dtpMatrix", b="missing"),
	  function(a, b, ...)
	  .Call("dtpMatrix_solve", a),
	  valueClass = "dtpMatrix")

setMethod("solve", signature(a = "dtpMatrix", b="matrix"),
	  function(a, b, ...)
	  .Call("dtpMatrix_matrix_solve", a, b),
	  valueClass = "matrix")

setMethod("t", signature(x = "dtpMatrix"),
          function(x) as(t(as(x, "dtrMatrix")), "dtpMatrix"),
          valueClass = "dtpMatrix")
###
