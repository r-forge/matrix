### Coercion and Methods for Triangular Triplet Matrices

setAs("dtTMatrix", "dtCMatrix",
      function(from) .Call("dtTMatrix_as_dtCMatrix", from))

setAs("dtTMatrix", "dgTMatrix",
      function(from) {
          d <- from@Dim
          if(uDiag <- from@diag == "U") # unit diagonal, need to add '1's
              uDiag <- (n <- d[1]) > 0
          new("dgTMatrix", Dim = d, Dimnames = from@Dimnames,
              i = c(from@i,from@j, if(uDiag) 1:n),
              j = c(from@j,from@i, if(uDiag) 1:n),
              x = c(from@x,from@x, if(uDiag) rep.int(1,n)))
      })

## Conversion to dense storage is first to a dtrMatrix
setAs("dtTMatrix", "dtrMatrix",
      function(from) .Call("dtTMatrix_as_dtrMatrix", from))

setAs("dtTMatrix", "matrix",
      function(from) as(as(from, "dtrMatrix"), "matrix"))

setAs("dtTMatrix", "dgeMatrix",
      function(from) as(as(from, "dtrMatrix"), "dgeMatrix"))

setMethod("t", signature(x = "dtTMatrix"),
          function(x)
          new("dtTMatrix", Dim = rev(x@Dim), diag = x@diag,
              i = x@j, j = x@i, x = x@x,
              uplo = if (x@uplo == "U") "L" else "U"),
          valueClass = "dtTMatrix")
