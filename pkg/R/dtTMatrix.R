setAs("dtTMatrix", "dtrMatrix",
      function(from) .Call("dtTMatrix_as_dtrMatrix", from))

setAs("dtTMatrix", "matrix",
      function(from) as(as(from, "dtrMatrix"), "matrix"))
      
setMethod("t", signature(x = "dtTMatrix"),
          function(x)
          new("dtTMatrix", Dim = rev(x@Dim), diag = x@diag,
              i = x@j, j = x@i, x = x@x,
              uplo = if (x@uplo == "U") "L" else "U"),
          valueClass = "dtTMatrix")
