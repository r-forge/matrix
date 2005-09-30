setMethod("t", signature(x = "dtCMatrix"),
          function(x) .Call("tsc_transpose", x),
          valueClass = "dtCMatrix")

setAs("dtCMatrix", "ltCMatrix",
      function(from) new("ltCMatrix", i = from@i, p = from@p,
                         uplo = from@uplo, diag = from@diag,
                         Dim = from@Dim, Dimnames = from@Dimnames))

setAs("dtCMatrix", "dgTMatrix",
      function(from) .Call("tsc_to_dgTMatrix", from))

setAs("dtCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))
