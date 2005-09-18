setMethod("t", signature(x = "dtCMatrix"),
          function(x) .Call("tsc_transpose", x),
          valueClass = "dtCMatrix")

setMethod("image", signature(x = "dtCMatrix"),
          function(x, ...) image(as(x, "dgTMatrix"), ...))

setAs("dtCMatrix", "dgTMatrix",
      function(from) .Call("tsc_to_dgTMatrix", from))

setAs("dtCMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))
