setMethod("t", signature(x = "tscMatrix"),
          function(x) .Call("tsc_transpose", x),
          valueClass = "tscMatrix")

setAs("tscMatrix", "dgTMatrix",
      function(from) .Call("tsc_to_dgTMatrix", from))

setAs("tscMatrix", "dgeMatrix",
      function(from) as(as(from, "dgTMatrix"), "dgeMatrix"))
