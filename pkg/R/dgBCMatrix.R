setAs("dgBCMatrix", "dgCMatrix",
      function(from) .Call("dgBCMatrix_to_dgCMatrix", from))

setAs("dgBCMatrix", "dgTMatrix",
      function(from) as(as(from, "dgCMatrix"), "dgTMatrix"))
