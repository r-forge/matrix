setAs("tripletMatrix", "cscMatrix",
      function(from)
      .Call("triplet_to_csc", from, PACKAGE = "Matrix")
      )
