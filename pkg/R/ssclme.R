setReplaceMethod("coef",
                 signature(object = "ssclme", value = "numeric"),
                 function(object, value) 
                 .Call("ssclme_coefGets", object, value, PACKAGE = "Matrix"))

setAs("ssclme", "sscMatrix",
      function(from)
      .Call("ssclme_asSscMatrix", from, PACKAGE = "Matrix"))

setAs("ssclme", "tscMatrix",
      function(from)
      .Call("ssclme_asTscMatrix", from, PACKAGE = "Matrix"))

