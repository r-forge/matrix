setReplaceMethod("coef", signature(object = "ssclme", value = "numeric"),
                 function(object, unconst = FALSE, ..., value)
                 .Call(ifelse(unconst, "ssclme_coefGetsUnc",
                              "ssclme_coefGets"),
                       object, as.double(value), PACKAGE = "Matrix"))

setAs("ssclme", "sscMatrix",
      function(from)
      .Call("ssclme_asSscMatrix", from, PACKAGE = "Matrix"))

setAs("ssclme", "tscMatrix",
      function(from)
      .Call("ssclme_asTscMatrix", from, PACKAGE = "Matrix"))

