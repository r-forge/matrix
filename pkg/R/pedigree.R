setAs("pedigree", "dtCMatrix", # representation as T^{-1}
      function(from) {
          sire <- from@sire
          n <- length(sire)
          animal <- seq(along = sire)
          j <- c(as.integer(sire), as.integer(from@dam))
          ind <- !is.na(j)
          as(new("dtTMatrix", i = rep.int(animal, 2)[ind] - 1:1,
                 j = j[ind] - 1:1, x = rep.int(-0.5, sum(ind)),
                 Dim = c(n,n), Dimnames = list(levels(sire), NULL),
                 uplo = "L", diag = "U"), "dtCMatrix")
      })

setMethod("show", signature(object = "pedigree"),
          function(object) print(as(object, "data.frame")))

setAs("pedigree", "data.frame",
      function(from)
      data.frame(sire = from@sire, dam = from@dam,
                 row.names = levels(from@sire)))

setMethod("head", "pedigree", function(x, ...)
          do.call("head", list(x = as(x, "data.frame"), ...)))

setMethod("tail", "pedigree", function(x, ...)
          do.call("tail", list(x = as(x, "data.frame"), ...)))
