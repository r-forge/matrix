setAs("matrix", "lgeMatrix",
      function(from) {
	  new("lgeMatrix",
	      x = as.logical(from),
	      Dim = as.integer(dim(from)),
	      Dimnames =
	      if(!is.null(dn <- dimnames(from))) dn else list(NULL,NULL)
	      )
      })

setAs("lgeMatrix", "matrix",
      function(from) {
	  array(from@x, dim = from@Dim, dimnames = from@Dimnames)
      })


## dense |-> compressed :
setAs("lgeMatrix", "lgTmatrix",
      function(from) {
          ##  cheap but not so efficient:
          ij <- which(as(from,"matrix"), arr.ind = TRUE)
          new("lgTMatrix", i = ij[,1], j = ij[,2],
              Dim = from@Dim, Dimnames = from@Dimnames)
      })
