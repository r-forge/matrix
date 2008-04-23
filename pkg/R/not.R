#### --- All method definitions for  "!" (not) ---

## Divert everything to  "lMatrix" and its subclasses :
setMethod("!", "Matrix", function(x) !as(x, "lMatrix"))

## -- diag ---

setMethod("!", "ldiMatrix", function(x) {
    if(x@diag == "N")
	x@x <- !x@x
    else { ## "U"
	x@diag <- "N"
	x@x <- rep.int(FALSE, x@Dim[1])
    }
    x
})

## -- lsparse --

setMethod("!", "lsparseMatrix",
          ## turns FALSE to TRUE --> dense matrix
          function(x) !as(x, "denseMatrix"))# was "lgeMatrix"

## Use "Matrix" method !as(. , "lMatrix")
## setMethod("!", "nsparseMatrix",
##           ## turns FALSE to TRUE --> dense matrix
##           function(x) !as(x, "ngeMatrix"))


## -- ldense ---

setMethod("!", "ltrMatrix",
	  function(x) {
	      x@x <- !x@x ## And now fill one triangle with '!FALSE' results :
	      ## TODO: the following should be .Call using
	      ##	a variation of make_array_triangular:
	      r <- as(x, "lgeMatrix")
	      n <- x@Dim[1]
	      if(x@diag == "U")
		  r@x[indDiag(n)] <- FALSE ## result has diagonal all FALSE
	      r@x[indTri(n, upper=x@uplo != "U")] <- TRUE
	      r
	  })

setMethod("!", "ltpMatrix", function(x) !as(x, "ltrMatrix"))

## for the other ldense* ones
setMethod("!", "lgeMatrix",
	  function(x) { x@x <- !x@x ; x })
## FIXME : this loses symmetry "lsy" and "lsp":
setMethod("!", "ldenseMatrix",
	  function(x) !as(x, "lgeMatrix"))

## -- ndense ---

setMethod("!", "ntrMatrix",
	  function(x) {
	      x@x <- !x@x
	      ## And now we must fill one triangle with '!FALSE' results :

	      ## TODO: the following should be .Call using
	      ##	a variation of make_array_triangular:
	      r <- as(x, "ngeMatrix")
	      n <- x@Dim[1]
	      coli <- rep(1:n, each=n)
	      rowi <- rep(1:n, n)
	      Udiag <- x@diag == "U"
	      log.i <-
		  if(x@uplo == "U") {
		      if(Udiag) rowi >= coli else rowi > coli
		  } else {
		      if(Udiag) rowi <= coli else rowi < coli
		  }
	      r@x[log.i] <- TRUE
	      r
	  })

setMethod("!", "ntpMatrix", function(x) !as(x, "ntrMatrix"))

## for the other ldense* ones
setMethod("!", "ngeMatrix",
          function(x) { x@x <- !x@x ; x })
## FIXME : this loses symmetry "nsy" and "nsp":
setMethod("!", "ndenseMatrix",
          function(x) !as(x, "ngeMatrix"))
