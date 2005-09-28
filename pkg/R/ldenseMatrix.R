l2d_Matrix <- function(from) {
    stopifnot(is(from, "lMatrix"))
    newCl <- sub("^l", "d", class(from))
    r <- new(newCl, x = as.double(from@x),
             Dim = from@Dim, Dimnames = from@Dimnames,
             factors = list()) ## FIXME: treat 'factors' smartly
    if(is(r, "triangularMatrix")) {
        r@uplo <- from@uplo
        r@diag <- from@diag
    } else if(is(r, "symmetricMatrix")) {
        r@uplo <- from@uplo
    }
    r
}

dummy_meth <- function(x) {
    cl <- class(x)
    as(callGeneric(as(x, sub("^l", "d", cl))), cl)
}

setAs("lgeMatrix", "dgeMatrix", l2d_Matrix)
setAs("ltrMatrix", "dtrMatrix", l2d_Matrix)
setAs("ltpMatrix", "dtpMatrix", l2d_Matrix)
setAs("lsyMatrix", "dsyMatrix", l2d_Matrix)
setAs("lspMatrix", "dspMatrix", l2d_Matrix)

setAs("lspMatrix", "lsyMatrix",
      function(from) .Call("lspMatrix_as_lsyMatrix", from) )
setAs("lsyMatrix", "lspMatrix",
      function(from) .Call("lsyMatrix_as_lspMatrix", from) )

setAs("ltpMatrix", "ltrMatrix",
      function(from) .Call("ltpMatrix_as_ltrMatrix", from) )
setAs("ltrMatrix", "ltpMatrix",
      function(from) .Call("ltrMatrix_as_ltpMatrix", from) )

setAs("ldenseMatrix", "matrix",
      function(from) as(as(from, sub("^l", "d", class(from))), "matrix"))

setAs("matrix", "ldenseMatrix",
      function(from) callGeneric(as(from, "lgeMatrix")))

setMethod("t", signature(x = "lgeMatrix"), t_geMatrix)
setMethod("t", signature(x = "ltrMatrix"), t_trMatrix)
setMethod("t", signature(x = "lsyMatrix"), t_trMatrix)
setMethod("t", signature(x = "ltpMatrix"),
          function(x) as(callGeneric(as(x, "ltrMatrix")), "ltpMatrix"))
setMethod("t", signature(x = "lspMatrix"),
          function(x) as(callGeneric(as(x, "lsyMatrix")), "lspMatrix"))
