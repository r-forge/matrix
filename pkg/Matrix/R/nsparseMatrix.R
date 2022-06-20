#### Superclass Methods for all sparse nonzero-pattern matrices

## MJ: no longer needed ... replacement in ./sparseMatrix.R
if(FALSE) {
.C2nC <- function(from, isTri = is(from, "triangularMatrix"))
    .Call(Csparse_to_nz_pattern, from, isTri)

setAs("CsparseMatrix", "nsparseMatrix", function(from) .C2nC(from))
setAs("CsparseMatrix", "nMatrix",       function(from) .C2nC(from))

setAs("nsparseMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))
} ## MJ

setMethod("is.na", signature(x = "nsparseMatrix"), is.na_nsp)

setMethod("image", "nsparseMatrix",
          function(x, ...) image(as(x, "dMatrix"), ...))
