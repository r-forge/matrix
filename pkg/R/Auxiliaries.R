#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

## For %*% (M = Matrix; n = numeric):
.M.n <- function(x, y) callGeneric(x, as.matrix(y))
.n.M <- function(x, y) callGeneric(rbind(x), y)

## chol() via "dpoMatrix"
cholMat <- function(x, pivot, LINPACK) {
    px <- as(x, "dpoMatrix")
    if(identical(TRUE, validObject(px, test=TRUE)))
        chol(px)
    else stop("'x' is not positive definite -- chol() undefined.")
}
