library(Matrix)

### Do all kinds of object creation and coercion


chk.matrix <- function(M) {
    ## check if "matrix" coercion looks ok
    cl <- class(M)
    d <- dim(M)
    cat("class ", dQuote(cl), " [",d[1]," x ",d[2],"]; slots (",
        paste(slotNames(M), collapse=","), ")\n", sep='')
    stopifnot(validObject(M),
              identical(dim(m <- as(M,"matrix")), dim(M))
              )
}

## "dge"
chk.matrix(m1 <- Matrix(1:6, ncol=2))
chk.matrix(m2 <- Matrix(1:7, ncol=3)) # a warning
## "dpo"
chk.matrix(cm <- crossprod(m1))
chk.matrix(as(cm, "dsyMatrix"))
chk.matrix(as(cm, "dgeMatrix"))
try( chk.matrix(as(cm, "Matrix")) ) # gives an error

## Cholesky
chk.matrix(ch <- chol(cm))
## FIXME:
try( chk.matrix(ch2 <- chol(as(cm, "dsyMatrix"))) ) # should not give an error
try( chk.matrix(ch3 <- chol(as(cm, "dgeMatrix"))) ) # nor that one
