library(Matrix)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

(t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2,
           Dim= as.integer(c(4,4))))
stopifnot(validObject(t1),
          validObject(t1c <- as(t1, "dtCMatrix")))
assert.EQ.mat(t1, as(t1c, "matrix"))

## from  0-diagonal to unit-diagonal {low-level step}:
tu <- t1 ; tu@diag <- "U"
tu
stopifnot(validObject(cu <- as(tu, "dtCMatrix")),
          validObject(t(cu)),
          validObject(t(tu)))
assert.EQ.mat(cu, as(tu,"matrix"), tol=0)
