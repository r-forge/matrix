library(Matrix)

data(mm)
stopifnot(is(mm) == c("dgCMatrix", "dMatrix", "Matrix"),
          dim(mm) == c(1850, 712),
          identical(dimnames(mm), list(NULL,NULL)))
str(mm)
tmm <- t(mm)
str(tmm)
stopifnot(validObject(tmm),
          identical(as(tmm, "matrix"), t(as(mm, "matrix"))))
