library(Matrix)

data(mm)
stopifnot(is(mm) == c("dgCMatrix", "dMatrix", "Matrix"),
          dim(mm) == (dm <- c(1850, 712)),
          identical(dimnames(mm), list(NULL,NULL)))
str(mm)
tmm <- t(mm)
str(tmm)

mTm <- crossprod(mm)
mmT  <- crossprod(tmm) ## should be the same as
if(FALSE) # << not yet available ! {missing  "csc_tcrossprod"}
mmT. <- tcrossprod(mm)
str(mTm)
str(mmT)# much larger (i.e less sparse)
stopifnot(validObject(tmm), dim(tmm) == dm[2:1],
          validObject(mTm), dim(mTm) == dm[c(2,2)],
          validObject(mmT), dim(mmT) == dm[c(1,1)],
          identical(as(tmm, "matrix"), t(as(mm, "matrix"))))
