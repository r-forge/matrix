library(Matrix)

## This is example(sp....) -- much extended

mEQ <- function(x,y) {
    ## first drop columns from y  which are all 0 :
    if(any(i0 <- colSums(abs(y)) == 0))
        y <- y[, !i0, drop=FALSE]
    isTRUE(all.equal(x,y, tol=0))
}

dd <- data.frame(a = gl(3,4), b = gl(4,1,12))# balanced 2-way
options("contrasts") # the default:  "contr.treatment"
sparse.model.matrix(~ a + b, dd)
sparse.model.matrix(~ a + b, dd, contrasts = list(a="contr.sum"))
sparse.model.matrix(~ a + b, dd, contrasts = list(b="contr.SAS"))

     ## Sparse method is equivalent to the traditional one :
stopifnot(mEQ(sparse.model.matrix(~ a + b, dd),
              Matrix(model.matrix(~ a + b, dd), sparse=TRUE)),
          mEQ(sparse.model.matrix(~ a + b, dd, contrasts = list(a="contr.sum")),
              Matrix(model.matrix(~ a + b, dd, contrasts = list(a="contr.sum")),
                     sparse=TRUE)),
          mEQ(sparse.model.matrix(~ a + b, dd, contrasts = list(a="contr.SAS")),
              Matrix(model.matrix(~ a + b, dd, contrasts = list(a="contr.SAS")),
                     sparse=TRUE))
          )

(dd3 <- cbind(dd, c = gl(2,6), d = gl(3,8)))
dd. <- dd3[- c(1, 13:15, 17), ]
##
sm <- sparse.model.matrix(~ a + b + c + d, dd.)
mm <- Matrix(model.matrix(~ a + b + c + d, dd.), sparse=TRUE)
stopifnot(mEQ(sm, mm)) ## ok

dim(mm <- Matrix(model.matrix(~ a + b:c + c + d, dd.), sparse=TRUE))
dim(sm <- sparse.model.matrix(~ a + b:c + c + d, dd.)) # wrong: 9 col. vs 12
stopifnot(mEQ(sm, mm))

set.seed(17)
dd4 <- cbind(dd, c = gl(2,6), d = gl(8,3))
dd4 <- cbind(dd4, x = round(rnorm(nrow(dd4)), 1))
dd4 <- dd4[- c(1, 13:15, 17), ]
##-> 'd' has unused levels
##
dim(mm <- Matrix(model.matrix(~ a + b + c + d, dd4), sparse=TRUE))
dim(sm <- sparse.model.matrix(~ a + b + c + d, dd4))
## dimension differ !!
stopifnot(mEQ(sm, mm)) ## but that's ok, since  mm has  all-0 column !
##
## look at this :
all(mm[,"d5"] == 0)  ## !!!! --- correct: a column of all 0  <--> dropped level!
stopifnot(all.equal(sm, mm[, - which("d5" == colnames(mm))])) ## indeed !
## i.e., sm has just dropped an all zero column --- which it should!
