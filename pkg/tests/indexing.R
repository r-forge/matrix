#### For both 'Extract' ("[") and 'Replace' ("[<-") Method testing

library(Matrix)

identical3 <- function(x,y,z)	identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d) identical(a,b) && identical3(b,c,d)

### Dense Matrices

m <- Matrix(1:28, nrow = 7)
validObject(m) ; m@x <- as.double(m@x) ; validObject(m)
stopifnot(identical(m, m[]),
          identical(m[2, 3],  16), # simple number
          identical(m[2, 3:4], c(16,23))) # simple numeric of length 2

m[2, 3:4, drop=FALSE] # sub matrix of class 'dgeMatrix'
m[-(4:7), 3:4]        # dito; the upper right corner of 'm'

## rows or columns only:
m[1,]     # first row, as simple numeric vector
m[,2]     # 2nd column
m[,1:2]   # sub matrix of first two columns
m[-(1:6),, drop=FALSE] # not the first 6 rows, i.e. only the 7th

## logical indexing
stopifnot(identical(m[2,3], m[(1:nrow(m)) == 2, (1:ncol(m)) == 3]),
          identical(m[2,], m[(1:nrow(m)) == 2, ]),
          identical(m[,3:4], m[, (1:4) >= 3]))

## dimnames index (TODO)

## TODO: more --- particularly once we have "m > 10" working!


### Sparse Matrices

m <- 1:800
set.seed(101) ; m[sample(800, 600)] <- 0
m <- Matrix(m, nrow = 40)
mm <- as(m, "matrix")
dimnames(mm) <- NULL ## << workaround: as(<sparse>, "matrix") has NULL dimnames
str(mC <- as(m, "dgCMatrix"))
str(mT <- as(m, "dgTMatrix"))
stopifnot(identical(mT, as(mC, "dgTMatrix")),
	  identical(mC, as(mT, "dgCMatrix")))

mC[,1]
mC[1:2,]
mC[7, drop = FALSE]

mT[,c(2,4)]
mT[1,]
mT[4, drop = FALSE]
stopifnot(identical3(mm[,1], mC[,1], mT[,1]),
	  identical3(mm[3,], mC[3,], mT[3,]),
	  identical3(mT[2,3], mC[2,3], 0),
	  identical(mT[], mT),
	  ## TODO: identical4() with  m[c(3,7), 2:4]
	  identical3(as(mC[c(3,7), 2:4],"matrix"), mm[c(3,7), 2:4],
		     as(mT[c(3,7), 2:4],"matrix")))


