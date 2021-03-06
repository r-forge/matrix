#### Currently a collection of simple tests
##	(since 'Matrix' takes long to load, rather have fewer source files!)

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

### Matrix() ''smartness''
(d4 <- Matrix(diag(4)))
(z4 <- Matrix(0*diag(4)))
(o4 <- Matrix(1+diag(4)))
(m4 <- Matrix(cbind(0,rbind(6*diag(3),0))))
dm4 <- Matrix(m4, sparse = FALSE)
class(mN <-  Matrix(NA, 3,4)) # NA *is* logical
stopifnot(validObject(d4), validObject(z4), validObject(o4),
          validObject(m4), validObject(dm4), validObject(mN))
assert.EQ.mat(dm4, as(m4, "matrix"))
if(FALSE)# assert.EQ.. fails for all NA
assert.EQ.mat(mN, matrix(NA, 3,4))
sN <- Matrix(, 3,4, sparse=TRUE)# warning: NA coerced to FALSE
stopifnot(length(sN@i) == 0, # all "FALSE"
          validObject(Matrix(c(NA,0), 4, 3, byrow = TRUE)),
          validObject(Matrix(c(NA,0), 4, 4)),
          is(Matrix(c(NA,0,0,0), 4, 4), "sparseMatrix"))

## large sparse ones: these now directly "go sparse":
str(m0 <- Matrix(0,     nrow=100, ncol = 1000))
str(l0 <- Matrix(FALSE, nrow=100, ncol = 200))

###--  Sparse Triangular :

(t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2,
           Dim= as.integer(c(4,4))))
stopifnot(validObject(t1),
          validObject(t1c <- as(t1, "dtCMatrix")))
assert.EQ.mat(t1, as(t1c, "matrix"))

## from  0-diagonal to unit-diagonal {low-level step}:
tu <- t1 ; tu@diag <- "U"
tu
stopifnot(validObject(cu <- as(tu, "dtCMatrix")),
	  validObject(tu. <- as(cu, "dtTMatrix")),
          ## NOT: identical(tu, tu.), # since T* is not unique!
	  identical(cu, as(tu., "dtCMatrix")),
	  validObject(t(cu)),
	  validObject(t(tu)))
assert.EQ.mat(cu, as(tu,"matrix"), tol=0)


###-- Numeric Dense: Crossprod & Solve

set.seed(123)
mm. <- mm <- Matrix(rnorm(500 * 150), nc = 150)
stopifnot(validObject(mm))
xpx <- crossprod(mm)
stopifnot(identical(mm, mm.),# once upon a time, mm was altered by crossprod()
          validObject(xpx))
str(mm) # 'dge*"
str(xpx)# 'dpo*"
xpy <- crossprod(mm, rnorm(500))
res <- solve(xpx, xpy)
str(xpx)# now with Cholesky factor
stopifnot(validObject(xpx),
          validObject(xpy),
          validObject(res))
stopifnot(all.equal(xpx %*% res, xpy, tol= 1e-12))

###-- more solve() methods  {was ./solve.R }

## first for "dgeMatrix" and all kinds of RHS :
(m6 <- 1 + as(diag(0:5), "dgeMatrix"))
rcond(m6)
I6 <- as(diag(6), "dgeMatrix")
stopifnot(all.equal(I6, m6 %*% solve(m6)),
          all.equal(I6, solve(m6) %*% m6) )

(i6 <- solve(m6, Matrix(1:6)))
stopifnot(identical(i6, as(cbind(c(-4, rep(1,5))), "dgeMatrix")),
          identical(i6, solve(m6, 1:6)),
          identical(i6, solve(m6, matrix(1:6))),
          identical(i6, solve(m6, matrix(c(1,2,3,4,5,6))))
          )

###-- row- and column operations  {was ./rowcolOps.R }

set.seed(321)
(m1 <- round(Matrix(rnorm(25), 5), 2))
m1k <- Matrix(round(rnorm(1000), 2), 50, 20)
m.m <- as(m1k, "matrix")
stopifnot(all.equal(colMeans(m1k), colMeans(m.m)),
          all.equal(colSums (m1k), colSums (m.m)),
          all.equal(rowMeans(m1k), rowMeans(m.m)),
          all.equal(rowSums (m1k), rowSums (m.m)))

###-- kronecker for nonsparse uses Matrix(.):
stopifnot(is(kr <- kronecker(m1, m6), "Matrix"))
assert.EQ.mat(kr,
              kronecker(as(m1, "matrix"),
                        as(m6, "matrix")),
              tol = 0)
## sparse:
(kt1 <- kronecker(t1, tu))
kt2 <- kronecker(t1c, cu)
ktf <- kronecker(as.matrix(t1), as.matrix(tu))

assert.EQ.mat(kt1, ktf, tol= 0)
assert.EQ.mat(kt2, ktf, tol= 0)
## but kt1 and kt2, both "dgT" are different since entries are not ordered!

## coercion from "dpo" or "dsy"
xx <- as(xpx, "dsyMatrix")
stopifnot(isSymmetric(xxS  <- as(xx,  "sparseMatrix")),
          isSymmetric(xpxS <- as(xpx, "sparseMatrix")))

tm <- matrix(0, 8,8)
tm[cbind(c(1,1,2,7,8),
         c(3,6,4,8,8))] <- c(2,-30,15,20,80)
(tM <- Matrix(tm))                ## dtC
(mM <- Matrix(m <- (tm + t(tm)))) ## dsC
mT <- as(mM, "dsTMatrix")
gC <- as(as(mT, "dgTMatrix"), "dgCMatrix")
## Check that 'mT' and gC print properly :
pr.mT <- capture.output(mT)
nn <- unlist(strsplit(gsub(" +\\.", "", sub("^....", "", pr.mT[-(1:2)])), " "))
stopifnot(as.numeric(nn[nn != ""]) == m[m != 0],
          capture.output(gC)[-1] == pr.mT[-1])
assert.EQ.mat(tM, tm, tol=0)
assert.EQ.mat(gC, m,  tol=0)
assert.EQ.mat(mT, m,  tol=0)
stopifnot(is(mM, "dsCMatrix"), is(tM, "dtCMatrix"),
          ## coercions  general <-> symmetric
          identical(as(as(mM, "dgCMatrix"), "dsCMatrix"), mM),
          identical(as(as(mM, "dgTMatrix"), "dsTMatrix"), mT),
          identical(as(as(tM, "dgCMatrix"), "dtCMatrix"), tM)
)
eM <- eigen(mM) # works thanks to base::as.matrix hack in ../R/zzz.R
stopifnot(all.equal(eM$values,
                { v <- c(162.462112512353, 30.0665927567458)
                  c(v, 15, 0, 0, 160-v[1], -15, -v[2])}, tol=1e-14))

##--- symmetric -> pos.def. needs valid test:
m5 <- Matrix(diag(5) - 1)
if(FALSE) # FIXME: this happily "works" but MM thinks it shouldn't:
assertError(as(m5, "dpoMatrix"))


###-- logical sparse : ----------

(lkt <- as(as(kt1, "dgCMatrix"), "lgCMatrix"))# ok
(clt <- crossprod(lkt))
if(FALSE) ## FIXME!
    crossprod(clt)
## CHOLMOD error: matrix cannot be symmetric


### "d" <-> "l"  for (symmetric) sparse :
data(KNex)
mm <- KNex$mm
xpx <- crossprod(mm)
lxpx <- as(xpx, "lsCMatrix")
if(FALSE)
    show(lxpx) ## gives error about "lsC" -> "lgT" coercion ..
## The bug is actually from *subsetting* the large matrix:
if(FALSE) ## FIXME
    r <- lxpx[1:2,]

lmm <- as(mm, "lgCMatrix")
xlx <- crossprod(lmm)
## now A = lxpx and B = xlx should be close, but not quite the same
## since <x,y> = 0 is well possible when x!=0 and y!=0 .
## However,  A[i,j] != 0 ==> B[i,j] != 0:
A <- as(as(lxpx, "lgCMatrix"), "lgTMatrix")
B <- as(as(xlx,  "lgCMatrix"), "lgTMatrix")
ij <- function(a) a@i + ncol(a) * a@j
stopifnot(all(ij(A) %in% ij(B)))




cat('Time elapsed: ', proc.time(),'\n') # "stats"
