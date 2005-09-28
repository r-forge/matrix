#### Currently a collection of simple tests
##	(since 'Matrix' takes long to load, rather have fewer source files!)

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

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
          validObject(t(cu)),
          validObject(t(tu)))
assert.EQ.mat(cu, as(tu,"matrix"), tol=0)


###-- Numeric Dense: Crossprod & Solve

set.seed(123)
mm <- Matrix(rnorm(500 * 150), nc = 150)
stopifnot(validObject(mm))
xpx <- crossprod(mm)# alters mm !
stopifnot(validObject(mm),
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
mm <- Matrix(round(rnorm(1000), 2), 50, 20)
m1 <- as(mm, "matrix")
stopifnot(all.equal(colMeans(mm), colMeans(m1)),
          all.equal(colSums(mm), colSums(m1)),
          all.equal(rowMeans(mm), rowMeans(m1)),
          all.equal(rowSums(mm), rowSums(m1)))

###-- Testing expansions of factorizations {was ./expand.R }

(m1 <- round(Matrix(rnorm(25), 5), 2))
(lul <- expand(lu(m1)))
stopifnot(all.equal(as(m1, "matrix"),
                    as(lul$P %*% (lul$L %*% lul$U), "matrix")))


proc.time()
