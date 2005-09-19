library(Matrix)

### Use ``non-unique'' versions of dgTMatrix objects

N <- 200
set.seed(1)
i <- as.integer(round(runif (N, 0, 100)))
j <- as.integer(3* rpois (N, lam=15))
x <- round(rnorm(N), 2)
which(duplicated(cbind(i,j))) # 8 index pairs are duplicated

m1 <- new("dgTMatrix", Dim = c(max(i)+1:1, max(j)+1:1), i = i, j = j, x = x)
mc <- as(m1, "dgCMatrix")
m2 <- as(mc, "dgTMatrix")## the same as 'm1' but without duplicates

stopifnot(!isTRUE(all.equal(m1, m2)),
          all.equal(as(m1,"matrix"), as(m2,"matrix"), tol=1e-15),
          all.equal(crossprod(m1), crossprod(m2), tol=1e-15),
          identical(mc, as(m2, "dgCMatrix")))

### -> uniq* functions now in ../R/Auxiliaries.R
(t2 <- system.time(um2 <- Matrix:::uniq(m1)))
