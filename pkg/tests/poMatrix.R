### Testing positive definite matrices

library(Matrix)

stopifnot(c(0,0) == dim(Hilbert(0)))

h9 <- Hilbert(9)
str(h9)
all.equal(determinant(h9)$modulus, -96.7369450737858, tol= 1e-15)
stopifnot(0 == length(h9@factorization))# nothing yet
str(f9 <- as(chol(h9), "trMatrix"))
## h9 now has factorization
stopifnot(names(h9@factorization) == "Cholesky")
rcond(h9)
rcond(f9)
str(h9)# has 'rcond' and 'factorization'
options(digits=4)
crossprod(f9)# looks the same as
h9 # but not internally!
##   i.e. this is all wrong : all.equal(h9, crossprod(f9))
