library(Matrix)
h9 <- Hilbert(9)
str(h9)
f9 <- as(chol(h9), "trMatrix")
str(h9)# now has factorization
if(is.null(h9@factorization))
  stop("didn't auto-update factorization slot")
str(f9)
rcond(h9)
rcond(f9)
options(digits=4)
crossprod(f9)# looks the same as
h9 # but not internally!
##   i.e. this is all wrong : all.equal(h9, crossprod(f9))
