### triangular packed
library(Matrix)
## round(): "Math2" group method
round(tp6 <- as(chol(Hilbert(6)),"dtpMatrix"), 3)
1/tp6 # "Arith" group : gives 'dgeMatrix'
str(tp6)
stopifnot(validObject(tp6),
          all.equal(tp6 %*% diag(6), as(tp6, "matrix")),
          class((tt6 <- t(tp6))) == "dtpMatrix",
          tp6@uplo == "U" && tt6@uplo == "L")
if(FALSE) # FIXME : still fails
diag(6) %*% tp6
(tr6 <- as(tp6, "dtrMatrix")) ## FIXME: prints horribly
D. <- determinant(tp6)
rc <- rcond(tp6)
stopifnot(all.equal(D.$modulus, -6.579251212),
          all.equal(rc, 1.791511257e-4),
          rc == tp6@rcond)
norm(tp6, "I")
norm(tp6, "1")
norm(tp6, "F")
object.size(tp6)
object.size(as(tp6, "dtrMatrix"))
object.size(as(tp6, "matrix"))

## larger case
set.seed(123)
rl <- new("dtpMatrix", uplo="L", diag="N", Dim = rep.int(1000:1000,2),
          x = rnorm(500*1001))
validObject(rl)
str(rl)
norm(rl, "I")
norm(rl, "1")
norm(rl, "F")
rcond(rl)# 0 !
all.equal(rl %*% diag(1000), as(rl, "matrix"))
object.size(rl)
object.size(as(rl, "dtrMatrix"))
object.size(as(rl, "matrix"))
determinant(rl)
q('no')
