library(Matrix)

### Matrix Products including  cross products

source(system.file("test-tools.R", package = "Matrix"))

m5 <- 1 + as(diag(-1:4)[-5,], "dgeMatrix")
## named dimnames:
dimnames(m5) <- list(Rows= LETTERS[1:5], paste("C", 1:6, sep=""))
stopifnot(dim(m5) == 5:6,
          class(cm5 <- crossprod(m5)) == "dpoMatrix")
assert.EQ.mat((c.m5 <- t(m5) %*% m5), as(cm5, "matrix"))
## classes differ; but the 'dimnames' are *both* missing -- FIXME
tc.m5 <- m5 %*% t(m5)
(tcm5 <- tcrossprod(m5)) # "dpo*"
assert.EQ.mat(tc.m5, mm5 <- as(tcm5, "matrix"))# missing dimnames - FIXME
## tcrossprod(x,y) :
assert.EQ.mat(tcrossprod(m5, m5), mm5)
assert.EQ.mat(tcrossprod(m5, as(m5,"matrix")), mm5)
assert.EQ.mat(tcrossprod(as(m5,"matrix"), m5), mm5)

## simple cases with 'scalars' treated as 1x1 matrices:
d <- Matrix(1:5)
d %*% 2
10 %*% t(d)
assertError(3 %*% d)             # must give an error , similar to
assertError(5 %*% as.matrix(d))  # -> error

## right and left "numeric" and "matrix" multiplication:
(p1 <- m5 %*% c(10, 2:6))
(p2 <- c(10, 2:5) %*% m5)
(pd1 <- m5 %*% diag(1:6))
(pd2 <- diag(10:6) %*% m5)
stopifnot(dim(crossprod(t(m5))) == c(5,5),
          c(class(p1),class(p2),class(pd1),class(pd2)) == "dgeMatrix"
          )
assert.EQ.mat(p1, cbind(c(20,30,33,38,54)))
assert.EQ.mat(pd1, as(m5,"matrix") %*% diag(1:6))
assert.EQ.mat(pd2, diag(10:6) %*% as(m5,"matrix"))


M <- mm[1:500, 1:200]
cpr <- t(mm) %*% mm
showMethods("%*%", class=class(M))

v1 <- rep(1, ncol(M))
str(r <-  M %*% Matrix(v1))
str(r. <- M %*% cbind(v1))
stopifnot(identical3(r, r., M %*% as(v1, "matrix")))

v2 <- rep(1,nrow(M))
r2 <- t(Matrix(v2)) %*% M
str(r2. <- v2 %*% M)
stopifnot(identical4(r2, r2., rbind(v2) %*% M, t(as(v2, "matrix")) %*% M))


###--- "logical" Matrices : ---------------------

## Robert's Example, a bit more readable
fromTo <- rbind(c(2,10),
                c(3, 9))
N <- 10
nrFT <- nrow(fromTo)
rowi <- rep.int(1:nrFT, fromTo[,2]-fromTo[,1] + 1) - 1:1
coli <- unlist(lapply(1:nrFT, function(x) fromTo[x,1]:fromTo[x,2])) - 1:1
sM  <- new("lgTMatrix", i = rowi, j=coli, Dim=as.integer(c(N,N)))
sM # nice

sm <- as(sM, "matrix")
sM %*% sM
assert.EQ.mat(sM %*% sM,        sm %*% sm)
assert.EQ.mat(t(sM) %*% sM,
              (t(sm) %*% sm) > 0, tol=0)
crossprod(sM)
tcrossprod(sM)
stopifnot(identical( crossprod(sM), t(sM) %*%   sM),
          identical(tcrossprod(sM),   sM  %*% t(sM)))

assert.EQ.mat( crossprod(sM),  crossprod(sm) > 0)
assert.EQ.mat(tcrossprod(sM), as(tcrossprod(sm),"matrix") > 0)

proc.time() # for ``statistical reasons''
