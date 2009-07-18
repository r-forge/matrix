#### Matrix Factorizations  --- of all kinds

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc


### "sparseQR" : Check consistency of methods
##   --------
data(KNex); mm <- KNex$mm; y <- KNex$y
stopifnot(is((Y <- Matrix(y)), "dgeMatrix"))
md <- as(mm, "matrix")                  # dense

system.time(mmq <- qr(mm))
system.time(mdq <- qr(md))# much (~ 150 x) slower

## qr.qy and qr.qty should be inverses
stopifnot(all.equal(qr.qy (mmq, qr.qty(mmq, y))@x, y),
          all.equal(qr.qty(mmq, qr.qy (mmq, y))@x, y),
          all.equal(qr.qty(mmq, y), qr.qty(mmq, Y)) )

## consistency of results dense and sparse
stopifnot(is.all.equal3(qr.coef  (mdq, y), qr.coef  (mmq,y)@x, qr.coef  (mmq,Y)@x) ,
          is.all.equal3(qr.resid (mdq, y), qr.resid (mmq,y)@x, qr.resid (mmq,Y)@x) ,
          is.all.equal3(qr.fitted(mdq, y), qr.fitted(mmq,y)@x, qr.fitted(mmq,Y)@x) )


### "denseLU"

## Testing expansions of factorizations {was ./expand.R, then in simple.R }

set.seed(1)
(m1 <- round(Matrix(rnorm(25), 5), 2))
str(lu1 <- lu(m1))
(luX <- expand(lu1))
stopifnot(all.equal(as(m1, "matrix"),
                    as(luX$P %*% (luX$L %*% luX$U), "matrix")))

### "sparseLU"
por1 <- readMM(system.file("external/pores_1.mtx", package = "Matrix"))
lu1 <- lu(por1)
pm <- as(por1, "CsparseMatrix")
(pmLU <- lu(pm)) # -> show(<MatrixFactorization>)
xp <- expand(pmLU)
## permute rows and columns of original matrix
ppm <- pm[pmLU@p + 1:1, pmLU@q + 1:1]
Ppm <- pmLU@L %*% pmLU@U
## identical only as long as we don't keep the original class info:
stopifnot(identical(lu1, pmLU),
	  identical(ppm, with(xp, P %*% pm %*% t(Q))),
	  sapply(xp, is, class="Matrix"))


## these two should be the same, and `are' in some ways:
assert.EQ.mat(ppm, as(Ppm, "matrix"), tol = 1e-14)
## *however*
length(ppm@x)# 180
length(Ppm@x)# 317 !
table(Ppm@x == 0)# (194, 123) - has 123 "zero" and 14 ``almost zero" entries


###________ Cholesky() ________

##--------  LDL' ---- small exact examples

set.seed(1)
for(n in c(5:12)) {
    cat("\nn = ",n,"\n-------\n")
    rr <- mkLDL(n)
    ##    -------- from 'test-tools.R'
    stopifnot(all(with(rr, A == as(L %*% D %*% t(L),
                           "symmetricMatrix"))))
    d <- rr$d.half
    A <- rr$A
    R <- chol(A)
    print(d. <- diag(R))
    D. <- Diagonal(x= d.^2)
    L. <- t(R) %*% Diagonal(x = 1/d.)
    stopifnot(all.equal(as.matrix(D.), as.matrix(rr$ D)),
              all.equal(as.matrix(L.), as.matrix(rr$ L)))
    ##
    CAp <- Cholesky(A)# perm=TRUE --> Permutation:
    p <- CAp@perm + 1L
    P <- as(p, "pMatrix")
    ## the inverse permutation:
    invP <- solve(P)@perm
    lDet <- sum(2* log(d))# the "true" value
    ldet <- Matrix:::.diag.dsC(Chx = CAp, res.kind = "sumLog")
    ##
    CA	<- Cholesky(A,perm=FALSE)
    ldet2 <- Matrix:::.diag.dsC(Chx = CA, res.kind = "sumLog")
    ## not printing CAp : ends up non-integer for n >= 11
    mCAp <- as(CAp,"sparseMatrix")
    print(mCA  <- drop0(as(CA, "sparseMatrix")))
    stopifnot(identical(A[p,p], as(P %*% A %*% t(P),
				   "symmetricMatrix")),
	      all.equal(lDet, sum(log(Matrix:::.diag.dsC(Chx= CAp,res.kind="diag")))),
	      relErr(d.^2, Matrix:::.diag.dsC(Chx= CA, res.kind="diag")) < 1e-14,
	      all.equal(lDet, ldet),
	      all.equal(lDet, ldet2),
	      relErr(A[p,p], tcrossprod(mCAp)) < 1e-14)
}## for()

set.seed(17)
(rr <- mkLDL(4))
(CA <- Cholesky(rr$A))
stopifnot(all.equal(determinant(rr$A),
		    determinant(as(rr$A, "matrix"))))
A12 <- mkLDL(12, 1/10)

(r12 <- allCholesky(A12$A))
aCh.hash <- r12$r.all %*% (2^(2:0))
if(FALSE)## if(require("sfsmisc"))
split(rownames(r12$r.all), Duplicated(aCh.hash))

## TODO: find cases for both choices when we leave it to CHOLMOD to chose
if(FALSE) ## not yet
for(n in 1:9) { ## # before seg.fault at n = 10 !
    cat(sprintf("n = %3d:\n", n))
    mkA <- mkLDL(1+rpois(1, 30), 1/10)
    r <- allCholesky(mkA$A)
    ## Compare .. apart from the NAs that happen from (perm=FALSE, super=TRUE)
    iNA <- apply(is.na(r$r.all), 1, any)
    stopifnot(aCh.hash[!iNA] == r$r.all[!iNA,] %*% (2^(2:0)))
    cat("--------\n")
}
### --- see seg.fault on  nb-mm {for the R version I can't 'dbg' ..}
## FIXME -- UNfinished

A. <-
    new("dsCMatrix", Dim = c(32L, 32L), uplo = "U"
        , i = c(0L, 1L, 2L, 3L, 4L, 2L, 5L, 5L, 6L, 0L, 7L, 8L, 8L, 9L, 3L,
          4L, 10L, 11L, 0L, 7L, 12L, 6L, 13L, 14L, 4L, 10L, 15L, 8L, 9L,
          16L, 17L, 1L, 2L, 5L, 18L, 6L, 9L, 13L, 15L, 19L, 12L, 20L, 0L,
          7L, 8L, 9L, 12L, 16L, 21L, 7L, 21L, 22L, 9L, 19L, 23L, 2L, 5L,
          18L, 19L, 21L, 22L, 24L, 10L, 11L, 18L, 25L, 0L, 6L, 7L, 12L,
          13L, 18L, 19L, 20L, 21L, 22L, 24L, 25L, 26L, 6L, 9L, 13L, 15L,
          16L, 19L, 23L, 26L, 27L, 9L, 11L, 19L, 23L, 25L, 27L, 28L, 1L,
          4L, 7L, 10L, 15L, 18L, 22L, 24L, 25L, 28L, 29L, 1L, 18L, 19L,
          24L, 29L, 30L, 14L, 21L, 22L, 24L, 26L, 31L)
        ##
        , p = c(0L, 1L, 2L, 3L, 4L, 5L, 7L, 9L, 11L, 12L, 14L, 17L, 18L, 21L,
          23L, 24L, 27L, 30L, 31L, 35L, 40L, 42L, 49L, 52L, 55L, 62L, 66L,
          79L, 88L, 95L, 106L, 112L, 118L)
        ##
        , x = c(25, 225, 441, 484, 529, 36603, 3038085, 828, 19333, 500, 10256,
          841, 52983, 3338253, 26136, 52371, 6596473, 121, 75, 1500, 1186,
          24276, 2039809, 144, 2116, 209484, 8825, 1682, 105966, 3428, 9, 21150,
          33516, 2781828, 4535397, 6358, 32724, 534072, 6859, 3576050, 49011,
          2500237, 2325, 46500, 47937, 3020031, 6975, 95874, 2949210, 1536,
          4032, 38340, 31752, 3206952, 3112720, 4851, 402633, 368676, 37908,
          17856, 124992, 2578213, 13200, 9317, 4455, 1398038, 200, 27455, 4000,
          79402, 2306220, 5832, 604010, 4018902, 56040, 262080, 1160640,
          320760, 11925262, 16184, 25596, 1401956, 23826, 3008, 3393938,
          2508408, 1537480, 7532296, 28512, 5445, 2879712, 2794176, 419505,
          2252448, 2768530, 225, 14812, 11264, 1466388, 59248, 21150, 67584,
          2500, 116, 6960, 976637, 3375, 317250, 74358, 3866616, 20231,
          9084758, 2592, 39744, 278208, 1232064, 2583360, 2789776)
        )
validObject(A.)

for(p in c(FALSE,TRUE)) # no NA here
    for(L in c(FALSE,TRUE, NA))
        for(s in c(FALSE,TRUE, NA)) {
            cat(sprintf("p,L,S = (%2d,%2d,%2d): ", p,L,s))
            r <- tryCatch(Cholesky(A., perm=p, LDL=L, super=s),
                          error = function(e)e)
            cat(if(inherits(r, "error")) " *** E ***" else
                sprintf("%3d", r@type),"\n", sep="")
        }
str(A., max=3) ## look at the 'factors'

facs <- A.@factors
names(facs) <- sub("Cholesky$", "", names(facs))
facs <- facs[order(names(facs))]

sapply(facs, class)
str(lapply(facs, slot, "type"))
## super = TRUE  currently always entails  LDL=FALSE :
## hence isLDL is TRUE for ("D" and not "S"):
sapply(facs, isLDL)

## --- now a "large" (712 x 712) real data example

data(KNex)
mtm <- with(KNex, crossprod(mm))
ld.3 <- .Call("dsCMatrix_LDL_D", mtm, perm=TRUE,  "sumLog")
stopifnot(names(mtm@factors) == "sPDCholesky")
ld.4 <- .Call("dsCMatrix_LDL_D", mtm, perm=FALSE, "sumLog")# clearly slower
stopifnot(names(mtm@factors) == paste(c("sPD", "spD"),"Cholesky", sep=''))
c2 <- Cholesky(mtm, super = TRUE)
stopifnot(names(mtm@factors) == paste(c("sPD", "spD", "SPd"),
               "Cholesky", sep=''))

## is now taken from cache
c1 <- Cholesky(mtm)

bv <- 1:nrow(mtm) # even integer
b <- matrix(bv)
## solve(c2, b) by default solves  Ax = b, where A = c2'c2 !
x <- solve(c2,b)
stopifnot(identical3(x, solve(c2, bv), solve(c2, b, system = "A")),
          all.equal(x, solve(mtm, b)))
for(sys in c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt")) {
    x <- solve(c2, b,  system = sys)
    cat(sys,":\n"); print(head(x))
    stopifnot(dim(x) == c(712, 1),
              identical(x, solve(c2, bv, system = sys)))
}

## log(|LL'|) - check if super = TRUE and simplicial give same determinant
ld1 <- .Call("CHMfactor_ldetL2", c1)
ld2 <- .Call("CHMfactor_ldetL2", c2)
(ld1. <- determinant(mtm))
## experimental
ld3 <- .Call("dsCMatrix_LDL_D", mtm, TRUE, "sumLog")
ld4 <- .Call("dsCMatrix_LDL_D", mtm, FALSE, "sumLog")
stopifnot(all.equal(ld1, ld2),
	  is.all.equal3(ld2, ld3, ld4),
	  all.equal(ld.3, ld3, tol = 1e-14),
	  all.equal(ld.4, ld4, tol = 1e-14),
	  all.equal(ld1, as.vector(ld1.$modulus), tol = 1e-14))

## Some timing measurements
mtm <- with(KNex, crossprod(mm))
I <- .symDiagonal(n=nrow(mtm))
set.seed(101); r <- runif(100)

system.time(D1 <- sapply(r, function(rho) Matrix:::ldet1.dsC(mtm + (1/rho) * I)))
## 0.842 on fast cmath-5
system.time(D2 <- sapply(r, function(rho) Matrix:::ldet2.dsC(mtm + (1/rho) * I)))
## 0.819
system.time(D3 <- sapply(r, function(rho) Matrix:::ldet3.dsC(mtm + (1/rho) * I)))
## 0.810
stopifnot(is.all.equal3(D1,D2,D3, tol = 1e-13))

## Updating LL'  should remain LL' and not become  LDL' :
if(FALSE) {
    data(Dyestuff, package = "lme4")
    Zt <- as(Dyestuff$Batch, "sparseMatrix")
} else {
    Zt <- new("dgCMatrix", Dim = c(6L, 30L), x = rep(1, 30),
              i = rep(0:5, each=5),
              p = 0:30, Dimnames = list(LETTERS[1:6], NULL))
}
Ut <- 0.78 * Zt
L <- Cholesky(tcrossprod(Ut), LDL = FALSE, Imult = 1)
L1 <- update(L, tcrossprod(Ut), mult = 1)
stopifnot(all.equal(L, L1))


## Schur() ----------------------
checkSchur <- function(A, SchurA = Schur(A), tol = 1e-14) {
    stopifnot(is(SchurA, "Schur"),
              isOrthogonal(Q <- SchurA@Q),
              all.equal(as.mat(A),
                        as.mat(Q %*% SchurA@T %*% t(Q)), tol = tol))
}

SH <- Schur(H5 <- Hilbert(5))
checkSchur(H5, SH)
checkSchur(Diagonal(x = 9:3))

p <- 4L
uTp <- new("dtpMatrix", x=c(2, 3, -1, 4:6, -2:1), Dim = c(p,p))
(uT <- as(uTp, "dtrMatrix"))
## Schur ( <general> )  <--> Schur( <triangular> )
Su <- Schur(uT) ;   checkSchur(uT, Su)
gT <- as(uT,"generalMatrix")
Sg <- Schur(gT) ;   checkSchur(gT, Sg)
Stg <- Schur(t(gT));checkSchur(t(gT), Stg)
Stu <- Schur(t(uT));checkSchur(t(uT), Stu)

stopifnot(identical3(Sg@T, uT, Su@T),
          identical(Sg@Q, as(diag(p), "dgeMatrix")),
          identical(Stg@T, as(t(gT[,p:1])[,p:1], "triangularMatrix")),
          identical(Stg@Q, as(diag(p)[,p:1], "dgeMatrix")),
          identical(Stu@T, Stg@T))
assert.EQ.mat(Stu@Q, as(Stg@Q,"matrix"), tol=0)

## the pedigreemm example where solve(.) failed:
p <- new("dtCMatrix", i = c(2L, 3L, 2L, 5L, 4L, 4:5), p = c(0L, 2L, 4:7, 7L),
	 Dim = c(6L, 6L), Dimnames = list(as.character(1:6), NULL),
	 x = rep.int(-0.5, 7), uplo = "L", diag = "U")
ip <- solve(p)
assert.EQ.mat(solve(ip), as(p,"matrix"))
