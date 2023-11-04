## For documenting performance +/-

stopifnot(requireNamespace("microbenchmark"))
BM <- microbenchmark::microbenchmark

library(Matrix)
set.seed(82919)
options(width = 200L)
MatrixVersion <- packageVersion("Matrix")

dgC <- rsparsematrix(1000L, 1000L, 0.01)
dgR <- as(dgC, "RsparseMatrix")
dgT <- as(dgC, "TsparseMatrix")
dsC <- forceSymmetric(dgC)
dsR <- forceSymmetric(dgR)
dsT <- forceSymmetric(dgT)
dtC <- triu(dgC)
dtR <- triu(dgR)
dtT <- triu(dgT)
ddC <- band(dgC, 0L, 0L)
ddR <- band(dgR, 0L, 0L)
ddT <- band(dgT, 0L, 0L)

dge <- new("dgeMatrix", Dim = c(1000L, 1000L), x = rnorm(1e+06L))
dsy <- forceSymmetric(dge)
dtr <- triu(dge)
dsp <- pack(dsy)
dtp <- pack(dtr)

if(MatrixVersion <= "1.4.1") {
    ## due to "bugs" in Matrix <= 1.4-1
    dsR <- as(dsR, "RsparseMatrix")
    dsT <- as(dsT, "TsparseMatrix")
    ddR <- as(ddR, "RsparseMatrix")
}


## TODO: many more improvements since 1.4-1 and yet more since 1.4-0 ...


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Getting diagonal from [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(diag(dgC), diag(dgR), diag(dgT),
   diag(dsC), diag(dsR), diag(dsT),
   diag(dtC), diag(dtR), diag(dtT),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Transpose of [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(t(dgC), t(dgR), t(dgT),
   t(dsC), t(dsR), t(dsT),
   t(dtC), t(dtR), t(dtT),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Coerce .g[CRT]Matrix to .[st][CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(forceSymmetric(dgC), triu(dgC), tril(dgC),
   forceSymmetric(dgR), triu(dgR), tril(dgR),
   forceSymmetric(dgT), triu(dgT), tril(dgT),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Test for symmetric .g[CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dsC.as.g <- as(dsC, "generalMatrix")
dsR.as.g <- as(dsR, "generalMatrix")
dsT.as.g <- as(dsT, "generalMatrix")

BM(isSymmetric(dsC.as.g),
   isSymmetric(dsR.as.g),
   isSymmetric(dsT.as.g),
   isSymmetric(dsC.as.g, tol = 0),
   isSymmetric(dsR.as.g, tol = 0),
   isSymmetric(dsT.as.g, tol = 0),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Test for triangular .g[CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dtC.as.g <- as(dtC, "generalMatrix")
dtR.as.g <- as(dtR, "generalMatrix")
dtT.as.g <- as(dtT, "generalMatrix")

BM(isTriangular(dtC.as.g, upper = TRUE),
   isTriangular(dtR.as.g, upper = TRUE),
   isTriangular(dtT.as.g, upper = TRUE),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Test for diagonal .g[CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ddC.as.g <- as(ddC, "generalMatrix")
ddR.as.g <- as(ddR, "generalMatrix")
ddT.as.g <- as(ddT, "generalMatrix")

BM(isDiagonal(ddC.as.g),
   isDiagonal(ddR.as.g),
   isDiagonal(ddT.as.g),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Determinant of positive definite ds[yp]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dsy.spd <- as(tcrossprod(dsy), "dsyMatrix")
dsp.spd <- pack(dsy.spd)

BM(determinant(dsy, TRUE), determinant(dsy.spd, TRUE),
   determinant(dsp, TRUE), determinant(dsp.spd, TRUE),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Getting (skew-)symmetric part of denseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(symmpart(dge), symmpart(dsy), symmpart(dtr), symmpart(dsp), symmpart(dtp),
   skewpart(dge), skewpart(dsy), skewpart(dtr), skewpart(dsp), skewpart(dtp),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Getting (skew-)symmetric part of [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(symmpart(dgC), symmpart(dsC), symmpart(dtC),
   symmpart(dgR), symmpart(dsR), symmpart(dtR),
   symmpart(dgT), symmpart(dsT), symmpart(dtT),
   skewpart(dgC), skewpart(dsC), skewpart(dtC),
   skewpart(dgR), skewpart(dsR), skewpart(dtR),
   skewpart(dgT), skewpart(dsT), skewpart(dtT),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Coercing between [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(as(dgC, "RsparseMatrix"), as(dgC, "TsparseMatrix"),
   as(dgR, "CsparseMatrix"), as(dgR, "TsparseMatrix"),
   as(dgT, "CsparseMatrix"), as(dgT, "RsparseMatrix"),
   as(dsC, "RsparseMatrix"), as(dsC, "TsparseMatrix"),
   as(dsR, "CsparseMatrix"), as(dsR, "TsparseMatrix"),
   as(dsT, "CsparseMatrix"), as(dsT, "RsparseMatrix"),
   as(dtC, "RsparseMatrix"), as(dtC, "TsparseMatrix"),
   as(dtR, "CsparseMatrix"), as(dtR, "TsparseMatrix"),
   as(dtT, "CsparseMatrix"), as(dtT, "RsparseMatrix"),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Coercing from denseMatrix to [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(as(dge, "CsparseMatrix"), as(dge, "RsparseMatrix"), as(dge, "TsparseMatrix"),
   as(dsy, "CsparseMatrix"), as(dsy, "RsparseMatrix"), as(dsy, "TsparseMatrix"),
   as(dtr, "CsparseMatrix"), as(dtr, "RsparseMatrix"), as(dtr, "TsparseMatrix"),
   as(dsp, "CsparseMatrix"), as(dsp, "RsparseMatrix"), as(dsp, "TsparseMatrix"),
   as(dtp, "CsparseMatrix"), as(dtp, "RsparseMatrix"), as(dtp, "TsparseMatrix"),
   unit = "microseconds")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Coercing from [CRT]sparseMatrix to denseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(MatrixVersion > "1.4.1") {
    .to.unpacked <- function(from) as(from, "unpackedMatrix")
    .to.packed   <- function(from) as(from, "packedMatrix")
} else {
    .to.unpacked <- function(from) as(from, "denseMatrix")
    .to.packed   <- function(from) pack(as(from, "denseMatrix"))
}

BM(.to.unpacked(dgC), .to.unpacked(dsC), .to.unpacked(dtC),
   .to.packed(dsC), .to.packed(dtC),
   .to.unpacked(dgR), .to.unpacked(dsR), .to.unpacked(dtR),
   .to.packed(dsR), .to.packed(dtR),
   .to.unpacked(dgT), .to.unpacked(dsT), .to.unpacked(dtT),
   .to.packed(dsT), .to.packed(dtT),
   unit = "microseconds")
