## Check consistency of methods for the sparseQR class

library(Matrix)

#source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

data(KNex); mm <- KNex$mm; y <- KNex$y
md <- as(mm, "matrix")                  # dense 

system.time(mmq <- qr(mm))
system.time(mdq <- qr(md))

## qr.qy and qr.qty should be inverses
all.equal(qr.qy(mmq, qr.qty(mmq, y))@x, y)
all.equal(qr.qty(mmq, qr.qy(mmq, y))@x, y)

## consistency of results dense and sparse
all.equal(qr.coef(mdq, y), qr.coef(mmq,y)@x)
all.equal(qr.resid(mdq, y), qr.resid(mmq,y)@x)
all.equal(qr.fitted(mdq, y), qr.fitted(mmq,y)@x)
