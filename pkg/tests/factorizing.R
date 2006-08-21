#### Matrix Factorizations  --- of all kinds

library(Matrix)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc


### "sparseQR" : Check consistency of methods
##   --------
data(KNex); mm <- KNex$mm; y <- KNex$y
stopifnot(is((Y <- Matrix(y)), "dgeMatrix"))
md <- as(mm, "matrix")                  # dense

system.time(mmq <- qr(mm))
system.time(mdq <- qr(md))# much slower

## qr.qy and qr.qty should be inverses
stopifnot(all.equal(qr.qy (mmq, qr.qty(mmq, y))@x, y),
          all.equal(qr.qty(mmq, qr.qy (mmq, y))@x, y),
          all.equal(qr.qty(mmq, y), qr.qty(mmq, Y)) )

## consistency of results dense and sparse
stopifnot(is.all.equal3(qr.coef  (mdq, y), qr.coef  (mmq,y)@x, qr.coef  (mmq,Y)@x) ,
          is.all.equal3(qr.resid (mdq, y), qr.resid (mmq,y)@x, qr.resid (mmq,Y)@x) ,
          is.all.equal3(qr.fitted(mdq, y), qr.fitted(mmq,y)@x, qr.fitted(mmq,Y)@x) )

