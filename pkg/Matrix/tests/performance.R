## For documenting performance +/-

stopifnot(requireNamespace("microbenchmark"))
BM <- microbenchmark::microbenchmark

library(Matrix)
set.seed(82919)
options(width = 77L)

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

if(packageVersion("Matrix") <= "1.4.1") {
    ## due to "bugs" in Matrix <= 1.4-1
    dsR <- as(dsR, "RsparseMatrix")
    dsT <- as(dsT, "TsparseMatrix")
    ddR <- as(ddR, "RsparseMatrix")
}

## TODO: many more improvements since 1.4-1 and yet more since 1.4-0 ...


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## GETTING DIAGONAL FROM [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(diag(dgC), diag(dgR), diag(dgT),
   diag(dsC), diag(dsR), diag(dsT),
   diag(dtC), diag(dtR), diag(dtT),
   unit = "microseconds")

## Matrix r3504

## Unit: microseconds
##       expr   min      lq     mean  median      uq     max neval
##  diag(dgC) 8.897  9.1840 11.23646 10.0040 10.4140 105.985   100
##  diag(dgR) 9.840 10.1270 11.25614 11.0905 11.3365  16.769   100
##  diag(dgT) 7.175  7.5235  8.16228  8.0360  8.3230  14.719   100
##  diag(dsC) 3.362  3.6490  6.40912  3.8950  4.5920 227.099   100
##  diag(dsR) 3.198  3.5260  5.17133  3.7720  4.5510  96.309   100
##  diag(dsT) 5.043  5.2480  6.93966  5.6990  6.2730 120.663   100
##  diag(dtC) 3.321  3.6080  5.35829  3.8950  4.5510 120.622   100
##  diag(dtR) 3.116  3.5465  5.02250  3.8950  4.5100  96.022   100
##  diag(dtT) 5.043  5.3095  6.92654  5.8015  6.3140  96.391   100

## Matrix 1.4-1

## Unit: microseconds
##       expr     min       lq      mean   median       uq      max neval
##  diag(dgC)  61.418  68.1830  74.01279  70.6635  74.7225  275.807   100
##  diag(dgR) 141.081 156.8660 164.16851 162.3805 167.6285  278.882   100
##  diag(dgT) 152.602 162.0320 172.02575 168.8380 176.2180  323.900   100
##  diag(dsC)  52.275  54.5505  62.54878  57.4000  60.2290  478.798   100
##  diag(dsR) 108.773 118.9205 154.25389 123.2050 128.6785 2502.189   100
##  diag(dsT) 115.374 122.3030 133.98636 126.1365 133.0450  669.038   100
##  diag(dtC)  49.364  50.9220  57.29668  53.6485  57.2770  277.160   100
##  diag(dtR) 105.903 115.4970 125.06599 119.3920 125.5010  554.361   100
##  diag(dtT) 110.700 115.8455 126.35544 120.3350 124.2710  593.393   100


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TRANSPOSE OF [CRT]sparseMatrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(t(dgC), t(dgR), t(dgT),
   t(dsC), t(dsR), t(dsT),
   t(dtC), t(dtR), t(dtT),
   unit = "microseconds")

## Matrix r3504

## Unit: microseconds
##    expr    min      lq     mean  median      uq     max neval
##  t(dgC) 26.609 30.1555 31.37115 31.5290 32.2875  39.975   100
##  t(dgR) 27.429 30.1145 31.46750 31.4060 32.9230  38.622   100
##  t(dgT)  3.649  4.1820  4.61865  4.3870  4.5920  28.823   100
##  t(dsC) 13.858 15.5595 17.68945 16.1950 17.7325 126.895   100
##  t(dsR) 14.719 15.9285 19.17283 16.4000 17.8555 246.697   100
##  t(dsT)  3.690  4.2640  5.71294  4.5305  4.7150 128.289   100
##  t(dtC) 14.063 15.6210 17.94734 16.3180 18.2450 114.718   100
##  t(dtR) 14.555 15.9695 17.92848 16.4615 18.2245 120.991   100
##  t(dtT)  3.649  4.1000  5.38822  4.3870  4.6740 106.518   100

## Matrix 1.4-1

## Unit: microseconds
##    expr     min       lq      mean   median       uq      max neval
##  t(dgC)  48.175  53.4845  57.34096  55.9445  61.2130   84.091   100
##  t(dgR) 215.414 232.0600 270.26913 239.5220 248.5420 2583.000   100
##  t(dgT)  66.051  67.4860  70.52123  69.2080  72.3035   99.384   100
##  t(dsC)  28.249  32.2670  34.84795  33.9275  36.8795   54.079   100
##  t(dsR) 207.009 214.5940 233.69385 220.6825 232.4495 1066.656   100
##  t(dsT) 191.962 199.0960 222.01582 203.4830 216.2340 1395.353   100
##  t(dtC)  32.718  37.5970  41.63755  39.6880  42.4350  179.949   100
##  t(dtR) 200.326 207.7265 222.73701 211.9085 222.3840  872.972   100
##  t(dtT) 192.495 200.2850 243.67981 204.3235 214.7170 2352.047   100


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FORCE SYMMETRIC|TRIANGULAR for .g[CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BM(forceSymmetric(dgC), triu(dgC), tril(dgC),
   forceSymmetric(dgR), triu(dgR), tril(dgR),
   forceSymmetric(dgT), triu(dgT), tril(dgT),
   unit = "microseconds")

## Matrix r3504

## Unit: microseconds
##                 expr    min      lq     mean  median      uq     max neval
##  forceSymmetric(dgC) 26.445 29.2740 32.66552 31.7955 34.6040  51.291   100
##            triu(dgC) 30.668 33.8660 36.99143 36.1415 39.2575  51.824   100
##            tril(dgC) 29.028 33.6610 37.81225 36.7770 39.2780 143.254   100
##  forceSymmetric(dgR) 20.828 25.2765 27.47451 26.9575 28.7820  45.961   100
##            triu(dgR) 38.376 41.5740 44.07049 43.2755 46.2275  53.136   100
##            tril(dgR) 36.818 41.3485 45.84128 43.2345 45.5100 259.448   100
##  forceSymmetric(dgT) 24.682 26.8345 29.05998 27.8390 29.8070  44.198   100
##            triu(dgT) 27.306 29.1100 31.58353 30.7910 33.1690  48.339   100
##            tril(dgT) 27.429 28.8025 32.31743 30.0940 32.2875 138.334   100

## Matrix 1.4-1

## Unit: microseconds
##                 expr     min       lq      mean   median       uq      max
##  forceSymmetric(dgC)  34.850  40.8155  45.01185  44.4645  48.9745   62.156
##            triu(dgC)  47.888  54.6735  58.23599  57.5025  61.1105   83.599
##            tril(dgC)  67.568  72.0780  75.33996  74.6610  77.7975   94.833
##  forceSymmetric(dgR) 126.772 140.6300 147.39459 146.2265 153.0325  180.974
##            triu(dgR) 257.685 282.6130 384.61362 289.5625 299.9560 8089.423
##            tril(dgR) 280.399 298.8900 391.82552 305.6960 313.1375 8231.283
##  forceSymmetric(dgT) 137.719 149.5475 157.96972 156.1280 162.6265  217.997
##            triu(dgT) 162.811 184.9510 193.10877 190.7115 199.3830  229.723
##            tril(dgT) 187.411 201.0230 211.25783 208.1365 216.5620  345.179
##  neval
##    100
##    100
##    100
##    100
##    100
##    100
##    100
##    100
##    100


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TEST FOR SYMMETRIC .g[CRT]Matrix
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

## Matrix r3504 ... FIXME: why slower for expressions 2 and 3?

## Unit: microseconds
##                            expr      min        lq       mean    median
##           isSymmetric(dsC.as.g)  892.488  942.7540 1053.34002  949.6420
##           isSymmetric(dsR.as.g) 1519.296 1607.7945 1780.47379 1622.4930
##           isSymmetric(dsT.as.g) 1465.996 1530.6325 1752.25103 1543.3425
##  isSymmetric(dsC.as.g, tol = 0)   12.259   13.5095   15.38976   15.1495
##  isSymmetric(dsR.as.g, tol = 0)   12.218   15.1495   16.02567   16.1540
##  isSymmetric(dsT.as.g, tol = 0)  116.153  125.4190  144.14575  133.6805
##         uq      max neval
##   963.2745 6053.035   100
##  1638.8110 7782.538   100
##  1574.1335 7605.910   100
##    16.1540   32.472   100
##    16.8100   21.156   100
##   141.4090 1097.201   100

## Matrix 1.4-1

## Unit: microseconds
##                            expr      min       lq     mean   median       uq
##           isSymmetric(dsC.as.g)  982.442 1032.462 1109.879 1044.249 1062.023
##           isSymmetric(dsR.as.g)  978.465 1034.676 1101.170 1050.030 1069.321
##           isSymmetric(dsT.as.g) 1216.962 1270.016 1362.696 1283.669 1304.764
##  isSymmetric(dsC.as.g, tol = 0) 1003.721 1034.943 1130.826 1044.885 1065.549
##  isSymmetric(dsR.as.g, tol = 0)  984.205 1033.918 1049.640 1047.181 1061.305
##  isSymmetric(dsT.as.g, tol = 0) 1204.047 1268.991 1337.705 1282.870 1295.784
##       max neval
##  4095.490   100
##  5530.818   100
##  4951.119   100
##  4825.003   100
##  1225.080   100
##  3924.028   100


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TEST FOR TRIANGULAR .g[CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dtC.as.g <- as(dtC, "generalMatrix")
dtR.as.g <- as(dtR, "generalMatrix")
dtT.as.g <- as(dtT, "generalMatrix")

BM(isTriangular(dtC.as.g, upper = TRUE),
   isTriangular(dtR.as.g, upper = TRUE),
   isTriangular(dtT.as.g, upper = TRUE),
   unit = "microseconds")

## Matrix r3504

## Unit: microseconds
##                                  expr   min    lq    mean median    uq
##  isTriangular(dtC.as.g, upper = TRUE) 3.690 3.731 3.97946 3.7925 3.854
##  isTriangular(dtR.as.g, upper = TRUE) 3.772 3.813 4.02169 3.8540 3.895
##  isTriangular(dtT.as.g, upper = TRUE) 2.870 2.911 3.08238 2.9520 2.993
##     max neval
##  13.202   100
##  12.054   100
##  12.177   100

## Matrix 1.4-1

## Unit: microseconds
##                                  expr     min      lq      mean   median
##  isTriangular(dtC.as.g, upper = TRUE) 856.367 878.507 938.48016 885.3950
##  isTriangular(dtR.as.g, upper = TRUE) 858.048 878.671 900.17099 889.1465
##  isTriangular(dtT.as.g, upper = TRUE)  16.646  17.835  18.72511  18.3270
##        uq      max neval
##  919.9785 4123.903   100
##  919.8145 1019.465   100
##   19.1675   37.638   100


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## TEST FOR DIAGONAL .g[CRT]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ddC.as.g <- as(ddC, "generalMatrix")
ddR.as.g <- as(ddR, "generalMatrix")
ddT.as.g <- as(ddT, "generalMatrix")

BM(isDiagonal(ddC.as.g),
   isDiagonal(ddR.as.g),
   isDiagonal(ddT.as.g),
   unit = "microseconds")

## Matrix r3504

## Unit: microseconds
##                  expr   min    lq    mean median    uq     max neval
##  isDiagonal(ddC.as.g) 1.599 1.599 2.55471  1.640 1.681  91.102   100
##  isDiagonal(ddR.as.g) 1.599 1.599 3.66212  1.640 1.681 202.171   100
##  isDiagonal(ddT.as.g) 0.820 0.861 1.96513  0.861 0.861 108.568   100

## Matrix 1.4-1

## Unit: microseconds
##                  expr    min      lq     mean  median      uq     max neval
##  isDiagonal(ddC.as.g) 13.120 14.0015 16.00230 14.4115 14.7190 179.006   100
##  isDiagonal(ddR.as.g) 21.894 23.4110 29.74509 23.7390 23.9850 593.639   100
##  isDiagonal(ddT.as.g)  1.845  2.0500  4.38454  2.2140  2.3985 217.177   100


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DETERMINANT OF POSITIVE DEFINITE ds[yp]Matrix
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dsy <- new("dsyMatrix", Dim = c(1000L, 1000L), x = rnorm(1e+06L))
dsp <- pack(dsy)
dsy.spd <- as(tcrossprod(dsy), "dsyMatrix")
dsp.spd <- pack(dsy.spd)

BM(determinant(dsy, TRUE), determinant(dsy.spd, TRUE),
   determinant(dsp, TRUE), determinant(dsp.spd, TRUE),
   unit = "microseconds")

## Matrix r3504

## Unit: microseconds
##                        expr        min          lq       mean      median
##      determinant(dsy, TRUE) 115013.200 115702.9840 118651.318 116976.1365
##  determinant(dsy.spd, TRUE)     13.776     15.2110   1624.231     28.8640
##      determinant(dsp, TRUE) 114556.911 115311.6595 117812.106 116318.0250
##  determinant(dsp.spd, TRUE)     13.981     14.8625   1689.809     23.4315
##          uq      max neval
##  118051.464 159599.3   100
##      43.706 159322.0   100
##  117379.289 159832.3   100
##      40.713 166178.1   100

## Matrix 1.4-1

## Unit: microseconds
##                        expr      min       lq     mean   median       uq
##      determinant(dsy, TRUE) 115330.9 117101.9 118136.5 117468.7 118258.4
##  determinant(dsy.spd, TRUE) 114112.4 116349.9 118703.5 116714.8 117830.7
##      determinant(dsp, TRUE) 116177.7 117460.8 118898.5 117771.1 118693.6
##  determinant(dsp.spd, TRUE) 114924.1 116493.3 118708.7 116991.7 117789.3
##       max neval
##  160809.6   100
##  160231.2   100
##  158672.1   100
##  159442.4   100
