### Assign Matrix-hidden utilities to .GlobalEnv

### -- this makes  ONLY  sense when debugging;
### i.e. evaluating fixed things in .GlobalEnv

is0       <- Matrix:::is0
isN0      <- Matrix:::isN0
.M.kind   <- Matrix:::.M.kind
.M.kindC  <- Matrix:::.M.kindC
.M.DN     <- Matrix:::.M.DN
##
isDiagonal   <- Matrix:::isDiagonal
isTriangular <- Matrix:::isTriangular
## .Call() arguments:
dense_to_Csparse <- Matrix:::dense_to_Csparse
dense_to_symmetric <- Matrix:::dense_to_symmetric
Csparse_diagU2N  <- Matrix:::Csparse_diagU2N
Csparse_general_to_symmetric <- Matrix:::Csparse_general_to_symmetric
