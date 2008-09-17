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
##
## .Call() arguments:
Csparse_diagU2N  <- Matrix:::Csparse_diagU2N
Csparse_drop     <- Matrix:::Csparse_drop
Csparse_general_to_symmetric <- Matrix:::Csparse_general_to_symmetric
Matrix_expand_pointers <- Matrix:::Matrix_expand_pointers
dense_to_Csparse   <- Matrix:::dense_to_Csparse
dense_to_symmetric <- Matrix:::dense_to_symmetric
Tsparse_to_Csparse <- Matrix:::Tsparse_to_Csparse
Tsparse_diagU2N    <- Matrix:::Tsparse_diagU2N
Csparse_symmetric_to_general <- Matrix:::Csparse_symmetric_to_general
diag_tC    <- Matrix:::diag_tC
dup_mMatrix_as_dgeMatrix <- Matrix:::dup_mMatrix_as_dgeMatrix
dup_mMatrix_as_geMatrix  <- Matrix:::dup_mMatrix_as_geMatrix
dsCMatrix_matrix_solve  <-  Matrix:::dsCMatrix_matrix_solve
dsyMatrix_matrix_solve  <-  Matrix:::dsyMatrix_matrix_solve
Csparse_to_nz_pattern   <-  Matrix:::Csparse_to_nz_pattern
Ops.x.x 	<- Matrix:::Ops.x.x
attr.all_Mat	<- Matrix:::attr.all_Mat
attrSlots    	<- Matrix:::attrSlots
attrSlotNames 	<- Matrix:::attrSlotNames
m_encodeInd 	<- Matrix:::m_encodeInd
m_encodeInd2 	<- Matrix:::m_encodeInd2
