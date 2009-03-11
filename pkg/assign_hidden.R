### Assign Matrix-hidden utilities to .GlobalEnv

### -- this makes  ONLY  sense when debugging;
### i.e. evaluating fixed things in .GlobalEnv

is0       <- Matrix:::is0
isN0      <- Matrix:::isN0
.M.kind   <- Matrix:::.M.kind
.M.kindC  <- Matrix:::.M.kindC
.M.shape   <- Matrix:::.M.shape
.M.DN     <- Matrix:::.M.DN
##
dimCheck  <- Matrix:::dimCheck

isDiagonal   <- Matrix:::isDiagonal
isTriangular <- Matrix:::isTriangular
##
## .Call() arguments:
dpoMatrix_chol <- Matrix:::dpoMatrix_chol
compressed_non_0_ij <- Matrix:::compressed_non_0_ij
Csparse_diagU2N  <- Matrix:::Csparse_diagU2N
Csparse_drop     <- Matrix:::Csparse_drop
Csparse_crossprod <- Matrix:::Csparse_crossprod
Csparse_transpose <- Matrix:::Csparse_transpose
Csparse_general_to_symmetric <- Matrix:::Csparse_general_to_symmetric
Matrix_expand_pointers <- Matrix:::Matrix_expand_pointers
dense_to_Csparse    <- Matrix:::dense_to_Csparse
dense_to_symmetric  <- Matrix:::dense_to_symmetric
compressed_to_TMatrix<- Matrix:::compressed_to_TMatrix
R_to_CMatrix  <- Matrix:::R_to_CMatrix
Tsparse_to_Csparse  <- Matrix:::Tsparse_to_Csparse
Tsparse_to_tCsparse <- Matrix:::Tsparse_to_tCsparse
Csparse_to_Tsparse  <- Matrix:::Csparse_to_Tsparse
Csparse_to_matrix   <- Matrix:::Csparse_to_matrix
dgTMatrix_to_matrix <- Matrix:::dgTMatrix_to_matrix
dgTMatrix_to_dgeMatrix <- Matrix:::dgTMatrix_to_dgeMatrix
Csparse_to_dense    <- Matrix:::Csparse_to_dense
Tsparse_diagU2N     <- Matrix:::Tsparse_diagU2N
Csparse_symmetric_to_general <- Matrix:::Csparse_symmetric_to_general
dense_band          <- Matrix:::dense_band
Csparse_band        <- Matrix:::Csparse_band
dsyMatrix_as_matrix <- Matrix:::dsyMatrix_as_matrix
dsyMatrix_as_dspMatrix <- Matrix:::dsyMatrix_as_dspMatrix
dspMatrix_as_dsyMatrix <- Matrix:::dspMatrix_as_dsyMatrix
dtrMatrix_as_matrix <- Matrix:::dtrMatrix_as_matrix
dsTMatrix_as_dsyMatrix <- Matrix:::dsTMatrix_as_dsyMatrix
dtTMatrix_as_dtrMatrix <- Matrix:::dtTMatrix_as_dtrMatrix
ltTMatrix_as_ltrMatrix <- Matrix:::ltTMatrix_as_ltrMatrix
lsyMatrix_as_lgeMatrix <- Matrix:::lsyMatrix_as_lgeMatrix
ltrMatrix_as_lgeMatrix <- Matrix:::ltrMatrix_as_lgeMatrix
ntTMatrix_as_ntrMatrix <- Matrix:::ntTMatrix_as_ntrMatrix

dMatrix_validate    <- Matrix:::dMatrix_validate
dgeMatrix_validate  <- Matrix:::dgeMatrix_validate
dsyMatrix_validate  <- Matrix:::dsyMatrix_validate
dpoMatrix_validate  <- Matrix:::dpoMatrix_validate
dppMatrix_validate  <- Matrix:::dppMatrix_validate
dtrMatrix_validate  <- Matrix:::dtrMatrix_validate
Csparse_validate    <- Matrix:::Csparse_validate
tCMatrix_validate   <- Matrix:::tCMatrix_validate
xCMatrix_validate   <- Matrix:::xCMatrix_validate
Tsparse_validate    <- Matrix:::Tsparse_validate
tTMatrix_validate   <- Matrix:::tTMatrix_validate
xTMatrix_validate   <- Matrix:::xTMatrix_validate
symmetricMatrix_validate  <- Matrix:::symmetricMatrix_validate
triangularMatrix_validate <- Matrix:::triangularMatrix_validate

dgeMatrix_getDiag <- Matrix:::dgeMatrix_getDiag
lgeMatrix_getDiag <- Matrix:::lgeMatrix_getDiag
diag_tC    <- Matrix:::diag_tC
dup_mMatrix_as_dgeMatrix <- Matrix:::dup_mMatrix_as_dgeMatrix
dup_mMatrix_as_geMatrix  <- Matrix:::dup_mMatrix_as_geMatrix
dsCMatrix_matrix_solve  <-  Matrix:::dsCMatrix_matrix_solve
dsyMatrix_matrix_solve  <-  Matrix:::dsyMatrix_matrix_solve
dtCMatrix_sparse_solve  <-  Matrix:::dtCMatrix_sparse_solve
Csparse_to_nz_pattern   <-  Matrix:::Csparse_to_nz_pattern

dgCMatrix_colSums <- Matrix:::dgCMatrix_colSums
lgCMatrix_colSums <- Matrix:::lgCMatrix_colSums
ngCMatrix_colSums <- Matrix:::ngCMatrix_colSums

## end{.Call} ------------------------------------------------------

Ops.x.x 	<- Matrix:::Ops.x.x
attr.all_Mat	<- Matrix:::attr.all_Mat
attrSlots    	<- Matrix:::attrSlots
attrSlotNames 	<- Matrix:::attrSlotNames
m_encodeInd 	<- Matrix:::m_encodeInd
m_encodeInd2 	<- Matrix:::m_encodeInd2
