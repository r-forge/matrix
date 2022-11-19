## METHODS FOR GENERIC: [<-                  UNFINISHED AND NOT-YET-USED
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## GOAL: automate method definitions and eventually replace ones in
##
##             ./Csparse.R
##             ./Matrix.R
##             ./Rsparse.R
##             ./Tsparse.R
##             ./denseMatrix.R
##             ./diagMatrix.R
##             ./indMatrix.R
##             ./sparseMatrix.R
##
##       need to write C-level functions
##
##             *_subassign_1ary    (x, i,    value)
##             *_subassign_1ary_mat(x, i,    value)
##             *_subassign_2ary    (x, i, j, value)
##
##       for * = unpackedMatrix,packedMatrix,
##               CsparseMatrix,RsparseMatrix,TsparseMatrix,
##               diagonalMatrix,indMatrix
