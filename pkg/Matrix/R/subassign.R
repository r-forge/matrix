## METHODS FOR GENERIC: [<-
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
##             *_subassign_1ary_vec(x, i,    value)
##             *_subassign_1ary_mat(x, i,    value)
##             *_subassign_2ary    (x, i, j, value)
##
##       for * = CRsparse,Tsparse,unpackedMatrix,packedMatrix
##       and where 'value' is allowed to be a vector _or_ sparseVector
##
##       diagonalMatrix and indMatrix should go via CsparseMatrix
