#### Logical Sparse Matrices in Compressed column-oriented format

### contains = "lsparseMatrix"

setMethod("%*%", signature(x = "lgCMatrix", y = "lgCMatrix"),
          function(x, y) .Call("lgCMatrix_lgCMatrix_mm", x, y),
          valueClass = "lgCMatrix")

setMethod("t", signature(x = "lgCMatrix"),
          function(x) .Call("lgCMatrix_trans", x),
          valueClass = "lgCMatrix")

