## METHODS ENHANCING PACKAGE: sparseM
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ COERCIONS FROM SparseM TO Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Natural pairs
setAs("matrix.csc", "dgCMatrix",
      function(from) new("dgCMatrix", Dim = from@dimension,
                         p = from@ia - 1L, i = from@ja - 1L, x = from@ra))
setAs("matrix.csr", "dgRMatrix",
      function(from) new("dgRMatrix", Dim = from@dimension,
                         p = from@ia - 1L, j = from@ja - 1L, x = from@ra))
setAs("matrix.coo", "dgTMatrix",
      function(from) new("dgTMatrix", Dim = from@dimension,
                         i = from@ia - 1L, j = from@ja - 1L, x = from@ra))

## Remaining matrix.c(sc|sr|oo) to dgCMatrix
setAs("matrix.csr", "dgCMatrix",
      function(from) as(as(from, "dgRMatrix"), "CsparseMatrix"))
setAs("matrix.coo", "dgCMatrix",
      function(from) as(as(from, "dgTMatrix"), "CsparseMatrix"))

## Each matrix.c(sc|sr|oo) to each [CRT]sparseMatrix
setAs("matrix.csc", "CsparseMatrix",
      function(from)    as(from, "dgCMatrix"))
setAs("matrix.csc", "RsparseMatrix",
      function(from) as(as(from, "dgCMatrix"), "RsparseMatrix"))
setAs("matrix.csc", "TsparseMatrix",
      function(from) as(as(from, "dgCMatrix"), "TsparseMatrix"))
setAs("matrix.csr", "CsparseMatrix",
      function(from) as(as(from, "dgRMatrix"), "CsparseMatrix"))
setAs("matrix.csr", "RsparseMatrix",
      function(from)    as(from, "dgRMatrix"))
setAs("matrix.csr", "TsparseMatrix",
      function(from) as(as(from, "dgRMatrix"), "TsparseMatrix"))
setAs("matrix.coo", "CsparseMatrix",
      function(from) as(as(from, "dgTMatrix"), "CsparseMatrix"))
setAs("matrix.coo", "RsparseMatrix",
      function(from) as(as(from, "dgTMatrix"), "RsparseMatrix"))
setAs("matrix.coo", "TsparseMatrix",
      function(from)    as(from, "dgTMatrix"))

## Each matrix.c(sc|sr|oo) to sparseMatrix, Matrix ("easy")
## NB: favouring column format over row format
setAs("matrix.csc", "sparseMatrix", function(from) as(from, "CsparseMatrix"))
setAs("matrix.csr", "sparseMatrix", function(from) as(from, "CsparseMatrix"))
setAs("matrix.coo", "sparseMatrix", function(from) as(from, "TsparseMatrix"))
setAs("matrix.csc",       "Matrix", function(from) as(from, "CsparseMatrix"))
setAs("matrix.csr",       "Matrix", function(from) as(from, "CsparseMatrix"))
setAs("matrix.coo",       "Matrix", function(from) as(from, "TsparseMatrix"))


## ~~~~ COERCIONS FROM Matrix TO SparseM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Natural pairs
setAs("dgCMatrix", "matrix.csc",
      function(from) new("matrix.csc", dimension = from@Dim,
                         ia = from@p + 1L, ja = from@i + 1L, ra = from@x))
setAs("dgRMatrix", "matrix.csr",
      function(from) new("matrix.csr", dimension = from@Dim,
                         ia = from@p + 1L, ja = from@j + 1L, ra = from@x))
setAs("dgTMatrix", "matrix.coo",
      function(from) new("matrix.coo", dimension = from@Dim,
                         ia = from@i + 1L, ja = from@j + 1L, ra = from@x))

## Everything else
setAs("Matrix", "matrix.csc",
      function(from) as(as(as(as(from, "dMatrix"), "generalMatrix"), "CsparseMatrix"), "matrix.csc"))
setAs("Matrix", "matrix.csr",
      function(from) as(as(as(as(from, "dMatrix"), "generalMatrix"), "RsparseMatrix"), "matrix.csr"))
setAs("Matrix", "matrix.coo",
      function(from) as(as(as(as(from, "dMatrix"), "generalMatrix"), "TsparseMatrix"), "matrix.coo"))
