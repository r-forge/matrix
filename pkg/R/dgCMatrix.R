#### Sparse Matrices in Compressed column-oriented format

### contains = "dsparseMatrix", "CsparseMatrix"

## Specific conversions, should they be necessary.  Better to convert as
## as(x, "TsparseMatrix") or as(x, "denseMatrix")

## Moved to ./Csparse.R :
## setAs("dgCMatrix", "dgTMatrix", ....
## setAs("dgCMatrix", "dgeMatrix", ....
## setAs("dgeMatrix", "dgCMatrix", ....

## rather use Csparse* to lsparse* in ./lsparseMatrix.R ,
## but this is for "back-compatibility" (have had tests for it..):

setAs("dgCMatrix", "ngCMatrix",
      function(from) .Call(Csparse_to_nz_pattern, from, FALSE))

setAs("dgCMatrix", "lgCMatrix",
      function(from) { ## FIXME use .Call() too!
	  r <- new("lgCMatrix")
	  r@x <- as.logical(from@x)
	  ## and copy the other slots
	  for(nm in c("i", "p", "Dim", "Dimnames"))
	      slot(r, nm) <- slot(from, nm)
	  r
      })

setMethod("image", "dgCMatrix",
	  function(x, ...) {
	      x <- as(x, "dgTMatrix")
	      callGeneric()
	  })

## Group Methods, see ?Arith (e.g.)
## -----
##
## "Arith" is now in ./Ops.R
##
## "Math" is up in ./Csparse.R
##
## "Math2" is up in ./dMatrix.R


###---- end {Group Methods} -----------------


## "[<-" methods { setReplaceMethod()s }  are now in ./Csparse.R

## setMethod("writeHB", signature(obj = "dgCMatrix"),
## 	  function(obj, file, ...) {
## 	      .Deprecated("writeMM")
## 	      .Call(Matrix_writeHarwellBoeing, obj,
## 		    as.character(file), "DGC")
## 	  })

##-> ./colSums.R  for colSums,... rowMeans

setMethod("determinant", signature(x = "dgCMatrix", logarithm = "logical"),
          detSparseLU)

## SPQR currently experimental, but possibly THE method with potential
## ----
## From spqr_user_guide.pdf :
## ordering : string describing what ordering method to use. Let [m2 n2]=size(S)
## where S is obtained by removing singletons from A. The singleton permutation places
## A*P in the form [A11 A12 ; 0 S] where A11 is upper triangular with diagonal entries
## all greater than tol.
## The default is to use COLAMD if m2<=2*n2; otherwise try AMD. Let f be the flops for
## chol((S*P)’*(S*P)) with the ordering P found by AMD. Then if f/nnz(R) >= 500
## and nnz(R)/nnz(S) >= 5 then try METIS, and take the best ordering found (AMD
## or METIS); otherwise use AMD without trying METIS. If METIS is not installed then
## the default ordering is to use COLAMD if m2<=2*n2 and to use AMD otherwise.
## The available orderings are:

## ’default’: the default ordering.
## ’amd’: use amd(S’*S).
## ’colamd’: use colamd(S).
## ’metis’: use metis(S’*S), only if METIS is installed.
## ’best’: try all three (AMD, COLAMD, METIS) and take the best.
## ’bestamd’: try AMD and COLAMD and take the best.
## ’fixed’: use P=I; this is the only option if P is not present in the output.
## ’natural’: singleton removal only.

#define SPQR_ORDERING_FIXED 0
#define SPQR_ORDERING_NATURAL 1
#define SPQR_ORDERING_COLAMD 2
#define SPQR_ORDERING_GIVEN 3       /* only used for C/C++ interface */
#define SPQR_ORDERING_CHOLMOD 4     /* CHOLMOD best-effort (COLAMD, METIS,...)*/
#define SPQR_ORDERING_AMD 5         /* AMD(A'*A) */
#define SPQR_ORDERING_METIS 6       /* metis(A'*A) */
#define SPQR_ORDERING_DEFAULT 7     /* SuiteSparseQR default ordering */
#define SPQR_ORDERING_BEST 8        /* try COLAMD, AMD, and METIS; pick best */
#define SPQR_ORDERING_BESTAMD 9     /* try COLAMD and AMD; pick best */

## Using spqr() is as Tim Davis has named his Matlab interface;
## but the function arguments will differ anyway
setGeneric("spqr", function(x, ...) standardGeneric("spqr"))
setMethod("spqr", signature(x = "dgCMatrix"),
	  function(x, econ = 0,
                   ordering = c("default", "fixed","natural","amd",
                                "colamd", "metis", "best", "bestamd"),
                   tol = -2)
      {
          ordering <- as.integer(
                                 c("fixed" = 0,
                                   "natural" = 1,
                                   "colamd" = 2,
                                   "given" = 3,
                                   "cholmod" = 4,
                                   "amd" = 5    ,
                                   "metis" = 6  ,
                                   "default" = 7,
                                   "best" = 8   ,
                                   "bestamd" = 9)[match.arg(ordering)])
          .Call(dgCMatrix_SPQR, x, ordering, econ)
      })

## an alternative would be (a version of) the following interface:
if(FALSE)
setMethod("qr", signature(x = "dgCMatrix"),
	  function(x, tol = 1e-07, LAPACK = FALSE, SPQR = FALSE)
          ## SPQR currently experimental, but possibly THE method with potential
          if(SPQR).Call(dgCMatrix_SPQR, x, ncol(x)) else
	  .Call(dgCMatrix_QR, x, TRUE))

setMethod("qr", signature(x = "dgCMatrix"),
	  function(x, tol = 1e-07, LAPACK = FALSE)
	  .Call(dgCMatrix_QR, x, TRUE))

setMethod("qr", signature(x = "sparseMatrix"),
	  function(x, ...)
	  qr(as(as(x, "CsparseMatrix"), "dsparseMatrix"), ...))

setMethod("lu", signature(x = "dgCMatrix"),
	  function(x, ...) {
	      .Call(dgCMatrix_LU, x,
			 TRUE, ## <- orderp
			 1) ## <- tol
	      })
setMethod("lu", signature(x = "sparseMatrix"),
	  function(x, ...) lu(as(as(x, "CsparseMatrix"), "dsparseMatrix"), ...))


setMethod("solve", signature(a = "dgCMatrix", b = "matrix"),
	  function(a, b, ...) .Call(dgCMatrix_matrix_solve, a, b),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgCMatrix", b = "ddenseMatrix"),
	  function(a, b, ...) .Call(dgCMatrix_matrix_solve, a, b),
	  valueClass = "dgeMatrix")
setMethod("solve", signature(a = "dgCMatrix", b = "dsparseMatrix"),
	  function(a, b, ...)
	  .Call(dgCMatrix_matrix_solve, a, as(b, "denseMatrix")),
	  valueClass = "dgeMatrix")

setMethod("solve", signature(a = "dgCMatrix", b = "missing"),
	  function(a, b, ...) .Call(dgCMatrix_matrix_solve, a, b=NULL),
	  valueClass = "dgeMatrix")
