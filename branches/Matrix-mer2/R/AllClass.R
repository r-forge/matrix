.onLoad <- function(lib, pkg) {
    if(is.null(getOption("max.print")))
	options(max.print = 10000)#-> show() of large matrices
}

## ------------- Virtual Classes ----------------------------------------

## Mother class of all Matrix objects
setClass("Matrix",
	 representation(Dim = "integer", Dimnames = "list", factors = "list",
			"VIRTUAL"),
	 prototype = prototype(Dim = integer(2), Dimnames = list(NULL,NULL)),
	 validity = function(object) {
	     Dim <- object@Dim
	     if (length(Dim) != 2)
		 return("Dim slot must be of length 2")
	     if (any(Dim < 0))
		 return("Dim slot must contain non-negative values")
	     Dn <- object@Dimnames
	     if (!is.list(Dn) || length(Dn) != 2)
		 return("'Dimnames' slot must be list of length 2")
	     ## 'else'	ok :
	     TRUE
	 })

setClass("symmetricMatrix",
	 representation(uplo = "character", "VIRTUAL"),
	 contains = "Matrix")

setClass("triangularMatrix",
	 representation(uplo = "character", diag = "character",
			"VIRTUAL"), contains = "Matrix")

## Virtual class of numeric matrices
setClass("dMatrix",
	 representation(x = "numeric", "VIRTUAL"), contains = "Matrix")

## Virtual class of integer matrices
setClass("iMatrix",
	 representation(x = "integer", "VIRTUAL"), contains = "Matrix")

## Virtual class of logical matrices
setClass("lMatrix", representation("VIRTUAL"), contains = "Matrix")
## Note that logical sparse matrices do not need an x slot so the x
## slot is part of the ldenseMatrix class

## Virtual class of complex matrices
setClass("zMatrix", # letter 'z' is as in the names of Lapack subroutines
	 representation(x = "complex", "VIRTUAL"), contains = "Matrix")

## Virtual class of dense matrices
setClass("denseMatrix", representation("VIRTUAL"),
         contains = "Matrix")

## Virtual class of dense, numeric matrices
setClass("ddenseMatrix", representation(rcond = "numeric", "VIRTUAL"),
	 contains = c("dMatrix", "denseMatrix"))

## Virtual class of dense, logical matrices
setClass("ldenseMatrix", representation(x = "logical", "VIRTUAL"),
	 contains = c("lMatrix", "denseMatrix"))

## virtual SPARSE ------------

setClass("sparseMatrix", representation("VIRTUAL"), contains = "Matrix")

## sparse matrices in Triplet representation (dgT, lgT, ..):
setClass("TsparseMatrix", representation(i = "integer", j = "integer", "VIRTUAL"),
         contains = "sparseMatrix")

setClass("CsparseMatrix", representation(i = "integer", p = "integer", "VIRTUAL"),
         contains = "sparseMatrix")

setClass("RsparseMatrix", representation(p = "integer", j = "integer", "VIRTUAL"),
         contains = "sparseMatrix")

setClass("dsparseMatrix", representation("VIRTUAL"),
	 contains = c("dMatrix", "sparseMatrix"))

setClass("lsparseMatrix", representation("VIRTUAL"),
	 contains = c("lMatrix", "sparseMatrix"))

## ------------------ Proper (non-virtual) Classes ----------------------------

##----------------------  DENSE	 -----------------------------------------

## numeric, dense, general matrices
setClass("dgeMatrix", contains = "ddenseMatrix",
	 ## checks that length( @ x) == prod( @ Dim):
	 validity =
         function(object) .Call("dgeMatrix_validate", object, PACKAGE = "Matrix")
	 )
## i.e. "dgeMatrix" cannot be packed, but "ddenseMatrix" can ..

## numeric, dense, non-packed, triangular matrices
setClass("dtrMatrix",
         ## FIXME?
         ##> 'ddense*' before 'dge*' so it can use d* or ddense* methods
         ##> WITHOUT a coerce to dge* (losing triangularity)
         ##> gives error from callNextMethod() in crossprod() dispatch {R bug?}
	 ##> contains = c("ddenseMatrix", "dgeMatrix", "triangularMatrix"),
	 contains = c("dgeMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
         function(object) .Call("dtrMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, dense, packed, triangular matrices
setClass("dtpMatrix",
	 contains = c("ddenseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
         function(object) .Call("dtpMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, dense, non-packed symmetric matrices
setClass("dsyMatrix",
         ## FIXME?
         ##> 'ddense*' before 'dge*' so it can use d* or ddense* methods
         ##> WITHOUT a coerce to dge* (losing triangularity)
         ##> gives error in crossprod() dispatch
	 ##> contains = c("ddenseMatrix", "dgeMatrix", "symmetricMatrix"),
	 contains = c("dgeMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U"),
	 validity =
         function(object) .Call("dsyMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, dense, packed symmetric matrices
setClass("dspMatrix",
	 prototype = prototype(uplo = "U"),
	 contains = c("ddenseMatrix", "symmetricMatrix"),
	 validity =
         function(object) .Call("dspMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, dense, non-packed, positive-definite, symmetric matrices
setClass("dpoMatrix", contains = "dsyMatrix",
	 validity =
         function(object) .Call("dpoMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, dense, packed, positive-definite, symmetric matrices
setClass("dppMatrix", contains = "dspMatrix",
	 validity =
         function(object) .Call("dppMatrix_validate", object, PACKAGE = "Matrix")
	 )

##----- logical dense Matrices -- e.g. as result of <ddenseMatrix>  COMPARISON

## numeric, dense, general matrices
setClass("lgeMatrix", contains = "ldenseMatrix",
	 ## checks that length( @ x) == prod( @ Dim):
	 validity = function(object) stopifnot(length(object@x) == prod(object@Dim))
	 )
## i.e. "lgeMatrix" cannot be packed, but "ldenseMatrix" can ..

## numeric, dense, non-packed, triangular matrices
setClass("ltrMatrix",
	 contains = c("lgeMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"))

## numeric, dense, packed, triangular matrices
setClass("ltpMatrix",
	 contains = c("ldenseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N")
	 ## validity: ldense*, triangular*  should suffice
	 )

## numeric, dense, non-packed symmetric matrices
setClass("lsyMatrix",
	 contains = c("lgeMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U")
         ##, validity = function(object) .Call("lsyMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, dense, packed symmetric matrices
setClass("lspMatrix",
	 prototype = prototype(uplo = "U"),
	 contains = c("ldenseMatrix", "symmetricMatrix"),
         validity = function(object)
         .Call("dspMatrix_validate", object, PACKAGE = "Matrix")
         ## "dsp" and "lsp" have the same validate
	 )

##-------------------- S P A R S E (non-virtual) --------------------------

##---------- numeric sparse matrix classes --------------------------------

## numeric, sparse, triplet general matrices
setClass("dgTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix"),
	 validity =
         function(object) .Call("dgTMatrix_validate", object, PACKAGE = "Matrix")
	 )

## Should not have dtTMatrix inherit from dgTMatrix because a dtTMatrix could
## be less than fully stored if diag = "U".  Methods for the dgTMatrix
## class would not produce correct results even though all the slots
## are present.

## numeric, sparse, triplet triangular matrices
setClass("dtTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 prototype = prototype(uplo = "U", diag = "N"),
	 validity =
         function(object) .Call("dtTMatrix_validate", object, PACKAGE = "Matrix")
	 )

## Should not have dsTMatrix inherit from dgTMatrix because a dsTMatrix
## is not fully stored.  Methods for the dgTMatrix class would not
## produce correct results even though all the slots are present.

## numeric, sparse, triplet symmetric matrices
setClass("dsTMatrix",
	 contains = c("TsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(uplo = "U"),
	 validity =
         function(object) .Call("dsTMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, sparse, sorted compressed sparse column-oriented general matrices
setClass("dgCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix"),
	 prototype = prototype(p = 0:0),# to be valid
	 validity =
         function(object) .Call("dgCMatrix_validate", object, PACKAGE = "Matrix")
	 )

## see comments for dtTMatrix above
## numeric, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("dtCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "triangularMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U", diag = "N"),# to be valid
	 validity =
         function(object) .Call("tsc_validate", object, PACKAGE = "Matrix")
	 )

## see comments for dsTMatrix above
## numeric, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("dsCMatrix",
	 contains = c("CsparseMatrix", "dsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U"),# to be valid
	 validity =
         function(object) .Call("dsCMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, sparse, sorted compressed sparse row-oriented general matrices
setClass("dgRMatrix",
	 representation(j = "integer", p = "integer"),
	 contains = "dsparseMatrix",
	 ##TODO: validity = function(object) .Call("dgRMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("dtRMatrix",
	 contains = c("dgRMatrix", "triangularMatrix"),
	 ##TODO: validity = function(object) .Call("dtRMatrix_validate", object, PACKAGE = "Matrix")
	 )

## numeric, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("dsRMatrix",
	 contains = c("dgRMatrix", "symmetricMatrix"),
	 ##TODO: validity = function(object) .Call("dsRMatrix_validate", object, PACKAGE = "Matrix")
	 )

##---------- logical sparse matrix classes --------------------------------

## these classes are used in symbolic analysis to determine the
## locations of non-zero entries

## logical, sparse, triplet general matrices
setClass("lgTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix"),
	 validity =
         function(object) .Call("lgTMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, triplet triangular matrices
setClass("ltTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "triangularMatrix"),
	 validity =
         function(object) .Call("ltTMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, triplet symmetric matrices
setClass("lsTMatrix",
	 contains = c("TsparseMatrix", "lsparseMatrix", "symmetricMatrix"),
	 validity =
         function(object) .Call("lsTMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, sorted compressed sparse column-oriented general matrices
setClass("lgCMatrix",
	 contains = c("lsparseMatrix", "CsparseMatrix"),
	 prototype = prototype(p = 0:0),# to be valid
	 validity =
         function(object) .Call("lgCMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, sorted compressed sparse column-oriented triangular matrices
setClass("ltCMatrix",
	 contains = c("lsparseMatrix", "CsparseMatrix", "triangularMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U", diag = "N"),# to be valid
	 validity =
         function(object) .Call("ltCMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, sorted compressed sparse column-oriented symmetric matrices
setClass("lsCMatrix",
	 contains = c("lsparseMatrix", "CsparseMatrix", "symmetricMatrix"),
	 prototype = prototype(p = 0:0, uplo = "U"),# to be valid
	 validity =
         function(object) .Call("lsCMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, sorted compressed sparse row-oriented general matrices
setClass("lgRMatrix",
	 representation(j = "integer", p = "integer"),
	 contains = "lsparseMatrix",
	 validity =
         function(object) .Call("lgRMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, sorted compressed sparse row-oriented triangular matrices
setClass("ltRMatrix",
	 contains = c("lgRMatrix", "triangularMatrix"),
	 validity =
         function(object) .Call("ltRMatrix_validate", object, PACKAGE = "Matrix")
	 )

## logical, sparse, sorted compressed sparse row-oriented symmetric matrices
setClass("lsRMatrix",
	 contains = c("lgRMatrix", "symmetricMatrix"),
	 validity =
         function(object) .Call("lsRMatrix_validate", object, PACKAGE = "Matrix")
	 )

## Factorization classes

setClass("Cholesky", contains = "dtrMatrix")

setClass("LDL", contains = "dtrMatrix")

setClass("correlation", representation(sd = "numeric"), contains = "dpoMatrix")

setClass("pCholesky", contains = "dtpMatrix")

setClass("BunchKaufman",
         representation(perm = "integer"),
         contains = "dtrMatrix",
	 validity =
         function(object) .Call("BunchKaufman_validate", object, PACKAGE = "Matrix")
         )

setClass("pBunchKaufman",
         representation(perm = "integer"),
         contains = "dtpMatrix",
	 validity =
         function(object) .Call("pBunchKaufman_validate", object, PACKAGE = "Matrix")
         )

setClass("dCholCMatrix",
	 representation(perm = "integer", Parent = "integer", D = "numeric"),
	 contains = "dtCMatrix",
	 validity =
         function(object) .Call("dCholCMatrix_validate", object, PACKAGE = "Matrix")
         )

setClass("lCholCMatrix",
	 representation(perm = "integer", Parent = "integer"),
	 contains = "ltCMatrix",
	 validity =
         function(object) .Call("lCholCMatrix_validate", object, PACKAGE = "Matrix")
         )

##-------------------- permutation ----------------------------------------

setClass("pMatrix", representation(perm = "integer"),
         contains = "sparseMatrix",
	 validity = function(object) {
	     d <- object@Dim
	     if (d[2] != (n <- d[1])) return("pMatrix must be square")
	     perm <- object@perm
	     if (length(perm) != n)
		 return(paste("length of 'perm' slot must be", n))
	     if(n > 0 &&
		!(all(range(perm) == c(1, n)) && length(unique(perm)) == n))
		 return("'perm' slot is not a valid permutation")
	     TRUE
	 })

### Class Union :  no inheritance, but is(*, <class>) :
setClassUnion("packedMatrix",
              members = c("dspMatrix", "dppMatrix", "dtpMatrix",
               "lspMatrix", "ltpMatrix"))


## --------------------- non-"Matrix" Classes --------------------------------

## --- "General" (not Matrix at all) ----

## for 'i' in x[i] or A[i,] :
setClassUnion("index", members =  c("numeric", "logical", "character"))


## --- Matrix - related ----

setClass("determinant",
	 representation(modulus = "numeric",
			logarithm = "logical",
			sign = "integer",
			call = "call"))

setClass("LU",
         representation(x = "numeric", perm = "integer"),
	 validity = function(object) .Call("LU_validate", object, PACKAGE = "Matrix")
         )

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")
setOldClass("externalptr")

## mixed effects representation
setClass("mer",
	 representation(## original data
			flist = "list", # list of grouping factors
                        Zt = "dgCMatrix",  # sparse representation of Z'
                        X = "matrix",      # X
                        y = "numeric",     # y
			method = "character", # parameter estimation method
                        useScale = "logical", # should scale factor be included
                        family = "family", # glm family
                        call = "call",     # call to model-fitting function
                        ## invariants derived from data structure
			cnames = "list",   # column names of model matrices
			nc = "integer",    # dimensions of blocks in Omega
			Gp = "integer",    # Pointers to groups of rows in Zt
			devComp = "numeric", # Components of deviance
			deviance = "numeric", # Current deviance (ML and REML)
			status = "logical", # factored? inverted?
                        ## quantities that vary when Z, X or y are updated
			XtX = "dpoMatrix", # X'X
			ZtZ = "dsCMatrix", # Z'Z
			ZtX = "dgeMatrix", # Z'X
                        Zty = "numeric",   # Z'y 
                        Xty = "numeric",   # X'y 
                        ## slots that vary during the optimization
			Omega = "list", # list of relative precision matrices
                        ## Cholesky factor of inflated [Z:X:y]'[Z:X:y]
                        L = "externalptr", # sparse Cholesky factor of Z'Z + Omega
			RZX = "dgeMatrix",
			RXX = "dtrMatrix",
                        rZy = "numeric",
                        rXy = "numeric",
                        ## lists of arrays
                        bVar = "list",
                        gradComp = "list"
			)
	)

## Representation of a linear or generalized linear mixed effects model
setClass("lmer",
	 representation(assign = "integer", fitted = "numeric",
			fixed = "numeric", frame = "data.frame",
			logLik = "logLik", residuals = "numeric",
			terms = "terms"),
	 contains = "mer")

setClass("lmer,ranef", contains = "list")


