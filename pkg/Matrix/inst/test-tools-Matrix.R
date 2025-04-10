#### Tools for Package Testing --- in Matrix, sourced by ./test-tools.R
#### -------------------------

### ------- Part III --  "Matrix" (classes) specific ----------------------

## lower.tri() and upper.tri()  -- masking  base definitions
##	R/src/library/base/R/lower.tri.R
##	R/src/library/base/R/upper.tri.R
## but we do __not__ want to coerce to "base R" 'matrix' via as.matrix():
##
lower.tri <- function(x, diag = FALSE) if(diag) row(x) >= col(x) else row(x) > col(x)
upper.tri <- function(x, diag = FALSE) if(diag) row(x) <= col(x) else row(x) < col(x)

lsM <- function(...) {
    for(n in ls(..., envir=parent.frame()))
        if(is((. <- get(n)),"Matrix"))
            cat(sprintf("%5s: '%s' [%d x %d]\n",n,class(.), nrow(.),ncol(.)))
}

asD <- function(m) { ## as "Dense"
    if(canCoerce(m, "denseMatrix")) as(m, "denseMatrix")
    else if(canCoerce(m, (cl <- paste(.M.kind(m), "denseMatrix", sep=''))))
        as(m, cl)
    else if(canCoerce(m, "dgeMatrix")) as(m, "dgeMatrix")
    else stop("cannot coerce to a typical dense Matrix")
}

## "normal" sparse Matrix: Csparse, no diag="U"
asCsp <- function(x) diagU2N(as(x, "CsparseMatrix"))

##' @title quasi-identical dimnames
Qidentical.DN <- function(dx, dy) {

    stopifnot(is.list(dx) || is.null(dx),
	      is.list(dy) || is.null(dy))
    ## "empty"
    (is.null.DN(dx) && is.null.DN(dy)) || identical(dx, dy)
}

##' quasi-identical()  for 'Matrix' matrices
Qidentical <- function(x,y, strictClass = TRUE) {
    if(!identical(class(x), cy <- class(y))) {
        if(strictClass || !is(x, cy))
           return(FALSE)
    } ## else try further
    ## if(identical(x, y))
    ##     return(TRUE)
    ## ## else try further
    slts <- slotNames(x) ## MJ: should be slotNames(y), since is(x, class(y)) ??
    if("Dimnames" %in% slts) { ## always (or we have no 'Matrix')
	slts <- slts[slts != "Dimnames"]
	if(!Qidentical.DN(x@Dimnames, y@Dimnames) &&
	   ## allow for "completion" of (NULL, <names>) dimnames of symmetricMatrix:
	   !Qidentical.DN(dimnames(x), dimnames(y)))
	    return(FALSE)
    }
    if("factors" %in% slts) { ## allow one empty and one non-empty 'factors'
        slts <- slts[slts != "factors"]
        ## if both are not empty, they must be the same:
        if(length(xf <- x@factors) && length(yf <- y@factors))
            if(!identical(xf, yf)) return(FALSE)
    }
    for(sl in slts)
        if(!identical(slot(x,sl), slot(y,sl)))
            return(FALSE)
    TRUE
}

## MJ: It seems intuitive to allow either of is(x, class(y))
##     and is(y, class(x)) when strictClass=FALSE ...
.MJ.Qidentical <- function(x, y, strictClass = TRUE, skipSlots = NULL) {
    isxy <- identical(cx <- class(x), cy <- class(y))
    if (!isxy) {
        if (strictClass)
            return(FALSE)
        isxy <- is(x, cy)
        if (!(isxy || is(y, cx)))
            return(FALSE)
        ## else try further
    }
    slts <- slotNames(if (isxy) y else x)
    if (length(skipSlots))
        slts <- setdiff(slts, skipSlots)
    if ("Dimnames" %in% slts) { ## always, or we have no "Matrix"
	## allow symmetrization of 'Dimnames' for "symmetricMatrix" :
	slts <- slts[slts != "Dimnames"]
        if(!Qidentical.DN(x@Dimnames, y@Dimnames) &&
	   !Qidentical.DN(dimnames(x), dimnames(y)))
	    return(FALSE)
    }
    if ("factors" %in% slts) {
        ## allow one empty and one non-empty 'factors' :
        slts <- slts[slts != "factors"]
        ## if both are not empty, they must be the same:
        if (length(xf <- x@factors) && length(yf <- y@factors) &&
           !identical(xf, yf))
            return(FALSE)
    }
    for (slt in slts)
        if (!identical(slot(x, slt), slot(y, slt)))
            return(FALSE)
    TRUE
}

##' quasi-identical()  for traditional ('matrix') matrices
mQidentical <- function(x,y, strictClass = TRUE) {
    if(!identical(class(x), cy <- class(y))) {
        if(strictClass || !is(x, cy))
            return(FALSE)
        ## else try further
    }
    if(!Qidentical.DN(dimnames(x), dimnames(y)))
        return(FALSE)
    identical(unname(x), unname(y))
}

Q.C.identical <- function(x,y, sparse = is(x,"sparseMatrix"),
                          checkClass = TRUE, strictClass = TRUE) {
    if(checkClass && class(x) != class(y)) {
        if(strictClass || !is(x, class(y)))
	   return(FALSE) ## else try further
    }
    if(sparse)
	Qidentical(as(x,"CsparseMatrix"), as(y,"CsparseMatrix"),
		   strictClass=strictClass)
    else Qidentical(x,y, strictClass=strictClass)
}

##' <description>
##'
##' <details>
##' @title  Quasi-equal for  'Matrix' matrices
##' @param x  Matrix
##' @param y  Matrix
##' @param superclasses  x and y must coincide in (not) extending these; set to empty,
##'  if no class/inheritance checks should happen.
##' @param dimnames.check  logical indicating if dimnames(.) much match
##' @param tol  NA (--> use "==") or numerical tolerance for all.equal()
##' @return   logical:  Are x and y (quasi) equal ?
Q.eq <- function(x, y,
		 superclasses =
		 c("sparseMatrix", "denseMatrix",
		   "dMatrix", "lMatrix", "nMatrix", "iMatrix", "zMatrix"),
		 dimnames.check = TRUE, tol = NA) {
    ## quasi-equal - for 'Matrix' matrices
    if(any(dim(x) != dim(y)))
	return(FALSE)
    if(dimnames.check &&
       !identical(dimnames(x),
		  dimnames(y))) return(FALSE)
    xcl <- getClassDef(class(x))
    ycl <- getClassDef(class(y))
    for(SC in superclasses) {
	if( extends(xcl, SC) &&
	   !extends(ycl, SC)) return(FALSE)
    }
    kind <- Matrix::: .Arith.kind(x, y, "+") # for asC() :
    asC <- ## asCommon
        if((isDense <- extends(xcl,"denseMatrix")))
            function(m) as(m, "matrix")
        else function(m) {
            as(as(as(m,"CsparseMatrix"), paste0(kind, "Matrix")), "generalMatrix") # => "{d,z,i}gC"
        }
    if(is.na(tol)) {
	if(isDense)
	    all(x == y | (is.na(x) & is.na(y)))
	else ## 'x == y' blows up for large sparse matrices:
	    isTRUE(all.equal(asC(x), asC(y), tolerance = 0.,
			     check.attributes = dimnames.check))
    }
    else if(is.numeric(tol) && tol >= 0) {
        isTRUE(all.equal(asC(x), asC(y), tolerance = tol,
                         check.attributes = dimnames.check))
    }
    else stop("'tol' must be NA or non-negative number")
} # Q.eq()

## Here, 'kind' (d, l, z, i, n) does not matter
Q.eq2 <- function(x, y,
		  superclasses = c("sparseMatrix", "denseMatrix"),
		  dimnames.check = FALSE, tol = NA)
    Q.eq(x,y, superclasses=superclasses,
         dimnames.check=dimnames.check, tol=tol)

##' <description>
##'
##' <details>
##' @title  Quasi-equality of  symmpart(m) + skewpart(m) with m
##' @param m  Matrix
##' @param tol  numerical tolerance for all.equal()
##' @return logical
##' @author Martin Maechler
Q.eq.symmpart <- function(m, tol = 8 * .Machine$double.eps)
{
    ss <- symmpart(m) + skewpart(m)
    if(hasNA <- any(iNA <- is.na(ss))) {
	## ss has the NA's symmetrically, but typically m has *not*
	iiNA <- which(iNA) # <- useful! -- this tests  which() methods!
	## assign NA's too -- using correct kind of NA:
	m[iiNA] <- as(NA, Matrix:::.type.kind[Matrix:::.M.kind(m)])
    }
    Q.eq2(m, ss, tol = tol)
}

##' sample.int(n, size, replace=FALSE) for really large n:
sampleL <- function(n, size) {
    if(n < .Machine$integer.max)
	sample.int(n, size)
    else {
	i <- unique(round(n * runif(1.8 * size)))
	while(length(i) < size) {
	    i <- unique(c(i, round(n * runif(size))))
	}
	i[seq_len(size)]
    }
}


## Useful Matrix constructors for testing:

##' @title Random Sparse Matrix
##' @param n
##' @param m number of columns; default (=n)  ==> square matrix
##' @param density the desired sparseness density:
##' @param nnz number of non-zero entries; default from \code{density}
##' @param repr character string specifying the sparseness kind of the result.
##' @param giveCsparse *deprecated* logical specifying if result should be CsparseMatrix
##' @return a [CTR]sparseMatrix,  n x m
##' @author Martin Maechler, Mar 2008; July 2020 ('repr' instead og 'giveCsparse')
rspMat <- function(n, m = n, density = 1/4, nnz = round(density * n*m),
                   repr = c("C","T","R"), giveCsparse)
{
    stopifnot(length(n) == 1, n == as.integer(n),
              length(m) == 1, m == as.integer(m),
              0 <= density, density <= 1,
              0 <= nnz,
	      nnz <= (N <- n*m))
    in0 <- sampleL(N, nnz)
    x <- sparseVector(i = in0, x = as.numeric(1L + seq_along(in0)), length = N)
    dim(x) <- c(n,m)#-> sparseMatrix
    ## silent, back compatible (not yet warning about 'giveCsparse' deprecation):
    repr <- if(missing(repr) && !missing(giveCsparse))
		if(giveCsparse) "C" else "T"
	    else match.arg(repr)
    switch(repr,
	   "C" = as(x, "CsparseMatrix"),
	   "T" =    x,# TsparseMatrix
	   "R" = as(x, "RsparseMatrix"))
}

## __DEPRECATED__ !!
rSparseMatrix <- function(nrow, ncol, nnz,
			  rand.x = function(n) round(rnorm(nnz), 2), ...)
{
    stopifnot((nnz <- as.integer(nnz)) >= 0,
	      nrow >= 0, ncol >= 0, nnz <= nrow * ncol)
    .Deprecated("rsparsematrix")
    ##=========
    sparseMatrix(i = sample(nrow, nnz, replace = TRUE),
		 j = sample(ncol, nnz, replace = TRUE),
		 x = rand.x(nnz), dims = c(nrow, ncol), ...)
}


rUnitTri <- function(n, upper = TRUE, ...)
{
    ## Purpose: random unit-triangular sparse Matrix .. built from rspMat()
    ## ----------------------------------------------------------------------
    ## Arguments:  n: matrix dimension
    ##         upper: logical indicating if upper or lower triangular
    ##         ...  : further arguments passed to rspMat(), eg. 'density'
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  5 Mar 2008, 11:35

    r <- (if(upper) triu else tril)(rspMat(n, ...))
    ## make sure the diagonal is empty
    diag(r) <- 0
    r <- drop0(r)
    r@diag <- "U"
    r
}

##' Construct a nice (with exact numbers) random artificial  \eqn{A = L D L'}
##' decomposition with a sparse \eqn{n \times n}{n x n} matrix \code{A} of
##' density \code{density} and square root \eqn{D} determined by \code{d0}.
##'
##' If one of \code{rcond} or \code{condest} is true, \code{A} must be
##' non-singular, both use an \eqn{LU} decomposition requiring
##' non-singularity.
##' @title Make Nice Artificial   A = L D L'  (With Exact Numbers) Decomposition
##' @param n matrix dimension \eqn{n \times n}{n x n}
##' @param density ratio of number of non-zero entries to total number
##' @param d0 The sqrt of the diagonal entries of D default \code{10}, to be
##' \dQuote{different} from \code{L} entries.  More generally these can be negative
##' @param rcond logical indicating if \code{\link{rcond}(A, useInv=TRUE)}
##' should be returned which requires non-singular A and D.
##' @param condest logical indicating if \code{\link{condest}(A)$est}
##' should be returned which requires non-singular A and D.
##' @return list with entries A, L, d.half, D, ..., where A inherits from
##' class \code{"\linkS4class{symmetricMatrix}"} and should be equal to
##' \code{as(L \%*\% D \%*\% t(L), "symmetricMatrix")}.
##' @author Martin Maechler, Date: 15 Mar 2008
mkLDL <- function(n, density = 1/3,
                  d0 = 10, d.half = d0 * sample.int(n), # random permutation
                  rcond = (n < 99), condest = (n >= 100))
{
    stopifnot(n == round(n), density <= 1)
    n <- as.integer(n)
    stopifnot(n >= 1, is.numeric(d.half),
              length(d.half) == n)# no longer (2023-05-24): d.half >= 0
    L <- Matrix(0, n,n)
    nnz <- round(density * n*n)
    L[sample(n*n, nnz)] <- seq_len(nnz)
    L <- tril(L, -1L)
    diag(L) <- 1
### FIXME: allow  *negative* d.half[] entries!
    dh2 <- d.half^2
    non.sing <- sum(dh2 > 0) == n
    D <- Diagonal(x = dh2)
    A <- tcrossprod(L * rep(d.half, each=n))
    ## = as(L %*% D %*% t(L), "symmetricMatrix")
    list(A = A, L = L, d.half = d.half, D = D,
	 rcond.A = if (rcond  && non.sing) rcond(A, useInv=TRUE),
	 cond.A  = if(condest && non.sing) condest(A)$est)
}

eqDeterminant <- function(m1, m2, d1 = determinant(m1), d2 = determinant(m2),
                          NA.Inf.ok=FALSE, tol=.Machine$double.eps^0.5, ...)
{
    ##  determinant(.) ## logarithm = TRUE -- large |modulus|  <==> rank deficient
    d1m <- as.vector(d1$modulus)# dropping attribute
    d2m <- as.vector(d2$modulus)
    if((identical(d1m, -Inf) && identical(d2m, -Inf)) ||
       ## <==> det(m1) == det(m2) == 0, then 'sign' may even differ !
       (is.na(d1m) && is.na(d2m)))
	## if both are NaN or NA, we "declare" that's fine here
	return(TRUE)
    else if(NA.Inf.ok && ## first can be NA, second infinite:
            ## wanted: base::determinant.matrix() sometimes gives -Inf instead
            ## of NA,e.g. for matrix(c(0,NA,0,0,NA,NA,0,NA,0,0,1,0,0,NA,0,1), 4,4))
            is.na(d1m) && is.infinite(d2m)) return(TRUE)
    ## else
    if(is.infinite(d1m)) d1$modulus <- sign(d1m)* .Machine$double.xmax
    if(is.infinite(d2m)) d2$modulus <- sign(d2m)* .Machine$double.xmax
    ## now they are finite or *one* of them is NA/NaN, and all.equal() will tell so:
    all.equal(d1, d2, tolerance=tol, ...)
}

##' @param A a non-negative definite sparseMatrix, typically "dsCMatrix"
##'
##' @return a list with components resulting from calling
##'    Cholesky(.,  perm = .P., LDL = .L., super = .S.)
##'
##'    for all 2*2*3 combinations of (.P., .L., .S.)
allCholesky <- function(A, verbose = FALSE, silentTry = FALSE)
{
    ## Author: Martin Maechler, Date: 16 Jul 2009

    ##' @param r   list of CHMfactor objects, typically with names() as '. | .'
    ##'
    ##' @return an is(perm,LDL,super) matrix with interesting and *named* rownames
    CHM_to_pLs <- function(r) {
        is.perm  <- function(.)
            if(inherits(., "try-error")) NA else .@ordering != 0L
        is.LDL   <- function(.)
            if(inherits(., "try-error")) NA else .hasSlot(x, "is_ll") && !x@is_ll
        is.super <- function(.)
            if(inherits(., "try-error")) NA else .hasSlot(x, "super")
        r.st <- cbind(perm  = sapply(r, is.perm),
                      LDL   = sapply(r, is.LDL),
                      super = sapply(r, is.super))
        names(dimnames(r.st)) <- list("  p L s", "")
        r.st
    }

    my.Cholesky <- {
	if(verbose)
	    function (A, perm = TRUE, LDL = !super, super = FALSE, Imult = 0, ...) {
		cat(sprintf("Chol..(*, perm= %1d, LDL= %1d, super=%1d):",
			    perm, LDL, super))
		r <- Cholesky(A, perm=perm, LDL=LDL, super=super, Imult=Imult, ...)
		cat(" [Ok]\n")
		r
	    }
	else Cholesky
    }
    logi <- c(FALSE, TRUE)
    d12 <- expand.grid(perm = logi, LDL = logi, super = c(logi,NA),
		       KEEP.OUT.ATTRS = FALSE)
    r1 <- lapply(seq_len(nrow(d12)),
		 function(i) try(do.call(my.Cholesky,
                                         c(list(A = A), as.list(d12[i,]))),
                                 silent=silentTry))
    names(r1) <- apply(d12, 1,
		       function(.) paste(symnum(.), collapse=" "))
    dup.r1 <- duplicated(r1)
    r.all <- CHM_to_pLs(r1)
    if(!identical(dup.r1, duplicated(r.all)))
        warning("duplicated( <pLs-matrix> ) differs from duplicated( <CHM-list> )",
                immediate. = TRUE)
    list(Chol.A = r1,
         dup.r.all = dup.r1,
	 r.all	= r.all,
	 r.uniq = CHM_to_pLs(r1[ ! dup.r1]))
}

##' Cheap  Boolean Arithmetic Matrix product
##' Should be equivalent to  %&%  which is faster [not for large dense!].
##' Consequently mainly used in  checkMatrix()
## The first version (up to Aug.2022)  -- possibly what we should use for dense case (!?)
boolProd0 <- function(x,y) as((abs(x) %*% abs(y)) > 0, "nMatrix")
## New since  Aug.13, 2022, ensuring that zeros are dropped
isCRT <- function(x, cl = getClass(class(x)))
    extends(cl, "CsparseMatrix") || extends(cl, "TsparseMatrix") || extends(cl, "RsparseMatrix")
boolProd <- function(x,y) {
    ## treat x & y, drop0() & coercing to "n" -- this treats  NA  <==>  1  (!)
    x <- if(isCRT(x)) .sparse2kind(x, kind="n", drop0=TRUE) else as(drop0(x), "nMatrix")
    y <- if(isCRT(y)) .sparse2kind(y, kind="n", drop0=TRUE) else as(drop0(y), "nMatrix")
    r <- (abs(x) %*% abs(y)) > 0
    if(isCRT(r))
        .sparse2kind(r, kind="n", drop0=TRUE)
    else # also for "sparseMatrix" cases "indMatrix" (incl "pMatrix") or "diagonalMatrix"
        ## NB: "diagonalMatrix already *does* drop0(.) when coerced to "nMatrix"
        as(r, "nMatrix")
}
.sparse2kind <- Matrix:::.sparse2kind # (FIXME -- a version of this should be exported!)

###----- Checking a "Matrix" -----------------------------------------

##' Check the compatibility of \pkg{Matrix} package Matrix with a
##' \dQuote{traditional} \R matrix and perform a host of internal consistency
##' checks.
##'
##' @title Check Compatibility of Matrix Package Matrix with Traditional R Matrices
##'
##' @param m   a "Matrix"
##' @param m.m as(m, "matrix")  {if 'do.matrix' }
##' @param do.matrix logical indicating if as(m, "matrix") should be applied;
##'    typically false for large sparse matrices
##' @param do.t  logical: is t(m) "feasible" ?
##' @param doNorm
##' @param doOps
##' @param doSummary
##' @param doCoerce
##' @param doCoerce2
##' @param do.prod
##' @param verbose logical indicating if "progress output" is produced.
##' @param catFUN (when 'verbose' is TRUE): function to be used as generalized cat()
##' @return TRUE (invisibly), unless an error is signalled
##' @author Martin Maechler, since 11 Apr 2008
checkMatrix <- function(m, m.m = if(do.matrix) as(m, "matrix"),
			do.matrix = !isSparse || prod(dim(m)) < 1e6,
			do.t = TRUE, doNorm = TRUE, doOps = TRUE,
                        doSummary = TRUE, doCoerce = TRUE,
			doCoerce2 = doCoerce && !isRsp, doDet = do.matrix,
			do.prod = do.t && do.matrix && !isRsp,
			verbose = TRUE, catFUN = cat,
                        MSG = if(interactive() || capabilities("long.double") ||
                                 isTRUE(get0("doExtras"))) message else function(...) {}
                        )
{
    ## is also called from  dotestMat()  in ../tests/Class+Meth.R

    stopifnot(is(m, "Matrix"))
    validObject(m) # or error(....)
    .isDense  <- Matrix ::: .isDense
    .isSparse <- Matrix ::: .isSparse

    clNam <- class(m)
    cld <- getClassDef(clNam) ## extends(cld, FOO) is faster than is(m, FOO)
    isGen <- extends(cld, "generalMatrix")
    isSym <- extends(cld, "symmetricMatrix")
    isTri <- extends(cld, "triangularMatrix")
    isCor <- isSym && (extends(cld, "corMatrix") || extends(cld, "copMatrix"))
    if(isSparse <- extends(cld, "sparseMatrix")) { # also true for these
        isCsp  <- extends(cld, "CsparseMatrix")
	isRsp  <- extends(cld, "RsparseMatrix")
        isTsp  <- extends(cld, "TsparseMatrix")
	isDiag <- extends(cld, "diagonalMatrix")
	isInd  <- extends(cld, "indMatrix")
	isPerm <- extends(cld, "pMatrix")
    } else isCsp <- isRsp <- isTsp <- isDiag <- isInd <- isPerm <- FALSE
    is.d     <- extends(cld, "dMatrix")
    is.l     <- extends(cld, "lMatrix")
    is.n     <- extends(cld, "nMatrix")
    is.z     <- extends(cld, "zMatrix")
    is.i     <- extends(cld, "iMatrix")
    nonMatr  <- clNam != (Mcl <- MatrixClass(clNam, cld))

    Cat	 <- function(...) if(verbose) cat(...)
    CatF <- function(...) if(verbose) catFUN(...)
    ## warnNow <- function(...) warning(..., call. = FALSE, immediate. = TRUE)

    DO.m <- function(expr) if(do.matrix) eval(expr) else TRUE

    vec <- function(x) {
	dim(x) <- c(length(x), 1L)
	dimnames(x) <- list(NULL,NULL)
	x
    }
    eps16 <- 16 * .Machine$double.eps

    ina <- is.na(m)
    if(do.matrix) {
	stopifnot(all(ina == is.na(m.m)),
                  all(is.nan(m) == is.nan(m.m)),
		  all(is.finite(m) == is.finite(m.m)),
		  all(is.infinite(m) == is.infinite(m.m)),
		  all(m == m | ina), ## check all() , "==" [Compare], "|" [Logic]
		  if(ncol(m) > 0) identical3(unname(m[,1]), unname(m.m[,1]),
					     as(m[,1,drop=FALSE], "vector"))
		  else identical(as(m, "vector"), as.vector(m.m)))
	if(any(m != m & !ina)) stop(" any (m != m) should not be true")
    } else {
	if(any(m != m)) stop(" any (m != m) should not be true")
        if(ncol(m) > 0)
             stopifnot(identical(unname(m[,1]), as(m[,1,drop=FALSE], "vector")))
        else stopifnot(identical(as(m, "vector"), as.vector(as(m, "matrix"))))
    }
    if(do.t) {
	tm <- t(m)
	if(isSym) ## check that t() swaps 'uplo'  L <--> U :
	    stopifnot(c("L","U") == sort(c(m@uplo, tm@uplo)))
	ttm <- t(tm)
        ## notInd: "pMatrix" ok, but others inheriting from "indMatrix" are not
        notInd <- (!isInd || isPerm)
	if(notInd && (isCsp || isGen || isDiag))
            stopifnot(Qidentical(m, ttm, strictClass = !nonMatr))
	else if(do.matrix) {
	    if(notInd) stopifnot(nonMatr || class(ttm) == clNam)
	    stopifnot(all(m == ttm | ina))
	    ## else : not testing
	}


	## crossprod()	%*%  etc
	if(do.prod) {
	    c.m <-  crossprod(m, boolArith = FALSE)
	    tcm <- tcrossprod(m, boolArith = FALSE)
            tolQ <- if(isSparse) NA else eps16
	    stopifnot(dim(c.m) == rep.int(ncol(m), 2),
		      dim(tcm) == rep.int(nrow(m), 2),
		      ## FIXME: %*% drops dimnames
		      Q.eq2(c.m, tm %*% m, tol = tolQ),
		      Q.eq2(tcm, m %*% tm, tol = tolQ),
                      ## should work with dimnames:
		      Q.eq(m %&% tm, boolProd(m, tm), superclasses=NULL, tol = 0)
                     ,
		      Q.eq(tm %&% m, boolProd(tm, m), superclasses=NULL, tol = 0)
                      )
	}
    }
    if(!do.matrix) {
	CatF(" will *not* coerce to 'matrix' since do.matrix is FALSE\n")
    } else if(doNorm) {
	CatF(sprintf(" norm(m [%d x %d]) :", nrow(m), ncol(m)))
	for(typ in c("1","I","F","M")) {
	    Cat('', typ, '')
	    stopifnot(all.equal(norm(m,typ), norm(m.m,typ)))
	}
	Cat(" ok\n")
    }
    if(do.matrix && doSummary) {
	summList <- lapply(getGroupMembers("Summary"), get,
			   envir = asNamespace("Matrix"))
	if(is.z) summList <- summList[match(getGroupMembers("Summary"), c("min", "max", "range"), 0L) == 0L]
	CatF(" Summary: ")
	for(f in summList) {
	    ## suppressWarnings():  e.g. any(<double>)	would warn here:
	    r <- suppressWarnings(identical(f(m), f(m.m)))
	    if(!isTRUE(r)) { ## typically for prod()
		f.nam <- sub("..$", '', sub("^\\.Primitive..", '', format(f)))
		## sum() and prod() are sensitive to order of f. p. operations
		## particularly on systems where sizeof(long double) == sizeof(double)
		(if(any(f.nam == c("sum", "prod"))) MSG else stop)(
		    sprintf("%s(m) [= %g] differs from %s(m.m) [= %g]",
			    f.nam, f(m), f.nam, f(m.m)))
	    }
	}
	if(verbose) cat(" ok\n")
    }

    ## and test 'dim()' as well:
    d <- dim(m)
    isSqr <- d[1] == d[2]
    if(do.t) stopifnot(identical(diag(m), diag(t(m))))
    ## rather FIXME (?)
    ## if(do.t) stopifnot(identical3(diag(m), diag(t(m)), diag(band(m, 0L,0L))))
    if(prod(d) < .Machine$integer.max && !extends(cld, "modelMatrix")) {
	vm <- vec(m)
	stopifnot(is(vm, "Matrix"), validObject(vm), dim(vm) == c(d[1]*d[2], 1))
    }

    if(!isInd)
        m.d <- local({ m. <- m
            diag(m.) <- diag(m) ## << *assigning* to 'm.' now typically annihilates @factor
            if(.hasSlot(m, "factors") && length(f <- m@factors))
                m.@factors <- f
            m. })
    if(do.matrix)
    stopifnot(identical(dim(m.m), dim(m)),
## now that "pMatrix" subsetting gives *LOGICAL*
## 	      if(isPerm) {
## 		  identical(as.integer(unname(diag(m))), unname(diag(m.m)))
## 	      } else
	      identical(diag(m), # base:: *and* Matrix diag()  now keep names
			diag(m.m)),## not for NA: diag(m) == diag(m.m),
	      identical(nnzero(m), sum(m.m != 0)),
	      identical(nnzero(m, na.counted = FALSE),
                        sum(m.m != 0, na.rm = TRUE)),
	      identical(nnzero(m, na.counted = TRUE),
                        sum(m.m != 0 | is.na(m.m)))
	      )

    if(isSparse) {
	n0m <- drop0(m) #==> n0m is Csparse
	has0 <- !Qidentical(n0m, as(m,"CsparseMatrix"))
    }

    if(isDiag)
        stopifnot(exprs = {
            .MJ.Qidentical(m, m.d, strictClass = FALSE,
                           skipSlots = if(m@diag != "N") c("diag", "x"))
            m@diag == "N" || (m.d@diag == "N" &&
                              identical(m.d@x, diag(m, names = FALSE)))
        })
    else if(isTri && m@diag != "N")
        stopifnot(exprs = {
            is(m.d, "triangularMatrix") && m.d@diag == "N"
            .MJ.Qidentical(m, m.d, strictClass = FALSE,
                           skipSlots = c("diag", "p", "i", "j", "x"))
            isSparse || all(m == m.d)
        })
    else if(!isInd && !(is.z && isSym && m@trans == "C"))
        stopifnot(.MJ.Qidentical(m, m.d, strictClass = FALSE,
                                 skipSlots =
                                     if(((isCsp || isRsp) && has0) || isTsp)
                                         c("p", "i", "j", "x")))

    ## use non-square matrix when "allowed":

    ## m12: sparse and may have 0s even if this is not: if(isSparse && has0)
    m12 <- as(as(  m, "lMatrix"),"CsparseMatrix")
    m12 <- drop0(m12)
    if(do.matrix) {
	## "!" should work (via as(*, "l...")) :
	m11 <- as(as(!!m,"CsparseMatrix"), "lMatrix")
	if(!Qidentical(m11, m12))
	    stopifnot(Qidentical(as(m11, "generalMatrix"),
				 as(m12, "generalMatrix")))
    }
    if(isSparse && !isDiag && !is.n) {
	## ensure that as(., "nMatrix") gives nz-pattern
	CatF("as(., \"nMatrix\") giving full nonzero-pattern: ")
	n1 <- as(m, "nMatrix")
	ns <- as(m, "nsparseMatrix")
	stopifnot(identical(n1,ns),
                  ## only testing [CR]sparseMatrix and indMatrix here ...
                  ## sum(<n.T>) excludes duplicated (i,j) pairs whereas
                  ## length(diagU2N(<[^n].T>)) includes them ...
                  isTsp ||
                  (if(isSym) length(if(.hasSlot(n1, "i")) n1@i else n1@j)
                   else sum(n1)) == length(if(isInd) m@perm else diagU2N(m)@x))
        Cat("ok\n")
    }

    if(doOps) {
	## makes sense with non-trivial m (!)
	CatF("2*m =?= m+m: ")
	if(identical(2L*m, m+m)) Cat("identical\n")
	else if(do.matrix) {
	    eq <- as(2L*m,"matrix") == as(m+m, "matrix") # but work for NA's:
	    stopifnot(all(eq | (is.na(m) & is.na(eq))))
	    Cat("ok\n")
	} else {# !do.matrix
	    stopifnot(identical(as(2L*m, "CsparseMatrix"),
                                as(m+m, "CsparseMatrix")))
	    Cat("ok\n")
	}
	if(do.matrix && !is.z) { # complex numbers are *not* ordered
	    ## m == m etc, now for all, see above
	    CatF("m >= m for all: "); stopifnot(all(m >= m | ina)); Cat("ok\n")
	}
	if(prod(d) > 0 && !is.z) {
	    CatF("m < m for none: ")
	    mlm <- m < m
	    if(!any(ina)) stopifnot(!any(mlm))
	    else if(do.matrix) stopifnot(!any(mlm & !ina))
	    else { ## !do.matrix & any(ina) :  !ina can *not* be used
		mlm[ina] <- FALSE
		stopifnot(!any(mlm))
	    }
	    Cat("ok\n")
	}

	if(isSqr) {
	    if(do.matrix) {
		## determinant(<dense>) "fails" for triangular with NA such as
		## (m <- matrix(c(1:0,NA,1), 2))
		CatF("symmpart(m) + skewpart(m) == m: ")
		Q.eq.symmpart(m)
		CatF("ok;  determinant(): ")
		if(!doDet || is.z || d[1] >= 100)
		    Cat(" skipped  ")
		else if(anyNA(m.m) && isTri)
		    Cat(" skipped: is triang. and has NA: ")
		else {
                    d1 <- determinant(m)
                    d2 <- determinant(m.m)
                    mD <- if(d[1]) max(abs(c(d1$modulus, d2$modulus)), na.rm = TRUE) / d[1] else 0
                    largeD <- mD > 4  # unfinished ... seen relatively large diff. for sparse / dense
                    stopifnot(eqDeterminant(m, m.m, d1=d1, d2=d2, NA.Inf.ok=TRUE,
                                            tol = if(largeD) 7/8 else 1e-4))
                    Cat("ok\n")
                }
	    }
	} else assertError(determinant(m))
    }# end{doOps}

    if(doCoerce && do.matrix && canCoerce("matrix", clNam)) {
	CatF("as(<matrix>, ",clNam,"): ", sep='')
	m3 <- as(m.m, clNam)
	Cat("valid:", validObject(m3), "\n")
	## m3 should ``ideally'' be identical to 'm'
    }

    if(doCoerce2 && do.matrix) { ## not for large m:  !m will be dense
        ## FIXME --- some stuff is for all 'kinds'  (e.g. the one for "n" ?)
	if(is.n) { # "nMatrix"
	    mM <- if(nonMatr) as(m, Mcl) else m
	    ## As currently all indexMatrix are nonzero pattern:
	    to. <- if(Mcl == "indMatrix") "indMatrix" else "nMatrix"
	    stopifnot(identical(mM, as(as(m, "dMatrix"),to.)),
		      identical(mM, as(as(m, "lMatrix"),to.)),
		      identical(which(m), which(m.m)))
	}
	else if(is.l) { ## "lMatrix"  -- should fulfill even with NA:
	    stopifnot(all(m | !m | ina), !any(!m & m & !ina))
	    if(isTsp) # allow modify, since at end here
		m <- asUniqueT(m, isT = TRUE)
	    stopifnot(identical(m, m & TRUE),
		      identical(m, FALSE | m))
	    ## also check the  coercions to [dln]Matrix
	    m. <- if(isSparse && has0) n0m else m
	    m1. <- m. # replace NA by 1 in m1. , carefully not changing class:
	    if(any(ina)) m1.@x[is.na(m1.@x)] <- TRUE
	    stopifnot(identical(m. , as(as(m. , "dMatrix"),"lMatrix")),
		      ## coercion to n* and back: only identical when no extra 0s:
		      identical(m1., as(as(m1., "nMatrix"),"lMatrix")),
		      identical(which(m), which(m.m)))
	}
	else if(is.d) { # "dMatrix"
	    m. <- if(isSparse && has0) n0m else m
	    m1 <- m1. <- (m. != 0)*1 # potentially dsparse*
            ## replace NA by 1 in m1. , carefully not changing class:
	    if(any(ina)) m1.@x[is.na(m1.@x)] <- 1
	    ## coercion to n* (nz-pattern!) and back: only identical when no extra 0s and no NAs:
            isSp <- isSparse || xor(.isDense(m1), .isDense(m.)) # only one is sparse
	    stopifnot(Q.C.identical(m1., as(as(m., "nMatrix"),"dMatrix"),
				    isSp, checkClass = FALSE),
		      Q.C.identical(m1 , as(as(m., "lMatrix"),"dMatrix"),
				    isSp, checkClass = FALSE))
	}
        else if(is.i) {
            ## __FIXME__  "iMatrix"
	}
        else if(is.z) {
            ## __FIXME__  "zMatrix"
        }
        maybeDense <- if(isSparse) identity else function(.) as(., "denseMatrix")
	if(isTri) {
	    mm. <- m
	    i0 <- if(m@uplo == "L")
		upper.tri(mm.) else lower.tri(mm.)
	    n.catchWarn <- if(is.n) suppressWarnings else identity
	    n.catchWarn( mm.[i0] <- 0 ) # ideally, mm. remained triangular, but can be dge*
            ## Aug.2022 - Coercion deprecations: No longer do as(*, clNam):
	    CatF("as(mm., \"triangularMatrix\"): ")
	    tm <- as(mm., "triangularMatrix")
	    Cat("valid:", validObject(tm), "\n")
	    if(m@uplo == tm@uplo) { ## otherwise, the matrix effectively was *diagonal*
                if(!isSparse && Matrix:::.isPacked(m)) m <- unpack(m) # to match tm
		## note that diagU2N(<dtr>) |-> dtC, now dtT:
		stopifnot(identical(tm, m) || Qidentical(tm, maybeDense(diagU2N(m))))
            }
	}
	else if(isDiag) {

	    ## TODO

	} else {

	    ## TODO
	}
    }# end {doCoerce2 && ..}

    if(doCoerce && isSparse) { ## coerce to sparseVector and back :
	v <- as(m, "sparseVector")
	stopifnot(length(v) == prod(d))
	dim(v) <- d
	stopifnot(Q.eq2(m, v))
    }

    invisible(TRUE)
} ## {checkMatrix}

### --- These use

##' Check QR-consistency of dense and sparse
chk.qr.D.S <- function(d., s., y, Y = Matrix(y), force = FALSE, tol = 1e-10) {
    stopifnot(is.qr(d.), is(s., "sparseQR"))
    cc <- qr.coef(d.,y)
    rank.def <- any(is.na(cc)) && d.$rank < length(d.$pivot)
    if(rank.def && force) cc <- mkNA.0(cc) ## set NA's to 0 .. ok, in some case

    ## when system is rank deficient, have differing cases, not always just NA <-> 0 coef
    ## FIXME though:  resid & fitted should be well determined
    if(force || !rank.def) stopifnot(
	is.all.equal3(	    cc	     , qr.coef  (s.,y), drop(qr.coef  (s.,Y)), tol=tol),
	is.all.equal3(qr.resid (d.,y), qr.resid (s.,y), drop(qr.resid (s.,Y)), tol=tol),
	is.all.equal3(qr.fitted(d.,y), qr.fitted(s.,y), drop(qr.fitted(s.,Y)), tol=tol)
	)
}

##' "Combi" calling chkQR() on both "(sparse)Matrix" and 'traditional' version
##' ------  and combine the two qr decompositions using chk.qr.D.S()
##'         [ chkQR() def. in >>>>> ./test-tools-1.R <<<<< ]
##'
##' @title check QR-decomposition, and compare sparse and dense one
##' @param A a 'Matrix' , typically 'sparseMatrix'
##' @param Qinv.chk
##' @param QtQ.chk
##' @param quiet
##' @return list with 'qA' (sparse QR) and 'qa' (traditional (dense) QR)
##' @author Martin Maechler
checkQR.DS.both <- function(A, Qinv.chk, QtQ.chk=NA,
                            quiet=FALSE, giveRE=TRUE, tol = 1e-13)
{
    stopifnot(is(A,"Matrix"))
    if(!quiet) cat("classical: ")
    qa <- chkQR(as(A, "matrix"), Qinv.chk=TRUE, QtQ.chk=TRUE, tol=tol, giveRE=giveRE)# works always
    if(!quiet) cat("[Ok] ---  sparse: ")
    qA <- chkQR(A, Qinv.chk=Qinv.chk, QtQ.chk=QtQ.chk, tol=tol, giveRE=giveRE)
    validObject(qA)
    if(!quiet) cat("[Ok]\n")
    chk.qr.D.S(qa, qA, y = 10 + 1:nrow(A), tol = 256*tol)# ok [not done in rank deficient case!]
    invisible(list(qA=qA, qa=qa))
}

non0.ij <- function(M) Matrix:::non0.i(as(M, "sparseMatrix"))

triuChk <- function(x, k) {
    ans <- triu(x, k)
    ij <- non0.ij(ans)
    stopifnot(identical(dim(x), dim(ans)), (ij %*% c(-1,1)) >= k)
    ans
}

trilChk <- function(x, k) {
    ans <- tril(x, k)
    ij <- non0.ij(ans)
    stopifnot(identical(dim(x), dim(ans)), (ij %*% c(-1,1)) <= k)
    ans
}
