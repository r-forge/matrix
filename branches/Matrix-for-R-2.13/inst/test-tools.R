#### Will be sourced by several R scripts in ../tests/

### ------- Part I --  unrelated to "Matrix" classes ---------------

paste0 <- function(...) paste(..., sep = '')

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)
identical5 <- function(a,b,c,d,e) identical(a,b) && identical4(b,c,d,e)

## Make sure errors are signaled
assertError <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, error = function(e) e)
    if(!inherits(t.res, "error"))
	stop(d.expr, "\n\t did not give an error", call. = FALSE)
    invisible(t.res)
}

isValid <- function(x, class) validObject(x, test=TRUE) && is(x, class)

is.all.equal3 <- function(x,y,z, tol = .Machine$double.eps^0.5)
    isTRUE(all.equal(x,y, tol=tol)) && isTRUE(all.equal(y,z, tol=tol))

is.all.equal4 <- function(x,y,z,u, tol = .Machine$double.eps^0.5)
    is.all.equal3(x,y,z, tol=tol) && isTRUE(all.equal(z,u, tol=tol))

## A version of all.equal() for the slots
all.slot.equal <- function(x,y, ...) {
    slts <- slotNames(x)
    for(sl in slts) {
        aeq <- all.equal(slot(x,sl), slot(y,sl), ...)
        if(!identical(TRUE, aeq))
            return(paste("slot '",sl,"': ", aeq, sep=''))
    }
    TRUE
}

## all.equal() for list-coercible objects -- apart from *some* components
all.equal.X <- function(x,y, except, ...)
{
    .trunc <- function(x) {
	ll <- as.list(x)
	ll[ - match(except, names(ll), nomatch = 0L)]
    }
    all.equal(.trunc(x), .trunc(y), ...)
}
## e.g. in lme4:
##  all.equal.X(env(m1), env(m2), except = c("call", "frame"))

## The relative error typically returned by all.equal:
relErr <- function(target, current) { ## make this work for 'Matrix'
    ## ==> no mean() ..
    n <- length(current)
    if(length(target) < n)
        target <- rep(target, length.out = n)
    sum(abs(target - current)) / sum(abs(target))
}

## is.R22 <- (paste(R.version$major, R.version$minor, sep=".") >= "2.2")

pkgRversion <- function(pkgname)
    substring(packageDescription(pkgname)[["Built"]], 3,5)

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

##' @title turn an S4 object (with slots) into a list with corresponding components
##' @param obj an R object with a formal class (aka "S4")
##' @return a list with named components where \code{obj} had slots
##' @author Martin Maechler
S4_2list <- function(obj) {
   sn <- slotNames(obj)
   structure(lapply(sn, slot, object = obj), .Names = sn)
}


### ------- Part II  -- related to matrices, but *not* "Matrix" -----------

add.simpleDimnames <- function(m) {
    stopifnot(length(d <- dim(m)) == 2)
    dimnames(m) <- list(if(d[1]) paste0("r", seq_len(d[1])),
                        if(d[2]) paste0("c", seq_len(d[2])))
    m
}

as.mat <- function(m) {
    ## as(., "matrix")	but with no extraneous empty dimnames
    m <- as(m, "matrix")
    if(identical(dimnames(m), list(NULL,NULL)))
	dimnames(m) <- NULL
    m
}

assert.EQ.mat <- function(M, m, tol = if(show) 0 else 1e-15, show=FALSE) {
    ## Purpose: check equality of  'Matrix' M with  'matrix' m
    ## ----------------------------------------------------------------------
    ## Arguments: M: is(., "Matrix") typically {but just needs working as(., "matrix")}
    ##            m: is(., "matrix")
    ##            show: if TRUE, return (and hence typically print) all.equal(...)
    validObject(M)
    MM <- as.mat(M)                     # as(M, "matrix")
    if(is.logical(MM) && is.numeric(m))
	storage.mode(MM) <- "integer"
    attr(MM, "dimnames") <- attr(m, "dimnames") <- NULL
    if(show) all.equal(MM, m, tol = tol)
    else if(!isTRUE(r <- all.equal(MM, m, tol = tol)))
	stop("all.equal() |->  ", r)
}

chk.matrix <- function(M) {
    ## check object; including coercion to "matrix" :
    cl <- class(M)
    cat("class ", dQuote(cl), " [",nrow(M)," x ",ncol(M),"]; slots (",
	paste(slotNames(M), collapse=","), ")\n", sep='')
    stopifnot(validObject(M),
	      dim(M) == c(nrow(M), ncol(M)),
	      identical(dim(m <- as(M, "matrix")), dim(M))
	      )
}

isOrthogonal <- function(x, tol = 1e-15) {
    all.equal(diag(as(zapsmall(crossprod(x)), "diagonalMatrix")),
              rep(1, ncol(x)), tol = tol)
}


### ------- Part III --  "Matrix" (classes) specific ----------------------

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
asCsp <- function(x) Matrix:::diagU2N(as(x, "CsparseMatrix"))

isNull.dimnames <- function(dn) {
    is.null(dn) || {
	dn <- unname(dn)
	ch0 <- character(0)
	identical(dn, list(NULL,NULL)) ||
	identical(dn, list(ch0, NULL)) ||
	identical(dn, list(NULL, ch0)) ||
	identical(dn, list(ch0, ch0))
    }
}

Qidentical.DN <- function(dx, dy) {
    ## quasi-identical dimnames
    stopifnot(is.list(dx) || is.null(dx),
	      is.list(dy) || is.null(dy))
    ## "empty"
    (isNull.dimnames(dx) && isNull.dimnames(dy)) || identical(dx, dy)
}

Qidentical <- function(x,y, strictClass = TRUE) {
    ## quasi-identical - for 'Matrix' matrices
    if(class(x) != class(y)) {
        if(strictClass || !is(x, class(y)))
           return(FALSE)
        ## else try further
    }
    slts <- slotNames(x)
    if("Dimnames" %in% slts) { ## always (or we have no 'Matrix')
	slts <- slts[slts != "Dimnames"]
	if(!(Qidentical.DN(x@Dimnames, y@Dimnames)))
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
##' @param superclasses   x and y must coincide in (not) extending these
##' @param dimnames.check  logical indicating if dimnames(.) much match
##' @param tol  NA (--> use "==") or numerical tolerance for all.equal()
##' @return   logical:  Are x and y (quasi) equal ?
Q.eq <- function(x, y,
		 superclasses =
		 c("sparseMatrix", "denseMatrix",
		   "dMatrix", "lMatrix", "nMatrix"),
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
    asC <- ## asCommon
        if((isDense <- extends(xcl,"denseMatrix")))
            function(m) as(m, "matrix")
        else function(m)
            as(as(as(m,"CsparseMatrix"), "dMatrix"), "dgCMatrix")
    if(is.na(tol)) {
	if(isDense)
	    all(x == y | (is.na(x) & is.na(y)))
	else ## 'x == y' blows up for large sparse matrices:
	    isTRUE(all.equal(asC(x), asC(y), tol = 0.,
			     check.attributes = dimnames.check))
    }
    else if(is.numeric(tol) && tol >= 0) {
        isTRUE(all.equal(asC(x), asC(y), tol = tol,
                         check.attributes = dimnames.check))
    }
    else stop("'tol' must be NA or non-negative number")
}

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


## Useful Matrix constructors for testing:

rspMat <- function(n, m = n, density = 1/4, nnz = round(density * n*m))
{
    ## Purpose: random sparse Matrix
    ## ----------------------------------------------------------------------
    ## Arguments: (n,m) : dimension [default m=n ==> *square* matrix}
    ##           density: the desired sparseness density:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 5 Mar 2008, 11:07
    stopifnot(length(n) == 1, n == as.integer(n),
              length(m) == 1, m == as.integer(m),
              0 <= density, density <= 1,
              0 <= nnz, nnz <= n*m)
    x <- numeric(n*m)
    ## entries 2 : (nnz+1) {so we can have '1' as 'special'}
    x[sample(n*m, nnz, replace=FALSE)] <- as.numeric(1L + seq_len(nnz))
    Matrix(x, n,m, sparse=TRUE)
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

## This is related to rUnitTri(), ver
mkLDL <- function(n, density = 1/3) {
    ## Purpose: make nice artificial   A = L D L'  (with exact numbers) decomp
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 15 Mar 2008
    stopifnot(n == round(n))
    n <- as.integer(n)
    L <- Matrix(0, n,n)

    nnz <- round(n*n * density)
    L[sample(n*n, nnz)] <- seq_len(nnz)
    L <- tril(L,-1)
    diag(L) <- 1
    d.half <- sample(10*(n:1))# random permutation ; use '10*' to be "different" from L entries
    D <- Diagonal(x = d.half * d.half)
    A <- tcrossprod(L * rep(d.half, each=n))
    ## = as(L %*% D %*% t(L), "symmetricMatrix")
    list(A = A, L = L, d.half = d.half, D = D)
}

eqDeterminant <- function(m1, m2, NA.Inf.ok=FALSE, ...) {
    d1 <- determinant(m1) ## logarithm = TRUE
    d2 <- determinant(m2)
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
    all.equal(d1, d2, ...)
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
        is.perm <- function(.)
            if(inherits(., "try-error")) NA else !all(.@perm == 0:(.@Dim[1]-1))
        is.LDL <- function(.)if(inherits(., "try-error")) NA else isLDL(.)
	r.st <-
	    cbind(perm	= sapply(r, is.perm),
		  LDL	= sapply(r, is.LDL),
		  super = sapply(r, class) == "dCHMsuper")
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
    list(dup.r.all = dup.r1,
	 r.all	= r.all,
	 r.uniq = CHM_to_pLs(r1[ ! dup.r1]))
}


###----- Checking a "Matrix" -----------------------------------------

##' <description>
##'
##' <details>
##' @title Compatibility tests "Matrix" <-> "traditional matrix"
##'        and many more consistency checks
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
			doCoerce2 = doCoerce && !extends(cld, "RsparseMatrix"),
			do.prod = do.t && do.matrix && !extends(cld, "RsparseMatrix"),
			verbose = TRUE, catFUN = cat)
{
    ## is also called from  dotestMat()  in ../tests/Class+Meth.R

    stopifnot(is(m, "Matrix"))
    validObject(m) # or error(....)

    clNam <- class(m)
    cld <- getClassDef(clNam) ## extends(cld, FOO) is faster than is(m, FOO)
    isCor    <- extends(cld, "corMatrix")
    isSparse <- extends(cld, "sparseMatrix") # also true for diagonalMatrix !
    isSym    <- extends(cld, "symmetricMatrix")
    isDiag   <- extends(cld, "diagonalMatrix")
    nonMatr  <- clNam != Matrix:::MatrixClass(clNam, cld)

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
	if(extends(cld, "CsparseMatrix") ||
	   extends(cld, "generalMatrix") || isDiag)
            stopifnot(Qidentical(m, ttm, strictClass = !nonMatr))
	else if(do.matrix)
	    stopifnot(nonMatr || class(ttm) == clNam,
                      all(m == ttm | ina))
        ## else : not testing

	## crossprod()	%*%  etc
	if(do.prod) {
	    c.m <-  crossprod(m)
	    tcm <- tcrossprod(m)
            tolQ <- if(isSparse) NA else eps16
	    stopifnot(dim(c.m) == rep.int(ncol(m), 2),
		      dim(tcm) == rep.int(nrow(m), 2),
		      ## FIXME: %*% drops dimnames
		      Q.eq2(c.m, tm %*% m, tol = tolQ),
		      Q.eq2(tcm, m %*% tm, tol = tolQ))
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
	CatF(" Summary: ")
	for(f in summList) {
	    ## suppressWarnings():  e.g. any(<double>)	would warn here:
	    r <- suppressWarnings(if(isCor) all.equal(f(m), f(m.m)) else
				  identical(f(m), f(m.m)))
	    if(!isTRUE(r)) {
		f.nam <- sub("..$", '', sub("^\\.Primitive..", '', format(f)))
		## prod() is delicate: NA or NaN can both happen
		(if(f.nam == "prod") message else stop)(
		    sprintf("%s(m) [= %g] differs from %s(m.m) [= %g]",
			    f.nam, f(m), f.nam, f(m.m)))
	    }
	}
	if(verbose) cat(" ok\n")
    }

    ## and test 'dim()' as well:
    d <- dim(m)
    isSqr <- d[1] == d[2]
    isPerm <- extends(cld, "pMatrix")
    if(do.t) stopifnot(identical(diag(m), diag(t(m))))
    ## TODO: also === diag(band(m,0,0))

    if(prod(d) < .Machine$integer.max && !extends(cld, "modelMatrix")) {
	vm <- vec(m)
	stopifnot(is(vm, "Matrix"), validObject(vm), dim(vm) == c(d[1]*d[2], 1))
    }

    if(do.matrix)
    stopifnot(identical(dim(m.m), dim(m)),
	      ## base::diag() keeps names [Matrix FIXME]
## now that "pMatrix" subsetting gives *LOGICAL*
## 	      if(isPerm) {
## 		  identical(as.integer(unname(diag(m))), unname(diag(m.m)))
## 	      } else
	      identical(unname(diag(m)),
			unname(diag(m.m))),## not for NA: diag(m) == diag(m.m),
	      identical(nnzero(m), sum(m.m != 0)),
	      identical(nnzero(m, na.= FALSE), sum(m.m != 0, na.rm = TRUE)),
	      identical(nnzero(m, na.= TRUE),  sum(m.m != 0 | is.na(m.m)))
	      )

    if(isSparse) {
	n0m <- drop0(m) #==> n0m is Csparse
	has0 <- !Qidentical(n0m, as(m,"CsparseMatrix"))
    }
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
    if(isSparse && !extends(cld,"nMatrix")) {
	## ensure that as(., "nMatrix") gives nz-pattern
	CatF("as(., \"nMatrix\") giving full nonzero-pattern: ")
	n1 <- as(m, "nMatrix")
	ns <- as(m, "nsparseMatrix")
	stopifnot(identical(n1,ns),
		  isDiag || ((if(isSym) Matrix:::nnzSparse else sum)(n1) ==
			     length(if(isPerm) m@perm else Matrix:::diagU2N(m)@x)))
        Cat("ok\n")
    }

    if(doOps) {
	## makes sense with non-trivial m (!)
	CatF("2*m =?= m+m: ")
	if(identical(2*m, m+m)) Cat("identical\n")
	else if(do.matrix) {
	    eq <- as(2*m,"matrix") == as(m+m, "matrix") # but work for NA's:
	    stopifnot(all(eq | (is.na(m) & is.na(eq))))
	    Cat("ok\n")
	} else {# !do.matrix
	    stopifnot(identical(as(2*m, "CsparseMatrix"),
                                as(m+m, "CsparseMatrix")))
	    Cat("ok\n")
	}
	if(do.matrix) {
	    ## m == m etc, now for all, see above
	    CatF("m >= m for all: "); stopifnot(all(m >= m | ina)); Cat("ok\n")
	}
	if(prod(d) > 0) {
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
		if(any(is.na(m.m)) && extends(cld, "triangularMatrix"))
		    Cat(" skipped: is triang. and has NA")
		else
		    stopifnot(eqDeterminant(m, m.m, NA.Inf.ok=TRUE))
		Cat("ok\n")
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

	if(extends(cld, "nMatrix")) {
	    stopifnot(identical(m, as(as(m, "dMatrix"),"nMatrix")),
		      identical(m, as(as(m, "lMatrix"),"nMatrix")),
		      identical(which(m), which(m.m)))
	}
	else if(extends(cld, "lMatrix")) { ## should fulfill even with NA:
	    stopifnot(all(m | !m | ina), !any(!m & m & !ina))
	    if(extends(cld, "TsparseMatrix")) # allow modify, since at end here
		m <- Matrix:::uniqTsparse(m, clNam)
	    stopifnot(identical(m, m & TRUE),
		      identical(m, FALSE | m))
	    ## also check the  coercions to [dln]Matrix
	    m. <- if(isSparse && has0) n0m else m
	    m1. <- m. # replace NA by 1 in m1. , carefully not changing class:
	    if(any(ina)) m1.@x[is.na(m1.@x)] <- TRUE
	    stopifnot(identical(m. , as(as(m. , "dMatrix"),"lMatrix")),
		      clNam == "ldiMatrix" || # <- there's no "ndiMatrix"
		      ## coercion to n* and back: only identical when no extra 0s:
		      identical(m1., as(as(m1., "nMatrix"),"lMatrix")),
		      identical(which(m), which(m.m)))
	}
	else if(extends(cld, "dMatrix")) {
	    m. <- if(isSparse && has0) n0m else m
	    m1 <- (m. != 0)*1
	    if(!isSparse && substr(clNam,1,3) == "dpp")
		## no "nppMatrix" possible
		m1 <- unpack(m1)

	    m1. <- m1 # replace NA by 1 in m1. , carefully not changing class:
	    if(any(ina)) m1.@x[is.na(m1.@x)] <- 1
	    ## coercion to n* (nz-pattern!) and back: only identical when no extra 0s and no NAs:
	    stopifnot(Q.C.identical(m1., as(as(m., "nMatrix"),"dMatrix"),
				    isSparse, checkClass = FALSE),
		      Q.C.identical(m1 , as(as(m., "lMatrix"),"dMatrix"),
				    isSparse, checkClass = FALSE))
	}

	if(extends(cld, "triangularMatrix")) {
	    mm. <- m
	    i0 <- if(m@uplo == "L")
		upper.tri(mm.) else lower.tri(mm.)
	    mm.[i0] <- 0 # ideally, mm. remained triangular, but can be dge*
	    CatF("as(<triangular (ge)matrix>, ",clNam,"): ", sep='')
	    tm <- as(as(mm., "triangularMatrix"), clNam)
	    Cat("valid:", validObject(tm), "\n")
	    if(m@uplo == tm@uplo) ## otherwise, the matrix effectively was *diagonal*
		## note that diagU2N(<dtr>) |-> dtC :
		stopifnot(Qidentical(tm, as(Matrix:::diagU2N(m), clNam)))
	}
	else if(extends(cld, "diagonalMatrix")) {

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
}

