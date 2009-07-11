####  Utilities  for  Sparse Model Matrices

## The "first" version {no longer used}:
fac2sparse <- function(from, to = c("d","i","l","n","z"), drop.unused.levels = TRUE)
{
    ## factor(-like) --> sparseMatrix {also works for integer, character}
    fact <- if (drop.unused.levels) factor(from) else as.factor(from)
    levs <- levels(fact)
    n <- length(fact)
    to <- match.arg(to)
    ## MM: using new() and then assigning slots has efficiency "advantage"
    ##     of *not* validity checking
    res <- new(paste(to, "gCMatrix", sep=''))
    res@i <- as.integer(fact) - 1L # 0-based
    res@p <- 0:n
    res@Dim <- c(length(levs), n)
    res@Dimnames <- list(levs, NULL)
    if(to != "n")
	res@x <- rep.int(switch(to,
				"d" = 1., "i" = 1L, "l" = TRUE, "z" = 1+0i),
			 n)
    res
}

## This version can deal with NA's -- but is less efficient (how much?) :
fac2sparse <- function(from, to = c("d","i","l","n","z"),
                       drop.unused.levels = TRUE)
{
    ## factor(-like) --> sparseMatrix {also works for integer, character}
    fact <- if (drop.unused.levels) factor(from) else as.factor(from)
    levs <- levels(fact)
    n <- length(fact)
    to <- match.arg(to)
    i <- as.integer(fact) - 1L                  # 0-based indices
    df <- data.frame(i = i, j = seq_len(n) - 1L)[!is.na(i),]
    if(to != "n")
	df$x <- rep.int(switch(to,
			       "d" = 1., "i" = 1L, "l" = TRUE, "z" = 1+0i),
			nrow(df))
    as(do.call("new", c(list(Class = paste(to, "gTMatrix", sep=''),
			     Dim = c(length(levs), n),
			     Dimnames = list(levs, names(fact))),
			df)),
       "CsparseMatrix")
}

setAs("factor", "sparseMatrix", function(from) fac2sparse(from, to = "d"))

## a version that uses contrasts --- *iff* contrasts.arg is not FALSE
fac2Sparse <- function(from, to = c("d","i","l","n","z"),
                       drop.unused.levels = TRUE, contrasts.arg = NULL)
{

    m <- fac2sparse(from, to=to,
                    drop.unused.levels=drop.unused.levels)
    if(identical(contrasts.arg, FALSE))
	return(m)
    ## else *do* use contrasts.arg
    if(is.null(contrasts.arg))
	contrasts.arg <- getOption("contrasts")[if(is.ordered(from))
						"ordered" else "unordered"]
    stopifnot(is.function(FUN <- get(contrasts.arg)))
    t(FUN(length(levels(from)), sparse = TRUE)) %*% m
}

if(getRversion() < "2.10.0" || R.version$`svn rev` < 48913) {
### Define contr.sum() etc  with 'sparse' argument :

.sparse.array <- function(x, dim, dimnames)
    Matrix::Matrix(x, nrow=dim[1], ncol=dim[2],
		   dimnames=dimnames, sparse=TRUE)

contr.helmert <-
    function (n, contrasts=TRUE, sparse=FALSE)
{
    if (length(n) <= 1) {
	if(is.numeric(n) && length(n) == 1 && n > 1) levels <- seq_len(n)
	else stop("not enough degrees of freedom to define contrasts")
    } else levels <- n
    llev <- length(levels <- as.character(levels))
    if(sparse) array <- .sparse.array
    if (contrasts) {
	cont <- array(-1, c(llev, llev-1L), list(levels, NULL))
	cont[col(cont) <= row(cont) - 2] <- 0
	cont[col(cont) == row(cont) - 1] <- 1L:(llev-1)
    } else {
	cont <- array(0, c(llev, llev), list(levels, levels))
	cont[col(cont) == row(cont)] <- 1
    }
    cont
}

contr.treatment <-
    function(n, base = 1, contrasts = TRUE, sparse = FALSE)
{
    if(is.numeric(n) && length(n) == 1) {
	if(n > 1) levels <- as.character(seq_len(n))
	else stop("not enough degrees of freedom to define contrasts")
    } else {
	levels <- as.character(n)
	n <- length(n)
    }

    if(sparse) array <- .sparse.array
    contr <- array(0, c(n, n), list(levels, levels))
    diag(contr) <- 1
    if(contrasts) {
	if(n < 2)
	    stop(gettextf("contrasts not defined for %d degrees of freedom",
                          n - 1), domain = NA)
	if (base < 1 | base > n)
	    stop("baseline group number out of range")
	contr <- contr[, -base, drop = FALSE]
    }
    contr
}

contr.sum <-
    function (n, contrasts=TRUE, sparse=FALSE)
{
    if (length(n) <= 1) {
	if (is.numeric(n) && length(n) == 1 && n > 1)
	    levels <- seq_len(n)
	else stop("not enough degrees of freedom to define contrasts")
    } else levels <- n
    llev <- length(levels <- as.character(levels))
    if(sparse) array <- .sparse.array
    if (contrasts) {
	cont <- array(0, c(llev, llev - 1L), list(levels, NULL))
	cont[col(cont) == row(cont)] <- 1
	cont[llev, ] <- -1
    } else {
	cont <- array(0, c(llev, llev), list(levels, levels))
	cont[col(cont) == row(cont)] <- 1
    }
    cont
}

contr.SAS <- function(n, contrasts = TRUE, sparse=FALSE)
{
    contr.treatment(n,
                    base = if (is.numeric(n) && length(n) == 1) n else length(n),
                    contrasts, sparse=sparse)
}

contr.poly <- function (n, scores = 1L:n, contrasts = TRUE, sparse = FALSE) {
    ## this is non-sense anyway
    as(stats::contr.poly(n, scores=scores, contrasts=contrasts),
       "sparseMatrix")
}


}## end if (<old R version>)


## Goal: an  model.sparseMatrix()
##      model.matrix(object, data = environment(object),
##                   contrasts.arg = NULL, xlev = NULL, ...)
## "FIXME": Rather should have
##    model.matrix(.......,  sparse = TRUE) ...
##
## This is Cut'n'Paste from model.matrix() ... just replacing one small part
sparse.model.matrix <- function(object, data = environment(object),
                                contrasts.arg = NULL, xlev = NULL, ...)
{
    t <- if(missing(data)) terms(object) else terms(object, data=data)
    if (is.null(attr(data, "terms")))
	data <- model.frame(object, data, xlev=xlev)
    else {
	reorder <- match(sapply(attr(t,"variables"),deparse,
                                width.cutoff=500)[-1L],
                         names(data))
	if (any(is.na(reorder)))
	    stop("model frame and formula mismatch in model.matrix()")
	if(!identical(reorder, seq_len(ncol(data))))
	    data <- data[,reorder, drop=FALSE]
    }
    int <- attr(t, "response")
    if(length(data)) {      # otherwise no rhs terms, so skip all this
        contr.funs <- as.character(getOption("contrasts"))
        namD <- names(data)
        ## turn any character columns into factors
        for(i in namD)
            if(is.character(data[[i]])) {
                data[[i]] <- factor(data[[i]])
                warning(gettextf("variable '%s' converted to a factor", i),
                        domain = NA)
            }
        isF <- sapply(data, function(x) is.factor(x) || is.logical(x) )
        isF[int] <- FALSE
        isOF <- sapply(data, is.ordered)
        for(nn in namD[isF])            # drop response
            if(is.null(attr(data[[nn]], "contrasts")))
                contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
        ## it might be safer to have numerical contrasts:
        ##	  get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
        if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
            if (is.null(namC <- names(contrasts.arg)))
                stop("invalid 'contrasts.arg' argument")
            for (nn in namC) {
                if (is.na(ni <- match(nn, namD)))
                    warning(gettextf("variable '%s' is absent, its contrast will be ignored", nn),
                            domain = NA)
                else {
                    ca <- contrasts.arg[[nn]]
                    if(is.matrix(ca)) contrasts(data[[ni]], ncol(ca)) <- ca
                    else contrasts(data[[ni]]) <- contrasts.arg[[nn]]
                }
            }
        }
    } else {               # internal model.matrix needs some variable
        isF <-  FALSE
        data <- list(x=rep(0, nrow(data)))
    }
    ## <Sparse> src/library/stats/R/models.R has
    ##    ans <- .Internal(model.matrix(t, data))
    ans <- model.spmatrix(t, data)
    ## </Sparse>

    cons <- if(any(isF))
	lapply(data[isF], function(x) attr(x,  "contrasts"))
    else NULL
    attr(ans, "contrasts") <- cons
    ans
}


sparseInt.r <- function(X, Y, do.names = TRUE)
{
    ## Produce the t(Z); Z = "design matrix" of (X : Y), where
    ##             --- t(Z) : aka rowwise -version : "r"

    ## X, Y either are numeric vectors (well matrix is possible)
    ##      or as(<factor>, sparseM):

    r <-
        if((nX <- is.numeric(X)) | (nY <- is.numeric(Y))) {
            if(nX) {
                if (nY) X * Y
                else {  ## numeric X,  sparseMatrix Y
                    r <- Y
                    r@x <- X * Y@x
                    r
                }
            }
            else { ## sparseMatrix X, numeric Y
                r <- X
                r@x <- Y * X@x
                r
            }

        }
        else { ## X & Y are both sparseMatrix
            do.call("rBind",
                    lapply(seq_len(nrow(Y)),
                           function(i) {
                               Yi <- Y[i, , drop = FALSE]
                               dimnames(Yi) <- list(NULL,NULL)
                               X * Yi[rep(1,nrow(X)), , drop = FALSE]
                           }))
        }

    if(do.names) {
        ## FIXME: This names business needs a good solution..
        ##        but maybe "up in the caller"
        NULL.if.0 <- function(.) if(length(.) > 0) . # else NULL
        if(!is.null(dim(r)) &&
           !is.null(nX <- rownames(X)) &&
           !is.null(nY <- rownames(Y)))
            rownames(r) <- outer(nX, nY, paste, sep = ":")
    }
    r
}


is.model.frame <- function(x)
{
  ## Purpose: check if x is a "valid" model.frame
  ## ------------------------------------------------------------
  ## Author: Martin Maechler, Date: 3 Jul 2009
    is.data.frame(x) &&
    !is.null(tms <- attr(x, "terms")) &&
    inherits(tms, "terms") && ## is.terms() would be better
    inherits(tms, "formula") &&
    is.matrix(attr(tms, "factors")) &&
    is.language(vv <- attr(tms, "variables")) &&
    vv[[1]] == as.symbol("list") &&
    all((vars <- sapply(as.list(vv[-1]), as.character)) %in% colnames(x))
    ## and we could go on testing vars
}

## This version uses  'rBind' and returns  X' { t(X) ] :
model.spmatrix <- function(trms, mf, transpose=FALSE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments: ...
    ##	transpose: if TRUE, return X' {which is faster}, otherwise (default): X
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  7 Jul 2009

    ## mf is a model frame or a "simple" data.frame [after reorder !]
    stopifnot(is.data.frame(mf))

##     cat("model.spmatrix(trms, mf): \n trms = \n")
##     str(trms)
##     cat("-------  mf = \n")
##     str(mf)

    hasInt <- attr(trms, "intercept") == 1
    if(!hasInt)
        stop("formula with*OUT* intercept not yet implemented")

    ## Create a sparse model matrix from a model frame.
    ## TH: I expect that at least one component is a factor or sparse


    termsFactors <- attr(trms, "factors")
    Names <- colnames(termsFactors)
    ## the degree of interaction:
    intOrder <- attr(trms, "order")
    isInteraction <- intOrder > 1L
    simpleNames <- Names[!isInteraction]
    namObj <- names(mf <- unclass(mf))
    attributes(mf) <- list(names = namObj) # i.e. drop all other attributes

    ## Convert character & factor to "Rowwise- sparseMatrix ("dummy"-matrix)
    ich <- sapply(mf, is.character)
    mf[ich] <- lapply(mf[ich], factor)
    is.f <- sapply(mf, is.factor)
    mf[is.f] <- mapply(function(f, nam) {
        r <- ## as(f, "sparseMatrix")
            fac2Sparse(f, to = "d",
                       drop.unused.levels = TRUE,
                       contrasts.arg = attr(f, "contrasts"))
        ## for some contrast {contr.sum}, the above *loses* rownames .. hmm ..
	rownames(r) <-
	    paste(nam, if(is.null(rownames(r))) seq_len(nrow(r)) else rownames(r),
		  sep="")
        ## if(nrow(r) > 1) seq_len(nrow(r)), sep="")
        r
    }, mf[is.f], namObj[is.f])

    ## mf: now a list of *numeric* vectors and sparseMatrix (ex-factors)

    result <- structure(vector("list", length = length(Names)),
                        names = Names)
    result[simpleNames] <- mf[simpleNames]

    ## Now handle interactions.
    ##  Use previous columns in result (especially for high-order interactions).
    if((maxDeg <- max(intOrder)) > 1)
        ## systematically do 2-way, 3-way {using 2-way}, etc etc
        for(deg in 2:maxDeg) {
            for(nm in Names[intOrder == deg]) {
                ## Split name into two pieces at the last ":"
                nmSplits <- strsplit(nm, ":", fixed=TRUE)[[1]]
                if(deg > 2) {
                    lengths <- nchar(nmSplits[-deg])
                    nmSplits <- c(substring(nm, 1L, sum(lengths)+deg-2L), nmSplits[deg])
                }
                if(getOption("verbose"))
                    cat(sprintf("interaction '%s' from   '%s' * '%s' \n",
                                nm, nmSplits[1],nmSplits[2]))
                result[[nm]] <- sparseInt.r(result[[nmSplits[1]]],
                                            result[[nmSplits[2]]])
            }
        }

    myFormula <- trms
    attributes(myFormula) <- NULL       # only need formula
    ## r :=  X^T
    r <- structure(do.call("rBind",
			   c(if(hasInt) list("(Intercept)" = 1), result)),
		   ## extra attributes added to the sparse Matrix:
                   ##  [do *NOT* cobble slots!]
		   termNames = names(result),
		   df = sapply(result, ncol),
		   call = match.call(),
		   formula = myFormula)
    if(!transpose) t(r) else r
}

### Keep this namespace-hidden: Would need to return a classed object

## FIXME: still test this function for both methods, since currently
## ----- both  dgCMatrix_cholsol and  dgCMatrix_qrsol are only called from here!
lm.fit.sparse <- function(x, y, offset = NULL, method = c("qr", "cholesky"),
                          tol = 1e-7, singular.ok = TRUE, order = NULL,
                          transpose = FALSE) ## NB: meaning of 'transpose'
                                        # is changed from original

### Fit a linear model, __ given __ a sparse model matrix 'x'
### using a sparse QR or a sparse Cholesky factorization
{
    cld <- getClass(class(x))
    stopifnot(extends(cld, "dsparseMatrix"))
## or     if(!is(x, "dsparseMatrix")) x <- as(x, "dsparseMatrix")
    yy <- as.numeric(y)
    if (!is.null(offset)) {
	stopifnot(length(offset) == length(y))
	yy <- yy - as.numeric(offset)
    }
    method <- match.arg(method)
    order <- {
        if(is.null(order)) ## recommended default depends on method :
            if(method == "qr") 3L else 1L
        else as.integer(order) }

    if(transpose) x <- t(x)
    ans <- switch(method,
		  cholesky =
		  .Call(dgCMatrix_cholsol,# has AS_CHM_SP(x)
			as(x, "CsparseMatrix"), yy),
		  qr =
		  .Call(dgCMatrix_qrsol, # has AS_CSP(): must be dgC or dtC:
			if(cld@className %in% c("dtCMatrix", "dgCMatrix")) x
			else as(x, "dgCMatrix"),
			yy, order),
		  ## otherwise:
		  stop("unknown method ", dQuote(method))
		  )
    ans
}

