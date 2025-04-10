#### Tools for Package Testing  --- in Matrix, sourced by ./test-tools.R
#### -------------------------
## to be used as, e.g.,
## source(system.file("test-tools-1.R", package="Matrix"), keep.source=FALSE)

### ------- Part I --  unrelated to "Matrix" classes ---------------

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)
identical5 <- function(a,b,c,d,e) identical(a,b) && identical4(b,c,d,e)
identical6 <- function(a,b,c,d,e,f)  identical(a,b) && identical5(b,c,d,e,f)
identical7 <- function(a,b,c,d,e,f,g)identical(a,b) && identical6(b,c,d,e,f,g)

require(tools)#-> assertError() and assertWarning()
assertWarningAtLeast <- function(expr, verbose=getOption("verbose"))
    tools::assertCondition(expr, "error", "warning", verbose=verbose)

##' [ from R's  demo(error.catching) ]
##' We want to catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' @title tryCatch both warnings and errors
##' @param expr
##' @return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' @author Martin Maechler
tryCatch.W.E <- function(expr)
{
    W <- NULL
    w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
				     warning = w.handler),
	 warning = W)
}


##' Is 'x' a valid object of class 'class' ?
isValid <- function(x, class) isTRUE(validObject(x, test=TRUE)) && is(x, class)

##' Signal an error (\code{\link{stop}}), if \code{x} is not a valid object
##' of class \code{class}.
##'
##' @title Stop if Not a Valid Object of Given Class
##' @param x any \R object
##' @param class character string specifying a class name
##' @return \emph{invisibly}, the value of \code{\link{validObject}(x)}, i.e.,
##'   \code{TRUE}; otherwise an error will have been signalled
##' @author Martin Maechler, March 2015
stopifnotValid <- function(x, class) {
    if(!is(x, class))
	stop(sprintf("%s is not of class \"%s\"",
		     deparse(substitute(x)), class), call. = FALSE)
    invisible(validObject(x))
}

## Some (sparse) Lin.Alg. algorithms return 0 instead of NA, e.g.
## qr.coef(<sparseQR>, y).
## For those cases, need to compare with a version where NA's are replaced by 0
mkNA.0 <- function(x) { x[is.na(x)] <- 0 ; x }

##' ... : further arguments passed to all.equal() such as 'check.attributes'
is.all.equal <- function(x,y, tol = .Machine$double.eps^0.5, ...)
    isTRUE(all.equal(x,y, tolerance=tol, ...))
is.all.equal3 <- function(x,y,z, tol = .Machine$double.eps^0.5, ...)
    is.all.equal(x,y, tol=tol, ...) && is.all.equal(y,z, tol=tol, ...)

is.all.equal4 <- function(x,y,z,u, tol = .Machine$double.eps^0.5, ...)
    is.all.equal3(x,y,z, tol=tol, ...) && isTRUE(all.equal(z,u, tolerance=tol, ...))

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
all.equal.X <- function(x,y, except, tol = .Machine$double.eps^0.5, ...)
{
    .trunc <- function(x) {
	ll <- as.list(x)
	ll[ - match(except, names(ll), nomatch = 0L)]
    }
    all.equal(.trunc(x), .trunc(y), tolerance = tol, ...)
}
## e.g. in lme4:
##  all.equal.X(env(m1), env(m2), except = c("call", "frame"))

##' @title The relative error typically returned by all.equal:
##' @note a copy of  sfsmisc::relErr()  unchanged since 2021-03-10
relErr <- function(target, current) { ## make this work, also for 'Matrix' ==> no mean() ..
    n <- length(current)
    if(length(target) < n) # (as we don't use mean())
        target <- rep(target, length.out = n)
    sum(abs(target - current)) / sum(abs(target))
}

##' Compute the signed relative error between target and current vector -- vectorized
##' @title Relative Error (:= 0 when absolute error == 0)
##' @param target  numeric, possibly scalar
##' @param current numeric of length() a multiple of length(target)
##' @return *vector* of the same length as current
##' @author Martin Maechler
##'
##' @note a copy of  sfsmisc::relErrV()  as of 2022-04-02
relErrV <- function(target, current, eps0 = .Machine$double.xmin) {
    n <- length(target <- as.vector(target))
    ## assert( <length current> is multiple of <length target>) :
    lc <- length(current)
    if(!n) {
	if(!lc) return(numeric()) # everything length 0
	else stop("length(target) == 0 differing from length(current)")
    } else if(!lc)
	stop("length(current) == 0 differing from length(target)")
    ## else n, lc  > 0
    if(lc %% n)
	stop("length(current) must be a multiple of length(target)")
    recycle <- (lc != n) # explicitly recycle
    R <- if(recycle)
	     target[rep(seq_len(n), length.out=lc)]
	 else
	     target # (possibly "mpfr")
    R[] <- 0
    ## use *absolute* error when target is zero {and deal with NAs}:
    t0 <- abs(target) < eps0 & !(na.t <- is.na(target))
    R[t0] <- current[t0]
    ## absolute error also when it is infinite, as (-Inf, Inf) would give NaN:
    dInf <- is.infinite(E <- current - target)
    R[dInf] <- E[dInf]
    useRE <- !dInf & !t0 & (na.t | is.na(current) | (current != target))
    R[useRE] <- (current/target)[useRE] - 1
    ## preserve {dim, dimnames, names}  from 'current' :
    if(!is.null(d <- dim(current)))
	array(R, dim=d, dimnames=dimnames(current))
    else if(!is.null(nm <- names(current)) && is.null(names(R))) # not needed for mpfr
	`names<-`(R, nm)
    else R
}


##' @title Number of correct digits: Based on relErrV(), recoding "Inf" to 'zeroDigs'
##' @param target  numeric vector of "true" values
##' @param current numeric vector of "approximate" values
##' @param zeroDigs how many correct digits should zero error give
##' @return basically   -log10 (| relErrV(target, current) | )
##' @author Martin Maechler, Summer 2011 (for 'copula')
nCorrDigits <- function(target, current, zeroDigs = 16) {
    stopifnot(zeroDigs >= -log10(.Machine$double.eps))# = 15.65
    RE <- relErrV(target, current)
    r <- -log10(abs(RE))
    r[RE == 0] <- zeroDigs
    r[is.na(RE) | r < 0] <- 0 # no correct digit, when relErr is NA
    r
}


## is.R22 <- (paste(R.version$major, R.version$minor, sep=".") >= "2.2")

pkgRversion <- function(pkgname)
    sub("^R ([0-9.]+).*", "\\1", packageDescription(pkgname)[["Built"]])

showSys.time <- function(expr, ...) {
    ## prepend 'Time' for R CMD Rdiff
    st <- system.time(expr, ...)
    writeLines(paste("Time", capture.output(print(st))))
    invisible(st)
}
showProc.time <- local({ ## function + 'pct' variable
    pct <- summary(proc.time())# length 3, shorter names
    function(final="\n", ind=TRUE) { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- summary(proc.time())
	delta <- (pct - ot)[ind]
	##  'Time' *not* to be translated:  tools::Rdiff() skips its lines!
	cat('Time', paste0("(",paste(names(delta),collapse=" "),"):"), delta, final)
    }
})

##' A version of sfsmisc::Sys.memGB() which should never give an error
##'  ( ~/R/Pkgs/sfsmisc/R/unix/Sys.ps.R  )
##' TODO: on Windows, with memory.size() & memory.limit() defunct, how do I get it ????
Sys.memGB <- function(kind = "MemTotal", ## "MemFree" is typically more relevant
                      NA.value = 2.10201) {
    if(!file.exists(pf <- "/proc/meminfo"))
	return(if(.Platform$OS.type == "windows")
                   NA.value ## memory.limit() / 1000  ## no longer with R 4.2.0
	       else
                   NA.value)
    mm <- tryCatch(drop(read.dcf(pf, fields=kind)),
                   error = function(e) NULL)
    if(is.null(mm) || any(is.na(mm)) || !all(grepl(" kB$", mm)))
        return(NA.value)
    ## return memory in giga bytes
    as.numeric(sub(" kB$", "", mm)) / (1000 * 1024)
}



##' @title turn an S4 object (with slots) into a list with corresponding components
##' @param obj an R object with a formal class (aka "S4")
##' @return a list with named components where \code{obj} had slots
##' @author Martin Maechler
S4_2list <- # <- "old" name (I like less: too hard to remember)
as.listS4 <- function(obj) {
   sn <- .slotNames(obj)
   ## structure(lapply(sn, slot, object = obj), .Names = sn)
   `names<-`(lapply(sn, slot, object = obj), sn)
}

assert.EQ <- function(target, current, tol = if(showOnly) 0 else 1e-15,
                      giveRE = FALSE, showOnly = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## showOnly: if TRUE, return (and hence typically print) all.equal(...)
    T <- isTRUE(ae <- all.equal(target, current, tolerance = tol, ...))
    if(showOnly) return(ae) else if(giveRE && T) { ## don't show if stop() later:
	ae0 <- if(tol == 0) ae else all.equal(target, current, tolerance = 0, ...)
	if(!isTRUE(ae0)) writeLines(ae0)
    }
    if(!T) stop("all.equal() |-> ", paste(ae, collapse=sprintf("%-19s","\n")))
    else if(giveRE) invisible(ae0)
}

##' a version with other "useful" defaults (tol, giveRE, check.attr..)
assert.EQ. <- function(target, current,
		       tol = if(showOnly) 0 else .Machine$double.eps^0.5,
		       giveRE = TRUE, showOnly = FALSE, ...) {
    assert.EQ(target, current, tol=tol, giveRE=giveRE, showOnly=showOnly,
	      check.attributes=FALSE, ...)
}

### ------- Part II  -- related to matrices, but *not* "Matrix" -----------

add.simpleDimnames <- function(m, named=FALSE) {
    stopifnot(length(d <- dim(m)) == 2)
    dimnames(m) <- setNames(list(if(d[1]) paste0("r", seq_len(d[1])),
				 if(d[2]) paste0("c", seq_len(d[2]))),
			    if(named) c("Row", "Col"))
    m
}

as.mat <- function(m) {
    ## as(., "matrix")	but with no extraneous empty dimnames
    d0 <- dim(m)
    m <- as(m, "matrix")
    if(!length(m) && is.null(d0)) dim(m) <- c(0L, 0L) # rather than (0, 1)
    if(identical(dimnames(m), list(NULL,NULL)))
	dimnames(m) <- NULL
    m
}

assert.EQ.mat <- function(M, m, tol = if(showOnly) 0 else 1e-15,
                          showOnly=FALSE, giveRE = FALSE, ...) {
    ## Purpose: check equality of  'Matrix' M with  'matrix' m
    ## ----------------------------------------------------------------------
    ## Arguments: M: is(., "Matrix") typically {but just needs working as(., "matrix")}
    ##            m: is(., "matrix")
    ##            showOnly: if TRUE, return (and hence typically print) all.equal(...)
    validObject(M)
    MM <- as.mat(M)                     # as(M, "matrix")
    if(is.logical(MM) && is.numeric(m))
	storage.mode(MM) <- "integer"
    attr(MM, "dimnames") <- attr(m, "dimnames") <- NULL
    assert.EQ(MM, m, tol=tol, showOnly=showOnly, giveRE=giveRE, ...)
}
## a short cut
assert.EQ.Mat <- function(M, M2, tol = if(showOnly) 0 else 1e-15,
                          showOnly=FALSE, giveRE = FALSE, ...)
    assert.EQ.mat(M, as.mat(M2), tol=tol, showOnly=showOnly, giveRE=giveRE, ...)

if(getRversion() <= "3.6.1" || R.version$`svn rev` < 77410)
    ## { methods::canCoerce() : use .class1(), not class() }
    canCoerce <- function(object, Class) {
        is(object, Class) ||
        !is.null(selectMethod("coerce", c(methods:::.class1(object), Class),
                              optional = TRUE,
                              useInherited = c(from=TRUE, to=FALSE)))
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
              rep(1, ncol(x)), tolerance = tol)
}

## .M.DN <- Matrix:::.M.DN -- but do *NOT* want to load Matrix namespace!
## from ../R/Auxiliaries.R :
`%||%` <- function(x, orElse) if(!is.null(x)) x else orElse
.M.DN <- function(x) dimnames(x) %||% list(NULL,NULL)
dnIdentical  <- function(x,y)   identical (.M.DN(x), .M.DN(y))
dnIdentical3 <- function(x,y,z) identical3(.M.DN(x), .M.DN(y), .M.DN(z))

##' @title Are two matrices practically equal - including dimnames
##' @param M1, M2: two matrices to be compared, maybe of _differing_ class
##' @param tol
##' @param dimnames logical indicating if dimnames must be equal
##' @param ... passed to all.equal(M1, M2)
##' @return TRUE or FALSE
is.EQ.mat <- function(M1, M2, tol = 1e-15, dimnames = TRUE, ...) {
    (if(dimnames) dnIdentical(M1,M2) else TRUE) &&
    is.all.equal(unname(as(M1, "matrix")),
		 unname(as(M2, "matrix")), tol=tol, ...)
}

##' see is.EQ.mat()
is.EQ.mat3 <- function(M1, M2, M3, tol = 1e-15, dimnames = TRUE, ...) {
    (if(dimnames) dnIdentical3(M1,M2,M3) else TRUE) &&
    is.all.equal3(unname(as(M1, "matrix")),
		  unname(as(M2, "matrix")),
		  unname(as(M3, "matrix")), tol=tol, ...)
}

##' here, as it also works for qr(<base matrix>)
chkQR <- function(a,
                  y = seq_len(nrow(a)),## RHS: made to contain no 0
                  a.qr = qr(a),
                  tol = 1e-11, # 1e-13 failing very rarely (interesting)
                  ##----------
                  Qinv.chk = !sp.rank.def,
                  QtQ.chk = !sp.rank.def,
                  verbose = getOption("Matrix.verbose", FALSE),
                  giveRE = verbose,
                  quiet = FALSE)
{
    d <- dim(a)
    stopifnot((n <- d[1]) >= (p <- d[2]), is.numeric(y))
    kind <- if(is.qr(a.qr)) "qr"
            else if(is(a.qr, "sparseQR")) "spQR"
            else stop("unknown qr() class: ", class(a.qr))
    if(!missing(verbose) && verbose) {
	op <- options(Matrix.verbose = verbose)
	on.exit(options(op))
    }
    ## rank.def <- switch(kind,
    ##     	       "qr"  = a.qr$rank < length(a.qr$pivot),
    ##     	       "spQR" = a.qr@V@Dim[1] > a.qr@Dim[1])
    sp.rank.def <- (kind == "spQR") && (a.qr@V@Dim[1] > a.qr@Dim[1])
    if(sp.rank.def && !quiet && (missing(Qinv.chk) || missing(QtQ.chk)))
	message("is sparse *structurally* rank deficient:  Qinv.chk=",
		Qinv.chk,", QtQ.chk=",QtQ.chk)
    if(is.na(QtQ.chk )) QtQ.chk  <- !sp.rank.def
    if(is.na(Qinv.chk)) Qinv.chk <- !sp.rank.def

    if(Qinv.chk) { ## qr.qy and qr.qty should be inverses,  Q'Q y = y = QQ' y :
        if(verbose) cat("Qinv.chk=TRUE: checking   Q'Q y = y = QQ' y :\n")
	## FIXME: Fails for structurally rank deficient sparse a's, but never for classical
	assert.EQ(drop(qr.qy (a.qr, qr.qty(a.qr, y))), y, giveRE=giveRE, tol = tol/64)
	assert.EQ(drop(qr.qty(a.qr, qr.qy (a.qr, y))), y, giveRE=giveRE, tol = tol/64)
    }

    piv <- switch(kind,
                  "qr" = a.qr$pivot,
                  "spQR" = 1L + a.qr@q)# 'q', not 'p' !!
    invP <- sort.list(piv)

    .ckQR <- function(cmpl) { ## local function, using parent's variables
        if(verbose) cat("complete = ",cmpl,": checking  X = Q R P*\n", sep="")
        Q <- qr.Q(a.qr, complete=cmpl) # NB: Q is already "back permuted"
        R <- qr.R(a.qr, complete=cmpl)
        rr <- if(cmpl) n else p
        stopifnot(dim(Q) == c(n,rr),
                  dim(R) == c(rr,p))
        assert.EQ.Mat(a, Q %*% R[, invP], giveRE=giveRE, tol=tol)
        ##            =  ===============
	if(QtQ.chk)
	    assert.EQ.mat(crossprod(Q), diag(rr), giveRE=giveRE, tol=tol)
        ##                ===========   ====
    }
    .ckQR(FALSE)
    .ckQR(TRUE)
    invisible(a.qr)
}## end{chkQR}
