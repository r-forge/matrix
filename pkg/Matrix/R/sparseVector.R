## METHODS FOR CLASS: sparseVector (virtual)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("mean", signature(x = "sparseVector"),
	  function(x, trim = 0, na.rm = FALSE, ...) {
              if(is.numeric(trim) && length(trim) == 1L && !is.na(trim) &&
                 trim == 0) {
                  ## Be fast in this special case :
                  if(isTRUE(na.rm))
                      ## FIXME: don't allocate !is.na(x)
                      x <- x[!is.na(x)]
                  sum(x) / length(x)
	      } else {
                  ## FIXME: don't allocate as.numeric(x); need 'sort' method
                  warning("suboptimally using as.numeric(x) to compute trimmed mean of sparseVector 'x'")
                  mean.default(as.numeric(x), trim = trim, na.rm = na.rm, ...)
              }
          })

setMethod("t", "sparseVector",
          function(x) spV2M(x, nrow = 1L, ncol = x@length, check = FALSE))

## MJ: unused
if(FALSE) {
## sort.default() does
##		x[order(x, na.last = na.last, decreasing = decreasing)]
## but that uses a *dense* integer order vector
## ==> need direct sort() method for "sparseVector" for mean(*,trim), median(),..
sortSparseV <- function(x, decreasing = FALSE, na.last = NA) {
    if(length(ina <- which(is.na(x)))) {
        if(is.na(na.last)) x <- x[-ina]
    }
    ## TODO
    .NotYetImplemented()
}

##' Uniquify sparceVectors, i.e., bring them in "regularized" from,
##' --- similar in spirit (and action!) as  uniqTsparse(.) for "TsparseMatrix"
##' __FIXME__ better name ??  , then export and document!  __TODO__
uniqSpVec <- function(x) {
    ii <- sort.list(x@i, method = "radix")
    x@i <- x@i[ii]
    x@x <- x@x[ii]
    x
}
} ## MJ

all.equal.sparseV <- function(target, current, ...)
{
    if(!is(target, "sparseVector") || !is(current, "sparseVector")) {
	return(paste0("target is ", data.class(target), ", current is ",
		      data.class(current)))
    }
    lt <- length(target)
    lc <- length(current)
    if(lt != lc) {
	return(paste0("sparseVector", ": lengths (", lt, ", ", lc, ") differ"))
    }

    t.has.x <- .hasSlot(target,  "x")## has "x" slot
    c.has.x <- .hasSlot(current, "x")## has "x" slot
    nz.t <- length(i.t <- target @i)
    nz.c <- length(i.c <- current@i)
    t.x <- if(t.has.x)	target@x else rep.int(TRUE, nz.t)
    c.x <- if(c.has.x) current@x else rep.int(TRUE, nz.c)
    if(nz.t != nz.c || any(i.t != i.c)) { ## "work" if indices are not the same
	i1.c <- setdiff(i.t, i.c)# those in i.t, not yet in i.c
	i1.t <- setdiff(i.c, i.t)
	if((n1t <- length(i1.t))) {
	    target@i <- i.t <- c(i.t, i1.t)
	    t.x <- c(t.x, rep.int(if(t.has.x) 0 else 0L, n1t))
	}
	if((n1c <- length(i1.c))) {
	    current@i <- i.c <- c(i.c, i1.c)
	    c.x <- c(c.x, rep.int(if(c.has.x) 0 else 0L, n1c))
	}
    }
    if(is.unsorted(i.t)) {  ## method="quick" {"radix" not ok for large range}
	ii <- sort.list(i.t, method = "quick", na.last=NA)
	target@i <- i.t <- i.t[ii]
	t.x <- t.x[ii]
    }
    if(is.unsorted(i.c)) {
	ii <- sort.list(i.c, method = "quick", na.last=NA)
	current@i <- i.c <- i.c[ii]
	c.x <- c.x[ii]
    }

    ## Now, we have extended both target and current
    ## *and* have sorted the respective i-slot, the i-slots should match!
    stopifnot(all(i.c == i.t))

    if(is.logical(t.x))
        all.equal.raw(t.x, c.x, ...)
    else
        all.equal.numeric(t.x, c.x, ...)
} ## all.equal.sparseV


## For these, we remain sparse:
setMethod("all.equal", c(target = "sparseVector", current = "sparseVector"),
	  all.equal.sparseV)
setMethod("all.equal", c(target = "sparseVector", current = "sparseMatrix"),
	  function(target, current, ...)
	  all.equal.sparseV(target, as(current, "sparseVector"), ...))
setMethod("all.equal", c(target = "sparseMatrix", current = "sparseVector"),
	  function(target, current, ...)
	  all.equal.sparseV(as(target, "sparseVector"), current, ...))
## For the others, where one is "dense", "go to" dense rather now than later:
setMethod("all.equal", c(target = "ANY", current = "sparseVector"),
	  function(target, current, ...)
	  all.equal(target, as.vector(current), ...))
setMethod("all.equal", c(target = "sparseVector", current = "ANY"),
	  function(target, current, ...)
	  all.equal(as.vector(target), current, ...))


## S3 method for 'c' [but only for dispatch on 1st arg, hence also exported as fn]
c.sparseVector <- function(...) {
    svl <- lapply(list(...), as, Class = "sparseVector")
    ## cls <- unique(unlist(lapply(svl, is)))
    ns <- vapply(svl, slot, 1, "length")
    if((N <- sum(ns)) < .Machine$integer.max) { # some purism ..
	ns <- as.integer(ns)
	N <- as.integer(N)
    }
    narg <- length(ns)
    iss <- lapply(svl, slot, "i")
    ## new 'i' slot:
    ii <- unlist(iss) + rep(cumsum(c(0L, ns[-narg])), lengths(iss))
    ## result must have 'x' slot if we have any
    has.x <- any(have.x <- vapply(svl, .hasSlot, logical(1L), name = "x"))
    if(has.x) {
	cls <- if     (any(vapply(svl, is, NA, "zsparseVector"))) "zsparseVector"
	       else if(any(vapply(svl, is, NA, "dsparseVector"))) "dsparseVector"
	       else if(any(vapply(svl, is, NA, "isparseVector"))) "isparseVector"
	       else "lsparseVector"
	if(!(all.x <- all(have.x)))
	    one <- if     (identical(cls, "lsparseVector")) TRUE
		   else if(identical(cls, "isparseVector")) 1L else 1.
	xx <- unlist(if(all.x) lapply(svl, slot, "x")
		     else lapply(seq_len(narg), function(i) {
			 if(have.x[[i]]) svl[[i]]@x
			 else rep_len(one, length(iss[[i]]))
		     }))
	new(cls, x = xx,     i = ii, length = N)
    } else ## no "x" slot
	new("nsparseVector", i = ii, length = N)
}

### rep(x, ...) -- rep() is primitive with internal default method with these args:
### -----------
### till R 2.3.1, it had  rep.default()  which we use as 'model' here.

repSpV <- function(x, times) {
    ## == rep.int(<sparseVector>, times)"
    times <- as.integer(times)# truncating as rep.default()
    n <- x@length
    has.x <- .hasSlot(x, "x")## has "x" slot
    ## just assign new correct slots:
    if(times <= 1) { ## be quick for {0, 1} times
        if(times < 0) stop("'times >= 0' is required")
        if(times == 0) {
            x@length <- 0L
            x@i <- integer(0)
            if(has.x) x@x <- rep.int(x@x, 0)
        }
        return(x)
    }
    n. <- as.double(n)
    if(n. * times >= .Machine$integer.max)
        n <- n. # so won't have overflow in subsequent multiplys
    x@length <- n * times
    x@i <- rep.int(x@i, times) + n * rep(0:(times-1L), each=length(x@i))
    ## := outer(x@i, 0:(times-1) * n, "+")   but a bit faster
    if(has.x) x@x <- rep.int(x@x, times)
    x
}

setMethod("rep", "sparseVector",
	  function(x, times, length.out, each, ...) {
	      if (length(x) == 0)
		  return(if(missing(length.out)) x else head(x, length.out))
	      if (!missing(each)) {
		  tm <- rep.int(each, length(x))
		  x <- rep(x, tm) # "recursively"
		  if(missing(length.out) && missing(times))
		      return(x)
	      } ## else :
	      if (!missing(length.out)) # takes precedence over times
		  times <- ceiling(length.out/length(x))
	      r <- repSpV(x, times)
	      if (!missing(length.out) && length(r) != length.out) {
		  if(length.out > 0) head(r, length.out) else r[integer(0)]
	      }
	      else r
	  })

##' indices of vector x[] to construct  Toeplitz matrix
##' FIXME: write in C, port to  R('stats' package), and use in stats::toeplitz()
ind4toeplitz <- function(n) {
    A <- matrix(raw(), n, n)
    abs(as.vector(col(A) - row(A))) + 1L
}

.toeplitz.spV <-  function(x, symmetric=TRUE, repr = c("C","T","R"), giveCsparse) {
    ## semantically "identical" to stats::toeplitz
    n <- length(x)
    r <- spV2M(x[ind4toeplitz(n)], n,n, symmetric = symmetric, check = FALSE)
    ##   ^^^^^ returning TsparseMatrix
    if(!missing(giveCsparse)) {
	if(missing(repr)) {
	    repr <- if(giveCsparse) "C" else "T"
	    warning(gettextf("'giveCsparse' has been deprecated; setting 'repr = \"%s\"' for you",
                             repr),
                    domain = NA)
	} else ## !missing(repr)
            Matrix.msg("'giveCsparse' has been deprecated; will use 'repr' instead")
    }
    switch(match.arg(repr), "C" = .M2C(r), "T" = r, "R" = .M2R(r))
}
setMethod("toeplitz", "sparseVector", .toeplitz.spV)
rm(.toeplitz.spV)
