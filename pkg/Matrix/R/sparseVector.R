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
