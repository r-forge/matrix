#### All Methods in relation with the sparseVector (sub)classes


# atomicVector : classUnion (logical,integer,double,....)
setAs("atomicVector", "sparseVector",
      function(from) {
	  n <- length(from)
	  r <- new(paste(.V.kind(from), "sparseVector", sep=''), length = n)
	  ii <- from != 0
	  r@x <- from[ii]
	  r@i <- seq_len(n)[ii]
	  r
      })


for(T in c("d","i","l","z")) {
    setAs("xsparseVector", paste(T, "sparseVector", sep=''),
          function(from) {
              from@x <- as(from@x, .type.kind[T])
              ## and now "the hack":
              class(from) <- paste(T, "sparseVector", sep='')
              from
          })
}

setAs("sparseVector", "nsparseVector",
      function(from) {
          if(any(is.na(from@x)))
              stop("cannot coerce 'NA's to \"nsparseVector\"")
          new("nsparseVector", i = from@i, length = from@length)
      })

sp2vec <- function(x, mode = .type.kind[substr(cl, 1,1)]) {
    ## sparseVector  ->  vector
    cl <- class(x)
    r <- vector(mode, x@length)
    r[x@i] <-
	if(cl != "nsparseVector") { # cheap test for 'has x slot'
	    if(is(x@x, mode)) x@x else as(x@x, mode)
	} else TRUE
    r
}

setAs("sparseVector", "vector", function(from) sp2vec(from))

setMethod("as.vector", signature(x = "sparseVector", mode = "missing"),
	  sp2vec)
setMethod("as.vector", signature(x = "sparseVector", mode = "character"),
	  sp2vec)

setMethod("as.numeric", "sparseVector", function(x) sp2vec(x, mode = "double"))

## the "catch all remaining" method:
setAs("ANY", "sparseVector",
      function(from) as(as.vector(from), "sparseVector"))

setAs("diagonalMatrix", "sparseVector",
      function(from) {
	  kind <- .M.kind(from) ## currently only "l" and "d" --> have 'x'
	  n <- nrow(from)
          n2 <- as.double(n) * n
	  if(n2 > .Machine$integer.max) { ## double (i, length)
	      ii <- seq(1, by = n+1, length.out = n) ## 1-based indexing
	  } else { # integer ok
	      n2 <- as.integer(n2)
	      ii <- as.integer(seq(1L, by = n+1L, length.out = n))
	  }
	  new(paste(kind, "sparseVector", sep=''),
	      length = n2, i = ii,
	      x = if(from@diag != "U") from@x else
		  rep.int(switch(kind, "d" = 1, "l" = TRUE, "i" = 1L, "z" = 1+0i), n))
	 })

setAs("sparseMatrix", "sparseVector",
      function(from) as(as(from, "TsparseMatrix"), "sparseVector"))

setAs("TsparseMatrix", "sparseVector",
      function(from) {
	  d <- dim(from)
	  n <- prod(d) # -> numeric, no integer overflow
	  kind <- .M.kind(from)
	  if(is_duplicatedT(from, di = d))
	      from <- uniqTsparse(from)
	  r <- new(paste(kind, "sparseVector", sep=''), length = n)
	  r@i <- if(n < .Machine$integer.max) {
	      1L + from@i + d[1] * from@j
	  } else {
	      1 + from@i + as.double(d[1]) * from@j
	  }
	  if(kind != "n") ## have 'x' slot
	      r@x <- from@x
	  r
      })



## TODO -- also want  (sparseVector, dim) |---> sparseMatrix
##  because of (nrow,ncol) specification can not (?)  use as(.).
##  Hence use  Matrix(.) ?  or my  spMatrix(.) ?

## For now, define this utility function:
spV2M <- function (x, nrow, ncol, byrow = FALSE)
{
    ## Purpose:	 sparseVector --> sparseMatrix	constructor
    ## ----------------------------------------------------------------------
    ## Arguments: x: "sparseVector" object
    ##		nrow, ncol, byrow: as for matrix() or Matrix()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 May 2007

    cx <- class(x)
    stopifnot(extends(cx, "sparseVector"))
    if(!missing(ncol)) { ncol <- as.integer(ncol)
			 if(ncol <= 0) stop("'ncol' must be >= 1") }
    if(!missing(nrow)) { nrow <- as.integer(nrow)
			 if(nrow <= 0) stop("'nrow' must be >= 1") }
    n <- length(x)
    if(missing(nrow)) {
	if(missing(ncol)) { ## both missing: --> (n x 1)
	    ncol <- 1L
	    nrow <- n
	} else {
	    if(n %% ncol != 0) warning("'ncol' is not a factor of length(x)")
	    nrow <- as.integer(ceiling(n / ncol))
	}
    } else {
	if(missing(ncol)) {
	    if(n %% nrow != 0) warning("'nrow' is not a factor of length(x)")
	    ncol <- as.integer(ceiling(n / nrow))
	} else { ## both nrow and ncol specified
	    n.n <- as.double(ncol) * nrow # no integer overflow
	    if(n.n <  n) stop("nrow * ncol < length(x)")
	    if(n.n != n) warning("nrow * ncol != length(x)")
	}
    }
    ## now nrow * ncol >= n
    ##	   ~~~~~~~~~~~~~~~~
    cld <- getClassDef(cx)
    kind <- .M.kindC(cld)		# "d", "n", "l", "i", "z", ...
    has.x <- kind != "n"
    ## "careful_new()" :
    cNam <- paste(kind, "gTMatrix", sep='')
    chngCl <- is.null(newCl <- getClass(cNam, .Force=TRUE))
    if(chngCl) { ## e.g. "igTMatrix" is not yet implemented
	if(substr(cNam,1,1) == "z")
	    stopifnot(sprintf("Class '%s' is not yet implemented", cNam))
	## coerce to "double":
	newCl <- getClass("dgTMatrix")
    }
    r <- new(newCl, Dim = c(nrow, ncol))
    ## now "compute"  the (i,j,x) slots given x@(i,x)
    i0 <- x@i - 1L
    if(byrow) { ## need as.integer(.) since <sparseVector> @ i can be double
	r@j <- as.integer(i0 %% ncol)
	r@i <- as.integer(i0 %/% ncol)
    } else {				# default{byrow = FALSE}
	r@i <- as.integer(i0 %% nrow)
	r@j <- as.integer(i0 %/% nrow)
    }
    if(has.x)
	r@x <- if(chngCl) as.numeric(x@x) else x@x
    r
}

## This is very similar to the 'x = "sparseMatrix"' method in ./sparseMatrix.R:
setMethod("dim<-", signature(x = "sparseVector", value = "ANY"),
	  function(x, value) {
	      if(!is.numeric(value) || length(value) != 2)
		  stop("dim(.) value must be numeric of length 2")
	      if(length(x) != prod(value <- as.integer(value)))
		  stop("dimensions don't match the number of cells")
	      spV2M(x, nrow=value[1], ncol=value[2])
	  })


setMethod("length", "sparseVector", function(x) x@length)

setMethod("show", signature(object = "sparseVector"),
   function(object) {
       n <- object@length
       cl <- class(object)
       cat(sprintf('sparse vector (nnz/length = %d/%.0f) of class "%s"\n',
		   length(object@i), as.double(n), cl))
       maxp <- max(1, getOption("max.print"))
       if(n <= maxp) {
	   prSpVector(object, maxp = maxp)
       } else { # n > maxp : will cut length of what we'll display :
	   ## cannot easily show head(.) & tail(.) because of "[1] .." printing of tail
	   prSpVector(object[seq_len(maxp)], maxp = maxp)
	   cat(" ............................",
	       "\n ........suppressing ", n - maxp,
	       " entries in show(); maybe adjust 'options(max.print= *)'",
	       "\n ............................\n\n", sep='')
       }
       invisible(object)
   })

prSpVector <- function(x, digits = getOption("digits"),
		    maxp = getOption("max.print"), zero.print = ".")
{
    cld <- getClassDef(cl <- class(x))
    stopifnot(extends(cld, "sparseVector"), maxp >= 1)
    if(is.logical(zero.print))
	zero.print <- if(zero.print) "0" else " "
##     kind <- .M.kindC(cld)
##     has.x <- kind != "n"
    n <- x@length
    if(n > 0) {
        if(n > maxp) { # n > maxp =: nn : will cut length of what we'll display :
            ## FIXME: very inefficient for very large maxp < n
            x <- x[seq_len(maxp)]       # need "[" to work ...
            n <- maxp
        }
        xi <- x@i
        logi <- extends(cld, "lsparseVector") || extends(cld, "nsparseVector")
        cx <- if(logi) rep.int("N", n) else character(n)
        cx[ -xi ] <- zero.print
        cx[  xi ] <- {
            if(logi) "|" else
            ## numeric (or --not yet-- complex): 'has.x' in any cases
            format(x@x, digits = digits)
        }
        ## right = TRUE : cheap attempt to get better "." alignment
        print(cx, quote = FALSE, right = TRUE, max = maxp)
    }
    invisible(x) # TODO? in case of n > maxp, "should" return original x
}

## This is a simplified intI() {-> ./Tsparse.R } -- for sparseVector indexing:
intIv <- function(i, n)
{
    ## Purpose: translate numeric | logical index     into  1-based integer
    ## --------------------------------------------------------------------
    ## Arguments: i: index vector (numeric | logical)
    ##		  n: array extent { ==	length(.) }
    if(missing(i))
	seq_len(n)
    else if(is(i, "numeric")) {
	storage.mode(i) <- "integer"
	if(any(i < 0L)) {
	    if(any(i > 0L))
		stop("you cannot mix negative and positive indices")
	    seq_len(n)[i]
	} else {
	    if(length(i) && max(i) > n)
		stop("indexing out of range 0:",n)
	    if(any(z <- i == 0))
		i <- i[!z]
	    i
	}
    }
    else if (is(i, "logical")) {
	seq_len(n)[i]
    } else stop("index must be numeric or logical for 'sparseVector' indexing")
}


setMethod("[", signature(x = "sparseVector", i = "index"),
	  function (x, i, j, ..., drop) {
	      cld <- getClassDef(class(x))
	      has.x <- !extends(cld, "nsparseVector")
	      n <- x@length
	      ii <- intIv(i, n)
	      anyDup <- any(iDup <- duplicated(ii))
	      m <- match(x@i, ii, nomatch = 0)
	      sel <- m > 0L
	      x@length <- length(ii)
	      x@i <- m[sel]
	      if(anyDup) {
		  i.i <- match(ii[iDup], ii)
		  jm <- lapply(i.i, function(.) which(. == m))
		  sel <- c(which(sel), unlist(jm))
		  x@i <- c(x@i, rep.int(which(iDup), sapply(jm, length)))
	      }
	      if (has.x)
		  x@x <- x@x[sel]
	      x
	  })

## This is much analogous to replTmat in ./Tsparse.R:
replSPvec <- function (x, i, value)
{
    n <- x@length
    ii <- intIv(i, n)
    lenRepl <- length(ii)
    lenV <- length(value)
    if(lenV == 0) {
	if(lenRepl != 0)
	    stop("nothing to replace with")
	else return(x)
    }
    ## else: lenV := length(value) > 0
    if(lenRepl %% lenV != 0)
	stop("number of items to replace is not a multiple of replacement length")
    anyDup <- any(duplicated(ii))
    if(anyDup) { ## multiple *replacement* indices: last one wins
	## TODO: in R 2.6.0 use	 duplicate(*, fromLast=TRUE)
	ir <- lenRepl:1
	keep <- match(ii, ii[ir]) == ir
	ii <- ii[keep]
	lenV <- length(value <- rep(value, length = lenRepl)[keep])
	lenRepl <- length(ii)
    }

    cld <- getClassDef(class(x))
    has.x <- !extends(cld, "nsparseVector")
    m <- match(x@i, ii, nomatch = 0)
    sel <- m > 0L

    ## the simplest case
    if(all0(value)) { ## just drop the non-zero entries
	if(any(sel)) { ## non-zero there
	    x@i <- x@i[!sel]
	    if(has.x)
		x@x <- x@x[!sel]
	}
	return(x)

    }
    ## else --	some( value != 0 ) --
    if(lenV > lenRepl)
	stop("too many replacement values")
    else if(lenV < lenRepl)
	value <- rep(value, length = lenRepl)
    ## now:  length(value) == lenRepl

    v0 <- is0(value)
    ## value[1:lenRepl]:  which are structural 0 now, which not?

    if(any(sel)) {
	## indices of non-zero entries -- WRT to subvector
	iN0 <- m[sel] ## == match(x@i[sel], ii)

	## 1a) replace those that are already non-zero with new val.
	vN0 <- !v0[iN0]
	if(any(vN0) && has.x)
	    x@x[sel][vN0] <- value[iN0[vN0]]

	## 1b) replace non-zeros with 0 --> drop entries
	if(any(!vN0)) {
	    i <- which(sel)[!vN0]
	    if(has.x)
		x@x <- x@x[-i]
	    x@i <- x@i[-i]
	}
	iI0 <- if(length(iN0) < lenRepl)
	    seq_len(lenRepl)[-iN0]
    } else iI0 <- seq_len(lenRepl)

    if(length(iI0) && any(vN0 <- !v0[iI0])) {
	## 2) add those that were structural 0 (where value != 0)
	ij0 <- iI0[vN0]
	x@i <- c(x@i, ii[ij0])
	if(has.x)
	    x@x <- c(x@x, value[ij0])
    }
    x

}

setReplaceMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
				value = "replValue"),
		 replSPvec)



## a "method" for c(<(sparse)Vector>, <(sparse)Vector>):
## FIXME: This is not exported, nor used (nor documented)
c2v <- function(x, y) {
    ## these as(., "sp..V..") check input implicitly:
    cx <- class(x <- as(x, "sparseVector"))
    cy <- class(y <- as(y, "sparseVector"))
    if(cx != cy) { ## find "common" class; result does have 'x' slot
        cxy <- c(cx,cy)
        commType <- {
            if(all(cxy %in% c("nsparseVector", "lsparseVector")))
                "lsparseVector"
            else { # ==> "numeric" ("integer") or "complex"
                xslot1 <- function(u, cl.u)
                    if(cl.u != "nsparseVector") u@x[1] else TRUE
                switch(typeof(xslot1(x, cx) + xslot1(y, cy)),
                       ## "integer", "double", or "complex"
                       "integer" = "isparseVector",
                       "double" = "dsparseVector",
                       "complex" = "zsparseVector")
            }
        }
        if(cx != commType) x <- as(x, commType)
        if(cy != commType) y <- as(y, commType)
        cx <- commType
    }
    ## now *have* common type -- transform 'x' into result:
    nx <- x@length
    x@length <- nx + y@length
    x@i <- c(x@i, nx + y@i)
    if(cx != "nsparseVector")
        x@x <- c(x@x, y@x)
    x
}

all.equal.sparseV <- function(target, current, ...)
{
    if(!is(target, "sparseVector") || !is(current, "sparseVector")) {
	return(paste("target is ", data.class(target), ", current is ",
		     data.class(current), sep = ""))
    }
    lt <- length(target)
    lc <- length(current)
    if(lt != lc) {
	return(paste("sparseVector", ": lengths (", lt, ", ", lc, ") differ",
		     sep = ""))
    }

    t.has.x <- class(target)  != "nsparseVector"
    c.has.x <- class(current) != "nsparseVector"
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

    all.equal.numeric(c.x, t.x, ...)
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


### rep(x, ...) -- rep() is primitive with internal default method with these args:
### -----------
### till R 2.3.1, it had  rep.default()  which we use as 'model' here.

repSpV <- function(x, times) {
    ## == rep.int(<sparseVector>, times)"
    times <- as.integer(times)# truncating as rep.default()
    n <- x@length
    has.x <- substr(class(x), 1,1) != "n" ## fast, but hackish
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
		  return(if(missing(length.out)) x else x[seq_len(length.out)])
	      if (!missing(each)) {
		  tm <- rep.int(each, length(x))
		  x <- rep(x, tm) # "recursively"
		  if(missing(length.out) && missing(times))
		      return(x)
	      } ## else :
	      if (!missing(length.out)) # takes precedence over times
		  times <- ceiling(length.out/length(x))
	      r <- repSpV(x, times)
	      if (!missing(length.out))
		  return(r[if(length.out > 0) 1:length.out else integer(0)])
	      return(r)
	  })


### Group Methods (!)

## o "Ops" , "Arith", "Compare"  :  ---> in ./Ops.R

