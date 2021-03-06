### Define Methods that can be inherited for all subclasses

### Idea: Coercion between *VIRTUAL* classes -- as() chooses "closest" classes
### ----  should also work e.g. for  dense-triangular --> sparse-triangular !

##-> see als ./dMatrix.R, ./ddenseMatrix.R  and  ./lMatrix.R

setAs("ANY", "sparseMatrix", function(from) as(from, "CsparseMatrix"))

setAs("sparseMatrix", "generalMatrix", as_gSparse)

setAs("sparseMatrix", "symmetricMatrix", as_sSparse)

setAs("sparseMatrix", "triangularMatrix", as_tSparse)

spMatrix <- function(nrow, ncol,
                     i = integer(), j = integer(), x = numeric())
{
    dim <- c(as.integer(nrow), as.integer(ncol))
    ## The conformability of (i,j,x) with itself and with 'dim'
    ## is checked automatically by internal "validObject()" inside new(.):
    kind <- .M.kind(x)
    new(paste(kind, "gTMatrix", sep=''), Dim = dim,
        x = if(kind == "d") as.double(x) else x,
        ## our "Tsparse" Matrices use  0-based indices :
        i = as.integer(i - 1L),
        j = as.integer(j - 1L))
}

sparseMatrix <- function(i, j, p, x, dims, dimnames, index0 = FALSE)
{
  ## Purpose: user-level substitute for most  new(<sparseMatrix>, ..) calls
  ## Author: Martin Maechler, Date: 6 Jan 2009, based on Doug Bates' idea
    if((m.i <- missing(i)) + (m.j <- missing(j)) + (m.p <- missing(p)) != 1)
        stop("exactly one of 'i', 'j', or 'p' must be missing from call")
    if(!m.p) {
        ## we *could* let  validObject() / .validateCsparse() check later ...
        p <- as.integer(p)
        if((lp <- length(p)) < 1 || p[1] != 0 || any((dp <- p[-1] - p[-lp]) < 0))
            stop("'p' must be a non-decreasing vector (0, ...)")
    }
    isPat <- missing(x) ## <-> patter"n" Matrix
    ## "minimal dimensions" from (i,j,p) :
    dims.min <- c(if(m.i) lp - 1L else max(i <- as.integer(i)),
                  if(m.j) lp - 1L else max(j <- as.integer(j)))
    if(any(is.na(dims.min))) stop("NA's in (i,j) are not allowed")
    if(missing(dims)) {
        dims <- dims.min
    } else { ## check dims
        stopifnot(all(dims >= dims.min))
        dims <- as.integer(dims)
    }
    kx <- if(isPat) "n" else .M.kind(x)
    if(m.j) { ## -> Csparse
        r <- new(paste(kx, "gCMatrix", sep=''))
        r@Dim <- dims
        r@p <- as.integer(p)
        if(!isPat) r@x <- if(kx == "d" && !is.double(x)) as.double(x) else x
        r@i <- as.integer(if(index0) i else i - 1L)
        vv <- .validateCsparse(r, sort.if.needed=TRUE)## modify 'r' in-place !!!
        if(is.character(vv)) stop(vv)
    }
    else if(m.i) { ## -> Rsparse
        stop("(j,p) --> RsparseMatrix :  not yet implemented")
    }
    else if(m.p) { ## -> Tsparse
        r <- new(paste(kx, "gTMatrix", sep=''))
        r@Dim <- dims
        if(!isPat) r@x <- if(kx == "d" && !is.double(x)) as.double(x) else x
        r@i <- as.integer(if(index0) i else i - 1L)
        r@j <- as.integer(if(index0) j else j - 1L)
        validObject(r)
        r
    }
    if(!missing(dimnames)) {
        ## FIXME: should we check here, or validObject(r) or ??
        r@Dimnames <- dimnames
    }
    r
}

sparseMatrix <- function(i = ep, j = ep, p, x, dims, dimnames, index1 = TRUE)
{
  ## Purpose: user-level substitute for most  new(<sparseMatrix>, ..) calls
  ## Author: Douglas Bates, Date: 12 Jan 2009, based on Martin's version
    if((m.i <- missing(i)) + (m.j <- missing(j)) + (m.p <- missing(p)) != 1)
        stop("exactly one of 'i', 'j', or 'p' must be missing from call")
    if(!m.p) {
        p <- as.integer(p)
        if((lp <- length(p)) < 1 || p[1] != 0 || any((dp <- p[-1] - p[-lp]) < 0))
            stop("'p' must be a non-decreasing vector (0, ...)")
        ep <- rep.int(seq_along(dp), dp)
    }
    ## i and j are now both defined.  Make them 1-based indices.
    i1 <- as.logical(index1)[1]
    i <- as.integer(i + !(m.i || i1))
    j <- as.integer(j + !(m.j || i1))

    ## "minimal dimensions" from (i,j,p); no warnings from empty i or j :
    dims.min <- suppressWarnings(c(max(i), max(j)))
    if(any(is.na(dims.min))) stop("NA's in (i,j) are not allowed")
    if(missing(dims)) {
        dims <- dims.min
    } else { ## check dims
        stopifnot(all(dims >= dims.min))
        dims <- as.integer(dims)
    }
    isPat <- missing(x) ## <-> patter"n" Matrix
    kx <- if(isPat) "n" else .M.kind(x)
    r <- new(paste(kx, "gTMatrix", sep=''))
    r@Dim <- dims
    if(!isPat) {
	if(kx == "d" && !is.double(x)) x <- as.double(x)
	if(length(x) != (n <- length(i))) { ## recycle
	    if(length(x) != 1 && n %% length(x) != 0)
		warning("length(i) is not a multiple of length(x)")
	    x <- rep(x, length.out = n)
	}
	r@x <- x
    }
    r@i <- i - 1L
    r@j <- j - 1L
    validObject(r)
    as(r, "CsparseMatrix")
}

## "graph" coercions -- this needs the graph package which is currently
##  -----               *not* required on purpose
## Note: 'undirected' graph <==> 'symmetric' matrix

## Add some utils that may no longer be needed in future versions of the 'graph' package
graph.has.weights <- function(g) "weight" %in% names(edgeDataDefaults(g))

graph.wgtMatrix <- function(g)
{
    ## Purpose: work around "graph" package's  as(g, "matrix") bug
    ## ----------------------------------------------------------------------
    ## Arguments: g: an object inheriting from (S4) class "graph"
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006

    ## MM: another buglet for the case of  "no edges":
    if(numEdges(g) == 0) {
      p <- length(nd <- nodes(g))
      return( matrix(0, p,p, dimnames = list(nd, nd)) )
    }
    ## Usual case, when there are edges:
    has.w <- "weight" %in% names(edgeDataDefaults(g))
    if(has.w) {
        w <- unlist(edgeData(g, attr = "weight"))
        has.w <- any(w != 1)
    } ## now 'has.w' is TRUE  iff  there are weights != 1
    m <- as(g, "matrix")
    ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
    if(has.w) { ## fix it if needed
        tm <- t(m)
        tm[tm != 0] <- w
        t(tm)
    }
    else m
}


setAs("graphAM", "sparseMatrix",
      function(from) {
	  symm <- edgemode(from) == "undirected" && isSymmetric(from@adjMat)
	  ## This is only ok if there are no weights...
	  if(graph.has.weights(from)) {
	      as(graph.wgtMatrix(from),
		 if(symm) "dsTMatrix" else "dgTMatrix")
	  }
	  else { ## no weights: 0/1 matrix -> logical
	      as(as(from, "matrix"),
		 if(symm) "nsTMatrix" else "ngTMatrix")
	  }
      })

setAs("graph", "CsparseMatrix",
      function(from) as(as(from, "graphNEL"), "CsparseMatrix"))

setAs("graphNEL", "CsparseMatrix",
      function(from) as(as(from, "TsparseMatrix"), "CsparseMatrix"))

setAs("graphNEL", "TsparseMatrix",
      function(from) {
          nd <- nodes(from)
          dm <- rep.int(length(nd), 2)
	  symm <- edgemode(from) == "undirected"

 	  if(graph.has.weights(from)) {
	      eWts <- edgeWeights(from)
	      lens <- unlist(lapply(eWts, length))
	      i <- rep.int(0:(dm[1]-1), lens) # column indices (0-based)
	      To <- unlist(lapply(eWts, names))
	      j <- as.integer(match(To,nd) - 1L) # row indices (0-based)
	      ## symm <- symm && <weights must also be symmetric>: improbable
	      ## if(symm) new("dsTMatrix", .....) else
	      new("dgTMatrix", i = i, j = j, x = unlist(eWts),
		  Dim = dm, Dimnames = list(nd, nd))
	  }
 	  else { ## no weights: 0/1 matrix -> logical
              edges <- lapply(from@edgeL[nd], "[[", "edges")
              lens <- unlist(lapply(edges, length))
              ## nnz <- sum(unlist(lens))  # number of non-zeros
              i <- rep.int(0:(dm[1]-1), lens) # column indices (0-based)
              j <- as.integer(unlist(edges) - 1) # row indices (0-based)
              if(symm) {            # symmetric: ensure upper triangle
                  tmp <- i
                  flip <- i > j
                  i[flip] <- j[flip]
                  j[flip] <- tmp[flip]
                  new("nsTMatrix", i = i, j = j, Dim = dm,
                      Dimnames = list(nd, nd), uplo = "U")
              } else {
                  new("ngTMatrix", i = i, j = j, Dim = dm,
                      Dimnames = list(nd, nd))
              }
          }
      })

setAs("sparseMatrix", "graph", function(from) as(from, "graphNEL"))
setAs("sparseMatrix", "graphNEL",
      ## since have specific method for Tsparse below, are Csparse or Rsparse,
      ## i.e. do not need to "uniquify" the T* matrix:
      function(from) Tsp2grNEL(as(from, "TsparseMatrix"), need.uniq=FALSE))

Tsp2grNEL <- function(from, need.uniq = is_not_uniqT(from)) {
    d <- dim(from)
    if(d[1] != d[2])
	stop("only square matrices can be used as incidence matrices for graphs")
    n <- d[1]
    if(n == 0) return(new("graphNEL"))
    if(is.null(rn <- dimnames(from)[[1]]))
	rn <- as.character(1:n)
    if(need.uniq) ## Need to 'uniquify' the triplets!
        from <- uniq(from)

    if(isSymmetric(from)) { # either "symmetricMatrix" or otherwise
	##-> undirected graph: every edge only once!
	if(!is(from, "symmetricMatrix")) {
	    ## a general matrix which happens to be symmetric
	    ## ==> remove the double indices
	    from <- tril(from)
	}
        eMode <- "undirected"
    } else {
        eMode <- "directed"
    }
    ## every edge is there only once, either upper or lower triangle
    ft1 <- cbind(rn[from@i + 1L], rn[from@j + 1L])
    ## not yet: graph::ftM2graphNEL(.........)
    ftM2graphNEL(ft1, W = from@x, V= rn, edgemode= eMode)

}
setAs("TsparseMatrix", "graphNEL", function(from) Tsp2grNEL(from))


### Subsetting -- basic things (drop = "missing") are done in ./Matrix.R

### FIXME : we defer to the "*gT" -- conveniently, but not efficient for gC !

## [dl]sparse -> [dl]gT   -- treat both in one via superclass
##                        -- more useful when have "z" (complex) and even more

setMethod("[", signature(x = "sparseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i,j, ..., drop) {
	      cld <- getClassDef(class(x))
##> why should this be needed; can still happen in <Tsparse>[..]:
##>	      if(!extends(cld, "generalMatrix")) x <- as(x, "generalMatrix")
##	      viaCl <- paste(.M.kind(x, cld), "gTMatrix", sep='')
	      x <- as(x, "TsparseMatrix")[i, , drop=drop]
##simpler than x <- callGeneric(x = as(x, "TsparseMatrix"), i=i, drop=drop)
	      ## try_as(x, c(cl, sub("T","C", viaCl)))
	      if(is(x, "Matrix") && extends(cld, "CsparseMatrix"))
		  as(x, "CsparseMatrix") else x
	  })

setMethod("[", signature(x = "sparseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x,i,j, ..., drop) {
	      cld <- getClassDef(class(x))
##> why should this be needed; can still happen in <Tsparse>[..]:
##>	      if(!extends(cld, "generalMatrix")) x <- as(x, "generalMatrix")
##	      viaCl <- paste(.M.kind(x, cld), "gTMatrix", sep='')

	      x <- as(x, "TsparseMatrix")[, j, drop=drop]
##simpler than x <- callGeneric(x = as(x, "TsparseMatrix"), j=j, drop=drop)
	      if(is(x, "Matrix") && extends(cld, "CsparseMatrix"))
		  as(x, "CsparseMatrix") else x
	  })

setMethod("[", signature(x = "sparseMatrix",
			 i = "index", j = "index", drop = "logical"),
	  function (x, i, j, ..., drop) {
	      cld <- getClassDef(class(x))
	      ## be smart to keep symmetric indexing of <symm.Mat.> symmetric:
##>	      doSym <- (extends(cld, "symmetricMatrix") &&
##>			length(i) == length(j) && all(i == j))
##> why should this be needed; can still happen in <Tsparse>[..]:
##>	      if(!doSym && !extends(cld, "generalMatrix"))
##>		  x <- as(x, "generalMatrix")
##	      viaCl <- paste(.M.kind(x, cld),
##			     if(doSym) "sTMatrix" else "gTMatrix", sep='')
	      x <- as(x, "TsparseMatrix")[i, j, drop=drop]
	      if(is(x, "Matrix") && extends(cld, "CsparseMatrix"))
		  as(x, "CsparseMatrix") else x
	  })


## setReplaceMethod("[", .........)
## -> ./Tsparse.R
## &  ./Csparse.R  & ./Rsparse.R {those go via Tsparse}
##
## Do not use as.vector() (see ./Matrix.R ) for sparse matrices :
setReplaceMethod("[", signature(x = "sparseMatrix", i = "missing", j = "ANY",
				value = "sparseMatrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, , j=j, value = as(value, "sparseVector")))

setReplaceMethod("[", signature(x = "sparseMatrix", i = "ANY", j = "missing",
				value = "sparseMatrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, , value = as(value, "sparseVector")))

setReplaceMethod("[", signature(x = "sparseMatrix", i = "ANY", j = "ANY",
				value = "sparseMatrix"),
		 function (x, i, j, ..., value)
		 callGeneric(x=x, i=i, j=j, value = as(value, "sparseVector")))



## Group Methods

setMethod("Math",
	  signature(x = "sparseMatrix"),
	  function(x) callGeneric(as(x, "CsparseMatrix")))

## further group methods -> see ./Ops.R {"Summary": ./dMatrix.R }



### --- show() method ---

## FIXME(?) -- ``merge this'' (at least ``synchronize'') with
## - - -   prMatrix() from ./Auxiliaries.R
## FIXME: prTriang() in ./Auxiliaries.R  should also get  align = "fancy"
##
printSpMatrix <- function(x, digits = getOption("digits"),
		       maxp = getOption("max.print"), zero.print = ".",
		       col.names, note.dropping.colnames = TRUE,
		       col.trailer = '', align = c("fancy", "right"))
{
    cl <- getClassDef(class(x))
    stopifnot(extends(cl, "sparseMatrix"))
    validObject(x) # have seen seg.faults for invalid objects
    d <- dim(x)
    if(prod(d) > maxp) { # "Large" => will be "cut"
        ## only coerce to dense that part which won't be cut :
        nr <- maxp %/% d[2]
	m <- as(x[1:max(1, nr), ,drop=FALSE], "Matrix")
    } else {
        m <- as(x, "matrix")
    }
    dn <- dimnames(m) ## will be === dimnames(cx)
    logi <- (extends(cl,"lsparseMatrix") || extends(cl,"nsparseMatrix") ||
	     extends(cl, "pMatrix"))
    if(logi)
	cx <- array("N", dim(m), dimnames=dn)
    else { ## numeric (or --not yet implemented-- complex):
	cx <- apply(m, 2, format)
	if(is.null(dim(cx))) {# e.g. in	1 x 1 case
	    dim(cx) <- dim(m)
	    dimnames(cx) <- dn
	}
    }
    if (missing(col.names))
        col.names <- {
            if(!is.null(cc <- getOption("sparse.colnames")))
                cc
            else if(is.null(dn[[2]]))
                FALSE
            else { # has column names == dn[[2]]
                ncol(x) < 10
            }
        }
    if(identical(col.names, FALSE))
        cx <- emptyColnames(cx, msg.if.not.empty = note.dropping.colnames)
    else if(is.character(col.names)) {
        stopifnot(length(col.names) == 1)
        cn <- col.names
        switch(substr(cn, 1,3),
               "abb" = {
                   iarg <- as.integer(sub("^[^0-9]*", '', cn))
                   colnames(cx) <- abbreviate(colnames(cx), minlength = iarg)
               },
               "sub" = {
                   iarg <- as.integer(sub("^[^0-9]*", '', cn))
                   colnames(cx) <- substr(colnames(cx), 1, iarg)
               },
               stop("invalid 'col.names' string: ", cn))
    }
    ## else: nothing to do for col.names == TRUE
    if(is.logical(zero.print))
	zero.print <- if(zero.print) "0" else " "
    if(logi) {
	cx[!m] <- zero.print
	cx[m] <- "|"
    } else { # non logical
	## show only "structural" zeros as 'zero.print', not all of them..
	## -> cannot use 'm'
        d <- dim(cx)
	ne <- length(iN0 <- 1L + .Call(m_encodeInd, non0ind(x, cl), di = d))
	if(0 < ne && ne < prod(d)) {
	    align <- match.arg(align)
	    if(align == "fancy" && !is.integer(m)) {
		fi <- apply(m, 2, format.info) ## fi[3,] == 0  <==> not expo.
		if(R.version$`svn rev` < 48021) {
		    ## work around format.info() bug R <= 2.8.1 :
		    if(any(L <- fi[2,] > .Machine$integer.max - 5))
			fi[2:3,L] <- 0
                }
		## now 'format' the zero.print by padding it with ' ' on the right:
		## case 1: non-exponent:  fi[2,] + as.logical(fi[2,] > 0)
		## the column numbers of all 'zero' entries -- (*large*)
		cols <- 1L + (0:(prod(d)-1L))[-iN0] %/% d[1]
		pad <-
		    ifelse(fi[3,] == 0,
			   fi[2,] + as.logical(fi[2,] > 0),
			   ## exponential:
			   fi[2,] + fi[3,] + 4)
                ## now be efficient ; sprintf() is relatively slow
                ## and pad is much smaller than 'cols'; instead of "simply"
		## zero.print <- sprintf("%-*s", pad[cols] + 1, zero.print)
		if(any(doP <- pad > 0)) {#
		    ## only pad those that need padding - *before* expanding
		    z.p.pad <- rep.int(zero.print, length(pad))
		    z.p.pad[doP] <- sprintf("%-*s", pad[doP] + 1, zero.print)
		    zero.print <- z.p.pad[cols]
		}
                else
                    zero.print <- rep.int(zero.print, length(cols))
	    } ## else "right" : nothing to do

	    cx[-iN0] <- zero.print
	} else if (ne == 0)# all zeroes
	    cx[] <- zero.print
    }
    if(col.trailer != '')
        cx <- cbind(cx, col.trailer, deparse.level = 0)
    ## right = TRUE : cheap attempt to get better "." alignment
    print(cx, quote = FALSE, right = TRUE, max = maxp)
    invisible(x)
} ## printSpMatrix()

printSpMatrix2 <- function(x, digits = getOption("digits"),
                           maxp = getOption("max.print"), zero.print = ".",
                           col.names, note.dropping.colnames = TRUE,
                           suppRows = NULL, suppCols = NULL,
                           col.trailer = if(suppCols) "......" else "",
                           align = c("fancy", "right"))
{
    d <- dim(x)
    if((identical(suppRows,FALSE) && identical(suppCols, FALSE)) ||
       (!isTRUE(suppRows) && !isTRUE(suppCols) && prod(d) <= maxp))
    {
        if(missing(col.trailer) && is.null(suppCols))
            suppCols <- FALSE # for 'col.trailer'
        printSpMatrix(x, digits=digits, maxp=maxp,
                      zero.print=zero.print, col.names=col.names,
                      note.dropping.colnames=note.dropping.colnames,
                      col.trailer=col.trailer, align=align)
    } else { ## d[1] > maxp / d[2] >= nr : -- this needs [,] working:
	validObject(x)
	nR <- d[1] ## nrow
	useW <- getOption("width") - (format.info(nR)[1] + 3+1)
	##			     space for "[<last>,] "

	## --> suppress rows and/or columns in printing ...

	if(is.null(suppCols)) suppCols <- (d[2] * 2 > useW)
	nc <- if(suppCols) (useW - (1 + nchar(col.trailer))) %/% 2 else d[2]
	nr <- maxp %/% nc
	if(is.null(suppRows)) suppRows <- (nr < nR)

	sTxt <- c("in show(); maybe adjust 'options(max.print= *)'",
		  "\n ..............................\n")
	if(suppRows) {
	    if(suppCols)
		x <- x[ , 1:nc, drop = FALSE]
	    n2 <- ceiling(nr / 2)
	    printSpMatrix(x[seq_len(min(nR, max(1, n2))), , drop=FALSE],
			  digits=digits, maxp=maxp,
			  zero.print=zero.print, col.names=col.names,
			  note.dropping.colnames=note.dropping.colnames,
			  col.trailer = col.trailer, align=align)
	    cat("\n ..............................",
		"\n ........suppressing rows ", sTxt, "\n", sep='')
	    ## tail() automagically uses "[..,]" rownames:
	    printSpMatrix(tail(x, max(1, nr-n2)),
			  digits=digits, maxp=maxp,
			  zero.print=zero.print, col.names=col.names,
			  note.dropping.colnames=note.dropping.colnames,
			  col.trailer = col.trailer, align=align)
	}
	else if(suppCols) {
	    printSpMatrix(x[ , 1:nc , drop = FALSE],
			  digits=digits, maxp=maxp,
			  zero.print=zero.print, col.names=col.names,
			  note.dropping.colnames=note.dropping.colnames,
			  col.trailer = col.trailer, align=align)
	    cat("\n .....suppressing columns ", sTxt, sep='')
	}
	else stop("logic programming error in printSpMatrix2(), please report")

	invisible(x)
    }
} ## printSpMatrix2 ()

setMethod("print", signature(x = "sparseMatrix"), printSpMatrix)

setMethod("show", signature(object = "sparseMatrix"),
          function(object) {
              d <- dim(object)
              cl <- class(object)
              cat(sprintf('%d x %d sparse Matrix of class "%s"\n',
                          d[1], d[2], cl))
              printSpMatrix2(object)
          })


## For very large and very sparse matrices,  the above show()
## is not really helpful;  Use  summary() as an alternative:

setMethod("summary", signature(object = "sparseMatrix"),
	  function(object, ...) {
	      d <- dim(object)
	      T <- as(object, "TsparseMatrix")
	      ## return a data frame (int, int,	 {double|logical|...})	:
	      r <- if(is(object,"nsparseMatrix"))
		  data.frame(i = T@i + 1L, j = T@j + 1L)
	      else data.frame(i = T@i + 1L, j = T@j + 1L, x = T@x)
	      attr(r, "header") <-
		  sprintf('%d x %d sparse Matrix of class "%s", with %d entries',
			  d[1], d[2], class(object), length(T@i))
	      ## use ole' S3 technology for such a simple case
	      class(r) <- c("sparseSummary", class(r))
	      r
	  })

print.sparseSummary <- function (x, ...) {
    cat(attr(x, "header"),"\n")
    print.data.frame(x, ...)
    invisible(x)
}

setMethod("isSymmetric", signature(object = "sparseMatrix"),
	  function(object, tol = 100*.Machine$double.eps, ...) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)

	      ## else slower test using t()  --

	      ## FIXME (for tol = 0): use cholmod_symmetry(A, 1, ...)
	      ##        for tol > 0   should modify  cholmod_symmetry(..) to work with tol

	      ## or slightly simpler, rename and export	 is_sym() in ../src/cs_utils.c


	      if (is(object, "dMatrix"))
		  ## use gC; "T" (triplet) is *not* unique!
		  isTRUE(all.equal(.as.dgC.0.factors(  object),
				   .as.dgC.0.factors(t(object)),
				   tol = tol, ...))
	      else if (is(object, "lMatrix"))
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(	 object,  "lgCMatrix"),
			    as(t(object), "lgCMatrix"))
	      else if (is(object, "nMatrix"))
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(	 object,  "ngCMatrix"),
			    as(t(object), "ngCMatrix"))
	      else stop("not yet implemented")
	  })


setMethod("isTriangular", signature(object = "CsparseMatrix"), isTriC)
setMethod("isTriangular", signature(object = "TsparseMatrix"), isTriT)

setMethod("isDiagonal", signature(object = "sparseMatrix"),
	  function(object) {
              d <- dim(object)
              if(d[1] != d[2]) return(FALSE)
              ## else
	      gT <- as(object, "TsparseMatrix")
	      all(gT@i == gT@j)
	  })


setMethod("determinant", signature(x = "sparseMatrix", logarithm = "missing"),
	  function(x, logarithm, ...)
	  determinant(x, logarithm = TRUE, ...))
setMethod("determinant", signature(x = "sparseMatrix", logarithm = "logical"),
	  function(x, logarithm = TRUE, ...)
	  determinant(as(x,"dsparseMatrix"), logarithm, ...))


setMethod("diag", signature(x = "sparseMatrix"),
	  function(x, nrow, ncol) diag(as(x, "CsparseMatrix")))

setMethod("dim<-", signature(x = "sparseMatrix", value = "ANY"),
	  function(x, value) {
	      if(!is.numeric(value) || length(value) != 2)
		  stop("dim(.) value must be numeric of length 2")
	      if(prod(dim(x)) != prod(value <- as.integer(value)))
		  stop("dimensions don't match the number of cells")
              ## be careful to keep things sparse
	      as(spV2M(as(x, "sparseVector"), nrow=value[1], ncol=value[2]),
		 class(x))
	  })


setMethod("norm", signature(x = "sparseMatrix", type = "character"),
	  function(x, type, ...) {
	      type <- toupper(substr(type[1], 1, 1))
	      switch(type,  ##  max(<empty>, 0)  |-->  0
		     "O" = ,
                     "1" = max(colSums(abs(x)), 0), ## One-norm (L_1)
		     "I" = max(rowSums(abs(x)), 0), ## L_Infinity
		     "F" = sqrt(sum(x^2)), ## Frobenius
		     "M" = max(abs(x), 0), ## Maximum modulus of all
		     ## otherwise:
		     stop("invalid 'type'"))
	  })

setMethod("rcond", signature(x = "sparseMatrix", norm = "character"),
	  function(x, norm, ...) {
	      d <- dim(x)
              ## FIXME: qr.R(qr(.)) warns about differing R (permutation!)
              ##        really fix qr.R() *or* go via dense in any cases
	      rcond(if(d[1] == d[2]) {
			warning("rcond(.) via sparse -> dense coercion")
			as(x, "denseMatrix")
		    } else if(d[1] > d[2]) qr.R(qr(x)) else qr.R(qr(t(x))),
		    norm = norm, ...)
	  })

setMethod("cov2cor", signature(V = "sparseMatrix"),
	  function(V) {
	      ## like stats::cov2cor() but making sure all matrices stay sparse
	      p <- (d <- dim(V))[1]
	      if (p != d[2])
		  stop("'V' is not a *square* matrix")
	      if(!is(V, "dMatrix"))
		  V <- as(V, "dMatrix")# actually "dsparseMatrix"
	      Is <- sqrt(1/diag(V))
	      if (any(!is.finite(Is))) ## original had 0 or NA
		  warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
	      ## TODO: if  <diagonal> %*% <sparse> was implemented more efficiently
	      ##       we'd rather use that!
	      Is <- as(Diagonal(x = Is), "sparseMatrix")
	      r <- Is %*% V %*% Is
	      r[cbind(1L:p,1L:p)] <- 1 # exact in diagonal
	      as(r, "symmetricMatrix")
	  })

setMethod("is.na", signature(x = "sparseMatrix"),## NB: nsparse* have own method!
	  function(x) {
	      if(any((inax <- is.na(x@x)))) {
		  cld <- getClassDef(class(x))
		  if(extends(cld, "triangularMatrix") && x@diag == "U")
		      inax <- is.na((x <- .diagU2N(x, cld))@x)
		  r <- as(x, "lMatrix") # will be "lsparseMatrix" - *has* x slot
		  r@x <- inax
                  if(!extends(cld, "CsparseMatrix"))
                      r <- as(r, "CsparseMatrix")
		  as(.Call(Csparse_drop, r, 0), "nMatrix") # a 'pattern matrix
	      }
	      else is.na_nsp(x)
	  })

## all.equal(): similar to all.equal_Mat() in ./Matrix.R ;
## -----------	eventually defer to  "sparseVector" methods:
setMethod("all.equal", c(target = "sparseMatrix", current = "sparseMatrix"),
	  function(target, current, check.attributes = TRUE, ...)
      {
	  msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
	  if(is.list(msg)) return(msg[[1]])
	  ## else
	  r <- all.equal(as(target, "sparseVector"), as(current, "sparseVector"),
			 check.attributes=check.attributes, ...)
	  if(is.null(msg) && (r.ok <- isTRUE(r))) TRUE else c(msg, if(!r.ok) r)
      })
setMethod("all.equal", c(target = "sparseMatrix", current = "ANY"),
	  function(target, current, check.attributes = TRUE, ...)
      {
	  msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
	  if(is.list(msg)) return(msg[[1]])
	  ## else
	  r <- all.equal(as(target, "sparseVector"), current,
			 check.attributes=check.attributes, ...)
	  if(is.null(msg) && (r.ok <- isTRUE(r))) TRUE else c(msg, if(!r.ok) r)
      })
setMethod("all.equal", c(target = "ANY", current = "sparseMatrix"),
	  function(target, current, check.attributes = TRUE, ...)
      {
	  msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
	  if(is.list(msg)) return(msg[[1]])
	  ## else
	  r <- all.equal(target, as(current, "sparseVector"),
			 check.attributes=check.attributes, ...)
	  if(is.null(msg) && (r.ok <- isTRUE(r))) TRUE else c(msg, if(!r.ok) r)
      })



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

## xtabs returning a sparse matrix.  This is cut'n'paste
## of xtabs() in <Rsrc>/src/library/stats/R/xtabs.R ;
## with the new argument 'sparse'
xtabs <- function(formula = ~., data = parent.frame(), subset, sparse = FALSE,
		  na.action, exclude = c(NA, NaN), drop.unused.levels = FALSE)
{
    if (missing(formula) && missing(data))
	stop("must supply either 'formula' or 'data'")
    if(!missing(formula)){
	## We need to coerce the formula argument now, but model.frame
	## will coerce the original version later.
	formula <- as.formula(formula)
	if (!inherits(formula, "formula"))
	    stop("'formula' missing or incorrect")
    }
    if (any(attr(terms(formula, data = data), "order") > 1))
	stop("interactions are not allowed")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
	m$data <- as.data.frame(data)
    m$... <- m$exclude <- m$drop.unused.levels <- m$sparse <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    if(length(formula) == 2) {
	by <- mf
	y <- NULL
    }
    else {
	i <- attr(attr(mf, "terms"), "response")
	by <- mf[-i]
	y <- mf[[i]]
    }
    by <- lapply(by, function(u) {
	if(!is.factor(u)) u <- factor(u, exclude = exclude)
	u[ , drop = drop.unused.levels]
    })
    if(!sparse) { ## identical to stats::xtabs
	x <-
	    if(is.null(y))
		do.call("table", by)
	    else if(NCOL(y) == 1)
		tapply(y, by, sum)
	    else {
		z <- lapply(as.data.frame(y), tapply, by, sum)
		array(unlist(z),
		      dim = c(dim(z[[1]]), length(z)),
		      dimnames = c(dimnames(z[[1]]), list(names(z))))
	    }
	x[is.na(x)] <- 0
	class(x) <- c("xtabs", "table")
	attr(x, "call") <- match.call()
	x

    } else { ## sparse
	if (length(by) != 2)
	    stop("xtabs(*, sparse=TRUE) applies only to two-way tables")
	rows <- by[[1]]
	cols <- by[[2]]
	rl <- levels(rows)
	cl <- levels(cols)
	if (is.null(y))
	    y <- rep.int(1, length(rows))
	as(new("dgTMatrix",
	       i = as.integer(rows) - 1L,
	       j = as.integer(cols) - 1L,
	       x = as.double(y),
	       Dim = c(length(rl), length(cl)),
	       Dimnames = list(rl, cl)), "CsparseMatrix")
    }
}
