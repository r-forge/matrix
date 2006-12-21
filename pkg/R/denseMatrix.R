### Simple fallback methods for all dense matrices
### These are "cheap" to program, but potentially far from efficient;
### Methods for specific subclasses will overwrite these:

setAs("ANY", "denseMatrix", function(from) Matrix(from, sparse=FALSE))


## dense to sparse:
setAs("denseMatrix", "dsparseMatrix",
## MM thought that  as() will take the ``closest'' match; but that fails!
##      function(from) as(as(from, "dgeMatrix"), "dsparseMatrix"))
      function(from) as(as(from, "dgeMatrix"), "dgCMatrix"))

setAs("denseMatrix", "CsparseMatrix",
      function(from) {
          cl <- class(from)
	  notGen <- !is(from, "generalMatrix")
	  if (notGen) { ## e.g. for triangular | symmetric
              ## FIXME: this is a *waste* in the case of packed matrices!
	      if     (extends(cl, "dMatrix")) from <- as(from, "dgeMatrix")
	      else if(extends(cl, "nMatrix")) from <- as(from, "ngeMatrix")
	      else if(extends(cl, "lMatrix")) from <- as(from, "lgeMatrix")
	      else if(extends(cl, "zMatrix")) from <- as(from, "zgeMatrix")
	      else stop("undefined method for class ", cl)
	  }
          ## FIXME: contrary to its name, this only works for "dge*" :
	  .Call(dense_to_Csparse, from)
      })

setAs("denseMatrix", "TsparseMatrix",
      function(from) as(as(from, "CsparseMatrix"), "TsparseMatrix"))


setMethod("show", signature(object = "denseMatrix"),
          function(object) prMatrix(object))
##- ## FIXME: The following is only for the "dMatrix" objects that are not
##- ##	      "dense" nor "sparse" -- i.e. "packed" ones :
##- ## But these could be printed better -- "." for structural zeros.
##- setMethod("show", signature(object = "dMatrix"), prMatrix)
##- ## and improve this as well:
##- setMethod("show", signature(object = "pMatrix"), prMatrix)
##- ## this should now be superfluous [keep for safety for the moment]:

## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

## use geClass() when 'i' or 'j' are missing:
## since  symmetric, triangular, .. will not be preserved anyway:
setMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, drop) {
	      r <- as(x, "matrix")[i, , drop=drop]
	      if(is.null(dim(r))) r else as(r, geClass(x))
	  })

setMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, j, drop) {
	      r <- as(x, "matrix")[, j, drop=drop]
	      if(is.null(dim(r))) r else as(r, geClass(x))
	  })

setMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
			 drop = "logical"),
	  function (x, i, j, drop) {
	      r <- callGeneric(x = as(x, "matrix"), i=i, j=j, drop=drop)
	      if(is.null(dim(r)))
		  r
	      else {
		  cl <- class(x)
		  if(extends(cl, "symmetricMatrix") &&
		     length(i) == length(j) && all(i == j))
		      as(r, cl) ## keep original symmetric class
		  else as_geClass(r, cl)
	      }
	  })

## Now the "[<-" ones --- see also those in ./Matrix.R
## It's recommended to use setReplaceMethod() rather than setMethod("[<-",.)
## even though the former is currently just a wrapper for the latter

## FIXME: 1) These are far from efficient
## -----  2) value = "numeric" is only ok for "ddense*"
setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
				value = "replValue"),
		 function (x, i, value) {
		     r <- as(x, "matrix")
		     r[i, ] <- value
		     as(r, geClass(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
				value = "replValue"),
		 function (x, j, value) {
		     r <- as(x, "matrix")
		     r[, j] <- value
		     as(r, geClass(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
				value = "replValue"),
		 function (x, i, j, value) {
		     r <- as(x, "matrix")
		     r[i, j] <- value
		     as_geClass(r, class(x)) ## was as(r, class(x))
		 })


setMethod("isSymmetric", signature(object = "denseMatrix"),
	  function(object, tol = 100*.Machine$double.eps) {
	      ## pretest: is it square?
	      d <- dim(object)
	      if(d[1] != d[2]) return(FALSE)
	      ## else slower test
	      if (is(object,"dMatrix"))
		  isTRUE(all.equal(as(object, "dgeMatrix"),
				   as(t(object), "dgeMatrix"), tol = tol))
	      else if (is(object, "nMatrix"))
		  identical(as(object, "ngeMatrix"),
			    as(t(object), "ngeMatrix"))
	      else if (is(object, "lMatrix"))# not possible currently
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(object, "lgeMatrix"),
			    as(t(object), "lgeMatrix"))
	      else if (is(object, "zMatrix"))
                  stop("'zMatrix' not yet implemented")
	      else if (is(object, "iMatrix"))
                  stop("'iMatrix' not yet implemented")
	  })

setMethod("isTriangular", signature(object = "triangularMatrix"),
	  function(object, ...) TRUE)

setMethod("isTriangular", signature(object = "denseMatrix"), isTriMat)

setMethod("isDiagonal", signature(object = "denseMatrix"), .is.diagonal)

.as.dge.Fun <- function(x, na.rm = FALSE, dims = 1) {
    x <- as(x, "dgeMatrix")
    callGeneric()
}
setMethod("colSums",  signature(x = "denseMatrix"), .as.dge.Fun)
setMethod("colMeans", signature(x = "denseMatrix"), .as.dge.Fun)
setMethod("rowSums",  signature(x = "denseMatrix"), .as.dge.Fun)
setMethod("rowMeans", signature(x = "denseMatrix"), .as.dge.Fun)

### cbind2
setMethod("cbind2", signature(x = "denseMatrix", y = "numeric"),
	  function(x, y) {
	      d <- dim(x); nr <- d[1]; nc <- d[2]
	      y <- rep(y, length.out = nr) # 'silent procrustes'
	      ## beware of (packed) triangular, symmetric, ...
	      x <- as(x, geClass(x))
	      x@x <- c(x@x, as.double(y))
	      x@Dim[2] <- nc + 1:1
	      if(is.character(dn <- x@Dimnames[[2]]))
		  x@Dimnames[[2]] <- c(dn, "")
	      x
	  })
## the same, (x,y) <-> (y,x):
setMethod("cbind2", signature(x = "numeric", y = "denseMatrix"),
	  function(x, y) {
	      d <- dim(y); nr <- d[1]; nc <- d[2]
	      x <- rep(x, length.out = nr)
	      y <- as(y, geClass(y))
	      y@x <- c(as.double(x), y@x)
	      y@Dim[2] <- nc + 1:1
	      if(is.character(dn <- y@Dimnames[[2]]))
		  y@Dimnames[[2]] <- c("", dn)
	      y
	  })


setMethod("cbind2", signature(x = "denseMatrix", y = "matrix"),
	  function(x, y) callGeneric(x, as(y, geClass(y))))
setMethod("cbind2", signature(x = "matrix", y = "denseMatrix"),
	  function(x, y) callGeneric(as(x, geClass(x)), y))

setMethod("cbind2", signature(x = "denseMatrix", y = "denseMatrix"),
	  function(x, y) {
	      nr <- rowCheck(x,y)
	      ncx <- x@Dim[2]
	      ncy <- y@Dim[2]
	      ## beware of (packed) triangular, symmetric, ...
	      hasDN <- !is.null(dnx <- dimnames(x)) |
	      !is.null(dny <- dimnames(y))
	      x <- as(x, geClass(x))
	      y <- as(y, geClass(y))
	      x@x <- c(x@x, y@x)
	      x@Dim[2] <- ncx + ncy
	      if(hasDN) {
		  ## R and S+ are different in which names they take
		  ## if they differ -- but there's no warning in any case
		  rn <- if(!is.null(dnx[[1]])) dnx[[1]] else dny[[1]]
		  cx <- dnx[[2]] ; cy <- dny[[2]]
		  cn <- if(is.null(cx) && is.null(cy)) NULL
		  else c(if(!is.null(cx)) cx else rep.int("", ncx),
			 if(!is.null(cy)) cy else rep.int("", ncy))
		  x@Dimnames <- list(rn, cn)
	      }
	      x
	  })

### rbind2 -- analogous to cbind2 --- more to do for @x though:

setMethod("rbind2", signature(x = "denseMatrix", y = "numeric"),
	  function(x, y) {
	      if(is.character(dn <- x@Dimnames[[1]])) dn <- c(dn, "")
	      y <- rbind2(as(x,"matrix"), y)
	      new(paste(.M.kind(y), "geMatrix", sep=''), x = c(y),
                  Dim = x@Dim + 1:0, Dimnames = list(dn, x@Dimnames[[2]]))
	  })
## the same, (x,y) <-> (y,x):
setMethod("rbind2", signature(x = "numeric", y = "denseMatrix"),
	  function(x, y) {
	      if(is.character(dn <- y@Dimnames[[1]])) dn <- c("", dn)
	      x <- rbind2(x, as(y,"matrix"))
	      new(paste(.M.kind(x), "geMatrix", sep=''), x = c(x),
                  Dim = y@Dim + 1:0, Dimnames = list(dn, y@Dimnames[[2]]))
	  })

setMethod("rbind2", signature(x = "denseMatrix", y = "matrix"),
	  function(x, y) callGeneric(x, as(y, geClass(y))))
setMethod("rbind2", signature(x = "matrix", y = "denseMatrix"),
	  function(x, y) callGeneric(as(x, geClass(x)), y))

setMethod("rbind2", signature(x = "denseMatrix", y = "denseMatrix"),
	  function(x, y) {
	      nc <- colCheck(x,y)
	      nrx <- x@Dim[1]
	      nry <- y@Dim[1]
	      dn <-
		  if(!is.null(dnx <- dimnames(x)) |
		     !is.null(dny <- dimnames(y))) {
		      ## R and S+ are different in which names they take
		      ## if they differ -- but there's no warning in any case
		      list(if(is.null(rx <- dnx[[1]]) && is.null(ry <- dny[[1]]))
			   NULL else
			   c(if(!is.null(rx)) rx else rep.int("", nrx),
			     if(!is.null(ry)) ry else rep.int("", nry)),
			   if(!is.null(dnx[[2]])) dnx[[2]] else dny[[2]])

		  } else list(NULL, NULL)
	      ## beware of (packed) triangular, symmetric, -> "cheap" (FIXME):
              x <- rbind2(as(x,"matrix"), as(y,"matrix"))
	      new(paste(.M.kind(x), "geMatrix", sep=''), x = c(x),
                  Dim = c(nrx + nry, nc), Dimnames = dn)
	  })
