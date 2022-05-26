## Methods for virtual class "denseMatrix" of dense matrices
##
## Some of these are merely fallback methods, which are "cheap" to implement
## but potentially far from efficient. More efficient methods for subclasses
## will overwrite these.

## ~~~~ COERCIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## README: Many of these coercions are for _subclasses_ of denseMatrix.
##         They have been centralized here quite on purpose, for easier
##         maintenance, to prevent infelicities among groups of similar
##         methods, and to avoid accidental gaps in implementation.

..dense2ge  <- function(from) .Call(R_dense_as_geMatrix, from, ".")
..dense2dge <- function(from) .Call(R_dense_as_geMatrix, from, "d")
..dense2lge <- function(from) .Call(R_dense_as_geMatrix, from, "l")
..dense2nge <- function(from) .Call(R_dense_as_geMatrix, from, "n")

..m2ge  <- function(from) .Call(R_matrix_as_geMatrix, from, ".")
..m2dge <- function(from) .Call(R_matrix_as_geMatrix, from, "d")
..m2lge <- function(from) .Call(R_matrix_as_geMatrix, from, "l")
..m2nge <- function(from) .Call(R_matrix_as_geMatrix, from, "n")

..dense2dkind <- function(from) .Call(R_dense_as_kind, from, "d")
..dense2lkind <- function(from) .Call(R_dense_as_kind, from, "l")
..dense2nkind <- function(from) .Call(R_dense_as_kind, from, "n")

..pack   <- function(from) pack(from)
..packT  <- function(from) pack(from, symmetric = FALSE)
..packS  <- function(from) pack(from, symmetric = TRUE)
..unpack <- function(from) unpack(from)

..pack3    <- function(from) .Call(unpackedMatrix_pack, from, FALSE, NA, NA)
..unpack3  <- function(from) .Call(packedMatrix_unpack, from, FALSE)

## To denseMatrix ..........................................

setAs("ANY", "denseMatrix",
      function(from) Matrix(from, sparse = FALSE, doDiag = FALSE))

## Keep synchronized with Matrix() in ./Matrix.R, where diagonal "matrix"
## (which are both symmetric and triangular) are coerced to "symmetricMatrix",
## _not_ "triangularMatrix"
setAs( "matrix", "denseMatrix", ..m2dense)
setAs("numLike", "denseMatrix", ..m2ge)

## To base matrix ..........................................

setAs("denseMatrix", "matrix", .dense2m)
.dsy2mat <- function(from, ## MJ: for backwards compatibility, until deprecated
                     keep.dimnames = TRUE) {
    to <- .dense2m(from)
    if (!keep.dimnames)
        dimnames(to) <- NULL
    to
}

## Faster:
setAs("dgeMatrix", "matrix", .ge2m)
setAs("lgeMatrix", "matrix", .ge2m)
setAs("ngeMatrix", "matrix", .ge2m)

## To packed ...............................................

setAs("denseMatrix", "packedMatrix", ..pack)
setAs(     "matrix", "packedMatrix", ..pack)
.dsy2dsp <- ..pack # MJ: for backwards compatibility, until deprecated

## To unpacked .............................................

setAs("denseMatrix", "unpackedMatrix", ..unpack)
setAs(     "matrix", "unpackedMatrix", ..m2ge)
setAs(    "numLike", "unpackedMatrix", ..m2ge)

## To general ..............................................

setAs("denseMatrix", "generalMatrix", ..dense2ge)
setAs(     "matrix", "generalMatrix", ..m2ge)
setAs(    "numLike", "generalMatrix", ..m2ge)

## MJ: no longer
if(FALSE) {
setAs("denseMatrix", "generalMatrix", as_geSimpl)
} ## MJ

## To symmetric ............................................

## setAs("denseMatrix", "symmetricMatrix", .) # inherited from Matrix
## setAs(     "matrix", "symmetricMatrix", .) # in ./symmetricMatrix.R

## To triangular ...........................................

## setAs("denseMatrix", "triangularMatrix", .) # inherited from Matrix
## setAs(     "matrix", "triangularMatrix", .) # in ./triangularMatrix.R

## To "kind" ...............................................

setAs("denseMatrix", "dMatrix", ..dense2dkind)
setAs("denseMatrix", "lMatrix", ..dense2lkind)
setAs("denseMatrix", "nMatrix", ..dense2nkind)
setAs("denseMatrix", "ddenseMatrix", ..dense2dkind)
setAs("denseMatrix", "ldenseMatrix", ..dense2lkind)
setAs("denseMatrix", "ndenseMatrix", ..dense2nkind)

setAs("matrix", "ddenseMatrix",
      function(from) {
          if (!is.double(from))
              storage.mode(from) <- "double"
          .m2dense(from)
      })
setAs("matrix", "ldenseMatrix",
      function(from) {
          if (!is.logical(from))
              storage.mode(from) <- "logical"
          .m2dense(from)
      })
setAs("matrix", "ndenseMatrix",
      function(from) as(as(from, "ldenseMatrix"), "nMatrix"))

setAs("numLike", "ddenseMatrix", ..m2dge)
setAs("numLike", "ldenseMatrix", ..m2lge)
setAs("numLike", "ndenseMatrix", ..m2nge)

## FIXME: Conceivably, we could also have
##
## setAs( "matrix", "[dln]Matrix", .)
## setAs("numLike", "[dln]Matrix", .)
##
## but those ought to return dense _or_ sparse.
## Note that there is already:
##
## setAs("matrix", "lMatrix", .)  [ in ./lMatrix.R   ]
## setAs("matrix", "nMatrix", .)  [ in ./ngTMatrix.R ]
##
## the former returning Matrix(`storage.mode<-`(from, "logical"))
## and the latter always returning an ngTMatrix (i.e.,
## the current implementation is neither complete nor consistent) ...

## More granular coercions .................................

.kinds <- c("d", "l", "n")
for (.kind in .kinds) {
    ## General to non-general, preserving kind
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "trMatrix"), ..M2tri)
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "syMatrix"), ..M2symm)
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "tpMatrix"), ..packT)
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "spMatrix"), ..packS)

    ## Non-general to general, preserving kind
    .f <- get(paste0("..dense2", .kind, "ge"),
              mode = "function", inherits = FALSE)
    for (.xx in c("tr", "sy", "tp", "sp"))
        setAs(paste0(.kind, .xx,   "Matrix"),
              paste0(.kind,      "geMatrix"), .f)

    ## Unpacked to packed, preserving kind and structure
    setAs(paste0(.kind, "trMatrix"), paste0(.kind, "tpMatrix"), ..pack3)
    setAs(paste0(.kind, "syMatrix"), paste0(.kind, "spMatrix"), ..pack3)

    ## Packed to unpacked, preserving kind and structure
    setAs(paste0(.kind, "tpMatrix"), paste0(.kind, "trMatrix"), ..unpack3)
    setAs(paste0(.kind, "spMatrix"), paste0(.kind, "syMatrix"), ..unpack3)

    ## Triangular to symmetric, preserving kind and storage
    setAs(paste0(.kind, "trMatrix"), paste0(.kind, "syMatrix"), ..M2symm)
    setAs(paste0(.kind, "tpMatrix"), paste0(.kind, "spMatrix"), ..M2symm)

    ## Symmetric to triangular, preserving kind and storage
    setAs(paste0(.kind, "syMatrix"), paste0(.kind, "trMatrix"), ..M2tri)
    setAs(paste0(.kind, "spMatrix"), paste0(.kind, "tpMatrix"), ..M2tri)

    ## Other kinds to this kind, preserving structure and storage
    .g <- get(paste0("..dense2", .kind, "kind"),
              mode = "function", inherits = FALSE)
    .otherkinds <- .kinds[.kinds != .kind]
    for (.otherkind in .otherkinds)
        for (.xx in c("ge", "tr", "sy", "tp", "sp"))
            setAs(paste0(.otherkind, .xx, "Matrix"),
                  paste0(     .kind, .xx, "Matrix"), .g)

    ## Rely on coercions to virtual classes when changing two or more
    ## of kind, structure, and storage at once:

    ##  NO: as(<dsyMatrix>, "dtpMatrix")
    ## YES: as(as(<dsyMatrix>, "triangularMatrix"), "packedMatrix")

    ##  NO: as(<lgeMatrix>, "nspMatrix")
    ## YES: as(as(<lgeMatrix>, "lspMatrix"), "nMatrix")
}
rm(.kind, .kinds, .otherkind, .otherkinds, .f, .g, .xx)

## From base matrix ........................................

## m->ge
setAs("matrix", "dgeMatrix", ..m2dge)
setAs("matrix", "lgeMatrix", ..m2lge)
setAs("matrix", "ngeMatrix", ..m2nge)

## m->tr
setAs("matrix", "dtrMatrix",
      function(from) {
          if(!is.double(from))
              storage.mode(from) <- "double"
          .M2tri(from)
      })
setAs("matrix", "ltrMatrix",
      function(from) {
          if(!is.logical(from))
              storage.mode(from) <- "logical"
          .M2tri(from)
      })
setAs("matrix", "ntrMatrix",
      function(from) as(as(from, "ltrMatrix"), "nMatrix"))

## m->sy
setAs("matrix", "dsyMatrix",
      function(from) {
          if(!is.double(from))
              storage.mode(from) <- "double"
          .M2symm(from)
      })
setAs("matrix", "lsyMatrix",
      function(from) {
          if(!is.logical(from))
              storage.mode(from) <- "logical"
          .M2symm(from)
      })
setAs("matrix", "nsyMatrix",
      function(from) {
          if(anyNA(from))
              from <- is.na(from) | from
          as(as(from, "lsyMatrix"), "nMatrix")
      })

## m->tp
setAs("matrix", "dtpMatrix",
      function(from) {
          if(!is.double(from))
              storage.mode(from) <- "double"
          pack(from, symmetric = FALSE)
      })
setAs("matrix", "ltpMatrix",
      function(from) {
          if(!is.logical(from))
              storage.mode(from) <- "logical"
          pack(from, symmetric = FALSE)
      })
setAs("matrix", "ntpMatrix",
      function(from) as(as(from, "ltpMatrix"), "nMatrix"))

## m->sp
setAs("matrix", "dspMatrix",
      function(from) {
          if(!is.double(from))
              storage.mode(from) <- "double"
          pack(from, symmetric = TRUE)
      })
setAs("matrix", "lspMatrix",
      function(from) {
          if(!is.logical(from))
              storage.mode(from) <- "logical"
          pack(from, symmetric = TRUE)
      })
setAs("matrix", "nspMatrix",
      function(from) {
          if(anyNA(from))
              from <- is.na(from) | from
          as(as(from, "lspMatrix"), "nMatrix")
      })

## From base vector ........................................

setAs("numLike", "dgeMatrix", ..m2dge)
setAs("numLike", "lgeMatrix", ..m2lge)
setAs("numLike", "ngeMatrix", ..m2nge)

rm(..dense2ge, ..dense2dge, ..dense2lge, ..dense2nge,
   ..m2ge, ..m2dge, ..m2lge, ..m2nge,
   ..dense2dkind, ..dense2lkind, ..dense2nkind,
   ..pack, ..packT, ..packS, ..unpack, ..pack3, ..unpack3)

## To sparse ...............................................

## dense to sparse:
## : if we do this, do it "right", i.e. preserve symmetric/triangular!
## setAs("denseMatrix", "dsparseMatrix",
## ## MM thought that  as() will take the ``closest'' match; but that fails!
## ##      function(from) as(as(from, "dgeMatrix"), "dsparseMatrix"))
##       function(from) as(as(from, "dgeMatrix"), "dgCMatrix"))

.dense2C <- function(from, kind = NA, uplo = "U", symDimnames = FALSE) {
    useK <- is.character(kind) && length(kind) == 1 &&
        kind %in% c("gen", "sym", "tri")
    if(!useK) {
        cl <- class(from)
        cld <- getClassDef(cl) ## get it once (speedup)
    }
    r <- .Call(dense_to_Csparse, from)# goes via "generalMatrix"
    ## FIXME: for symmetric / triangular matrices, this is a waste, notably if packed
    if (useK && kind == "gen"  ||  !useK && extends(cld, "generalMatrix"))
	r
    else if(useK && kind == "sym" || !useK && extends(cld, "symmetricMatrix"))
	forceCspSymmetric(r, uplo, isTri = FALSE, symDimnames=symDimnames)
    else if(!useK && extends(cld, "diagonalMatrix"))
	stop("diagonalMatrix in .dense2C() -- should never happen, please report!")
    else { ## we have "triangular" :
        if(useK) {
            cl <- class(from)
            cld <- getClassDef(cl) ## get it once (speedup)
        }
	if	(extends(cld,"dMatrix")) as(r, "dtCMatrix")
        else if (extends(cld,"lMatrix")) as(r, "ltCMatrix")
        else if (extends(cld,"nMatrix")) as(r, "ntCMatrix")
        else if (extends(cld,"zMatrix")) as(r, "ztCMatrix")
	else stop(gettextf("undefined method for class %s", dQuote(cl)), domain=NA)
    }
}

setAs("denseMatrix", "CsparseMatrix", function(from) .dense2C(from))

## This sometimes fails (eg. for "lsyMatrix"), and we really want to
## use the generic ``go via Csparse'' (top of ./sparseMatrix.R) instead
## setAs("denseMatrix",  "sparseMatrix",
##       function(from) {
## 	  cl <- class(from)
## 	  cld <- getClassDef(cl)
## 	  if (extends(cld, "generalMatrix"))
## 	      .Call(dense_to_Csparse, from)
## 	  else ## i.e. triangular | symmetric
## 	      as_Csparse(from, cld)
##       })

setAs("denseMatrix", "TsparseMatrix",
      function(from) as(.dense2C(from), "TsparseMatrix"))


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dense.band <- function(x, k1, k2, ...) .Call(R_dense_band, x, k1, k2)
.dense.triu <- function(x, k = 0,  ...) .Call(R_dense_band, x, k, NULL)
.dense.tril <- function(x, k = 0,  ...) .Call(R_dense_band, x, NULL, k)
for (.cl in c("denseMatrix", "matrix")) {
    setMethod("band", signature(x = .cl), .dense.band)
    setMethod("triu", signature(x = .cl), .dense.triu)
    setMethod("tril", signature(x = .cl), .dense.tril)
}
rm(.dense.band, .dense.triu, .dense.tril, .cl)

setMethod("show", signature(object = "denseMatrix"),
          function(object) prMatrix(object))
##- ## FIXME: The following is only for the "dMatrix" objects that are not
##- ##	      "dense" nor "sparse" -- i.e. "packed" ones :
##- ## But these could be printed better -- "." for structural zeros.
##- setMethod("show", signature(object = "dMatrix"), prMatrix)
##- ## and improve this as well:
##- setMethod("show", signature(object = "pMatrix"), prMatrix)
##- ## this should now be superfluous [keep for safety for the moment]:

setMethod("dim<-", signature(x = "denseMatrix", value = "ANY"),
	  function(x, value) {
	      if(!is.numeric(value) || length(value) != 2)
		  stop("dim(.) value must be numeric of length 2")
	      if(prod(dim(x)) != prod(value <- as.integer(value)))
		  stop("dimensions don't match the number of cells")
	      clx <- as.character(MatrixClass(class(x))) # as.*(): drop attr
	      if(substring(clx,2) == "geMatrix") {
		  x@Dim <- value
		  if(length(x@factors) > 0)
		      x@factors <- list()
		  x
	      } else { ## other "denseMatrix"
		  x <- as_geSimpl2(x, clx)
		  dim(x) <- value
                  x
	      }
          })



## Using "index" for indices should allow
## integer (numeric), logical, or character (names!) indices :

## use geClass() when 'i' or 'j' are missing:
## since  symmetric, triangular, .. will not be preserved anyway:
setMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop) {
	      if((na <- nargs()) == 3)
		  r <- as(x, "matrix")[i, drop=drop]
	      else if(na == 4)
		  r <- as(x, "matrix")[i, , drop=drop]
	      else stop(gettextf("invalid nargs()= %d", na), domain=NA)
	      if(is.null(dim(r))) r else as(r, geClass(x))
	  })

setMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
			 drop = "logical"),
	  function (x, i, j, ..., drop) {
	      r <- as(x, "matrix")[, j, drop=drop]
	      if(is.null(dim(r))) r else as(r, geClass(x))
	  })

setMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
			 drop = "logical"),
	  function (x, i, j, ..., drop) {
	      r <- callGeneric(x = as(x, "matrix"), i=i, j=j, drop=drop)
	      if(is.null(dim(r)))
		  r
	      else {
		  cld <- getClassDef(cl <- class(x))
		  if(extends(cld, "symmetricMatrix") &&
		     length(i) == length(j) && isTRUE(all(i == j)))
                      ## keep original symmetric class (but not "dpo")
                      as(r, class2(cl, .M.kindC(cld)))

		  else as_smartClass(r, cl)
	      }
	  })

.dense.sub.i.2col <- function(x, i, j, ..., drop) {
    r <- as(x, "matrix")[ i ]
    if(is.null(dim(r))) r else as(r, geClass(x))
}
setMethod("[", signature(x = "denseMatrix", i = "matrix", j = "missing"),#drop="ANY"
	  .dense.sub.i.2col)
setMethod("[", signature(x = "denseMatrix", i = "matrix", j = "missing", drop="missing"),
	  .dense.sub.i.2col)


## Now the "[<-" ones --- see also those in ./Matrix.R
## It's recommended to use setReplaceMethod() rather than setMethod("[<-",.)
## even though the former is currently just a wrapper for the latter

## x[] <- value :
setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "missing",
				value = "ANY"),## double/logical/...
	  function (x, value) {
	      x <- as(x, "generalMatrix")
	      x@x[] <- value
	      validObject(x)# check if type and lengths above match
	      x
	  })

## FIXME: 1) These are far from efficient
## -----
setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "missing",
				value = "replValue"),
		 function (x, i, j, ..., value) {
		     r <- as(x, "matrix")
## 		     message("`[<-` with nargs()= ",nargs())
		     if((na <- nargs()) == 3)
			 r[i] <- value
		     else if(na == 4)
			 r[i, ] <- value
		     else stop(gettextf("invalid nargs()= %d", na), domain=NA)
		     as(r, geClass(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "missing", j = "index",
				value = "replValue"),
		 function (x, i, j, ..., value) {
		     r <- as(x, "matrix")
		     r[, j] <- value
		     as(r, geClass(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "index", j = "index",
				value = "replValue"),
		 function (x, i, j, ..., value) {
		     r <- as(x, "matrix")
		     r[i, j] <- value
		     as_smartClass(r, class(x)) ## was as(r, class(x))
		 })

setReplaceMethod("[", signature(x = "denseMatrix", i = "matrix",  # 2-col.matrix
				j = "missing", value = "replValue"),
		 function(x, i, j, ..., value) {
		     r <- as(x, "matrix")
		     r[ i ] <- value
		     as(r, geClass(x))
		 })

## MJ: no longer needed ... replacement in ./(un)?packedMatrix.R
if(FALSE) {
setMethod("isSymmetric", signature(object = "denseMatrix"),
	  function(object, tol = 100*.Machine$double.eps, tol1 = 8*tol, ...) {
	      ## pretest: is it square?
	      d <- dim(object)
              if((n <- d[1L]) != d[2L]) return(FALSE)
	      if(n <= 1L) return(TRUE)
	      ## else: square (n x n) matrix, n >= 2 :
              is.z <- is(object, "zMatrix")
	      ## initial tests, fast for large non-symmetric:
	      if(length(tol1)) {
		  ## initial pre-tests, fast for large non-symmetric:
		  Cj <- if(is.z) Conj else identity
		  for(i in unique(c(1L, 2L, n-1L, n)))
		      if(is.character(all.equal(object[i, ], Cj(object[, i]),
						tolerance = tol1, ...))) return(FALSE)
	      }
	      ## else slower test
	      if (is(object,"dMatrix"))
		  isTRUE(all.equal(as(  object,  "dgeMatrix"),
				   as(t(object), "dgeMatrix"), tolerance = tol, ...))
	      else if (is(object, "nMatrix"))
		  identical(as(  object,  "ngeMatrix"),
			    as(t(object), "ngeMatrix"))
	      else if (is(object, "lMatrix"))# not possible currently
		  ## test for exact equality; FIXME(?): identical() too strict?
		  identical(as(  object,  "lgeMatrix"),
			    as(t(object), "lgeMatrix"))
	      else if (is.z) ## will error out here
		  isTRUE(all.equal(as(       object,   "zgeMatrix"),
				   as(Conj(t(object)), "zgeMatrix"),
				   tolerance = tol, ...))
	      else if (is(object, "iMatrix")) ## will error out here
		  identical(as(object, "igeMatrix"),
			    as(t(object), "igeMatrix"))
	  })
setMethod("isTriangular", signature(object = "denseMatrix"), isTriMat)
setMethod("isDiagonal", signature(object = "denseMatrix"), .is.diagonal)
} ## MJ

setMethod("rcond", signature(x = "denseMatrix", norm = "character"),
	  function(x, norm, ...)
	  rcond(as(as(x, "dMatrix"), "dgeMatrix"), norm=norm, ...))

setMethod("symmpart", signature(x = "denseMatrix"),
	  function(x) symmpart(as(x, "dMatrix")))
setMethod("skewpart", signature(x = "denseMatrix"),
	  function(x) skewpart(as(x, "dMatrix")))

setMethod("is.na", signature(x = "denseMatrix"),
	  function(x) {
	      if(any((inax <- is.na(x@x)))) {
		  r <- as(x, "lMatrix")#-> logical x-slot
		  r@x <- inax
		  as(r, "nMatrix")
	      } else {
		  d <- x@Dim
		  new("ngCMatrix", Dim = d, Dimnames = dimnames(x),
		      i = integer(0), p = rep.int(0L, d[2]+1L))
	      }
	  })

if(.Matrix.avoiding.as.matrix) {
setMethod("qr", signature(x = "ddenseMatrix"),
	  function(x, ...) qr.default(.dense2m(x), ...))
setMethod("qr", signature(x = "denseMatrix"),
	  function(x, ...) qr(as(x, "ddenseMatrix"), ...))
}

