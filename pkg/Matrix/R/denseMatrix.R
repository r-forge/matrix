## Methods for virtual class "denseMatrix" of dense matrices
## ... and many for base matrices and base vectors, too

## ~~~~ COERCIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## README: Many of these coercions are for _subclasses_ of denseMatrix.
##         They have been centralized here quite on purpose, for easier
##         maintenance, to prevent infelicities among groups of similar
##         methods, and to avoid accidental gaps in implementation.

..dense2g  <- function(from) .Call(R_dense_as_general, from, ".")
..dense2dg <- function(from) .Call(R_dense_as_general, from, "d")
..dense2lg <- function(from) .Call(R_dense_as_general, from, "l")
..dense2ng <- function(from) .Call(R_dense_as_general, from, "n")

..dense2ddense <- function(from) .Call(R_dense_as_kind, from, "d")
..dense2ldense <- function(from) .Call(R_dense_as_kind, from, "l")
..dense2ndense <- function(from) .Call(R_dense_as_kind, from, "n")

..dense2dsparse <- function(from)
    .Call(R_dense_as_sparse, from, "d.C", NULL, NULL)
..dense2lsparse <- function(from)
    .Call(R_dense_as_sparse, from, "l.C", NULL, NULL)
..dense2nsparse <- function(from)
    .Call(R_dense_as_sparse, from, "n.C", NULL, NULL)

..dense2Csparse <- function(from)
    .Call(R_dense_as_sparse, from, "..C", NULL, NULL)
..dense2Rsparse <- function(from)
    .Call(R_dense_as_sparse, from, "..R", NULL, NULL)
..dense2Tsparse <- function(from)
    .Call(R_dense_as_sparse, from, "..T", NULL, NULL)

..m2ge  <- function(from) .Call(R_matrix_as_geMatrix, from, ".")
..m2dge <- function(from) .Call(R_matrix_as_geMatrix, from, "d")
..m2lge <- function(from) .Call(R_matrix_as_geMatrix, from, "l")
..m2nge <- function(from) .Call(R_matrix_as_geMatrix, from, "n")

..m2dense  <- function(from) .m2dense(from, ".")
..m2ddense <- function(from) .m2dense(from, "d")
..m2ldense <- function(from) .m2dense(from, "l")
..m2ndense <- function(from) .m2dense(from, "n")

..m2dsparse <- function(from) .m2sparse(from, "d", "C")
..m2lsparse <- function(from) .m2sparse(from, "l", "C")
..m2nsparse <- function(from) .m2sparse(from, "n", "C")

..m2Csparse <- function(from) .m2sparse(from, ".", "C")
..m2Rsparse <- function(from) .m2sparse(from, ".", "R")
..m2Tsparse <- function(from) .m2sparse(from, ".", "T")

..pack    <- function(from) pack(from)
..packT   <- function(from) pack(from, symmetric = FALSE)
..packS   <- function(from) pack(from, symmetric = TRUE)
..unpack  <- function(from) unpack(from)
## returning .(ge|tr|sy|tp|sp)Matrix, never a subclass:
..pack3   <- function(from) .Call(unpackedMatrix_pack, from, FALSE, NA, NA)
..unpack3 <- function(from) .Call(packedMatrix_unpack, from, FALSE)

## To denseMatrix ..........................................

setAs("ANY", "denseMatrix",
      function(from) Matrix(from, sparse = FALSE, doDiag = FALSE))

setAs( "matrix", "denseMatrix", ..m2dense)
setAs("numLike", "denseMatrix", ..m2ge)

## To sparse ...............................................

setAs("denseMatrix",  "sparseMatrix", ..dense2Csparse)
setAs("denseMatrix", "CsparseMatrix", ..dense2Csparse)
setAs("denseMatrix", "RsparseMatrix", ..dense2Rsparse)
setAs("denseMatrix", "TsparseMatrix", ..dense2Tsparse)

setAs("matrix",  "sparseMatrix", ..m2Csparse)
setAs("matrix", "CsparseMatrix", ..m2Csparse)
setAs("matrix", "RsparseMatrix", ..m2Rsparse)
setAs("matrix", "TsparseMatrix", ..m2Tsparse)

setAs("numLike",  "sparseMatrix", ..dense2Csparse)
setAs("numLike", "CsparseMatrix", ..dense2Csparse)
setAs("numLike", "RsparseMatrix", ..dense2Rsparse)
setAs("numLike", "TsparseMatrix", ..dense2Tsparse)

## To base matrix, base vector .............................

setAs("denseMatrix", "matrix", .dense2m)
setAs("denseMatrix", "vector", .dense2v)

setMethod("as.vector", signature(x = "denseMatrix"),
          function(x, mode) as.vector(.dense2v(x), mode))

setMethod("as.numeric", signature(x = "denseMatrix"),
          function(x, ...) as.double(.dense2v(x)))
setMethod("as.numeric", signature(x = "ddenseMatrix"),
          function(x, ...) .dense2v(x))

setMethod("as.logical", signature(x = "denseMatrix"),
          function(x, ...) as.logical(.dense2v(x)))
setMethod("as.logical", signature(x = "ldenseMatrix"),
          function(x, ...) .dense2v(x))
setMethod("as.logical", signature(x = "ndenseMatrix"),
          function(x, ...) .dense2v(x))

## Faster:
for (.from in paste0(c("d", "l", "n"), "geMatrix")) {
    setAs(.from, "matrix", .ge2m)
    setAs(.from, "vector", .ge2v)

    setMethod("as.vector", signature(x = .from),
              function(x, mode) as.vector(x@x, mode))
}
rm(.from)

## To "kind" ...............................................

setAs("denseMatrix", "dMatrix", ..dense2ddense)
setAs("denseMatrix", "lMatrix", ..dense2ldense)
setAs("denseMatrix", "nMatrix", ..dense2ndense)

## With a dense _or_ sparse result:
for (.from in c("matrix", "numLike")) {
    setAs(.from, "dMatrix",
          function(from) {
              storage.mode(from) <- "double"
              Matrix(from)
          })
    setAs(.from, "lMatrix",
          function(from) {
              storage.mode(from) <- "logical"
              Matrix(from)
          })
    setAs(.from, "nMatrix",
          function(from) {
              storage.mode(from) <- "logical"
              if(anyNA(from))
                  from[is.na(from)] <- TRUE # needed before symmetry test!
              as(Matrix(from), "nMatrix")
          })
}
rm(.from)

setAs("denseMatrix", "ddenseMatrix", ..dense2ddense)
setAs("denseMatrix", "ldenseMatrix", ..dense2ldense)
setAs("denseMatrix", "ndenseMatrix", ..dense2ndense)

setAs("matrix", "ddenseMatrix", ..m2ddense)
setAs("matrix", "ldenseMatrix", ..m2ldense)
setAs("matrix", "ndenseMatrix", ..m2ndense)

setAs("numLike", "ddenseMatrix", ..m2dge)
setAs("numLike", "ldenseMatrix", ..m2lge)
setAs("numLike", "ndenseMatrix", ..m2nge)

setAs("denseMatrix", "dsparseMatrix", ..dense2dsparse)
setAs("denseMatrix", "lsparseMatrix", ..dense2lsparse)
setAs("denseMatrix", "nsparseMatrix", ..dense2nsparse)

setAs("matrix", "dsparseMatrix", ..m2dsparse)
setAs("matrix", "lsparseMatrix", ..m2lsparse)
setAs("matrix", "nsparseMatrix", ..m2nsparse)

setAs("numLike", "dsparseMatrix", ..dense2dsparse)
setAs("numLike", "lsparseMatrix", ..dense2lsparse)
setAs("numLike", "nsparseMatrix", ..dense2nsparse)

## To general ..............................................

setAs("denseMatrix", "generalMatrix", ..dense2g)
setAs(     "matrix", "generalMatrix", ..m2ge)
setAs(    "numLike", "generalMatrix", ..m2ge)

## To symmetric ............................................

## setAs("denseMatrix", "symmetricMatrix", .) # inherited from Matrix
## setAs(     "matrix", "symmetricMatrix", .) # in ./symmetricMatrix.R

## To triangular ...........................................

## setAs("denseMatrix", "triangularMatrix", .) # inherited from Matrix
## setAs(     "matrix", "triangularMatrix", .) # in ./triangularMatrix.R

## To unpacked .............................................

setAs("denseMatrix", "unpackedMatrix", ..unpack)
setAs(     "matrix", "unpackedMatrix", ..m2dense)
setAs(    "numLike", "unpackedMatrix", ..m2ge)

## To packed ...............................................

setAs("denseMatrix", "packedMatrix", ..pack)
setAs(     "matrix", "packedMatrix", ..pack)

## More granular coercions .................................

## DEPRECATED IN 1.4-2; see ./zzz.R
if(FALSE) {
.kinds <- c("d", "l", "n")
for (.kind in .kinds) {
    ## General to non-general, preserving kind
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "trMatrix"), ..M2tri)
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "syMatrix"), ..M2symm)
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "tpMatrix"), ..packT)
    setAs(paste0(.kind, "geMatrix"), paste0(.kind, "spMatrix"), ..packS)

    ## Non-general to general, preserving kind
    for (.xy in c("tr", "sy", "tp", "sp"))
        setAs(paste0(.kind, .xy,   "Matrix"),
              paste0(.kind,      "geMatrix"),
              get(sprintf("..dense2%sg", .kind),
                  mode = "function", inherits = FALSE))

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
    .otherkinds <- .kinds[.kinds != .kind]
    for (.otherkind in .otherkinds)
        for (.xy in c("ge", "tr", "sy", "tp", "sp"))
            setAs(paste0(.otherkind, .xy, "Matrix"),
                  paste0(     .kind, .xy, "Matrix"),
                  get(sprintf("..dense2%sdense", .kind),
                      mode = "function", inherits = FALSE))

    ## Dense to sparse, preserving kind and structure
    for (.xy in c("ge", "tr", "sy", "tp", "sp"))
        for (.repr in c("C", "R", "T"))
            setAs(paste0(.kind,            .xy,                 "Matrix"),
                  paste0(.kind, `substr<-`(.xy, 2L, 2L, .repr), "Matrix"),
                  get(sprintf("..dense2%ssparse", .repr),
                      mode = "function", inherits = FALSE))
}
rm(.kind, .kinds, .otherkind, .otherkinds, .xy, .repr)

## For whatever reason, we also have these granular ones in Matrix 1.4-1:
setAs("dgeMatrix", "dsTMatrix",
      function(from) .dense2sparse(.M2symm(from), "..T", NULL, NULL))
} ## DEPRECATED IN 1.4-2; see ./zzz.R

## From base matrix, base vector ...........................

## DEPRECATED IN 1.4-2; see ./zzz.R
if(FALSE) {
## m|v->ge
for (.from in c("matrix", "numLike")) {
    setAs(.from, "dgeMatrix", ..m2dge)
    setAs(.from, "lgeMatrix", ..m2lge)
    setAs(.from, "ngeMatrix", ..m2nge)
}
rm(.from)

## m->tr
setAs("matrix", "dtrMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .M2tri(from)
      })
setAs("matrix", "ltrMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          .M2tri(from)
      })
setAs("matrix", "ntrMatrix",
      function(from) as(as(from, "ltrMatrix"), "nMatrix"))

## m->sy
setAs("matrix", "dsyMatrix",
      function(from) {
          storage.mode(from) <- "double"
          .M2symm(from)
      })
setAs("matrix", "lsyMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          .M2symm(from)
      })
setAs("matrix", "nsyMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          if(anyNA(from))
              from[is.na(from)] <- TRUE # needed before symmetry test!
          as(as(from, "lsyMatrix"), "nMatrix")
      })

## m->tp
setAs("matrix", "dtpMatrix",
      function(from) {
          storage.mode(from) <- "double"
          pack(from, symmetric = FALSE)
      })
setAs("matrix", "ltpMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          pack(from, symmetric = FALSE)
      })
setAs("matrix", "ntpMatrix",
      function(from) as(as(from, "ltpMatrix"), "nMatrix"))

## m->sp
setAs("matrix", "dspMatrix",
      function(from) {
          storage.mode(from) <- "double"
          pack(from, symmetric = TRUE)
      })
setAs("matrix", "lspMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          pack(from, symmetric = TRUE)
      })
setAs("matrix", "nspMatrix",
      function(from) {
          storage.mode(from) <- "logical"
          if(anyNA(from))
              from[is.na(from)] <- TRUE # needed before symmetry test!
          as(as(from, "lspMatrix"), "nMatrix")
      })

## m|v->g[CRT]
.to <- ".g.Matrix"
.def <- function(from) {
    .dense2sparse(from, .CODE, NULL, NULL)
}
.b <- body(.def)
for (.kind in c("d", "l", "n")) {
    for (.repr in c("C", "R", "T")) {
        substr(.to, 1L, 1L) <- .kind
        substr(.to, 3L, 3L) <- .repr
        body(.def) <-
            do.call(substitute, list(.b, list(.CODE = substr(.to, 1L, 3L))))
        for (.from in c("matrix", "numLike"))
            setAs(.from, .to, .def)
    }
}

## m->t[CRT]
.to <- ".t.Matrix"
.def <- function(from) {
    if(it <- isTriangular(from))
        .dense2sparse(from, .CODE, attr(it, "kind"), "N")
    else
        stop("matrix is not triangular; consider triu() or tril()")
}
.b <- body(.def)
for (.kind in c("d", "l", "n")) {
    for (.repr in c("C", "R", "T")) {
        substr(.to, 1L, 1L) <- .kind
        substr(.to, 3L, 3L) <- .repr
        body(.def) <-
            do.call(substitute, list(.b, list(.CODE = substr(.to, 1L, 3L))))
        setAs("matrix", .to, .def)
    }
}

## m->s[CRT]
.to <- ".s.Matrix"
.def <- function(from) {
    if(isSymmetric(from))
        .dense2sparse(from, .CODE, "U", NULL)
    else
        stop("matrix is not symmetric; consider forceSymmetric() or symmpart()")
}
.b <- body(.def)
for (.kind in c("d", "l", "n")) {
    for (.repr in c("C", "R", "T")) {
        substr(.to, 1L, 1L) <- .kind
        substr(.to, 3L, 3L) <- .repr
        body(.def) <-
            do.call(substitute, list(.b, list(.CODE = substr(.to, 1L, 3L))))
        setAs("matrix", .to, .def)
    }
}
rm(.from, .to, .def, .b, .kind, .repr)
} ## DEPRECATED IN 1.4-2; see ./zzz.R

## Exported functions, now just aliases or wrappers ........
## (some or all could be made deprecated) ..................

..2dge <- ..dense2dg
.dense2sy <- .M2symm
.dsy2dsp <- ..pack3
.dsy2mat <- function(from, keep.dimnames = TRUE) {
    to <- .dense2m(from)
    if (!keep.dimnames)
        dimnames(to) <- NULL
    to
}
.m2ngCn <- function(from, na.is.not.0 = FALSE) {
    if (!na.is.not.0 && anyNA(from))
        stop("attempt to coerce matrix with NA to \"ngCMatrix\"")
    .dense2sparse(from, "ngC", NULL, NULL)
}
.m2ngTn <- function(from, na.is.not.0 = FALSE) {
    if (!na.is.not.0 && anyNA(from))
        stop("attempt to coerce matrix with NA to \"ngTMatrix\"")
    .dense2sparse(from, "ngT", NULL, NULL)
}
.m2dgC <- function(from) .dense2sparse(from, "dgC", NULL, NULL)
.m2lgC <- function(from) .dense2sparse(from, "lgC", NULL, NULL)
.m2ngC <- function(from) .m2ngCn(from)

rm(..dense2g, ..dense2dg, ..dense2lg, ..dense2ng,
   ..dense2ddense, ..dense2ldense, ..dense2ndense,
   ..dense2dsparse, ..dense2lsparse, ..dense2nsparse,
   ..dense2Csparse, ..dense2Rsparse, ..dense2Tsparse,
   ..m2ge, ..m2dge, ..m2lge, ..m2nge,
   ..m2dense, ..m2ddense, ..m2ldense, ..m2ndense,
   ..m2dsparse, ..m2lsparse, ..m2nsparse,
   ..m2Csparse, ..m2Rsparse, ..m2Tsparse,
   ..pack, ..packT, ..packS, ..unpack,
   ..pack3, ..unpack3)

## MJ: no longer ... replacement above
if(FALSE) {
setAs("denseMatrix", "generalMatrix", as_geSimpl)

## dense to sparse:
## : if we do this, do it "right", i.e. preserve symmetric/triangular!
## setAs("denseMatrix", "dsparseMatrix",
## ## MM thought that  as() will take the ``closest'' match; but that fails!
## ##      function(from) as(as(from, "dgeMatrix"), "dsparseMatrix"))
##       function(from) as(as(from, "dgeMatrix"), "dgCMatrix"))

.dense2C <- function(from, kind = NA, uplo = "U") {
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
	forceSymmetricCsparse(r, uplo)
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
} ## MJ


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

