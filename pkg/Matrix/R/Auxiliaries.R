#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

`%||%` <- function(L, R) if (is.null(L)) R else L # stats:::`%||%`

is0  <- function(x) !(is.na(x) | x)
isN0 <- function(x)   is.na(x) | x
is1  <- function(x)  !is.na(x) & x == 1
isN1 <- function(x)   is.na(x) | x != 1
isT  <- function(x)  !is.na(x) &  x
isNT <- function(x)   is.na(x) | !x

##allFalse <- function(x) !any(x) && !any(is.na(x))## ~= all0, but allFalse(NULL) = TRUE w/warning
##all0 <- function(x) !any(is.na(x)) && all(!x) ## ~= allFalse
allFalse <- function(x) if(is.atomic(x)) .Call(R_all0, x) else !any(x) && !any(is.na(x))
all0     <- function(x) if(is.atomic(x)) .Call(R_all0, x) else all(!x) && !any(is.na(x))

##anyFalse <- function(x) isTRUE(any(!x))		 ## ~= any0
## any0 <- function(x) isTRUE(any(x == 0))	      ## ~= anyFalse
anyFalse <-
any0 <- function(x) if(is.atomic(x)) .Call(R_any0, x) else isTRUE(any(!x))

## These work "identically" for  1 ('==' TRUE)  and 0 ('==' FALSE)
##	(but give a warning for "double"  1 or 0)
## TODO: C versions of these would be faster
allTrue  <- function(x) all(x) && !anyNA(x)

## Note that mode(<integer>) = "numeric" -- as0(), as1() return "double"
## which is good *AS LONG AS* we do not really have i..Matrix integer matrices
as1 <- function(x, mod=mode(x))
    switch(mod, "integer"= 1L, "double"=, "numeric"= 1, "logical"= TRUE,
	   "complex"= 1+0i, stop(gettextf("invalid 'mod': %s", mod), domain = NA))
as0 <- function(x, mod=mode(x))
    switch(mod, "integer"= 0L, "double"=, "numeric"= 0, "logical"= FALSE,
	   "complex"= 0+0i, stop(gettextf("invalid 'mod': %s", mod), domain = NA))

## MJ: no longer used
if(FALSE) {
.bail.out.1 <- function(fun, cl) {
    stop(gettextf(
     'not-yet-implemented method for %s(<%s>).\n ->>  Ask the package authors to implement the missing feature.',
		  fun, cl[1L]), call. = FALSE, domain=NA)
}
} ## MJ

.bail.out.2 <- function(fun, cl1, cl2) {
    stop(gettextf(
     'not-yet-implemented method for %s(<%s>, <%s>).\n ->>  Ask the package authors to implement the missing feature.',
		  fun, cl1[1L], cl2[1L]), call. = FALSE, domain=NA)
}

Matrix.verbose <- function()
    getOption("Matrix.verbose", .MatrixEnv[["verbose"]])

Matrix.warn <- function()
    getOption("Matrix.warn", .MatrixEnv[["warn"]])

## MJ: or maybe a general one?
if(FALSE) {
Matrix.getOption <- function(x, default = .Matrix.Env[[x]])
    getOption(paste0("Matrix.", x), default)
} ## MJ

Matrix.msg <- function(..., .M.level = 1, call. = FALSE, domain = NULL) {
    if(Matrix.verbose() >= .M.level) {
        m <-
            if((w <- Matrix.warn()) < 1)
                function(..., call., domain) message(..., domain = domain)
            else if(w < 2)
                warning
            else stop
        m(..., call. = call., domain = domain)
    }
}

mmultDim <- function(d.a, d.b, type = 1L) {
    ## Return the 'dim' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b    {crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    ## after asserting that ncol(<left operand>) == nrow(<right operand>)
    i.a <- 1L + (type != 2L)
    i.b <- 1L + (type == 3L)
    if(d.a[i.a] != d.b[i.b])
        stop(gettextf("non-conformable matrix dimensions in %s",
                      deparse(sys.call(sys.parent()))),
             call. = FALSE, domain = NA)
    c(d.a[-i.a], d.b[-i.b])
}

mmultDimnames <- function(dn.a, dn.b, type = 1L) {
    ## Return the 'dimnames' of the product indicated by 'type':
    ##     type 1:    a  %*%   b
    ##          2:  t(a) %*%   b    {crossprod}
    ##          3:    a  %*% t(b)  {tcrossprod}
    c(if(is.null(dn.a)) list(NULL) else dn.a[2L - (type != 2L)],
      if(is.null(dn.b)) list(NULL) else dn.b[2L - (type == 3L)])
}

## MJ: Implement forceTriangular() and export this and that?
## MJ: Notably this provides a model for (maybe, in the future) allowing
## forceSymmetric(<non-square>) ... by truncating the "too long" dimension.
forceDiagonal <- function(x, diag = NA_character_) {
    y <- diag(x, names = FALSE) # FIXME? don't allocate if diag == "U"
    cl <- switch(typeof(y),
                 logical = {
                     if(is.na(diag))
                         diag <- if(allTrue(y)) "U" else "N"
                     "ldiMatrix" },
                 integer = {
                     if(is.na(diag))
                         diag <- if(allTrue(y == 1L)) "U" else "N"
                     storage.mode(y) <- "double"
                     "ddiMatrix" },
                 ## integer = {
                 ##     if(is.na(diag))
                 ##         diag <- if(allTrue(y == 1L)) "U" else "N"
                 ##     "idiMatrix" },
                 double = {
                     if(is.na(diag))
                         diag <- if(allTrue(y == 1)) "U" else "N"
                     "ddiMatrix" },
                 complex =
                     stop("complex \"diagonalMatrix\" not yet implemented"),
                 ## complex = {
                 ##     if(is.na(diag))
                 ##         diag <- if(allTrue(y == 1+0i)) "U" else "N"
                 ##     "zdiMatrix" },
                 stop(gettextf("cannot coerce matrix of type \"%s\" to \"diagonalMatrix\"",
                               typeof(y)),
                      domain = NA))
    n <- length(y)
    d <- dim(x)
    dn <- dimnames(x) %||% list(NULL, NULL)
    if(any(d > n)) {
        d <- c(n, n)
        w <- if(d[1L] > n) 1L else 2L
        if(!is.null(dnw <- dn[[w]]))
            dn[[w]] <- dnw[seq_len(n)]
    }
    new(cl, Dim = d, Dimnames = dn, diag = diag,
        x = if(diag == "N") y else y[FALSE])
}

.tCRT <- function(x, lazy = TRUE) .Call(R_sparse_transpose, x, lazy)


.drop0 <- function(x, tol = 0)
    .Call(R_sparse_drop0, x, tol)

drop0 <- function(x, tol = 0, is.Csparse = NA, give.Csparse = TRUE) {
    tryCoerce <-
        if(give.Csparse)
            is.na(is.Csparse) || !is.Csparse
        else if(.M.kind(x) != "n")
            !.isCRT(x)
        else if(all(.M.repr(x) != c("C", "R", "T", "i")))
            TRUE
        else return(x) # n[gst][CRT]Matrix or indMatrix
    if(tryCoerce)
        x <- if(isS4(x)) .M2C(x) else .m2sparse.checking(x, ".", "C")
    .Call(R_sparse_drop0, x, as.double(tol))
}

indDiag <- function(n, upper = TRUE, packed = FALSE)
    .Call(R_index_diagonal, n, packed, upper)

indTri <- function(n, upper = TRUE, diag = FALSE, packed = FALSE)
    .Call(R_index_triangle, n, packed, upper, diag)

## For sparseness handling, return a
## 2-column (i,j) matrix of 0-based indices of non-zero entries:

##' the workhorse for non0ind(.), but occasionally used directly
non0.i <- function(M, cM = class(M), uniqT = TRUE) {
    cld <- getClassDef(cM)
    if(extends(cld, "CsparseMatrix"))
        .Call(compressed_non_0_ij, M, TRUE)
    else if(extends(cld, "RsparseMatrix"))
        .Call(compressed_non_0_ij, M, FALSE)
    else if(extends(cld, "TsparseMatrix")) {
        if(uniqT && is_not_uniqT(M))
	    .Call(compressed_non_0_ij, .M2C(M), TRUE)
	else cbind(M@i, M@j, deparse.level = 0L)
    } else if(extends(cld, "diagonalMatrix")) {
        i <- seq.int(from = 0L, length.out = M@Dim[1L])
        if(M@diag == "N")
	    i <- i[isN0(M@x)]
	cbind(i, i, deparse.level = 0L)
    } else if(extends(cld, "indMatrix")) {
        perm <- M@perm
        i <- seq.int(from = 0L, length.out = length(perm))
        if(M@margin == 1L)
            cbind(i, perm - 1L, deparse.level = 0L)
        else cbind(perm - 1L, i, deparse.level = 0L)
    } else
        stop(gettextf("non0.i() not yet implemented for class %s", dQuote(cM)),
             domain = NA)
}

##' the "more versatile / user" function (still not exported):
non0ind <- function(x, cld = getClassDef(class(x)),
		    uniqT = TRUE, xtendSymm = TRUE, check.Udiag = TRUE)
{
    if(is.numeric(x)) {
        n <- length(x)
        return(if(n == 0L)
                   integer(0L)
               else if(is.matrix(x))
                   arrayInd(seq_len(n)[isN0(x)], dim(x)) - 1L
               else (0:(n-1L))[isN0(x)])
    }
    stopifnot(extends(cld, "sparseMatrix"))

    ij <- non0.i(x, cld, uniqT = uniqT)
    if(xtendSymm && extends(cld, "symmetricMatrix")) {
        ## also get "other" triangle, but not the diagonal again
	notdiag <- ij[, 1L] != ij[, 2L]
	rbind(ij, ij[notdiag, 2:1], deparse.level = 0L)
    } else if(check.Udiag && extends(cld, "triangularMatrix") &&
              x@diag == "U") {
        i <- seq.int(from = 0L, length.out = x@Dim[1L])
        rbind(ij, cbind(i, i, deparse.level = 0L), deparse.level = 0L)
    } else ij
}

##' Decode "encoded" (i,j) indices back to  cbind(i,j)
##' This is the inverse of encodeInd(.)
##'
##' @title Decode "Encoded" (i,j) Indices
##' @param code integer in 0:((n x m - 1)  <==> encodeInd(.) result
##' @param nr the number of rows
##' @return
##' @author Martin Maechler
decodeInd <- function(code, nr)
    cbind(as.integer(code %% nr), as.integer(code %/% nr),
          deparse.level = 0L)

complementInd <- function(ij, dim, orig1=FALSE, checkBnds=FALSE) {
    ## Purpose: Compute the complement of the 2-column 0-based ij-matrix
    ##		but as 1-based indices
    n <- prod(dim)
    if(n == 0L)
        return(integer(0L))
    seq_len(n)[-(1L + .Call(m_encodeInd, ij, dim, orig1, checkBnds))]
}

WhichintersectInd <- function(ij1, ij2, di, orig1=FALSE, checkBnds=FALSE) {
    ## from 2-column (i,j) matrices where i \in {0,.., nrow-1},
    ## find *where*  common entries are in ij1 & ij2
    m1 <- match(.Call(m_encodeInd, ij1, di, orig1, checkBnds),
                .Call(m_encodeInd, ij2, di, orig1, checkBnds))
    ni <- !is.na(m1)
    list(which(ni), m1[ni])
}

### There is a test on this in ../tests/dgTMatrix.R !

uniqTsparse <- function(x, class.x = c(class(x))) {
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## -----------------------------------------------------------
    ## The following is not quite efficient, but easy to program,
    ## and much based on C code
    ##
    ## TODO: faster for the case where 'x' is already 'uniq'?  if(anyDuplicatedT(.))
    if(!extends(class.x, "TsparseMatrix"))
	stop(gettextf("not yet implemented for class \"%s\"", dQuote(class.x)),
	     domain = NA)
    .M2T(.M2C(x))
}

##' non-exported version with*OUT* check -- called often only  if(anyDuplicatedT(.))
.uniqTsparse <- function(x) .M2T(.M2C(x))

asTuniq <- function(x) {
    if(is(x, "TsparseMatrix"))
        uniqTsparse(x)
    else as(x, "TsparseMatrix")
}

## is 'x' a uniq Tsparse Matrix ?
is_not_uniqT <- function(x, di = dim(x))
    is.unsorted(x@j) || anyDuplicatedT(x, di)

## is 'x' a TsparseMatrix with duplicated entries (to be *added* for uniq):
anyDuplicatedT <- function(x, di = dim(x))
    anyDuplicated(.Call(m_encodeInd2, x@i, x@j, di, FALSE, FALSE))

mat2triplet <- function(x, uniqT = FALSE) {
    T <- as(x, "TsparseMatrix")
    if(uniqT && anyDuplicatedT(T)) T <- .uniqTsparse(T)
    if(is(T, "nsparseMatrix"))
         list(i = T@i + 1L, j = T@j + 1L)
    else list(i = T@i + 1L, j = T@j + 1L, x = T@x)
}

.validateCsparse <- function(x, sort.if.needed = FALSE)
    .Call(Csparse_validate2, x, sort.if.needed)

dmperm <- function(x, nAns = 6L, seed = 0L) {
    stopifnot(length(nAns <- as.integer(nAns)) == 1L, nAns %in% c(2L, 4L, 6L),
              length(seed <- as.integer(seed)) == 1L, seed %in% -1:1)
    if(isS4(x)) {
        cld <- getClassDef(class(x))
        if(!extends(cld, "CsparseMatrix"))
            cld <- getClassDef(class(x <- as(x, "CsparseMatrix")))
        if(extends(cld, "symmetricMatrix"))
            cld <- getClassDef(class(x <- .M2gen(x)))
        if(!(extends(cld, "dMatrix") || extends(cld, "nMatrix")))
            x <- .M2kind(x, "d")
    } else { # typically a traditional matrix
        x <- .m2sparse(x, "dgC", NULL, NULL)
    }
    .Call(Csparse_dmperm, x, seed, nAns) # tolerating only [dn][gt]CMatrix 'x'
}

## (matrix|denseMatrix)->denseMatrix as similar as possible to "target"
as_denseClass <- function(x, cl, cld = getClassDef(cl)) {
    kind <- .M.kind(x)
    cl <- .M.nonvirtual(new(cld))
    if(cl == "indMatrix")
        cl <- "ngeMatrix"
    symmetric <- substr(cl, 2L, 2L) == "s" && isSymmetric(x)
    triangular <- !symmetric &&
        substr(cl, 2L, 2L) == "t" && (it <- isTriangular(x))
    packed <- (symmetric || triangular) && substr(cl, 3L, 3L) == "p"
    if(isS4(x)) {
        r <- if(symmetric || triangular)
                 .M2kind(if(symmetric)
                             forceSymmetric(x)
                         else if(attr(it, "kind") == "U")
                             triu(x)
                         else tril(x),
                         kind)
             else .M2gen(x, kind)
        if(packed) pack(r) else r
    }
    else if(symmetric)
        .m2dense(x, paste0(kind, "s", if(packed) "p" else "y"), "U", NULL)
    else if(triangular)
        .m2dense(x, paste0(kind, "t", if(packed) "p" else "r"), attr(it, "kind"), "N")
    else .m2dense(x, paste0(kind, "ge"))
}

## (matrix|sparseMatrix)->CsparseMatrix as similar as possible to "target"
as_CspClass <- function(x, cl, cld = getClassDef(cl)) {
    kind <- .M.kind(x)
    cl <- .M.nonvirtual(new(cld))
    if(cl == "indMatrix")
        cl <- "ngeMatrix"
    symmetric <- substr(cl, 2L, 2L) == "s" && isSymmetric(x)
    triangular <- !symmetric &&
        substr(cl, 2L, 2L) == "t" && (it <- isTriangular(x))
    if(isS4(x)) {
        r <- if(symmetric || triangular)
                 .M2kind(if(symmetric)
                             forceSymmetric(x)
                         else if(attr(it, "kind") == "U")
                             triu(x)
                         else tril(x),
                         kind)
             else .M2gen(x, kind)
        .M2C(r)
    }
    else if(symmetric)
        .m2sparse(x, paste0(kind, "sC"), "U", NULL)
    else if(triangular)
        .m2sparse(x, paste0(kind, "tC"), attr(it, "kind"), "N")
    else .m2sparse(x, paste0(kind, "gC"))
}

## as(<[Mm]atrix>, <non-unit triangular CsparseMatrix>)
asCspN <- function(x) .Call(R_sparse_diag_U2N, as(x, "CsparseMatrix"))

diagU2N <- function (x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "triangularMatrix") && x@diag == "U")
        .diagU2N(x, cl = cl, checkDense = checkDense)
    else x
}

.diagU2N <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(!checkDense && extends(cl, "denseMatrix"))
        x <- as(x, "CsparseMatrix")
    ..diagU2N(x)
}

..diagU2N <- function(x) {
    diag(x) <- TRUE
    x
}

diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(extends(cl, "triangularMatrix") && x@diag == "N")
        .diagN2U(x, cl = cl, checkDense = checkDense)
    else x
}

.diagN2U <- function(x, cl = getClassDef(class(x)), checkDense = FALSE) {
    if(!checkDense & (isDense <- extends(cl, "denseMatrix")))
        ..diagN2U(as(x, "CsparseMatrix"), sparse = TRUE)
    else ..diagN2U(x, sparse = !isDense)
}

..diagN2U <- function(x, sparse) {
    if(sparse && x@Dim[1L] > 0L)
        x <- switch(x@uplo,
                    U = .Call(R_sparse_band, x, 1L, NULL),
                    L = .Call(R_sparse_band, x, NULL, -1L),
                    stop("invalid 'uplo'"))
    x@diag <- "U"
    x
}

## Caches 'value' in the 'factors' slot of 'x' { NOT a copy of 'x' ! }
## and returns 'value'
.set.factor <- function(x, name, value, warn.no.slot = FALSE)
    .Call(R_set_factor, x, name, value, warn.no.slot)

## Empties 'factors' slot of 'x' { NOT a copy of 'x' ! }
## and returns TRUE if 'x' was modified and FALSE if not
.empty.factors <- function(x, warn.no.slot = FALSE)
    .Call(R_empty_factors, x, warn.no.slot)

## all-0 matrix from 'x' which must inherit from Matrix
.setZero <- function(x, kind = .M.kind(x)) {
    cl <- if(.hasSlot(x, "diag"))
              ".tCMatrix"
          else if(.hasSlot(x, "uplo"))
              ".sCMatrix"
          else ".gCMatrix"
    substr(cl, 1L, 1L) <- kind
    d <- x@Dim
    new(cl, Dim = d, Dimnames = x@Dimnames, p = rep.int(0L, d[2L] + 1))
}

##' Compute the three "parts" of two sets:
.setparts <- function(x, y) {
    n1 <- length(m1 <- match(x, y, 0L))
    n2 <- length(m2 <- match(y, x, 0L))
    ix <- which(m1 == 0L)
    iy <- which(m2 == 0L)
    list(x.only = x[ix], ix.only = ix, mx = m1,
         y.only = y[iy], iy.only = iy, my = m2,
         int = if(n1 < n2) y[m1] else x[m2])
}

##' *Only* to be used as function in
##'    setMethod.("Compare", ...., .Cmp.swap)  -->  ./Ops.R  & ./diagMatrix.R
.Cmp.swap <- function(e1, e2) {
    ## "swap RHS and LHS" and use the method below:
    switch(.Generic,
	   "==" =,
           "!=" = callGeneric(e2, e1),
	   "<"	= e2 >	e1,
	   "<=" = e2 >= e1,
	   ">"	= e2 <	e1,
	   ">=" = e2 <= e1)
}

### These two are very similar, the first one has the advantage
### to be applicable to 'Chx' directly:

## FIXME:  kind = "diagBack" is not yet implemented
##	would be much more efficient, but there's no CHOLMOD UI (?)

.diag.dsC <- function(x, Chx = Cholesky(x, LDL = TRUE), res.kind = "diag") {
    force(Chx)
    if(!missing(Chx))
        stopifnot(.CHM.is.LDL(Chx), is.integer(Chx@p), is.double(Chx@x))
    .Call(diag_tC, Chx, res.kind)
    ##    ^^^^^^^ from ../src/Csparse.c
    ## => res.kind in ("trace", "sumLog", "prod", "min", "max", "range", "diag", "diagBack")
}

dimScale <- function(x, d1 = sqrt(1/diag(x, names = FALSE)), d2 = d1) {
    dim.x <- dim(x)
    D1 <- Diagonal(n = dim.x[1L], x = d1)
    D2 <- if(missing(d2)) D1 else Diagonal(n = dim.x[2L], x = d2)
    y <- D1 %*% x %*% D2 # inefficient for symmetricMatrix 'x', but "general"
    if(isS4(x) && is(x, "symmetricMatrix") && identical(d1, d2))
        y <- forceSymmetric(y, x@uplo)
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}

rowScale <- function(x, d) {
    y <- Diagonal(n = nrow(x), x = d) %*% x
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}

colScale <- function(x, d) {
    y <- x %*% Diagonal(n = ncol(x), x = d)
    if(is.list(dn <- dimnames(x)))
        y@Dimnames <- dn
    y
}
