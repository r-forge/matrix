## METHODS FOR CLASS: sparseMatrix (virtual)
## sparse matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("mean", signature(x = "sparseMatrix"),
          function(x, ...) mean(as(x, "sparseVector"), ...))

## For very large and very sparse matrices,  the above show()
## is not really helpful;  Use  summary() showing "triplet" as an alternative:

mat2triplet <- function(x, uniqT = FALSE) {
    T <- as(x, "TsparseMatrix")
    if(uniqT && anyDuplicatedT(T)) T <- .uniqTsparse(T)
    if(is(T, "nsparseMatrix"))
         list(i = T@i + 1L, j = T@j + 1L)
    else list(i = T@i + 1L, j = T@j + 1L, x = T@x)
}

setMethod("dim<-", signature(x = "sparseMatrix"),
          function(x, value) {
              if(!is.numeric(value) || length(value) != 2L)
                  stop("dimensions must be numeric of length 2")
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0))
                  stop("dimensions cannot contain negative values")
              if(!is.integer(value)) {
                  if(any(value > .Machine$integer.max))
                      stop("dimensions cannot exceed 2^31-1")
                  value <- as.integer(value)
              }
              if(all(value == (d <- x@Dim)))
                  return(x)
              if((pv <- prod(value)) != (pd <- prod(d)))
                  stop(gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                                pv, pd, domain = NA))
              r <- spV2M(as(x, "sparseVector"),
                         nrow = value[1L], ncol = value[2L])
              ## 'r' is a TsparseMatrix
              if(extends(cd <- getClassDef(class(x)) , "CsparseMatrix"))
                  as(r, "CsparseMatrix")
              else if(extends(cd, "RsparseMatrix"))
                  as(r, "RsparseMatrix")
              else r
          })

setMethod("rep", "sparseMatrix",
          function(x, ...) rep(as(x, "sparseVector"), ...))

if(FALSE) ### FIXME: This would *NOT* be needed, if    as.matrix(<sparseMatrix>) was a no-op ;
          ### -----  and then,  base::scale() -> base::scale.default() would work "magically" already..
## scale() is S3 generic in base
scale.sparseMatrix <- function(x, center = FALSE, scale = TRUE) {
    if(center) warning("a sparseMatrix should rarely be centered: will not be sparse anymore")
    ## x <- as.matrix(x)

    ## This rest is *identically*  == base :: scale.default :
    nc <- ncol(x)
    if (is.logical(center)) {
        if (center) {
            center <- colMeans(x, na.rm=TRUE)
            x <- sweep(x, 2L, center, check.margin=FALSE)
        }
    }
    else if (is.numeric(center) && (length(center) == nc))
        x <- sweep(x, 2L, center, check.margin=FALSE)
    else
        stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
        if (scale) {
            f <- function(v) {
                v <- v[!is.na(v)]
                sqrt(sum(v^2) / max(1, length(v) - 1L))
            }
            scale <- apply(x, 2L, f)
            x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
        }
    }
    else if (is.numeric(scale) && length(scale) == nc)
        x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
        stop("length of 'scale' must equal the number of columns of 'x'")
    if(is.numeric(center)) attr(x, "scaled:center") <- center
    if(is.numeric(scale)) attr(x, "scaled:scale") <- scale
    x
}

.sparse.diag.get <- function(x, nrow, ncol, names = TRUE)
    .Call(R_sparse_diag_get, x, names)
.sparse.diag.set <- function(x, value)
    .Call(R_sparse_diag_set, x, value)
.sparse.band <- function(x, k1, k2, ...) .Call(R_sparse_band, x, k1, k2)
.sparse.triu <- function(x, k = 0L, ...) .Call(R_sparse_band, x, k, NULL)
.sparse.tril <- function(x, k = 0L, ...) .Call(R_sparse_band, x, NULL, k)
.sparse.t    <- function(x)              .Call(R_sparse_transpose, x, FALSE)
.sparse.fS1  <- function(x, uplo) .Call(R_sparse_force_symmetric, x, NULL)
.sparse.fS2  <- function(x, uplo) .Call(R_sparse_force_symmetric, x, uplo)
.sparse.symmpart <- function(x) .Call(R_sparse_symmpart, x)
.sparse.skewpart <- function(x) .Call(R_sparse_skewpart, x)

.C.is.di <- function(object)
    .Call(Csparse_is_diagonal, object)
.R.is.di <- function(object)
    .Call(Rsparse_is_diagonal, object)
.T.is.di <- function(object)
    .Call(Tsparse_is_diagonal, object)

.C.is.tr <- function(object, upper = NA, ...)
    .Call(Csparse_is_triangular, object, upper)
.R.is.tr <- function(object, upper = NA, ...)
    .Call(Rsparse_is_triangular, object, upper)
.T.is.tr <- function(object, upper = NA, ...)
    .Call(Tsparse_is_triangular, object, upper)

.C.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(Csparse_is_symmetric, object, checkDN)
}
.R.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(Rsparse_is_symmetric, object, checkDN)
}
.T.is.sy <- function(object, checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    .Call(Csparse_is_symmetric, as(object, "CsparseMatrix"), checkDN)
}
.sparse.is.sy.dz <- function(object, tol = 100 * .Machine$double.eps,
                             checkDN = TRUE, ...) {
    ## backwards compatibility: don't check DN if check.attributes=FALSE
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    ## be very fast when requiring exact symmetry
    if(tol <= 0) {
        if(!.hasSlot(object, "p"))
            return(.Call(Csparse_is_symmetric, as(object, "CsparseMatrix"), checkDN))
        else if(.hasSlot(object, "i"))
            return(.Call(Csparse_is_symmetric, object, checkDN))
        else
            return(.Call(Rsparse_is_symmetric, object, checkDN))
    }
    ## pretest: is it square?
    d <- object@Dim
    if((n <- d[1L]) != d[2L])
        return(FALSE)
    ## pretest: are DN symmetric in the sense of validObject(<symmetricMatrix>)?
    if(checkDN && !isSymmetricDN(object@Dimnames))
        return(FALSE)
    if(n <= 1L)
        return(TRUE)
    ## now handling an n-by-n [CRT]sparseMatrix, n >= 2:
    x  <- as(  object,  "sparseVector")
    tx <- as(t(object), "sparseVector")
    if(is(tx, "zsparseVector"))
        tx@x <- Conj(tx@x)
    ae <- function(check.attributes, ...) {
        ## discarding possible user-supplied check.attributes:
        all.equal(..., check.attributes = FALSE)
    }
    isTRUE(ae(target = x, current = tx, tolerance = tol, ...))
}

.sparse.subclasses <- names(getClassDef("sparseMatrix")@subclasses)

for (.cl in grep("^[CRT]sparseMatrix$", .sparse.subclasses, value = TRUE)) {
    setMethod("diag",   signature(x = .cl), .sparse.diag.get)
    setMethod("diag<-", signature(x = .cl), .sparse.diag.set)
    setMethod("band", signature(x = .cl), .sparse.band)
    setMethod("triu", signature(x = .cl), .sparse.triu)
    setMethod("tril", signature(x = .cl), .sparse.tril)
    setMethod("t",    signature(x = .cl), .sparse.t)
    setMethod("forceSymmetric", signature(x = .cl, uplo = "missing"),
              .sparse.fS1)
    setMethod("forceSymmetric", signature(x = .cl, uplo = "character"),
              .sparse.fS2)
    setMethod("symmpart", signature(x = .cl), .sparse.symmpart)
    setMethod("skewpart", signature(x = .cl), .sparse.skewpart)
    setMethod("isDiagonal", signature(object = .cl),
              get(paste0(".", substr(.cl, 1L, 1L), ".is.di"),
                  mode = "function", inherits = FALSE))
}

for (.cl in grep("^.g[CRT]Matrix$", .sparse.subclasses, value = TRUE))
    setMethod("isTriangular", signature(object = .cl),
              get(paste0(".", substr(.cl, 3L, 3L), ".is.tr"),
                  mode = "function", inherits = FALSE))
for (.cl in grep("^[lni]g[CRT]Matrix$", .sparse.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl),
              get(paste0(".", substr(.cl, 3L, 3L), ".is.sy"),
                  mode = "function", inherits = FALSE))
for (.cl in grep("^[dz][gt][CRT]Matrix$", .sparse.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .sparse.is.sy.dz)

rm(.cl, .sparse.subclasses, .sparse.is.sy.dz,
   list = c(grep("^[.]sparse[.](band|tri[ul]|t|fS[21]|symmpart|skewpart)$",
                 ls(all.names = TRUE), value = TRUE),
            grep("^[.][CRT][.]is[.](di|tr|sy)$",
                 ls(all.names = TRUE), value = TRUE)))
