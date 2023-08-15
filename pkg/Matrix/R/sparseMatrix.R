## METHODS FOR CLASS: sparseMatrix (virtual)
## sparse matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~ CONSTRUCTORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

spMatrix <- function(nrow, ncol, i = integer(), j = integer(), x = double())
    new(paste0(.M.kind(x), "gTMatrix"), # rely on new() to check validity
        Dim = c(as.integer(nrow), as.integer(ncol)),
        i = as.integer(i) - 1L,
        j = as.integer(j) - 1L,
        x = if(is.integer(x)) as.double(x) else x)

sparseMatrix <- function(i, j, p, x, dims, dimnames,
                         symmetric = FALSE, triangular = FALSE, index1 = TRUE,
                         repr = c("C", "R", "T"), giveCsparse,
                         check = TRUE, use.last.ij = FALSE)
{
    if((m.i <- missing(i)) + (m.j <- missing(j)) + (m.p <- missing(p)) != 1L)
        stop("exactly one of 'i', 'j', and 'p' must be missing from call")
    if(symmetric && triangular)
        stop("use Diagonal() to construct diagonal (symmetric && triangular) sparse matrices")
    index1 <- as.logical(index1) # allowing {0,1}

    repr <-
        ## NB: prior to 2020-05, we had 'giveCsparse' {T->"C" [default], F->"T"}
        ##     but no 'repr' ... the following is to remain backwards compatible
        if(missing(giveCsparse))
            match.arg(repr)
        else if(!missing(repr)) {
            warning("'giveCsparse' is deprecated; using 'repr' instead")
            match.arg(repr)
        ## } else {
        ##     repr <- if(giveCsparse) "C" else "T"
        ##     warning(gettextf("'giveCsparse' is deprecated; setting repr=\"%s\" for you", repr),
        ##             domain = NA)
        ## }
        } else if(giveCsparse) {
            ## NOT YET:
            ## warning("'giveCsparse' is deprecated; setting repr=\"C\" for you")
            "C"
        } else {
            warning("'giveCsparse' is deprecated; setting repr=\"T\" for you")
            "T"
        }

    if(!m.p) {
        p <- as.integer(p)
        if((n.p <- length(p)) == 0L || anyNA(p) || p[1L] != 0L ||
           any((dp <- p[-1L] - p[-n.p]) < 0L))
            stop("'p' must be a nondecreasing vector c(0, ...)")
        if((n.dp <- length(dp)) > .Machine$integer.max)
            stop("dimensions cannot exceed 2^31-1")
        i. <- rep.int(seq.int(from = 0L, length.out = n.dp), dp)
        if(m.i) i <- i. else j <- i.
    }

    if(!m.i)
        i <- if(index1) as.integer(i) - 1L else as.integer(i) # need 0-index
    if(!m.j)
        j <- if(index1) as.integer(j) - 1L else as.integer(j) # need 0-index

    rij <- cbind(if(n.i <- length(i)) range(i) else 0:-1,
                 if(n.j <- length(j)) range(j) else 0:-1,
                 deparse.level = 0L)
    if(anyNA(rij))
        stop("'i' and 'j' must not contain NA") # and not overflow
    if(any(rij[1L, ] < 0L))
        stop("'i' and 'j' must be ", if(index1) "positive" else "non-negative")
    dims <-
        if(!missing(dims)) {
            if(length(dims) != 2L ||
               any(is.na(dims) | dims < 0L | dims >= .Machine$integer.max + 1))
                stop("invalid 'dims'")
            if(any(dims - 1L < rij[2L, ]))
                stop("'dims' must contain all (i,j) pairs")
            as.integer(dims)
        } else if(symmetric || triangular)
            rep.int(max(rij), 2L) + 1L
        else rij[2L, ] + 1L

    kind <- if(m.x <- missing(x)) "n" else .M.kind(x)
    shape <-
        if(symmetric) {
            if(dims[1L] != dims[2L])
                stop("symmetric matrix must be square")
            "s"
        } else if(triangular) {
            if(dims[1L] != dims[2L])
                stop("triangular matrix must be square")
            "t"
        } else "g"

    r <- new(paste0(kind, shape, "TMatrix"))
    r@Dim <- dims
    if(!missing(dimnames) && !is.null(dimnames))
        r@Dimnames <-
            if(is.character(validDN(dimnames, dims)))
                dimnames
            else fixupDN(dimnames) # needs a valid argument
    if((symmetric || triangular) && all(i >= j))
        r@uplo <- "L" # else "U", the prototype
    if(!m.x) {
        if(is.integer(x))
            x <- as.double(x)
        if((n.x <- length(x)) > 0L && n.x != n.i) {
            if(n.x < n.i) {
                if(n.i %% n.x != 0L)
                    warning(if(m.i) "p[length(p)] " else "length(i) ",
                            "is not an integer multiple of length(x)")
                x <- rep_len(x, n.i) # recycle
            } else if(n.x == 1L)
                x <- x[0L] # tolerate length(i) = 0, length(x) = 1
            else stop("length(x) must not exceed ",
                      if(m.i) "p[length(p)]" else "length(i)")
        }
        if(use.last.ij && n.i == n.j &&
           anyDuplicated.matrix(ij <- cbind(i, j, deparse.level = 0L),
                                fromLast = TRUE)) {
            which.not.dup <- which(!duplicated(ij, fromLast = TRUE))
            i <- i[which.not.dup]
            j <- j[which.not.dup]
            x <- x[which.not.dup]
        }
        r@x <- x
    }
    r@i <- i
    r@j <- j

    if(check)
        validObject(r)
    switch(repr, "C" = .M2C(r), "T" = r, "R" = .M2R(r),
           ## should never happen:
           stop("invalid 'repr'; must be \"C\", \"R\", or \"T\""))
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
              Is <- Diagonal(x = Is)# , names = TRUE
              r <- Is %*% V %*% Is
              r[cbind(1:p,1:p)] <- 1 # exact in diagonal
              as(`dimnames<-`(r, symmDN(dimnames(V))), "symmetricMatrix")
              ## as(r, "symmetricMatrix")
          })

## all.equal(): similar to all.equal_Mat() in ./Matrix.R ;
## -----------	eventually defer to  "sparseVector" methods:
setMethod("all.equal", c(target = "sparseMatrix", current = "sparseMatrix"),
          function(target, current, check.attributes = TRUE, ...) {
              msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
              if(is.list(msg)) msg[[1]]
              else .a.e.comb(msg,
                             all.equal(as(target, "sparseVector"), as(current, "sparseVector"),
                                       check.attributes=check.attributes, ...))
          })
setMethod("all.equal", c(target = "sparseMatrix", current = "ANY"),
          function(target, current, check.attributes = TRUE, ...) {
              msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
              if(is.list(msg)) msg[[1]]
              else .a.e.comb(msg,
                             all.equal(as(target, "sparseVector"), current,
                                       check.attributes=check.attributes, ...))
          })
setMethod("all.equal", c(target = "ANY", current = "sparseMatrix"),
          function(target, current, check.attributes = TRUE, ...) {
              msg <- attr.all_Mat(target, current, check.attributes=check.attributes, ...)
              if(is.list(msg)) msg[[1]]
              else .a.e.comb(msg,
                             all.equal(target, as(current, "sparseVector"),
                                       check.attributes=check.attributes, ...))
          })


setMethod("writeMM", "sparseMatrix",
          function(obj, file, ...)
              writeMM(as(obj, "CsparseMatrix"), as.character(file), ...))

### --- sparse model matrix,  fac2sparse, etc ----> ./spModels.R

###  xtabs(*, sparse = TRUE) ---> part of standard package 'stats' since R 2.10.0

##' @title Random Sparse Matrix
##' @param nrow,
##' @param ncol number of rows and columns, i.e., the matrix dimension
##' @param nnz number of non-zero entries
##' @param rand.x random number generator for 'x' slot
##' @param ... optionally further arguments passed to sparseMatrix()
##' @return a sparseMatrix of dimension (nrow, ncol)
##' @author Martin Maechler
##' @examples M1 <- rsparsematrix(1000, 20, nnz = 200)
##'           summary(M1)
if(FALSE) ## better version below
rsparsematrix <- function(nrow, ncol, nnz,
                          rand.x = function(n) signif(rnorm(nnz), 2),
                          warn.nnz = TRUE, ...) {
    maxi.sample <- 2^31 # maximum n+1 for which sample(n) returns integer
    stopifnot((nnz <- as.integer(nnz)) >= 0,
              nrow >= 0, ncol >= 0, nnz <= nrow * ncol,
              nrow < maxi.sample, ncol < maxi.sample)
    ## to ensure that nnz is strictly followed, must act on duplicated (i,j):
    i <- sample.int(nrow, nnz, replace = TRUE)
    j <- sample.int(ncol, nnz, replace = TRUE)
    dim <- c(nrow, ncol)
    it <- 0
    while((it <- it+1) < 100 &&
          anyDuplicated(n.ij <- encodeInd2(i, j, dim, checkBnds = FALSE))) {
        m <- length(k.dup <- which(duplicated(n.ij)))
        Matrix.msg(sprintf("%3g duplicated (i,j) pairs", m), .M.level = 2)
        if(runif(1) <= 1/2)
            i[k.dup] <- sample.int(nrow, m, replace = TRUE)
        else
            j[k.dup] <- sample.int(ncol, m, replace = TRUE)
    }
    if(warn.nnz && it == 100 &&
       anyDuplicated(encodeInd2(i, j, dim, checkBnds = FALSE)))
        warning("number of non zeros is smaller than 'nnz' because of duplicated (i,j)s")
    sparseMatrix(i = i, j = j, x = rand.x(nnz), dims = dim, ...)
}

## No warn.nnz needed, as we sample the encoded (i,j) with*out* replacement:
rsparsematrix <- function(nrow, ncol, density,
                          nnz = round(density * maxE), symmetric = FALSE,
                          rand.x = function(n) signif(rnorm(n), 2), ...) {
    maxE <- if(symmetric) nrow*(nrow+1)/2 else nrow*ncol
    stopifnot((nnz <- as.integer(nnz)) >= 0,
              nrow >= 0, ncol >= 0, nnz <= maxE)
    ## sampling with*out* replacement (replace=FALSE !):
    ijI <- -1L +
        if(symmetric) sample(indTri(nrow, diag=TRUE), nnz)
        else sample.int(maxE, nnz)
    ## i,j below correspond to  ij <- decodeInd(code, nr) :
    if(is.null(rand.x))
        sparseMatrix(i = ijI  %% nrow,
                     j = ijI %/% nrow,
                     index1 = FALSE,
                     symmetric = symmetric,
                     dims = c(nrow, ncol), ...)
    else
        sparseMatrix(i = ijI  %% nrow,
                     j = ijI %/% nrow,
                     x = rand.x(nnz),
                     index1 = FALSE,
                     symmetric = symmetric,
                     dims = c(nrow, ncol), ...)
}

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
