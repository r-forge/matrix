## METHODS FOR CLASS: sparseMatrix, [CRT]sparseMatrix (virtual)
## sparse matrices, in some cases restricted to CSC, CSR, triplet
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.sparse.band <- function(x, k1, k2, ...)
    .Call(R_sparse_band, x, k1, k2)
.sparse.triu <- function(x, k = 0L, ...)
    .Call(R_sparse_band, x, k, NULL)
.sparse.tril <- function(x, k = 0L, ...)
    .Call(R_sparse_band, x, NULL, k)
.sparse.diag.get <- function(x = 1, nrow, ncol, names = TRUE)
    .Call(R_sparse_diag_get, x, names)
.sparse.diag.set <- function(x, value)
    .Call(R_sparse_diag_set, x, value)
.sparse.t <- function(x)
    .Call(R_sparse_transpose, x, "T", FALSE)
.sparse.ct <- function(x)
    .Call(R_sparse_transpose, x, "C", FALSE)
.sparse.fS  <- function(x, uplo = NULL, trans = "C", ...)
    .Call(R_sparse_force_symmetric, x, uplo, trans)
.sparse.symmpart <- function(x, trans = "C", ...)
    .Call(R_sparse_symmpart, x, trans)
.sparse.skewpart <- function(x, trans = "C", ...)
    .Call(R_sparse_skewpart, x, trans)
.sparse.is.sy <- function(object,
                          tol = 100 * .Machine$double.eps,
                          trans = "C", checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    stopifnot(is.numeric(tol), length(tol) == 1L, !is.na(tol))
    ans <- .Call(R_sparse_is_symmetric, object, trans, tol <= 0, checkDN)
    if(!is.na(ans))
        return(ans)
    ## 'object' is an n-by-n [dz]sparseMatrix, n >= 1
    ae <- function(target, current, tolerance, scale = NULL, ...)
        .V.a.e(target = target, current = current,
               tolerance = tolerance, scale = scale,
               check.attributes = FALSE, check.class = FALSE)
    conjugate <- is.complex(object@x) && identical(trans, "C")
    op <- if(conjugate) ct else t
    isTRUE(ae(target = .M2V(object), current = .M2V(op(object)),
              tolerance = tol, ...))
}
.sparse.is.tr <- function(object, upper = NA)
    .Call(R_sparse_is_triangular, object, upper)
.sparse.is.di <- function(object)
    .Call(R_sparse_is_diagonal, object)

setMethod("diff", c(x = "sparseMatrix"),
          ## Mostly cut and paste of base::diff.default :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) != 1L ||
                  lag < 1L || differences < 1L)
                  stop(gettextf("'%s' and '%s' must be positive integers",
                                "lag", "differences"),
                       domain = NA)
              if(lag * differences >= x@Dim[1L])
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences)) {
                  m <- x@Dim[1L]
                  x <- x[i1, , drop = FALSE] -
                      x[-m:-(m - lag + 1L), , drop = FALSE]
              }
              x
          })

setMethod("mean", c(x = "sparseMatrix"),
          function(x, ...) mean(as(x, "sparseVector"), ...))

setMethod("rep", c(x = "sparseMatrix"),
          function(x, ...)  rep(as(x, "sparseVector"), ...))

for(.cl in paste0(c("C", "R", "T"), "sparseMatrix")) {
setMethod("band"  , c(x = .cl), .sparse.band)
setMethod("triu"  , c(x = .cl), .sparse.triu)
setMethod("tril"  , c(x = .cl), .sparse.tril)
setMethod("diag"  , c(x = .cl), .sparse.diag.get)
setMethod("diag<-", c(x = .cl), .sparse.diag.set)
setMethod("t"     , c(x = .cl), .sparse.t)
setMethod("ct"    , c(x = .cl), .sparse.ct)
setMethod("forceSymmetric", c(x = .cl), .sparse.fS)
setMethod("symmpart", c(x = .cl), .sparse.symmpart)
setMethod("skewpart", c(x = .cl), .sparse.skewpart)
setMethod("isSymmetric" , c(object = .cl), .sparse.is.sy)
setMethod("isTriangular", c(object = .cl), .sparse.is.tr)
setMethod("isDiagonal"  , c(object = .cl), .sparse.is.di)
}
rm(.cl)

rm(list = c(grep("^[.]sparse[.](band|tri[ul]|diag[.](get|set)|c?t|fS|symmpart|skewpart|is[.](sy|tr|di))$",
                 ls(all.names = TRUE), value = TRUE)))


if(FALSE) ### FIXME: This would *NOT* be needed, if    as.matrix(<sparseMatrix>) was a no-op ;
         ### -----  and then,  base::scale() -> base::scale.default() would work "magically" already..
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
