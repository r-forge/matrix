## METHODS FOR CLASS: sparseMatrix, [CRT]sparseMatrix (virtual)
## sparse matrices, in some cases restricted to CSC, CSR, triplet
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.sparse.fS  <- function(x, uplo = NULL, trans = "C", ...)
    .Call(R_sparse_force_symmetric, x, uplo, trans)
.sparse.symmpart <- function(x, uplo = "U", trans = "C", ...)
    .Call(R_sparse_symmpart, x, uplo, trans)
.sparse.skewpart <- function(x, trans = "C", ...)
    .Call(R_sparse_skewpart, x, trans)

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
setMethod("forceSymmetric", c(x = .cl), .sparse.fS)
setMethod("symmpart", c(x = .cl), .sparse.symmpart)
setMethod("skewpart", c(x = .cl), .sparse.skewpart)
}
rm(.cl)

rm(list = c(grep("^[.]sparse[.](fS|symmpart|skewpart)$",
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
