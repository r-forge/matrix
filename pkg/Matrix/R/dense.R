## METHODS FOR CLASS: denseMatrix (virtual)
## dense matrices with unpacked _or_ packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.dense.band <- function(x, k1, k2, ...)
    .Call(R_dense_band, x, k1, k2)
.dense.triu <- function(x, k = 0L, ...)
    .Call(R_dense_band, x, k, NULL)
.dense.tril <- function(x, k = 0L, ...)
    .Call(R_dense_band, x, NULL, k)
.dense.diag.get <- function(x = 1, nrow, ncol, names = TRUE)
    .Call(R_dense_diag_get, x, names)
.dense.diag.set <- function(x, value)
    .Call(R_dense_diag_set, x, value)
.dense.t <- function(x)
    .Call(R_dense_transpose, x, "T")
.dense.ct <- function(x)
    .Call(R_dense_transpose, x, "C")
.dense.fS  <- function(x, uplo = NULL, trans = "C", ...)
    .Call(R_dense_force_symmetric, x, uplo, trans)
.dense.symmpart <- function(x, trans = "C", ...)
    .Call(R_dense_symmpart, x, "U", trans)
.dense.skewpart <- function(x, trans = "C", ...)
    .Call(R_dense_skewpart, x, trans)
.dense.is.sy <- function(object,
                         tol = 100 * .Machine$double.eps, tol1 = 8 * tol,
                         trans = "C", checkDN = TRUE, ...) {
    if(checkDN) {
        ca <- function(check.attributes = TRUE, ...) check.attributes
        checkDN <- ca(...)
    }
    stopifnot(is.numeric(tol), length(tol) == 1L, !is.na(tol))
    ans <- .Call(R_dense_is_symmetric, object, trans, tol <= 0, checkDN)
    if(!is.na(ans))
        return(ans)
    ## 'object' is an n-by-n [dz]denseMatrix, n >= 1
    ae <- function(target, current, tolerance, scale = NULL, ...)
        all.equal.numeric(target = target, current = current,
                          tolerance = tolerance, scale = scale,
                          check.attributes = FALSE, check.class = FALSE)
    conjugate <- is.complex(object@x) && identical(trans, "C")
    if(length(tol1) && (n <- object@Dim[1L]) > 1L) {
        op <- if(conjugate) Conj else identity
        for(i in if(n > 4L) c(1L, 2L, n - 1L, n) else 1L:n)
            if(is.character(ae(target = object[i, ], current = op(object[, i]),
                               tolerance = tol1, ...)))
                return(FALSE)
    }
    object <- .M2gen(object)
    op <- if(conjugate) ct else t
    isTRUE(ae(target = object@x, current = op(object)@x,
              tolerance = tol, ...))
}
.dense.is.tr <- function(object, upper = NA)
    .Call(R_dense_is_triangular, object, upper)
.dense.is.di <- function(object)
    .Call(R_dense_is_diagonal, object)

setMethod("diff", c(x = "denseMatrix"),
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

setMethod("mean"  , c(x = "denseMatrix"),
          function(x, ...) mean.default(.M2v(x), ...))

setMethod("rep"   , c(x = "denseMatrix"),
          function(x, ...)          rep(.M2v(x), ...))

setMethod("band"  , c(x = "denseMatrix"), .dense.band)

setMethod("triu"  , c(x = "denseMatrix"), .dense.triu)

setMethod("tril"  , c(x = "denseMatrix"), .dense.tril)

setMethod("diag"  , c(x = "denseMatrix"), .dense.diag.get)

setMethod("diag<-", c(x = "denseMatrix"), .dense.diag.set)

setMethod("t"     , c(x = "denseMatrix"), .dense.t)

setMethod("ct"    , c(x = "denseMatrix"), .dense.ct)

setMethod("forceSymmetric", c(x = "denseMatrix"), .dense.fS)

setMethod("symmpart", c(x = "denseMatrix"), .dense.symmpart)

setMethod("skewpart", c(x = "denseMatrix"), .dense.skewpart)

setMethod("isSymmetric" , c(object = "denseMatrix"), .dense.is.sy)

setMethod("isTriangular", c(object = "denseMatrix"), .dense.is.tr)

setMethod("isDiagonal"  , c(object = "denseMatrix"), .dense.is.di)


## METHODS FOR CLASS: unpackedMatrix (virtual)
## dense matrices with unpacked storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("unpack", c(x = "packedMatrix"),
          function(x, ...) .Call(R_dense_as_unpacked, x))

setMethod("pack", c(x = "packedMatrix"),
          function(x, ...) x)


## METHODS FOR CLASS: packedMatrix (virtual)
## dense matrices with packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.uM.pack <-
function(x, ...) .Call(R_dense_as_packed, x, NULL, NULL, NULL)

.uM.pack.ge <-
function(x, symmetric = NA, upperTri = NA, ...) {
    if(((sna <- is.na(symmetric)) || symmetric) &&
       ({ trans <- "C"; isSymmetric(x, trans = trans, ...) } ||
        (is.complex(x@x) &&
         { trans <- "T"; isSymmetric(x, trans = trans, ...) })))
        .Call(R_dense_as_packed, x, "U", trans, NULL)
    else if((sna || !symmetric) &&
            (it <- isTriangular(x, upper = upperTri))) {
        uplo <-
            if(is.na(upperTri))
                attr(it, "kind")
            else if(upperTri)
                "U"
            else "L"
        .Call(R_dense_as_packed, x, uplo, NULL, "N")
    } else {
        if(sna)
            stop("matrix is not symmetric or triangular")
        else if(symmetric)
            stop("matrix is not symmetric")
        else stop("matrix is not triangular")
    }
}

setMethod("unpack", c(x = "unpackedMatrix"),
          function(x, ...) x)

setMethod("pack", c(x = "unpackedMatrix"), .uM.pack)

.uM.subclasses <- names(getClassDef("unpackedMatrix")@subclasses)
for(.cl in grep("^[nlidz]geMatrix$", .uM.subclasses, value = TRUE))
setMethod("pack", c(x = .cl), .uM.pack.ge)
rm(.cl, .uM.subclasses)


## METHODS FOR CLASS: matrix
## traditional matrices, which really are "dense"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.m.pack <- .uM.pack.ge
body(.m.pack)[[2L]][[2L]][[3L]] <-
    quote(isSymmetric(x, ...))
body(.m.pack)[[2L]][[3L]] <-
    quote(.m2dense(x, ".sp"))
body(.m.pack)[[2L]][[4L]][[3L]][[3L]] <-
    quote(.m2dense(x, ".tp", uplo = uplo))

setMethod("unpack", c(x = "matrix"),
          function(x, ...) .m2dense.checking(x, "."))
setMethod("pack", c(x = "matrix"), .m.pack)
setMethod("band", c(x = "matrix"), .dense.band)
setMethod("triu", c(x = "matrix"), .dense.triu)
setMethod("tril", c(x = "matrix"), .dense.tril)
setMethod("ct", c(x = "matrix"),
          function(x) if(is.complex(x)) Conj(t(x)) else t(x))
setMethod("forceSymmetric", c(x = "matrix"),
          function(x, uplo = "U", trans = "C", ...)
              .m2dense(x, ".sy", uplo = uplo, trans = trans))
setMethod("symmpart", c(x = "matrix"),
          function(x, trans = "C", ...) {
              op <- if(is.complex(x) && identical(trans, "C")) Conj else identity
              symmetrizeDN(0.5 * (x + op(t(x))))
          })
setMethod("skewpart", c(x = "matrix"),
          function(x, trans = "C", ...) {
              op <- if(is.complex(x) && identical(trans, "C")) Conj else identity
              symmetrizeDN(0.5 * (x - op(t(x))))
          })
setMethod("isTriangular", c(object = "matrix"), .dense.is.tr)
setMethod("isDiagonal"  , c(object = "matrix"), .dense.is.di)

rm(.uM.pack, .uM.pack.ge, .m.pack,
   list = c(grep("^[.]dense[.](band|tri[ul]|diag[.](get|set)|c?t|fS|symmpart|skewpart|is[.](sy|tr|di))$",
                 ls(all.names = TRUE), value = TRUE)))
