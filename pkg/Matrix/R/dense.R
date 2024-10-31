## METHODS FOR CLASS: denseMatrix (virtual)
## dense matrices with unpacked _or_ packed storage
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

rm(.uM.pack, .uM.pack.ge, .m.pack)
