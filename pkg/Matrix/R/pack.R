## METHODS FOR GENERIC: pack
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.M.pack <- .m.pack <-
function(x, symmetric = NA, upperTri = NA, ...) {
    if (((sna <- is.na(symmetric)) || symmetric) &&
        ({ trans <- "C"; isSymmetric(x, trans = trans, ...) } ||
         (is.complex(x@x) &&
          { trans <- "T"; isSymmetric(x, trans = trans, ...) })))
        .Call(R_dense_as_packed, x, "U", trans, NULL)
    else if ((sna || !symmetric) &&
             (it <- isTriangular(x, upper = upperTri))) {
        uplo <-
            if (is.na(upperTri))
                attr(it, "kind")
            else if (upperTri)
                "U"
            else "L"
        .Call(R_dense_as_packed, x, uplo, NULL, "N")
    } else {
        if (sna)
            stop("matrix is not symmetric or triangular")
        else if (symmetric)
            stop("matrix is not symmetric")
        else stop("matrix is not triangular")
    }
}
body(.m.pack)[[2L]][[2L]][[3L]] <-
    quote(isSymmetric(x, ...))
body(.m.pack)[[2L]][[3L]] <-
    quote(.m2dense(x, ".sp"))
body(.m.pack)[[2L]][[4L]][[3L]][[3L]] <-
    quote(.m2dense(x, ".tp", uplo = uplo))

setMethod("pack", c(x = "matrix"), .m.pack)

setMethod("pack", c(x = "unpackedMatrix"),
          function(x, ...)
              .Call(R_dense_as_packed, x, NULL, NULL, NULL))

for (.cl in grep("^[nlidz]geMatrix$",
                 names(getClassDef("unpackedMatrix")@subclasses),
                 value = TRUE))
setMethod("pack", c(x = .cl), .M.pack)

setMethod("pack", c(x = "packedMatrix"),
          function(x, ...) x)

rm(.cl, .M.pack, .m.pack)


## METHODS FOR GENERIC: unpack
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("unpack", c(x = "matrix"),
          function(x, ...)
              .m2dense.checking(x, "."))

setMethod("unpack", c(x = "unpackedMatrix"),
          function(x, ...) x)

setMethod("unpack", c(x = "packedMatrix"),
          function(x, ...)
              .Call(R_dense_as_unpacked, x))
