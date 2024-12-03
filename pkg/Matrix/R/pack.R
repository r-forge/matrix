## METHODS FOR GENERIC: pack
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("pack", c(x = "matrix"),
          function(x, symmetric = NA, upperTri = NA, ...) {
              if (((sna <- is.na(symmetric)) || symmetric) &&
                  isSymmetric(x, ...))
                  .m2dense(x, ".sp")
              else if ((sna || !symmetric) &&
                       (it <- isTriangular(x, upper = upperTri))) {
                  uplo <-
                      if (is.na(upperTri))
                          attr(it, "kind")
                      else if (upperTri)
                          "U"
                      else "L"
                  .m2dense(x, ".tp", uplo = uplo)
              }
              else if (sna)
                  stop("matrix is not symmetric or triangular")
              else if (symmetric)
                  stop("matrix is not symmetric")
              else stop("matrix is not triangular")
          })

setMethod("pack", c(x = "unpackedMatrix"),
          function(x, symmetric = NA, upperTri = NA, ...) {
              if (.M.shape(x) != "g")
                  .Call(R_dense_as_packed, x, NULL, NULL, NULL)
              else if (((sna <- is.na(symmetric)) || symmetric) &&
                       ({ trans <- "C"; isSymmetric(x, trans = trans, ...) } ||
                        (.M.kind(x) == "z" &&
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
              }
              else if (sna)
                  stop("matrix is not symmetric or triangular")
              else if (symmetric)
                  stop("matrix is not symmetric")
              else stop("matrix is not triangular")
          })

setMethod("pack", c(x = "packedMatrix"),
          function(x, ...) x)


## METHODS FOR GENERIC: unpack
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("unpack", c(x = "matrix"),
          function(x, ...) .m2dense.checking(x, "."))

setMethod("unpack", c(x = "unpackedMatrix"),
          function(x, ...) x)

setMethod("unpack", c(x = "packedMatrix"),
          function(x, ...) .Call(R_dense_as_unpacked, x))
