## NB: Here, we do _not_ define generic versions of functions that:
##     * are already (S3) generic in their package of origin
##       (e.g., isSymmetric), because our first setMethod call
##       will define a correct (S4) generic function automatically
##     * already have _implicit_ generic definitions in package
##       methods (e.g., qr.R); see methods/R/makeBasicFunsList.R

setGeneric("%&%",
           function(x, y)
               standardGeneric("%&%"))

setGeneric("BunchKaufman",
           function(x, ...)
               standardGeneric("BunchKaufman"))

setGeneric("Cholesky",
           function(A, ...)
               standardGeneric("Cholesky"))

setGeneric("Schur",
           function(x, vectors = TRUE, ...)
               standardGeneric("Schur"),
           signature = "x")

setGeneric("band",
           function(x, k1, k2, ...)
               standardGeneric("band"),
           signature = "x")

setGeneric("ct",
           function(x)
               standardGeneric("ct"))

setGeneric("expand",
           function(x, ...)
               standardGeneric("expand"))

setGeneric("expand1",
           function(x, which, ...)
               standardGeneric("expand1"),
           signature = "x")

setGeneric("expand2",
           function(x, ...)
               standardGeneric("expand2"))

setGeneric("expm",
           function(x)
               standardGeneric("expm"))

setGeneric("facmul",
           function(x, factor, y, trans = FALSE, left = TRUE, ...)
               standardGeneric("facmul"),
           signature = c("x", "y"))

setGeneric("forceCanonical",
           function(x, ...)
               standardGeneric("forceCanonical"))

setGeneric("forceDiagonal",
           function(x, diag = NULL, ...)
               standardGeneric("forceDiagonal"),
           signature = "x")

setGeneric("forceSymmetric",
           function(x, ...)
               standardGeneric("forceSymmetric"))

setGeneric("forceTriangular",
           function(x, ...)
               standardGeneric("forceTriangular"))

setGeneric("isCanonical",
           function(object, ...)
               standardGeneric("isCanonical"))

setGeneric("isDiagonal",
           function(object, ...)
               standardGeneric("isDiagonal"))

setGeneric("isTriangular",
           function(object, upper = NA, ...)
               standardGeneric("isTriangular"),
           signature = "object")

setGeneric("lu",
           function(x, ...)
               standardGeneric("lu"))

setGeneric("nnzero",
           function(x, ...)
               standardGeneric("nnzero"))

setGeneric("pack",
           function(x, ...)
               standardGeneric("pack"))

setGeneric("skewpart",
           function(x, ...)
               standardGeneric("skewpart"))

setGeneric("symmpart",
           function(x, ...)
               standardGeneric("symmpart"))

setGeneric("tril",
           function(x, k = 0L, ...)
               standardGeneric("tril"),
           signature = "x")

setGeneric("triu",
           function(x, k = 0L, ...)
               standardGeneric("triu"),
           signature = "x")

setGeneric("unpack",
           function(x, ...)
               standardGeneric("unpack"))

setGeneric("updown",
           function(update, C, L)
               standardGeneric("updown"))

setGeneric("writeMM",
           function(obj, file, ...)
               standardGeneric("writeMM"))
