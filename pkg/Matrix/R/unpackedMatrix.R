.upM.subclasses <- names(getClass("unpackedMatrix")@subclasses)

.upM.pack <- function(x, ...) {
    .Call(unpackedMatrix_pack, x, NA, NA)
}
.upM.pack.ge <- .m.pack <- function(x, symmetric = NA, upperTri = NA, ...) {
    if(((sna <- is.na(symmetric)) || symmetric) && isSymmetric(x)) {
        .Call(unpackedMatrix_pack, x, FALSE, TRUE)
    } else if((sna || !symmetric) &&
               (tt <- isTriangular(x, upper = upperTri))) {
        upper <- if(is.null(kind <- attr(tt, "kind"))) upperTri else kind == "U"
        .Call(unpackedMatrix_pack, x, TRUE, upper)
    } else {
        kind <- if(is.na(upperTri)) "" else if(upperTri) "upper " else "lower "
        if(sna)
            stop("'x' is not symmetric or ", kind, "triangular")
        else if(symmetric)
            stop("'x' is not symmetric")
        else
            stop("'x' is not ", kind, "triangular")
    }
}
body(.m.pack)[[2L]][[3L]][[2L]][[2L]] <-
body(.m.pack)[[2L]][[4L]][[3L]][[3L]][[2L]] <- quote(matrix_pack)

setMethod("unpack", "unpackedMatrix", function(x, ...) x)
setMethod("pack", "unpackedMatrix", .upM.pack)

for (.cl in grep("^.geMatrix$", .upM.subclasses, value = TRUE))
    setMethod("pack", signature(x = .cl), .upM.pack.ge)

setMethod("pack", signature(x = "matrix"), .m.pack)

.upM.is.sy <- function(object, ...) {
    .Call(unpackedMatrix_is_symmetric, object) # requiring exact symmetry
}
.upM.is.sy.dz <- function(object, tol = 100 * .Machine$double.eps,
                          tol1 = 8 * tol, ...) {
    ## be very fast when requiring exact symmetry
    if (tol <= 0)
        return(.Call(unpackedMatrix_is_symmetric, object))
    ## pretest: is it square?
    d <- object@Dim
    if((n <- d[1L]) != d[2L]) return(FALSE)
    if(n <= 1L) return(TRUE)
    ## now handling n-by-n [dz]..Matrix, n >= 2 :
    if(is(object, "zMatrix")) {
        ge <- "zgeMatrix"
        Cj <- Conj
    } else {
        ge <- "dgeMatrix"
        Cj <- identity
    }
    ## pretest: outermost rows ~= outermost columns? (fast for large asymmetric)
    if(length(tol1)) {
        i. <- if (n <= 4L) 1:n else c(1L, 2L, n-1L, n)
        for(i in i.)
            if(!isTRUE(all.equal(object[i, ], Cj(object[, i]),
                                 tolerance = tol1, ...))) return(FALSE)
    }
    ## followed by slower test
    isTRUE(all.equal(as(     object  , ge),
                     as(Cj(t(object)), ge), tolerance = tol, ...))
}
.upM.is.tr <- function(object, upper = NA, ...) {
    .Call(unpackedMatrix_is_triangular, object, upper)
}
.upM.is.di <- function(object) {
    .Call(unpackedMatrix_is_diagonal, object)
}
.m.is.sy <- function(object, tol = 100 * .Machine$double.eps,
                     tol1 = 8 * tol, ...) {
    if (is.logical(object) || is.integer(object) || tol <= 0)
        .Call(matrix_is_symmetric, object) # requiring exact symmetry
    else
        isSymmetric.matrix(object, tol = tol, tol1 = tol1, ...)
}
.m.is.tr <- function(object, upper = NA, ...) {
    .Call(matrix_is_triangular, object, upper)
}
.m.is.di <- function(object) {
    .Call(matrix_is_diagonal, object)
}

## methods for .syMatrix in ./symmetricMatrix.R
for (.cl in c("unpackedMatrix",
              grep("^[lni](ge|tr)Matrix$", .upM.subclasses, value = TRUE)))
    setMethod("isSymmetric", signature(object = .cl), .upM.is.sy)

for (.cl in grep("^[dz](ge|tr)Matrix$", .upM.subclasses, value = TRUE))
    setMethod("isSymmetric", signature(object = .cl), .upM.is.sy.dz)

## methods for .trMatrix in ./triangularMatrix.R
for (.cl in c("unpackedMatrix",
              grep("^.(ge|sy)Matrix$", .upM.subclasses, value = TRUE)))
    setMethod("isTriangular", signature(object = .cl), .upM.is.tr)

setMethod("isDiagonal", signature(object = "unpackedMatrix"), .upM.is.di)

setMethod("isSymmetric", signature(object = "matrix"), .m.is.sy)
setMethod("isTriangular", signature(object = "matrix"), .m.is.tr)
setMethod("isDiagonal", signature(object = "matrix"), .m.is.di)

setMethod("t", signature(x = "unpackedMatrix"),
          function(x) .Call(unpackedMatrix_t, x))
setMethod("diag", signature(x = "unpackedMatrix"),
          function(x, nrow, ncol, names) .Call(unpackedMatrix_diag_get, x, names))
setMethod("diag<-", signature(x = "unpackedMatrix"),
          function(x, value) .Call(unpackedMatrix_diag_set, x, value))

rm(.upM.pack, .upM.pack.ge, .upM.is.sy, .upM.is.sy.dz, .upM.is.tr, .upM.is.di,
   .m.pack, .m.is.sy, .m.is.tr, .m.is.di, .cl, .upM.subclasses)
