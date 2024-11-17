## METHODS FOR GENERIC: coerce, as.*
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.M2kind <- function(from, kind = ".", sparse = NA)
    .Call(R_Matrix_as_kind, from, kind, sparse)

.M2gen <- function(from, kind = ".")
    .Call(R_Matrix_as_general, from, kind)
..M2gen <- function(from) # for setAs()
    .Call(R_Matrix_as_general, from, ".")

.M2sym <- function(from, ...) {
    if(isSymmetric(from, ...))
        forceSymmetric(from)
    else
        stop("matrix is not symmetric; consider forceSymmetric(.) or symmpart(.)")
}
..M2sym <- .M2sym # for setAs()
formals(..M2sym) <- formals(..M2sym)[-2L]
body(..M2sym)[[2L]][[2L]] <-
    body(..M2sym)[[2L]][[2L]][-3L]

.M2tri <- function(from, ...) {
    if(!(it <- isTriangular(from, ...)))
        stop("matrix is not triangular; consider triu(.) or tril(.)")
    else if(attr(it, "kind") == "U")
        triu(from)
    else
        tril(from)
}
..M2tri <- .M2tri # for setAs()
formals(..M2tri) <- formals(..M2tri)[-2L]
body(..M2tri)[[2L]][[2L]][[2L]][[2L]][[3L]] <-
    body(..M2tri)[[2L]][[2L]][[2L]][[2L]][[3L]][-3L]

.M2diag <- function(from) {
    if (!isDiagonal(from))
        stop("matrix is not diagonal; consider Diagonal(x=diag(.))")
    forceDiagonal.backcomp(from)
}

.M2v <- function(from)
    .Call(R_Matrix_as_vector, from)

.M2m <- function(from)
    .Call(R_Matrix_as_matrix, from)

.M2unpacked <- function(from)
    .Call(R_Matrix_as_unpacked, from)

.M2packed <- function(from)
    .Call(R_Matrix_as_packed, from)

.M2C <- function(from)
    .Call(R_Matrix_as_Csparse, from)

.M2R <- function(from)
    .Call(R_Matrix_as_Rsparse, from)

.M2T <- function(from)
    .Call(R_Matrix_as_Tsparse, from)

.M2V <- function(from)
    .Call(R_Matrix_as_Vector, from)

.sparse2dense <- function(from, packed = FALSE)
    .Call(R_sparse_as_dense, from, packed)

.diag2dense <- function(from, kind = ".", shape = "t", packed = FALSE,
                        uplo = "U", trans = "T")
    .Call(R_diagonal_as_dense, from, kind, shape, packed, uplo, trans)

.ind2dense <- function(from, kind = "n")
    .Call(R_index_as_dense, from, kind)

.dense2sparse <- function(from, repr = "C")
    .Call(R_dense_as_sparse, from, repr)

.diag2sparse <- function(from, kind = ".", shape = "t", repr = "C",
                         uplo = "U", trans = "T")
    .Call(R_diagonal_as_sparse, from, kind, shape, repr, uplo, trans)

.ind2sparse <- function(from, kind = "n", repr = ".")
    .Call(R_index_as_sparse, from, kind, repr)

.v2dense <- function(from, class = ".ge",
                     uplo = "U", trans = "C", diag = "N",
                     nrow = 1L, ncol = 1L, byrow = FALSE,
                     dimnames = NULL)
    .Call(R_vector_as_dense, from, class, uplo, trans, diag,
          nrow, ncol, byrow, dimnames)

.m2dense <- function(from, class = ".ge",
                     uplo = "U", trans = "C", diag = "N", margin = 2L)
    .Call(R_matrix_as_dense, from, class, uplo, trans, diag, margin)

.m2dense.checking <- function(from, kind = ".", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("invalid type \"%s\" in '%s'",
                         typeof(from), ".m2dense.checking"),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid %s=\"%s\" to '%s'",
                                 "kind", kind, ".m2dense.checking"),
                        domain = NA))
        if(kind == "n" && anyNA(from))
            from[is.na(from)] <- TRUE
    }
    if(isSymmetric(from, ...))
        .m2dense(from, paste0(kind, "sy"),
                 uplo = "U", trans = "C")
    else if(it <- isTriangular(from))
        .m2dense(from, paste0(kind, "tr"),
                 uplo = attr(it, "kind"), diag = "N")
    else
        .m2dense(from, paste0(kind, "ge"))
}

.V2sparse <- function(from, class = ".gC",
                      uplo = "U", trans = "C", diag = "N",
                      nrow = 1L, ncol = 1L, byrow = FALSE,
                      dimnames = NULL)
    .Call(R_Vector_as_sparse, from, class, uplo, trans, diag,
          nrow, ncol, byrow, dimnames)

.m2sparse <- function(from, class = ".gC",
                      uplo = "U", trans = "C", diag = "N", margin = 2L)
    .Call(R_matrix_as_sparse, from, class, uplo, trans, diag, margin)

.m2sparse.checking <- function(from, kind = ".", repr = "C", ...) {
    switch(typeof(from), logical =, integer =, double = NULL,
           stop(gettextf("invalid type \"%s\" in '%s'",
                         typeof(from), ".m2sparse.checking"),
                domain = NA))
    if(kind != ".") {
        ## These must happen before isSymmetric() call
        storage.mode(from) <-
            switch(kind, n =, l = "logical", d = "double",
                   stop(gettextf("invalid %s=\"%s\" to '%s'",
                                 "kind", kind, ".m2sparse.checking"),
                        domain = NA))
        if(kind == "n" && anyNA(from))
            from[is.na(from)] <- TRUE
    }
    if(isSymmetric(from, ...))
        .m2sparse(from, paste0(kind, "s", repr),
                  uplo = "U", trans = "C")
    else if(it <- isTriangular(from))
        .m2sparse(from, paste0(kind, "t", repr),
                  uplo = attr(it, "kind"), diag = "N")
    else
        .m2sparse(from, paste0(kind, "g", repr))
}

.m2V <- function(from, kind = ".")
    .Call(R_vector_as_Vector, from, kind)

.V2kind <- function(from, kind = ".") {
    if(kind == ".")
        return(from)
    kind. <- .M.kind(from)
    if(kind == ",")
        kind <- if(kind. == "z") "z" else "d"
    if(kind == kind.)
        return(from)
    to <- new(paste0(kind, "sparseVector"))
    to@length <- from@length
    to@i <- from@i
    if(kind != "n")
        to@x <-
            if(kind. == "n")
                rep.int(switch(kind, "l" = TRUE, "i" = 1L, "d" = 1, "z" = 1+0i), length(from@i))
            else as.vector(from@x, typeof(to@x))
    to
}

.V2v <- function(from) {
    if(.M.kind(from) != "n") {
        to <- vector(typeof(from@x), from@length)
        to[from@i] <- from@x
    } else {
        to <- logical(from@length)
        to[from@i] <- TRUE
    }
    to
}

.V2m <- function(from) {
    if(is.double(m <- length(from)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    to <- .V2v(from)
    dim(to) <- c(m, 1L)
    to
}

.V2a <- function(from) {
    if(is.double(m <- length(from)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    to <- .V2v(from)
    dim(to) <- m
    to
}

.V2unpacked <- function(from) {
    if(is.double(m <- length(from)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "geMatrix"))
    to@Dim <- c(m, 1L)
    to@x <- replace(vector(typeof(to@x), m), from@i,
                    if(kind == "n") TRUE else from@x)
    to
}

.V2C <- function(from) {
    if(is.double(m <- length(from)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "gCMatrix"))
    to@Dim <- c(m, 1L)
    to@p <- c(0L, length(from@i))
    to@i <- as.integer(from@i) - 1L
    if(kind != "n")
        to@x <- if(kind == "i") as.double(from@x) else from@x
    to
}

.V2R <- function(from) {
    if(is.double(m <- length(from)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "gRMatrix"))
    to@Dim <- c(m, 1L)
    to@p <- c(0L, cumsum(replace(logical(m), from@i, TRUE)))
    to@j <- integer(length(from@i))
    if(kind != "n")
        to@x <- if(kind == "i") as.double(from@x) else from@x
    to
}

.V2T <- function(from) {
    if(is.double(m <- length(from)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    kind <- .M.kind(from)
    to <- new(paste0(if(kind == "i") "d" else kind, "gTMatrix"))
    to@Dim <- c(m, 1L)
    to@i <- as.integer(from@i) - 1L
    to@j <- integer(length(from@i))
    if(kind != "n")
        to@x <- if(kind == "i") as.double(from@x) else from@x
    to
}


## ==== To vector ======================================================

## Need 'base' functions calling as.*() to dispatch to our S4 methods:
if (FALSE) {
## 2023-08-10: breaks iGraphMatch, mcmcsae, mcompanion
## which define proper subclasses of Matrix not extending
## any of _our_ proper subclasses of Matrix
as.matrix.Matrix <- function(x, ...) .M2m(x)
 as.array.Matrix <- function(x, ...) .M2m(x)
} else {
as.matrix.Matrix <- function(x, ...) as(x, "matrix")
 as.array.Matrix <- function(x, ...) as(x, "matrix")
setAs("Matrix", "matrix", .M2m)
}
as.matrix.sparseVector <- function(x, ...) .V2m(x)
 as.array.sparseVector <- function(x, ...) .V2a(x)

setMethod("as.vector" , c(x = "Matrix"),
          function(x, mode = "any") as.vector(.M2v(x), mode))
setMethod("as.matrix" , c(x = "Matrix"),
          as.matrix.Matrix)
setMethod("as.array"  , c(x = "Matrix"),
          as.array.Matrix)
setMethod("as.logical", c(x = "Matrix"),
          function(x, ...) as.logical(.M2v(x)))
setMethod("as.integer", c(x = "Matrix"),
          function(x, ...) as.integer(.M2v(x)))
setMethod("as.numeric", c(x = "Matrix"),
          function(x, ...) as.numeric(.M2v(x)))
setMethod("as.complex", c(x = "Matrix"),
          function(x, ...) as.complex(.M2v(x)))

setMethod("as.vector" , c(x = "sparseVector"),
          function(x, mode = "any") as.vector(.V2v(x), mode))
setMethod("as.matrix" , c(x = "sparseVector"),
          as.matrix.sparseVector)
setMethod("as.array"  , c(x = "sparseVector"),
          as.array.sparseVector)
setMethod("as.logical", c(x = "sparseVector"),
          function(x, ...) as.logical(.V2v(x)))
setMethod("as.integer", c(x = "sparseVector"),
          function(x, ...) as.integer(.V2v(x)))
setMethod("as.numeric", c(x = "sparseVector"),
          function(x, ...) as.numeric(.V2v(x)))
setMethod("as.complex", c(x = "sparseVector"),
          function(x, ...) as.complex(.V2v(x)))


## ==== To Matrix ======================================================

setAs("sparseVector", "Matrix",
      .V2C)
setAs("matrix", "Matrix",
      function(from) {
          if(isDiagonal(from))
              forceDiagonal.backcomp(from)
          else if(.sparseDefault(from))
              .m2sparse.checking(from, ".", "C")
          else .m2dense.checking(from, ".")
      })
setAs("vector", "Matrix",
      function(from) {
          if(is.object(from) && length(dim(from)) == 2L) # e.g., data.frame
              as(as.matrix(from), "Matrix")
          else if(.sparseDefault(from))
              .m2sparse(from, ".gC")
          else .m2dense(from, ".ge")
      })
setAs(   "ANY", "Matrix",
      function(from) as(as(from, "matrix"), "Matrix"))

if(FALSE) {
## MJ: not yet ... must first make as(<sparseCholesky>, "Matrix") defunct
setAs("MatrixFactorization", "Matrix",
      function(from) {
          n <- length(x <- expand2(from))
          to <- x[[1L]]
          if(n >= 2L) for(i in 2L:n) to <- to %*% x[[i]]
          to
      })
}


## ==== To sparseVector ================================================

setAs("Matrix", "sparseVector",
      function(from) .M2V(from))
setAs("matrix", "sparseVector",
      function(from) .m2V(from))
setAs("vector", "sparseVector",
      function(from) .m2V(from))
setAs(   "ANY", "sparseVector",
      function(from) .m2V(as.vector(from)))


## ==== To "kind" ======================================================

setAs("Matrix", "nMatrix",
      function(from) .M2kind(from, "n",    NA))
setAs("Matrix", "lMatrix",
      function(from) .M2kind(from, "l",    NA))
setAs("Matrix", "dMatrix",
      function(from) .M2kind(from, "d",    NA))

setAs("Matrix", "ndenseMatrix",
      function(from) .M2kind(from, "n", FALSE))
setAs("Matrix", "ldenseMatrix",
      function(from) .M2kind(from, "l", FALSE))
setAs("Matrix", "ddenseMatrix",
      function(from) .M2kind(from, "d", FALSE))

setAs("Matrix", "nsparseMatrix",
      function(from) .M2kind(from, "n",  TRUE))
setAs("Matrix", "lsparseMatrix",
      function(from) .M2kind(from, "l",  TRUE))
setAs("Matrix", "dsparseMatrix",
      function(from) .M2kind(from, "d",  TRUE))

setAs("matrix", "nMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse.checking(from, "n", "C")
          else .m2dense.checking(from, "n")
      })
setAs("matrix", "lMatrix",
      function(from) {
          if(isDiagonal(from))
              forceDiagonal.backcomp(`storage.mode<-`(from, "logical"))
          else if(.sparseDefault(from))
              .m2sparse.checking(from, "l", "C")
          else .m2dense.checking(from, "l")
      })
setAs("matrix", "dMatrix",
      function(from) {
          if(isDiagonal(from))
              forceDiagonal.backcomp(`storage.mode<-`(from, "double"))
          else if(.sparseDefault(from))
              .m2sparse.checking(from, "d", "C")
          else .m2dense.checking(from, "d")
      })

setAs("matrix", "ndenseMatrix",
      function(from) .m2dense.checking(from, "n"))
setAs("matrix", "ldenseMatrix",
      function(from) .m2dense.checking(from, "l"))
setAs("matrix", "ddenseMatrix",
      function(from) .m2dense.checking(from, "d"))

setAs("matrix", "nsparseMatrix",
      function(from) .m2sparse.checking(from, "n", "C"))
setAs("matrix", "lsparseMatrix",
      function(from) .m2sparse.checking(from, "l", "C"))
setAs("matrix", "dsparseMatrix",
      function(from) .m2sparse.checking(from, "d", "C"))

setAs("vector", "nMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse(from, "ngC")
          else .m2dense(from, "nge")
      })
setAs("vector", "lMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse(from, "lgC")
          else .m2dense(from, "lge")
      })
setAs("vector", "dMatrix",
      function(from) {
          if(.sparseDefault(from))
              .m2sparse(from, "dgC")
          else .m2dense(from, "dge")
      })

setAs("vector", "ndenseMatrix",
      function(from) .m2dense(from, "nge"))
setAs("vector", "ldenseMatrix",
      function(from) .m2dense(from, "lge"))
setAs("vector", "ddenseMatrix",
      function(from) .m2dense(from, "dge"))

setAs("vector", "nsparseMatrix",
      function(from) .m2sparse(from, "ngC"))
setAs("vector", "lsparseMatrix",
      function(from) .m2sparse(from, "lgC"))
setAs("vector", "dsparseMatrix",
      function(from) .m2sparse(from, "dgC"))

setAs("sparseVector", "nsparseVector",
      function(from) .V2kind(from, "n"))
setAs("sparseVector", "lsparseVector",
      function(from) .V2kind(from, "l"))
setAs("sparseVector", "isparseVector",
      function(from) .V2kind(from, "i"))
setAs("sparseVector", "dsparseVector",
      function(from) .V2kind(from, "d"))
setAs("sparseVector", "zsparseVector",
      function(from) .V2kind(from, "z"))

setAs("vector", "nsparseVector",
      function(from) .m2V(from, "n"))
setAs("vector", "lsparseVector",
      function(from) .m2V(from, "l"))
setAs("vector", "isparseVector",
      function(from) .m2V(from, "i"))
setAs("vector", "dsparseVector",
      function(from) .m2V(from, "d"))
setAs("vector", "zsparseVector",
      function(from) .m2V(from, "z"))


## ==== To "shape" =====================================================

..m2gen <- function(from) .m2dense(from, ".ge")

setAs(      "Matrix", "generalMatrix", ..M2gen)
setAs(      "matrix", "generalMatrix", ..m2gen)
setAs(      "vector", "generalMatrix", ..m2gen)
setAs("sparseVector", "generalMatrix",  .V2C)

setAs("Matrix",  "symmetricMatrix", ..M2sym)
setAs("matrix",  "symmetricMatrix", ..M2sym)

setAs("Matrix", "triangularMatrix", ..M2tri)
setAs("matrix", "triangularMatrix", ..M2tri)

rm(..m2gen)

setAs("diagonalMatrix",  "symmetricMatrix",
      function(from) {
          if(!isSymmetricDN(from@Dimnames))
              stop("matrix is not symmetric; consider forceSymmetric(.) or symmpart(.)")
          .diag2sparse(from, ".", "s", "C")
      })

setAs("diagonalMatrix", "triangularMatrix",
      function(from)
          .diag2sparse(from, ".", "t", "C"))


## ==== To "representation" ============================================

setAs("Matrix",    "denseMatrix", .M2unpacked)
setAs("Matrix", "unpackedMatrix", .M2unpacked)
setAs("Matrix",   "packedMatrix", .M2packed)
setAs("Matrix",   "sparseMatrix", .M2C)
setAs("Matrix",  "CsparseMatrix", .M2C)
setAs("Matrix",  "RsparseMatrix", .M2R)
setAs("Matrix",  "TsparseMatrix", .M2T)

## Do test for structure:
## FIXME: wrongly assumes that methods are defined for pack(<sparseMatrix>) ...
setAs("generalMatrix", "packedMatrix", function(from) pack(from))

setAs("matrix",    "denseMatrix",
      function(from)  .m2dense.checking(from, "."))
setAs("matrix", "unpackedMatrix",
      function(from)  .m2dense.checking(from, "."))
setAs("matrix",   "packedMatrix",
      function(from) pack(from))
setAs("matrix",   "sparseMatrix",
      function(from) .m2sparse.checking(from, ".", "C"))
setAs("matrix",  "CsparseMatrix",
      function(from) .m2sparse.checking(from, ".", "C"))
setAs("matrix",  "RsparseMatrix",
      function(from) .m2sparse.checking(from, ".", "R"))
setAs("matrix",  "TsparseMatrix",
      function(from) .m2sparse.checking(from, ".", "T"))

## Many people want this coercion to be available and fast:
setAs("matrix", "dgCMatrix",
      function(from) .m2sparse(from, "dgC"))

setAs("vector",    "denseMatrix",
      function(from)
          if(is.object(from) && length(dim(from)) == 2L) # e.g., data.frame
              as(as.matrix(from),  "denseMatrix")
          else .m2dense(from, ".ge"))
setAs("vector", "unpackedMatrix",
      function(from)  .m2dense(from, ".ge"))
setAs("vector",   "sparseMatrix",
      function(from)
          if(is.object(from) && length(dim(from)) == 2L) # e.g., data.frame
              as(as.matrix(from), "sparseMatrix")
          else .m2sparse(from, ".gC"))
setAs("vector",  "CsparseMatrix",
      function(from) .m2sparse(from, ".gC"))
setAs("vector",  "RsparseMatrix",
      function(from) .m2sparse(from, ".gR"))
setAs("vector",  "TsparseMatrix",
      function(from) .m2sparse(from, ".gT"))

setAs("ANY",  "denseMatrix",
      function(from) as(as(from, "matrix"),  "denseMatrix"))
setAs("ANY", "sparseMatrix",
      function(from) as(as(from, "matrix"), "sparseMatrix"))

setAs("sparseVector",    "denseMatrix", .V2unpacked)
setAs("sparseVector", "unpackedMatrix", .V2unpacked)
setAs("sparseVector",   "sparseMatrix", .V2C)
setAs("sparseVector",  "CsparseMatrix", .V2C)
setAs("sparseVector",  "RsparseMatrix", .V2R)
setAs("sparseVector",  "TsparseMatrix", .V2T)

setAs("Matrix", "diagonalMatrix", .M2diag)
setAs("matrix", "diagonalMatrix", .M2diag)

setAs("Matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))
setAs("matrix", "indMatrix",
      function(from) as(as(from, "nsparseMatrix"), "indMatrix"))

setAs("Matrix",   "pMatrix",
      function(from) as(as(from, "nsparseMatrix"),   "pMatrix"))
setAs("matrix",   "pMatrix",
      function(from) as(as(from, "nsparseMatrix"),   "pMatrix"))


## ==== More to indMatrix ==============================================

.perm2ind <-
function(perm, n, margin = 1L, check.p = 0L) {
    if (mode(perm) != "numeric")
        stop(gettextf("'%s' is not of type \"%s\" or \"%s\"",
                      "perm", "integer", "double"),
             domain = NA)
    else if ((m <- length(perm)) == 0L)
        perm <- integer(0L)
    else if (anyNA(r <- range(perm)))
        stop(gettextf("'%s' contains NA",
                      "perm"),
             domain = NA)
    else if (r[1L] < 1L)
        stop(gettextf("'%s' has elements less than %d",
                      "perm", 1L),
             domain = NA)
    else if (m > .Machine$integer.max ||
            (is.double(perm) && trunc(r[2L]) > .Machine$integer.max))
        stop(gettextf("dimensions cannot exceed %s",
                      "2^31-1"),
             domain = NA)
    else { perm <- as.integer(perm); r <- as.integer(r) }

    if (m.n <- missing(n))
        n <- if (m == 0L) 0L else r[2L]
    else if (mode(n) != "numeric" ||
             length(n) != 1L || is.na(n) || n < 0L)
        stop(gettextf("'%s' is not a non-negative number",
                      "n"),
             domain = NA)
    else if (is.double(n) && trunc(n) > .Machine$integer.max)
        stop(gettextf("dimensions cannot exceed %s",
                      "2^31-1"),
             domain = NA)
    else if (r[2L] > as.integer(n))
        stop(gettextf("'%s' has elements exceeding '%s'",
                      "perm", "n"),
             domain = NA)
    else n <- as.integer(n)

    if (mode(margin) != "numeric" ||
        length(margin) != 1L || is.na(margin) ||
        (margin != 1L && margin != 2L))
        stop(gettextf("'%s' is not %d or %d",
                      "margin", 1L, 2L),
             domain = NA)

    give.p <- check.p >= 1L && m == n &&
        (m == 0L || (all(r == c(1L, m)) && !anyDuplicated.default(perm)))
    if (check.p >= 2L && !give.p)
        stop(gettextf("'%s' is not a permutation of seq_len(%s)",
                      "perm", if (m.n) "max(perm, 0)" else "n"),
             domain = NA)

    J <- new(if (give.p) "pMatrix" else "indMatrix")
    nms <- names(perm)
    if (margin == 1L) {
        J@Dim <- c(m, n)
        J@Dimnames = list(nms, if (give.p) nms)
    } else {
        J@Dim <- c(n, m)
        J@Dimnames = list(if (give.p) nms, nms)
        J@margin <- 2L
    }
    J@perm <- perm
    J
}

setAs("numeric", "indMatrix",
      function(from)
          .perm2ind(from))

## FIXME: deprecate this method and export '.perm2ind'
setAs("list", "indMatrix",
      function(from)
          do.call(.perm2ind, unname(from)))

setAs("nsparseMatrix", "indMatrix",
      function(from) {
          from <- .M2gen(from)
          J <- new("indMatrix")
          J@Dim <- from@Dim
          J@Dimnames <- from@Dimnames
          from. <- .M2R(from)
          p <- from.@p
          m <- length(p) - 1L
          if (all(p == 0:m)) {
              J@perm <- from.@j + 1L
              return(J)
          }
          from. <- .M2C(from)
          p <- from.@p
          n <- length(p) - 1L
          if (all(p == 0:n)) {
              J@perm <- from.@i + 1L
              J@margin <- 2L
              return(J)
          }
          stop("matrix must have exactly one entry in each row or column")
      })


## ==== More to pMatrix ================================================

## MJ: could export without dot
.changeMargin <-
function(x) {
    x@margin <- if (x@margin == 1L) 2L else 1L
    x@perm <- invertPerm(x@perm)
    x
}

setAs("numeric", "pMatrix",
      function(from)
          .perm2ind(from, check.p = 2L))

setAs("nsparseMatrix", "pMatrix",
      function(from) {
          d <- from@Dim
          if ((n <- d[1L]) != d[2L])
              stop(gettextf("attempt to coerce non-square matrix to %s",
                            "pMatrix"),
                   domain = NA)
          from <- .M2gen(from)
          J <- new("pMatrix")
          J@Dim <- d
          J@Dimnames <- from@Dimnames
          from. <- .M2R(from)
          p <- from.@p
          m <- length(p) - 1L
          if (all(p == 0:m) && !anyDuplicated.default(j <- from.@j)) {
              J@perm <- j + 1L
              return(J)
          }
          from. <- .M2C(from)
          p <- from.@p
          n <- length(p) - 1L
          if (all(p == 0:n) && !anyDuplicated.default(i <- from.@i)) {
              J@perm <- i + 1L
              J@margin <- 2L
              return(J)
          }
          stop("matrix must have exactly one entry in each row and column")
      })

setAs("indMatrix", "pMatrix",
      function(from)
          new("pMatrix", from))


## ==== More from MatrixFactorization ==================================

setAs("denseSchur", "dgeMatrix",
      function(from) {
          to <- new("dgeMatrix")
          to@Dim <- from@Dim
          to@x <- from@x
          to
      })

setAs("denseLU", "dgeMatrix",
      function(from) {
          to <- new("dgeMatrix")
          to@Dim <- from@Dim
          to@x <- from@x
          to
      })

setAs("denseBunchKaufman", "dtrMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if (!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if (!packed) to else unpack(to)
      })

setAs("denseBunchKaufman", "dtpMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if (!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if (!packed) pack(to) else to
      })

setAs("denseCholesky", "dtrMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if (!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if (!packed) to else unpack(to)
      })

setAs("denseCholesky", "dtpMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if (!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if (!packed) pack(to) else to
      })

setAs("simplicialCholesky", "dtCMatrix",
      function(from) {
          nz <- from@nz
          k <- sequence.default(nz, from@p[seq_along(nz)] + 1L)

          to <- new("dtCMatrix")
          to@Dim  <- from@Dim
          to@uplo <- "L"
          to@p <- c(0L, cumsum(nz))
          to@i <- from@i[k]
          to@x <- from@x[k]
          to
      })

setAs("supernodalCholesky", "dgCMatrix",
      function(from) {
          super <- from@super
          pi <- from@pi
          b <- length(super)

          nr <- pi[-1L] - pi[-b]
          nc <- super[-1L] - super[-b]
          dp <- rep.int(nr, nc)

          to <- new("dgCMatrix")
          to@Dim <- from@Dim
          to@p <- c(0L, cumsum(dp))
          to@i <- from@s[sequence.default(dp, rep.int(pi[-b] + 1L, nc))]
          to@x <- from@x
          to
      })


## ==== More from/to positive semidefinite =============================

## Operations such as rounding can lose positive semidefiniteness
## but not symmetry, hence:
.indefinite <- function(x) as(x, .M.class(x, 6L))

setAs("Matrix", "posdefMatrix",
      function(from) {
          if (!isSymmetric(from))
              stop("matrix is not Hermitian")
          from <- .M2kind(forceSymmetric(from), ",")
          repr <- .M.repr(from)
          tryCatch(switch(repr,
                          "C" =, "R" =, "T" =
                              suppressWarnings(
                              Cholesky(from, perm = FALSE, LDL = FALSE, super = FALSE)
                              ),
                          "n" =, "p" =
                              Cholesky(from, perm = FALSE)),
                   error = function(e) stop("matrix is not positive definite"))
          z <- is.complex(y <- from@x)
          to <- new(paste0(if (z) "zp" else "dp",
                           if (repr == "n") "o" else repr,
                           "Matrix"))
          to@Dim <- from@Dim
          to@Dimnames <- from@Dimnames
          to@uplo <- from@uplo
          to@factors <- from@factors
          switch(repr,
                 "C" = { to@p <- from@p; to@i <- from@i },
                 "R" = { to@p <- from@p; to@j <- from@j },
                 "T" = { to@i <- from@i; to@j <- from@j })
          to@x <- y
          to
      })

setAs("matrix", "posdefMatrix",
      function(from) {
          if (!isSymmetric(from))
              stop("matrix is not Hermitian")
          if (length(from) > 0L)
          tryCatch(chol(from),
                   error = function(e) stop("matrix is not positive definite"))
          z <- is.complex(from)
          to <- new(if (z) "zpoMatrix" else "dpoMatrix")
          to@Dim <- dim(from)
          if (!is.null(dn <- dimnames(from)))
              to@Dimnames <- symDN(dn)
          to@x <- if (z) as.complex(from) else as.double(from)
          to
      })

setAs("Matrix", "corMatrix",
      function(from) as(as(as(as(from, "dMatrix", strict = FALSE), "posdefMatrix", strict = FALSE), "unpackedMatrix", strict = FALSE), "corMatrix"))

setAs("Matrix", "copMatrix",
      function(from) as(as(as(as(from, "dMatrix", strict = FALSE), "posdefMatrix", strict = FALSE),   "packedMatrix", strict = FALSE), "copMatrix"))

setAs("matrix", "corMatrix",
      function(from) as(as(as(`storage.mode<-`(from, "double"), "posdefMatrix", strict = FALSE), "unpackedMatrix", strict = FALSE), "corMatrix"))

setAs("matrix", "copMatrix",
      function(from) as(as(as(`storage.mode<-`(from, "double"), "posdefMatrix", strict = FALSE),   "packedMatrix", strict = FALSE), "copMatrix"))

setAs("dsyMatrix", "corMatrix",
      function(from) as(as(from, "posdefMatrix"), "corMatrix"))

setAs("dspMatrix", "copMatrix",
      function(from) as(as(from, "posdefMatrix"), "copMatrix"))

setAs("dpoMatrix", "corMatrix",
      function(from) {
          to <- new("corMatrix")
          to@Dim <- d <- from@Dim
          to@Dimnames <- from@Dimnames
          to@uplo <- from@uplo
          if((n <- d[1L]) > 0L) {
          x <- from@x
          k <- seq.int(from = 1L, by = n + 1, length.out = n)
          to@sd <- sd <- sqrt(x[k])
          x <- x / sd / rep(sd, each = n)
          x[k] <- 1
          to@x <- x
          }
          to
      })

setAs("dppMatrix", "copMatrix",
      function(from) {
          to <- new("copMatrix")
          to@Dim <- d <- from@Dim
          to@Dimnames <- from@Dimnames
          to@uplo <- uplo <- from@uplo
          if((n <- d[1L]) > 0L) {
          x <- from@x
          if(uplo == "U") {
              r <- 1L:n
              s <- 1L
              k <- cumsum(r)
          } else {
              r <- n:1L
              s <- 1L:n
              k <- cumsum(c(1L, r[-n]))
          }
          to@sd <- sd <- sqrt(x[k])
          x <- x / rep.int(sd, r) / sd[sequence.default(r, s)]
          x[k] <- 1
          to@x <- x
          }
          to
      })
