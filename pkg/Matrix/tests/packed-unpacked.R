## These are tests related to the replacement of many methods for the
## proper subclasses of 'denseMatrix' with methods for the (new, more
## general) virtual subclasses '(un)?packedMatrix'.

library(Matrix)
set.seed(145206)

if (interactive()) {
    options(Matrix.verbose = TRUE, warn = 1, error = recover)
} else {
    options(Matrix.verbose = TRUE, warn = 1)
}

.className <- function(Class) {
    if (is.character(Class))
        Class
    else if (is(Class, "classRepresentation"))
        Class@className
    else
        stop("'Class' is not a class name or a class definition")
}

stopifnotExtends1 <- function(Class1, Class2) {
    if (!extends(Class1, Class2) || Class1@virtual)
        stop(sprintf("'%s' is not a proper subclass of \"%s\"",
                     deparse(substitute(Class1)), className(Class2)))
}

stopifnotExtends2 <- function(Class1, Class2, Class3, Class4) {
    if (!extends(Class3, Class4))
        stop(sprintf("'%s' extends \"%s\" but '%s' does not extend \"%s\"",
                     deparse(substitute(Class1)), className(Class2),
                     deparse(substitute(Class3)), className(Class4)))
}

denseMatrix      <- getClassDef("denseMatrix")
packedMatrix     <- getClassDef("packedMatrix")
unpackedMatrix   <- getClassDef("unpackedMatrix")
generalMatrix    <- getClassDef("generalMatrix")
symmetricMatrix  <- getClassDef("symmetricMatrix")
triangularMatrix <- getClassDef("triangularMatrix")
dMatrix          <- getClassDef("dMatrix")
lMatrix          <- getClassDef("lMatrix")
nMatrix          <- getClassDef("nMatrix")

.PMAP <- c(dgeMatrix = NA,
           lgeMatrix = NA,
           ngeMatrix = NA,
           dsyMatrix = "dspMatrix",
           lsyMatrix = "lspMatrix",
           nsyMatrix = "nspMatrix",
           dtrMatrix = "dtpMatrix",
           ltrMatrix = "ltpMatrix",
           ntrMatrix = "ntpMatrix",
           dpoMatrix = "dppMatrix",
           corMatrix = "dppMatrix",
           Cholesky     = "pCholesky",
           BunchKaufman = "pBunchKaufman")
.UPMAP <- .PMAP[-c(1:3, 11L)]
.UPMAP <- `names<-`(names(.UPMAP), .UPMAP)

U <- function(x) x[upper.tri(x, TRUE)]
L <- function(x) x[lower.tri(x, TRUE)]
`U<-` <- function(x, value) { x[upper.tri(x, TRUE)] <- value; x }
`L<-` <- function(x, value) { x[lower.tri(x, TRUE)] <- value; x }

mkDN <- function(Dim) list(A = paste0("a", seq_len(Dim[1L])),
                           B = paste0("b", seq_len(Dim[2L])))

testDenseClass <- function(Class1, Class2, n) {
    if (!is(Class1, "classRepresentation"))
        Class1 <- getClassDef(Class1)
    stopifnotExtends1(Class1, denseMatrix)

    is.p  <- extends(Class1, packedMatrix)
    is.ge <- extends(Class1, generalMatrix)
    is.sy <- !is.ge && extends(Class1, symmetricMatrix)
    is.d  <- extends(Class1, dMatrix)
    is.l  <- !is.d && extends(Class1, lMatrix)

    if (is.ge)
        Class2 <- NULL
    if (!is.null(Class2)) {
        if (!is(Class2, "classRepresentation"))
            Class2 <- getClassDef(Class2)
        stopifnotExtends1(Class1, denseMatrix)
        if (!is.ge)
            if (is.p)
                stopifnotExtends2(Class1, packedMatrix, Class2, unpackedMatrix)
            else
                stopifnotExtends2(Class1, unpackedMatrix, Class2, packedMatrix)
        if (is.d)
            stopifnotExtends2(Class1, dMatrix, Class2, dMatrix)
        else if (is.l)
            stopifnotExtends2(Class1, lMatrix, Class2, lMatrix)
        else
            stopifnotExtends2(Class1, nMatrix, Class2, nMatrix)
    }

    .mkX <- if (is.d)
                function(n) rlnorm(n) # not rnorm(n), for validObject(<dp[op]M>)
            else
                function(n) sample(c(NA, FALSE, TRUE), n, TRUE)
    mkX <- if (is.p)
               function(Dim) .mkX(choose(Dim[1L]+1L, 2L))
           else
               function(Dim) .mkX(prod(Dim))

    if (is.ge)
        .Dim <- list(c(n, n), c(n+1L, n), c(n, n+1L)) # non-square general
    else
        .Dim <- list(c(n, n))

    .a <- list(Dim = seq_along(.Dim))
    if (!is.ge)
        .a <- c(.a,
                list(uplo = c("U", "L")),
                if (!is.sy) list(diag = c("N", "U")))
    .a <- c(.a, list(stringsAsFactors = FALSE))

    newargs <- do.call(expand.grid, .a)
    newargs[["Dim"]] <- .Dim[.i <- newargs[["Dim"]]]
    newargs[["Dimnames"]] <- lapply(.Dim, mkDN)[.i]
    newargs[["x"]] <- lapply(.Dim, mkX)[.i]
    all(unlist(.mapply(testDenseMatrix,
                       newargs,
                       list(Class1 = Class1, Class2 = Class2))))
}

testDenseMatrix <- function(Class1, Class2, ...) {
    is.p  <- extends(Class1, packedMatrix)
    is.ge <- extends(Class1, generalMatrix)
    is.sy <- !is.ge && extends(Class1, symmetricMatrix)
    is.tr <- !is.ge && !is.sy
    is.bk <- is.tr &&
        extends(Class1, if (is.p) "pBunchKaufman" else "BunchKaufman")
    is.cr <- is.sy && !is.p && extends(Class1, "corMatrix")
    is.d  <- extends(Class1, dMatrix)
    is.l  <- !is.d && extends(Class1, lMatrix)
    is.n  <- !is.d && !is.l

    args <- list(Class = Class1, ...)
    if (is.cr)
        args[["sd"]] <- rep.int(1, args[["Dim"]][1L])
    M <- do.call(new, args)
    cN <- Class1@className

    if (!is.ge)
        if (M@uplo == "U") {
            tri1 <- U; `tri1<-` <- `U<-`
            tri0 <- L; `tri0<-` <- `L<-`
        } else {
            tri1 <- L; `tri1<-` <- `L<-`
            tri0 <- U; `tri0<-` <- `U<-`
        }

    if (is.p) {
        m1 <- array(as.vector(0, typeof(M@x)),
                   dim = M@Dim, dimnames = M@Dimnames)
        tri1(m1) <- M@x
        m2 <- m1
        if (!is.ge)
            if (is.sy) {
                tri0(m2) <- tri0(t(m2))
                dimnames(m2) <- Matrix:::symmDN(M@Dimnames)
            } else if (M@diag == "U") {
                diag(m2) <- as.vector(1, typeof(M@x))
            }

        pM <- M
        upM <- do.call(new, replace(args, c("Class", "x"),
                                    list(.UPMAP[[cN]], as.vector(m1))))
        tM <- do.call(new, replace(args, c("Class", "Dimnames", "x", "uplo"),
                                   list(if (is.bk) "dtpMatrix" else cN,
                                        M@Dimnames[if (is.sy) 1:2 else 2:1],
                                        as.vector(tri0(t(m1))),
                                        switch(M@uplo, U = "L", "U"))))
    } else {
        m1 <- m2 <- array(M@x, dim = M@Dim, dimnames = M@Dimnames)
        if (!is.ge)
            if (is.sy) {
                tri0(m2) <- tri0(t(m2))
                dimnames(m2) <- Matrix:::symmDN(M@Dimnames)
            } else if (M@diag == "U") {
                diag(m2) <- as.vector(1, typeof(M@x))
            }

        upM <- M
        if (!is.ge)
            pM <- do.call(new, replace(if (is.cr) args[names(args) != "sd"]
                                       else args,
                                       c("Class", "x"),
                                       list(.PMAP[[cN]], tri1(m1))))
        tM <- do.call(new, replace(args,
                                   c("Class", "Dim", "Dimnames", "x",
                                     if (!is.ge) "uplo"),
                                   c(list(if (is.bk) "dtrMatrix" else cN,
                                          M@Dim[2:1],
                                          M@Dimnames[if (is.sy) 1:2 else 2:1],
                                          as.vector(t(m1))),
                                     if (!is.ge) switch(M@uplo, U = "L", "U"))))
    }

    stopifnot(identical(unpack(M), upM),
              is.ge || identical(pack(M), pM),
              identical(t(M), tM),
              identical(diag(M, names = FALSE), diag(m2, names = FALSE)),
              identical(diag(M, names = TRUE),  diag(m2, names = TRUE)))
    TRUE
}

.dense.subclasses <- c(names(getClassDef("packedMatrix")@subclasses),
                       names(getClassDef("unpackedMatrix")@subclasses))

mapply(testDenseClass,
       Class1 = .dense.subclasses,
       Class2 = c(.PMAP, .UPMAP)[.dense.subclasses],
       n = 4L)
