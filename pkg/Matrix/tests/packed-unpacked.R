## These are tests related to the replacement of many methods for the
## proper subclasses of 'denseMatrix' with methods for the (new, more
## general) virtual subclasses '(un)?packedMatrix'.

## for R_DEFAULT_PACKAGES=NULL :
library(stats)

library(Matrix)
set.seed(145206)

if (interactive()) {
    options(Matrix.verbose = TRUE, warn = 1, error = recover)
} else options(Matrix.verbose = TRUE, warn = 1)

U <- function(x, diag = FALSE) x[upper.tri(x, diag)]
L <- function(x, diag = FALSE) x[lower.tri(x, diag)]
`U<-` <- function(x, diag = FALSE, value) { x[upper.tri(x, diag)] <- value; x }
`L<-` <- function(x, diag = FALSE, value) { x[lower.tri(x, diag)] <- value; x }

mkDN <- function(Dim) list(A = paste0("a", seq_len(Dim[1L])),
                           B = paste0("b", seq_len(Dim[2L])))

denseMatrix      <- getClassDef("denseMatrix")
packedMatrix     <- getClassDef("packedMatrix")
unpackedMatrix   <- getClassDef("unpackedMatrix")
generalMatrix    <- getClassDef("generalMatrix")
symmetricMatrix  <- getClassDef("symmetricMatrix")
triangularMatrix <- getClassDef("triangularMatrix")
dMatrix          <- getClassDef("dMatrix")
lMatrix          <- getClassDef("lMatrix")
nMatrix          <- getClassDef("nMatrix")

## FIXME: Implement in C??
unpackedClass <- function(packedClass) {
    getClassDef(c(dspMatrix = "dsyMatrix",
                  lspMatrix = "lsyMatrix",
                  nspMatrix = "nsyMatrix",
                  dtpMatrix = "dtrMatrix",
                  ltpMatrix = "ltrMatrix",
                  ntpMatrix = "ntrMatrix",
                  dppMatrix = "dpoMatrix",
                  copMatrix = "corMatrix")[[packedClass@className]])
}
packedClass <- function(unpackedClass) {
    getClassDef(c(dgeMatrix = NA,
                  lgeMatrix = NA,
                  ngeMatrix = NA,
                  dsyMatrix = "dspMatrix",
                  lsyMatrix = "lspMatrix",
                  nsyMatrix = "nspMatrix",
                  dtrMatrix = "dtpMatrix",
                  ltrMatrix = "ltpMatrix",
                  ntrMatrix = "ntpMatrix",
                  dpoMatrix = "dppMatrix",
                  corMatrix = "copMatrix")[[unpackedClass@className]])
}
...Class <- function(denseClass) {
    cl <- "...Matrix"
    substr(cl, 1L, 1L) <-
        if (extends(denseClass, dMatrix))
            "d"
        else if (extends(denseClass, lMatrix))
            "l"
        else
            "n"
    substr(cl, 2L, 3L) <-
        if (extends(denseClass, generalMatrix) -> g)
            "ge"
        else if (extends(denseClass, symmetricMatrix))
            "sy"
        else
            "tr"
    if (!g && extends(denseClass, packedMatrix))
        substr(cl, 3L, 3L) <- "p"
    getClassDef(cl)
}

## Tests methods for packed (unpacked) class 'Class'
## using randomly generated matrices of size 'n'
testDenseClass <- function(Class, n) {
    if (!is(Class, "classRepresentation"))
        Class <- getClassDef(Class)
    stopifnot(extends(Class, denseMatrix), !isVirtualClass(Class))

    is.p  <- extends(Class, packedMatrix)
    is.ge <- extends(Class, generalMatrix)
    is.sy <- !is.ge && extends(Class, symmetricMatrix)
    is.tr <- !is.ge && !is.sy
    is.d  <- extends(Class, dMatrix)
    is.l  <- !is.d && extends(Class, lMatrix)
    is.n  <- !is.d && !is.l

    ## For randomly generating matrix data
    .mkX <- if (is.d)
                function(n) rlnorm(n) # not rnorm(n), for validObject(<dp[op]M>)
            else
                function(n) sample(c(NA, FALSE, TRUE), n, TRUE)
    mkX <- if (is.p)
               function(Dim) .mkX((Dim[1L] * (Dim[1L] + 1L)) %/% 2L)
           else
               function(Dim) .mkX(prod(Dim))

    ## "Full factorial" design with varying 'Dim', 'uplo', 'diag' slots
    .Dim <- c(list(c(n, n)), if (is.ge) list(c(n+1L, n), c(n, n+1L)))
    .a <- c(list(Dim = seq_along(.Dim)),
            if (!is.ge) list(uplo = c("U", "L")),
            if ( is.tr) list(diag = c("N", "U")),
            list(stringsAsFactors = FALSE))
    newargs <- do.call(expand.grid, .a)
    newargs[["Dim"]] <- .Dim[.i <- newargs[["Dim"]]]
    newargs[["Dimnames"]] <- lapply(.Dim, mkDN)[.i]
    newargs[["x"]] <- lapply(.Dim, mkX)[.i]

    if (is.sy) {
        iD <- Matrix:::indDiag
        if (extends(Class, "corMatrix"))
            newargs[["x"]] <-
                lapply(newargs[["x"]], replace, iD(n), 1)
        else if (extends(Class, "copMatrix"))
            newargs[["x"]] <-
                Map(function(x, upper) replace(x, iD(n, upper, TRUE), 1),
                    newargs[["x"]], newargs[["uplo"]] == "U")
    }

    ## Test the matrices generated by each set of arguments to 'new'
    all(unlist(.mapply(testDenseMatrix, newargs, list(Class = Class))))
}

testDenseMatrix <- function(Class, ...) {
    is.p  <- extends(Class, packedMatrix)
    is.ge <- extends(Class, generalMatrix)
    is.sy <- !is.ge && extends(Class, symmetricMatrix)
    is.tr <- !is.ge && !is.sy
    is.d  <- extends(Class, dMatrix)
    is.l  <- !is.d  && extends(Class, lMatrix)
    is.n  <- !is.d  && !is.l

    ## This class needs special care because it has an additional 'sd' slot
    is.cr <- is.sy &&
        (extends(Class, "corMatrix") || extends(Class, "copMatrix"))

    newargs <- list(Class = Class, ...)
    if (is.cr)
        newargs[["sd"]] <- rep.int(1, newargs[["Dim"]][2L])
    .M <- M <- do.call(new, newargs)

    m <- M@Dim[1L]
    n <- M@Dim[2L]
    r <- min(m, n)
    p0 <- (n * (n - 1L)) %/% 2L
    p1 <- (n * (n + 1L)) %/% 2L
    .ZERO  <- as.vector(0, typeof(M@x))
    .ONE  <- as.vector(1, typeof(M@x))
    .NA <- as.vector(NA, typeof(M@x))
    loup <- if (is.ge) NA_character_ else if (M@uplo == "U") "L" else "U"

    ## For conveniently getting and setting (non)trivial triangles
    ## of _traditional_ symmetric or triangular matrices
    if (!is.ge) {
        if (M@uplo == "U") {
            tri1 <- U; `tri1<-` <- `U<-`
            tri0 <- L; `tri0<-` <- `L<-`
        } else {
            tri1 <- L; `tri1<-` <- `L<-`
            tri0 <- U; `tri0<-` <- `U<-`
        }
    }

    .m <- m1 <- m2 <-
        if (is.p)
            `tri1<-`(array(.ZERO, dim = M@Dim, dimnames = M@Dimnames),
                     diag = TRUE, M@x)
        else
            array(M@x, dim = M@Dim, dimnames = M@Dimnames)

    diag(.M) <- diag(.m) <-
        if (is.d)
            rlnorm(r)
        else
            sample(c(.NA, .ZERO, .ONE), r, TRUE)

    if (is.n && anyNA(m2))
        m2[is.na(m2)] <- TRUE
    if (is.sy) {
        tri0(m2, diag = TRUE) <- tri0(t(m2), diag = TRUE)
        dimnames(m2) <- Matrix:::symDN(M@Dimnames)
    }
    if (is.tr && M@diag == "U")
        diag(m2) <- .ONE

    pM <-
        if (is.p)
            M
        else if (!is.ge)
            do.call(new, replace(newargs, c("Class", "x"),
                                 list(packedClass(Class),
                                      tri1(m1, diag = TRUE))))

    upM <-
        if (!is.p)
            M
        else
            do.call(new, replace(newargs, c("Class", "x"),
                                 list(unpackedClass(Class), as.vector(m1))))

    tM <- do.call(new, replace(newargs,
                               c("Dim", "Dimnames", "x",
                                 if (!is.ge) "uplo"),
                               c(list(M@Dim[2:1],
                                      M@Dimnames[if (is.sy) 1:2 else 2:1],
                                      if (is.p)
                                          tri0(t(m1), diag = TRUE)
                                      else
                                          as.vector(t(m1))),
                                 if (!is.ge) list(loup))))

    dM <- do.call(new, replace(newargs[names(newargs) != "sd"],
                               c("Class", "x", if (is.tr) "diag"),
                               c(list(...Class(Class),
                                      if (is.p)
                                          tri1(.m, diag = TRUE)
                                      else
                                          as.vector(.m)),
                                 if (is.tr) list("N"))))

    stopifnot(is.ge || identical(pack(M), pM),
              identical(unpack(M), upM),
              identical(t(M), tM),
              identical(diag(M, names = FALSE), diag(m2, names = FALSE)),
              identical(diag(M, names = TRUE),  diag(m2, names = TRUE)),
              identical(.M, dM))

    if (is.ge) {
        if (m == n) {
            ## Not symmetric and not triangular
            U(m2) <- if (is.d)  rlnorm(p0) else rep_len(c(.NA, .ZERO, .ONE), p0)
            L(m2) <- if (is.d) -rlnorm(p0) else rep_len(c(.ZERO, .ONE, .NA), p0)
            M@x <- as.vector(m2)
            stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE) == (n <= 1L),
                      isTriangular(M, upper = NA) == (n <= 1L),
                      isDiagonal(M) == (n <= 1L))

            ## Symmetric but not triangular
            L(m2) <- L(t(m2))
            M@x <- as.vector(m2)
            stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE),
                      isTriangular(M, upper = NA) == (n <= 1L),
                      isDiagonal(M) == (n <= 1L))

            ## Not symmetric but triangular
            L(m2) <- .ZERO
            M@x <- as.vector(m2)
            stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE) == (n <= 1L),
                      identical(isTriangular(M, upper = TRUE), TRUE),
                      identical(isTriangular(M, upper = FALSE), FALSE),
                      identical(isTriangular(M, upper = NA),
                                `attr<-`(TRUE, "kind", "U")),
                      isDiagonal(M) == (n <= 1L))

            ## Symmetric and triangular
            U(m2) <- .ZERO
            M@x <- as.vector(m2)
            stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE),
                      identical(isTriangular(M, upper = TRUE), TRUE),
                      identical(isTriangular(M, upper = FALSE), TRUE),
                      identical(isTriangular(M, upper = NA),
                                `attr<-`(TRUE, "kind", "U")),
                      isDiagonal(M))
        } else {
            ## Non-square ... _never_ symmetric, triangular, or diagonal
            stopifnot(!isSymmetric(M, tol = 0, checkDN = FALSE),
                      !isTriangular(M, upper = NA),
                      !isDiagonal(M))
        }

    } else if (is.sy) {
        ## Not triangular
        tri1(m2) <- if (is.d) rlnorm(p0) else rep_len(c(.NA, .ZERO, .ONE), p0)
        M@x <- if (is.p) tri1(m2, diag = TRUE) else as.vector(m2)
        stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE),
                  isTriangular(M, upper = NA) == (n <= 1L),
                  isDiagonal(M) == (n <= 1L))

        ## Triangular
        tri1(m2) <- .ZERO
        M@x <- if (is.p) tri1(m2, diag = TRUE) else as.vector(m2)
        stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE),
                  identical(isTriangular(M, upper = TRUE), TRUE),
                  identical(isTriangular(M, upper = FALSE), TRUE),
                  identical(isTriangular(M, upper = NA),
                            `attr<-`(TRUE, "kind", M@uplo)),
                  isDiagonal(M))
    } else {
        ## Not symmetric
        tri1(m2) <- if (is.d) rlnorm(p0) else rep_len(c(.NA, .ZERO, .ONE), p0)
        M@x <- if (is.p) tri1(m2, diag = TRUE) else as.vector(m2)
        stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE) == (n <= 1L),
                  identical(isTriangular(M, upper = M@uplo == "U"), TRUE),
                  identical(isTriangular(M, upper = M@uplo != "U"), n <= 1L),
                  identical(isTriangular(M, upper = NA),
                            `attr<-`(TRUE, "kind", M@uplo)),
                  isDiagonal(M) == (n <= 1L))

        ## Symmetric
        tri1(m2) <- .ZERO
        M@x <- if (is.p) tri1(m2, diag = TRUE) else as.vector(m2)
        stopifnot(isSymmetric(M, tol = 0, checkDN = FALSE),
                  identical(isTriangular(M, upper = M@uplo == "U"), TRUE),
                  identical(isTriangular(M, upper = M@uplo != "U"), TRUE),
                  identical(isTriangular(M, upper = NA),
                            `attr<-`(TRUE, "kind", M@uplo)),
                  isDiagonal(M))
    }

    TRUE
}

.dense.subclasses <- c(names(getClassDef("packedMatrix")@subclasses),
                       names(getClassDef("unpackedMatrix")@subclasses))
## For now:
.dense.subclasses <- grep("^[nld]", .dense.subclasses, value = TRUE)
stopifnot(all(vapply(.dense.subclasses, testDenseClass, NA, n = 4L)))


## diag(<non-square .geMatrix>, names = TRUE) preserves names
## if head of longer character vector matches shorter character vector
n <- 4L
m <- array(rnorm(n * (n + 1L)), dim = c(n, n + 1L))
M <- new("dgeMatrix", x = as.vector(m), Dim = dim(m))
rn <- letters[seq_len(n)]
cn <- letters[seq_len(n + 1L)]
ldn <- list(list(rn, cn), list(rn, replace(cn, 1L, "")))
for (dn in ldn)
    stopifnot(identical(diag(`slot<-`(M, "Dimnames", TRUE, dn), names = TRUE),
                        diag(`dimnames<-`(m, dn), names = TRUE)))

cat("Time elapsed:", proc.time(), "\n") # "stats"
