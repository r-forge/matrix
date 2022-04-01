library(Matrix)

if (interactive()) {
    options(Matrix.verbose = TRUE, warn = 1, error = recover)
} else {
    options(Matrix.verbose = TRUE, warn = 1)
}

## For getting and setting '[dD]imnames' on '[mM]atrix'
DN <- function(x) {
    if (is(x, "Matrix")) {
        x@Dimnames
    } else {
        dimnames(x)
    }
}
`DN<-` <- function(x, value) {
    if (is(x, "Matrix")) {
        x@Dimnames <- value
    } else {
        dimnames(x) <- value
    }
    x
}

## SDN1(dn) is documented to behave as SDN2(dn, NULL)
SDN1 <- Matrix:::symmDN
SDN2 <- function(dn, uplo = NULL) {
    J <-
        if (is.null(uplo)) {
            if (!is.null(dn[[1L]]) && is.null(dn[[2L]])) 1L else 2L
        } else {
            if (uplo == "U") 2L else 1L
        }
    rep(dn[J], 2L)
}

## Various possible (a)symmetries of 'Dimnames'
n <- 4L
rn <- letters[seq_len(n)]
cn <- LETTERS[seq_len(n)]
ldn <- list(list(rn, cn),
            list(rn, NULL),
            list(NULL, cn),
            list(NULL, NULL),
            list(x = rn, y = cn),
            list(x = rn, y = NULL),
            list(x = NULL, y = cn),
            list(x = NULL, y = NULL))

## 'matrix' and _most_ 'd..Matrix' ...
## zero matrices are fine for the purpose of testing handling of 'Dimnames'
lM <- c(list(matrix(0, n, n),
             new("ddiMatrix", x = double(n), Dim = c(n, n)),
             new("dgeMatrix", x = double(n * n), Dim = c(n, n))),
        .mapply(new,
                expand.grid(Class = c("dsyMatrix", "dtrMatrix"),
                            uplo = c("U", "L"),
                            stringsAsFactors = FALSE),
                list(x = double(n * n), Dim = c(n, n))),
        .mapply(new,
                expand.grid(Class = c("dspMatrix", "dtpMatrix"),
                            uplo = c("U", "L"),
                            stringsAsFactors = FALSE),
                list(x = double((n * (n + 1L)) %/% 2L), Dim = c(n, n))),
        list(new("dgCMatrix", x = double(0L), Dim = c(n, n),
                 i = integer(0L), p = rep.int(0L, n + 1L))),
        .mapply(new,
                expand.grid(Class = c("dsCMatrix", "dtCMatrix"),
                            uplo = c("U", "L"),
                            stringsAsFactors = FALSE),
                list(x = double(0L), Dim = c(n, n),
                     i = integer(0L), p = rep.int(0L, n + 1L))))

for (dn in ldn) {
    stopifnot(identical(sdn <- SDN1(dn), SDN2(dn)))
    for (M in lM) {
        DN(M) <- dn
        if (is.s <- is(M, "symmetricMatrix")) {
            ## 'dimnames' should symmetrize
            stopifnot(identical(dimnames(M), sdn))
        }

        if (is.s && !identical(dn[1L], dn[2L])) {
            ## Methods for 'symmetricMatrix' assume symmetric 'Dimnames'
            ## for efficiency ... should they?
            next
        }

        ## FIXME: Conditions for symmetrization are inconsistent
        ## between 'Csparse_general_to_symmetric' from ../src/Csparse.c
        ## and 'dense_to_symmetric' from ../src/dense.c
        fsdn <-
            if (is.s) {
                sdn
            } else if (is(M, "sparseMatrix") &&
                       (is.null(dn[[1L]]) || is.null(dn[[2L]])) &&
                       is.null(names(dn))) {
                dn
            } else if (is(M, "triangularMatrix")) {
                SDN2(dn, M@uplo)
            } else {
                SDN2(dn, "U")
            }
        stopifnot(identical(DN(forceSymmetric(M)), fsdn))

        if (is(M, "diagonalMatrix")) {
            ## Methods for 'diagonalMatrix' assume symmetric 'Dimnames'
            ## for efficiency (or maybe by accident?) ... should they?
            next
        }

        stopifnot(identical(DN(symmpart(M)), sdn),
                  identical(DN(skewpart(M)), sdn))
        ## others?
    }
}

cat("Time elapsed:", proc.time(), "\n") # "stats"
