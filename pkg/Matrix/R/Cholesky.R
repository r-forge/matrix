## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", c(x = "generalMatrix"),
          function(x, uplo = "U", ...) {
              ch <- chol(forceSymmetric(x, uplo), ...)
              ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("chol", c(x = "symmetricMatrix"),
          function(x, ...)
              chol(.M2kind(x, ","), ...))

setMethod("chol", c(x = "triangularMatrix"),
          function(x, uplo = "U", ...) {
              if(identical(uplo, x@uplo)) {
                  ch <- chol(forceSymmetric(x, uplo), ...)
                  ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
                  ch
              } else chol(forceDiagonal(x, x@diag), ...)
          })

setMethod("chol", c(x = "diagonalMatrix"),
          function(x, ...)
              chol(.M2kind(x, ","), ...))

setMethod("chol", c(x = "dsyMatrix"),
          function(x, pivot = FALSE, tol = -1, ...) {
              ch <- as(Cholesky(x, perm = pivot, tol = tol), "dtrMatrix")
              ch@Dimnames <- dimnames(x)
              if(ch@uplo != "U") t(ch) else ch
          })

setMethod("chol", c(x = "dspMatrix"),
          function(x, ...) {
              ch <- as(Cholesky(x), "dtpMatrix")
              ch@Dimnames <- dimnames(x)
              if(ch@uplo != "U") t(ch) else ch
          })

for(.cl in paste0("ds", c("C", "R", "T"), "Matrix"))
setMethod("chol", c(x = .cl),
          function(x, pivot = FALSE, ...) {
              ch <- t(as(Cholesky(x, perm = pivot, LDL = FALSE, super = FALSE),
                         "dtCMatrix")) # FIXME? give dtRMatrix, dtTMatrix?
              ch@Dimnames <- dimnames(x)
              ch
          })
rm(.cl)

setMethod("chol", c(x = "ddiMatrix"),
          function(x, ...) {
              if(length(y <- x@x)) {
                  if(is.na(min.y <- min(y)) || min.y < 0)
                      stop(gettextf("%1$s(%2$s) is undefined: '%2$s' is not positive semidefinite",
                                    "chol", "x"),
                           domain = NA)
                  x@x <- sqrt(y)
              }
              x
          })


## METHODS FOR GENERIC: Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Cholesky", c(A = "generalMatrix"),
          function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              ch@Dimnames <- A@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("Cholesky", c(A = "symmetricMatrix"),
          function(A, ...)
              Cholesky(.M2kind(A, ","), ...))

setMethod("Cholesky", c(A = "triangularMatrix"),
          function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              ch@Dimnames <- A@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("Cholesky", c(A = "diagonalMatrix"),
          function(A, ...)
              Cholesky(.M2kind(A, ","), ...))

setMethod("Cholesky", c(A = "dsyMatrix"),
          function(A, perm = TRUE, tol = -1, ...)
              .Call(poMatrix_trf, A, if(perm) 1L else 2L, perm, tol))

setMethod("Cholesky", c(A = "dspMatrix"),
          function(A, ...)
              .Call(ppMatrix_trf, A, 2L))

setMethod("Cholesky", c(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = !super, super = FALSE,
                   Imult = 0, ...)
              .Call(pCMatrix_trf, A, 2L, perm, !LDL, super, Imult))

setMethod("Cholesky", c(A = "dsRMatrix"),
          function(A, ...)
              Cholesky(.tCRT(A), ...))

setMethod("Cholesky", c(A = "dsTMatrix"),
          function(A, ...)
              Cholesky(.M2C(A), ...))

setMethod("Cholesky", c(A = "ddiMatrix"),
          function(A, ...) {
              if(length(y <- A@x) && (is.na(min.y <- min(y)) || min.y < 0))
                  stop(gettextf("%1$s(%2$s) is undefined: '%2$s' is not positive semidefinite",
                                "Cholesky", "x"),
                       domain = NA)
              n <- (d <- A@Dim)[1L]
              r <- new("dsimplicialCholesky")
              r@Dim <- d
              r@Dimnames <- A@Dimnames
              r@minor <- n
              r@colcount <- r@nz <- rep.int(1L, n)
              r@`next` <- c(seq_len(n), -1L, 0L)
              r@ prev  <- c(n + 1L, s, -1L)
              r@p <- 0:n
              r@i <- s <- seq.int(0L, length.out = n)
              r@x <- if(length(y)) y else rep.int(1, n)
              r@ordering <- 0L
              r@is_ll <- FALSE
              r@is_monotonic <- TRUE
              r
          })

setMethod("Cholesky", c(A = "matrix"),
          function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              if(!is.null(dn <- dimnames(A)))
                  ch@Dimnames <- dn # restore asymmetric 'Dimnames'
              ch
          })


## METHODS FOR GENERIC: chol2inv
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol2inv", c(x = "generalMatrix"),
          function(x, uplo = "U", ...) {
              d <- x@Dim
              if(d[1L] != d[2L])
                  stop("matrix is not square")
              chol2inv((if(  uplo == "U") triu else tril)(x), ...)
          })

setMethod("chol2inv", c(x = "symmetricMatrix"),
          function(x, ...)
              chol2inv((if(x@uplo == "U") triu else tril)(x), ...))

setMethod("chol2inv", c(x = "triangularMatrix"),
          function(x, ...)
              chol2inv(.M2kind(x, ","), ...))

setMethod("chol2inv", c(x = "diagonalMatrix"),
          function(x, ...)
              chol2inv(.M2kind(x, ","), ...))

for(.cl in paste0("dt", c("r", "p"), "Matrix"))
setMethod("chol2inv", c(x = .cl),
          function(x, ...) {
              if(x@diag != "N")
                  x <- ..diagU2N(x)
              r <- .Call(denseCholesky_solve, x, NULL)
              i <- if(x@uplo == "U") 2L else 1L
              r@Dimnames <- x@Dimnames[c(i, i)]
              r
          })
rm(.cl)

for(.cl in paste0("dt", c("C", "R", "T"), "Matrix"))
setMethod("chol2inv", c(x = .cl),
          function(x, ...)
              (if(x@uplo == "U") tcrossprod else crossprod)(solve(x)))
rm(.cl)

## 'uplo' can affect the 'Dimnames' of the result here :
setMethod("chol2inv", c(x = "ddiMatrix"),
          function(x, uplo = "U", ...)
              (if(  uplo == "U") tcrossprod else crossprod)(solve(x)))


## METHODS FOR CLASS: denseCholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("denseCholesky", "dtrMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if(!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if(!packed) to else as(to, "unpackedMatrix")
      })

setAs("denseCholesky", "dtpMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if(!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if(!packed) as(to, "packedMatrix") else to
      })

setMethod("diag", c(x = "denseCholesky"),
          function(x = 1, nrow, ncol, names = TRUE) {
              packed <- length(x@x) != prod(x@Dim)
              d <- diag(as(x, if(!packed) "dtrMatrix" else "dtpMatrix"),
                        names = FALSE)
              d * d
          })

.f1 <- function(x, which, ...) {
    d <- x@Dim
    switch(which,
           "P1" =, "P1." =
               {
                   r <- new("pMatrix")
                   r@Dim <- d
                   r@perm <- if(length(x@perm)) x@perm else seq_len(d[1L])
                   if(which == "P1.")
                       r@margin <- 2L
                   r
               },
           "L" =, "L." =, "L1" =, "L1." =
               {
                   packed <- length(x@x) != prod(x@Dim)
                   r <- as(x, if(!packed) "dtrMatrix" else "dtpMatrix")
                   uplo <- x@uplo
                   if(which == "L1" || which == "L1.") {
                       r.ii <- diag(r, names = FALSE)
                       r@x <- r@x /
                           if(uplo == "U") {
                               if(!packed)
                                   r.ii
                               else r.ii[sequence.default(seq_len(d[1L]))]
                           } else {
                               if(!packed)
                                   rep(r.ii, each = d[1L])
                               else rep.int(r.ii, seq.int(to = 1L, by = -1L, length.out = d[1L]))
                           }
                       r@diag <- "U"
                   }
                   if((which == "L." || which == "L1.") == (uplo == "U"))
                       r
                   else t(r)
               },
           "D" =
               {
                   r <- new("ddiMatrix")
                   r@Dim <- d
                   r@x <- diag(x, names = FALSE)
                   r
               },
           stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%3$s\", \"%3$s.\", \"%3$s1\", \"%3$s1.\", or \"%4$s\"",
                         "which", "P", "L", "D"),
                domain = NA))
}

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
.f2 <- function(x, LDL = TRUE, ...) {
    packed <- length(x@x) != prod(x@Dim)
    d <- x@Dim
    dn <- x@Dimnames
    uplo <- x@uplo
    perm <- x@perm

    P <- new("pMatrix")
    P@Dim <- d
    P@Dimnames <- c(list(NULL), dn[2L])
    P@margin <- 2L
    P@perm <- if(length(perm)) invertPerm(perm) else seq_len(d[1L])

    P. <- P
    P.@Dimnames <- c(dn[1L], list(NULL))
    P.@margin <- 1L

    X <- as(x, if(!packed) "dtrMatrix" else "dtpMatrix")
    if(LDL) {
        L.ii <- diag(X, names = FALSE)
        X@x <- X@x /
            if(uplo == "U") {
                if(!packed)
                    L.ii
                else L.ii[sequence.default(seq_len(d[1L]))]
            } else {
                if(!packed)
                    rep(L.ii, each = d[1L])
                else rep.int(L.ii, seq.int(to = 1L, by = -1L, length.out = d[1L]))
            }
        X@diag <- "U"
    }
    L  <- if(uplo == "U") t(X) else   X
    L. <- if(uplo == "U")   X  else t(X)
    if(LDL) {
        D <- new("ddiMatrix")
        D@Dim <- d
        D@x <- L.ii * L.ii
        list(P1. = P., L1 = L, D = D, L1. = L., P1 = P)
    } else list(P1. = P., L = L, L. = L., P1 = P)
}


setMethod("expand1", c(x = "denseCholesky"), .f1)
setMethod("expand2", c(x = "denseCholesky"), .f2)

rm(.f1, .f2)


## METHODS FOR CLASS: sparseCholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.CHF.is.perm <- function(x)
    x@ordering != 0L
.CHF.is.LDL <- function(x)
    .hasSlot(x, "is_ll") && !x@is_ll
.CHF.is.super <- function(x)
    .hasSlot(x, "super")

# Exported:
isLDL <- function(x) {
    if(is(x, "sparseCholesky"))
        .CHF.is.LDL(x)
    else stop(gettextf("'%s' does not inherit from virtual class %s",
                       "x", "sparseCholesky"),
              domain = NA)
}

## Exported:
.updateCHMfactor <- function(object, parent, mult = 0)
    .Call(sparseCholesky_update, object, parent, mult)

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

setMethod("diag", c(x = "sparseCholesky"),
          function(x = 1, nrow, ncol, names = TRUE)
              .Call(sparseCholesky_diag_get, x, TRUE))

setMethod("expand1", c(x = "simplicialCholesky"),
          function(x, which, ...) {
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- d <- x@Dim
                         r@perm <- if(length(x@perm)) x@perm + 1L else seq_len(d[1L])
                         if(which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "L" =, "L." =, "L1" =, "L1." = {
                         r <- as(x, "dtCMatrix")
                         LDL. <- .CHF.is.LDL(x)
                         if(which == "L1" || which == "L1.") {
                             if(!LDL.) {
                                 r.ii <- diag(r, names = FALSE)
                                 r.p <- r@p
                                 r@x <- r@x / rep.int(r.ii , r.p[-1L] - r.p[-length(r.p)])
                             }
                             diag(r) <- 1
                         } else if(LDL.) {
                             r.ii <- diag(r, names = FALSE)
                             r.p <- r@p
                             if(anyNA(r.ii))
                                 stop(gettextf("D[i,i] is NA, i=%d",
                                               which.max(is.na(r.ii))),
                                      domain = NA)
                             if(min(r.ii) < 0)
                                 stop(gettextf("D[i,i] is negative, i=%d",
                                               which.max(r.ii < 0)),
                                      domain = NA)
                             diag(r) <- 1
                             r@x <- r@x * rep.int(sqrt(r.ii), r.p[-1L] - r.p[-length(r.p)])
                         }
                         if(which == "L" || which == "L1")
                             r
                         else t(r)
                     },
                     "D" = {
                         r <- new("ddiMatrix")
                         r@Dim <- x@Dim
                         r@x <- diag(x, names = FALSE)
                         r
                     },
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%3$s\", \"%3$s.\", \"%3$s1\", \"%3$s1.\", or \"%4$s\"",
                                   "which", "P", "L", "D"),
                          domain = NA))
          })

setMethod("expand1", c(x = "supernodalCholesky"),
          function(x, which, ...) {
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- d <- x@Dim
                         r@perm <- if(length(x@perm)) x@perm + 1L else seq_len(d[1L])
                         if(which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "L" =, "L." =, "L1" =, "L1." = {
                         r <- as(x, "dgCMatrix")
                         if(which == "L1" || which == "L1.") {
                             r.ii <- diag(r, names = FALSE)
                             r.p <- r@p
                             r@x <- r@x / rep.int(r.ii, r.p[-1L] - r.p[-length(r.p)])
                             diag(r) <- 1
                         }
                         if(which == "L" || which == "L1")
                             r
                         else t(r)
                     },
                     "D" = {
                         r <- new("ddiMatrix")
                         r@Dim <- x@Dim
                         r@x <- diag(x, names = FALSE)
                         r
                     },
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%3$s\", \"%3$s.\", \"%3$s1\", \"%3$s1.\", or \"%4$s\"",
                                   "which", "P", "L", "D"),
                          domain = NA))
          })

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
setMethod("expand2", c(x = "simplicialCholesky"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              perm <- x@perm

              P <- new("pMatrix")
              P@Dim <- d
              P@Dimnames <- c(list(NULL), dn[2L])
              P@margin <- 2L
              P@perm <- if(length(perm))
                            invertPerm(perm, 0L, 1L)
                        else seq_len(d[1L])

              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L

              L <- as(x, "dtCMatrix")
              LDL. <- .CHF.is.LDL(x)
              if(!LDL && !LDL.)
                  return(list(P1. = P., L = L, L. = t(L), P1 = P))
              L.ii <- diag(L, names = FALSE)
              L.p <- L@p
              if(!LDL) {
                  if(anyNA(L.ii))
                      stop(gettextf("D[i,i] is NA, i=%d",
                                    which.max(is.na(L.ii))),
                           domain = NA)
                  if(min(L.ii) < 0)
                      stop(gettextf("D[i,i] is negative, i=%d",
                                    which.max(L.ii < 0)),
                           domain = NA)
                  diag(L) <- 1
                  L@x <- L@x * rep.int(sqrt(L.ii), L.p[-1L] - L.p[-length(L.p)])
                  return(list(P1. = P., L = L, L. = t(L), P1 = P))
              }
              D <- new("ddiMatrix")
              D@Dim <- d
              if(LDL.) {
                  diag(L) <- 1
                  D@x <- L.ii
              } else {
                  L@x <- L@x / rep.int(L.ii , L.p[-1L] - L.p[-length(L.p)])
                  D@x <- L.ii * L.ii
              }
              list(P1. = P., L1 = L, D = D, L1. = t(L), P1 = P)
          })

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
setMethod("expand2", c(x = "supernodalCholesky"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              perm <- x@perm

              P <- new("pMatrix")
              P@Dim <- d
              P@Dimnames <- c(list(NULL), dn[2L])
              P@margin <- 2L
              P@perm <- if(length(perm))
                            invertPerm(perm, 0L, 1L)
                        else seq_len(d[1L])

              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L

              L <- as(x, "dgCMatrix")
              if(!LDL)
                  return(list(P1. = P., L = L, L. = t(L), P1 = P))
              L.ii <- diag(L, names = FALSE)
              L.p <- L@p
              L@x <- L@x / rep.int(L.ii, L.p[-1L] - L.p[-length(L.p)])
              diag(L) <- 1
              D <- new("ddiMatrix")
              D@Dim <- d
              D@x <- L.ii * L.ii
              list(P1. = P., L1 = L, D = D, L1. = t(L), P1 = P)
          })

## returning list(P, L), where A = P' L L' P
## MJ: for backwards compatibility
setMethod("expand", c(x = "sparseCholesky"),
          function(x, ...)
              list(P = expand1(x, "P1"), L = expand1(x, "L")))

setMethod("update", c(object = "sparseCholesky"),
          function(object, parent, mult = 0, ...) {
              parent <- .M2kind(.M2C(parent), ",")
              if((shape <- .M.shape(parent)) != "s") {
                  Matrix.message(gettextf("'%1$s' is not formally symmetric; factorizing tcrossprod(%1$s)",
                                          "parent"),
                                 domain = NA)
                  if(shape == "t" && parent@diag != "N")
                      parent <- ..diagU2N(parent)
              }
              .Call(sparseCholesky_update, object, parent, mult)
          })

setMethod("updown",
          c(update = "character", C = "ANY", L = "ANY"),
          function(update, C, L)
              updown(identical(update, "+"), C, L))

setMethod("updown",
          c(update = "logical", C = "Matrix", L = "sparseCholesky"),
          function(update, C, L)
              updown(update, .M2kind(.M2C(C), ","), L))

setMethod("updown",
          c(update = "logical", C = "matrix", L = "sparseCholesky"),
          function(update, C, L)
              updown(update, .m2sparse(C, ",gC"), L))

for(.cl in c("dgCMatrix", "dsCMatrix"))
setMethod("updown",
          c(update = "logical", C = .cl, L = "sparseCholesky"),
          function(update, C, L) {
              if(length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .Call(sparseCholesky_updown, L, C, update)
          })
rm(.cl)

setMethod("updown",
          c(update = "logical", C = "dtCMatrix", L = "sparseCholesky"),
          function(update, C, L) {
              if(C@diag != "N")
                  C <- ..diagU2N(C)
              if(length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .Call(sparseCholesky_updown, L, C, update)
          })
