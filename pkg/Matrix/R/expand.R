## METHODS FOR GENERIC: expand
## NB: only for backwards compatibility; prefer expand1, expand2 (below)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(L, U, P), where A = P L U
setMethod("expand", c(x = "denseLU"),
          function(x, ...) {
              r <- expand2(x)[c(2L, 3L, 1L)]
              names(r) <- c("L", "U", "P")
              NN <- list(NULL, NULL)
              for (i in 1:3)
                  r[[i]]@Dimnames <- NN
              r
          })

## returning list(P, L, U, Q), where A = P' L U Q
setMethod("expand", c(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              p <- x@p + 1L
              q <- x@q + 1L
              if (length(p) && !is.null(rn <- dn[[1L]]))
                  dn[[1L]] <- rn[p]
              if (length(q) && !is.null(cn <- dn[[2L]]))
                  dn[[2L]] <- cn[q]
              P <- new("pMatrix")
              P@Dim <- d
              P@perm <- p
              Q <- new("pMatrix")
              Q@Dim <- d
              Q@perm <- if (length(q)) q else seq_len(d[1L])
              L <- x@L
              L@Dimnames <- c(dn[1L], list(NULL))
              U <- x@U
              U@Dimnames <- c(list(NULL), dn[2L])
              list(P = P, L = L, U = U, Q = Q)
          })

## returning list(P, L), where A = P' L L' P
setMethod("expand", c(x = "sparseCholesky"),
          function(x, ...)
              list(P = expand1(x, "P1"), L = expand1(x, "L")))


## METHODS FOR GENERIC: expand1
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("expand1", c(x = "denseSchur"),
          function(x, which, ...) {
              n <- (d <- x@Dim)[1L]
              switch(which,
                     "T" =
                         {
                             y <- x@x
                             cl <-
                                 if (n > 0L && length(y) == 0L) {
                                     y <- x@values
                                     "ddiMatrix"
                                 }
                                 else if (is.complex(y))
                                     "ztrMatrix"
                                 else if (is.double(x@values))
                                     "dtrMatrix"
                                 else "dgeMatrix"
                             T <- new(cl)
                             T@Dim <- d
                             T@x <- y
                             T
                         },
                     "Q" =, "Q." =
                         {
                             y <- x@vectors
                             cl <-
                                 if (n > 0L && length(y) == 0L)
                                     stop("missing requisite Schur vectors")
                                 else if (is.complex(y))
                                     "zgeMatrix"
                                 else "dgeMatrix"
                             Q <- new(cl)
                             Q@Dim <- d
                             Q@x <- y
                             switch(which, "Q" = Q, "Q." = ct(Q))
                         },
                     stop(gettextf("'%s' is not \"%s\", \"%s\", or \"%s\"",
                                   "which", "Q", "T", "Q."),
                          domain = NA))
          })

.expand.mkL <-
function(x, m, n) {
    if (m == n)
        x
    else if (m < n)
        x[seq_len(prod(m, m))]
    else {
        x[seq.int(from = 1L, by = m + 1, length.out = n)] <- 1
        if (n > 1L) {
            nvec <- seq_len(n - 1L)
            from <- seq.int(from = m + 1, by = m, length.out = n - 1L)
            x[sequence.default(nvec, from)] <- 0
        }
        x
    }
}

.expand.mkU <-
function(x, m, n) {
    if (m == n)
        x
    else if (m > n) {
        nvec <- rep.int(n, n)
        from <- seq.int(from = 1L, by = m, length.out = n)
        x[sequence.default(nvec, from)]
    }
    else {
        if (m > 1L) {
            nvec <- seq.int(to = 1L, by = -1L, length.out = m - 1L)
            from <- seq.int(from = 2L, by = m + 1, length.out = m - 1L)
            x[sequence.default(nvec, from)] <- 0
        }
        x
    }
}

setMethod("expand1", c(x = "denseLU"),
          function(x, which, ...) {
              d <- x@Dim
              m <- d[1L]
              n <- d[2L]
              switch(which,
                     "P1" =, "P1." =
                         {
                             P1 <- new("pMatrix")
                             P1@Dim <- c(m, m)
                             P1@perm <- asPerm(x@perm, n = m)
                             if (which == "P1.")
                                 P1@margin <- 2L
                             P1
                         },
                     "L" =
                         {
                             z <- is.complex(y <- x@x)
                             if (m <= n) {
                                 L <- new(if (z) "ztrMatrix" else "dtrMatrix")
                                 L@Dim <- c(m, m)
                                 L@uplo <- "L"
                                 L@diag <- "U"
                             } else {
                                 L <- new(if (z) "zgeMatrix" else "dgeMatrix")
                                 L@Dim <- d
                             }
                             L@x <- .expand.mkL(y, m, n)
                             L
                         },
                     "U" =
                         {
                             z <- is.complex(y <- x@x)
                             if (m >= n) {
                                 U <- new(if (z) "ztrMatrix" else "dtrMatrix")
                                 U@Dim <- c(n, n)
                             } else {
                                 U <- new(if (z) "zgeMatrix" else "dgeMatrix")
                                 U@Dim <- d
                             }
                             U@x <- .expand.mkU(y, m, n)
                             U
                         },
                     stop(gettextf("'%s' is not \"%s\", \"%s\", \"%s\", or \"%s\"",
                                   "which", "P1", "P1.", "L", "U"),
                          domain = NA))
          })

setMethod("expand1", c(x = "denseBunchKaufman"),
          function(x, which, ...) {
              l <- .Call(denseBunchKaufman_expand, x)
              b <- length(l) - 1L
              switch(which,
                     "DU" =, "DL" =
                         {
                             if (!endsWith(which, x@uplo))
                                 stop(gettextf("%s=\"%s\" invalid for %s=\"%s\"",
                                               "which", which, "x@uplo", x@uplo),
                                      domain = NA)
                             l[[b + 1L]]
                         },
                     "U" =, "U." =, "L" =, "L." =
                         {
                             if (!startsWith(which, x@uplo))
                                 stop(gettextf("%s=\"%s\" invalid for %s=\"%s\"",
                                               "which", which, "x@uplo", x@uplo),
                                      domain = NA)
                             z <- is.complex(x@x)
                             if (b > 0L) {
                                 L <- l[[b]]
                                 if (b > 1L)
                                     for (i in (b - 1L):1L)
                                         L <- l[[i]] %*% L
                                 if (!endsWith(which, "."))
                                     L
                                 else if (!z || x@trans == "C")
                                     ct(L)
                                 else t(L)
                             } else {
                                 D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                                 D@Dim <- x@Dim
                                 D@diag <- "U"
                                 D
                             }
                         },
                     stop(gettextf("'%s' is not \"%s\", \"%s\", or \"%s\"",
                                   "which", x@uplo, paste0("D", x@uplo), paste0(x@uplo, ".")),
                          domain = NA))
          })

setMethod("expand1", c(x = "denseCholesky"),
          function(x, which, ...) {
              d <- x@Dim
              switch(which,
                     "P1" =, "P1." =
                         {
                             P1 <- new("pMatrix")
                             P1@Dim <- d
                             P1@perm <- if (length(perm <- x@perm))
                                            perm
                                        else seq_len(d[1L])
                             if (which == "P1.")
                                 P1@margin <- 2L
                             P1
                         },
                     "L" =, "L." =, "L1" =, "L1." =
                         {
                             z <- is.complex(y <- x@x)
                             packed <- length(y) != prod(d)
                             L <- new(paste0(if (z) "z" else "d", if (!packed) "trMatrix" else "tpMatrix"))
                             L@Dim <- d
                             L@uplo <- x@uplo
                             L@x <- y
                             if (which == "L1" || which == "L1.") {
                                 L.jj <- diag(L, names = FALSE)
                                 L@x <- y /
                                     if (x@uplo == "U") {
                                         if (!packed)
                                             L.jj
                                         else L.jj[sequence.default(seq_len(d[1L]))]
                                     } else {
                                         if (!packed)
                                             rep(L.jj, each = d[1L])
                                         else rep.int(L.jj, seq.int(to = 1L, by = -1L, length.out = d[1L]))
                                     }
                                 L@diag <- "U"
                             }
                             if ((which == "L" || which == "L1") == (x@uplo != "U"))
                                 L
                             else ct(L)
                         },
                     "D" =
                         {
                             z <- is.complex(y <- diag(x, names = FALSE))
                             D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                             D@Dim <- d
                             D@x <- y
                             D
                         },
                     stop(gettextf("'%s' is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\"",
                                   "which", "P1", "P1.", "L", "L.", "L1", "L1.", "D"),
                          domain = NA))
          })

setMethod("expand1", c(x = "sparseQR"),
          function(x, which, ...) {
              .qr.rank.def.warn(x)
              R <- x@R
              d <- R@Dim
              m <- d[1L]
              n <- d[2L]
              switch(which,
                     "P1" =, "P1." =
                         {
                             P1 <- new("pMatrix")
                             P1@Dim <- c(m, m)
                             P1@perm <- x@p + 1L
                             if (which == "P1.")
                                 P1@margin <- 2L
                             P1
                         },
                     "P2" =, "P2." =
                         {
                             P2 <- new("pMatrix")
                             P2@Dim <- c(n, n)
                             P2@perm <- if (length(perm <- x@q))
                                            perm + 1L
                                        else seq_len(n)
                             if (which == "P2")
                                 P2@margin <- 2L
                             P2
                         },
                     "Q"  = .Call(sparseQR_matmult, x, NULL, 6L,  TRUE, NULL),
                     "Q1" = .Call(sparseQR_matmult, x, NULL, 6L, FALSE, NULL),
                     "R"  = R,
                     "R1" = triu(if (m == n) R else R[seq_len(n), , drop = FALSE]),
                     stop(gettextf("'%s' is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\"",
                                   "which", "P1", "P1.", "P2", "P2.", "Q", "Q1", "R", "R1"),
                          domain = NA))
          })

setMethod("expand1", c(x = "sparseLU"),
          function(x, which, ...) {
              d <- x@Dim
              switch(which,
                     "P1" =, "P1." =
                         {
                             P1 <- new("pMatrix")
                             P1@Dim <- d
                             P1@perm <- x@p + 1L
                             if (which == "P1.")
                                 P1@margin <- 2L
                             P1
                         },
                     "P2" =, "P2." =
                         {
                             P2 <- new("pMatrix")
                             P2@Dim <- d
                             P2@perm <- if (length(perm <- x@q))
                                            perm + 1L
                                        else seq_len(d[1L])
                             if (which == "P2")
                                 P2@margin <- 2L
                             P2
                         },
                     "L" = x@L,
                     "U" = x@U,
                     stop(gettextf("'%s' is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\"",
                                   "which", "P1", "P1.", "P2", "P2.", "L", "U"),
                          domain = NA))
          })

setMethod("expand1", c(x = "simplicialCholesky"),
          function(x, which, ...) {
              d <- x@Dim
              switch(which,
                     "P1" =, "P1." =
                         {
                             P1 <- new("pMatrix")
                             P1@Dim <- d
                             P1@perm <- if (length(perm <- x@perm))
                                            perm + 1L
                                        else seq_len(d[1L])
                             if (which == "P1.")
                                 P1@margin <- 2L
                             P1
                         },
                     "L" =, "L." =, "L1" =, "L1." =
                         {
                             nz <- x@nz

                             z <- is.complex(y <- x@x)
                             L <- new(if (z) "ztCMatrix" else "dtCMatrix")
                             L@Dim <- d
                             L@uplo <- "L"
                             L@p <- L.p <- c(0L, cumsum(nz))
                             k <- sequence.default(nz, x@p[seq_along(nz)] + 1L)
                             L@i <- x@i[k]
                             L@x <- y <- y[k]
                             if (which == "L1" || which == "L1.") {
                                 if (x@is_ll) {
                                 L.jj <- diag(L, names = FALSE)
                                 L@x <- y / rep.int(L.jj, L.p[-1L] - L.p[-length(L.p)])
                                 }
                                 diag(L) <- 1
                             } else if (!x@is_ll) {
                                 L.jj <- diag(L, names = FALSE)
                                 if (is.na(tmp <- min(0, L.jj)) || tmp < 0)
                                     stop(gettextf("D[j,j] is NaN or negative, j=%d",
                                                   which.max(is.na(L.jj) | L.jj < 0)),
                                          domain = NA)
                                 diag(L) <- 1
                                 L@x <- L@x * rep.int(sqrt(L.jj), L.p[-1L] - L.p[-length(L.p)])
                             }
                             if (which == "L" || which == "L1")
                                 L
                             else ct(L)
                         },
                     "D" =
                         {
                             z <- is.complex(y <- diag(x, names = FALSE))
                             D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                             D@Dim <- d
                             D@x <- y
                             D
                         },
                     stop(gettextf("'%s' is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\"",
                                   "which", "P1", "P1.", "L", "L.", "L1", "L1.", "D"),
                          domain = NA))
          })

setMethod("expand1", c(x = "supernodalCholesky"),
          function(x, which, ...) {
              d <- x@Dim
              switch(which,
                     "P1" =, "P1." =
                         {
                             P1 <- new("pMatrix")
                             P1@Dim <- d
                             P1@perm <- if (length(perm <- x@perm))
                                            perm + 1L
                                        else seq_len(d[1L])
                             if (which == "P1.")
                                 P1@margin <- 2L
                             P1
                         },
                     "L" =, "L." =, "L1" =, "L1." =
                         {
                             super <- x@super
                             pi <- x@pi
                             b <- length(super)

                             nr <- pi[-1L] - pi[-b]
                             nc <- super[-1L] - super[-b]
                             nz <- rep.int(nr, nc)

                             z <- is.complex(y <- x@x)
                             L <- new(if (z) "zgCMatrix" else "dgCMatrix")
                             L@Dim <- d
                             L@p <- L.p <- c(0L, cumsum(nz))
                             L@i <- x@s[sequence.default(nz, rep.int(pi[-b] + 1L, nc))]
                             L@x <- y
                             if (which == "L1" || which == "L1.") {
                                 L.jj <- diag(L, names = FALSE)
                                 L@x <- y / rep.int(L.jj, L.p[-1L] - L.p[-length(L.p)])
                             }
                             if (which == "L" || which == "L1")
                                 L
                             else ct(L)
                         },
                     "D" =
                         {
                             z <- is.complex(y <- diag(x, names = FALSE))
                             D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                             D@Dim <- d
                             D@x <- y
                             D
                         },
                     stop(gettextf("'%s' is not \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", or \"%s\"",
                                   "which", "P1", "P1.", "L", "L.", "L1", "L1.", "D"),
                          domain = NA))
           })


## METHODS FOR GENERIC: expand2
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(Q, T, Q'), where A = Q T Q'
setMethod("expand2", c(x = "denseSchur"),
          function(x, ...) {
              T <- expand1(x, "T")
              Q <- expand1(x, "Q")
              Q. <- ct(Q)
              dn <- x@Dimnames
              Q @Dimnames <- c(dn[1L], list(NULL))
              Q.@Dimnames <- c(list(NULL), dn[2L])
              list(Q = Q, T = T, Q. = Q.)
          })

## returning list(P1', L, U), where A = P1' L U
setMethod("expand2", c(x = "denseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              m <- d[1L]
              n <- d[2L]
              z <- is.complex(y <- x@x)

              P1. <- new("pMatrix")
              P1.@Dim <- c(m, m)
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@perm <- invertPerm(asPerm(x@perm, n = m))

              if (m <= n) {
                  L <- new(if (z) "ztrMatrix" else "dtrMatrix")
                  L@Dim <- c(m, m)
                  L@uplo <- "L"
                  L@diag <- "U"
              } else {
                  L <- new(if (z) "zgeMatrix" else "dgeMatrix")
                  L@Dim <- d
              }
              L@x <- .expand.mkL(y, m, n)

              if (m >= n) {
                  U <- new(if (z) "ztrMatrix" else "dtrMatrix")
                  U@Dim <- c(n, n)
              } else {
                  U <- new(if (z) "zgeMatrix" else "dgeMatrix")
                  U@Dim <- d
              }
              U@Dimnames <- c(list(NULL), dn[2L])
              U@x <- .expand.mkU(y, m, n)

              list(P1. = P1., L = L, U = U)
          })

## returning
## list(U, DU, U') where A = U DU U' and U = P[b] U[b] ... P[1] U[1]
## OR
## list(L, DL, L') where A = L DL L' and L = P[1] L[1] ... P[b] L[b]
setMethod("expand2", c(x = "denseBunchKaufman"),
          function(x, complete = FALSE, ...) {
              l <- .Call(denseBunchKaufman_expand, x)
              b <- length(l) - 1L
              z <- is.complex(x@x)
              if (complete) {
                  if (b > 0L)
                      l <- c(l, lapply(l[b:1L], if (!z || x@trans == "C") ct else t))
              } else {
                  if (b > 0L) {
                      L <- l[[b]]
                      if (b > 1L)
                          for (i in (b - 1L):1L)
                              L <- l[[i]] %*% L
                      l <- list(L, l[[b + 1L]], if (!z || x@trans == "C") ct(L) else t(L))
                  } else {
                      D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                      D@Dim <- x@Dim
                      D@diag <- "U"
                      l <- list(D, l[[1L]], D)
                  }
                  names(l) <- if (x@uplo == "U") c("U", "DU", "U.") else c("L", "DL", "L.")
              }
              dn <- x@Dimnames
              if (length(l) == 1L)
                  l[[1L]]@Dimnames <- dn
              else {
                  l[[1L]]@Dimnames <- c(dn[1L], list(NULL))
                  l[[length(l)]]@Dimnames <- c(list(NULL), dn[2L])
              }
              l
          })

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
setMethod("expand2", c(x = "denseCholesky"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              z <- is.complex(y <- x@x)
              packed <- length(y) != prod(d)

              P1 <- new("pMatrix")
              P1@Dim <- d
              P1@Dimnames <- c(list(NULL), dn[2L])
              P1@margin <- 2L
              P1@perm <- if (length(perm <- x@perm))
                             invertPerm(perm)
                         else seq_len(d[1L])

              P1. <- P1
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@margin <- 1L

              L <- new(paste0(if (z) "z" else "d", if (!packed) "trMatrix" else "tpMatrix"))
              L@Dim <- d
              L@uplo <- x@uplo
              L@x <- y
              if (LDL) {
                  L.jj <- diag(L, names = FALSE)
                  L@x <- y /
                      if (x@uplo == "U") {
                          if (!packed)
                              L.jj
                          else L.jj[sequence.default(seq_len(d[1L]))]
                      } else {
                          if (!packed)
                              rep(L.jj, each = d[1L])
                          else rep.int(L.jj, seq.int(to = 1L, by = -1L, length.out = d[1L]))
                      }
                  L@diag <- "U"
              }
              if (x@uplo == "U")
                  L <- ct(L. <- L)
              else L. <- ct(L)
              if (LDL) {
                  D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                  D@Dim <- d
                  D@x <- L.jj * L.jj
                  list(P1. = P1., L1 = L, D = D, L1. = L., P1 = P1)
              } else list(P1. = P1., L = L, L. = L., P1 = P1)
          })

## returning list(P1', Q, R, P2'), where A = P1' Q R P2'
setMethod("expand2", c(x = "sparseQR"),
          function(x, complete = FALSE, ...) {
              m0 <- .qr.rank.def.warn(x)
              R <- x@R
              d <- R@Dim
              dn <- x@Dimnames
              m <- d[1L]
              n <- d[2L]

              if (m0 && !is.null(dn[[1L]]))
                  length(dn[[1L]]) <- m
              Q <- .Call(sparseQR_matmult, x, NULL, 6L, complete, NULL)
              if (!complete && n < m)
                  R <- R[seq_len(n), , drop = FALSE]

              P1. <- new("pMatrix")
              P1.@Dim <- c(m, m)
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@margin <- 1L
              P1.@perm <- invertPerm(x@p, 0L, 1L)

              P2. <- new("pMatrix")
              P2.@Dim <- c(n, n)
              P2.@Dimnames <- c(list(NULL), dn[2L])
              P2.@margin <- 2L
              P2.@perm <- if (length(perm <- x@q))
                              invertPerm(perm, 0L, 1L)
                          else seq_len(n)

              if (complete)
                  list(P1. = P1., Q = Q, R = R, P2. = P2.)
              else list(P1. = P1., Q1 = Q, R1 = triu(R), P2. = P2.)
          })

## returning list(P1', L, U, P2'), where A = P1' L U P2'
setMethod("expand2", c(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames

              P1. <- new("pMatrix")
              P1.@Dim <- d
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@perm <- invertPerm(x@p, 0L, 1L)

              P2. <- new("pMatrix")
              P2.@Dim <- d
              P2.@Dimnames <- c(list(NULL), dn[2L])
              P2.@margin <- 2L
              P2.@perm <- if (length(perm <- x@q))
                              invertPerm(perm, 0L, 1L)
                          else seq_len(d[1L])

              list(P1. = P1., L = x@L, U = x@U, P2. = P2.)
          })

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
setMethod("expand2", c(x = "simplicialCholesky"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              z <- is.complex(y <- x@x)

              P1 <- new("pMatrix")
              P1@Dim <- d
              P1@Dimnames <- c(list(NULL), dn[2L])
              P1@margin <- 2L
              P1@perm <- if (length(perm <- x@perm))
                             invertPerm(perm, 0L, 1L)
                         else seq_len(d[1L])

              P1. <- P1
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@margin <- 1L

              nz <- x@nz

              L <- new(if (z) "ztCMatrix" else "dtCMatrix")
              L@Dim <- d
              L@uplo <- "L"
              L@p <- L.p <- c(0L, cumsum(nz))
              k <- sequence.default(nz, x@p[seq_along(nz)] + 1L)
              L@i <- x@i[k]
              L@x <- y <- y[k]
              if (!LDL && x@is_ll)
                  return(list(P1. = P1., L = L, L. = ct(L), P1 = P1))

              L.jj <- diag(L, names = FALSE)
              if (LDL) {
                  D <- new(if (z) "zdiMatrix" else "ddiMatrix")
                  D@Dim <- d
                  if (x@is_ll) {
                      L@x <- y / rep.int(L.jj, L.p[-1L] - L.p[-length(L.p)])
                      D@x <- L.jj * L.jj
                  } else {
                      diag(L) <- 1
                      D@x <- L.jj
                  }
                  list(P1. = P1., L1 = L, D = D, L1. = ct(L), P1 = P1)
              } else {
                  if (is.na(tmp <- min(0, L.jj)) || tmp < 0)
                      stop(gettextf("D[i,i] is NaN or negative, i=%d",
                                    which.max(is.na(L.jj) | L.jj < 0)),
                           domain = NA)
                  diag(L) <- 1
                  L@x <- L@x * rep.int(sqrt(L.jj), L.p[-1L] - L.p[-length(L.p)])
                  list(P1. = P1., L = L, L. = ct(L), P1 = P1)
              }
          })

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
setMethod("expand2", c(x = "supernodalCholesky"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              z <- is.complex(y <- x@x)

              P1 <- new("pMatrix")
              P1@Dim <- d
              P1@Dimnames <- c(list(NULL), dn[2L])
              P1@margin <- 2L
              P1@perm <- if (length(perm <- x@perm))
                             invertPerm(perm, 0L, 1L)
                         else seq_len(d[1L])

              P1. <- P1
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@margin <- 1L

              super <- x@super
              pi <- x@pi
              b <- length(super)

              nr <- pi[-1L] - pi[-b]
              nc <- super[-1L] - super[-b]
              nz <- rep.int(nr, nc)

              z <- is.complex(y <- x@x)
              L <- new(if (z) "zgCMatrix" else "dgCMatrix")
              L@Dim <- d
              L@p <- L.p <- c(0L, cumsum(nz))
              L@i <- x@s[sequence.default(nz, rep.int(pi[-b] + 1L, nc))]
              L@x <- y
              if (!LDL)
                  return(list(P1. = P1., L = L, L. = ct(L), P1 = P1))

              L.jj <- diag(L, names = FALSE)
              L@x <- y / rep.int(L.jj, L.p[-1L] - L.p[-length(L.p)])

              D <- new(if (z) "zdiMatrix" else "ddiMatrix")
              D@Dim <- d
              D@x <- L.jj * L.jj

              list(P1. = P1., L1 = L, D = D, L1. = ct(L), P1 = P1)
          })
