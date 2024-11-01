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
              d <- x@Dim
              switch(which,
                     "T" =
                         {
                             values <- x@values
                             x <- x@x
                             if (d[1L] > 0L && length(x) == 0L)
                                 new("ddiMatrix", Dim = d, x = values)
                             else if (is.double(values))
                                 new("dtrMatrix", Dim = d, x = x)
                             else new("dgeMatrix", Dim = d, x = x)
                         },
                     "Q" =, "Q." =
                         {
                             vectors <- x@vectors
                             if (d[1L] > 0L && length(vectors) == 0L)
                                 stop("missing requisite Schur vectors")
                             Q <- new("dgeMatrix", Dim = d, x = vectors)
                             switch(which, "Q" = Q, "Q." = t(Q))
                         },
                     stop(gettextf("'%1$s' is not \"%2$s\", \"%3$s\", or \"%2$s.\"",
                                   "which", "Q", "T"),
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
                             r <- new("pMatrix")
                             r@Dim <- c(m, m)
                             r@perm <- x@p + 1L
                             if (which == "P1.")
                                 r@margin <- 2L
                             r
                         },
                     "P2" =, "P2." =
                         {
                             r <- new("pMatrix")
                             r@Dim <- c(n, n)
                             r@perm <- if (length(x@q)) x@q + 1L else seq_len(n)
                             if (which == "P2")
                                 r@margin <- 2L
                             r
                         },
                     "Q"  = .Call(sparseQR_matmult, x, NULL, 6L,  TRUE, NULL),
                     "Q1" = .Call(sparseQR_matmult, x, NULL, 6L, FALSE, NULL),
                     "R"  = R,
                     "R1" = triu(if (m == n) R else R[seq_len(n), , drop = FALSE]),
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%2$s2\", \"%2$s2.\", \"%3$s\", \"%3$s1\", \"%4$s\", or \"%4$s1\"",
                                   "which", "P", "Q", "R"),
                          domain = NA))
          })

.mkL <-
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

.mkU <-
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
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- c(m, m)
                         r@perm <- asPerm(x@perm, n = m)
                         if (which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "L" = {
                         if (m <= n) {
                             r <- new("dtrMatrix")
                             r@Dim <- c(m, m)
                             r@uplo <- "L"
                             r@diag <- "U"
                         } else {
                             r <- new("dgeMatrix")
                             r@Dim <- d
                         }
                         r@x <- .mkL(x@x, m, n)
                         r
                     },
                     "U" = {
                         if (m >= n) {
                             r <- new("dtrMatrix")
                             r@Dim <- c(n, n)
                         } else {
                             r <- new("dgeMatrix")
                             r@Dim <- d
                         }
                         r@x <- .mkU(x@x, m, n)
                         r
                     },
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%3$s\", or \"%4$s\"",
                                   "which", "P", "L", "U"),
                          domain = NA))
          })

setMethod("expand1", c(x = "sparseLU"),
          function(x, which, ...) {
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- x@Dim
                         r@perm <- x@p + 1L
                         if (which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "P2" =, "P2." = {
                         r <- new("pMatrix")
                         r@Dim <- d <- x@Dim
                         r@perm <- if (length(x@q)) x@q + 1L else seq_len(d[1L])
                         if (which == "P2")
                             r@margin <- 2L
                         r
                     },
                     "L" = x@L,
                     "U" = x@U,
                     stop(gettextf("'%1$s' is not \"%2$s1\", \"%2$s1.\", \"%2$s2\", \"%2$s2.\", \"%3$s\", or \"%4$s\"",
                                   "which", "P", "L", "U"),
                          domain = NA))
          })

setMethod("expand1", c(x = "denseBunchKaufman"),
          function(x, which, ...) {
              r <- .Call(denseBunchKaufman_expand, x)
              b <- length(r) - 1L
              switch(which,
                     "DU" =, "DL" =
                         {
                             if (!endsWith(which, x@uplo))
                                 stop(gettextf("%s=\"%s\" invalid for %s@uplo=\"%s\"",
                                               "which", which, "x", x@uplo),
                                      domain = NA)
                             r[[b + 1L]]
                         },
                     "U" =, "U." =, "L" =, "L." =
                         {
                             if (!startsWith(which, x@uplo))
                                 stop(gettextf("%s=\"%s\" invalid for %s@uplo=\"%s\"",
                                               "which", which, "x", x@uplo),
                                      domain = NA)
                             if (b > 0L) {
                                 m <- r[[b]]
                                 if (b > 1L)
                                     for (i in (b - 1L):1L)
                                         m <- r[[i]] %*% m
                                 if (endsWith(which, ".")) t(m) else m
                             } else {
                                 m <- new("ddiMatrix")
                                 m@Dim <- x@Dim
                                 m@diag <- "U"
                                 m
                             }
                         },
                     stop(gettextf("'%s' is not \"%1$s\", \"D%1$s\", or \"%1$s.\"",
                                   "which", x@uplo),
                          domain = NA))
          })

setMethod("expand1", c(x = "denseCholesky"),
          function(x, which, ...) {
              d <- x@Dim
              switch(which,
                     "P1" =, "P1." =
                         {
                             r <- new("pMatrix")
                             r@Dim <- d
                             r@perm <- if (length(x@perm)) x@perm else seq_len(d[1L])
                             if (which == "P1.")
                                 r@margin <- 2L
                             r
                         },
                     "L" =, "L." =, "L1" =, "L1." =
                         {
                             packed <- length(x@x) != prod(x@Dim)
                             r <- as(x, if (!packed) "dtrMatrix" else "dtpMatrix")
                             uplo <- x@uplo
                             if (which == "L1" || which == "L1.") {
                                 r.ii <- diag(r, names = FALSE)
                                 r@x <- r@x /
                                     if (uplo == "U") {
                                         if (!packed)
                                             r.ii
                                         else r.ii[sequence.default(seq_len(d[1L]))]
                                     } else {
                                         if (!packed)
                                             rep(r.ii, each = d[1L])
                                         else rep.int(r.ii, seq.int(to = 1L, by = -1L, length.out = d[1L]))
                                     }
                                 r@diag <- "U"
                             }
                             if ((which == "L." || which == "L1.") == (uplo == "U"))
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
          })

setMethod("expand1", c(x = "simplicialCholesky"),
          function(x, which, ...) {
              switch(which,
                     "P1" =, "P1." = {
                         r <- new("pMatrix")
                         r@Dim <- d <- x@Dim
                         r@perm <- if (length(x@perm)) x@perm + 1L else seq_len(d[1L])
                         if (which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "L" =, "L." =, "L1" =, "L1." = {
                         r <- as(x, "dtCMatrix")
                         LDL. <- .CHF.is.LDL(x)
                         if (which == "L1" || which == "L1.") {
                             if (!LDL.) {
                                 r.ii <- diag(r, names = FALSE)
                                 r.p <- r@p
                                 r@x <- r@x / rep.int(r.ii , r.p[-1L] - r.p[-length(r.p)])
                             }
                             diag(r) <- 1
                         } else if (LDL.) {
                             r.ii <- diag(r, names = FALSE)
                             r.p <- r@p
                             if (anyNA(r.ii))
                                 stop(gettextf("D[i,i] is NA, i=%d",
                                               which.max(is.na(r.ii))),
                                      domain = NA)
                             if (min(r.ii) < 0)
                                 stop(gettextf("D[i,i] is negative, i=%d",
                                               which.max(r.ii < 0)),
                                      domain = NA)
                             diag(r) <- 1
                             r@x <- r@x * rep.int(sqrt(r.ii), r.p[-1L] - r.p[-length(r.p)])
                         }
                         if (which == "L" || which == "L1")
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
                         r@perm <- if (length(x@perm)) x@perm + 1L else seq_len(d[1L])
                         if (which == "P1.")
                             r@margin <- 2L
                         r
                     },
                     "L" =, "L." =, "L1" =, "L1." = {
                         r <- as(x, "dgCMatrix")
                         if (which == "L1" || which == "L1.") {
                             r.ii <- diag(r, names = FALSE)
                             r.p <- r@p
                             r@x <- r@x / rep.int(r.ii, r.p[-1L] - r.p[-length(r.p)])
                             diag(r) <- 1
                         }
                         if (which == "L" || which == "L1")
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


## METHODS FOR GENERIC: expand2
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(Q, T, Q'), where A = Q T Q'
setMethod("expand2", c(x = "denseSchur"),
          function(x, ...) {
              T <- expand1(x, "T")
              Q <- expand1(x, "Q")
              Q. <- t(Q)
              dn <- x@Dimnames
              Q @Dimnames <- c(dn[1L], list(NULL))
              Q.@Dimnames <- c(list(NULL), dn[2L])
              list(Q = Q, T = T, Q. = Q.)
          })

## returning list(P1', Q, R, P2'), where A = P1' Q R P2'
setMethod("expand2", c(x = "sparseQR"),
          function(x, complete = FALSE, ...) {
              m0 <- .qr.rank.def.warn(x)
              R <- x@R
              d <- R@Dim
              m <- d[1L]
              n <- d[2L]
              dn <- x@Dimnames
              if (m0 && !is.null(dn[[1L]]))
                  length(dn[[1L]]) <- m
              Q <- .Call(sparseQR_matmult, x, NULL, 6L, complete, NULL)
              if (!complete && n < m)
                  R <- R[seq_len(n), , drop = FALSE]
              p1 <- x@p
              p2 <- x@q
              P1. <- new("pMatrix",
                         Dim = c(m, m),
                         Dimnames = c(dn[1L], list(NULL)),
                         margin = 1L,
                         perm = invertPerm(p1, 0L, 1L))
              P2. <- new("pMatrix",
                         Dim = c(n, n),
                         Dimnames = c(list(NULL), dn[2L]),
                         margin = 2L,
                         perm = if (length(p2)) invertPerm(p2, 0L, 1L) else seq_len(n))
              if (complete)
                  list(P1. = P1., Q = Q, R = R, P2. = P2.)
              else
                  list(P1. = P1., Q1 = Q, R1 = triu(R), P2. = P2.)
          })

## returning list(P1', L, U), where A = P1' L U
setMethod("expand2", c(x = "denseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              m <- d[1L]
              n <- d[2L]

              P1. <- new("pMatrix")
              P1.@Dim <- c(m, m)
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@perm <- invertPerm(asPerm(x@perm, n = m))

              if (m <= n) {
                  L <- new("dtrMatrix")
                  L@Dim <- c(m, m)
                  L@uplo <- "L"
                  L@diag <- "U"
              } else {
                  L <- new("dgeMatrix")
                  L@Dim <- d
              }
              L@x <- .mkL(x@x, m, n)

              if (m >= n) {
                  U <- new("dtrMatrix")
                  U@Dim <- c(n, n)
              } else {
                  U <- new("dgeMatrix")
                  U@Dim <- d
              }
              U@Dimnames <- c(list(NULL), dn[2L])
              U@x <- .mkU(x@x, m, n)

              list(P1. = P1., L = L, U = U)
          })

## returning list(P1', L, U, P2'), where A = P1' L U P2'
setMethod("expand2", c(x = "sparseLU"),
          function(x, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              p1 <- x@p
              p2 <- x@q

              P1. <- new("pMatrix")
              P1.@Dim <- d
              P1.@Dimnames <- c(dn[1L], list(NULL))
              P1.@perm <- invertPerm(p1, 0L, 1L)

              P2. <- new("pMatrix")
              P2.@Dim <- d
              P2.@Dimnames <- c(list(NULL), dn[2L])
              P2.@margin <- 2L
              P2.@perm <-
                  if (length(p2))
                      invertPerm(p2, 0L, 1L)
                  else seq_len(d[1L])

              list(P1. = P1., L = x@L, U = x@U, P2. = P2.)
          })

## returning
## list(U, DU, U') where A = U DU U' and U = P[b] U[b] ... P[1] U[1]
## OR
## list(L, DL, L') where A = L DL L' and L = P[1] L[1] ... P[b] L[b]
setMethod("expand2", c(x = "denseBunchKaufman"),
          function(x, complete = FALSE, ...) {
              r <- .Call(denseBunchKaufman_expand, x)
              b <- length(r) - 1L
              if (complete) {
                  if (b > 0L)
                      r <- c(r, lapply(r[b:1L], t))
              } else {
                  if (b > 0L) {
                      m <- r[[b]]
                      if (b > 1L)
                          for (i in (b - 1L):1L)
                              m <- r[[i]] %*% m
                      r <- list(m, r[[b + 1L]], t(m))
                  } else {
                      m <- new("ddiMatrix")
                      m@Dim <- x@Dim
                      m@diag <- "U"
                      r <- list(m, r[[1L]], m)
                  }
                  names(r) <- if (x@uplo == "U") c("U", "DU", "U.") else c("L", "DL", "L.")
              }
              dn <- x@Dimnames
              if (length(r) == 1L)
                  r[[1L]]@Dimnames <- dn
              else {
                  r[[1L]]@Dimnames <- c(dn[1L], list(NULL))
                  r[[length(r)]]@Dimnames <- c(list(NULL), dn[2L])
              }
              r
          })

## returning list(P1', L1, D, L1', P1) or list(P1', L, L', P1),
## where  A = P1' L1 D L1' P1 = P1' L L' P1  and  L = L1 sqrt(D)
setMethod("expand2", c(x = "denseCholesky"),
          function(x, LDL = TRUE, ...) {
              packed <- length(x@x) != prod(x@Dim)
              d <- x@Dim
              dn <- x@Dimnames
              uplo <- x@uplo
              perm <- x@perm

              P <- new("pMatrix")
              P@Dim <- d
              P@Dimnames <- c(list(NULL), dn[2L])
              P@margin <- 2L
              P@perm <- if (length(perm)) invertPerm(perm) else seq_len(d[1L])

              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L

              X <- as(x, if (!packed) "dtrMatrix" else "dtpMatrix")
              if (LDL) {
                  L.ii <- diag(X, names = FALSE)
                  X@x <- X@x /
                      if (uplo == "U") {
                          if (!packed)
                              L.ii
                          else L.ii[sequence.default(seq_len(d[1L]))]
                      } else {
                          if (!packed)
                              rep(L.ii, each = d[1L])
                          else rep.int(L.ii, seq.int(to = 1L, by = -1L, length.out = d[1L]))
                      }
                  X@diag <- "U"
              }
              L  <- if (uplo == "U") t(X) else   X
              L. <- if (uplo == "U")   X  else t(X)
              if (LDL) {
                  D <- new("ddiMatrix")
                  D@Dim <- d
                  D@x <- L.ii * L.ii
                  list(P1. = P., L1 = L, D = D, L1. = L., P1 = P)
              } else list(P1. = P., L = L, L. = L., P1 = P)
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
              P@perm <- if (length(perm))
                            invertPerm(perm, 0L, 1L)
                        else seq_len(d[1L])

              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L

              L <- as(x, "dtCMatrix")
              LDL. <- .CHF.is.LDL(x)
              if (!LDL && !LDL.)
                  return(list(P1. = P., L = L, L. = t(L), P1 = P))
              L.ii <- diag(L, names = FALSE)
              L.p <- L@p
              if (!LDL) {
                  if (anyNA(L.ii))
                      stop(gettextf("D[i,i] is NA, i=%d",
                                    which.max(is.na(L.ii))),
                           domain = NA)
                  if (min(L.ii) < 0)
                      stop(gettextf("D[i,i] is negative, i=%d",
                                    which.max(L.ii < 0)),
                           domain = NA)
                  diag(L) <- 1
                  L@x <- L@x * rep.int(sqrt(L.ii), L.p[-1L] - L.p[-length(L.p)])
                  return(list(P1. = P., L = L, L. = t(L), P1 = P))
              }
              D <- new("ddiMatrix")
              D@Dim <- d
              if (LDL.) {
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
              P@perm <- if (length(perm))
                            invertPerm(perm, 0L, 1L)
                        else seq_len(d[1L])

              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L

              L <- as(x, "dgCMatrix")
              if (!LDL)
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
