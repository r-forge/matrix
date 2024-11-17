## METHODS FOR GENERIC: solve
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.solve.checkDim1 <-
function(nrow.a, ncol.a) {
    if (nrow.a != ncol.a)
        stop(gettextf("'%s' is not square",
                      "a"),
             domain = NA)
}

.solve.checkDim2 <-
function(nrow.a, nrow.b) {
    if (nrow.a != nrow.b)
        stop(gettextf("dimensions of '%s' and '%s' are inconsistent",
                      "a", "b"),
             domain = NA)
}

.solve.checkCondition <-
function(a, tol, rcond.a = rcond(a)) {
    if (tol > 0 && a@Dim[1L] > 0L && rcond.a < tol)
        stop(gettextf("'%s' is computationally singular, rcond(%s)=%g",
                      "a", "a", rcond.a),
             domain = NA)
}

## Factorize  A  as  P1 A P2 = L U .  Then :
##
##     kappa(A) <= kappa(P1') * kappa(L) * kappa(U) * kappa(P2')
##              ==          1 * kappa(L) * kappa(U) * 1
##
## If  U  is diagonally dominant, i.e., if there exists  a > 1
## such that
##
##     abs(U[i, i]) >= a * sum(abs(U[i, -i]))        i = 1, ... , n
##
## then :
##
##     kappa(U) <= (a + 1)/(a - 1) * max(abs(diag(U)))/min(abs(diag(U)))
##
## The bound contracts as  a --> Inf  and in the limit we have:
##
##     kappa(A)/kappa(L) <= kappa(U) <= max(abs(diag(U)))/min(abs(diag(U)))

.solve.checkConditionBound <-
function(u, tol, rad.u = range(abs(diag(u, names = FALSE)))) {
    if (tol > 0 && u@Dim[1L] > 0L) {
        r <- rad.u[1L] / rad.u[2L]
        if (r < tol)
            stop(gettextf("'%s' is computationally singular, min(d)/max(d)=%g, d=abs(diag(U))",
                          "a", r),
                 domain = NA)
    }
}

.solve.checkDiagonal <-
function(diag.a) {
    if (any(zero <- !diag.a, na.rm = TRUE))
        stop(gettextf("'%s' is exactly singular, D[i,i]=0, i=%d",
                      "a", which.max(zero)),
             domain = NA)
}

.solve.checkIndex <-
function(perm, margin, n) {
    if (anyDuplicated.default(perm)) {
        i <- seq_len(n)[-perm][1L]
        if (margin == 1L)
        stop(gettextf("'%s' is exactly singular, J[,j]=0, j=%d",
                      "a", i),
             domain = NA)
        else
        stop(gettextf("'%s' is exactly singular, J[i,]=0, i=%d",
                      "a", i),
             domain = NA)
    }
}

for (.cl.a in c("MatrixFactorization", "Matrix")) {
setMethod("solve", c(a = .cl.a, b = "vector"),
          function(a, b, ...)
         drop(solve(a, .m2dense(b, ",ge"), ...)))

setMethod("solve", c(a = .cl.a, b = "matrix"),
          function(a, b, ...)
              solve(a, .m2dense(b, ",ge"), ...))

setMethod("solve", c(a = .cl.a, b = "denseMatrix"),
          function(a, b, ...)
              solve(a, .dense2sparse(b, "C"), ...))

setMethod("solve", c(a = .cl.a, b = "CsparseMatrix"),
          function(a, b, ...)
              solve(a, .sparse2dense(b, FALSE), ...))

setMethod("solve", c(a = .cl.a, b = "RsparseMatrix"),
          function(a, b, ...)
              solve(a, .M2C(b), ...))

setMethod("solve", c(a = .cl.a, b = "TsparseMatrix"),
          function(a, b, ...)
              solve(a, .M2C(b), ...))

setMethod("solve", c(a = .cl.a, b = "diagonalMatrix"),
          function(a, b, ...)
              solve(a, .diag2sparse(b, ",", "g", "C"), ...))

setMethod("solve", c(a = .cl.a, b = "indMatrix"),
          function(a, b, ...)
              solve(a, .ind2sparse(b, ",", "C"), ...))

setMethod("solve", c(a = .cl.a, b = "sparseVector"),
          function(a, b, ...)
         drop(solve(a, .V2C(b), ...)))
}

for (.cl.b in c("Matrix", "sparseVector"))
setMethod("solve", c(a = "matrix", b = .cl.b),
          function(a, b, ...)
              solve(.m2dense(a, ",ge"), b, ...))

setMethod("solve", c(a = "denseSchur", b = "missing"),
          function(a, b, ...) {
              adim <- a@Dim
              if (adim[1L] > 0L && length(a@vectors) == 0L)
                  stop("missing requisite Schur vectors")
              Q <- expand1(a, "Q")
              T <- expand1(a, "T")
              r <- Q %*% solve(T, ct(Q))
              r@Dimnames <- a@Dimnames[2:1]
              r
          })

setMethod("solve", c(a = "denseSchur", b = "denseMatrix"),
          function(a, b, ...) {
              adim <- a@Dim
              if (adim[1L] > 0L && length(a@vectors) == 0L)
                  stop("missing requisite Schur vectors")
              .solve.checkDim2(adim[1L], b@Dim[1L])
              Q <- expand1(a, "Q")
              T <- expand1(a, "T")
              r <- Q %*% solve(T, crossprod(Q, b, trans = "C"))
              r@Dimnames <- c(a@Dimnames[2L], dimnames(b)[2L])
              r
          })

setMethod("solve", c(a = "denseLU", b = "missing"),
          function(a, b, ...)
              .Call(denseLU_solve, a, NULL))

setMethod("solve", c(a = "denseLU", b = "denseMatrix"),
          function(a, b, ...)
              .Call(denseLU_solve, a, .M2gen(b, ",")))

setMethod("solve", c(a = "denseBunchKaufman", b = "missing"),
          function(a, b, ...)
              .Call(denseBunchKaufman_solve, a, NULL))

setMethod("solve", c(a = "denseBunchKaufman", b = "denseMatrix"),
          function(a, b, ...)
              .Call(denseBunchKaufman_solve, a, .M2gen(b, ",")))

setMethod("solve", c(a = "denseCholesky", b = "missing"),
          function(a, b, ...)
              .Call(denseCholesky_solve, a, NULL))

setMethod("solve", c(a = "denseCholesky", b = "denseMatrix"),
          function(a, b, ...)
              .Call(denseCholesky_solve, a, .M2gen(b, ",")))

setMethod("solve", c(a = "sparseQR", b = "missing"),
          function(a, b, sparse = TRUE, ...) {
              z <- is.complex(a@R@x)
              b <- new(if (z) "zgeMatrix" else "dgeMatrix")
              b@Dim <- { n <- a@Dim[2L]; c(n, n) }
              b@x <- as.vector(diag(if (z) 1+0i else 1, n))
              r <- qr.coef(a, b)
              if (is.na(sparse) || sparse) .dense2sparse(r, "C") else r
          })

setMethod("solve", c(a = "sparseQR", b = "denseMatrix"),
          function(a, b, sparse = FALSE, ...) {
              r <- qr.coef(a, .M2gen(b, ","))
              if (is.na(sparse) || sparse) .dense2sparse(r, "C") else r
          })

setMethod("solve", c(a = "sparseQR", b = "CsparseMatrix"),
          function(a, b, sparse = TRUE, ...) {
              r <- qr.coef(a, .sparse2dense(.M2gen(b, ",")))
              if (is.na(sparse) || sparse) .dense2sparse(r, "C") else r
          })

setMethod("solve", c(a = "sparseLU", b = "missing"),
          function(a, b, sparse = TRUE,
                   tol = .Machine$double.eps, ...) {
              .solve.checkConditionBound(a@U, tol)
              .Call(sparseLU_solve, a, NULL, sparse)
          })

setMethod("solve", c(a = "sparseLU", b = "denseMatrix"),
          function(a, b, sparse = FALSE,
                   tol = .Machine$double.eps, ...) {
              .solve.checkConditionBound(a@U, tol)
              if (is.na(sparse) || sparse)
                  b <- .dense2sparse(b, "C")
              .Call(sparseLU_solve, a, .M2gen(b, ","), sparse)
          })

setMethod("solve", c(a = "sparseLU", b = "CsparseMatrix"),
          function(a, b, sparse = TRUE,
                   tol = .Machine$double.eps, ...) {
              .solve.checkConditionBound(a@U, tol)
              if (!(is.na(sparse) || sparse))
                  b <- .sparse2dense(b, FALSE)
              .Call(sparseLU_solve, a, .M2gen(b, ","), sparse)
          })

setMethod("solve", c(a = "sparseCholesky", b = "missing"),
          function(a, b, sparse = TRUE,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt",
                              "D", "P", "Pt"),
                   ...) {
              if (!missing(system)) {
                  ## Do the "right" thing :
                  if (identical(system, "D")) {
                      z <- is.complex(a@x)
                      r <- new(if (z) "zdiMatrix" else "ddiMatrix")
                      r@Dim <- a@Dim
                      r@Dimnames <- a@Dimnames[2:1]
                      if (.CHF.is.LDL(a))
                          r@x <- a@x[a@p + 1L]
                      else r@diag <- "U"
                      return(r)
                  } else if (identical(system, "P") ||
                             identical(system, "Pt")) {
                      perm <- a@perm
                      r <- new("pMatrix")
                      r@Dim <- d <- a@Dim
                      r@Dimnames <- a@Dimnames[2:1]
                      if (system == "Pt")
                          r@margin <- 2L
                      r@perm <- if (length(perm))
                                    perm + 1L
                                else seq_len(d[1L])
                      return(r)
                  }
              }
              .Call(sparseCholesky_solve, a, NULL, sparse, system)
          })

setMethod("solve", c(a = "sparseCholesky", b = "denseMatrix"),
          function(a, b, sparse = FALSE,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt",
                              "D", "P", "Pt"),
                   ...) {
              if (is.na(sparse) || sparse)
                  b <- .dense2sparse(b, "C")
              .Call(sparseCholesky_solve, a, .M2gen(b, ","), sparse, system)
          })

setMethod("solve", c(a = "sparseCholesky", b = "CsparseMatrix"),
          function(a, b, sparse = TRUE,
                   system = c("A", "LDLt", "LD", "DLt", "L", "Lt",
                              "D", "P", "Pt"),
                   ...) {
              if (!(is.na(sparse) || sparse))
                  b <- .sparse2dense(b, FALSE)
              .Call(sparseCholesky_solve, a, .M2gen(b, ","), sparse, system)
          })

setMethod("solve", c(a = "denseMatrix", b = "missing"),
          function(a, b, tol = .Machine$double.eps, ...) {
              adim <- a@Dim
              .solve.checkDim1(adim[1L], adim[2L])
              a <- .M2kind(a, ",")
              .solve.checkCondition(a, tol)
              trf <-
              switch(.M.shape(a, 2L),
                     "g" = lu(a, warnSing = FALSE),
                     "s" = BunchKaufman(a, warnSing = FALSE),
                     "p" = Cholesky(a, perm = FALSE),
                     "t" = return(.Call(trMatrix_solve, a, NULL)),
                     stop("should never happen ..."))
              solve(trf, ...)
          })

setMethod("solve", c(a = "denseMatrix", b = "denseMatrix"),
          function(a, b, tol = .Machine$double.eps, ...) {
              adim <- a@Dim
              .solve.checkDim1(adim[1L], adim[2L])
              a <- .M2kind(a, ",")
              .solve.checkCondition(a, tol)
              trf <-
              switch(.M.shape(a, 2L),
                     "g" = lu(a, warnSing = FALSE),
                     "s" = BunchKaufman(a, warnSing = FALSE),
                     "p" = Cholesky(a, perm = FALSE),
                     "t" =
                         {
                             r <- .Call(trMatrix_solve, .M2kind(a, ","), .M2gen(b, ","))
                             if (.M.shape(b) == "t" && a@uplo == b@uplo) {
                                 r <- if (a@uplo == "U") triu(r) else tril(r)
                                 if (a@diag != "N" && b@diag != "N")
                                     r <- ..diagN2U(r)
                             }
                             return(r)
                         },
                     stop("should never happen ..."))
              solve(trf, b, ...)
          })

setMethod("solve", c(a = "CsparseMatrix", b = "missing"),
          function(a, b, sparse = TRUE, ...) {
              adim <- a@Dim
              .solve.checkDim1(adim[1L], adim[2L])
              a <- .M2kind(a, ",")
              trf <-
              switch(.M.shape(a, 2L),
                     "g" = lu(a, errSing = TRUE),
                     "s" = tryCatch(Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                                    error = function(e) lu(a, errSing = TRUE)),
                     "p" = Cholesky(a, perm = TRUE, LDL = FALSE, super = FALSE),
                     "t" =
                         {
                             if (a@diag != "N")
                                 a <- ..diagU2N(a)
                             return(.Call(tCMatrix_solve, a, NULL, sparse))
                         },
                     stop("should never happen ..."))
              solve(trf, sparse = sparse, ...)
          })

setMethod("solve", c(a = "CsparseMatrix", b = "denseMatrix"),
          function(a, b, sparse = FALSE, ...) {
              adim <- a@Dim
              .solve.checkDim1(adim[1L], adim[2L])
              a <- .M2kind(a, ",")
              trf <-
              switch(.M.shape(a, 2L),
                     "g" = lu(a, errSing = TRUE),
                     "s" = tryCatch(Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                                    error = function(e) lu(a, errSing = TRUE)),
                     "p" = Cholesky(a, perm = TRUE, LDL = FALSE, super = FALSE),
                     "t" =
                         {
                             if (a@diag != "N")
                                 a <- ..diagU2N(a)
                             if (is.na(sparse) || sparse)
                                 b <- .dense2sparse(b, "C")
                             r <- .Call(tCMatrix_solve, .M2kind(a, ","), .M2gen(b, ","), sparse)
                             if (.M.shape(b) == "t" && a@uplo == b@uplo) {
                                 r <- if (a@uplo == "U") triu(r) else tril(r)
                                 if (a@diag != "N" && b@diag != "N")
                                     r <- ..diagN2U(r)
                             }
                             return(r)
                         },
                     stop("should never happen ..."))
              solve(trf, b, sparse = sparse, ...)
          })

setMethod("solve", c(a = "CsparseMatrix", b = "CsparseMatrix"),
          function(a, b, sparse = TRUE, ...) {
              adim <- a@Dim
              .solve.checkDim1(adim[1L], adim[2L])
              a <- .M2kind(a, ",")
              trf <-
              switch(.M.shape(a, 2L),
                     "g" = lu(a, errSing = TRUE),
                     "s" = tryCatch(Cholesky(a, perm = TRUE, LDL = TRUE, super = FALSE),
                                    error = function(e) lu(a, errSing = TRUE)),
                     "p" = Cholesky(a, perm = TRUE, LDL = FALSE, super = FALSE),
                     "t" =
                         {
                             if (a@diag != "N")
                                 a <- ..diagU2N(a)
                             if (!(is.na(sparse) || sparse))
                                 b <- .sparse2dense(b, FALSE)
                             r <- .Call(tCMatrix_solve, .M2kind(a, ","), .M2gen(b, ","), sparse)
                             if (.M.shape(b) == "t" && a@uplo == b@uplo) {
                                 r <- if (a@uplo == "U") triu(r) else tril(r)
                                 if (a@diag != "N" && b@diag != "N")
                                     r <- ..diagN2U(r)
                             }
                             return(r)
                         },
                     stop("should never happen ..."))
              solve(trf, b, sparse = sparse, ...)
          })

## TODO: implement triangular solver for [dz]tRMatrix, so that we can
##       handle  A = <[dz]gRMatrix>  and  A' = .tCRT(A)  like so:
##
##                   P1 A' P2 = L U
##       A x = b  <==================>  x = P1' inv(L') inv(U') P2' b

setMethod("solve", c(a = "RsparseMatrix", b = "missing"),
          function(a, b, ...)
              solve(.M2C(a), ...))

setMethod("solve", c(a = "RsparseMatrix", b = "denseMatrix"),
          function(a, b, ...)
              solve(.M2C(a), b, ...))

setMethod("solve", c(a = "RsparseMatrix", b = "CsparseMatrix"),
          function(a, b, ...)
              solve(.M2C(a), b, ...))

setMethod("solve", c(a = "TsparseMatrix", b = "missing"),
          function(a, b, ...)
              solve(.M2C(a), ...))

setMethod("solve", c(a = "TsparseMatrix", b = "denseMatrix"),
          function(a, b, ...)
              solve(.M2C(a), b, ...))

setMethod("solve", c(a = "TsparseMatrix", b = "CsparseMatrix"),
          function(a, b, ...)
              solve(.M2C(a), b, ...))

setMethod("solve", c(a = "diagonalMatrix", b = "missing"),
          function(a, b, ...) {
              ax <- a@x
              z <- is.complex(ax)
              r <- new(if (z) "zdiMatrix" else "ddiMatrix")
              r@Dim <- a@Dim
              r@Dimnames <- a@Dimnames[2:1]
              if (a@diag == "N") {
                  if (.M.kind(a) == "n")
                      ax <- ax | is.na(ax)
                  .solve.checkDiagonal(ax)
                  r@x <- 1 / ax
              }
              r
          })

setMethod("solve", c(a = "diagonalMatrix", b = "vector"),
          function(a, b, ...) {
              .solve.checkDim2(a@Dim[1L], length(b))
              ax <- a@x
              z <- is.complex(ax) || is.complex(b)
              as. <- if (z) as.complex else as.double
              r <-
              if (a@diag == "N") {
                  if (.M.kind(a) == "n")
                      ax <- ax | is.na(ax)
                  .solve.checkDiagonal(ax)
                  as.(b) / ax
              } else as.(b)
              names(r) <- a@Dimnames[[2L]]
              r
          })

setMethod("solve", c(a = "diagonalMatrix", b = "matrix"),
          function(a, b, ...) {
              .solve.checkDim2(a@Dim[1L], dim(b)[1L])
              solve(a) %*% b
          })

setMethod("solve", c(a = "indMatrix", b = "missing"),
          function(a, b, ...) {
              adim <- a@Dim
              aperm <- a@perm
              amargin <- a@margin
              if (.M.class(a, 0L) != "pMatrix") {
                  .solve.checkDim1(adim[1L], adim[2L])
                  .solve.checkIndex(aperm, amargin, adim[1L])
              }
              r <- new("pMatrix")
              r@Dim <- adim
              r@Dimnames <- a@Dimnames[2:1]
              r@perm <- aperm
              r@margin <- if (amargin == 1L) 2L else 1L
              r
          })

setMethod("solve", c(a = "indMatrix", b = "vector"),
          function(a, b, ...) {
              adim <- a@Dim
              aperm <- a@perm
              amargin <- a@margin
              if (.M.class(a, 0L) != "pMatrix") {
                  .solve.checkDim1(adim[1L], adim[2L])
                  .solve.checkIndex(aperm, amargin, adim[1L])
              }
              .solve.checkDim2(adim[1L], length(b))
              q <- if (amargin == 1L) invertPerm(aperm) else aperm
              z <- is.complex(b)
              as. <- if (z) as.complex else as.double
              r <- as.(b)[q]
              names(r) <- a@Dimnames[[2L]]
              r
          })

setMethod("solve", c(a = "indMatrix", b = "matrix"),
          function(a, b, ...) {
              .solve.checkDim2(a@Dim[1L], dim(b)[1L])
              solve(a) %*% b
          })

for (.cl.a in c("diagonalMatrix", "indMatrix"))
for (.cl.b in c("denseMatrix", "CsparseMatrix", "RsparseMatrix",
                "TsparseMatrix", "diagonalMatrix", "indMatrix"))
setMethod("solve", c(a = .cl.a, b = .cl.b),
          function(a, b, ...) {
              .solve.checkDim2(a@Dim[1L], b@Dim[1L])
              solve(a) %*% b
          })

rm(.cl.a, .cl.b)


## These are exported and used only (?) by MatrixModels :

## a=dgCMatrix
## b=vector or 1-column matrix
## x=double vector
.solve.dgC.qr <- function(a, b, order = 3L, check = TRUE) {
    if (check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .Call(dgCMatrix_qrsol, a, b, order)
}

## a=dgCMatrix
## b=vector, matrix, or Matrix
## x=dgCMatrix or dgeMatrix
.solve.dgC.lu <- function(a, b, tol = .Machine$double.eps, check = TRUE) {
    if (check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    trf <- lu(a, errSing = TRUE)
    solve(trf, b, tol = tol)
}

## a=dgCMatrix
## b=vector or 1-column matrix
## x=list(L, coef, Xty, resid)
.solve.dgC.chol <- function(a, b, check = TRUE) {
    if (check && !is(a, "dgCMatrix"))
        a <- as(as(as(a, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    .Call(dgCMatrix_cholsol, a, b)
}
