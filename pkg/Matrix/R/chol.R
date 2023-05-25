## METHODS FOR GENERIC: chol
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol", signature(x = "generalMatrix"),
	  function(x, uplo = "U", ...) {
              ch <- chol(forceSymmetric(x, uplo), ...)
              ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("chol", signature(x = "triangularMatrix"),
	  function(x, uplo = "U", ...) {
              if(identical(uplo, x@uplo)) {
                  ch <- chol(forceSymmetric(x, uplo), ...)
                  ch@Dimnames <- x@Dimnames # restore asymmetric 'Dimnames'
                  ch
              } else chol(forceDiagonal(x, x@diag), ...)
          })

setMethod("chol", signature(x = "symmetricMatrix"),
	  function(x, ...)
              chol(as(x, "dMatrix"), ...))

setMethod("chol", signature(x = "diagonalMatrix"),
	  function(x, ...)
              chol(..diag2d(x), ...))

setMethod("chol", signature(x = "dsyMatrix"),
          function(x, ...) {
              ch <- as(Cholesky(x), "dtrMatrix")
              if(ch@uplo != "U") t(ch) else ch
          })

setMethod("chol", signature(x = "dspMatrix"),
          function(x, ...) {
              ch <- as(Cholesky(x), "dtpMatrix")
              if(ch@uplo != "U") t(ch) else ch
          })

for(.cl in paste0("ds", c("C", "R", "T"), "Matrix"))
setMethod("chol", signature(x = .cl),
          function(x, pivot = FALSE, ...) {
              ch <- Cholesky(x, perm = pivot, LDL = FALSE, super = FALSE)
              t(as(ch, "dtCMatrix")) # FIXME? give dtRMatrix, dtTMatrix?
          })

setMethod("chol", signature(x = "ddiMatrix"),
          function(x, ...) {
              if(length(y <- x@x)) {
                  if(is.na(min.y <- min(y)) || min.y <= 0)
                      stop("chol(x) is undefined: 'x' is not positive definite")
                  x@x <- sqrt(y)
              }
              x
          })


## METHODS FOR GENERIC: Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("Cholesky", signature(A = "generalMatrix"),
	  function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              ch@Dimnames <- A@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("Cholesky", signature(A = "triangularMatrix"),
	  function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              ch@Dimnames <- A@Dimnames # restore asymmetric 'Dimnames'
              ch
          })

setMethod("Cholesky", signature(A = "symmetricMatrix"),
	  function(A, ...)
              Cholesky(as(A, "dMatrix"), ...))

setMethod("Cholesky", signature(A = "diagonalMatrix"),
	  function(A, ...)
              Cholesky(..diag2d(A), ...))

setMethod("Cholesky", signature(A = "dsyMatrix"),
          function(A, ...)
              .Call(dpoMatrix_trf, A, 2L))

setMethod("Cholesky", signature(A = "dspMatrix"),
          function(A, ...)
              .Call(dppMatrix_trf, A, 2L))

setMethod("Cholesky", signature(A = "dsCMatrix"),
          function(A, perm = TRUE, LDL = !super, super = FALSE,
                   Imult = 0, ...)
              .Call(dpCMatrix_trf, A, perm, LDL, super, Imult))

setMethod("Cholesky", signature(A = "dsRMatrix"),
          function(A, ...)
              Cholesky(.tCR2RC(A), ...))

setMethod("Cholesky", signature(A = "dsTMatrix"),
          function(A, ...)
              Cholesky(.T2C(A), ...))

setMethod("Cholesky", signature(A = "ddiMatrix"),
          function(A, ...) {
              if(length(y <- A@x) && (is.na(min.y <- min(y)) || min.y <= 0))
                  stop("Cholesky(A) is undefined: 'A' is not positive definite")
              n <- (d <- A@Dim)[1L]
              r <- new("dCHMsimpl")
              r@Dim <- d
              r@Dimnames <- A@Dimnames
              r@colcount <- r@nz <- rep.int(1L, n)
              r@type <- c(0L, 0L, 0L, 1L, 0L, 0L)
              r@p <- 0:n
              r@i <- s <- seq.int(0L, length.out = n)
              r@x <- if(length(y)) y else rep.int(1, n)
              r@nxt <- c(seq_len(n), -1L, 0L)
              r@prv <- c(n + 1L, s, -1L) # @<- will error if n + 1L overflows
              r
          })

setMethod("Cholesky", signature(A = "matrix"),
	  function(A, uplo = "U", ...) {
              ch <- Cholesky(forceSymmetric(A, uplo), ...)
              if(!is.null(dn <- dimnames(A)))
                  ch@Dimnames <- dn # restore asymmetric 'Dimnames'
              ch
          })


## METHODS FOR GENERIC: chol2inv
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("chol2inv", signature(x = "denseMatrix"), # ->dtrMatrix
	  function(x, ...)
              chol2inv(unpack(.M2tri(..dense2d(x))), ...))

setMethod("chol2inv", signature(x = "dtrMatrix"),
	  function (x, ...) {
	      if(x@diag != "N")
                  x <- ..diagU2N(x)
	      .Call(dtrMatrix_chol2inv, x)
	  })

setMethod("chol2inv", signature(x = "sparseMatrix"),
	  function (x, ...)
	      tcrossprod(solve(.M2tri(x))))

setMethod("chol2inv", signature(x = "diagonalMatrix"),
	  function (x, ...)
	      tcrossprod(solve(x)))

setMethod("chol2inv", signature(x = "CHMfactor"),
	  function (x, ...)
	      solve(x, system = "A"))


## METHODS FOR CLASS: p?Cholesky
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## returning list(L, L') or list(L1, D, L1'),
## where  A = L L' = L1 D L1'  and  L = L1 sqrt(D)
.def.unpacked <- .def.packed <- function(x, LDL = TRUE, ...) {
    x <- as(x, .CL)
    dn <- x@Dimnames
    up <- x@uplo == "U"
    nu <- x@diag == "N"
    if(LDL && nu) {
        L.ii <- diag(x, names = FALSE)
        x@x <- x@x / if(up) .UP else .LO
        x@diag <- "U"
    }
    L  <- if(up) t(x) else   x
    L. <- if(up)   x  else t(x)
    L @Dimnames <- c(dn[1L], list(NULL))
    L.@Dimnames <- c(list(NULL), dn[2L])
    if(LDL) {
        D <- new("ddiMatrix")
        D@Dim <- x@Dim
        if(nu)
            D@x <- L.ii * L.ii
        else D@diag <- "U"
        list(L1 = L, D = D, L1. = L.)
    } else list(L = L, L. = L.)
}

body(.def.unpacked) <-
    do.call(substitute,
            list(body(.def.unpacked),
                 list(.CL = "dtrMatrix",
                      .UP = quote(L.ii),
                      .LO = quote(rep(L.ii, each = x@Dim[1L])))))
body(.def.packed) <-
    do.call(substitute,
            list(body(.def.packed),
                 list(.CL = "dtpMatrix",
                      .UP = quote(L.ii[sequence.default(seq_len(x@Dim[1L]))]),
                      .LO = quote(rep.int(L.ii, seq.int(to = 1L, by = -1L, length.out = x@Dim[1L]))))))

## returning list(L, L') or list(L1, D, L1'), where A = L L' = L1 D L1'
setMethod("expand2", signature(x =  "Cholesky"), .def.unpacked)
setMethod("expand2", signature(x = "pCholesky"), .def.packed)

rm(.def.unpacked, .def.packed)


## METHODS FOR CLASS: CHMfactor
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diag", signature(x = "CHMfactor"),
          function(x, nrow, ncol, names)
              .Call(CHMfactor_diag_get, x, TRUE))

.CHM.is.perm <- function(x)
    !as.logical(x@type[1L])
.CHM.is.LDL <- function(x)
    !as.logical(x@type[2L])
.CHM.is.super <- function(x)
     as.logical(x@type[3L])

# Exported:
isLDL <- function(x) {
    if(is(x, "CHMfactor"))
        .CHM.is.LDL(x)
    else stop("'x' does not inherit from virtual class CHMfactor")
}

setAs("CHMsimpl", "dtCMatrix",
      function(from) {
          nz <- from@nz
          k <- sequence.default(nz, from@p[seq_along(nz)] + 1L)

          to <- new("dtCMatrix")
          to@Dim  <- from@Dim
          to@Dimnames <- from@Dimnames
          to@uplo <- "L"
          to@p <- c(0L, cumsum(nz))
          to@i <- from@i[k]
          to@x <- from@x[k]
          to
      })

setAs("CHMsimpl", "CsparseMatrix",
      function(from) {
          to <- as(from, "dtCMatrix")
          if(.CHM.is.LDL(from)) {
              L.ii <- sqrt(diag(to, names = FALSE))
              L.p <- to@p
              diag(to) <- 1
              to@x <- to@x * rep.int(L.ii, L.p[-1L] - L.p[-length(L.p)])
          }
          to
      })

setAs("CHMsuper", "dgCMatrix",
      function(from) {
          super <- from@super
          pi <- from@pi
          b <- length(super)

          nr <- pi[-1L] - pi[-b]
          nc <- super[-1L] - super[-b]
          dp <- rep.int(nr, nc)

          to <- new("dgCMatrix")
          to@Dim <- from@Dim
          to@Dimnames <- from@Dimnames
          to@p <- c(0L, cumsum(dp))
          to@i <- from@s[sequence.default(dp, rep.int(pi[-b] + 1L, nc))]
          to@x <- from@x
          to
      })

setAs("CHMsuper", "CsparseMatrix",
      function(from)
          as(from, "dgCMatrix"))

setAs("CHMfactor", "pMatrix",
      function(from) {
          r <- new("pMatrix")
          r@Dim <- d <- from@Dim
          r@perm <- if(length(perm <- from@perm)) perm + 1L else seq_len(d[1L])
          r
      })

## returning list(P1', L, L', P1) or list(P1', L1, D, L1', P1),
## where  A = P1' L L' P1 = P1' L1 D L1' P1  and  L = L1 sqrt(D)
setMethod("expand2", signature(x = "CHMsimpl"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              perm <- x@perm
              P <- new("pMatrix",
                       Dim = d,
                       Dimnames = c(list(NULL), dn[2L]),
                       margin = 2L,
                       perm = if(length(perm))
                                  invertPerm(perm, 0L, 1L)
                              else seq_len(d[1L]))
              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L
              L <- as(x, "dtCMatrix")
              LDL. <- .CHM.is.LDL(x)
              if(!LDL && !LDL.)
                  return(list(P1. = P., L = L, L. = t(L), P1 = P))
              L.ii <- diag(L, names = FALSE)
              L.p <- L@p
              if(!LDL) {
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

## returning list(P1', L, L', P1) or list(P1', L1, D, L1', P1),
## where  A = P1' L L' P1 = P1' L1 D L1' P1  and  L = L1 sqrt(D)
setMethod("expand2", signature(x = "CHMsuper"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              perm <- x@perm
              P <- new("pMatrix",
                       Dim = d,
                       Dimnames = c(list(NULL), dn[2L]),
                       margin = 2L,
                       perm = if(length(perm))
                                  invertPerm(perm, 0L, 1L)
                              else seq_len(d[1L]))
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
setMethod("expand", signature(x = "CHMfactor"),
          function(x, ...)
              list(P = as(x, "pMatrix"), L = as(x, "CsparseMatrix")))

.updateCHMfactor <- function(object, parent, mult = 0)
    .Call(CHMfactor_update, object, parent, mult)

setMethod("update", signature(object = "CHMfactor"),
	  function(object, parent, mult = 0, ...) {
              s <- .M.repr(parent)
              if(!nzchar(s))
                  stop("'parent' is not formally sparse")
              if(s != "C")
                  parent <- as(parent, "CsparseMatrix")
              s <- .M.shape(parent)
              if(s != "s") {
                  Matrix.msg("'parent' is not formally symmetric; factorizing tcrossprod(parent)")
                  if(s == "t" && parent@diag != "N")
                      parent <- ..diagU2N(parent)
              }
              s <- .M.kind(parent)
              if(s != "d")
                  parent <- .sparse2kind(parent, "d")
              .updateCHMfactor(object, parent, mult)
          })

.updownCHMfactor <- function(update, C, L)
    .Call(CHMfactor_updown, L, C, update)

for(.cl in c("Matrix", "matrix")) {
setMethod("updown",
          signature(update = "character", C = .cl, L = "CHMfactor"),
          function(update, C, L)
              updown(!identical(update, "-"), C, L))

setMethod("updown",
          signature(update = "logical", C = .cl, L = "CHMfactor"),
          function(update, C, L)
              updown(update, as(as(C, "CsparseMatrix"), "dMatrix"), L))
}
rm(.cl)

for(.cl in c("dgCMatrix", "dsCMatrix"))
setMethod("updown",
          signature(update = "logical", C = .cl, L = "CHMfactor"),
          function(update, C, L) {
              if(length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .updownCHMfactor(update, C, L)
          })
rm(.cl)

setMethod("updown",
          signature(update = "logical", C = "dtCMatrix", L = "CHMfactor"),
          function(update, C, L) {
              if(C@diag != "N")
                  C <- ..diagU2N(C)
              if(length(perm <- L@perm))
                  C <- C[perm + 1L, , drop = FALSE]
              .updownCHMfactor(update, C, L)
          })
