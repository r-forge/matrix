## METHODS FOR GENERIC: determinant
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.mkDet <-
function(modulus = sum(log(abs(x))),
         logarithm = TRUE,
         sign =
             if (is.complex(x))
                 NA_integer_
             else if (prod(x) < 0)
                 -1L
             else 1L,
         x) {
    if (!logarithm)
        modulus <- exp(modulus)
    attr(modulus, "logarithm") <- logarithm
    val <- list(modulus = modulus, sign = sign)
    class(val) <- "det"
    val
}

## base::det calls base::determinant, which is not S4 generic,
## so we export own our 'det' calling Matrix::determinant
det <- base::det
environment(det) <- environment() # the Matrix namespace

for (.cl in c("Matrix", "MatrixFactorization"))
setMethod("determinant",
          c(x = .cl, logarithm = "missing"),
          function(x, logarithm = TRUE, ...)
              determinant(x, TRUE, ...))

setMethod("determinant",
          c(x = "denseSchur", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              if (is.complex(values <- x@values) == is.complex(x@x))
                  .mkDet(x = values, logarithm = logarithm)
              else if (is.complex(values))
                  .mkDet(x = values, logarithm = logarithm,
                         sign = if (Re(prod(values)) < 0) -1L else 1L)
              else stop("should never happen ...")
          })

setMethod("determinant",
          c(x = "denseLU", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(denseLU_determinant, x, logarithm))

setMethod("determinant",
          c(x = "denseBunchKaufman", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(denseBunchKaufman_determinant, x, logarithm))

setMethod("determinant",
          c(x = "denseCholesky", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(denseCholesky_determinant, x, logarithm))

setMethod("determinant",
          c(x = "sparseQR", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(sparseQR_determinant, x, logarithm))

setMethod("determinant",
          c(x = "sparseLU", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .Call(sparseLU_determinant, x, logarithm))

setMethod("determinant",
          c(x = "sparseCholesky", logarithm = "logical"),
          function(x, logarithm = TRUE, sqrt = TRUE, ...) {
              if (missing(sqrt)) {
                  w <- getOption("Matrix.warnSqrtDefault",
                                 .MatrixEnv[["warnSqrtDefault"]])
                  if (is.atomic(w) && length(w) == 1L &&
                      ((w.na <- is.na(w <- as.integer(w))) || w > 0L)) {
                      if (w.na)
                          on.exit(options(Matrix.warnSqrtDefault = 0L))
                      else if (w > 1L) {
                          oop <- options(warn = 2L)
                          on.exit(options(oop))
                      }
                      warning(gettextf("the default value of argument '%s' of the method for '%s' with signature c(%s=\"%s\", %s=\"%s\") may change from %s to %s as soon as the next release of Matrix; set '%s' when programming",
                                       "sqrt", "determinant", "x", "sparseCholesky", "logarithm", "logical", "TRUE", "FALSE", "sqrt"),
                              domain = NA)
                  }
              }
              .Call(sparseCholesky_determinant, x, logarithm, sqrt)
          })

setMethod("determinant",
          c(x = "generalMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if (d[1L] != d[2L])
                  stop("determinant of non-square matrix is undefined")
              repr <- .M.repr(x)
              switch(repr, "R" = x <- .tCRT(x), "T" = x <- .M2C(x))
              trf <-
                  if (repr == "n")
                      lu(x, warnSing = FALSE)
                  else lu(x, errSing = FALSE)
              if (isS4(trf))
                  determinant(trf, logarithm, ...)
              else .mkDet(if (anyNA(x)) NaN else -Inf, logarithm, 1L)
          })

setMethod("determinant",
          c(x = "symmetricMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              posdef <- substr(.M.class(x, 2L), 2L, 2L) == "p"
              repr <- .M.repr(x)
              switch(repr, "R" = x <- .tCRT(x), "T" = x <- .M2C(x))
              trf <-
                  if (repr == "n" || repr == "p") {
                      if (posdef)
                          tryCatch(Cholesky(x, perm = FALSE),
                                   error = function(e) BunchKaufman(x, warnSing = FALSE))
                      else BunchKaufman(x, warnSing = FALSE)
                  } else {
                      if (posdef)
                          tryCatch(Cholesky(x, perm = TRUE, LDL = FALSE, super = FALSE),
                                   error = function(e) NULL)
                      else tryCatch(Cholesky(x, perm = TRUE, LDL = TRUE, super = FALSE),
                                    error = function(e) lu(x, errSing = FALSE))
                  }
              if (isS4(trf))
                  determinant(trf, logarithm, ...)
              else .mkDet(if (anyNA(x)) NaN else -Inf, logarithm, 1L)
          })

for (.cl in c("triangularMatrix", "diagonalMatrix"))
setMethod("determinant",
          c(x = .cl, logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              if (x@diag == "N")
                  .mkDet(x = diag(x, names = FALSE), logarithm = logarithm)
              else .mkDet(0, logarithm, 1L))

setMethod("determinant",
          c(x = "indMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...) {
              d <- x@Dim
              if (d[1L] != d[2L])
                  stop("determinant of non-square matrix is undefined")
              if (anyDuplicated.default(perm <- x@perm))
                  .mkDet(-Inf, logarithm, 1L)
              else .mkDet(0, logarithm, signPerm(perm))
          })

setMethod("determinant",
          c(x = "pMatrix", logarithm = "logical"),
          function(x, logarithm = TRUE, ...)
              .mkDet(0, logarithm, signPerm(x@perm)))

rm(.cl)
