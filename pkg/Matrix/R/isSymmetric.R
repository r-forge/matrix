## METHODS FOR GENERIC: isSymmetric
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("isSymmetric" , c(object = "denseMatrix"),
          function(object,
                   tol = 100 * .Machine$double.eps, tol1 = 8 * tol,
                   trans = "C", checkDN = TRUE, ...) {
              if (checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  checkDN <- ca(...)
              }
              stopifnot(is.numeric(tol), length(tol) == 1L, !is.na(tol))
              ans <- .Call(R_dense_is_symmetric, object,
                           trans, tol <= 0, checkDN)
              if (!is.na(ans))
                  return(ans)
              ## 'object' is an n-by-n [dz]denseMatrix, n >= 1
              ae <- function(target, current,
                             tolerance, scale = NULL, ...)
                  all.equal.numeric(target = target, current = current,
                                    tolerance = tolerance,
                                    scale = scale,
                                    check.attributes = FALSE,
                                    check.class = FALSE)
              conjugate <- is.complex(object@x) && identical(trans, "C")
              if (length(tol1) && (n <- object@Dim[1L]) > 1L) {
                  op <- if (conjugate) Conj else identity
                  for (i in if (n > 4L) c(1L, 2L, n - 1L, n) else 1L:n)
                      if (is.character(ae(target = object[i, ],
                                          current = op(object[, i]),
                                          tolerance = tol1, ...)))
                          return(FALSE)
              }
              object <- .M2gen(object)
              op <- if (conjugate) ct else t
              isTRUE(ae(target = object@x,
                        current = op(object)@x,
                        tolerance = tol, ...))
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("isSymmetric" , c(object = .cl),
          function(object,
                   tol = 100 * .Machine$double.eps,
                   trans = "C", checkDN = TRUE, ...) {
              if (checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  checkDN <- ca(...)
              }
              stopifnot(is.numeric(tol), length(tol) == 1L, !is.na(tol))
              ans <- .Call(R_sparse_is_symmetric, object,
                           trans, tol <= 0, checkDN)
              if (!is.na(ans))
                  return(ans)
              ## 'object' is an n-by-n [dz]sparseMatrix, n >= 1
              ae <- function(target, current,
                             tolerance, scale = NULL, ...)
                  .V.a.e(target = target, current = current,
                         tolerance = tolerance, scale = scale,
                         check.attributes = FALSE, check.class = FALSE)
              conjugate <- is.complex(object@x) && identical(trans, "C")
              op <- if (conjugate) ct else t
              isTRUE(ae(target = .M2V(object),
                        current = .M2V(op(object)),
                        tolerance = tol, ...))
          })
rm(.cl)

setMethod("isSymmetric", c(object = "diagonalMatrix"),
          ## needs care w/ dimnames() and for <complex> :
          function(object,
                   tol = 100 * .Machine$double.eps,
                   trans = "C", checkDN = TRUE, ...) {
              if (checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  if (ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              !is.complex(x <- object@x) || # check more only if <complex>
                  !identical(trans, "C") ||
                  object@diag != "N" ||
                  {
                      ae <- function(target, current,
                                     tolerance, scale = NULL, ...)
                          all.equal.numeric(target = target,
                                            current = current,
                                            tolerance = tolerance,
                                            scale = scale,
                                            check.attributes = FALSE,
                                            check.class = FALSE)
                      isTRUE(ae(x, Conj(x), tolerance = tol, ...))
                  }
          })

setMethod("isSymmetric", c(object = "indMatrix"),
          function(object, checkDN = TRUE, ...) {
              d <- object@Dim
              if ((n <- d[2L]) != d[1L])
                  return(FALSE)
              if (checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  if (ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              perm <- object@perm
              all(perm[perm] == seq_len(n))
          })
