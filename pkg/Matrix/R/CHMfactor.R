for(.cl in c("dtCMatrix", "dsparseMatrix", "triangularMatrix", "CsparseMatrix",
             "dMatrix", "sparseMatrix", "Matrix"))
setAs("CHMfactor", .cl,
      function(from) .Call(CHMfactor_to_sparse, from))
rm(.cl)

setAs("CHMfactor", "RsparseMatrix",
      function(from) .CR2RC(.Call(CHMfactor_to_sparse, from)))
setAs("CHMfactor", "TsparseMatrix",
      function(from) .CR2T(.Call(CHMfactor_to_sparse, from)))
setAs("CHMfactor", "pMatrix",
      function(from) {
          r <- new("pMatrix")
          r@Dim <- d <- from@Dim
          r@perm <- if(length(perm <- from@perm)) perm + 1L else seq_len(d[1L])
          r
      })

## returning list(P1', L, L', P1) or list(P1', L1, D, L1', P1),
## where  A = P1' L L' P1 = P1' L1 D L1' P1  and  L = L1 sqrt(D)
setMethod("expand2", signature(x = "CHMfactor"),
          function(x, LDL = TRUE, ...) {
              d <- x@Dim
              dn <- x@Dimnames
              perm <- x@perm
              perm <- if(length(perm))
                          invertPerm(perm, 0L, 1L)
                      else seq_len(d[1L])
              P <- new("pMatrix",
                       Dim = d,
                       Dimnames = c(list(NULL), dn[2L]),
                       margin = 2L,
                       perm = perm)
              P. <- P
              P.@Dimnames <- c(dn[1L], list(NULL))
              P.@margin <- 1L
              L <- .Call(CHMfactor_to_sparse, x) # nonunit lower triangular
              if(LDL) {
                  ## FIXME: be faster if isLDL(x) is TRUE
                  L.ii <- diag(L, names = FALSE)
                  L.p <- L@p
                  L@x <- L@x / rep.int(L.ii, L.p[-1L] - L.p[-length(L.p)])
                  L <- ..diagN2U(L, sparse = TRUE)
                  D <- new("ddiMatrix")
                  D@Dim <- d
                  D@x <- L.ii * L.ii
                  list(P1. = P., L1 = L, D = D, L1. = t(L), P1 = P)
              } else list(P1. = P., L = L, L. = t(L), P1 = P)
          })

## returning list(P, L), where A = P' L L' P
## MJ: for backwards compatibility
setMethod("expand", signature(x = "CHMfactor"),
          function(x, ...)
              list(P = as(x, "pMatrix"), L = as(x, "dtCMatrix")))

setMethod("update", signature(object = "CHMfactor"),
	  function(object, parent, mult = 0, ...) {
              cld <- getClassDef(class(parent))
              stopifnot(extends(cld, "sparseMatrix"))
              if(!extends(cld, "CsparseMatrix")) {
                  parent <- as(parent, "CsparseMatrix")
                  cld <- getClassDef(class(parent))
              }
              if(!extends(cld, "dMatrix")) {
                  parent <- ..sparse2d(parent)
                  cld <- getClassDef(class(parent))
              }
              if(!extends(cld, "symmetricMatrix") && Matrix.verbose() >= 1)
                  message("'parent' is not formally symmetric, will be handled as 'tcrossprod(parent)'")
              chkDots(..., which.call = -2L)
              .Call(CHMfactor_update, object, parent, mult)
          })

## Exported fast version, somewhat hidden;
## here 'parent' _must_ inherit from d[gs]CMatrix
.updateCHMfactor <- function(object, parent, mult = 0)
    .Call(CHMfactor_update, object, parent, mult)

setMethod("updown",
          signature(update = "logical", C = "mMatrix", L = "CHMfactor"),
          function(update, C, L) {
              if(length(update) != 1L || is.na(update))
                  stop("'update' must be TRUE, FALSE, \"+\", or \"-\"")
              bnew <- as(L, "pMatrix") %*% C
              .Call(CHMfactor_updown,
                    update, as(bnew, "CsparseMatrix"), L)
          })

setMethod("updown",
          signature(update = "character", C = "mMatrix", L = "CHMfactor"),
          function(update, C, L){
              if(length(update) != 1L || is.na(update) ||
                 (update != "+" && update != "-"))
                  stop("'update' must be TRUE, FALSE, \"+\", or \"-\"")
              bnew <- as(L, "pMatrix") %*% C
              .Call(CHMfactor_updown,
                    update == "+", as(bnew, "CsparseMatrix"), L)
          })

## "Fallback" giving a "good" error message
setMethod("updown", signature(update = "ANY", C = "ANY", L = "ANY"),
	  function(update, C, L) stop("'update' must be TRUE, FALSE, \"+\", or \"-\"; 'C' a [mM]atrix; and 'L' a CHMfactor"))

##' Test whether a CHMfactor object is LDL or LL
##' @param x a CHMfactor object
##' @return TRUE if 'x' is LDL, otherwise FALSE
isLDL <- function(x)
{
    if(!is(x, "CHMfactor"))
        stop("'x' does not inherit from virtual class CHMfactor")
    !x@type[2L] # '!'<=>'not' as type[2L] := (cholmod_factor)->is_ll
}
.isLDL <- function(x)
    !x@type[2L]

## Currently not exported, but called from some examples with ':::'
ldetL2up <- function(x, parent, Imult)
{
    ## Purpose: compute  log Det |A + m*I|  for many values of m
    ## ----------------------------------------------------------------------
    ## Arguments: x: CHMfactor to be updated
    ##      parent : CsparseMatrix M; for symmetric M, A = M, otherwise A = MM'
    ##       Imult : a numeric *vector* of 'm's (= I multipliers)
    ## ----------------------------------------------------------------------
    ## Author: Doug Bates, Date: 19 Mar 2008

    stopifnot(is(x, "CHMfactor"), is(parent, "CsparseMatrix"),
              nrow(x) == nrow(parent))
    .Call(CHMfactor_ldetL2up, x, parent, as.double(Imult))
}

## MJ: unused
if(FALSE) {
##' Update a sparse Cholesky factorization in place
##' @param L A sparse Cholesky factor that inherits from CHMfactor
##' @param parent a sparse matrix for updating the factor.  Either a
##'   dsCMatrix, in which case L is updated to the Cholesky
##'   factorization of parent, or a dgCMatrix, in which case L is
##'   updated to the Cholesky factorization of tcrossprod(parent)
##' @param Imult an optional positive scalar to be added to the
##'   diagonal before factorization,
##' @return NULL.  This function always returns NULL.  It is called
##'   for its side-effect of updating L in place.
##' @note This function violates the functional language semantics of
##'   R in that it updates its argument L in place (i.e. without copying).
##'   This is intentional but it means the function should be used
##'   with caution.  If the preceding sentences do not make sense to
##'   you, you should not use this function,.
destructive_Chol_update <- function(L, parent, Imult = 1)
{
    stopifnot(is(L, "CHMfactor"), is(parent, "sparseMatrix"))
    .Call(destructive_CHM_update, L, parent, Imult)
}
} ## MJ
