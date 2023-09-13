## METHODS FOR CLASS: indMatrix
## index matrices, i.e., matrices with standard unit vectors
## for all rows _or_ all columns
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.perm2ind <- function(perm, n, margin = 1L, check.p = FALSE) {
    perm.i <- perm
    if(!is.numeric(perm))
        stop("'perm' must be numeric")
    else if(anyNA(r <- range(perm)) || r[1L] < 1L ||
            (is.double(perm) && any(perm != (perm.i <- as.integer(perm)))))
        stop("elements of 'perm' must be positive integers")
    else if((m.i <- length(perm)) > (M <- .Machine$integer.max) || r[2L] > M)
        stop("dimensions cannot exceed 2^31-1")

    if(missing(n))
        n.i <- as.integer(r[2L])
    else {
        n.i <- n
        if(!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0L ||
           (is.double(n) && n != (n.i <- as.integer(n))))
            stop("'n' must be a non-negative integer")
        else if(n > M)
            stop("dimensions cannot exceed 2^31-1")
        else if(r[2L] > n)
            stop("elements of 'perm' cannot exceed 'n'")
    }

    if(!is.numeric(margin) || length(margin) != 1L || is.na(margin) ||
       (margin != 1L && margin != 2L))
        stop("'margin' must be 1 or 2")
    margin.i <- as.integer(margin)

    give.p <- check.p && m.i == n.i &&
        (m.i == 0L || (all(r == c(1L, m.i)) && !anyDuplicated.default(perm.i)))

    J <- new(if(give.p) "pMatrix" else "indMatrix")
    nms <- names(perm)
    if(margin.i == 1L) {
        J@Dim <- c(m.i, n.i)
        J@Dimnames = list(nms, if(give.p) nms)
    } else {
        J@Dim <- c(n.i, m.i)
        J@Dimnames = list(if(give.p) nms, nms)
        J@margin <- 2L
    }
    J@perm <- perm.i
    J
}

setAs("numeric", "indMatrix",
      function(from) {
          J <- new("indMatrix")
          if((m <- length(from)) == 0L)
              return(J)
          from.i <- from
          if(anyNA(r <- range(from)) || r[1L] < 1L ||
             (is.double(from) && any(from != (from.i <- as.integer(from)))))
              stop("elements of 'perm' slot must be positive integers")
          if(m > (M <- .Machine$integer.max) || r[2L] > M)
              stop("dimensions cannot exceed 2^31-1")
          J@Dim <- c(m, as.integer(r[2L]))
          J@Dimnames <- list(names(from), NULL)
          J@perm <- from.i
          J
      })

## FIXME: deprecate this method and export more general function .perm2ind
setAs("list", "indMatrix",
      function(from) do.call(.perm2ind, unname(from)))

setAs("nsparseMatrix", "indMatrix",
      function(from) {
          from <- .M2gen(from)
          J <- new("indMatrix")
          J@Dim <- from@Dim
          J@Dimnames <- from@Dimnames
          from. <- .M2R(from)
          p <- from.@p
          m <- length(p) - 1L
          if(all(p == 0:m)) {
              J@perm <- from.@j + 1L
              return(J)
          }
          from. <- .M2C(from)
          p <- from.@p
          n <- length(p) - 1L
          if(all(p == 0:n)) {
              J@perm <- from.@i + 1L
              J@margin <- 2L
              return(J)
          }
          stop("matrix must have exactly one nonzero element in each row or column")
      })

setMethod("band", signature(x = "indMatrix"),
          function(x, k1, k2, ...) band(.M2kind(x, "n"), k1, k2))

setMethod("triu", signature(x = "indMatrix"),
          function(x, k = 0L, ...) triu(.M2kind(x, "n")))

setMethod("tril", signature(x = "indMatrix"),
          function(x, k = 0L, ...) tril(.M2kind(x, "n")))

setMethod("diag", signature(x = "indMatrix"),
          function(x, nrow, ncol, names = TRUE) {
              if((m <- min(x@Dim)) == 0L)
                  return(logical(0L))
              i <- seq_len(m)
              r <- x@perm[i] == i
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 identical(nms <- dn[[1L]][i], dn[[2L]][i]))
                  names(r) <- nms
              r
          })

setMethod("diag<-", signature(x = "indMatrix"),
          function(x, value) `diag<-`(.M2kind(x, "n"), value))

setMethod("t", signature(x = "indMatrix"),
          function(x) {
              r <- new("indMatrix")
              r@Dim <- x@Dim[2:1]
              r@Dimnames = x@Dimnames[2:1]
              r@perm <- x@perm
              if(x@margin == 1L)
                  r@margin <- 2L
              r
          })

setMethod("forceSymmetric", signature(x = "indMatrix", uplo = "missing"),
          function(x, uplo) forceSymmetric(.M2kind(x, "n")))

setMethod("forceSymmetric", signature(x = "indMatrix", uplo = "character"),
          function(x, uplo) forceSymmetric(.M2kind(x, "n"), uplo))

setMethod("symmpart", signature(x = "indMatrix"),
          function(x) symmpart(.M2kind(x, "d")))

setMethod("skewpart", signature(x = "indMatrix"),
          function(x) skewpart(.M2kind(x, "d")))

setMethod("isDiagonal", signature(object = "indMatrix"),
          function(object) {
              d <- object@Dim
              if((n <- d[2L]) != d[1L])
                  return(FALSE)
              all(object@perm == seq_len(n))
          })

setMethod("isTriangular", signature(object = "indMatrix"),
          function(object, upper = NA, ...) {
              d <- object@Dim
              if((n <- d[2L]) != d[1L])
                  return(FALSE)
              if(object@margin == 1L) {
                  i <- seq_len(n)
                  j <- object@perm
              } else {
                  i <- object@perm
                  j <- seq_len(n)
              }
              if(is.na(upper)) {
                  if(all(j >= i))
                      return(`attr<-`(TRUE, "kind", "U"))
                  if(all(i <= j))
                      return(`attr<-`(TRUE, "kind", "L"))
                  FALSE
              } else if(upper) {
                  all(j >= i)
              } else {
                  all(i <= j)
              }
          })

setMethod("isSymmetric", signature(object = "indMatrix"),
          function(object, checkDN = TRUE, ...) {
              d <- object@Dim
              if((n <- d[2L]) != d[1L])
                  return(FALSE)
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...) check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              perm <- object@perm
              all(perm[perm] == seq_len(n))
          })



## METHODS FOR CLASS: pMatrix
## permutation matrices, i.e., matrices with standard unit vectors
## for all rows _and_ all columns
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MJ: could export without dot
.changeMargin <- function(x) {
    x@margin <- if(x@margin == 1L) 2L else 1L
    x@perm <- invertPerm(x@perm)
    x
}

setAs("numeric", "pMatrix",
      function(from) {
          J <- new("pMatrix")
          if((m <- length(from)) == 0L)
              return(J)
          if(m > .Machine$integer.max)
              stop("dimensions cannot exceed 2^31-1")
          from.i <- from
          if(anyNA(r <- range(from)) || any(r != c(1L, m)) ||
             (is.double(from) && any(from != (from.i <- as.integer(from)))) ||
             anyDuplicated.default(from.i))
              stop("'perm' slot must be a permutation of seq_along(perm)")
          nms <- names(from)
          J@Dim <- c(m, m)
          J@Dimnames <- list(nms, nms)
          J@perm <- from.i
          J
      })

setAs("nsparseMatrix", "pMatrix",
      function(from) {
          d <- from@Dim
          if((n <- d[1L]) != d[2L])
              stop("attempt to coerce non-square matrix to pMatrix")
          from <- .M2gen(from)
          J <- new("pMatrix")
          J@Dim <- d
          J@Dimnames <- from@Dimnames
          from. <- .M2R(from)
          p <- from.@p
          m <- length(p) - 1L
          if(all(p == 0:m) && !anyDuplicated.default(j <- from.@j)) {
              J@perm <- j + 1L
              return(J)
          }
          from. <- .M2C(from)
          p <- from.@p
          n <- length(p) - 1L
          if(all(p == 0:n) && !anyDuplicated.default(i <- from.@i)) {
              J@perm <- i + 1L
              J@margin <- 2L
              return(J)
          }
          stop("matrix must have exactly one nonzero element in each row and column")
      })

setAs("indMatrix", "pMatrix",
      function(from) new("pMatrix", from))

setMethod("t", signature(x = "pMatrix"),
          function(x) {
              r <- new("pMatrix")
              r@Dim <- x@Dim
              r@Dimnames = x@Dimnames[2:1]
              r@perm <- x@perm
              if(x@margin == 1L)
                  r@margin <- 2L
              r
          })
