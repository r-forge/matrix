## METHODS FOR CLASS: indMatrix
## index matrices, i.e., matrices with standard unit vectors
## for all rows _or_ all columns
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.perm2ind <- function(perm, n, margin = 1L, check.p = 0L) {
    if(mode(perm) != "numeric")
        stop(gettextf("'%s' is not of type \"%s\" or \"%s\"",
                      "perm", "integer", "double"),
             domain = NA)
    else if((m <- length(perm)) == 0L)
        perm <- integer(0L)
    else if(anyNA(r <- range(perm)))
        stop(gettextf("'%s' contains NA", "perm"),
             domain = NA)
    else if(r[1L] < 1L)
        stop(gettextf("'%s' has elements less than %d", "perm", 1L),
             domain = NA)
    else if(m > .Machine$integer.max ||
            (is.double(perm) && trunc(r[2L]) > .Machine$integer.max))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    else { perm <- as.integer(perm); r <- as.integer(r) }

    if(m.n <- missing(n))
        n <- if(m == 0L) 0L else r[2L]
    else if(mode(n) != "numeric" || length(n) != 1L || is.na(n) || n < 0L)
        stop(gettextf("'%s' is not a non-negative number", "n"),
             domain = NA)
    else if(is.double(n) && trunc(n) > .Machine$integer.max)
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    else if(r[2L] > as.integer(n))
        stop(gettextf("'%s' has elements exceeding '%s'", "perm", "n"),
             domain = NA)
    else n <- as.integer(n)

    if(mode(margin) != "numeric" || length(margin) != 1L || is.na(margin) ||
       (margin != 1L && margin != 2L))
        stop(gettextf("'%s' is not %d or %d", "margin", 1L, 2L),
             domain = NA)

    give.p <- check.p >= 1L && m == n &&
        (m == 0L || (all(r == c(1L, m)) && !anyDuplicated.default(perm)))
    if(check.p >= 2L && !give.p)
        stop(gettextf("'%s' is not a permutation of seq_len(%s)",
                      "perm", if(m.n) "max(perm, 0)" else "n"),
             domain = NA)

    J <- new(if(give.p) "pMatrix" else "indMatrix")
    nms <- names(perm)
    if(margin == 1L) {
        J@Dim <- c(m, n)
        J@Dimnames = list(nms, if(give.p) nms)
    } else {
        J@Dim <- c(n, m)
        J@Dimnames = list(if(give.p) nms, nms)
        J@margin <- 2L
    }
    J@perm <- perm
    J
}

setAs("numeric", "indMatrix",
      function(from) .perm2ind(from))

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
          stop("matrix must have exactly one entry in each row or column")
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
      function(from) .perm2ind(from, check.p = 2L))

setAs("nsparseMatrix", "pMatrix",
      function(from) {
          d <- from@Dim
          if((n <- d[1L]) != d[2L])
              stop(gettextf("attempt to coerce non-square matrix to %s",
                            "pMatrix"),
                   domain = NA)
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
          stop("matrix must have exactly one entry in each row and column")
      })

setAs("indMatrix", "pMatrix",
      function(from) new("pMatrix", from))
