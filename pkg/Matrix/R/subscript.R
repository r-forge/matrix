## METHODS FOR GENERIC: [
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.subscript.error.ist <- function(i) {
    if(isS4(i))
        gettextf("invalid subscript class \"%s\"", class(i))
    else gettextf("invalid subscript type \"%s\"", typeof(i))
}

.subscript.recycle <- function(i, mn, pattern) {
    ## Return integer vector corresponding to [nl]sparseVector 'i'
    ## recycled to length 'mn' :
    if(length(i.i <- i@i) == 0L)
        integer(0L)
    else if((i.length <- i@length) >= mn) {
        if(i.length > mn) {
            if(mn < 0x1p+53) {
                if(i.i[length(i.i)] >= mn + 1)
                    i.i[i.i >= mn + 1] <- NA
            } else {
                if(i.i[length(i.i)] > mn)
                    i.i[i.i > mn] <- NA
            }
        }
        if(pattern) i.i else i.i[i@x]
    } else {
        r <- ceiling(mn / i.length)
        mn. <- r * i.length
        i.i <-
            if(mn. <= .Machine$integer.max)
                rep.int(as.integer(i.i), r) +
                    rep(seq.int(from = 0L,
                                by = as.integer(i.length),
                                length.out = r),
                        each = length(i.i))
            else if(i.i[length(i.i)] + (r - 1) * i.length <=
                    0x1p+53)
                rep.int(as.double(i.i), r) +
                    rep(seq.int(from = 0,
                                by = as.double(i.length),
                                length.out = r),
                        each = length(i.i))
            else stop("recycled [nl]sparseVector would have maximal index exceeding 2^53")
        if(pattern) {
            if(mn. > mn) i.i[      i.i <= mn] else i.i
        } else {
            if(mn. > mn) i.i[i@x & i.i <= mn] else i.i[i@x]
        }
    }
}

## x[i] where 'i' is NULL or any vector or sparseVector
.subscript.1ary <- function(x, i) {
    mn <- prod(x@Dim)
    if(is.null(i))
        i <- integer(0L)
    else if(isS4(i)) {
        if(!.isVector(i))
            stop(.subscript.error.ist(i), domain = NA)
        kind <- .V.kind(i)
        if((pattern <- kind == "n") || kind == "l") {
            ## [nl]sparseVector
            i <- .subscript.recycle(i, mn, pattern)
            return(..subscript.1ary(x, i, unsorted = !pattern && anyNA(i)))
        }
        i <- i@x
    }
    switch(typeof(i),
           double =
               {
                   r <- min(1, i, na.rm = TRUE)
                   if(r < 1)
                       i <- if(r <= -1)
                                seq_len(mn)[i] # FIXME
                            else i[i >= 1]
                   ..subscript.1ary(x, i)
               },
           integer =
               {
                   r <- min(1L, i, na.rm = TRUE)
                   if(r < 1L)
                       i <- if(r <= -1L)
                                seq_len(mn)[i] # FIXME
                            else i[i >= 1L]
                   ..subscript.1ary(x, i)
               },
           logical =
               {
                   if(length(i) && !is.na(a <- all(i)) && a) {
                       if((len <- length(i)) <= mn)
                           return(as.vector(x))
                       else return(c(as.vector(x), rep.int(NA, len - mn)))
                   }
                   .subscript.1ary(x, as(i, "sparseVector")) # recursively
               },
           character =
               {
                   rep.int(if(.hasSlot(x, "x")) x@x[NA_integer_] else NA,
                           length(i))
               },
           stop(.subscript.error.ist(i), domain = NA))
}

## x[i] where 'i' is vector of type "integer" or "double"
## with elements greater than or equal to 1 (or NA)
..subscript.1ary <- function(x, i,
                             shape = .M.shape(x),
                             repr = .M.repr(x),
                             unsorted = is.unsorted(i)) {
    if(!any(repr == c("C", "R", "T")))
        return(.Call(R_subscript_1ary, x, i))
    if(shape == "t" && x@diag != "N")
        x <- ..diagU2N(x)
    if(shape == "s" || repr == "R") {
        mn <- prod(d <- x@Dim)
        if(mn < 0x1p+53) {
            r <- max(mn, i, na.rm = TRUE)
            if(r >= mn + 1)
                i[i >= mn + 1] <- NA
        } else if(is.double(i)) {
            r <- max(0x1p+53, i, na.rm = TRUE)
            if(r > 0x1p+53) {
                if(any(i > 0x1p+53 && i <= mn, na.rm = TRUE))
                    ## could be avoided in C, which has 64-bit integers :
                    warning("subscripts exceeding 2^53 replaced with NA")
                i[i > 0x1p+53] <- NA
            }
        }
        i1s <- i - 1L
        m <- d[1L]
        i. <- as.integer(i1s %% m)
        if(shape == "s") {
            j. <- as.integer(i1s %/% m)
            op <- if(x@uplo == "U") `>` else `<`
            if(length(w <- which(op(i., j.)))) {
                i.. <- i.[w]
                j.. <- j.[w]
                if(repr == "R")
                    i.[w] <- j..
                if(mn > .Machine$integer.max)
                    m <- as.double(m)
                i[w] <- m * i.. + j.. + 1L
            }
        }
    }
    o <-
        if(repr == "R")
            order(i., i)
        else if(is.na(unsorted) || unsorted)
            sort.list(i)
        else return(.Call(R_subscript_1ary, x, i))
    if(is.unsorted(o)) {
        s <- .Call(R_subscript_1ary, x, i[o])
        s[o] <- s
        s
    } else .Call(R_subscript_1ary, x, i)
}

## x[i] where 'i' is any array or Matrix
.subscript.1ary.mat <- function(x, i) {
    if(isS4(i)) {
        if(!.isMatrix(i))
            stop(.subscript.error.ist(i), domain = NA)
        if((logic <- any(.M.kind(i) == c("n", "l"))) || i@Dim[2L] != 2L) {
            if(logic && all(i@Dim) && !is.na(a <- all(i)) && a) {
                x <- as.vector(x)
                if((len <- prod(i@Dim)) <= (mn <- length(x)))
                    return(x)
                else return(c(x, rep.int(NA, len - mn)))
            }
            v <- if(.isDense(i)) "vector" else "sparseVector"
            return(.subscript.1ary(x, as(i, v)))
        }
        i <- as(i, "matrix")
    } else if(is.logical(i) || length(di <- dim(i)) != 2L || di[2L] != 2L)
        return(.subscript.1ary(x, i))
    switch(typeof(i),
           double =,
           integer =
               {
                   storage.mode(i) <- "integer"
                   if(min(1L, i, na.rm = TRUE) < 1L)
                       stop("negative values are not allowed in a matrix subscript")
                   dx <- x@Dim
                   m <- dx[1L]
                   n <- dx[2L]
                   if(m == n) {
                       if(max(n, i, na.rm = TRUE) > n)
                           stop("subscript out of bounds")
                   } else {
                       if(max(m, i[, 1L], na.rm = TRUE) > m ||
                          max(n, i[, 2L], na.rm = TRUE) > n)
                           stop("subscript out of bounds")
                   }
                   ## * rows containing 0 are deleted
                   ## * rows containing NA are kept
                   ## * rows containing both 0 and NA are handled
                   ##   according to value in first column
                   if(is.na(a <- all(i. <- i[, 1L])) || !a)
                       i <- i[i. > 0L, , drop = FALSE]
                   if(!all(j. <- i[, 2L], na.rm = TRUE))
                       i <- i[j. > 0L, , drop = FALSE]
                   ..subscript.1ary.mat(x, i)
               },
           character =
               {
                   dnx <- dimnames(x)
                   m <- c(match(i[, 1L], dnx[[1L]]), match(i[, 2L], dnx[[2L]]))
                   dim(m) <- di
                   if(any(!rowSums(is.na(i)) & rowSums(is.na(m))))
                       ## error if character row contains zero NA and
                       ## integer row contains at least one NA,
                       ## indicating non-match that should not be ignored
                       stop("subscript out of bounds")
                   ..subscript.1ary.mat(x, m)
               },
           stop(.subscript.error.ist(i), domain = NA))
}

## x[i] where 'i' is a 2-column matrix of type "integer"
## with i[, 1L] in 1:m (or NA) and i[, 2L] in 1:n (or NA)
..subscript.1ary.mat <- function(x, i,
                                 shape = .M.shape(x),
                                 repr = .M.repr(x)) {
    if(!any(repr == c("C", "R", "T")))
        return(.Call(R_subscript_1ary_mat, x, i))
    if(shape == "t" && x@diag != "N")
        x <- ..diagU2N(x)
    i. <- i[, 1L]
    j. <- i[, 2L]
    if(shape == "s") {
        op <- if(x@uplo == "U") `>` else `<`
        if(length(w <- which(op(i., j.)))) {
            i[w, ] <- i[w, 2:1]
            i. <- i[, 1L]
            j. <- i[, 2L]
        }
    }
    o <-
        if(repr == "R") {
            if(anyNA(j.))
                i.[is.na(j.)] <- NA
            order(i., j.)
        } else {
            if(anyNA(i.))
                j.[is.na(i.)] <- NA
            order(j., i.)
        }
    if(is.unsorted(o)) {
        s <- .Call(R_subscript_1ary_mat, x, i[o, , drop = FALSE])
        s[o] <- s
        s
    } else .Call(R_subscript_1ary_mat, x, i)
}

## x[i, j, drop] where 'i' and 'j' are NULL or any vector
.subscript.2ary <- function(x, i, j, drop) {
    d <- x@Dim
    l <- list(if(missing(i)) NULL else if(is.null(i)) integer(0L) else i,
              if(missing(j)) NULL else if(is.null(j)) integer(0L) else j)
    for(pos in 1:2) {
        if(!is.null(k <- l[[pos]])) {
            l[pos] <- list(
                switch(typeof(k),
                       double =
                           {
                               r <- d[pos]
                               if(max(r, k, na.rm = TRUE) >= r + 1)
                                   stop("subscript out of bounds")
                               if(min(1, k, na.rm = TRUE) < 1)
                                   seq_len(r)[k]
                               else as.integer(k)
                           },
                       integer =
                           {
                               r <- d[pos]
                               if(max(r, k, na.rm = TRUE) > r)
                                   stop("subscript out of bounds")
                               if(min(1L, k, na.rm = TRUE) < 1L)
                                   seq_len(r)[k]
                               else k
                           },
                       logical =
                           {
                               r <- d[pos]
                               if(length(k) > r)
                                   stop("logical subscript too long")
                               if(length(k) && !is.na(a <- all(k)) && a)
                                   NULL
                               else seq_len(r)[k]
                           },
                       character =
                           {
                               if(length(k) == 0L)
                                   integer(0L)
                               else if(is.null(nms <- dimnames(x)[[pos]]) ||
                                       anyNA(k <- match(k, nms)))
                                   stop("subscript out of bounds")
                               else k
                           },
                       stop(.subscript.error.ist(k), domain = NA)))
        }
    }
    if(is.double(lengths(l, use.names = FALSE)))
        stop("dimensions cannot exceed 2^31-1")
    ..subscript.2ary(x, l[[1L]], l[[2L]], drop = drop[1L])
}

## x[i, j, drop] where 'i' and 'j' are vectors of type "integer"
## of length not exceeding 2^31-1 with 'i' in 1:m (or NA) and 'j'
## in 1:n (or NA) ... NULL => missing
..subscript.2ary <- function(x, i, j, drop) {
    if(is.null(i) && is.null(j))
        r <- x
    else {
        r <- .Call(R_subscript_2ary, x, i, j)
        dn <- dimnames(x)
        if(!(is.null(i) || is.null(rn <- dn[[1L]])))
            dn[1L] <- list(if(length(i)) rn[i] else NULL)
        if(!(is.null(j) || is.null(cn <- dn[[2L]])))
            dn[2L] <- list(if(length(j)) cn[j] else NULL)
        r@Dimnames <- dn
    }
    if((is.na(drop) || drop) && any(r@Dim == 1L)) drop(as(r, "matrix")) else r
}

setMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 ".", ".", ".", na),
                         .M.level = 2)
              if(na == 2L) {
                  ## x[]
                  x
              } else if(na == 3L) {
                  ## x[, ]
                  drop(x)
              } else {
                  ## x[, , ], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 ".", ".", "l", na),
                         .M.level = 2)
              if(na < 4L) {
                  ## x[drop=], x[, drop=], x[drop=, ]
                  x
              } else if(na == 4L) {
                  ## x[, , drop=], x[, drop=, ], x[drop=, , ]
                  if(is.na(drop <- drop[1L]) || drop) drop(x) else x
              } else {
                  ## x[, , , drop=], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "i", ".", ".", na),
                         .M.level = 2)
              if(na == 2L) {
                  ## x[i=]
                  .subscript.1ary(x, i)
              } else if(na == 3L) {
                  ## x[i=, ], x[, i=]
                  .subscript.2ary(x, i, , drop = TRUE)
              } else {
                  ## x[i=, , ], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "i", ".", "l", na),
                         .M.level = 2)
              if(na == 3L) {
                  ## x[i=, drop=]
                  .subscript.1ary(x, i)
              } else if(na == 4L) {
                  ## x[i=, , drop=], x[, i=, drop=]
                  .subscript.2ary(x, i, , drop = drop)
              } else {
                  ## x[i=, , , drop=], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 ".", "i", ".", na),
                         .M.level = 2)
              if(na == 2L) {
                  ## x[j=]
                  .subscript.1ary(x, j)
              } else if(na == 3L) {
                  ## x[j=, ], x[, j=]
                  .subscript.2ary(x, , j, drop = TRUE)
              } else {
                  ## x[, j=, ], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 ".", "i", "l", na),
                         .M.level = 2)
              if(na == 3L) {
                  ## x[j=, drop=]
                  .subscript.1ary(x, j)
              } else if(na == 4L) {
                  ## x[j=, , drop=], x[, j=, drop=]
                  .subscript.2ary(x, , j, drop = drop)
              } else {
                  ## x[, j=, , drop=], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "i", "i", ".", na),
                         .M.level = 2)
              if(na == 3L) {
                  ## x[i=, j=], x[j=, i=]
                  .subscript.2ary(x, i, j, drop = TRUE)
              } else {
                  ## x[i=, j=, ], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "i", "i", "l", na),
                         .M.level = 2)
              if(na == 4L) {
                  ## x[i=, j=, drop=], x[j=, i=, drop=]
                  .subscript.2ary(x, i, j, drop = drop)
              } else {
                  ## x[i=, j=, , drop=], etc.
                  stop("incorrect number of dimensions")
              }
          })

for(.cl in c("matrix", "nMatrix", "lMatrix"))
setMethod("[", signature(x = "Matrix", i = .cl, j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "m", ".", ".", na),
                         .M.level = 2)
              if(na == 2L) {
                  ## x[i=]
                  .subscript.1ary.mat(x, i)
              } else if(na == 3L) {
                  ## x[i=, ], x[, i=]
                  .subscript.2ary(x, i, , drop = TRUE)
              } else {
                  ## x[i=, , ], etc.
                  stop("incorrect number of dimensions")
              }
          })

setMethod("[", signature(x = "Matrix", i = "NULL", j = "ANY",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              i <- integer(0L)
              callGeneric()
          })

setMethod("[", signature(x = "Matrix", i = "ANY", j = "NULL",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              j <- integer(0L)
              callGeneric()
          })

setMethod("[", signature(x = "Matrix", i = "NULL", j = "NULL",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              i <- integer(0L)
              j <- integer(0L)
              callGeneric()
          })

if(FALSE) {
## MJ: unfinished
setMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if(nargs() != 2L)
                  stop("incorrect number of dimensions")
              mn <- x@length
              switch(typeof(i),
                     double = {
                         . # TODO
                     },
                     integer = {
                         . # TODO
                     },
                     logical = {
                         if(length(i) && !is.na(a <- all(i)) && a) {
                             if((len <- length(i)) <= mn)
                                 return(x)
                             else return(c(x, rep.int(NA, len - mn)))
                         }
                         x[.subscript.recycle(as(i, "sparseVector"), mn, FALSE)] # recursively
                     },
                     stop(.subscript.error.ist(i), domain = NA))
          })
} else {
setMethod("[", signature(x = "sparseVector", i = "index", j = "missing",
                         drop = "missing"),
          function (x, i, j, ..., drop) {
              has.x <- .hasSlot(x, "x") ## has "x" slot
              n <- x@length
              if(extends(cl.i <- getClass(class(i)), "numeric") && any(i < 0)) {
                  if(any(i > 0))
                      stop("you cannot mix negative and positive indices")
                  if(any(z <- i == 0))
                      i <- i[!z]
                  ## all (i < 0), negative indices:
                  ## want to remain sparse --> *not* using intIv()
                  ##
                  ## TODO: more efficient solution would use C ..
                  i <- unique(sort(-i)) # so we need to drop the 'i's
                  nom <- is.na(m <- match(x@i, i))
                  ## eliminate those non-0 which do match:
                  x@i <- x@i[nom]
                  if(has.x)
                      x@x <- x@x[nom]
                  ## now all x@i "appear in 'i' but must be adjusted
                  ## for the removals:
                  x@i <- x@i - findInterval(x@i, i)
                  x@length <- n - length(i)
              } else { ## i >= 0  or  non-numeric  'i'
                  ii <- intIv(i, n, cl.i=cl.i)
                  m <- match(x@i, ii, nomatch = 0)
                  sel <- m > 0L
                  x@length <- length(ii)
                  x@i <- m[sel]
                  if(any(iDup <- duplicated(ii))) {
                      i.i <- match(ii[iDup], ii)
                      jm <- lapply(i.i, function(.) which(. == m))
                      if (has.x)
                          sel <- c(which(sel), unlist(jm))
                      x@i <- c(x@i, rep.int(which(iDup), lengths(jm)))
                  }
                  if(doSort <- is.unsorted(x@i)) {
                      io <- order(x@i, method="radix")
                      x@i <- x@i[io]
                  }
                  if (has.x)
                      x@x <- if(doSort) x@x[sel][io] else x@x[sel]
              }
              x
          })
}

setMethod("[", signature(x = "sparseVector", i = "sparseVector", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if(nargs() != 2L)
                  stop("incorrect number of dimensions")
              kind <- .V.kind(i)
              if((pattern <- kind == "n") || kind == "l")
                  x[.subscript.recycle(i, x@length, pattern)]
              else x[i@x]
          })

setMethod("[", signature(x = "sparseVector", i = "NULL", j = "ANY",
                         drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              i <- integer(0L)
              callGeneric()
          })


## METHODS FOR GENERIC: head
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("head", signature(x = "Matrix"),
          head.matrix)

setMethod("head", signature(x = "sparseVector"),
          function(x, n = 6L, ...) {
              stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
              len <- x@length
              n <- if(n < 0L) max(len + n, 0L) else min(n, len)
              if(n >= len)
                  return(x)
              nnz <- length(i <- x@i)
              x@length <- n <- if(is.integer(i)) as.integer(n) else trunc(n)
              if(nnz > 0L && i[nnz] > n) {
                  if(i[1L] > n) {
                      x@i <- integer(0L)
                      if(.hasSlot(x, "x"))
                          x@x <- vector(typeof(x@x), 0L)
                  } else {
                      ii <- 1L:(which.max(i > n) - 1L)
                      x@i <- i[ii]
                      if(.hasSlot(x, "x"))
                          x@x <- x@x[ii]
                  }
              }
              x
          })


## METHODS FOR GENERIC: tail
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("tail", signature(x = "Matrix"),
          tail.matrix)

setMethod("tail", signature(x = "sparseVector"),
          function(x, n = 6L, ...) {
              stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
              len <- x@length
              n <- if(n < 0L) max(len + n, 0L) else min(n, len)
              if(n >= len)
                  return(x)
              nnz <- length(i <- x@i)
              x@length <- n <- if(is.integer(i)) as.integer(n) else trunc(n)
              if(nnz > 0L && i[1L] <= (k <- len - n)) {
                  if(i[nnz] <= k) {
                      x@i <- integer(0L)
                      if(.hasSlot(x, "x"))
                          x@x <- vector(typeof(x@x), 0L)
                  } else {
                      ii <- which.min(i <= k):nnz
                      x@i <- i[ii] - k
                      if(.hasSlot(x, "x"))
                          x@x <- x@x[ii]
                  }
              }
              x
          })
