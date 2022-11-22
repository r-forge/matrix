## METHODS FOR GENERIC: [                    UNFINISHED AND NOT-YET-USED
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## GOAL: automate method definitions and eventually replace ones in
##
##             ./Csparse.R
##             ./Matrix.R
##             ./Tsparse.R
##             ./denseMatrix.R
##             ./diagMatrix.R
##             ./indMatrix.R
##             ./packedMatrix.R
##             ./sparseMatrix.R
##
##       need to write C-level functions
##
##             *_subscript_1ary    (x, i         )
##             *_subscript_1ary_mat(x, i         )
##             *_subscript_2ary    (x, i, j, drop)
##
##       for * = unpackedMatrix,packedMatrix,
##               CsparseMatrix,RsparseMatrix,TsparseMatrix,
##               diagonalMatrix,indMatrix

.subscript.error.dim <- "incorrect number of dimensions"
.subscript.error.oob <- "subscript out of bounds"
.subscript.error.neg <- "negative values are not allowed in a matrix subscript"
.subscript.error.lng <- "logical subscript too long"
.subscript.error.ist <- function(i) {
    if(isS4(i))
        gettextf("invalid subscript class \"%s\"", class(i))
    else
        gettextf("invalid subscript type \"%s\"", typeof(i))
}

## x[i] where 'i' is NULL or any vector or sparseVector
.subscript.1ary <- function(x, i) {
    mn <- prod(x@Dim)
    if(is.null(i))
        i <- integer(0L)
    else if(isS4(i)) {
        if(!is(i, "sparseVector"))
            stop(.subscript.error.ist(i), domain = NA)
        if(!(has.x <- .hasSlot(i, "x")) || is.logical(i.x <- i@x)) {
            ## [nl]sparseVector
            i.length <- i@length
            i.i <- i@i
            i <-
                if(i.length >= mn) {
                    if(i.length > mn && i.i[length(i.i)] >= mn + 1)
                        i.i[i.i >= mn + 1] <- NA
                    if(has.x) i.i[i.x] else i.i
                } else {
                    r <- ceiling(mn / i.length)
                    i.i <-
                        if(r * i.length <= .Machine$integer.max)
                            rep.int(as.integer(i.i), r) +
                                rep(seq.int(from = 0L,
                                            by = as.integer(i.length),
                                            length.out = r),
                                    each = length(i.i))
                        else {
                            rep.int(as.double(i.i), r) +
                                rep(seq.int(from = 0,
                                            by = as.double(i.length),
                                            length.out = r),
                                    each = length(i.i))
                        }
                    if(has.x) {
                        if(r * i.length > mn) i.i[i.x & i.i <= mn] else i.i[i.x]
                    } else {
                        if(r * i.length > mn) i.i[      i.i <= mn] else i.i
                    }
                }
            return(..subscript.1ary(x, i, unsorted = FALSE))
        }
        i <- i.x
    }
    switch(typeof(i),
           double =
               {
                   m <- min(1, i, na.rm = TRUE)
                   if(m < 1)
                       i <- if(m <= -1)
                                seq_len(mn)[i] # FIXME
                            else {
                                if(is.object(i))
                                    i <- as.double(i)
                                i[i >= 1]
                            }
                   ..subscript.1ary(x, i)
               },
           integer =
               {
                   m <- min(1L, i, na.rm = TRUE)
                   if(m < 1L)
                       i <- if(m <= -1L)
                                seq_len(mn)[i] # FIXME
                            else {
                                if(is.object(i))
                                    i <- as.integer(i)
                                i[i >= 1L]
                            }
                   ..subscript.1ary(x, i)
               },
           logical =
               {
                   if(!is.na(a <- all(i)) && a) {
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
    if(!nzchar(repr))
        return(.Call(R_subscript_1ary, x, i))
    if(repr == "R" || shape == "s") {
        d <- x@Dim
        m <- d[1L]
        mn. <- min(prod(d), 0x1p+52) # R_XLEN_T_MAX
        if(max(mn., i, na.rm = TRUE) >= mn. + 1)
            i[i >= mn. + 1] <- NA
        i1s <- i - 1L
        i. <- as.integer(i1s %% m)
        if(shape == "s") {
            j. <- as.integer(i1s %/% m)
            op <- if(x@uplo == "U") `>` else `<`
            if(length(w <- which(op(i., j.)))) {
                tmp <- i.[w]
                i.[w] <- j.[w]
                j.[w] <- tmp
            }
            if (mn. > .Machine$integer.max)
                m <- as.double(m)
            i <- m * j. + i. + 1L
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
        s[o] <- tmp
        s
    } else .Call(R_subscript_1ary, x, i)
}

## x[i] where 'i' is any array or Matrix
.subscript.1ary.mat <- function(x, i) {
    if(isS4(i)) {
        cld <- getClassDef(class(i))
        if(!extends(cld, "Matrix"))
            stop(.subscript.error.ist(i), domain = NA)
        logic <-
            extends(cld, "nMatrix") ||
            extends(cld, "lMatrix") ||
            extends(cld, "indMatrix")
        if(logic || i@Dim[2L] != 2L) {
            if(logic && !is.na(a <- all(i)) && a) {
                x <- as.vector(x)
                if((len <- prod(i@Dim)) <= (mn <- prod(x@Dim)))
                    return(x)
                else return(c(x, rep.int(NA, len - mn)))
            }
            v <- if(extends(cld, "denseMatrix")) "vector" else "sparseVector"
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
                       stop(.subscript.error.neg)
                   dx <- x@Dim
                   m <- dx[1L]
                   n <- dx[2L]
                   if(m == n) {
                       if(max(n, i, na.rm = TRUE) > n)
                           stop(.subscript.error.oob)
                   } else {
                       if(max(m, i[, 1L], na.rm = TRUE) > m ||
                          max(n, i[, 2L], na.rm = TRUE) > n)
                           stop(.subscript.error.oob)
                   }
                   ## * rows containing 0 are deleted
                   ## * rows containing NA are kept
                   ## * rows containing both 0 and NA are handled
                   ##   according to value in first column
                   if(is.na(a <- all(i. <- i[, 1L])) || !a)
                       i <- i[!i., , drop = FALSE]
                   if(!all(j. <- i[, 2L], na.rm = TRUE))
                       i <- i[!j., , drop = FALSE]
                   .Call(R_subscript_1ary_mat, x, i)
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
                       stop(.subscript.error.oob)
                   .Call(R_subscript_1ary_mat, x, m)
               },
           stop(.subscript.error.ist(i), domain = NA))
}

## x[i] where 'i' is a 2-column matrix of type "integer"
## with i[, 1L] in 1:m (or NA) and i[, 2L] in 1:n (or NA)
..subscript.1ary.mat <- function(x, i,
                                 shape = .M.shape(x),
                                 repr = .M.repr(x)) {
    if(!nzchar(repr))
        return(.Call(R_subscript_1ary_mat, x, i))
    i. <- i[, 1L]
    j. <- i[, 2L]
    if(shape == "s") {
        op <- if(x@uplo == "U") `>` else `<`
        if(length(w <- which(op(i., j.)))) {
            tmp <- i.[w]
            i.[w] <- j.[w]
            j.[w] <- tmp
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
    i <- list(if(missing(i)) NULL else if(is.null(i)) integer(0L) else i,
              if(missing(j)) NULL else if(is.null(j)) integer(0L) else j)
    for(pos in 1:2) {
        if(!is.null(k <- i[[pos]])) {
            i[[pos]] <-
                switch(typeof(k),
                       double =
                           {
                               r <- d[pos]
                               if(is.object(k))
                                   k <- as.double(k)
                               if(max(r, k, na.rm = TRUE) >= r + 1)
                                   stop(.subscript.error.oob)
                               seq_len(r)[k]
                           },
                       integer =
                           {
                               r <- d[pos]
                               if(is.object(k))
                                   k <- as.integer(k)
                               if(max(r, k, na.rm = TRUE) > r)
                                   stop(.subscript.error.oob)
                               if(min(1L, k, na.rm = TRUE) < 1L)
                                   seq_len(r)[k]
                               else r
                           },
                       logical =
                           {
                               r <- d[pos]
                               if(length(k) > r)
                                   stop(.subscript.error.lng)
                               if(!is.na(a <- all(k)) && a)
                                   NULL
                               else seq_len(r)[k]
                           },
                       character =
                           {
                               nms <- dimnames(x)[[pos]]
                               if(is.null(nms) || anyNA(k <- match(k, nms)))
                                   stop(.subscript.error.oob)
                               k
                           },
                       stop(.subscript.error.ist(k), domain = NA))
        }
    }
    ..subscript.2ary(x, i, drop = is.na(drop <- drop[1L]) && drop)
}

## x[i[[1L]], i[[2L]], drop] where 'i' is a length-2 list
## of vectors of type "integer" with i[[1L]] in 1:m (or NA)
## and i[[2L]] in 1:n (or NA) ... NULL indicates missingness
..subscript.2ary <- function(x, i, drop) {
    mi. <- is.null(i. <- i[[1L]])
    mj. <- is.null(j. <- j[[2L]])
    if(mi. && mj.) {

    } else if (mi.) {

    } else if (mj.) {

    } else {

    }
}

setMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "missing",
                         drop = "logical"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "missing",
                         drop = "logical"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
                         drop = "missing"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "missing", j = "index",
                         drop = "logical"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "missing"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "i", "i", ".", na),
                         .M.level = 2)
	      if(na == 3L) {
                  ## x[i=, j=], x[j=, i=]
                  .subscript.2ary(x, i, j, drop = TRUE)
              } else {
                  ## x[i=, j=, ], etc.
                  stop(.subscript.error.dim)
              }
          })

setMethod("[", signature(x = "Matrix", i = "index", j = "index",
                         drop = "logical"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("M[%s%s%s] : nargs() = %d",
                                 "i", "i", "l", na),
                         .M.level = 2)
              if(na == 4L) {
                  ## x[i=, j=, drop=], x[j=, i=, drop=]
                  .subscript.2ary(x, i, j, drop = drop)
              } else {
                  ## x[i=, j=, , drop=], etc.
                  stop(.subscript.error.dim)
              }
          })

for(.cl in c("matrix", "nMatrix", "lMatrix"))
setMethod("[", signature(x = "Matrix", i = .cl, j = "missing",
                         drop = "missing"),
          function(x, i, j, ..., drop) {
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
                  stop(.subscript.error.dim)
              }
          })
