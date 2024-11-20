## METHODS FOR GENERIC: [
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.subscript.invalid <-
function(i) {
    if (is.object(i))
        gettextf("invalid subscript class \"%s\"", class(i)[1L])
    else gettextf("invalid subscript type \"%s\"", typeof(i))
}

.subscript.recycle <-
function(i, n, pattern) {
    ## Return what would be the result of seq_len(n)[as.vector(i)] for
    ## 'i' of class nsparseVector (pattern = TRUE) or lsparseVector;
    ## length(i) < n = 0,...,(2^31-1)^2.
    i.length <- length(i)
    if (i.length >= n)
        stop("should never happen ...")
    i.i <- i@i
    if (length(i.i) == 0L)
        return(integer(0L))
    r <- ceiling(n / i.length)
    as. <- if (r * i.length - 1 < .Machine[["integer.max"]])
               as.integer
           else as.double
    J <- as.(i.i) - 1L # go to 0-based
    J <- rep.int(J, r) +
        rep(as.(seq.int(from = 0L, by = i.length, length.out = r)),
            each = length(J))
    J <- if (pattern)
             (if (J[length(J)] >= n) J[      J < n] else J     )
         else
             (if (J[length(J)] >= n) J[i@x & J < n] else J[i@x])
    if (J[length(J)] >= 0x1p+53) {
        ## Elements of 'J' may be the result of rounding a
        ## non-representable number; no obvious way (in R)
        ## to handle this possibility, hence:
        warning(gettextf("subscripts exceeding %s replaced with NA", "2^53"),
                domain = NA)
        J[J >= 0x1p+53] <- NA
    }
    J + 1L # return to 1-based
}

## x[i]
.subscript.1ary <-
function(x, i) {
    x.length <- length(x)
    if (is.null(i))
        i <- integer(0L)
    else if (typeof(i) == "S4") {
        .isM <- .isMatrix(i)
        .isV <- !.isM && .isVector(i)
        if ((!.isM && !.isV) || (kind <- .M.kind(i)) == "z")
            stop(.subscript.invalid(i), domain = NA)
        if (.isM)
            i <-
                if ((kind == "i" || kind == "d") && i@Dim[2L] == 2L)
                    .M2m(i)
                else if (.isDense(i))
                    .M2v(i)
                else {
                    .isV <- TRUE
                    .M2V(i)
                }
        if (.isV)
            i <-
                if (kind == "i" || kind == "d")
                    i@x
                else if (length(i) < x.length)
                    .subscript.recycle(i, x.length, kind == "n")
                else if (kind == "n")
                    i@i
                else i@i[i@x]
    }
    else if (is.object(i))
        i <- switch(typeof(i),
                    double =,
                    integer =,
                    character =
                        if (length(d <- dim(i)) == 2L && d[2L] == 2L)
                            `dim<-`(as.vector(i), d)
                        else as.vector(i),
                    logical = as.vector(i),
                    i)
    switch(typeof(i),
           double =,
           integer =
               {
                   if (length(d <- dim(i)) == 2L && d[2L] == 2L) {

                   if (min(0L, i, na.rm = TRUE) <= -1L)
                       stop("negative values are not allowed in a matrix subscript")
                   d <- x@Dim
                   m <- d[1L]
                   n <- d[2L]
                   i. <- i[, 1L]
                   j. <- i[, 2L]
                   if (max(m, i., na.rm = TRUE) - 1L >= m ||
                       max(n, j., na.rm = TRUE) - 1L >= n)
                       stop("subscript out of bounds")
                   if (is.double(i))
                       storage.mode(i.) <- storage.mode(j.) <-
                           storage.mode(i) <- "integer"
                   ## rows containing 0 are deleted;
                   ## rows containing NA are kept;
                   ## rows containing both 0 and NA are handled
                   ## according to value in first column
                   if (is.na(a <- all(i.)) || !a)
                       i <- i[i. >= 1L, , drop = FALSE]
                   if (!all(j., na.rm = TRUE))
                       i <- i[j. >= 1L, , drop = FALSE]
                   ..subscript.1ary.2col(x, i)

                   } else {

                   r <- min(1L, i, na.rm = TRUE)
                   if (r < 1L)
                       i <- if (r <= -1L)
                                seq_len(x.length)[i] # FIXME
                            else i[i >= 1L]
                   ..subscript.1ary(x, i)

                   }
               },
           logical =
               {
                   i.length <- length(i)
                   if (i.length > 0L && !is.na(a <- all(i)) && a) {
                       if (i.length <= x.length)
                           as.vector(x)
                       else c(as.vector(x), rep.int(NA, i.length - x.length))
                   }
                   else .subscript.1ary(x, .m2V(i)) # recursively
               },
           character =
               {
                   if (length(d <- dim(i)) == 2L && d[2L] == 2L) {

                   dn <- dimnames(x)
                   m <- c(match(i[, 1L], dn[[1L]]),
                          match(i[, 2L], dn[[2L]]))
                   dim(m) <- d
                   if (any(!rowSums(is.na(i)) & rowSums(is.na(m))))
                       ## error if character row contains zero NA
                       ## and integer row contains at least one NA,
                       ## indicating non-match that should not be
                       ## ignored
                       stop("subscript out of bounds")
                   ..subscript.1ary.2col(x, m)

                   } else {

                   rep.int(if (.M.kind(x) == "n") NA else x@x[NA_integer_],
                           length(i))

                   }
               },
           stop(.subscript.invalid(i), domain = NA))
}

## x[i] where 'i' is a vector of type "integer" or "double"
## with elements greater than or equal to 1 (or NA)
..subscript.1ary <-
function(x, i,
         shape = .M.shape(x),
         repr = .M.repr(x),
         unsorted = is.unsorted(i)) {
    if (repr == "R" || ((repr == "C" || repr == "T") && shape == "s")) {
        x.length <- length(x)
        if (x.length <= 0x1p+53) {
            r <- max(x.length, i, na.rm = TRUE)
            if (r - 1L >= x.length)
                i[i - 1L >= x.length] <- NA
        }
        else if (is.double(i)) {
            r <- max(0x1p+53, i, na.rm = TRUE)
            if (r > 0x1p+53) {
                if (any(i > 0x1p+53 & i <= x.length, na.rm = TRUE))
                    ## could be avoided in C, which has 64-bit integers :
                    warning(gettextf("subscripts exceeding %s replaced with NA", "2^53"),
                            domain = NA)
                i[i > 0x1p+53] <- NA
            }
        }
        i1s <- i - 1L
        m <- x@Dim[1L]
        i. <- as.integer(i1s %% m)
        j. <- as.integer(i1s %/% m)
        if (shape == "s") {
            op <- if (x@uplo == "U") `>` else `<`
            if (length(w <- which(op(i., j.)))) {
                i.. <- i.[w]
                j.. <- j.[w]
                i.[w] <- j..
                j.[w] <- i..
                unsorted <- TRUE # Bug #6839
            }
        }
    }
    o <-
    if (repr == "R")
        order(i., j.)
    else if ((repr == "C" || repr == "T") &&
             (is.na(unsorted) || unsorted))
        (if (shape == "s") order(j., i.) else sort.list(i))
    .Call(R_subscript_1ary, x, i, o)
}

## x[i] where 'i' is a 2-column matrix of type "integer"
## with i[, 1L] in 1:m (or NA) and i[, 2L] in 1:n (or NA)
..subscript.1ary.2col <-
function(x, i,
         shape = .M.shape(x),
         repr = .M.repr(x)) {
    o <-
    if (repr == "C" || repr == "R" || repr == "T") {
        i. <- i[, 1L]
        j. <- i[, 2L]
        if (shape == "s") {
            op <- if (x@uplo == "U") `>` else `<`
            if (length(w <- which(op(i., j.)))) {
                i. <- i[, 1L]
                j. <- i[, 2L]
            }
        }
        if (repr == "R") {
            if (anyNA(j.))
                i.[is.na(j.)] <- NA
            order(i., j.)
        } else {
            if (anyNA(i.))
                j.[is.na(i.)] <- NA
            order(j., i.)
        }
    }
    .Call(R_subscript_1ary_2col, x, i, o)
}

## x[i, j, drop]
.subscript.2ary <-
function(x, i, j, drop) {
    l <-
    list(if (missing(i)) NULL else if (is.null(i)) integer(0L) else i,
         if (missing(j)) NULL else if (is.null(j)) integer(0L) else j)
    d <- x@Dim
    for (pos in 1:2)
    if (!is.null(k <- l[[pos]])) {
        r <- d[pos]
        if (typeof(k) == "S4") {
            .isM <- .isMatrix(k)
            .isV <- !.isM && .isVector(k)
            if ((!.isM && !.isV) || (kind <- .M.kind(k)) == "z")
                stop(.subscript.invalid(k), domain = NA)
            if (.isM)
                k <-
                    if (.isDense(k))
                        .M2v(k)
                    else {
                        .isV <- TRUE
                        .M2V(k)
                    }
            if (.isV)
                k <-
                    if (kind == "i" || kind == "d")
                        k@x
                    else if (length(k) > r)
                        stop("logical subscript too long")
                    else if (length(k) < r)
                        .subscript.recycle(k, r, kind == "n")
                    else if (kind == "n")
                        k@i
                    else k@i[k@x]
        }
        else if (is.object(k))
            k <- switch(typeof(i),
                        double =,
                        integer =,
                        logical =,
                        character = as.vector(k),
                        k)
        l[pos] <-
        list(switch(typeof(k),
                    double =,
                    integer =
                        {
                            if (max(r, k, na.rm = TRUE) - 1L >= r)
                                stop("subscript out of bounds")
                            if (min(1L, k, na.rm = TRUE) < 1L)
                                seq_len(r)[k]
                            else as.integer(k)
                        },
                    logical =
                        {
                            if (length(k) > r)
                                stop("logical subscript too long")
                            if (length(k) > 0L && !is.na(a <- all(k)) && a)
                                NULL
                            else seq_len(r)[k]
                        },
                    character =
                        {
                            if (is.null(nms <- dimnames(x)[[pos]]) ||
                                anyNA(k <- match(k, nms)))
                                stop("subscript out of bounds")
                            else k
                        },
                    stop(.subscript.invalid(k), domain = NA)))
    }
    if (is.double(lengths(l, use.names = FALSE)))
        stop(gettextf("dimensions cannot exceed %s", "2^31-1"),
             domain = NA)
    ..subscript.2ary(x, l[[1L]], l[[2L]], drop = drop[1L])
}

## x[i, j, drop] where 'i' and 'j' are vectors of type "integer"
## of length not exceeding 2^31-1 with 'i' in 1:m (or NA) and 'j'
## in 1:n (or NA) ... NULL => missing
..subscript.2ary <-
function(x, i, j, drop) {
    if (is.null(i) && is.null(j))
        r <- x
    else {
        r <- .Call(R_subscript_2ary, x, i, j)
        dn <- dimnames(x)
        if (!(is.null(i) || is.null(rn <- dn[[1L]])))
            dn[1L] <- list(if (length(i)) rn[i] else NULL)
        if (!(is.null(j) || is.null(cn <- dn[[2L]])))
            dn[2L] <- list(if (length(j)) cn[j] else NULL)
        r@Dimnames <- dn
    }
    if ((is.na(drop) || drop) && any(r@Dim == 1L))
        drop(.M2m(r))
    else r
}

setMethod("[",
          c(x = "Matrix", i = "missing", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 2L)
                  ## x[]
                  x
              else if (na == 3L)
                  ## x[, ]
                  drop(x)
              else
                  ## x[, , ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "missing", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na < 4L)
                  ## x[drop=], x[, drop=], x[drop=, ]
                  x
              else if (na == 4L)
                  ## x[, , drop=], x[, drop=, ], x[drop=, , ]
                  if (is.na(drop <- drop[1L]) || drop) drop(x) else x
              else
                  ## x[, , , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "ANY", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 2L)
                  ## x[i=]
                  .subscript.1ary(x, i)
              else if (na == 3L)
                  ## x[i=, ], x[, i=]
                  .subscript.2ary(x, i, , drop = TRUE)
              else
                  ## x[i=, , ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "ANY", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 3L)
                  ## x[i=, drop=]
                  .subscript.1ary(x, i)
              else if (na == 4L)
                  ## x[i=, , drop=], x[, i=, drop=]
                  .subscript.2ary(x, i, , drop = drop)
              else
                  ## x[i=, , , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "missing", j = "ANY", drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 2L)
                  ## x[j=]
                  .subscript.1ary(x, j)
              else if (na == 3L)
                  ## x[j=, ], x[, j=]
                  .subscript.2ary(x, , j, drop = TRUE)
              else
                  ## x[, j=, ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "missing", j = "ANY", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 3L)
                  ## x[j=, drop=]
                  .subscript.1ary(x, j)
              else if (na == 4L)
                  ## x[j=, , drop=], x[, j=, drop=]
                  .subscript.2ary(x, , j, drop = drop)
              else
                  ## x[, j=, , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "ANY", j = "ANY", drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 3L)
                  ## x[i=, j=], x[j=, i=]
                  .subscript.2ary(x, i, j, drop = TRUE)
              else
                  ## x[i=, j=, ], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "Matrix", i = "ANY", j = "ANY", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              na <- nargs()
              if (na == 4L)
                  ## x[i=, j=, drop=], x[j=, i=, drop=]
                  .subscript.2ary(x, i, j, drop = drop)
              else
                  ## x[i=, j=, , drop=], etc.
                  stop("incorrect number of dimensions")
          })

setMethod("[",
          c(x = "sparseVector", i = "missing", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if (nargs() != 2L)
                  stop("incorrect number of dimensions")
              x
          })

setMethod("[",
          c(x = "sparseVector", i = "ANY", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop = TRUE) {
              if (nargs() != 2L)
                  stop("incorrect number of dimensions")
              x.length <- length(x)
              if (is.null(i))
                  i <- integer(0L)
              else if (typeof(i) == "S4") {
                  .isM <- .isMatrix(i)
                  .isV <- !.isM && .isVector(i)
                  if ((!.isM && !.isV) || (kind <- .M.kind(i)) == "z")
                      stop(.subscript.invalid(i), domain = NA)
                  if (.isM)
                      i <-
                          if (.isDense(i))
                              .M2v(i)
                          else {
                              .isV <- TRUE
                              .M2V(i)
                          }
                  if (.isV)
                      i <-
                          if (kind == "i" || kind == "d")
                              i@x
                          else if (length(i) < x.length)
                              .subscript.recycle(i, x.length, kind == "n")
                          else if (kind == "n")
                              i@i
                          else i@i[i@x]
              }
              else if (is.object(i))
                  i <- switch(typeof(i),
                              double =,
                              integer =,
                              logical =,
                              character = as.vector(i),
                              i)
              kind <- .M.kind(x)
              switch(typeof(i),
                     double =,
                     integer =
                         {
                             trunc. <-
                             function(x)
                                 if (is.double(x)) trunc(x) else x
                             x.i <- trunc.(x@i)
                             r <- min(1L, i, na.rm = TRUE)
                             if (r <= -1L) {
                                 if (r <= -x.length - 1L)
                                     i <- i[i > -x.length - 1L]
                                 r <- max(-1L, i)
                                 if (is.na(r) || r >= 1L)
                                     stop("only zeros may be mixed with negative subscripts")
                                 if (r > -1L)
                                     i <- i[i <= -1L]
                                 d <- unique.default(sort.int(-trunc.(i)))
                                 m <- match(x.i, d, 0L)
                                 y <- new(paste0(kind, "sparseVector"))
                                 y@length <- length(x) - length(d)
                                 y@i <-
                                     {
                                         tmp <- x.i[m == 0L]
                                         tmp - findInterval(tmp, d) # !!
                                     }
                                 if (kind != "n")
                                     y@x <- x@x[m == 0L]
                                 y
                             } else {
                                 if (r < 1L)
                                     i <- i[i >= 1L]
                                 if (max(0L, i, na.rm = TRUE) - 1L >= x.length)
                                     i[i - 1L >= x.length] <- NA
                                 a <- anyNA(i)
                                 if (a && kind == "n") {
                                     x <- .V2kind(x, "l")
                                     kind <- "l"
                                 }
                                 d <- trunc.(i)
                                 m <- match(d, x.i, 0L)
                                 y <- new(paste0(kind, "sparseVector"))
                                 y@length <- length(i)
                                 y@i <-
                                     if (!a)
                                         which(m != 0L)
                                     else {
                                         i. <- is.na(i)
                                         m[i.] <- NA
                                         which(m != 0L | i.)
                                     }
                                 if (kind != "n")
                                     y@x <- x@x[m]
                                 y
                             }
                         },
                     logical =
                         {
                             i.length <- length(i)
                             if (i.length > 0L && !is.na(a <- all(i)) && a) {
                                 if (i.length <= x.length)
                                     x
                                 else {
                                     if (kind == "n") {
                                         x <- .V2kind(x, "l")
                                         kind <- "l"
                                     }
                                     y <- new(paste0(kind, "sparseVector"))
                                     y@length <- i.length
                                     y@i <- c(x@i, (x.length + 1):i.length)
                                     if (kind != "n")
                                         y@x <- c(x@x, rep.int(NA, i.length - x.length))
                                     y
                                 }
                             }
                             else x[.m2V(i)] # recursively
                         },
                     character =
                         {
                             i.length <- length(i)
                             y <- new(paste0(if (kind == "n") "l" else kind, "sparseVector"))
                             y@length <- i.length
                             y@i <- seq_len(i.length)
                             y@x <- rep.int(x@x[NA], i.length)
                             y
                         },
                     stop(.subscript.invalid(i), domain = NA))
          })


## METHODS FOR GENERIC: head
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("head", c(x = "Matrix"),
          head.matrix)

setMethod("head", c(x = "sparseVector"),
          function(x, n = 6L, ...) {
              stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
              len <- length(x)
              n <- if (n < 0L) max(len + n, 0L) else min(n, len)
              if (n >= len)
                  return(x)
              nnz <- length(i <- x@i)
              x@length <- n <- if (is.integer(i)) as.integer(n) else trunc(n)
              if (nnz > 0L && i[nnz] > n) {
                  pattern <- .M.kind(x) == "n"
                  if (i[1L] > n) {
                      x@i <- integer(0L)
                      if (!pattern)
                          x@x <- x@x[0L]
                  } else {
                      ii <- 1L:(which.max(i > n) - 1L)
                      x@i <- i[ii]
                      if (!pattern)
                          x@x <- x@x[ii]
                  }
              }
              x
          })


## METHODS FOR GENERIC: tail
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("tail", c(x = "Matrix"),
          tail.matrix)

setMethod("tail", c(x = "sparseVector"),
          function(x, n = 6L, ...) {
              stopifnot(is.numeric(n), length(n) == 1L, !is.na(n))
              len <- length(x)
              n <- if (n < 0L) max(len + n, 0L) else min(n, len)
              if (n >= len)
                  return(x)
              nnz <- length(i <- x@i)
              x@length <- n <- if (is.integer(i)) as.integer(n) else trunc(n)
              if (nnz > 0L && i[1L] <= (k <- len - n)) {
                  pattern <- .M.kind(x) == "n"
                  if (i[nnz] <= k) {
                      x@i <- integer(0L)
                      if (!pattern)
                          x@x <- x@x[0L]
                  } else {
                      ii <- which.min(i <= k):nnz
                      x@i <- i[ii] - k
                      if (!pattern)
                          x@x <- x@x[ii]
                  }
              }
              x
          })
