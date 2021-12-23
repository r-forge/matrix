## TODO/IDEAS:
## * we would potentially benefit from reusing some of the machinery
##   in R_MAIN_DIR/subscript.c to convert supplied subscripts to
##   "nice" in-bounds integer subscripts, which are necessary to compute
##   "triangular" subscripts ... though none of those utilities appear
##   to be "public" ...
##
## * methods for x[i, j] and x[i, j, drop=] with neither i nor j missing
## * haven't given much thought to long vector support
## * ditto subassignment
## * ditto future '[iz]Matrix' classes
##
## * an efficient 'apply' analogue for "packedMatrix" should be relatively
##   easy to implement now that row/column extraction are fast
##
## NTS:
## * need to use dimnames(x) instead of x@Dimnames to get symmetric
##   dimnames when dealing with "symmetricMatrix"


## UTILITIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NOTATION:
## utility .pM.sub<code>.<class> subsets "packedMatrix"
## code=
## 0: 1-ary indexing, as in x[i]
## 1: 2-ary indexing with one of 'i', 'j' missing, as in x[i, , drop=]
## 2: 2-ary indexing with neither 'i' nor 'j' missing, as in x[i, j, drop=]
## class=
##        index: vector index [dispatches to chr,logi,num]
## chr,logi,num: vector index
##          mat: matrix, lMatrix index
##         null: NULL index

.pM.error.oob <- function() {
    stop("subscript out of bounds")
}
.pM.error.ist <- function(i) {
    if (isS4(i)) {
        stop(sprintf("invalid subscript (S4) class '%s%'", class(i)))
    } else {
        stop(sprintf("invalid subscript type '%s%'", typeof(i)))
    }
}
.pM.error.dim <- function() {
    stop("incorrect number of dimensions")
}

## to avoid defining methods for each "index" subclass separately,
## at the cost of having to dispatch on index type ourselves ...
.pM.sub0.index <- function(x, i) {
    switch(mode(i),
           numeric =
               {
                   if (!is.numeric(i)) {
                       ## use numeric representation of factors, etc.
                       class(i) <- NULL
                   }
                   .pM.sub0.num(x, i)
               },
           logical = .pM.sub0.logi(x, i),
           character = .pM.sub0.chr(x, i),
           .pM.error.ist(i))
}
.pM.sub1.index <- function(x, i, drop, col) {
    switch(mode(i),
           numeric =
               {
                   if (!is.numeric(i)) {
                       ## use numeric representation of factors, etc.
                       class(i) <- NULL
                   }
                   .pM.sub1.num(x, i, drop = drop, col = col)
               },
           logical = .pM.sub1.logi(x, i, drop = drop, col = col),
           character = .pM.sub1.chr(x, i, drop = drop, col = col),
           .pM.error.ist(i))
}

## emulating 'stringSubscript' in R_MAIN_DIR/subscript.c,
## though the first case, 'x[<character but not rank-2 matrix>]',
## is quite pathological and really need not be supported ...
.pM.sub0.chr <- function(x, i) {
    rep.int(x@x[1L][NA], length(i))
}
.pM.sub1.chr <- function(x, i, drop, col) {
    if (length(i) > 0L) {
        nms <- dimnames(x)[[1L + col]]
        if (is.null(nms) || anyNA(i <- match(i, nms))) {
            .pM.error.oob()
        }
        .pM.sub1.num(x, i, drop = drop, col = col)
    } else {
        .pM.sub1.num(x, integer(0L), drop = drop, col = col)
    }
}

.pM.sub0.logi <- function(x, i) {
    ni <- length(i)
    if (ni == 0L) {
        return(x@x[0L])
    }
    n <- x@Dim[1L]
    ## optimize when identical(x@x, as(x, "[dln]geMatrix")@x)
    if (n <= 1L) {
        return(x@x[i])
    }
    ## FIXME: not sure how useful these optimizations are in practice ...
    ## they could be deleted or performed only if 'length(i)' is small
    if (anyNA(i)) {
        ## optimize x[NA], etc.
        if (all(is.na(i))) {
            return(rep.int(x@x[1L][NA], max(n * n, length(i))))
        }
    } else {
        ## optimize x[FALSE], etc.
        if (!any(i)) {
            return(x@x[0L])
        }
        ## optimize x[TRUE], etc.
        if (all(i)) {
            x <- as(x, geClass(x))@x
            if (ni > n * n) {
                length(x) <- ni
            }
            return(x)
        }
    }
    ## dispatch
    ## FIXME: inefficient... is R_MAIN_DIR/subscript.c machinery better?
    ## though notably still an improvement over 'as(x, "matrix")[i]'
    ## which allocates a double vector of length n*n when 'x' is a "dMatrix"
    .pM.sub0.num(x, seq_len(n * n)[i])
}
.pM.sub1.logi <- function(x, i, drop, col) {
    p <- 1L + col # subset dimension [1=i, 2=j]
    d <- x@Dim
    dn <- dimnames(x)
    empty <- function() {
        d[p] <- 0L
        dn[p] <- list(NULL)
        new(geClass(x), x = x@x[0L], Dim = d, Dimnames = dn)
    }
    ni <- length(i)
    if (ni == 0L) {
        return(empty())
    }
    n <- d[1L]
    if (ni > n) {
        stop("logical subscript too long")
    }
    ## FIXME: not sure how useful these optimizations are in practice ...
    ## they could be deleted or performed only if 'length(i)' is small
    if (anyNA(i)) {
        ## optimize x[NA, ], etc.
        if (all(is.na(i))) {
            cl <- geClass(x)
            x <- rep.int(x@x[1L][NA], n * n)
            if (!is.null(dn[[p]])) {
                dn[[p]][] <- NA
            }
            return(new(cl, x = x, Dim = d, Dimnames = dn))
        }
    } else {
        ## optimize x[FALSE, ], etc.
        if (!any(i)) {
            return(empty())
        }
        ## optimize x[TRUE, ], etc.
        if (all(i)) {
            if (n == 1L && !isFALSE(drop[1L])) {
                pp <- 1L + !col # marginal dimension [1=i, 2=j]
                x <- x@x
                if (!is.null(dn[[pp]])) {
                    names(x) <- dn[[pp]]
                }
            }
            return(x)
        }
    }
    ## dispatch
    .pM.sub1.num(x, seq_len(n)[i], drop = drop, col = col)
}

.pM.sub0.num <- function(x, i) {
    if (length(i) == 0L) {
        return(x@x[0L])
    }
    n <- x@Dim[1L]
    ## optimize when identical(x@x, as(x, "[dln]geMatrix")@x)
    if (n <= 1L) {
        return(x@x[i])
    }
    ## FIXME: inefficient... is R_MAIN_DIR/subscript.c machinery better?
    if (any(i < 0, na.rm = TRUE)) {
        i <- seq_len(n * n)[i]
    } else {
        if (is.double(i)) {
            i <- as.integer(i)
        }
        i <- i[i > 0L]
        i[i > n * n] <- NA
    }
    ## Expects "nice" indices: in-bounds integers and NA only
    .Call(packedMatrix_sub0_1ary, x, i)
}
.pM.sub1.num <- function(x, i, drop, col) {
    n <- x@Dim[1L]
    if (any(i >= n + 1L, na.rm = TRUE)) {
        .pM.error.oob()
    }
    ## Expects "nice" indices: in-bounds integers and NA only
    .Call(packedMatrix_sub1, x, seq_len(n)[i], drop, col)
}

## Could be adapted to support indexing with "[dn]Matrix" and "array" ...
## leaving out for now
.pM.sub0.mat <- function(x, i) {
    if (is(i, "lMatrix")) {
        return(.pM.sub0.logi(x, as.vector(i)))
    }
    if (is.logical(i)) {
        dim(i) <- NULL # avoids a copy if !identical(i, as.vector(i))
        return(.pM.sub0.logi(x, i))
    }
    d <- dim(i)
    if (length(d) == 2L && d[2L] == 2L) {
        if (is.numeric(i)) {
            if (is.double(i)) {
                i <- as.integer(i)
                dim(i) <- d
            }
            ## rows containing 0 are deleted, rows containing NA result in NA,
            ## rows containing both are handled according to the first column
            i <- i[i[, 1L] != 0L, , drop = FALSE] # NA,j -> NA,NA
            i <- i[i[, 2L] != 0L, , drop = FALSE]
            if (dim(i)[1L] == 0L) {
                return(x@x[0L])
            }
            if (any(i < 1L, na.rm = TRUE)) {
                stop("negative values are not allowed in a matrix subscript")
            }
            if (any(i > x@Dim[1L], na.rm = TRUE)) {
                .pM.error.oob()
            }
            ## Expects "nice" indices: in-bounds integers and NA only
            .Call(packedMatrix_sub0_2ary, x, i)
        } else if (is.character(i)) {
            if (d[1L] == 0L) {
                return(x@x[0L])
            }
            dn <- dimnames(x)
            m <- c(match(i[, 1L], dn[[1L]]), match(i[, 2L], dn[[2L]]))
            dim(m) <- d
            ## error if character row contains zero NA but integer row
            ## contains at least one NA, indicating nonmatch that should
            ## not be ignored
            if (any(rowSums(is.na(i)) == 0L & rowSums(is.na(m)) > 0L)) {
                .pM.error.oob()
            }
            ## Expects "nice" indices: in-bounds integers and NA only
            .Call(packedMatrix_sub0_2ary, x, m)
        } else {
            .pM.error.ist(i)
        }
    } else {
        dim(i) <- NULL
        if (is.numeric(i)) {
            .pM.sub0.num(x, i)
        } else if (is.character(i)) {
            .pM.sub0.chr(x, i)
        } else {
            .pM.error.ist(i)
        }
    }
}
.pM.sub1.mat <- function(x, i, drop, col) {
    dim(i) <- NULL # avoids a copy if !identical(i, as.vector(i))
    if (is.numeric(i)) {
        .pM.sub1.num(x, i, drop = drop, col = col)
    } else if (is.logical(i)) {
        .pM.sub1.logi(x, i, drop = drop, col = col)
    } else if (is.character(i)) {
        .pM.sub1.chr(x, i, drop = drop, col = col)
    } else {
        .pM.error.ist(i)
    }
}

## unused argument 'i' included for consistency with other utilities
.pM.sub0.null <- function(x, i) {
    .pM.sub0.num(x, integer(0L))
}
.pM.sub1.null <- function(x, i, drop, col) {
    .pM.sub1.num(x, integer(0L), drop = drop, col = col)
}


## METHOD DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("[", signature(x = "packedMatrix", i = "missing", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("pM[m, m, m] : nargs() = %d", na), .M.level = 2)
	      if (na == 2L) { # x[]
                  x
              } else if (na == 3L) { # x[, ]
                  ## drop=TRUE implicit
                  if (x@Dim[1L] == 1L) x@x else x
              } else { # x[, , ], and so on
                  .pM.error.dim()
              }
          })

setMethod("[", signature(x = "packedMatrix", i = "missing", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("pM[m, m, %s] : nargs() = %d", drop, na), .M.level = 2)
              if (na == 4L) { # x[, , drop=], x[, drop=, ], x[drop=, , ]
                  if (!isFALSE(drop[1L]) && x@Dim[1L] == 1L) x@x else x
              } else if (na < 4L) { # x[drop=], x[, drop=], x[drop=, ]
                  x
              } else { # x[, , , drop=], and so on
                  .pM.error.dim()
              }
          })

## Generalizing here to avoid repetition ...
.j.map <- c(index = "index", matrix = "mat", `NULL` = "null")
.i.map <- c(.j.map, lMatrix = "mat")
.op <- options(keep.source = FALSE)
for (.i.cl in names(.i.map)) {
    .nm0 <- as.name(paste0(".pM.sub0.", .i.map[[.i.cl]]))
    .nm1 <- as.name(paste0(".pM.sub1.", .i.map[[.i.cl]]))
    eval(bquote({
        setMethod("[", signature(x = "packedMatrix", i = .i.cl, j = "missing", drop = "missing"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[%s, m, m] : nargs() = %d", .(.i.cl), na), .M.level = 2)
                      if (na == 2L) { # x[i]
                          .(.nm0)(x, i)
                      } else if (na == 3L) { # x[i, ]
                          .(.nm1)(x, i, drop = TRUE, col = FALSE)
                      } else { # x[i, , ], etc.
                          .pM.error.dim()
                      }
                  })
        setMethod("[", signature(x = "packedMatrix", i = .i.cl, j = "missing", drop = "logical"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[%s, m, %s] : nargs() = %d", .(.i.cl), drop, na), .M.level = 2)
                      if (na == 3L) { # x[i, drop=]
                          .(.nm0)(x, i)
                      } else if (na == 4L) { # x[i, , drop=]
                          .(.nm1)(x, i, drop = drop, col = FALSE)
                      } else { # x[i, , , drop=], etc.
                          .pM.error.dim()
                      }
                  })
    })) # eval(bquote({
}
for (.j.cl in names(.j.map)) {
    .nm0 <- as.name(paste0(".pM.sub0.", .j.map[[.j.cl]]))
    .nm1 <- as.name(paste0(".pM.sub1.", .j.map[[.j.cl]]))
    eval(bquote({
        setMethod("[", signature(x = "packedMatrix", i = "missing", j = .j.cl, drop = "missing"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[m, %s, m] : nargs() = %d", .(.j.cl), na), .M.level = 2)
                      if (na == 2L) { # x[j=]
                          .(.nm0)(x, j)
                      } else if (na == 3L) { # x[, j]
                          .(.nm1)(x, j, drop = TRUE, col = TRUE)
                      } else { # x[, j, ], etc.
                          .pM.error.dim()
                      }
                  })
        setMethod("[", signature(x = "packedMatrix", i = "missing", j = .j.cl, drop = "logical"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[m, %s, %s] : nargs() = %d", .(.j.cl), drop, na), .M.level = 2)
                      if (na == 3L) { # x[j=, drop=]
                          .(.nm0)(x, j)
                      } else if (na == 4L) { # x[, j, drop=]
                          .(.nm1)(x, j, drop = drop, col = TRUE)
                      } else { # x[, j, , drop=], etc.
                          .pM.error.dim()
                      }
                  })
    })) # eval(bquote({
}
options(.op)
rm(.i.map, .j.map, .i.cl, .j.cl, .op, .nm0, .nm1)

setMethod("t", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_t, x))
setMethod("diag", signature(x = "packedMatrix"),
          function(x, nrow, ncol, names) .Call(packedMatrix_diag_get, x, names))
setMethod("diag<-", signature(x = "packedMatrix"),
          function(x, value) .Call(packedMatrix_diag_set, x, value))
