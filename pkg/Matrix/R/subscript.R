## METHODS FOR GENERIC: [
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
##             *_subscript_1ary_vec(x, i         )
##             *_subscript_1ary_mat(x, i         )
##             *_subscript_2ary    (x, i, j, drop)
##
##       for * = CRsparse,Tsparse,unpackedMatrix,packedMatrix
##
##       diagonalMatrix and indMatrix should go via CsparseMatrix,
##       except for row indexing of indMatrix

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

.subscript.1ary.vec <- function(x, i, .NAME) { # .NAME = *_subscript_1ary_vec
    switch(typeof(i),
           double =
               {
                   m <- min(1, i, na.rm = TRUE)
                   if(m < 1)
                       i <- if(m <= -1)
                                seq_len(prod(x@Dim))[i]
                            else {
                                if(is.object(i))
                                    i <- as.double(i)
                                i[i >= 1]
                            }
                   .Call(.NAME, x, i)
               },
           integer =
               {
                   m <- min(1L, i, na.rm = TRUE)
                   if(m < 1L)
                       i <- if(m <= -1L)
                                seq_len(prod(x@Dim))[i]
                            else {
                                if(is.object(i))
                                    i <- as.integer(i)
                                i[i >= 1L]
                            }
                   .Call(.NAME, x, i)
               },
           logical =
               {
                   i <- seq_len(prod(x@Dim))[i]
                   .Call(.NAME, x, i)
               },
           character =
               {
                   rep.int(if(.hasSlot(x, "x")) x@x[NA_integer_] else NA,
                           length(i))
               },
           stop(.subscript.error.ist(i), domain = NA))
}

.subscript.1ary.mat <- function(x, i, .NAME) { # .NAME = *_subscript_1ary_mat
    if(isS4(i)) {
        cld <- getClassDef(class(i))
        if(!extends(cld,  "Matrix"))
            stop(.subscript.error.ist(i), domain = NA)
        if(!extends(cld, "dMatrix") || i@Dim[2L] != 2L)
            return(.subscript.1ary.vec(x, as.vector(i))) # FIXME: sparseVector
        i <- as(i, "matrix")
    } else if(is.logical(i) || length(di <- dim(i)) != 2L || di[2L] != 2L)
        return(.subscript.1ary.vec(x, i))
    switch(typeof(i),
           double =,
           integer =
               {
                   ## * rows containing 0 are deleted
                   ## * rows containing NA result in NA
                   ## * rows containing both 0 and NA are handled
                   ##   according to value in first column
                   storage.mode(i) <- "integer"
                   i <- i[i[, 1L] != 0L, , drop = FALSE] # [NA,j] -> [NA,NA]
                   i <- i[i[, 2L] != 0L, , drop = FALSE]
                   if(min(1L, i, na.rm = TRUE) < 1L)
                       stop(.subscript.error.neg)
                   d <- x@Dim
                   m <- d[1L]
                   n <- d[2L]
                   if(m == n) {
                       if(max(n, i, na.rm = TRUE) > n)
                           stop(.subscript.error.oob)
                   } else {
                       if(max(m, i[, 1L], na.rm = TRUE) > m ||
                          max(n, i[, 2L], na.rm = TRUE) > n)
                           stop(.subscript.error.oob)
                   }
                   .Call(.NAME, x, i)
               },
           character =
               {
                   dn <- dimnames(x)
                   m <- c(match(i[, 1L], dn[[1L]]), match(i[, 2L], dn[[2L]]))
                   dim(m) <- di
                   if(any(rowSums(is.na(i)) == 0L & rowSums(is.na(m)) > 0L))
                       ## error if character row contains zero NA and
                       ## integer row contains at least one NA,
                       ## indicating non-match that cannot be ignored
                       stop(.subscript.error.oob)
                   .Call(.NAME, x, m)
               },
           stop(.subscript.error.ist(i), domain = NA))
}

.subscript.2ary <- function(x, i, j, drop, .NAME) { # .NAME = *_suscriptb_2ary
    d <- x@Dim
    index <- list(if(missing(i)) NULL else if(is.null(i)) integer(0L) else i,
                  if(missing(j)) NULL else if(is.null(j)) integer(0L) else j)
    for(pos in 1:2) {
        if(!is.null(k <- index[[pos]])) {
            index[[pos]] <-
                switch(typeof(k),
                       double =
                           {
                               r <- d[pos]
                               if(is.object(k))
                                   k <- as.double(k)
                               if(max(r + 1, k, na.rm = TRUE) > r + 1)
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
                               seq_len(r)[k]
                           },
                       logical =
                           {
                               r <- d[pos]
                               if(length(k) > r)
                                   stop(.subscript.error.lng)
                               seq_len(r)[k]
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
    .Call(.NAME, x, index[[1L]], index[[2L]], if(missing(drop)) TRUE else drop)
}
