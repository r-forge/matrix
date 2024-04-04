## METHODS FOR GENERIC: BunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("BunchKaufman", c(x = "dsyMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(syMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", c(x = "dspMatrix"),
          function(x, warnSing = TRUE, ...)
              .Call(spMatrix_trf, x, as.logical(warnSing)))

setMethod("BunchKaufman", c(x = "matrix"),
          function(x, uplo = "U", ...)
              BunchKaufman(.m2dense(x, ",sy", uplo), ...))


## METHODS FOR CLASS: denseBunchKaufman
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setAs("denseBunchKaufman", "dtrMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if(!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if(!packed) to else unpack(to)
      })

setAs("denseBunchKaufman", "dtpMatrix",
      function(from) {
          packed <- length(from@x) != prod(from@Dim)
          to <- new(if(!packed) "dtrMatrix" else "dtpMatrix")
          to@Dim <- from@Dim
          to@uplo <- from@uplo
          to@x <- from@x
          if(!packed) pack(to) else to
      })

.f1 <- function(x, which, ...) {
    r <- .Call(denseBunchKaufman_expand, x)
    b <- length(r) - 1L
    switch(which,
           "DU" =, "DL" =
               {
                   if(!endsWith(which, x@uplo))
                       stop(gettextf("%s=\"%s\" invalid for %s@uplo=\"%s\"",
                                     "which", which, "x", x@uplo),
                            domain = NA)
                   r[[b + 1L]]
               },
           "U" =, "U." =, "L" =, "L." =
               {
                   if(!startsWith(which, x@uplo))
                       stop(gettextf("%s=\"%s\" invalid for %s@uplo=\"%s\"",
                                     "which", which, "x", x@uplo),
                            domain = NA)
                   if(b > 0L) {
                       m <- r[[b]]
                       if(b > 1L)
                           for(i in (b - 1L):1L)
                               m <- r[[i]] %*% m
                       if(endsWith(which, ".")) t(m) else m
                   } else {
                       m <- new("ddiMatrix")
                       m@Dim <- x@Dim
                       m@diag <- "U"
                       m
                   }
               },
           stop(gettextf("'%s' is not \"%1$s\", \"D%1$s\", or \"%1$s.\"",
                         "which", x@uplo),
                domain = NA))
}

## returning
## list(U, DU, U') where A = U DU U' and U = P[b] U[b] ... P[1] U[1]
## OR
## list(L, DL, L') where A = L DL L' and L = P[1] L[1] ... P[b] L[b]
.f2 <- function(x, complete = FALSE, ...) {
    r <- .Call(denseBunchKaufman_expand, x)
    b <- length(r) - 1L
    if(complete) {
        if(b > 0L)
            r <- c(r, lapply(r[b:1L], t))
    } else {
        if(b > 0L) {
            m <- r[[b]]
            if(b > 1L)
                for(i in (b - 1L):1L)
                    m <- r[[i]] %*% m
            r <- list(m, r[[b + 1L]], t(m))
        } else {
            m <- new("ddiMatrix")
            m@Dim <- x@Dim
            m@diag <- "U"
            r <- list(m, r[[1L]], m)
        }
        names(r) <- if(x@uplo == "U") c("U", "DU", "U.") else c("L", "DL", "L.")
    }
    dn <- x@Dimnames
    if(length(r) == 1L)
        r[[1L]]@Dimnames <- dn
    else {
        r[[1L]]@Dimnames <- c(dn[1L], list(NULL))
        r[[length(r)]]@Dimnames <- c(list(NULL), dn[2L])
    }
    r
}

setMethod("expand1", c(x = "denseBunchKaufman"), .f1)
setMethod("expand2", c(x = "denseBunchKaufman"), .f2)

rm(.f1, .f2)
