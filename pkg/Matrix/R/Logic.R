## METHODS FOR GENERIC: Logic (group)                   WORK IN PROGRESS
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Logic") # excluding unary "!" -> ./not.R
## [1] "&" "|"

.Logic.swap <- function(e1, e2)
    switch(.Generic, "&" = e2 & e1, "|" = e2 | e1,
           stop(gettextf("unexpected .Generic=\"%s\" in '%s' method",
                         .Generic, "Logic"),
                domain = NA))


## .... denseMatrix ....................................................

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "denseMatrix", e2 = "vector"),
          function(e1, e2) {

          })


## .... sparseMatrix ...................................................

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "denseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "sparseVector"),
          function(e1, e2) {

          })

setMethod("Logic", signature(e1 = "sparseMatrix", e2 = "vector"),
          function(e1, e2) {

          })


## .... sparseVector ...................................................

setMethod("Logic", signature(e1 = "sparseVector", e2 = "Matrix"),
          .Logic.swap)

setMethod("Logic", signature(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) {
              if(e1.longer <- (n <- length(e1)) > (N <- length(e2))) {
                  tmp <- n
                  n <- N
                  N <- tmp
              }
              if(n > 0L && n < N) {
                  if(N %% n != 0L)
                      warning("longer object length is not a multiple of shorter object length")
                  if(e1.longer)
                      e2 <- rep(e2, length.out = N)
                  else e1 <- rep(e1, length.out = N)
              }
              k1 <- .V.kind(e1)
              k2 <- .V.kind(e2)
              pattern <- k1 == "n" && k2 == "n"
              r <- new(if(pattern) "nsparseVector" else "lsparseVector")
              r@length <- if(n == 0L) n else N
              if(.Generic == "&") {
                  if(n > 0L &&
                     length(i1 <- e1@i) > 0L &&
                     length(i2 <- e2@i) > 0L &&
                     length(i3 <- if(k1 == "n")
                                      i2[m2 <- match(i2, i1, 0L) > 0L]
                                  else i1[m1 <- match(i1, i2, 0L) > 0L]) > 0L) {
                      r@i <- i3
                      if(!pattern)
                          r@x <- if(k1 == "n")
                                     as.logical(e2@x[m2])
                                 else if(k2 == "n")
                                     as.logical(e1@x[m1])
                                 else e1@x[m1] & e2@x[match(i2, i1, 0L) > 0L]
                  }
              } else if(.Generic == "|") {
                  if(n > 0L &&
                     (length(i1 <- e1@i) > 0L | length(i2 <- e2@i) > 0L)) {
                      r@i <- i3 <- sort.int(unique.default(c(i1, i2)))
                      if(!pattern) {
                          m1 <- match(i3, i1, 0L) > 0L
                          m2 <- match(i3, i2, 0L) > 0L
                          x3 <- logical(length(i3))
                          if(k1 == "n") {
                              x3[m1] <- TRUE
                              x3[m2] <- x3[m2] | e2@x
                          } else if(k2 == "n") {
                              x3[m2] <- TRUE
                              x3[m1] <- x3[m1] | e1@x
                          } else {
                              x3[m1] <- as.logical(e1@x)
                              x3[m2] <- x3[m2] | e2@x
                          }
                          r@x <- x3
                      }
                  }
              } else stop(gettextf("unexpected .Generic=\"%s\" in '%s' method",
                                   .Generic, "Logic"),
                          domain = NA)
              r
          })

setMethod("Logic", signature(e1 = "sparseVector", e2 = "vector"),
          function(e1, e2) {
              if(e1.longer <- (n <- length(e1)) > (N <- length(e2))) {
                  tmp <- n
                  n <- N
                  N <- tmp
              }
              if(n > 0L && n < N) {
                  if(N %% n != 0L)
                      warning("longer object length is not a multiple of shorter object length")
                  if(e1.longer)
                      e2 <- rep(e2, length.out = N)
                  else e1 <- rep(e1, length.out = N)
              }
              k1 <- .V.kind(e1)
              if(.Generic == "&") {
                  r <- new("lsparseVector")
                  r@length <- if(n == 0L) n else N
                  if(n > 0L && length(i1 <- e1@i) > 0L) {
                      r@i <- i1
                      r@x <- if(k1 == "n") as.logical(e2[i1]) else e1@x & e2[i1]
                  }
              } else if(.Generic == "|") {
                  r <- if(n == 0L) logical(0L) else as.logical(e2)
                  if(n > 0L && length(i1 <- e1@i))
                      r[i1] <- if(k1 == "n") TRUE else e1@x | r[i1]
              } else stop(gettextf("unexpected .Generic=\"%s\" in '%s' method",
                                   .Generic, "Logic"),
                          domain = NA)
              r
          })


## .... vector .........................................................

setMethod("Logic", signature(e1 = "vector", e2 = "Matrix"),
          .Logic.swap)

setMethod("Logic", signature(e1 = "vector", e2 = "sparseVector"),
          .Logic.swap)
