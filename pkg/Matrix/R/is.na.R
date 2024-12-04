## METHODS FOR GENERIC: anyNA
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("anyNA", c(x = "denseMatrix"),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && anyNA(forceCanonical(x)@x))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("anyNA", c(x = .cl),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && anyNA(forceCanonical(x)@x))

setMethod("anyNA", c(x = "diagonalMatrix"),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && length(y <- x@x) > 0L && anyNA(y))

setMethod("anyNA", c(x = "indMatrix"),
          function(x, recursive = FALSE)
              FALSE)

setMethod("anyNA", c(x = "sparseVector"),
          function(x, recursive = FALSE)
              .M.kind(x) != "n" && anyNA(x@x))

rm(.cl)


## METHODS FOR GENERIC: is.na
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.na", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- kind != "n"
              substr(cl, 1L, 1L) <- "n"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g") {
                  r@uplo <- x@uplo
                  if (if (shape == "s") kind == "z" && x@trans == "C" else maybe && x@diag != "N")
                      x <- forceCanonical(x)
              }
              r@x <- if (maybe) is.na(x@x) else logical(length(x@x))
              r
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("is.na", c(x = .cl),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- kind != "n"
              substr(cl, 1L, 1L) <- if (maybe) "l" else "n"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g") {
                  r@uplo <- x@uplo
                  if (shape == "s" && kind == "z" && x@trans == "C")
                      x <- forceCanonical(x)
              }
              if (maybe) {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { r@i <- x@i; r@j <- x@j })
                  r@x <- is.na(x@x)
                  .M2kind(.drop0(r), "n")
              } else {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- integer(d[2L] + 1) },
                         "R" = { r@p <- integer(d[1L] + 1) })
                  r
              }
          })

setMethod("is.na", c(x = "diagonalMatrix"),
          function(x) {
              maybe <- .M.kind(x) != "n" && x@diag == "N"
              r <- new("ndiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if (maybe) is.na(x@x) else logical(d[1L])
              r
          })

setMethod("is.na", c(x = "indMatrix"),
          function(x) {
              m <- x@margin
              r <- new(if (m == 1L) "ngRMatrix" else "ngCMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@p <- integer(d[m] + 1)
              r
          })

setMethod("is.na", c(x = "sparseVector"),
          function(x) {
              maybe <- .M.kind(x) != "n"
              r <- new("nsparseVector")
              r@length <- x@length
              if (maybe)
                  r@i <- x@i[is.na(x@x)]
              r
          })

rm(.cl)


## METHODS FOR GENERIC: is.nan
## NB: mostly parallel to is.na, completely parallel to is.infinite
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.nan", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- any(kind == c("d", "z"))
              substr(cl, 1L, 1L) <- "n"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g") {
                  r@uplo <- x@uplo
                  if (if (shape == "s") kind == "z" && x@trans == "C" else maybe && x@diag != "N")
                      x <- forceCanonical(x)
              }
              r@x <- if (maybe) is.nan(x@x) else logical(length(x@x))
              r
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("is.nan", c(x = .cl),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- any(kind == c("d", "z"))
              substr(cl, 1L, 1L) <- if (maybe) "l" else "n"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g") {
                  r@uplo <- x@uplo
                  if (shape == "s" && kind == "z" && x@trans == "C")
                      x <- forceCanonical(x)
              }
              if (maybe) {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { r@i <- x@i; r@j <- x@j })
                  r@x <- is.nan(x@x)
                  .M2kind(.drop0(r), "n")
              } else {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- integer(d[2L] + 1) },
                         "R" = { r@p <- integer(d[1L] + 1) })
                  r
              }
          })

setMethod("is.nan", c(x = "diagonalMatrix"),
          function(x) {
              maybe <- any(.M.kind(x) == c("d", "z")) && x@diag == "N"
              r <- new("ndiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if (maybe) is.nan(x@x) else logical(d[1L])
              r
          })

setMethod("is.nan", c(x = "indMatrix"),
          function(x) {
              m <- x@margin
              r <- new(if (m == 1L) "ngRMatrix" else "ngCMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@p <- integer(d[m] + 1)
              r
          })

setMethod("is.nan", c(x = "sparseVector"),
          function(x) {
              maybe <- any(.M.kind(x) == c("d", "z"))
              r <- new("nsparseVector")
              r@length <- x@length
              if (maybe)
                  r@i <- x@i[is.nan(x@x)]
              r
          })

rm(.cl)


## METHODS FOR GENERIC: is.infinite
## NB: mostly parallel to is.na, completely parallel to is.nan
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.infinite", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- any(kind == c("d", "z"))
              substr(cl, 1L, 1L) <- "n"
              r <- new(cl)
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g") {
                  r@uplo <- x@uplo
                  if (if (shape == "s") kind == "z" && x@trans == "C" else maybe && x@diag != "N")
                      x <- forceCanonical(x)
              }
              r@x <- if (maybe) is.infinite(x@x) else logical(length(x@x))
              r
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("is.infinite", c(x = .cl),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- any(kind == c("d", "z"))
              substr(cl, 1L, 1L) <- if (maybe) "l" else "n"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape != "g") {
                  r@uplo <- x@uplo
                  if (shape == "s" && kind == "z" && x@trans == "C")
                      x <- forceCanonical(x)
              }
              if (maybe) {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- x@p; r@i <- x@i },
                         "R" = { r@p <- x@p; r@j <- x@j },
                         "T" = { r@i <- x@i; r@j <- x@j })
                  r@x <- is.infinite(x@x)
                  .M2kind(.drop0(r), "n")
              } else {
                  switch(substr(cl, 3L, 3L),
                         "C" = { r@p <- integer(d[2L] + 1) },
                         "R" = { r@p <- integer(d[1L] + 1) })
                  r
              }
          })

setMethod("is.infinite", c(x = "diagonalMatrix"),
          function(x) {
              maybe <- any(.M.kind(x) == c("d", "z")) && x@diag == "N"
              r <- new("ndiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- if (maybe) is.infinite(x@x) else logical(d[1L])
              r
          })

setMethod("is.infinite", c(x = "indMatrix"),
          function(x) {
              m <- x@margin
              r <- new(if (m == 1L) "ngRMatrix" else "ngCMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@p <- integer(d[m] + 1)
              r
          })

setMethod("is.infinite", c(x = "sparseVector"),
          function(x) {
              maybe <- any(.M.kind(x) == c("d", "z"))
              r <- new("nsparseVector")
              r@length <- x@length
              if (maybe)
                  r@i <- x@i[is.infinite(x@x)]
              r
          })

rm(.cl)


## METHODS FOR GENERIC: is.finite
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("is.finite", c(x = "denseMatrix"),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- kind != "n"
              if (shape == "s")
                  substr(cl, 1L, 1L) <- "n"
              else substr(cl, 1L, 3L) <- "nge"
              r <- new(cl)
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape == "s") {
                  r@uplo <- x@uplo
                  if (kind == "z" && x@trans == "C")
                      x <- forceCanonical(x)
              }
              else if (maybe && shape == "t")
                  x <- .M2gen(x)
              r@x <- if (maybe) is.finite(x@x) else rep.int(TRUE, if (shape == "s") length(x@x) else prod(d))
              r
          })

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix"))
setMethod("is.finite", c(x = .cl),
          function(x) {
              cl <- .M.class(x)
              kind  <- substr(cl, 1L, 1L)
              shape <- substr(cl, 2L, 2L)
              maybe <- kind != "n"
              r <- new(if (shape == "s") "nsyMatrix" else "ngeMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape == "s") {
                  r@uplo <- x@uplo
                  if (kind == "z" && x@trans == "C")
                      x <- forceCanonical(x)
              }
              tmp <- rep.int(TRUE, prod(d))
              if (maybe && !all(k <- is.finite(x@x))) {
                  if (substr(cl, 3L, 3L) != "T") {
                      x <- .M2T(x)
                      if (length(k) > length(x@x)) # was overallocated
                          k <- is.finite(x@x)
                  }
                  if (prod(d) > .Machine[["integer.max"]])
                      d <- as.double(d)
                  i <- x@j * d[1L] + x@i + 1L
                  tmp[i] <- k
              }
              r@x <- tmp
              r
          })

setMethod("is.finite", c(x = "diagonalMatrix"),
          function(x) {
              maybe <- .M.kind(x) != "n" && x@diag == "N"
              r <- new("nsyMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              tmp <- rep.int(TRUE, prod(d))
              if (maybe && !all(k <- is.finite(x@x))) {
                  i <- seq.int(from = 1L, by = d[1L] + 1, length.out = d[1L])
                  tmp[i] <- k
              }
              r@x <- tmp
              r
          })

setMethod("is.finite", c(x = "indMatrix"),
          function(x) {
              r <- new("ngeMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- rep.int(TRUE, prod(d))
              r
          })

setMethod("is.finite", c(x = "sparseVector"),
          function(x) {
              maybe <- .M.kind(x) != "n"
              r <- rep.int(TRUE, x@length)
              if (maybe)
                  r[x@i[!is.finite(x@x)]] <- FALSE
              r
          })

rm(.cl)
