## METHODS FOR GENERIC: ! (not)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("!", c(x = "Matrix"),
          function(x) !.M2kind(x, "l"))

setMethod("!", c(x = "sparseVector"),
          function(x) !.V2kind(x, "l"))

setMethod("!", c(x = "ndenseMatrix"),
          function(x) {
              if ((shape <- .M.shape(x)) == "t")
                  x <- .M2gen(x)
              r <- new(.M.class(x))
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape == "s")
                  r@uplo <- x@uplo
              r@x <- { y <- x@x; !(y | is.na(y)) }
              r
          })

setMethod("!", c(x = "ldenseMatrix"),
          function(x) {
              if ((shape <- .M.shape(x)) == "t")
                  x <- .M2gen(x)
              r <- new(.M.class(x))
              r@Dim <- x@Dim
              r@Dimnames <- x@Dimnames
              if (shape == "s")
                  r@uplo <- x@uplo
              r@x <- !x@x
              r
          })

setMethod("!", c(x = "nsparseMatrix"),
          function(x) {
              if (.M.shape(x) == "t")
                  x <- .M2gen(x)
              x <- .sparse2dense(x)
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "lsparseMatrix"),
          function(x) {
              if (.M.shape(x) == "t")
                  x <- .M2gen(x)
              x <- .sparse2dense(x)
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "ndiMatrix"),
          function(x) {
              if (x@diag == "N" && anyNA(y <- x@x))
                  x@x <- y | is.na(y)
              x <- .diag2dense(x, ".", "g")
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "ldiMatrix"),
          function(x) {
              x <- .diag2dense(x, ".", "g")
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "indMatrix"),
          function(x) {
              x <- .ind2dense(x)
              x@x <- !x@x
              x
          })

setMethod("!", c(x = "nsparseVector"),
          function(x) !.V2v(x))

setMethod("!", c(x = "lsparseVector"),
          function(x) !.V2v(x))
