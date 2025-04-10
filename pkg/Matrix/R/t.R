## METHODS FOR GENERIC: t, ct
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("ct", c(x = "matrix"),
          function(x)
              if (is.complex(x)) Conj(t.default(x)) else t.default(x))

setMethod( "t", c(x = "denseMatrix"),
          function(x)
              .Call(R_dense_transpose, x, "T"))

setMethod("ct", c(x = "denseMatrix"),
          function(x)
              .Call(R_dense_transpose, x, "C"))

for (.cl in paste0(c("C", "R", "T"), "sparseMatrix")) {
setMethod( "t", c(x = .cl),
          function(x)
              .Call(R_sparse_transpose, x, "T", FALSE))

setMethod("ct", c(x = .cl),
          function(x)
              .Call(R_sparse_transpose, x, "C", FALSE))
}

setMethod( "t", c(x = "diagonalMatrix"),
          function(x) {
              x@Dimnames <- x@Dimnames[2:1]
              x
          })

setMethod("ct", c(x = "diagonalMatrix"),
          function(x) {
              x@Dimnames <- x@Dimnames[2:1]
              if (is.complex(y <- x@x))
                  x@x <- Conj(y)
              x
          })

.fn <-
function(x) {
    r <- new("indMatrix")
    r@Dim <- x@Dim[2:1]
    r@Dimnames = x@Dimnames[2:1]
    r@perm <- x@perm
    if (x@margin == 1L)
        r@margin <- 2L
    r
}

setMethod( "t", c(x = "indMatrix"), .fn)
setMethod("ct", c(x = "indMatrix"), .fn)

.fn <-
function(x) {
    r <- new("pMatrix")
    r@Dim <- x@Dim
    r@Dimnames = x@Dimnames[2:1]
    r@perm <- x@perm
    if (x@margin == 1L)
        r@margin <- 2L
    r
}

setMethod( "t", c(x = "pMatrix"), .fn)
setMethod("ct", c(x = "pMatrix"), .fn)

setMethod( "t", c(x = "sparseVector"),
          function(x) .tCRT(.V2C(x), trans = "T"))
setMethod("ct", c(x = "sparseVector"),
          function(x) .tCRT(.V2C(x), trans = "C"))

rm(.cl, .fn)
