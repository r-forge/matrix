#### Index Matrices -- Coercion and Methods

## The typical   'constructor' : coerce from  'perm'
setAs("integer", "indMatrix",
      function(from) {
          n <- length(from)
          d <- length(unique(from))
          nn <- names(from)
          new("indMatrix", Dim = c(n, d), Dimnames = list(nn, NULL),
              perm = from)
      })

setAs("numeric", "indMatrix",
      function(from)
          if(all(from == (i <- as.integer(from)))) as(i, "indMatrix")
      else stop("coercion to \"indMatrix\" only works from integer numeric"))

setAs("indMatrix", "matrix",
      function(from) {
          fp <- from@perm
          r <- ldiag(n = from@Dim[2])[fp,]
          if(.has.DN(from)) dimnames(r) <- from@Dimnames
          r
      })


## coerce to 0/1 sparse matrix, i.e. sparse pattern
setAs("indMatrix", "ngTMatrix",
      function(from) {
          d <- from@Dim
          new("ngTMatrix", i = seq_len(d[1]) - 1L, j = from@perm - 1L,
              Dim = d, Dimnames = from@Dimnames)
      })

setAs("indMatrix", "TsparseMatrix", function(from) as(from, "ngTMatrix"))
setAs("indMatrix", "nMatrix",      function(from) as(from, "ngTMatrix"))
setAs("indMatrix", "lMatrix", function(from) as(as(from, "nMatrix"), "lMatrix"))
setAs("indMatrix", "dMatrix", function(from) as(as(from, "nMatrix"), "dMatrix"))
setAs("indMatrix", "dsparseMatrix", function(from) as(from, "dMatrix"))
setAs("indMatrix", "nsparseMatrix", function(from) as(from, "nMatrix"))
setAs("indMatrix", "CsparseMatrix",
      function(from) as(as(from, "ngTMatrix"), "CsparseMatrix"))
setAs("indMatrix", "ngeMatrix", function(from) as(as(from, "ngTMatrix"),"ngeMatrix"))

setAs("nMatrix", "indMatrix",
      function(from) {
          from <- as(as(from, "TsparseMatrix"), "ngTMatrix")
          n <- (d <- from@Dim)[1]
          if(n < d[2]) stop("not a skinny matrix")
          if(length(i <- from@i) != n)
              stop("the number of non-zero entries differs from nrow(.)")
          if((need.sort <- is.unsorted(i))) {
              ii <- sort.list(i)
              i <- i[ii]
          }
          if(n >= 1 && !identical(i, 0:(n - 1)))
              stop("must have exactly one non-zero entry per row")
          
          new("indMatrix", ## validity checking checks the 'perm' slot:
              perm = 1L + if(need.sort) from@j[ii] else from@j,
              Dim = d, Dimnames = from@Dimnames)
      })

setAs("matrix", "indMatrix", function(from) as(as(from, "nMatrix"), "indMatrix"))

setAs("indMatrix", "matrix", function(from) as(as(from, "nMatrix"), "matrix"))

setAs("sparseMatrix", "indMatrix", function(from)
    as(as(from, "nsparseMatrix"), "indMatrix"))
setMethod("is.na", signature(x = "indMatrix"), is.na_nsp)
setMethod("is.infinite", signature(x = "indMatrix"), is.na_nsp)
setMethod("is.finite", signature(x = "indMatrix"), allTrueMatrix)

setMethod("t", signature(x = "indMatrix"), function(x) t(as(x, "ngTMatrix")))



setMethod("%*%", signature(x = "matrix", y = "indMatrix"),
          function(x, y) {
              mmultCheck(x,y)
              x%*%as(y, "CsparseMatrix") 
          })
setMethod("%*%", signature(x = "Matrix", y = "indMatrix"),
          function(x, y) {
              mmultCheck(x,y)
              x%*%as(y, "CsparseMatrix") 
          })

setMethod("%*%", signature(x = "indMatrix", y = "matrix"),
          function(x, y) { mmultCheck(x,y); y[x@perm ,] })
setMethod("%*%", signature(x = "indMatrix", y = "Matrix"),
          function(x, y) { mmultCheck(x,y); y[x@perm ,] })


setMethod("crossprod", signature(x = "indMatrix", y = "matrix"),
          function(x, y) { mmultCheck(x,y, 2L); t(y%*%as(x, "CsparseMatrix"))})
setMethod("crossprod", signature(x = "indMatrix", y = "Matrix"),
          function(x, y) { mmultCheck(x,y, 2L); t(y%*%as(x, "CsparseMatrix"))})
setMethod("crossprod", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y) {
              mmultCheck(x,y, 2L)
              Matrix(data=tabulate(interaction(x@perm, y@perm)), 
                     nrow=x@Dim[2], 
                     ncol=y@Dim[2])
          })

setMethod("tcrossprod", signature(x = "matrix", y = "indMatrix"),
          function(x, y) { mmultCheck(x,y, 3L); x[, y@perm] })
setMethod("tcrossprod", signature(x = "Matrix", y = "indMatrix"),
          function(x, y) { mmultCheck(x,y, 3L); x[, y@perm] })
setMethod("tcrossprod", signature(x = "indMatrix", y = "indMatrix"),
          function(x, y) {
              mmultCheck(x,y, 3L)
              x[,y@perm]
          })


setMethod("crossprod", signature(x = "indMatrix", y = "missing"),
          function(x, y=NULL) Diagonal(x=as.numeric(table(x@perm))))
setMethod("tcrossprod", signature(x = "indMatrix", y = "missing"),
          function(x, y=NULL) x[,x@perm] )


setMethod("kronecker", signature(X = "indMatrix", Y = "indMatrix"),
          function (X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              if (FUN != "*") stop("kronecker method must use default 'FUN'")
              perm <- as.integer(interaction(rep(X@perm, each=Y@Dim[1]), 
                                   rep(Y@perm, times=X@Dim[1]), lex.order=TRUE))
              new("indMatrix", perm=perm, Dim=X@Dim*Y@Dim)
          })          


.indMat.nosense <- function (x, i, j, ..., value)
    stop('replacing "indMatrix" entries is not allowed, as rarely sensible')
setReplaceMethod("[", signature(x = "indMatrix", i = "perm"), .indMat.nosense)
setReplaceMethod("[", signature(x = "indMatrix", i = "missing", j = "perm"),
                 .indMat.nosense) ##   explicit  ^^^^^^^^^^^^ for disambiguation
setReplaceMethod("[", signature(x = "indMatrix", i = "missing", j = "missing"),
                 .indMat.nosense)