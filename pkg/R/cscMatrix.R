setMethod("crossprod", signature(x = "cscMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_crossprod", x, PACKAGE = "Matrix"))

setMethod("crossprod", signature(x = "cscMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, y, PACKAGE = "Matrix"))

setMethod("crossprod", signature(x = "cscMatrix", y = "numeric"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", x, as.matrix(y), PACKAGE = "Matrix"))

setMethod("dim", signature(x = "cscMatrix"),
          function(x) x@Dim, valueClass = "integer")

setMethod("diag", signature(x = "cscMatrix"),
          function(x = 1, nrow, ncol = n)
          .Call("csc_getDiag", x, PACKAGE = "Matrix"))

setAs("cscMatrix", "tripletMatrix",
      function(from)
      .Call("csc_to_triplet", from, PACKAGE = "Matrix")
      )

setAs("cscMatrix", "matrix",
      function(from)
      .Call("csc_to_matrix", from, PACKAGE = "Matrix"))

setAs("cscMatrix", "geMatrix",
      function(from)
      .Call("csc_to_geMatrix", from, PACKAGE = "Matrix"))

setAs("matrix", "cscMatrix",
      function(from) {
          storage.mode(from) = "double"
          .Call("matrix_to_csc", from, PACKAGE = "Matrix")
      })

setMethod("t", signature(x = "cscMatrix"),
          function(x) .Call("csc_transpose", x, PACKAGE = "Matrix"),
          valueClass = "cscMatrix")

setMethod("image", "cscMatrix",
          function(x, ...) {
              dots = list(...)
              col = c("white", "black")
              if (!is.null(tmp <- dots$col)) col = tmp
              xlab = "column"
              if (!is.null(tmp <- dots$xlab)) xlab = tmp
              ylab = "row"
              if (!is.null(tmp <- dots$ylab)) ylab = tmp
              dd = dim(x)
              nr = dd[1]
              nc = dd[2]
              opar = par(las = 1)
              image(x = 1:nc, y = 1:nr,
                    z = .Call("csc_to_imagemat", x, PACKAGE = "Matrix"),
                    zlim = 0:1, axes = FALSE, col = col, las = 1,
                    xlab = xlab, ylab = ylab, ...)
              axis(1, pretty(c(1,nc)))
              axis(2, 1+nr-pretty(c(1,nr)), labels = pretty(c(1,nr)))
              box()
              par(opar)
          })

setMethod("image", "tripletMatrix",
          function(x,
                   xlim = c(0, matdim[2] + 1),
                   ylim = c(matdim[1] + 1, 0),
                   sub = sprintf("Dimensions: %d x %d", matdim[1], matdim[2]),
                   xlab = "Column", ylab = "Row",
                   cuts = 20,
                   col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
                   ...) {
              require(lattice)
              
              matdim <- x@Dim
              levelplot(abs(trip@x) ~ trip@j * trip@i,
                        sub = sub,
                        xlab = xlab, ylab = ylab,
                        xlim = xlim, ylim = ylim,
                        col.regions = col.regions,
                        par.settings = list(background = list(col = "transparent")),
                        ...)
          })
