## METHODS FOR GENERIC: cov2cor
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("cov2cor", c(V = "Matrix"),
          function(V) {
              if (.M.kind(V) == "z")
                  stop("matrix is not numeric")
              V <- as(forceSymmetric(V), "posdefMatrix")
              V@factors <- list()
              Vjj <- diag(V, names = FALSE)
              if ((n <- length(Vjj)) > 0L) {
              if (min(1, Vjj, na.rm = TRUE) <= 0)
                  stop("matrix is not positive definite")
              sd <- sqrt(Vjj)
              repr <- .M.repr(V)
              switch(repr,
                     "n" =
                         {
                             k <- seq.int(from = 1L, by = n + 1, length.out = n)
                             y <- V@x / sd / rep(sd, each = n)
                         },
                     "p" =
                         {
                             if (V@uplo == "U") {
                                 r <- 1L:n
                                 s <- 1L
                                 k <- cumsum(r)
                             } else {
                                 r <- n:1L
                                 s <- 1L:n
                                 k <- cumsum(c(1L, r[-n]))
                             }
                             y <- V@x / rep.int(sd, r) / sd[sequence.default(r, s)]
                         },
                     "C" =, "R" =
                         {
                             T <- .M2T(V)
                             y <- T@x / sd[T@i + 1L] / sd[T@j + 1L]
                             k <- T@i == T@j
                             ## V may be overallocated; we need to
                             ## ensure that the 'x' and 'i' or 'j'
                             ## slots continue to have equal length.
                             if (repr == "C")
                                 V@i <- T@i
                             else V@j <- T@j
                         },
                     "T" =
                         {
                             V <- aggregateT(V)
                             y <- V@x / sd[V@i + 1L] / sd[V@j + 1L]
                             k <- T@i == T@j
                         })
              }
              y[k] <- 1
              V@x <- y
              V
          })
