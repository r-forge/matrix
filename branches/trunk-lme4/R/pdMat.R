# Additional methods for pdMat objects

setMethod("LMEgradient",
          signature(x="pdDiag", A="matrix", nlev="numeric"),
          function(x, A, nlev) {
              .Call("pdDiag_LMEgradient", x, A, nlev, PACKAGE="Matrix")
          })

setMethod("LMEhessian",
          signature(x="pdDiag", A="matrix", H="array",
                    nlev="numeric"),
          function(x, A, H, nlev)
      {
          theta <- x@param
          q <- length(theta)

          ## the part involving D_iD_j
          ans <-
              if (q > 1)
                  diag(-exp(2*theta)*colSums(A*A))
              else as.matrix(-exp(2*theta)*A*A)

          ## add the part not involving D_iD_j
          tmp <- matrix(0, nrow=q, ncol=q)
          for (j in 1:q)
              tmp[, j] <- H[seq(from = 1+(j-1)*q*q*(q+1), length=q, by = q+1)]
          ans <- ans + tmp*outer(theta, theta,
                                 function(v, w) exp(2*(v+w)))
          nm <- names(x)
          ans <- ans+t(ans)
          if (!is.null(nm))
              dimnames(ans) <- list(nm, nm)
          ans
      })

setReplaceMethod("EMupdate",
                 signature(x="pdDiag", nlev="numeric", value="matrix"),
                 function(x, nlev, value) {
                     .Call("pdDiag_EMupdate", x, nlev, value, PACKAGE="Matrix")
                 })

setMethod("LMEgradient",
          signature(x="pdIdent", A="matrix", nlev="numeric"),
          function(x, A, nlev) {
              .Call("pdIdent_gradient", x, A, nlev, PACKAGE="lme4")
          })

setReplaceMethod("EMupdate",
                 signature(x="pdIdent", nlev="numeric", value="matrix"),
                 function(x, nlev, value) {
                     .Call("pdIdent_EMupdate", x, nlev, value, PACKAGE="Matrix")
                 })

setMethod("pdgradient", "pdIdent",
          function(x) {
              mat <- as(x, "pdmatrix")
              dn <- dimnames(mat)
              if (!is.null(dn)) dn <- c(list(dimnames(mat)), NULL)
              array(2. * c(mat), dim = c(dim(mat), 1), dimnames = dn)
          })

setMethod("LMEgradient",
          signature(x="pdLogChol", A="matrix", nlev="numeric"),
          function(x, A, nlev)
          .Call("pdLogChol_LMEgradient", x, A, nlev, PACKAGE="Matrix")
          )

setMethod("LMEhessian",
          signature(x="pdLogChol", A="matrix", H="array",
                    nlev="numeric"),
          function(x, A, H, nlev)
      {
          .Call("pdLogChol_LMEhessian", x, A, H, nlev, PACKAGE="Matrix")
      })

setReplaceMethod("EMupdate",
                 signature(x="pdLogChol", nlev="numeric", value="matrix"),
                 function(x, nlev, value) {
                     .Call("pdLogChol_EMupdate", x, nlev, value, PACKAGE="Matrix")
                 })

setMethod("pdgradient", "pdLogChol",
          function(x) {
              .Call("pdLogChol_pdgradient", x, PACKAGE="Matrix")
          })

