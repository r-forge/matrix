setReplaceMethod("LMEoptimize", signature(x="ssclme", value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1) return(x)
                 st = coef(x, unconstr = TRUE) # starting values
                 if (value$optimizer == "optim") {
                     optimRes =
                         if (value$analyticGradient) {
                             optim(st,
                                   fn = function(pars) {
                                       coef(x, unconstr = TRUE) = pars
                                       deviance(x, REML = value$REML)
                                   },
                                   gr = function(pars) {
                                       coef(x, unconstr = TRUE) = pars
                                       gradient(x, REML = value$REML)
                                   },
                                   method = "BFGS",
                                   control = list(trace = value$msVerbose,
                                                  reltol = value$msTol,
                                                  maxit = value$msMaxIter))
                         } else {
                             optim(st,
                                   fn = function(pars) {
                                       coef(x, unconstr = TRUE) = pars
                                       deviance(x, REML = value$REML)
                                   },
                                   method = "BFGS",
                                   control = list(trace = value$msVerbose,
                                                  reltol = value$msTol,
                                                  maxit = value$msMaxIter))
                         }
                     if (optimRes$convergence != 0) {
                         warning("optim failed to converge")
                     }
                     coef(x, unconstr = TRUE) = optimRes$par
                 } else {
                     typsize <- rep(1.0, length(st))
                     if (is.null(value$nlmStepMax))
                         value$nlmStepMax <-
                             max(100 * sqrt(sum((st/typsize)^2)), 100)
                     nlmRes =
                         nlm(f = if (value$analyticGradient) {
                             function(pars) {
                                 coef(x, unconstr = TRUE) = pars
                                 ans = deviance(x, REML = value$REML)
                                 attr(ans, "gradient") =
                                     gradient(x, REML = value$REML)
                                 ans
                             }
                         } else {
                             function(pars)
                             {
                                 coef(x, unconstr = TRUE) = pars
                                 deviance(x, REML = value$REML)
                             }
                         },
                             p = st,
                             print.level = if (value$msVerbose) 2 else 0,
                             steptol = value$msTol,
                             gradtol = value$msTol,
                             stepmax = value$nlmStepMax,
                             typsize=typsize,
                             iterlim = value$msMaxIter)
                     coef(x, unconstr = TRUE) = nlmRes$estimate
                 }
                 return(x)
             })

setMethod("deviance", signature(object = "ssclme"),
          function(object, REML = FALSE, ...) {
              .Call("ssclme_factor", object, PACKAGE = "Matrix")
              object@deviance[ifelse(REML, 2, 1)]
          })

setMethod("coef", signature(object = "ssclme"),
          function(object, unconst = FALSE, ...) {
              .Call(ifelse(as(unconst, "logical")[1],
                           "ssclme_coefUnc", "ssclme_coef"),
                    object, PACKAGE = "Matrix")
          })

setMethod("ranef", signature(object = "ssclme"),
          function(object, ...) {
              val = .Call("ssclme_ranef", object, PACKAGE = "Matrix")
              bv = object@bVar
              names(val) = names(bv)
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = dimnames(bv[[i]])[-1]
              }
              lapply(val, t)
          })


setMethod("fixef", signature(object = "ssclme"),
          function(object, ...) {
              val = .Call("ssclme_fixef", object, PACKAGE = "Matrix")
              names(val) = dimnames(object@XtX)[[2]][seq(along = val)]
              val
          })

setMethod("vcov", signature(object = "ssclme"),
          function(object, ...) {
              sigma = .Call("ssclme_sigma", object, PACKAGE = "Matrix")
              rr = object@RXX
              nr = nrow(rr)
              rr = rr[-nr, -nr, drop = FALSE]
              sigma^2 * rr %*% t(rr)
          })

setMethod("VarCorr", signature(x = "ssclme"),
          function(x) {
              val = .Call("ssclme_variances", x, TRUE, PACKAGE = "Matrix")
              bVar = x@bVar
              for (i in seq(along = val[1]))
                  dimnames(val[1][[i]]) = dimnames(bVar[[i]][1:2])
              new("VarCorr", scale = val[[2]], reSumry = val[1],
                   useScale = TRUE)
          })
