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
