setMethod("GLMM", signature(formula = "missing"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          resp = getResponseFormula(data)[[2]]
          cov = getCovariateFormula(data)[[2]]
          nCall$formula = eval(substitute(resp ~ cov))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("GLMM", signature(formula = "formula",
                            data = "groupedData", random = "missing"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("GLMM", signature(random = "formula"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          cov = getCovariateFormula(random)
          nCall$random <- lapply(getGroupsFormula(random, asList = TRUE),
                                 function(f) cov)
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("GLMM", signature(formula = "formula",
                            data = "groupedData",
                            random = "list"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          nCall = mCall = match.call()
          nCall$data <- data@data
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })





make.glm.call <- 
    function (mf, frm) 
{
    m <- match(c("formula", "family", "data", "weights", "subset", "na.action", 
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("glm")
    ## environment(frm) = environment(formula) ???
    mf$formula = frm
    mf$model = FALSE
    mf$x = FALSE
    mf$y = TRUE
    mf
}










setMethod("GLMM",
          signature(formula = "formula",
                    random = "list"),
          function(formula, family, data, random,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          random <-
              lapply(random,
                     get("formula", pos = parent.frame(), mode = "function"))
          controlvals <- do.call("lmeControl", control)
          controlvals$REML <- FALSE
          method <- match.arg(method)

          ## BEGIN glm fit without random effects

          ## several arguments are handled at this point.
          ## what about

          ##   subset, na.action (arguments to model.frame) 

          ## Are these  already handled ?
          ## The rest could be part of ..., why have them explicitly ?
          ## Do we want possibility of supplying control in glm() ?

          glm.fit <- eval(make.glm.call(match.call(expand.dots = TRUE),
                                        formula), parent.frame())
          family <- glm.fit$family
          offset <- if (is.null(glm.fit$offset)) 0 else glm.fit$offset
          weights <- sqrt(abs(glm.fit$prior.weights))
          ## initial 'fitted' values on linear scale
          eta <- glm.fit$linear.predictors
          etaold <- eta

          ## END using glm fit results
          ## Note: offset is on the linear scale

          ## Not clear how offset works






          data <- eval(make.mf.call(match.call(expand.dots = FALSE),
                                    formula, random), parent.frame())
          facs <- lapply(names(random),
                         function(x) eval(as.name(x), envir = data))
          names(facs) <- names(random)
          facs <-              # order in decreasing number of levels
              facs[rev(order(unlist(lapply(facs,
                                           function(fac)
                                           length(levels(fac))))))]
          ## creates model matrices
          mmats.unadjusted <-
              c(lapply(random,
                       function(x) model.matrix(formula(x), data = data)),
                list(.Xy =
                     cbind(model.matrix(formula, data = data),
                           .response = glm.fit$y))) #  WAS: model.response(data)
          responseIndex <- ncol(mmats.unadjusted$.Xy)
          obj <-  ## creates ssclme structure
              .Call("ssclme_create", facs,
                    unlist(lapply(mmats.unadjusted, ncol)),
                    as.integer(controlvals$maxLIsize),
                    PACKAGE = "Matrix")
          facs = facshuffle(obj, facs)
          obj = obj[[1]]
          mmats <- mmats.unadjusted
          ## the next line is to force a copy of mmats, because we are
          ## going to use both mmats and mmats.unadjusted as arguments
          ## in a .Call where one of them will be modified (don't want
          ## the other to be modified as well)
          mmats[[1]][1,1] <- 1 
          conv <- FALSE
          firstIter <- TRUE

          for (iter in seq(length = controlvals$glmmMaxIter))
          {
              mu <- family$linkinv(eta)
              dmu.deta <- family$mu.eta(eta)
              ## adjusted response
              z <- eta + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta - offset
              ## weights (note: weights is already square-rooted)
              w <- weights * dmu.deta / sqrt(family$variance(mu))
              .Call("nlme_weight_matrix_list",
                    mmats.unadjusted, w, z, mmats, PACKAGE="lme4")
              .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
              if (firstIter) .Call("ssclme_initial", obj, PACKAGE="Matrix")
              .Call("ssclme_EMsteps", obj,
                    controlvals$niterEM,
                    FALSE, #controlvals$REML,
                    controlvals$EMverbose,
                    PACKAGE = "Matrix")
              LMEoptimize(obj) = controlvals
              eta[] <-
                  .Call("ssclme_fitted", obj, facs,
                        mmats.unadjusted, PACKAGE = "Matrix")
              cat(paste("Iteration", iter,"Termination Criterion:",
                        format(max(abs(eta - etaold)) /
                        (0.1 + max(abs(eta)))), "\n"))
              ## use this to determine convergence
              if (max(abs(eta - etaold)) <
                  (0.1 + max(abs(eta))) * controlvals$tolerance)
              {
                  conv <- TRUE
                  break
              }
              etaold[] <- eta

              ## Changing number of iterations on second and
              ## subsequent iterations.
              if (firstIter)
              {
                  controlvals$niterEM <- 2
                  controlvals$msMaxIter <- 10
                  firstIter <- FALSE
              }
          }
          if (!conv) warning("IRLS iterations for glmm did not converge")







          ## Need to optimize L(theta, beta) using Laplace approximation

          ## Things needed for that:
          ##
          ## 1. reduced ssclme object, offset, weighted model matrices
          ## 2. facs, reduced model matrices

          ## Of these, those in 2 will be fixed given theta and beta,
          ## and can be thought of arguments to the L(theta, beta)
          ## function. However, the ones in 1 will have the same
          ## structure. So the plan is to pre-allocate them and pass
          ## them in too so they can be used without creating/copying
          ## them more than once


          ## reduced ssclme

          reducedObj <- .Call("ssclme_collapse", obj, PACKAGE = "Matrix")
          reducedMmats.unadjusted <- mmats.unadjusted
          reducedMmats.unadjusted$.Xy <-
              reducedMmats.unadjusted$.Xy[, responseIndex, drop = FALSE]
          reducedMmats <- mmats
          reducedMmats$.Xy <-
              reducedMmats$.Xy[, responseIndex, drop = FALSE]
          .Call("ssclme_update_mm", reducedObj, facs, reducedMmats, PACKAGE="Matrix")

          ## make obj comparable
          ## .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")






          ## define function that calculates bhats given theta and beta 

          bhat <- 
              function(pars = NULL) # 1:(responseIndex-1) - beta, rest - theta
              {
                  if (is.null(pars))
                  {
                      off <- drop(mmats.unadjusted$.Xy %*% c(fixef(obj), 0))
                  }
                  else
                  {
#                      cat(pars)
#                      cat("\n")
#                      print(reducedObj@Omega) 
                      coef(reducedObj) <- pars[responseIndex:length(pars)]
#                      print(reducedObj@Omega)
#                      print("-----------------------")
                      off <- drop(mmats.unadjusted$.Xy %*%
                                  c(pars[1:(responseIndex-1)], 0))
                  }

                  niter <- 20
                  conv <- FALSE

                  eta <- off + 
                      .Call("ssclme_fitted", reducedObj, facs,
                            reducedMmats.unadjusted, PACKAGE = "Matrix")
#                  eta <-
#                      .Call("ssclme_fitted", obj, facs,
#                            mmats.unadjusted, PACKAGE = "Matrix")
#                  print(all.equal(eta, eta.check)) # not really sure why not TRUE
#                  plot(eta, eta.check)
                  etaold <- eta


                  for (iter in seq(length = niter))
                  {
                      mu <- family$linkinv(eta)
                      dmu.deta <- family$mu.eta(eta)
                      z <- eta + (reducedMmats.unadjusted$.Xy[, 1] - mu) /
                          dmu.deta - offset
                      w <- weights * dmu.deta / sqrt(family$variance(mu))
                      .Call("nlme_weight_matrix_list",
                            reducedMmats.unadjusted, w, z, reducedMmats, PACKAGE="lme4")
                      .Call("ssclme_update_mm", reducedObj, facs, reducedMmats, PACKAGE="Matrix")
                      eta[] <- off + 
                          .Call("ssclme_fitted", reducedObj, facs,
                                reducedMmats.unadjusted, PACKAGE = "Matrix")
                      ##cat(paste("bhat Criterion:", max(abs(eta - etaold)) /
                      ##          (0.1 + max(abs(eta))), "\n"))
                      ## use this to determine convergence
                      if (max(abs(eta - etaold)) <
                          (0.1 + max(abs(eta))) * controlvals$tolerance)
                      {
                          conv <- TRUE
                          break
                      }
                      etaold[] <- eta

                  }
                  if (!conv) warning("iterations for bhat did not converge")

                  ## bhat doesn't really need to return anything, we
                  ## just want the side-effect of modifying reducedObj
                  ## In particular, We are interested in
                  ## ranef(reducedObj) and reducedObj@bVar (?). But
                  ## the mu-scale response will be useful for log-lik
                  ## calculations later, so return them anyway

                  invisible(family$linkinv(eta)) 
              }


          ## function that calculates log likelihood (the thing that
          ## needs to get evaluated at each Gauss-Hermite location)

          ## log scale ? worry about details later, get the pieces in place

          ## this is for the Laplace approximation only. GH is more
          ## complicated 

          loglikLaplace <- function(pars = NULL)
          {
              mu <- bhat(pars = pars)  ## gets correct values of bhat and bvars

              ## GLM family log likelihood (upto constant ?)(log scale)
              ## FIXME: need to adjust for sigma^2 for appropriate models (not trivial!)
              ans <- ## log lik from observations given fixed and random effects
                  sum(family$dev.resids(y = mmats.unadjusted$.Xy[, responseIndex],
                                        mu = mu,
                                        wt = weights^2))
              ranefs <- ranef(reducedObj)
              Omega <- obj@Omega
              jacobian <- numeric(length(ranefs))
              for (i in seq(along = ranefs))
              {
                  Omega[[i]] <- Omega[[i]] + t(Omega[[i]])
                  diag(Omega[[i]]) <- diag(Omega[[i]]) / 2
                  ## want log of `const det(Omega) exp(-1/2 b' Omega b )`
                  ## i.e., const + log det(Omega) - .5 * (b' Omega b)
                  ## FIXME: need to adjust for sigma^2 for appropriate models (easy)
                  ## these are all the b'Omega b, summed as they eventually need to be
                  ## think of this as sum(rowSums((ranefs[[i]] %*% Omega[[i]]) * ranefs[[i]]))
                  ans <- ans - sum((ranefs[[i]] %*% Omega[[i]]) * ranefs[[i]]) +
                      nrow(ranefs[[i]]) * determinant(Omega[[i]], logarithm = TRUE)$modulus
                  ## Jacobian adjustment
                  ## log.jacobian[i] <- sum(log(apply(bVar[[i]], 3, function(x) sum(diag(x)))))
                  ans <- ans + sum(log(apply(reducedObj@bVar[[i]], 3, function(x) sum(diag(x)))))
              }
              ans ## this is (upto some constant) log of the Laplacian approximation of the likelihood
          }

          if (method == "Laplace")
          {
              cat(paste("Using optimizer", controlvals$optim), "\n")

              if (controlvals$optimizer == "optim") {
                  optimRes =
                      ## if (controlvals$analyticGradient)  ??
                      optim(fn = loglikLaplace,
                            par = c(fixef(obj), coef(obj)),
                            ## hessian = TRUE,
                            method = "BFGS",
                            control = list(trace = TRUE, #controlvals$msVerbose,
                            reltol = controlvals$msTol,
                            ##fnscale = -xval,
                            ##parscale = 1/controlvals$msScale(coef(obj)),
                            maxit = controlvals$msMaxIter))
                  if (optimRes$convergence != 0) {
                      warning("optim failed to converge")
                  }
                  ##fixef(obj) <- optimRes$par[seq(length = responseIndex - 1)]
                  print(fixef(obj))
                  print(optimRes$par[seq(length = responseIndex - 1)])
                  print(coef(obj))
                  coef(obj) <- optimRes$par[responseIndex:length(optimRes$par)]
                  print(coef(obj))
                  ## need to calculate likelihood
              }
              else if (controlvals$optimizer == "nlm") {
                  ## typsize <- 1/controlvals$msScale(coef(obj))
                  typsize <- rep(1.0, length(coef(obj)) + responseIndex - 1)
                  if (is.null(controlvals$nlmStepMax))
                      controlvals$nlmStepMax <-
                          max(100 * sqrt(sum((c(fixef(obj),
                                                coef(obj))/typsize)^2)), 100)
                  nlmRes =
                      nlm(f = loglikLaplace, ## if (controlvals$analyticGradient)  ??
                          p = c(fixef(obj), coef(obj)),
                          ## hessian = TRUE,
                          print.level = 2, #if (controlvals$msVerbose) 2 else 0,
                          steptol = controlvals$msTol,
                          gradtol = controlvals$msTol,
                          stepmax = controlvals$nlmStepMax,
                          typsize=typsize,
                          ## fscale=xval,
                          iterlim = controlvals$msMaxIter)
                  ##fixef(obj) <- nlmRes$estimate[seq(length = responseIndex - 1)]
                  print(fixef(obj))
                  print(nlmRes$estimate[seq(length = responseIndex - 1)])
                  print(coef(obj))
                  coef(obj) <- nlmRes$estimate[responseIndex:length(nlmRes$estimate)]
                  print(coef(obj))
                  ## need to calculate likelihood
              }

          }



          
          ## Get updated ranef estimates
#          bhat()

#          str(ranef(reducedObj))


          new("lme", call = match.call(), facs = facs,
              x = if (x) mmats else list(),
              model = if(model) data else data.frame(list()),
              REML = FALSE, rep = obj, fitted = eta)
      })








### Local variables:
### mode: R
### End:
