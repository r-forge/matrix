setMethod("GLMM", signature(formula = "missing"),
          function(formula, family, data, random, ...)
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
          function(formula, family, data, random, ...)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("GLMM", signature(random = "formula"),
          function(formula, family, data, random, ...)
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
          function(formula, family, data, random, ...)
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
                   method = match.arg(c("PQL", "Laplace")),
                   control = list(),
                   subset,
                   weights,
                   na.action,
                   offset,
                   model = TRUE, x = FALSE, ...)
      {
          random <-
              lapply(random,
                     get("formula", pos = parent.frame(), mode = "function"))
          controlvals <- do.call("lmeControl", control)
          controlvals$REML <- FALSE


          ## problems with ..1, ..2
          ## print(match.call(expand.dots = FALSE))




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
          .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")




#          print(str(ranef(reducedObj)))






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
                      coef(reducedObj) <- pars[responseIndex:length(pars)]
                      off <- drop(mmats.unadjusted$.Xy %*%
                                  c(pars[1:(responseIndex-1)], 0))
                  }

                  niter <- 20
                  conv <- FALSE

                  eta <-
                      .Call("ssclme_fitted", reducedObj, facs,
                            reducedMmats.unadjusted, PACKAGE = "Matrix") + off
#                  eta.check <-
#                      .Call("ssclme_fitted", obj, facs,
#                            mmats.unadjusted, PACKAGE = "Matrix")
#                  print(all.equal(eta, eta.check))
#                  plot(eta, eta.check)
#                  etaold <- eta


                  for (iter in seq(length = niter))
                  {
                      mu <- family$linkinv(eta)
                      dmu.deta <- family$mu.eta(eta)
                      z <- eta + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta - offset
                      w <- weights * dmu.deta / sqrt(family$variance(mu))
                      .Call("nlme_weight_matrix_list",
                            reducedMmats.unadjusted, w, z, reducedMmats, PACKAGE="lme4")
                      .Call("ssclme_update_mm", reducedObj, facs, reducedMmats, PACKAGE="Matrix")
                      eta[] <-
                          .Call("ssclme_fitted", reducedObj, facs,
                                reducedMmats.unadjusted, PACKAGE = "Matrix") + off
                      cat(paste("bhat Criterion:", max(abs(eta - etaold)) /
                                (0.1 + max(abs(eta))), "\n"))
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
                  ## bhat doesn't return anything, it just modifies reducedObj
                  ## In particular, We are interested in ranef(reducedObj) and 
                  ## reducedObj@bVar (?)
                  invisible()
              }



          ## Get updated ranef estimates
          bhat()

          print(str(ranef(reducedObj)))


          new("lme", call = match.call(), facs = facs,
              x = if (x) mmats else list(),
              model = if(model) data else data.frame(list()),
              REML = FALSE, rep = obj, fitted = eta)
      })










### Local variables:
### mode: R
### End:
