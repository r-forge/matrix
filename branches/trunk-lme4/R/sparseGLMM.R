


testfun <- function()
{
    ## simulation: 300 obs, 30 students, 10 obs per student, one covariate

    dat <- data.frame(id = gl(30, 10),
                      x = rnorm(300))

    dat$resp <-
        with(dat,
             rbinom(300,
                    size = 1,
                    prob = binomial()$linkinv( x + rnorm(30, sd = .5)[id] + rnorm(300, sd = .1) )))


    fm.bin <-
        sparseGLMM(resp ~ x, data = dat, family = binomial(),
                   random = list(id = ~1),
                   control = lmeControl(EMverbose = F))


    return (fm.bin)
}





setGeneric("sparseGLMM",
           function(formula, family, data, random,

                    subset, method, na.action, control, model, x, ...)

           standardGeneric("sparseGLMM"))



# if (!isGeneric("sparseGLMM")) {
#     setGeneric("sparseGLMM",
#                function(formula, family, data, random,

#                         control, niter, method, verbose, ...)


#                standardGeneric("sparseGLMM"))
# }






## lme for reference (remove later)

setMethod("sparseGLMM", signature(formula = "formula", family = "family", random = "list"),
          function(formula, family, data, random, subset,
                   method, na.action, control, model, x, ...)
      {

          ##print("got here")
          
          if (missing(model))
              model = TRUE
          if (missing(x))
              x = FALSE
          random = lapply(as(random, "list"),
                   get("formula", pos = parent.frame(), mode = "function"))
                   #lapply(random, function(x)
                          #if(inherits(x, "formula")) pdLogChol(x) else x)
          method = if (missing(method)) "REML" else
                   match.arg(method, c("REML", "ML"))
          controlvals <- if (missing(control)) lmeControl() else
                            do.call("lmeControl", control)
          mCall <- match.call(expand.dots = FALSE)


          mCall[[1]] <- as.name("model.frame")
          names(mCall)[2] <- "formula"
          mCall$family <- mCall$random <- mCall$method <-
              mCall$control <- mCall$model <- mCall$x <- NULL
          form <- formula
          form[[3]] <- (~a+b)[[2]]
          form[[3]][[2]] <- formula[[3]]
          form[[3]][[3]] <-
              as.formula((parse(text=paste("~",
                                paste(names(random),
                                      collapse = "+")))[[1]]))[[2]]
          for (pdm in random) {
              tmp <- form
              tmp[[3]] <- (~a+b)[[2]]
              tmp[[3]][[2]] <- form[[3]]
              tmp[[3]][[3]] <- formula(pdm)[[2]]
              form <- tmp
          }
          environment(form) <- environment(formula)
          mCall$formula <- form
          mCall$drop.unused.levels <- TRUE

          data <- eval(mCall, parent.frame())
          facs <- lapply(names(random),
                         function(x) eval(as.name(x), envir = data))
          names(facs) <- names(random)


          ## creates model matrices
          mmats.unadjusted <-
              c(lapply(random,
                       function(x) model.matrix(formula(x), data = data)),
                list(.Xy =
                     cbind(model.matrix(formula, data = data),
                           .response = model.response(data))))
          responseIndex <- ncol(mmats.unadjusted$.Xy)


          ## creates ssclme structure
          obj <- .Call("ssclme_create", facs, unlist(lapply(mmats.unadjusted, ncol)),
                       as.integer(2e5), PACKAGE = "Matrix")
          facs = facshuffle(obj, facs)
          obj = obj[[1]]




          ## get initial estimates
          fm.glm <- glm(formula, family, data)
          coefFixed <- c(coef(fm.glm), 0)
          ## what happens for more than one random effect ?
          coefRanef <- rep(0, nlevels(facs[[1]]))





#          print(str(mmats.unadjusted))





          mmats <- mmats.unadjusted
          niter <- 20

          for (iter in seq(length = niter))
          {
              eta <- mmats.unadjusted$.Xy %*% coefFixed
              for (facname in names(facs))
              {
                  ## FIXME: coefRanef depends on facname
                  eta <- eta + mmats.unadjusted[[facname]] * 
                      coefRanef[facs[[facname]]]
              }
              mu <- family$linkinv(eta)
              deta.dmu <- 1 / family$mu.eta(eta)

              ## adjusted response
              z <- eta + (mmats.unadjusted$.Xy[, responseIndex] - mu) * deta.dmu
              ## weights (needs sqrt ?)
              w <- 1 / (deta.dmu^2 * family$variance(mu))

              plot(z, mmats$.Xy[, responseIndex])

              ## simple option:  mmats <- mmats.unadjusted
              ## Does this prevent overwriting of components ?
              for (facname in names(facs))
                  mmats[[facname]][,] <- mmats.unadjusted[[facname]] * w
              mmats$.Xy[,] <- mmats.unadjusted$.Xy
              mmats$.Xy[, responseIndex] <- z
              mmats$.Xy[,] <- mmats$.Xy * w


              .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
              .Call("ssclme_initial", obj, PACKAGE="Matrix")
              .Call("ssclme_factor", obj, PACKAGE = "Matrix")
              .Call("ssclme_EMsteps", obj, controlvals$niterEM,
                    method == "REML", controlvals$EMverbose, PACKAGE = "Matrix")

              coefFixed <- c(.Call("ssclme_fixef", obj, PACKAGE = "Matrix"), 0)
              ## what happens for more than one random effect ?
              coefRanef <- .Call("ssclme_ranef", obj, PACKAGE = "Matrix")
          }

          new("lme", call = match.call(), facs = facs,
              x = if(x) mmats else list(),
              model = if(model) data else data.frame(list()),
              REML = method == "REML", rep = obj)
      })













## trying to modify

# setMethod("sparseGLMMold", signature(formula = "formula", random = "list"),
#           function(formula, family, data, random, control, niter,
#                    method, verbose, nEM.IRLS, model, x, ...)
#       {
#           if (missing(nEM.IRLS))
#               nEM.IRLS <- 1
#           if (missing(verbose)) verbose = FALSE
#           m <- Call <- match.call()
#           method <- if(missing(method)) "PQL" else
#                     match.arg(method, c("PQL", "Laplace"))
#           nm <- names(m)[-1]
#           dontkeep <-
#               is.element(nm, c("correlation", "control", "niter",
#                                "verbose", "nEM.IRLS", "method"))
#           for(i in nm[dontkeep]) m[[i]] <- NULL
#           m[[1]] <- as.name("glmmStruct")
#           m$nextraCols <- 1
#           m$method <- method
#           fit <- eval(m, parent.frame())

#           off <- fit@reStruct@offset
#           w <-  fit@prior.weights
#           origy <- fit@origy
#           fam <- fit@family

#           ## We always do the PQL steps. For that, we extract the reStruct
#           ## slot from fit
#           fit <- as(fit, "reStruct")
#           control = if (missing(control)) lmeControl() else
#                     do.call("lmeControl", control)
#           EMsteps(fit) <- control

#           control$niterEM <- nEM.IRLS
#           converged <- FALSE
#           eta <- .Call("nlme_reStruct_fitted", fit, NULL, PACKAGE="lme4")

#           for(i in seq(length=if(missing(niter)) 20 else niter)) {
#               ##update zz and wz
#               mu <- fam$linkinv(eta)
#               mu.eta.val <- fam$mu.eta(eta)
#               zz <- eta + (origy - mu)/mu.eta.val  - off
#               wz <- w * mu.eta.val^2 / fam$variance(mu)

#               response(fit) <- zz
#               weighted(fit) <- sqrt(abs(wz))
#               EMsteps(fit) <- control
#               LMEoptimize(fit) <- control
#               if(verbose) {
#                   cat("iteration", i, "\n")
# ###             class(fit) <- "glmmStruct"
# ###             fit@logLik <- as.numeric(NA)
# ###             cat("Approximate logLik:",
# ###                 .Call("nlme_glmmLa2_logLikelihood",
# ###                     fit, NULL, PACKAGE="lme4"), "\n")
# ###            class(fit) <- "reStruct"
# ###            fit@logLik <- as.numeric(NA)
#                   cat("Parameters:", coef(fit), "\n")
#                   cat("Fixed Effects:", fixef(fit), "\n")
#               }
#               etaold <- eta
#               eta <- .Call("nlme_reStruct_fitted", fit, NULL, PACKAGE="lme4")
#               if(sum((eta-etaold)^2) < 1e-6*sum(eta^2)) {
#                   converged <- TRUE
#                   break
#               }
#           }
#           if (control$msMaxIter > 0 && !converged)
#               stop("IRLS iterations in sparseGLMM failed to converge")
#           ## We recreate the glmm object and set the reStruct slot
#           .Call("nlme_replaceSlot", fit, "logLik",
#                 as.numeric(NA), PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, "dontCopy", TRUE, PACKAGE = "lme4")
#           fit <- .Call("nlme_replaceSlot", eval(m, parent.frame()),
#                        "reStruct", fit, PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, c("reStruct", "dontCopy"),
#                 TRUE, PACKAGE = "lme4")
#           fit <- .Call("nlme_glmmLaplace_solveOnly", fit,
#                        500, 1, PACKAGE="lme4")
#           ## dirtyStored is FALSE but should be TRUE
# #          .Call("nlme_replaceSlot", fit, c("reStruct", "dirtyStored"),
# #                TRUE, PACKAGE = "lme4")
# #          .Call("nlme_replaceSlot", fit, "reStruct",
# #                .Call("nlme_commonDecompose", fit@reStruct,
# #                      NULL, PACKAGE="lme4"), PACKAGE = "lme4")
# #          .Call("nlme_replaceSlot", fit, c("reStruct", "dontCopy"),
# #                FALSE, PACKAGE = "lme4")

#           if (method != "PQL") {
#               ## Do the 2nd order Laplace fit here
#               LMEoptimize(fit) <- control
#           }
#           ## dirtyStored is FALSE but should be TRUE
#           .Call("nlme_replaceSlot", fit, c("reStruct", "dirtyStored"),
#                 TRUE, PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, "reStruct",
#                 .Call("nlme_commonDecompose", fit@reStruct,
#                       NULL, PACKAGE="lme4"), PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, c("reStruct", "dontCopy"),
#                 FALSE, PACKAGE = "lme4")
#           ## zero some of the matrix slots
#           if (!missing(x) && x == FALSE)
#               .Call("nlme_replaceSlot", fit, c("reStruct", "original"),
#                     matrix(0.0, nrow = 0, ncol = 0), PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, c("reStruct", "decomposed"),
#                 matrix(0.0, nrow = 0, ncol = 0), PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, c("reStruct", "weighted"),
#                 matrix(0.0, nrow = 0, ncol = 0), PACKAGE = "lme4")

#           if (!missing(model) && model == FALSE)
#               .Call("nlme_replaceSlot", fit, "frame",
#                     data.frame(), PACKAGE = "lme4")
#           .Call("nlme_replaceSlot", fit, "fitted",
#                 fam$linkinv(if (is.null(fit@na.action)) {
#                     fitted(fit@reStruct)[fit@reStruct@reverseOrder]
#                 } else {
#                     napredict(attr(data, "na.action"),
#                               fitted(fit@reStruct)[fit@reStruct@reverseOrder])
#                 }), PACKAGE = "lme4")
#           fit
#       })


### Local variables:
### mode: R
### End:
