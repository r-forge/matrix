


testfun <- function(...) # test function
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
                   control = lmeControl(EMverbose = F), ...)


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



          debug <- TRUE ## check if fitted() works. Remove all such code later




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
          controlvals$REML <- method == "REML"



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
          .Call("ssclme_initial", obj, PACKAGE="Matrix")




          ## get initial estimates
          fm.glm <- glm(formula, family, data)
          coefFixed <- c(coef(fm.glm), 0)


          if (debug)
          {
              coefRanef <-
                  lapply(facs,
                         function(x) numeric(nlevels(x)))
          }


          mmats <- mmats.unadjusted
          niter <- 20
          conv <- FALSE
          firstIter <- TRUE
          devold <- 0 ## deviance for convergence check

          ## initial 'fitted' values on linear scale
          eta <- drop(mmats.unadjusted$.Xy %*% coefFixed)
          etaold <- eta + 1


          for (iter in seq(length = niter))
          {
              if (debug)
              {
                  eta.check <- mmats.unadjusted$.Xy %*% coefFixed
                  for (facname in names(facs))
                  {
                      eta.check <- eta.check + mmats.unadjusted[[facname]] * 
                          coefRanef[[facname]][facs[[facname]]]
                  }

                  if (any(eta.check != eta)) {
                      warning("fitted() does not match calculation, diff: ",
                              sum(((eta.check - eta)^2)))
                  }
              }


              mu <- family$linkinv(eta)
              dmu.deta <- family$mu.eta(eta)
              ## adjusted response
              z <- eta + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta
              ## weights
              w <- dmu.deta / sqrt(family$variance(mu))

plot(z, mmats$.Xy[, responseIndex])

              ## Does this prevent overwriting of components ?
              for (facname in names(facs))
                  mmats[[facname]][] <- mmats.unadjusted[[facname]] * w
              mmats$.Xy[] <- mmats.unadjusted$.Xy
              mmats$.Xy[, responseIndex] <- z
              mmats$.Xy[] <- mmats$.Xy * w

              .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
              ## ssclme_initial should only be called on the first iteration
              if (firstIter) .Call("ssclme_initial", obj, PACKAGE="Matrix")
              .Call("ssclme_EMsteps", obj, controlvals$niterEM,
                    method == "REML", controlvals$EMverbose, PACKAGE = "Matrix")
              LMEoptimize(obj) = controlvals
              eta[] <- .Call("ssclme_fitted", obj, facs, mmats.unadjusted, PACKAGE = "Matrix")

print(sqrt(max((eta - etaold)^2)) /
                  (0.1 + sqrt(max(eta^2))))

              ## use this to determine convergence
              if (sqrt(max((eta - etaold)^2)) /
                  (0.1 + sqrt(max(eta^2))) <
                  controlvals$tolerance) {
                  conv <- TRUE
                  break
              }
              etaold[] <- eta





              if (debug)
              {
                  coefFixed <- c(.Call("ssclme_fixef", obj, PACKAGE = "Matrix"), 0)
                  coefRanef <- .Call("ssclme_ranef", obj, PACKAGE = "Matrix")
                  if (!is.null(names(coefRanef))) print("time to modify code here")
                  names(coefRanef) <- names(facs)
              }



### May want to adjust number of EMiterations on second and subsequent
### iterations.  Probably one or two iterations will suffice

### Do you want to call LMEoptimize(obj) = controlvals here?  Probably
### set number of iterations low.

              if (firstIter) {
                  controlvals$niterEM <- 2
                  controlvals$msMaxIter <- 10
                  firstIter <- FALSE
              }

          }
          if (!conv) warning("IRLS iterations for glmm did not converge")

          new("lme", call = match.call(), facs = facs,
              x = if(x) mmats else list(),
              model = if(model) data else data.frame(list()),
              REML = method == "REML", rep = obj, fitted = eta)
      })











## what happens in glm.fit:


#             eta <- drop(x %*% start)
#             mu <- linkinv(eta <- eta + offset)
#             dev <- sum(dev.resids(y, mu, weights))
#             if (control$trace)
#                 cat("Deviance =", dev, "Iterations -", iter, "\n")



#             ## check for divergence
#             boundary <- FALSE
#             if (!is.finite(dev)) {
#                 if(is.null(coefold))
#                     stop("no valid set of coefficients has been found:please supply starting values", call. = FALSE)
#                 warning("Step size truncated due to divergence", call. = FALSE)
#                 ii <- 1
#                 while (!is.finite(dev)) {
#                     if (ii > control$maxit)
#                         stop("inner loop 1; can't correct step size")
#                     ii <- ii + 1
#                     start <- (start + coefold)/2
#                     eta <- drop(x %*% start)
#                     mu <- linkinv(eta <- eta + offset)
#                     dev <- sum(dev.resids(y, mu, weights))
#                 }
#                 boundary <- TRUE
#                 if (control$trace)
#                     cat("Step halved: new deviance =", dev, "\n")
#             }



#             ## check for fitted values outside domain.
#             if (!(valideta(eta) && validmu(mu))) {
#                 warning("Step size truncated: out of bounds", call. = FALSE)
#                 ii <- 1
#                 while (!(valideta(eta) && validmu(mu))) {
#                     if (ii > control$maxit)
#                         stop("inner loop 2; can't correct step size")
#                     ii <- ii + 1
#                     start <- (start + coefold)/2
#                     eta <- drop(x %*% start)
#                     mu <- linkinv(eta <- eta + offset)
#                 }
#                 boundary <- TRUE
#                 dev <- sum(dev.resids(y, mu, weights))
#                 if (control$trace)
#                     cat("Step halved: new deviance =", dev, "\n")
#             }




#             ## check for convergence
#             if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
#                 conv <- TRUE
#                 coef <- start
#                 break
#             } else {
#                 devold <- dev
#                 coef <- coefold <- start
#             }
#         } ##-------------- end IRLS iteration -------------------------------

#         if (!conv) warning("Algorithm did not converge")
#         if (boundary) warning("Algorithm stopped at boundary value")
#         eps <- 10*.Machine$double.eps
#         if (family$family == "binomial") {
#             if (any(mu > 1 - eps) || any(mu < eps))
#                 warning("fitted probabilities numerically 0 or 1 occurred")
#         }
#         if (family$family == "poisson") {
#             if (any(mu < eps))
#                 warning("fitted rates numerically 0 occurred")
#         }
#         ## If X matrix was not full rank then columns were pivoted,
#         ## hence we need to re-label the names ...
#         ## Original code changed as suggested by BDR---give NA rather
#         ## than 0 for non-estimable parameters
#         if (fit$rank < nvars) coef[fit$pivot][seq(fit$rank+1, nvars)] <- NA
#         xxnames <- xnames[fit$pivot]
#         residuals <- rep.int(NA, nobs)
#         residuals[good] <- z - (eta - offset)[good] # z does not have offset in.
#         fit$qr <- as.matrix(fit$qr)
#         nr <- min(sum(good), nvars)
#         if (nr < nvars) {
#             Rmat <- diag(nvars)
#             Rmat[1:nr, 1:nvars] <- fit$qr[1:nr, 1:nvars]
#         }
#         else Rmat <- fit$qr[1:nvars, 1:nvars]
#         Rmat <- as.matrix(Rmat)
#         Rmat[row(Rmat) > col(Rmat)] <- 0
#         names(coef) <- xnames
#         colnames(fit$qr) <- xxnames
#         dimnames(Rmat) <- list(xxnames, xxnames)


















### Local variables:
### mode: R
### End:
