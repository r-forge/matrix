


testfun <- function(...) # test function
{
    ## simulation: 300 obs, 30 students, 10 obs per student, one covariate

    dat <- data.frame(id1 = gl(30, 10),
                      id2 = gl(60, 5),
                      x = rnorm(300))

    dat$resp <-
        with(dat,
             rbinom(300,
                    size = 1,
                    prob = binomial()$linkinv( x + rnorm(30, sd = .5)[id1])))


#    fm.GLMM <- GLMM(formula = resp ~ x, data = dat, family = binomial(), random = list(id = ~1))


    fm.sparse1 <-
        sparseGLMM(resp ~ x, data = dat, family = binomial(),
                   random = list(id1 = ~1, id2 = ~1),
                   control = lmeControl(...))
    fm.sparse2 <-
        sparseGLMM(resp ~ x, data = dat, family = binomial(),
                   random = list(id2 = ~1, id1 = ~1),
                   control = lmeControl(...))
    list(fm.sparse1, fm.sparse2)
}



testfun2 <- function(...)
{
    data(guImmun)
    guImmun$numimmun <- as.numeric(guImmun$immun) - 1
    fm1.gui <-
        sparseGLMM(numimmun ~ kid2p + mom25p + ord + ethn +
                   momEd + husEd + momWork + rural + pcInd81,
                   data = guImmun, family = binomial(),
                   random = list(mom = ~1, comm = ~1),
                   control = lmeControl(...))
    fm2.gui <-
        sparseGLMM(numimmun ~ kid2p + mom25p + ord + ethn +
                   momEd + husEd + momWork + rural + pcInd81,
                   data = guImmun, family = binomial(),
                   random = list(comm = ~1, mom = ~1),
                   control = lmeControl(...))
    list(fm1 = fm1.gui, fm2 = fm2.gui)
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

setMethod("sparseGLMM",
          signature(formula = "formula",
                    family = "family",
                    random = "list"),
          function(formula, family, data, random, subset,
                   method, na.action, control, model, x, ...)
      {
          if (missing(method)) method <- "PQL" # FIXME: see glmm.R 

          debug <- FALSE ## check if fitted() works. Remove all such code later



          if (missing(model))
              model = TRUE
          if (missing(x))
              x = FALSE
          random = lapply(as(random, "list"),
                   get("formula", pos = parent.frame(), mode = "function"))
          controlvals <- if (missing(control)) lmeControl() else
                            do.call("lmeControl", control)
          controlvals$REML <- FALSE
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
          facs =                        # order in decreasing number of levels
              facs[rev(order(unlist(lapply(facs,
                                           function(fac)
                                           length(levels(fac))))))]





          ## creates model matrices
          mmats.unadjusted <-
              c(lapply(random,
                       function(x) model.matrix(formula(x), data = data)),
                list(.Xy =
                     cbind(model.matrix(formula, data = data),
                           .response = model.response(data))))
          responseIndex <- ncol(mmats.unadjusted$.Xy)


          ## creates ssclme structure
          obj <- .Call("ssclme_create", facs,
                       unlist(lapply(mmats.unadjusted, ncol)),
                       as.integer(2e5), PACKAGE = "Matrix")


          ## FIXME: names of facs lost, but may be useful later
          facs = facshuffle(obj, facs)
          print(names(facs)) ## = names(random)

          obj = obj[[1]]


          ## get initial estimates from glm
          coefFixef <- c(coef(glm(formula, family, data)), 0)


#           ## debugging code, remove later *******************
#           if (debug)
#           {
#               coefRanef <-
#                   lapply(facs,
#                          function(x) numeric(nlevels(x)))
#               for (facname in names(facs))
#               {
#                   coefRanef[[facname]] <-
#                       matrix(rep(coefRanef[[facname]], ncol(mmats.unadjusted[[facname]])),
#                              ncol = ncol(mmats.unadjusted[[facname]]))
#               }
#           }
#           ## **********************************************



          mmats <- mmats.unadjusted
          conv <- FALSE
          firstIter <- TRUE

          ## initial 'fitted' values on linear scale
          eta <- drop(mmats.unadjusted$.Xy %*% coefFixef)
          etaold <- eta

          for (iter in seq(length = controlvals$glmmMaxIter))
          {

              cat(paste("Iteration", iter, "\n"))


#               ## debugging code, remove later *******************
#               if (debug)
#               {
#                   eta.check <- mmats.unadjusted$.Xy %*% coefFixef
#                   for (facname in names(facs))
#                   {
#                       eta.check <- eta.check +
#                           mmats.unadjusted[[facname]] * coefRanef[[facname]][facs[[facname]],]
#                   }
#                   eta.check <- drop(eta.check)


#                   if (!is.logical(all.equal(eta.check, eta))) {
#                       warning("fitted() does not match calculation, diff: ",
#                               sum(((eta.check - eta)^2)))
#                   }
#                   else print("fitted values match")
#                   #eta <- drop(eta.check)
#               }
#               ## *************************************************





              mu <- family$linkinv(eta)
              dmu.deta <- family$mu.eta(eta)



              ## adjusted response
              z <- eta + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta
              ## weights
              w <- dmu.deta / sqrt(family$variance(mu))

              ## Does this prevent overwriting of components ?
              for (facname in names(facs))
                  mmats[[facname]][] <- mmats.unadjusted[[facname]] * w
              mmats$.Xy[] <- mmats.unadjusted$.Xy
              mmats$.Xy[, responseIndex] <- z
              mmats$.Xy[] <- mmats$.Xy * w







              .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
              if (firstIter) .Call("ssclme_initial", obj, PACKAGE="Matrix")
print(str(obj@Omega))
              .Call("ssclme_EMsteps", obj,
                    controlvals$niterEM,
                    FALSE, #controlvals$REML,
                    controlvals$EMverbose,
                    PACKAGE = "Matrix")
              LMEoptimize(obj) = controlvals
              eta[] <-
                  .Call("ssclme_fitted", obj, facs,
                        mmats.unadjusted, PACKAGE = "Matrix")


              cat(paste("\tTermination Criteria:", max(abs(eta - etaold)) /
                        (0.1 + max(abs(eta))), "\n"))

              ## use this to determine convergence
              if (max(abs(eta - etaold)) <
                  (0.1 + max(abs(eta))) * controlvals$tolerance)
              {
                  conv <- TRUE
                  break
              }
              etaold[] <- eta





#               ## debugging code, remove later *******************
#               if (debug)
#               {
#                   coefFixef <- c(fixef(obj), 0)
#                   coefRanef <- ranef(obj)
#               }
#               ## *********************************************


              ## Changing number of iterations on second and
              ## subsequent iterations.

              if (firstIter)
              {
                  controlvals$niterEM <- 2 # FIXME, should be 2
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

          offset <- mmats.unadjusted$.Xy %*% c(fixef(obj), 0)
          reducedObj <- .Call("ssclme_collapse", obj, PACKAGE = "Matrix")
          reducedMmats.unadjusted <- mmats.unadjusted
          reducedMmats.unadjusted$.Xy <-
              reducedMmats.unadjusted$.Xy[, responseIndex, drop = FALSE]
          reducedMmats <- reducedMmats.unadjusted

          laplik <-
              LaplaceLikelihood(obj = reducedObj,
                                facs = facs, offset = offset,
                                mm.unadjusted = reducedMmats.unadjusted,
                                mm = reducedMmats,
                                family = family)




          new("lme", call = match.call(), facs = facs,
              x = if (x) mmats else list(),
              model = if(model) data else data.frame(list()),
              REML = FALSE, rep = obj, fitted = eta)
      })





LaplaceLikelihood <-
    function(obj, facs, offset, mmats.unadjusted, mmats, family)
{

    niter <- 20
    conv <- FALSE

    ## will this work ?
    eta <-
        .Call("ssclme_fitted", obj, facs,
              mmats.unadjusted, PACKAGE = "Matrix")
    etaold <- eta


    for (iter in seq(length = niter))
    {

        mu <- family$linkinv(eta)
        dmu.deta <- family$mu.eta(eta)

        ## adjusted response (how does offset get involved ?)
        z <- eta + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta
        ## weights
        w <- dmu.deta / sqrt(family$variance(mu))

        ## Does this prevent overwriting of components ?
        for (facname in names(facs))
            mmats[[facname]][] <- mmats.unadjusted[[facname]] * w
        mmats$.Xy[] <- mmats.unadjusted$.Xy
        mmats$.Xy[, responseIndex] <- z
        mmats$.Xy[] <- mmats$.Xy * w


        ## the following needs to be changed, not sure how

        .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")



        if (firstIter) .Call("ssclme_initial", obj, PACKAGE="Matrix")
        etaold <- eta
        eta[] <-
            .Call("ssclme_fitted", obj, facs,
                  mmats.unadjusted, PACKAGE = "Matrix") + offset


        if (FALSE) ## converged
        {
            conv <- TRUE
            break
        }


    }


}







### Local variables:
### mode: R
### End:
