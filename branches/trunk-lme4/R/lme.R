facshuffle = function(sslm, facs)       # unexported utility
{
    if (!length(sslm[[2]])) return(facs)
    s1 = sslm[[1]]
    s2 = sslm[[2]]
    lens = diff(s1@Gp)
    ff = vector("list", length(facs))
    for (i in seq(along = lens)) {
        sq = seq(lens[i])
        perm = 1 + s2[sq]
        s2 = s2[-sq] - lens[i]
        fi = facs[[i]]
        fip = factor(perm[as.integer(fi)])
        levels(fip)[perm] = levels(fi)
        ff[[i]] = fip
    }
    ff
}

lmeControl <-
  ## Control parameters for lme
  function(maxIter = 50, msMaxIter = 50, tolerance =
           sqrt((.Machine$double.eps)), niterEM = 25,
           msTol = sqrt(.Machine$double.eps), msScale, msVerbose = FALSE,
           glmmMaxIter = 20,
           returnObject = FALSE, gradHess = TRUE, apVar = TRUE,
           .relStep = (.Machine$double.eps)^(1/3), minAbsParApVar = 0.05,
           nlmStepMax = NULL,
           natural = TRUE, optimizer="nlm", EMverbose=FALSE,
           analyticGradient = FALSE,
           analyticHessian=FALSE)
{
    if (missing(msScale)) msScale = function(start) {
        scale <- abs(start)
        nonzero <- scale > 0
        if (any(nonzero)) {
            scale[nonzero] <- 1/scale[nonzero]
            scale[!nonzero] <- median(scale[nonzero])
        }
        else {
            scale <- rep(1, length(scale))
        }
        scale
    }
    list(maxIter = maxIter, msMaxIter = msMaxIter, tolerance = tolerance,
         niterEM = niterEM, msTol = msTol, msScale = msScale,
         msVerbose = msVerbose,
         glmmMaxIter = glmmMaxIter,
         returnObject = returnObject,
         gradHess = gradHess , apVar = apVar, .relStep = .relStep,
         nlmStepMax = nlmStepMax,
         minAbsParApVar = minAbsParApVar, natural = natural,
         optimizer=optimizer, EMverbose=EMverbose,
         analyticHessian=analyticHessian,
         analyticGradient=analyticGradient)
}

setMethod("lme", signature(data = "missing"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          nCall$data = list()
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("lme", signature(formula = "missing", data = "groupedData"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          resp = getResponseFormula(data)[[2]]
          cov = getCovariateFormula(data)[[2]]
          nCall$formula = eval(substitute(resp ~ cov))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("lme", signature(formula = "formula", data = "groupedData",
                           random = "missing"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          cov = formula[[3]]
          grps = getGroupsFormula(data)[[2]]
          nCall$random = eval(substitute(~ cov | grps))
          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })


setMethod("lme", signature(formula = "formula", random = "formula"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
          nCall = mCall = match.call()
          nCall$random = lapply(getGroupsFormula(random, asList = TRUE),
                                function(x, form) form,
                                form = pdLogChol(getCovariateFormula(random)))

          nCall$data <- as(data, "data.frame")

          .Call("nlme_replaceSlot", eval(nCall, parent.frame()), "call",
                mCall, PACKAGE = "lme4")
      })

setMethod("lme", signature(formula = "formula", random = "list"),
          function(formula, data, random, correlation, weights, subset,
                   method, na.action, control, model, x)
      {
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
          controlvals$REML = method == "REML"
          mCall <- match.call(expand.dots = FALSE)
          mCall[[1]] <- as.name("model.frame")
          names(mCall)[2] <- "formula"
          mCall$random <- mCall$correlation <- mCall$method <-
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
          environment(form) = environment(formula)
          mCall$formula = form
          mCall$drop.unused.levels = TRUE
          data = eval(mCall, parent.frame())
          facs = lapply(names(random),
                         function(x) eval(as.name(x), envir = data))
          names(facs) = names(random)
          mmats <- c(lapply(random,
                            function(x) model.matrix(formula(x), data = data)),
                     list(.Xy = cbind(model.matrix(formula, data = data),
                          .response = model.response(data))))
          obj = .Call("ssclme_create", facs, unlist(lapply(mmats, ncol)),
                       as.integer(2e5), PACKAGE = "Matrix")
          facs = facshuffle(obj, facs)
          obj = obj[[1]]
          .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
          .Call("ssclme_initial", obj, PACKAGE="Matrix")
          .Call("ssclme_EMsteps", obj, controlvals$niterEM,
                controlvals$REML, controlvals$EMverbose, PACKAGE = "Matrix")
          LMEoptimize(obj) = controlvals
          fitted = .Call("ssclme_fitted", obj, facs, mmats, PACKAGE = "Matrix")
          new("lme", call = match.call(), facs = facs,
              x = if(x) mmats else list(),
              model = if(model) data else data.frame(list()),
              REML = method == "REML", rep = obj, fitted = fitted)
      })

setMethod("fitted", signature=c(object="lme"),
          function(object, ...)
      {
          object@fitted
      })

setMethod("residuals", signature=c(object="lme"),
          function(object, ...) NULL)

setMethod("logLik", signature(object="lme"),
          function(object, REML = FALSE, ...)
          -deviance(object@rep, REML = REML)/2)

setMethod("deviance", signature(object="lme"),
          function(object, REML, ...)
          deviance(object@rep,
                   REML = ifelse(missing(REML), object@REML, REML))
          )

#setMethod("summary", signature(object="lme"),
#          function(object, ...) {
#              llik <- logLik(object)    # has an oldClass
#              resd <- residuals(object, type="pearson")
#              if (length(resd) > 5) {
#                  resd <- quantile(resd)
#                  names(resd) <- c("Min","Q1","Med","Q3","Max")
#              }
#              new("summary.lme",
#                  call = object@call,
#                  logLik = llik,
#                  AIC = AIC(llik),
#                  BIC = BIC(llik),
#                  re = summary(as(object, "reStruct")),
#                  residuals = resd)
#          })

#setMethod("show", "summary.lme",
#          function(object) {
#              rdig <- 5
#              cat("Linear mixed-effects model fit by ")
#              cat(ifelse(object@re@REML, "REML\n", "maximum likelihood\n") )
#              cat(" Data:", deparse( object@call$data ), "\n")
#              if (!is.null(object@call$subset)) {
#                  cat("  Subset:",
#                      deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
#              }
#              print(data.frame(AIC = object@AIC, BIC = object@BIC,
#                               logLik = c(object@logLik), row.names = ""))
#              cat("\n")
#              object@re@useScale = TRUE
#              object@re@showCorrelation = TRUE
#              show(object@re)
#              ## Should this be part of the show method for summary.reStruct?
#              cat("\nNumber of Observations:", object@re@nobs)
#              cat("\nNumber of Groups: ")
#              ngrps <- object@re@ngrps
#              if ((length(ngrps)) == 1) {
#                  cat(ngrps,"\n")
#              } else {				# multiple nesting
#                  cat("\n")
#                  print(ngrps)
#              }
#              invisible(object)
#          })

setMethod("show", "lme",
          function(object)
      {
          #sumry = summary(object)
          rdig <- 5
          cat("Linear mixed-effects model\n")
          cat(" Data:", deparse( object@call$data ), "\n")
          if (!is.null(object@call$subset)) {
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
          }
          cat(paste(" log-", ifelse(object@REML, "restricted-", ""),
                    "likelihood: ", sep = ''), logLik(object), "\n")
          show(fixef(object))
          show(c(object@rep@Omega,
                 sigmaSq = .Call("ssclme_sigma",
                                 object@rep, PACKAGE="Matrix")^2))
          nc = object@rep@nc
          cat("\nNumber of Observations:", nc[length(nc)], "\n")
          invisible(object)
      })

setMethod("anova", signature(object = "lme"),
          function(object, ...)
          cat("anova method for lme not yet implemented\n"))

setMethod("fixef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@rep
          callGeneric()
      })

setMethod("formula", "lme", function(x, ...) x@call$formula)

setMethod("intervals", signature(object = "lme", level = "ANY"),
          function(object, level = 0.95, ...)
          cat("intervals method for lme not yet implemented\n"))

setMethod("plot", signature(x = "lme"),
          function(x, y, ...)
          cat("plot method for lme not yet implemented\n"))

setMethod("ranef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@rep
          callGeneric()
      })

setMethod("coef", signature(object = "lme"),
          function(object, ...)
      {
          object = object@rep
          callGeneric()
      })

setMethod("update", signature(object = "lme"),
          function(object, formula., ..., evaluate = TRUE)
      {
          call <- object@call
          if (is.null(call))
              stop("need an object with call component")
          extras <- match.call(expand.dots = FALSE)$...
          if (!missing(formula.))
              call$formula <- update.formula(formula(object), formula.)
          if (length(extras) > 0) {
              existing <- !is.na(match(names(extras), names(call)))
              for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
              if (any(!existing)) {
                  call <- c(as.list(call), extras[!existing])
                  call <- as.call(call)
              }
          }
          if (evaluate)
              eval(call, parent.frame())
          else call
      })

setMethod("vcov", signature(object = "lme"),
          function(object, ...) {
              object = object@rep
              callGeneric()
          })

setMethod("VarCorr", signature(x = "lme"),
          function(x) {
              x = x@rep
              callGeneric()
          })
