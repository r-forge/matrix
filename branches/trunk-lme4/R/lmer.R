contr.SAS <- function(n, contrasts = TRUE)
{
    if (is.numeric(n) && length(n) == 1) contr.treatment(n, n, contrasts)
    else contr.treatment(n, length(n), contrasts)
}

lmerControl <-                            # Control parameters for lmer
  function(maxIter = 50,
           msMaxIter = 50,
           tolerance = sqrt((.Machine$double.eps)),
           niterEM = 20,
           msTol = sqrt(.Machine$double.eps),
           msVerbose = getOption("verbose"),
           PQLmaxIt = 20,
           .relStep = (.Machine$double.eps)^(1/3),
           EMverbose = getOption("verbose"),
           analyticGradient = TRUE,
           analyticHessian=FALSE)
{
    list(maxIter = maxIter,
         msMaxIter = msMaxIter,
         tolerance = tolerance,
         niterEM = niterEM,
         msTol = msTol,
         msVerbose = msVerbose,
         PQLmaxIt = PQLmaxIt,
         .relStep = .relStep,
         EMverbose=EMverbose,
         analyticHessian=analyticHessian,
         analyticGradient=analyticGradient)
}

setMethod("lmer", signature(formula = "formula"),
          function(formula, data,
                   method = c("REML", "ML"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
                                        # match and check parameters
          method <- match.arg(method)
          controlvals <- do.call("lmerControl", control)
          controlvals$REML <- method == "REML"
          if (length(formula) < 3) stop("formula must be a two-sided formula")
                                        # create the model frame as frm
          mf <- match.call()
          m <- match(c("data", "subset", "weights", "na.action", "offset"),
                     names(mf), 0)
          mf <- mf[c(1, m)]
          mf[[1]] <- as.name("model.frame")
          frame.form <- subbars(formula)
          environment(frame.form) <- environment(formula)
          mf$formula <- frame.form
          mf$drop.unused.levels <- TRUE
          frm <- eval(mf, parent.frame())

          ## grouping factors and model matrices for random effects
          bars <- findbars(formula[[3]])
          random <-
              lapply(bars,
                     function(x) list(model.matrix(eval(substitute(~term,
                                                                   list(term=x[[2]]))),
                                                   frm),
                                      eval(substitute(as.factor(fac),
                                                      list(fac = x[[3]])), frm)))
          names(random) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

          ## order factor list by decreasing number of levels
          ford <- rev(order(sapply(random, function(x) length(levels(x[[2]])))))
          if (any(ford != seq(a = random))) { # re-order both facs and random
              random <- random[ford]
          }
          mmats <- c(lapply(random, "[[", 1),
                     .fixed = list(cbind(model.matrix(nobars(formula), frm),
                     .response = model.response(frm))))
          obj <- .Call("lmer_create", lapply(random, "[[", 2), mmats, PACKAGE = "Matrix")
          obj@call <- match.call()
          obj@REML <- method == "REML"
          .Call("lmer_initial", obj, PACKAGE="Matrix")
          .Call("lmer_ECMEsteps", obj, 
                controlvals$niterEM,
                controlvals$REML,
                controlvals$EMverbose,
                PACKAGE = "Matrix")
          LMEoptimize(obj) <- controlvals
          #fitted = .Call("ssclme_fitted", obj, facs, mmats, TRUE, PACKAGE = "Matrix")
          #residuals = mmats$.Xy[,".response"] - fitted
          #if (as.logical(x)[1]) x = mmats else x = list()
          #rm(mmats)
          obj
      })

setReplaceMethod("LMEoptimize", signature(x="lmer", value="list"),
                 function(x, value)
             {
                 if (value$msMaxIter < 1) return(x)
                 st <- ccoef(x)         # starting values
                 nc <- x@nc
                 nc <- nc[1:(length(nc) - 2)]
                 constr <- unlist(lapply(nc, function(k) 1:((k*(k+1))/2) <= k))
                 fn <- function(pars) {
                     ccoef(x) <- pars
                     deviance(x, REML = value$REML)
                 }
                 gr <- if (value$analyticGradient)
                     function(pars) {
                         ccoef(x) <- pars
                         grad <- lme4:::gradient(x, REML = value$REML, unconst = TRUE)
                         grad[constr] <- -grad[constr]/pars[constr]
                         grad
                     } else NULL
                 optimRes <- optim(st, fn, gr,
                                   method = "L-BFGS-B",
                                   lower = ifelse(constr, 1e-10, -Inf),
                                   control = list(maxit = value$msMaxIter,
                                   trace = as.integer(value$msVerbose)))
                 if (optimRes$convergence != 0) {
                     warning(paste("optim returned message",optimRes$message,"\n"))
                 }
                 ccoef(x) <- optimRes$par
                 return(x)
             })

setMethod("deviance", signature(object = "lmer"),
          function(object, REML = FALSE, ...) {
              .Call("lmer_factor", object, PACKAGE = "Matrix")
              object@deviance[ifelse(REML, 2, 1)]
          })

setMethod("ranef", signature(object = "lmer"),
          function(object, ...) {
              .Call("lmer_ranef", object, PACKAGE = "Matrix")
          })

setMethod("fixef", signature(object = "lmer"),
          function(object, ...) {
              val = .Call("lmer_fixef", object, PACKAGE = "Matrix")
              names(val) = object@cnames[[".fixed"]]
              val[-length(val)]
          })

setMethod("vcov", signature(object = "lmer"),
          function(object, REML = TRUE, useScale = TRUE,...) {
              ## force an "lmer_invert"
              sc <- .Call("lmer_sigma", object, REML, PACKAGE = "Matrix")
              rr <- object@RXX
              nms <- object@cnames[[".fixed"]]
              dimnames(rr) <- list(nms, nms)
              nr <- nrow(rr)
              rr <- rr[-nr, -nr, drop = FALSE]
              rr <- rr %*% t(rr)
              if (useScale) {
                  rr = sc^2 * rr
              }
              rr
          })

setMethod("VarCorr", signature(x = "lmer"),
          function(x, REML = TRUE, useScale = TRUE, ...) {
              val = .Call("lmer_variances", x, PACKAGE = "Matrix")
              for (i in seq(along = val)) {
                  dimnames(val[[i]]) = list(x@cnames[[i]], x@cnames[[i]])
                  val[[i]] = as(as(val[[i]], "pdmatrix"), "corrmatrix")
              }
              new("VarCorr",
                  scale = .Call("lmer_sigma", x, REML),
                  reSumry = val,
                  useScale = useScale)
          })

setMethod("gradient", signature(x = "lmer"),
          function(x, REML, unconst, ...)
          .Call("lmer_gradient", x, REML, unconst))

setMethod("summary", "lmer",
          function(object, REML = TRUE, useScale = TRUE, ...) {
              fcoef <- fixef(object)
              corF <- as(as(vcov(object, REML, useScale), "pdmatrix"),
                         "corrmatrix")
              DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@stdDev, DF)
              nc <- object@nc
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
              new("summary.lmer",
                  coefficients = as.matrix(coefs),
                  scale = .Call("lmer_sigma", object, REML),
                  denomDF = as.integer(DF),
                  REML = REML,
                  ngrps = unlist(lapply(object@flist,
                                        function(x) length(levels(x)))),
                  nobs = nc[length(nc)],
                  corFixed = corF,
                  VarCorr = VarCorr(object, REML, useScale),
                  useScale = useScale,
                  showCorrelation = FALSE)
          })

setMethod("show", "lmer",
          function(object) {
              fcoef <- fixef(object)
              corF <- as(as(vcov(object, REML = TRUE, useScale = TRUE), "pdmatrix"),
                         "corrmatrix")
              DF <- getFixDF(object)
              coefs <- cbind(fcoef, corF@stdDev, DF)
              nc <- object@nc
              dimnames(coefs) <-
                  list(names(fcoef), c("Estimate", "Std. Error", "DF"))
              new("summary.ssclme",
                  coefficients = as.matrix(coefs),
                  scale = .Call("lmer_sigma", object, REML = TRUE),
                  denomDF = as.integer(DF),
                  REML = TRUE,
                  ngrps = unlist(lapply(object@flist,
                                        function(x) length(levels(x)))),
                  nobs = nc[length(nc)],
                  corFixed = corF,
                  VarCorr = VarCorr(object, REML = TRUE, useScale = TRUE),
                  useScale = TRUE,
                  showCorrelation = FALSE)
          })

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

setMethod("getFixDF", signature(object="lmer"),
          function(object, ...)
      {
          nc <- object@nc[-seq(along = object@Omega)]
          p <- nc[1] - 1
          n <- nc[2]
          rep(n-p, p)
      })

setMethod("fitted", signature(object="lmer"),
          function(object, ...)
      {
          object@fitted
      })

setMethod("residuals", signature(object="lmer"),
          function(object, ...) object@residuals )

setMethod("logLik", signature(object="lmer"),
          function(object, REML = object@REML, ...) {
              val <- -deviance(object, REML = REML)/2
              nc <- object@nc[-seq(a = object@Omega)]
              attr(val, "nall") <- attr(val, "nobs") <- nc[2]
              attr(val, "df") <- nc[1] + length(coef(object))
              attr(val, "REML") <- REML 
              class(val) <- "logLik"
              val
          })

# setMethod("summary", signature(object="lme"),
#           function(object, ...) {
#               llik <- logLik(object)
#               resd <- residuals(object, type="pearson")
#               if (length(resd) > 5) {
#                   resd <- quantile(resd)
#                   names(resd) <- c("Min","Q1","Med","Q3","Max")
#               }
#               new("summary.lme",
#                   call = object@call,
#                   logLik = llik,
#                   re = summary(object@rep, REML = object@REML,
#                                useScale = TRUE),
#                   residuals = resd)
#           })

# setMethod("show", signature(object = "summary.lme"),
#           function(object)
#       {
#           rdig <- 5
#           cat("Linear mixed-effects model fit by ")
#           cat(ifelse(object@re@REML, "REML\n", "maximum likelihood\n") )
#           if (!is.null(object@call$formula)) {
#               cat("Fixed:", deparse(object@call$formula),"\n")
#           }
#           if (!is.null(object@call$data)) {
#               cat(" Data:", deparse(object@call$data), "\n")
#           }
#           if (!is.null(object@call$subset)) {
#               cat(" Subset:",
#                   deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")
#           }
#           llik <- object@logLik
#           print(data.frame(AIC = AIC(llik), BIC = BIC(llik),
#                            logLik = c(object@logLik), row.names = ""))
#           cat("\n")
#           object@re@useScale <- TRUE
#           object@re@showCorrelation <- TRUE
#           show(object@re)
#           invisible(object)
#       })

setMethod("anova", signature(object = "lmer"),
          function(object, ...)
      {
          mCall <- match.call(expand.dots = TRUE)
          dots <- list(...)
          modp <- logical(0)
          if (length(dots))
              modp <- sapply(dots, inherits, "lmer") | sapply(dots, inherits, "lm")
          if (any(modp)) {              # multiple models - form table
              opts <- dots[!modp]
              mods <- c(list(object), dots[modp])
              names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)], as.character)
              mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE), attr, "df"))]
              calls <- lapply(mods, slot, "call")
              data <- lapply(calls, "[[", "data")
              if (any(data != data[[1]])) stop("all models must be fit to the same data object")
              header <- paste("Data:", data[[1]])
              subset <- lapply(calls, "[[", "subset")
              if (any(subset != subset[[1]])) stop("all models must use the same subset")
              if (!is.null(subset[[1]]))
                  header <-
                      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
              llks <- lapply(mods, logLik, REML = FALSE)
              Df <- sapply(llks, attr, "df")
              llk <- unlist(llks)
              chisq <- 2 * pmax(0, c(NA, diff(llk)))
              dfChisq <- c(NA, diff(Df))
              val <- data.frame(Df = Df,
                                AIC = sapply(llks, AIC),
                                BIC = sapply(llks, BIC),
                                logLik = llk,
                                "Chisq" = chisq,
                                "Chi Df" = dfChisq,
                                "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
                                check.names = FALSE)
              class(val) <- c("anova", class(val))
              attr(val, "heading") <-
                  c(header, "", "Models:",
                    paste(names(mods),
                          unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
                          sep = ": "),"")
              return(val)
          } else {
#               mCall$terms <- object@terms
#               mCall$assign <- object@assign
          foo <- object
          foo@status["factored"] <- FALSE
          .Call("lmer_factor", foo, PACKAGE="Matrix")
          dfr <- lme4:::getFixDF(foo)
          rcol <- ncol(foo@RXX)
          ss <- foo@RXX[ , rcol]^2
          ssr <- ss[[rcol]]
          ss <- ss[seq(along = dfr)]
          names(ss) <- object@cnames[[".fixed"]][seq(along = dfr)]
          # FIXME: This only gives single degree of freedom tests
          ms <- ss
          df <- rep(1, length(ms))
          f <- ms/(ssr/dfr)
          P <- pf(f, df, dfr, lower.tail = FALSE)
          table <- data.frame(df, ss, ms, dfr, f, P)
          dimnames(table) <-
              list(names(ss),
                   c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
          if (any(match(names(ss), "(Intercept)", nomatch = 0)))
              table <- table[-1,]
          attr(table, "heading") <- "Analysis of Variance Table"
          class(table) <- c("anova", "data.frame")
          table

          }
      })

setMethod("formula", "lmer", function(x, ...) x@call$formula)

setMethod("plot", signature(x = "lmer"),
          function(x, y, ...)
          cat("plot method for lmer not yet implemented\n"))

setMethod("update", signature(object = "lmer"),
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


setMethod("confint", signature(object = "lmer"),
          function (object, parm, level = 0.95, ...) 
      {
          cf <- fixef(object)
          pnames <- names(cf)
          if (missing(parm)) 
              parm <- seq(along = pnames)
          else if (is.character(parm)) 
              parm <- match(parm, pnames, nomatch = 0)
          a <- (1 - level)/2
          a <- c(a, 1 - a)
          pct <- paste(round(100 * a, 1), "%")
          ci <- array(NA, dim = c(length(parm), 2),
                      dimnames = list(pnames[parm], pct))
          ses <- sqrt(diag(vcov(object)))[parm]
          ci[] <- cf[parm] + ses * t(outer(a, getFixDF(object@rep)[parm], qt))
          ci
      })




