setMethod("GLMM1", signature(formula = "formula"),
          function(formula, family, data,
                   method = c("PQL", "Laplace"),
                   control = list(),
                   subset, weights, na.action, offset,
                   model = TRUE, x = FALSE, y = FALSE, ...)
      {
          gVerb <- getOption("verbose")
                                        # match and check parameters
          controlvals <- do.call("lmerControl", control)
          controlvals$REML <- FALSE
          if (length(formula) < 3) stop("formula must be a two-sided formula")

          ## initial glm fit and evaluation of model frame
          mf <- match.call()            
          m <- match(c("family", "data", "subset", "weights",
                       "na.action", "offset"),
                     names(mf), 0)
          mf <- mf[c(1, m)]
          mf[[1]] <- as.name("glm")
          fixed.form <- nobars(formula)
          if (!inherits(fixed.form, "formula")) fixed.form <- ~ 1 # default formula
          environment(fixed.form) <- environment(formula)
          mf$formula <- fixed.form
          mf$x <- mf$model <- mf$y <- TRUE
          glm.fit <- eval(mf, parent.frame())
          family <- glm.fit$family
          offset <- glm.fit$offset
          if (is.null(offset)) offset <- 0
          weights <- sqrt(abs(glm.fit$prior.weights))
          ## initial 'fitted' values on linear scale
          eta <- glm.fit$linear.predictors
          etaold <- eta
          ## Note: offset is on the linear scale
          
          mf$x <- mf$model <- mf$y <- mf$family <- NULL
          mf$drop.unused.levels <- TRUE
          this.form <- subbars(formula)
          environment(this.form) <- environment(formula)
          mf$formula <- this.form
          mf[[1]] <- as.name("model.frame")
          frm <- eval(mf, parent.frame())
          
          ## grouping factors and model matrices for random effects
          bars <- findbars(formula[[3]])
          random <-
              lapply(bars,
                     function(x) list(model.matrix(eval(substitute(~term,
                                                                   list(term=x[[2]]))),
                                                   frm),
                                      eval(substitute(as.factor(fac)[,drop = TRUE],
                                                      list(fac = x[[3]])), frm)))
          names(random) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
          
          ## order factor list by decreasing number of levels
          nlev <- sapply(random, function(x) length(levels(x[[2]])))
          if (any(diff(nlev) > 0)) {
              random <- random[rev(order(nlev))]
          }
          mmats <- c(lapply(random, "[[", 1),
                     .fixed = list(cbind(glm.fit$x, .response = glm.fit$y)))
          ## FIXME: Use Xfrm and Xmat to get the terms and assign
          ## slots, pass these to lmer_create, then destroy Xfrm, Xmat, etc.
          obj <- .Call("lmer_create", lapply(random, "[[", 2), mmats, PACKAGE = "Matrix")
          obj@terms <- attr(glm.fit$model, "terms")
          obj@assign <- attr(glm.fit$x, "assign")
          obj@call <- match.call()
          obj@REML <- FALSE
          rm(glm.fit)
          .Call("lmer_initial", obj, PACKAGE="Matrix")
          mmats.unadjusted <- mmats
          conv <- FALSE
          firstIter <- TRUE
          msMaxIter.orig <- controlvals$msMaxIter

          for (iter in seq(length = controlvals$PQLmaxIt))
          {
              mu <- family$linkinv(eta)
              dmu.deta <- family$mu.eta(eta)
              ## weights (note: weights is already square-rooted)
              w <- weights * dmu.deta / sqrt(family$variance(mu))
              ## adjusted response (should be comparable to X \beta, not including offset
              z <- eta - offset + (mmats.unadjusted$.Xy[, responseIndex] - mu) / dmu.deta
              .Call("nlme_weight_matrix_list",
                    mmats.unadjusted, w, z, mmats, PACKAGE="Matrix")
              .Call("lmer_update_mm", obj, facs, mmats, PACKAGE="Matrix")
              if (firstIter) {
                  .Call("lmer_initial", obj, PACKAGE="Matrix")
                  if (gVerb) cat(" PQL iterations convergence criterion\n")
              }
              .Call("lmer_ECMEsteps", obj, 
                    controlvals$niterEM,
                    FALSE,
                    controlvals$EMverbose,
                    PACKAGE = "Matrix")
              LMEoptimize(obj) <- controlvals
          }
          obj
      })
