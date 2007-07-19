# lmer, glmer and nlmer plus methods and utilities

### FIXME: Move this function to the stats package
rWishart <- function(n, df, invScal)
### Random sample from a Wishart distribution
    .Call(lme4_rWishart, n, df, invScal)

### Utilities for parsing the mixed model formula

findbars <- function(term)
### Return the pairs of expressions that separated by vertical bars
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

nobars <- function(term)
### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

subnms <- function(term, nlist)
### Substitute any names from nlist in term with 1
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

slashTerms <- function(x)
### Return the list of '/'-separated terms in an expression that
### contains slashes
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
### expand any slashes in the grouping factors returned by findbars
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

### Utilities used in lmer, glmer and nlmer

lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
### Create the model frame, X, Y, wts, offset and terms

### mc - matched call of calling function
### formula - two-sided formula
### contrasts - contrasts argument
### vnms - names of variables to be included in the model frame
{
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    ## The model formula for evaluation of the model frame.  It looks
    ## like a linear model formula but includes any random effects
    ## terms and any names of parameters used in a nonlinear mixed model.
    frame.form <- subbars(formula)      # substitute `+' for `|'
    if (length(vnms) > 0)               # add the variables names for nlmer
        frame.form[[3]] <-
            substitute(foo + bar,
                       list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                            bar = frame.form[[3]]))

    ## The model formula for the fixed-effects terms only.
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (inherits(fixed.form, "name"))   # RHS is empty - use `y ~ 1'
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))

    ## attach the correct environment
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf                            # save a copy of the call
    mf <- eval(mf, parent.frame(2))

    ## evaluate the terms for the fixed-effects only (used in anova)
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2)) # allow model.frame to update them
    mt <- attr(fe, "terms")

    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }

    ## Extract X checking for a null model. This check shouldn't be
    ## needed because an empty formula is changed to ~ 1 but it can't hurt.
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
    storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical

    ## Extract the weights and offset.  For S4 classes we want the
    ## `not used' condition to be numeric(0) instead of NULL
    wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
    off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

    ## check weights and offset
    if (any(wts < 0))
        stop(gettextf("negative weights not allowed"), domain = NA)
    if(length(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), NROW(Y), domain = NA))

    ## remove the terms attribute from mf
    attr(mf, "terms") <- NULL
### FIXME: Should I instead attach mt as the terms attribute for mf?
### Will it break anything to have more variables in the frame than
### are used in the terms?  I think any function using the terms
### should be defined here (e.g. the anova method) or use the terms()
### extractor.
    list(Y = Y, X = X, wts = wts, off = off, mt = mt, mf = mf)
}

lmerFactorList <- function(formula, mf, rmInt, drop)
### Create the list of grouping factors and the corresponding
### transposed model matrices.
### rmInt is a logical scalar indicating if the `(Intercept)` column
### should be removed before creating Zt
### drop is a logical scalar indicating if 0 should be dropped
{
    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
             {
                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf)
                 im <- as(ff, "sparseMatrix") # transpose of indicators
                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                    list(expr = x[[2]]))),
                                    mf)
                 if (rmInt) {
                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                     if (ncol(mm) < 2)
                         stop("lhs of a random-effects term cannot be an intercept only")
                     mm <- mm[ , -icol , drop = FALSE]
                 }
                 ans <- list(f = ff,
                             Zt = do.call(rBind,
                             lapply(seq_len(ncol(mm)),
                                    function(j) {im@x <- mm[,j]; im})),
                             cnames = colnames(mm))
                 if (drop) ans$Zt <- drop0(ans$Zt)
                 ans
             })

    ## order factor list by decreasing number of levels but don't
    ## change the order if this is already true
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    ord <- seq_along(nlev)
    if (any(diff(nlev) < 0)) ord <- rev(order(nlev))
    fl[ord]
### FIXME: 1. Detect and collapse repeated factors. Use an integer
### index vector to map random-effects terms to factors.
### 2. Check to see if the factors form a nested sequence after
### removing repeated factors (for information only).
}

checkSTform <- function(ST, STnew)
### Check that the 'STnew' argument matches the form of ST.
{
    stopifnot(is.list(STnew), length(STnew) == length(ST),
              all.equal(names(ST), names(STnew)))
    lapply(seq_along(STnew), function (i)
           stopifnot(class(STnew[[i]]) == class(ST[[i]]),
                     all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
    all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

lmerControl <- function(msVerbose = getOption("verbose"))
### Control parameters for lmer, glmer and nlmer
{
    list(
### FIXME: Should the user have control of maxIter and tolerance? If
### so, how should they be passed to the lmer_optimize C function?
         ## maxIter = as.integer(maxIter),
         ## tolerance = as.double(tolerance),
	 msVerbose = as.integer(msVerbose))# "integer" on purpose
}

VecFromNames <- function(nms, mode = "numeric")
### Generate a named vector of the given mode
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans
}

mkdims <- function(fr, FL, start)
### Create the standard versions of flist, Zt, Gp, ST, cnames and dd
{
    flist <- lapply(FL, get("[["), "f")
    Ztl <- lapply(FL, get("[["), "Zt")
    cnames <- lapply(FL, get("[["), "cnames")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    Gp <- unname(c(0L, cumsum(sapply(Ztl, nrow))))
    rm(Ztl, FL)                         # because they could be large
    nc <- sapply(cnames, length)        # # of columns in els of ST
    ST <- lapply(nc, function(n) matrix(0, n, n))
    .Call(ST_initialize, ST, Gp, Zt)
    if (!is.null(start) && checkSTform(ST, start)) ST <- start
    Vt <- .Call(mer_create_Vt, Zt, ST, Gp)

    ## record dimensions and algorithm settings
    dd <-
        VecFromNames(c("nf", "n", "p", "q", "s", "np", "REML", "ftyp", "nest"),
                       "integer")
    dd["nf"] <- length(cnames)          # number of random-effects terms
    dd["n"] <- nrow(fr$mf)              # number of observations
    dd["p"] <- ncol(fr$X)               # number of fixed-effects coefficients
    dd["q"] <- nrow(Zt)                 # number of random effects
    dd["s"] <- 1L                       # always 1 except in nlmer
    nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
### FIXME: Check number of variance components versus number of
### levels in the factor for each term. Warn or stop as appropriate
    dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
    dd["REML"] <- 0L                    # glmer and nlmer don't use REML

    list(Gp = Gp, ST = ST, Vt = Vt, Zt = Zt,
         cnames = cnames, dd = dd,
         dev = VecFromNames(
         c("ML", "REML", "ldL2", "ldRX2", "lpdisc", "disc", "bqd"),
         "numeric"), flist = flist)
}

mkFltype <- function(family)
### check for predefined families
{
    fltype <- 0                         # not a predefined type
    if (family$family == "gaussian" && family$link == "identity") fltype <- -1
    if (family$family == "binomial" && family$link == "logit") fltype <- 1
    if (family$family == "binomial" && family$link == "probit") fltype <- 2
    if (family$family == "poisson" && family$link == "log") fltype <- 3
    as.integer(fltype)
}

mkFamilyEnv <- function(glmFit)
### Create and populate the family function evaluation environment.
{
    env <- new.env()
    assign("devResid", unname(resid(glmFit, type = "deviance")), env = env)
    assign("eta", unname(glmFit$linear.predictors), env = env)
    assign("mu", unname(glmFit$fitted.values), env = env)
    assign("weights", unname(glmFit$weights), env = env)
    assign("z", unname(glmFit$residuals), env = env)
### FIXME: install the family functions and create evaluation expressions in the environment
    env
}

### The main event

lmer <-
    function(formula, data, family = NULL, method = c("REML", "ML"),
             control = list(), start = NULL, verbose, subset,
             weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
### Linear Mixed-Effects in R
{
    mc <- match.call()
    if (!is.null(family)) {             # call glmer
        mc[[1]] <- as.name("glmer")
        return(eval(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lmerFactorList(formula, fr$mf, 0L, 0L) # flist, Zt, cnames
    Y <- as.double(fr$Y)
    dm <- mkdims(fr, FL, start)
    stopifnot(length(levels(dm$flist[[1]])) < length(Y))
### FIXME: A kinder, gentler error message may be in order.
### This checks that the number of levels in a grouping factor < n
### Only need to check the first factor because it is the one with
### the most levels.
    dm$dd["REML"] <- match.arg(method) == "REML"
    dm$dd["ftyp"] <- -1L              # gaussian family, identity link

    ## Create the dense matrices to be used in the deviance evaluation
### FIXME: incorporate weights and offset in the creation of ZtXy et al.
    Xy <- cbind(fr$X, Y)
    ZtXy <- as(dm$Zt %*% Xy, "matrix")
    rownames(fr$X) <- dimnames(ZtXy) <- dimnames(Xy) <- NULL
    RXy <- XytXy <- crossprod(Xy)
    RXy[] <- 0
    fixef <- numeric(dm$dd["p"])
    names(fixef) <- colnames(fr$X)
    dimnames(fr$X) <- NULL

    ans <- new("lmer",
               frame = if (model) fr$mf else fr$mf[0,],
               call = mc, terms = fr$mt, flist = dm$flist,
               Zt = dm$Zt, X = if (x) fr$X else fr$X[0,],
### FIXME: Should y retain its names? As it stands any row names in the
### frame are dropped.  Really?  Are they part of the frame slot (if not
### reduced to 0 rows)?
               y = unname(Y), ZtXy = ZtXy, XytXy = XytXy,
               weights = unname(fr$wts), offset = unname(fr$off),
               cnames = unname(dm$cnames), Gp = unname(dm$Gp),
               dims = dm$dd, ST = dm$ST, Vt = dm$Vt,
               L = .Call(mer_create_L, dm$Vt),
               RXy = RXy, RVXy = ZtXy,
               deviance = dm$dev, fixef = fixef,
               ranef = numeric(dm$dd["q"]),
               uvec = numeric(dm$dd["q"]))
    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
    .Call(mer_optimize, ans, verbose, 0)
    .Call(lmer_update_effects, ans)
    ans
}

glmer <-
function(formula, data, family = gaussian, method = c("Laplace", "AGQ"),
         control = list(), start = NULL, verbose, subset, weights,
         na.action, offset, contrasts = NULL, model = TRUE, ...)
### Fit a generalized linear mixed model
{
    warning("glmer is in development.  Results reported here are wrong.")
    mc <- match.call()
                                        # Evaluate and check the family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if(is.null(family$family)) stop("'family' not recognized")
    if(family$family == "gaussian" && family$link == "identity") {
        mc[[1]] <- "lmer"               # use lmer not glmer
        mc$family <- NULL
        return(eval.parent(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    glmFit <- glm.fit(fr$X, fr$Y, weights = fr$weights, # glm on f.e.
                      offset = fr$offset, family = family,
                      intercept = attr(fr$mt, "intercept") > 0)
    FL <- lmerFactorList(formula, fr$mf, 0L, 0L) # flist, Zt, cnames
    dm <- mkdims(fr, FL, start)
    dm$dd["ftyp"] <- mkFltype(glmFit$family)
    dimnames(fr$X) <- NULL

    ans <- new("glmer",
               env = mkFamilyEnv(glmFit),
               famName = unlist(glmFit$family[c("family", "link")]),
               frame = if (model) fr$mf else fr$mf[0,],
               call = mc, terms = fr$mt, flist = dm$flist,
               Zt = dm$Zt, X = fr$X,
               y = unname(as.double(glmFit$y)),
               weights = unname(glmFit$prior.weights),
               offset = unname(fr$off),
               cnames = unname(dm$cnames), Gp = unname(dm$Gp),
               dims = dm$dd, ST = dm$ST, Vt = dm$Vt,
               L = .Call(mer_create_L, dm$Vt),
               deviance = dm$dev,
               fixef = coef(glmFit),
               ranef = numeric(dm$dd["q"]),
               uvec = numeric(dm$dd["q"]))
    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
#    .Call(mer_optimize, ans, verbose, 2)
    ans
}

nlmer <- function(formula, data, control = list(), start = NULL,
                  verbose = FALSE, subset, weights, na.action,
                  contrasts = NULL, model = TRUE, ...)
### Fit a nonlinear mixed-effects model
{
    ##    warning("nlmer is in development.  Results reported here are wrong.")
    mc <- match.call()
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
    if (length(nlform) < 3)
        stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
    if (is.numeric(start)) start <- list(fixed = start)
    s <- length(pnames <- names(start$fixed))
    stopifnot(length(start$fixed) > 0, s > 0,
              inherits(data, "data.frame"), nrow(data) > 1)
### FIXME: Allow for a situation where data is not specified.  What
### should it default to?
    if (any(pnames %in% names(data)))
        stop("parameter names must be distinct from names of the variables in data")
    anms <- all.vars(nlmod)
    if (!all(pnames %in% anms))
        stop("not all parameter names are used in the nonlinear model expression")

    if (!length(vnms <- setdiff(anms, pnames)))
        stop("there are no variables used in the nonlinear model expression")
    ## create a frame in which to evaluate the factor list
    fr <- lmerFrames(mc,
                     eval(substitute(foo ~ bar,
                                     list(foo = nlform[[2]],
                                          bar = subnms(formula[[3]], lapply(pnames, as.name))))),
                     contrasts, vnms)
    mf <- fr$mf
    env <- new.env()
    lapply(names(mf), function(nm) assign(nm, env = env, mf[[nm]]))
    n <- nrow(mf)
    lapply(pnames, function(nm) assign(nm, env = env, rep(start$fixed[[nm]], length.out = n)))

    n <- nrow(mf)
    mf <- mf[rep(seq_len(n), s), ]
    row.names(mf) <- seq_len(nrow(mf))
    ss <- rep.int(n, s)
    for (nm in pnames) mf[[nm]] <- rep.int(as.numeric(nm == pnames), ss)
                                        # factor list and model matrices
    FL <- lmerFactorList(substitute(foo ~ bar, list(foo = nlform[[2]], bar = formula[[3]])),
                         mf, TRUE, TRUE)
    dm <- mkdims(fr, FL, start$STpars)
    dm$dd["s"] <- s
    dm$dd["ftyp"] <- -1L              # gaussian family, identity link
    dm$dd["p"] <- length(start$fixed)

### FIXME: Probably store X as a dense matrix.  Sparse storage may actually be larger.
### Then the X slot can be moved from the lmer and glmer representations to mer.
    X <- as.matrix(mf[,pnames])
    rownames(X) <- NULL

### FIXME: The only times there would be additional columns in the
### fixed effects would be as interactions with parameter names and
### they must be constructed differently
    xnms <- colnames(fr$X)
    if (!is.na(icol <- match("(Intercept)",xnms))) xnms <- xnms[-icol]
    if (length(xnms) > 0)
        Xt <- rBind(Xt, t(Matrix(fr$X[rep.int(seq_len(n), s), xnms, drop = FALSE])))

    Mt <- .Call(nlmer_create_Mt, dm$Vt, s)
    ans <- new("nlmer",
               env = env, model = nlmod, pnames = pnames,
               mu = numeric(dm$dd["n"]), Mt = Mt,
               frame = if (model) fr$mf else fr$mf[0,],
               call = mc, terms = fr$mt, flist = dm$flist, X = X,
               Zt = dm$Zt, Vt = dm$Vt, y = unname(as.double(fr$Y)),
               weights = unname(fr$wts), 
               cnames = unname(dm$cnames), Gp = unname(dm$Gp),
               dims = dm$dd, ST = dm$ST,
               L = .Call(mer_create_L, Mt),
               deviance = dm$dev, fixef = start$fixed,
               ranef = numeric(dm$dd["q"]),
               uvec = numeric(dm$dd["q"]))
    .Call(nlmer_eval_model, ans)
    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
#    .Call(mer_optimize, ans, verbose, 1L)
    ans
}

#### Extractors specific to mixed-effects models

setMethod("coef", signature(object = "mer"),
	  function(object, ...)
      {
          if (length(list(...)))
              warning(paste('arguments named "',
                            paste(names(list(...)), collapse = ", "),
                                  '" ignored', sep = ''))
          fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
          ref <- ranef(object)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
          for (i in seq(a = val)) {
              refi <- ref[[i]]
              row.names(val[[i]]) <- row.names(refi)
              nmsi <- colnames(refi)
              if (!all(nmsi %in% names(fef)))
                  stop("unable to align random and fixed effects")
              for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
          }
          new("coef.lmer", val)
       })

setAs("mer", "dtCMatrix", function(from)
### Extract the L matrix
      as(from@L, "sparseMatrix"))

setMethod("fixef", signature(object = "mer"),
          function(object, ...)
### Extract the fixed effects
          object@fixef)

setMethod("ranef", signature(object = "mer"),
	  function(object, postVar = FALSE, ...)
### Extract the random effects
### FIXME: This will need to be modified if flist is collapsed
      {
          if (postVar) .NotYetImplemented()
          Gp <- object@Gp
          ii <- lapply(diff(Gp), seq_len)
          rr <- object@ranef
          cn <- object@cnames
          rn <- lapply(object@flist, levels)
          ans <-
              lapply(seq_len(length(ii)),
                     function(i)
                     data.frame(matrix(rr[ii[[i]] + Gp[i]],
                                       nc = length(cn[[i]]),
                                       dimnames = list(rn[[i]], cn[[i]])),
                                check.names = FALSE))
          names(ans) <- names(object@ST)
##               if (postVar) {
##                   pV <- .Call(mer_postVar, object)
##                   for (i in seq(along = ans))
##                       attr(ans[[i]], "postVar") <- pV[[i]]
##               }
          ans
      })

setMethod("VarCorr", signature(x = "mer"),
	  function(x, REML = NULL, ...)
### Create the VarCorr object of variances and covariances
### FIXME: Modify this to handle glmer models without a scale parameter
      {
	  sc <- .Call(mer_sigma, x, REML)
	  cnames <- x@cnames
	  ans <- x@ST
          for (i in seq(along = ans)) {
              ai <- ans[[i]]
              dm <- dim(ai)
              if (dm[1] < 2) {
                  el <- (sc * ai)^2
              } else {
                  dd <- diag(ai)
                  diag(ai) <- rep(1, dm[1])
                  el <- sc^2 * crossprod(dd * t(ai))
              }
              el <- as(el, "dpoMatrix")
              el@Dimnames <- list(cnames[[i]], cnames[[i]])
	      el@factors$correlation <- as(el, "corMatrix")
	      ans[[i]] <- el
	  }
	  attr(ans, "sc") <- sc
	  ans
      })

#### Methods for standard extractors for fitted models

setMethod("anova", signature(object = "lmer"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "lmer") | sapply(dots, is, "lm") else logical(0)
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
					attr, "df"))]
	      calls <- lapply(mods, slot, "call")
	      data <- lapply(calls, "[[", "data")
	      if (any(data != data[[1]]))
		  stop("all models must be fit to the same data object")
	      header <- paste("Data:", data[[1]])
	      subset <- lapply(calls, "[[", "subset")
	      if (any(subset != subset[[1]]))
		  stop("all models must use the same subset")
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
		  c(header, "Models:",
		    paste(names(mods),
			  unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
			  sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
              p <- object@dims["p"]
	      ss <- (object@RXy[seq_len(p), p + 1L, drop = TRUE])^2
	      names(ss) <- names(object@fixef)
	      asgn <- attr(object@X, "assign")
	      terms <- terms(object)
	      nmeffects <- attr(terms, "term.labels")
	      if ("(Intercept)" %in% names(ss))
		  nmeffects <- c("(Intercept)", nmeffects)
	      ss <- unlist(lapply(split(ss, asgn), sum))
	      df <- unlist(lapply(split(asgn,  asgn), length))
	      #dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	      ms <- ss/df
	      f <- ms/(.Call(mer_sigma, object, 0L)^2)
	      #P <- pf(f, df, dfr, lower.tail = FALSE)
	      #table <- data.frame(df, ss, ms, dfr, f, P)
	      table <- data.frame(df, ss, ms, f)
	      dimnames(table) <-
		  list(nmeffects,
#			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		       c("Df", "Sum Sq", "Mean Sq", "F value"))
	      if ("(Intercept)" %in% nmeffects)
		  table <- table[-match("(Intercept)", nmeffects), ]
	      attr(table, "heading") <- "Analysis of Variance Table"
	      class(table) <- c("anova", "data.frame")
	      table
	  }
      })

setMethod("confint", signature(object = "lmer"),
	  function(object, parm, level = 0.95, ...)
	  .NotYetImplemented()
	  )

setMethod("deviance", signature(object="mer"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- object@dims["REML"]
          object@deviance[ifelse(REML, "REML", "ML")]
      })

setMethod("fitted", signature(object = "lmer"),
	  function(object, ...)
	  .NotYetImplemented()
	  )

setMethod("formula", signature(x = "mer"),
	  function(x, ...)
	  x@call$formula
	  )

setMethod("logLik", signature(object="mer"),
	  function(object, REML = NULL, ...)
### Extract the log-likelihood or restricted log-likelihood
      {
          dims <- object@dims
          if (is.null(REML) || is.na(REML[1]))
              REML <- object@dims["REML"]
          val <- -deviance(object, REML = REML)/2
          attr(val, "nall") <- attr(val, "nobs") <- dims["n"]
          attr(val, "df") <-
              dims["p"] + length(.Call(ST_getPars, object))
          attr(val, "REML") <-  as.logical(REML)
          class(val) <- "logLik"
          val
      })

setMethod("residuals", signature(object = "glmer"),
	  function(object, ...) .NotYetImplemented())

setMethod("residuals", signature(object = "lmer"),
	  function(object, ...) object@y - fitted(object))

### FIXME: There should not be two identical methods like this but I'm not
##        sure how to pass the ... argument to a method for another generic
##        cleanly.
setMethod("resid", signature(object = "glmer"),
	  function(object, ...) .NotYetImplemented())

setMethod("resid", signature(object = "lmer"),
	  function(object, ...) object@y - fitted(object))

setMethod("simulate", "lmer",
          function(object, nsim = 1, seed = NULL, ...)
      {
          .NotYetImplemented()
	  if(!is.null(seed)) set.seed(seed)
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1)		     # initialize the RNG if necessary
          RNGstate <- .Random.seed
          dims <- object@dims
          ans <- array(0, c(dims["n"], nsim))
          fe <- fixef(object)
          re <- ranef(object)
          nc <- sapply(re, ncol)
          nr <- sapply(re, nrow)
          sz <- nc * nr
          vc <- VarCorr(object)

          cholL <- lapply(vc, chol)
          n <- object@dims["n"]
          for (i in seq_len(nsim))
              ans[, 1] <- crossprod(object@ZXyt,
                                    c(unlist(lapply(seq_along(re), function(k)
                                                    (t(cholL[[k]]) %*%
                                                     matrix(rnorm(sz[k]),
                                                            nc = nr[k]))@x)),
                                      fe, 0))@x
          ans + rnorm(prod(dim(ans)), sd = attr(vc, "sc"))
      })

setMethod("summary", signature(object = "mer"),
	  function(object, ...)
      {
          fcoef <- fixef(object)
          vcov <- vcov(object)
          corF <- vcov@factors$correlation
          dims <- object@dims
          ## DF <- getFixDF(object)
          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
          REML <- object@dims["REML"]
          llik <- logLik(object, REML)
          dev <- object@deviance

          glz <- is(object, "glmer")
          methTitle <-
              if (glz)
                  "Generalized linear mixed model fit using Laplace"
##                   paste("Generalized linear mixed model fit using ",
##                         switch(object@status["glmm"],
##                                "PQL", "Laplace", "AGQ"))
              else paste("Linear mixed-effects model fit by",
                         if(REML) "REML" else "maximum likelihood")

          AICframe <- {
              if (glz)
                  data.frame(AIC = AIC(llik), BIC = BIC(llik),
                             logLik = c(llik),
                             deviance = -2*llik,
                             row.names = "")
              else
                  data.frame(AIC = AIC(llik), BIC = BIC(llik),
                             logLik = c(llik),
                             MLdeviance = dev["ML"],
                             REMLdeviance = dev["REML"],
                             row.names = "")
          }
          varcor <- VarCorr(object)
          REmat <- formatVC(varcor)
          if (is.na(attr(varcor, "sc")))
              REmat <- REmat[-nrow(REmat), , drop = FALSE]

          if (nrow(coefs) > 0) {
              if (dims["ftyp"] >= 0) {
                  coefs <- coefs[, 1:2, drop = FALSE]
                  stat <- coefs[,1]/coefs[,2]
                  pval <- 2*pnorm(abs(stat), lower = FALSE)
                  coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
              } else {
                  stat <- coefs[,1]/coefs[,2]
                  ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                  coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
              }
          } ## else : append columns to 0-row matrix ...
          oclass <- ""
          if (is(object, "glmer")) oclass <- "summary.glmer"
          if (is(object, "lmer")) oclass <- "summary.lmer"
          if (is(object, "nlmer")) oclass <- "summary.nlmer"
          if (!nchar(oclass))
              stop(gettextf("unrecognized class %s", class(object)))
          new(oclass,
              object,
              methTitle = methTitle,
              logLik = llik,
              ngrps = sapply(object@flist, function(x) length(levels(x))),
              sigma = .Call(mer_sigma, object, REML),
              coefs = coefs,
              vcov = vcov,
              REmat = REmat,
              AICtab= AICframe
              )
      })## summary()

setMethod("terms", signature(x = "mer"),
	  function(x, ...) x@terms)

setMethod("update", signature(object = "mer"),
	  function(object, formula., ..., evaluate = TRUE)
      {
	  call <- object@call
	  if (is.null(call))
	      stop("need an object with call slot")
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

setMethod("vcov", signature(object = "lmer"),
	  function(object, REML = 0, ...)
### Extract the conditional variance-covariance matrix of the fixed effects
      {
          rr <- as(.Call(mer_sigma, object, REML)^2 *
                   chol2inv(object@RXy, size = object@dims['p']), "dpoMatrix")
          nms <- colnames(object@X)
          dimnames(rr) <- list(nms, nms)
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })

setMethod("vcov", "mer",                # placeholder
          function (object, ...)
      {
          rr <- as(diag(length(object@fixef)), "dpoMatrix")
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })


setMethod("with", signature(data = "lmer"),
	  function(data, expr, ...) {
	      dat <- eval(data@call$data)
	      if (!is.null(na.act <- attr(data@frame, "na.action")))
		  dat <- dat[-na.act, ]
	      lst <- c(list(. = data), data@flist, data@frame, dat)
	      eval(substitute(expr), lst[unique(names(lst))])
	  })

### Show and print methods and utilities for them

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{  
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, function(el) el@factors$correlation)
    reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}

## This is modeled a bit after  print.summary.lm :
printMer <- function(x, digits = max(3, getOption("digits") - 3),
                      correlation = TRUE, symbolic.cor = FALSE,
                      signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    REML <- so@dims["REML"]
    llik <- so@logLik
    dev <- so@deviance
    dims <- x@dims

    cat(so@methTitle, "\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x@call$subset)[[2]]),"\n")
    if (inherits(x, "glmer"))
        cat(" Family: ", so@famName["family"], "(",
            so@famName["link"], " link)\n", sep = "")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    cat(sprintf("Number of obs: %d, groups: ", dims["n"]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (is.na(so@sigma))
	cat("\nEstimated scale (compare to 1):",
            sqrt(exp(so@deviance["lr2"])/so@dims["n"]), "\n")
    if (nrow(so@coefs) > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(correlation) {
	    rn <- rownames(so@coefs)
	    corF <- so@vcov@factors$correlation
	    if (!is.null(corF)) {
		p <- ncol(corF)
		if (p > 1) {
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			print(symnum(as(corF, "matrix"), abbr.col = NULL))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p)
			dimnames(corf) <- list(abbreviate(rn, minlen=11),
					       abbreviate(rn, minlen=6))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}

setMethod("print", "lmer", printMer)
setMethod("print", "glmer", printMer)
setMethod("show", "lmer", function(object) printMer(object))
setMethod("show", "glmer", function(object) printMer(object))

setMethod("show", "nlmer", function(object)
### FIXME: Does nlmer need a separate show method?
      {
          dims <- object@dims
          cat("Nonlinear mixed model fit by Laplace\n")
          if (!is.null(object@call$formula))
              cat("Formula:", deparse(object@call$formula),"\n")
          if (!is.null(object@call$data))
              cat("   Data:", deparse(object@call$data), "\n")
          if (!is.null(object@call$subset))
              cat(" Subset:",
                  deparse(asOneSidedFormula(object@call$subset)[[2]]),"\n")

          cat("Random effects:\n")
          print(formatVC(VarCorr(object)), quote = FALSE,
                digits = max(3, getOption("digits") - 3))

          cat(sprintf("Number of obs: %d, groups: ", dims["n"]))
          ngrps <- sapply(object@flist, function(x) length(levels(x)))
          cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
          cat("\n")
          cat("\nFixed effects:\n")
          print(object@fixef)
          invisible(object)
      })

#### Methods for secondary, derived classes

setMethod("deviance", signature(object = "summary.mer"), function(object) object@deviance)
setMethod("logLik", signature(object = "summary.mer"), function(object) object@logLik)
setMethod("vcov", signature(object = "summary.mer"), function(object) object@vcov)
setMethod("summary", signature(object = "summary.mer"), function(object) object)

#### Methods to produce specific plots

setMethod("plot", signature(x = "coef.lmer"),
          function(x, y, ...)
      {
          varying <- unique(do.call("c",
                                    lapply(x, function(el)
                                           names(el)[sapply(el,
                                                            function(col)
                                                            any(col != col[1]))])))
          gf <- do.call("rBind", lapply(x, "[", j = varying))
          gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
          switch(min(length(varying), 3),
                 qqmath(eval(substitute(~ x | .grp,
                                        list(x = as.name(varying[1])))), gf, ...),
                 xyplot(eval(substitute(y ~ x | .grp,
                                        list(y = as.name(varying[1]),
                                             x = as.name(varying[2])))), gf, ...),
                 splom(~ gf | .grp, ...))
      })

setMethod("plot", signature(x = "ranef.lmer"),
	  function(x, y, ...)
      {
	  lapply(x, function(x) {
	      cn <- lapply(colnames(x), as.name)
	      switch(min(ncol(x), 3),
		     qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
		     xyplot(eval(substitute(y ~ x,
					    list(y = cn[[1]],
						 x = cn[[2]]))), x, ...),
		     splom(~ x, ...))
	  })
      })

setMethod("qqmath", signature(x = "ranef.lmer"),
          function(x, data, ...) {
              prepanel.ci <- function(x, y, se, subscripts, ...) {
                  y <- as.numeric(y)
                  se <- as.numeric(se[subscripts])
                  hw <- 1.96 * se
                  list(ylim = range(y - hw, y + hw, finite = TRUE))
              }
              panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
                  panel.grid(h = -1,v = -1)
                  panel.abline(h = 0)
                  x <- as.numeric(x)
                  y <- as.numeric(y)
                  se <- as.numeric(se[subscripts])
                  ly <- y - 1.96 * se
                  uy <- y + 1.96 * se
                  panel.segments(x, y - 1.96*se, x, y + 1.96 * se,
                                 col = 'black')
                  panel.xyplot(x, y, pch = pch, ...)
              }
              f <- function(x) {
                  if (!is.null(pv <- attr(x, "postVar"))) {
                      cols <- 1:(dim(pv)[1])
                      se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
                      nr <- nrow(x)
                      nc <- ncol(x)
                      ord <- unlist(lapply(x, order)) +
                          rep((0:(nc - 1)) * nr, each = nr)
                      rr <- 1:nr
                      ind <- gl(ncol(x), nrow(x), labels = names(x))
                      xyplot(unlist(x)[ord] ~
                             rep(qnorm((rr - 0.5)/nr), ncol(x)) | ind[ord],
                             se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                             scales = list(y = list(relation = "free")),
                             xlab = "Standard normal quantiles",
                             ylab = NULL, aspect = 1, ...)
                  } else {
                      qqmath(~values|ind, stack(x),
                             scales = list(y = list(relation = "free")),
                             xlab = "Standard normal quantiles",
                             ylab = NULL, ...)
                  }
              }
              lapply(x, f)
          })


#### Creating and displaying a Markov Chain Monte Carlo sample from
#### the posterior distribution of the parameters

setMethod("mcmcsamp", signature(object = "lmer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
		   trans = TRUE, deviance = FALSE, ...)
### Generate a Markov chain Monte Carlo sample from the posterior distribution
### of the parameters in a linear mixed model
      {
          .NotYetImplemented()
          ans <- t(.Call(lmer_MCMCsamp, object, saveb, n,
                         trans, verbose, deviance))
	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
	  class(ans) <- "mcmc"
	  mcmccompnames(ans, object, saveb, trans,
			glmer=FALSE, deviance=deviance)
      })

abbrvNms <- function(gnm, cnms)
### Abbreviate names of columns in grouping factors
### gnm - group name
### cnms - column names
{
    ans <- paste(abbreviate(gnm), abbreviate(cnms), sep = '.')
    if (length(cnms) > 1) {
	anms <- lapply(cnms, abbreviate, minlength = 3)
	nmmat <- outer(anms, anms, paste, sep = '.')
	ans <- c(ans, paste(abbreviate(gnm, minlength = 3),
			    nmmat[upper.tri(nmmat)], sep = '.'))
    }
    ans
}

mcmccompnames <- function(ans, object, saveb, trans, glmer, deviance)
### Mangle the names of the columns of the mcmcsamp result ans
### This operation is common to the methods for "lmer" and "glmer"
{
    .NotYetImplemented()
    gnms <- names(object@flist)
    cnms <- object@cnames
    ff <- fixef(object)
    colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                unlist(lapply(seq(along = gnms),
                              function(i)
                              abbrvNms(gnms[i],cnms[[i]]))))
    if (trans) {
        ## parameter type: 0 => fixed effect, 1 => variance,
        ##		 2 => covariance
        ptyp <- c(integer(length(ff)), if (glmer) integer(0) else 1:1,
                  unlist(lapply(seq(along = gnms),
                                function(i)
                            {
                                k <- length(cnms[[i]])
                                rep(1:2, c(k, (k*(k-1))/2))
                            })))
        colnms[ptyp == 1] <-
            paste("log(", colnms[ptyp == 1], ")", sep = "")
        colnms[ptyp == 2] <-
            paste("atanh(", colnms[ptyp == 2], ")", sep = "")
    }
    if (deviance) colnms <- c(colnms, "deviance")
### FIXME: this will fail for a mer2 object
    if(saveb) {## maybe better colnames, "RE.1","RE.2", ... ?
        rZy <- object@rZy
        colnms <- c(colnms,
                    paste("b", sprintf(paste("%0",
                                             1+floor(log(length(rZy),10)),
                                             "d", sep = ''),
                                       seq(along = rZy)),
                          sep = '.'))
    }
    colnames(ans) <- colnms
    ans
}

#### Odds and ends

## simulestimate <- function(x, FUN, nsim = 1, seed = NULL, control = list())
## {
##     FUN <- match.fun(FUN)
##     stopifnot((nsim <- as.integer(nsim[1])) > 0,
## 	      inherits(x, "lmer"))
##     if (!is.null(seed)) set.seed(seed)
##     ## simulate the linear predictors
##     lpred <- .Call(mer_simulate, x, nsim)
##     sc <- abs(x@devComp[8])
##     ## add fixed-effects contribution and per-observation noise term
##     lpred <- lpred + drop(x@X %*% fixef(x)) + rnorm(prod(dim(lpred)), sd = sc)

##     cv <- do.call(lmerControl, control)
##     Omega <- x@Omega
##     x@wrkres <- x@y <- lpred[,1]
##     .Call(mer_update_ZXy, x)
##     LMEoptimize(x) <- cv
##     template <- FUN(x)
##     if (!is.numeric(template))
##         stop("simulestimate currently only handles functions that return numeric vectors")
##     ans <- matrix(template, nr = nsim, nc = length(template), byrow = TRUE)
##     colnames(ans) <- names(template)
##     for (i in 1:nsim) {
##         x@wrkres <- x@y <- lpred[,i]
##         x@Omega <- Omega
##         .Call(mer_update_ZXy, x)
##         LMEoptimize(x) <- cv
##         foo <- try(FUN(x))
##         ans[i,] <- if (inherits(foo, "try-error")) NA else foo
##     }
##     ans
## }

hatTrace <- function(x)
{
    .NotYetImplemented()
    stopifnot(is(x, "mer"))
##     .Call(mer_hat_trace2, x)
}

ST2Omega <- function(ST)
### Temporary function to convert the ST representation of the
### relative variance-covariance matrix returned by lmer into the
### Omega representation required by lmer
{
    if (nrow(ST) == 1) return(as(1/ST^2, "dpoMatrix"))
    dd <- diag(ST)
    T <- as(ST, "dtrMatrix")
    T@diag <- "U"
    crossprod(solve(T)/dd)
}


## Extract the random effects
## setMethod("ranef", signature(object = "mer"),
## 	  function(object, postVar = FALSE, ...) {
## 	      ans <- new("ranef.lmer",
##                          lapply(.Call(mer_ranef, object),
##                                 data.frame, check.names = FALSE))
##               names(ans) <- names(object@flist)
##               if (postVar) {
##                   pV <- .Call(mer_postVar, object)
##                   for (i in seq(along = ans))
##                       attr(ans[[i]], "postVar") <- pV[[i]]
##               }
##               ans
## 	  })

## setMethod("mcmcsamp", signature(object = "lmer"),
## 	  function(object, n = 1, verbose = FALSE, saveb = FALSE,
## 		   trans = TRUE, deviance = FALSE, ...)
##       {
##           ans <- t(.Call(mer_MCMCsamp, object, saveb, n, trans, verbose, deviance))
## 	  attr(ans, "mcpar") <- as.integer(c(1, n, 1))
## 	  class(ans) <- "mcmc"
## 	  mcmccompnames(ans, object, saveb, trans,
## 			glmer=FALSE, deviance=deviance)
##       })

## setMethod("simulate", signature(object = "mer"),
## 	  function(object, nsim = 1, seed = NULL, ...)
##       {
## 	  if(!exists(".Random.seed", envir = .GlobalEnv))
## 	      runif(1)		     # initialize the RNG if necessary
## 	  if(is.null(seed))
## 	      RNGstate <- .Random.seed
## 	  else {
## 	      R.seed <- .Random.seed
## 	      set.seed(seed)
## 	      RNGstate <- structure(seed, kind = as.list(RNGkind()))
## 	      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
## 	  }

##           stopifnot((nsim <- as.integer(nsim[1])) > 0,
##                     inherits(object, "lmer"))
## 	  ## similate the linear predictors
## 	  lpred <- .Call(mer_simulate, object, nsim)
## 	  sc <- abs(object@devComp[8])

## 	  ## add fixed-effects contribution and per-observation noise term
## 	  lpred <- as.data.frame(lpred + drop(object@X %*% fixef(object)) +
## 				 rnorm(prod(dim(lpred)), sd = sc))
## 	  ## save the seed
## 	  attr(lpred, "seed") <- RNGstate
## 	  lpred
##       })

## We need to define an S4 print method, since using an S3 print
## method fails as soon as you call print() explicitly, e.g. when
## wanting to specify options.

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

## setMethod("getFixDF", signature(object="mer"),
## 	  function(object, ...) {
## 	      devc <- object@devComp
## 	      rep(as.integer(devc[1]- devc[2]), devc[2])
## 	  })

## simss <- function(fm0, fma, nsim)
## {
##     ysim <- simulate(fm0, nsim)
##     cv <- list(gradient = FALSE, msMaxIter = 200:200,
## 	       msVerbose = 0:0)
##     sapply(ysim, function(yy) {
## 	.Call(mer_update_y, fm0, yy)
## 	LMEoptimize(fm0) <- cv
## 	.Call(mer_update_y, fma, yy)
## 	LMEoptimize(fma) <- cv
## 	exp(c(H0 = fm0@devComp[["logryy2"]],
## 	      Ha = fma@devComp[["logryy2"]]))
##     })
## }

## setMethod("denomDF", "mer",
##           function(x, ...)
##       {
##           mm <- x@X
##           aa <- attr(mm, "assign")
##           tt <- x@terms
##           if (!isNested(x))
##               return(list(coef = as.numeric(rep(NA, length(x@fixef))),
##                           terms = as.numeric(rep(NA,
##                           length(attr(tt, "order"))))))
##           hasintercept <- attr(tt, "intercept") > 0
##           ## check which variables vary within levels of grouping factors
##           vars <- eval(attr(tt, "variables"), x@frame)
##           fl <- x@flist
##           vv <- matrix(0:0, nrow = length(vars), ncol = length(fl),
##                         dimnames = list(NULL, names(fl)))
##           ## replace this loop by C code.
##           for (i in 1:nrow(ans))        # check if variables vary within factors
##               for (j in 1:ncol(ans))
##                   ans[i,j] <- all(tapply(vars[[i]], fl[[j]],
##                                          function(x) length(unique(x)) == 1))
##           ## which terms vary within levels of which grouping factors?
##           tv <- crossprod(attr(tt, "factors"), !ans)
##           ## maximum level at which the term is constant
##           ml <- apply(tv, 1, function(rr) max(0, which(as.logical(rr))))
##           ## unravel assignment applied to terms
##           ll <- attr(tt, "term.labels")
##           if (hasintercept)
##               ll <- c("(Intercept)", ll)
##           aaa <- factor(aa, labels = ll)
##           asgn <- split(order(aa), aaa)
##           nco <- lapply(asgn, length)   # number of coefficients per term
##           nlev <- lapply(fl, function(x) length(levels(x)))
##           if (hasintercept) asgn$"(Intercept)" <- NULL
##           list(ml = ml, nco = nco, nlev = nlev)
##       })




