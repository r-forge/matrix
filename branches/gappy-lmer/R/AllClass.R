## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("lmList.confint", contains = "array")

## -------------------- lmer-related Classes --------------------------------

setOldClass("data.frame")
setOldClass("family")
setOldClass("logLik")
setOldClass("terms")

### FIXME: the nlmer class has 'mu' and 'resid' slots.  Should those
### be part of the mer class?  Should eta, the linear predictor, be a
### slot in the mer class?

setClass("mer", ## Slots common to all three types of mixed models
	 representation(## original data
                        frame = "data.frame", # model frame (or empty frame)
                        call = "call",      # matched call
                        terms = "terms",    # terms for fixed-effects
			flist = "list",     # list of grouping factors
                        X = "matrix",       # fixed effects model matrix
                                            # (may have 0 rows in lmer)
			Zt = "dgCMatrix",   # sparse form of Z'
			weights = "numeric", #  length 0 -> constant wts
                        offset = "numeric", # length 0 -> no offset
                        y = "numeric",      # response vector
			cnames = "list",    # row/column names of els of ST
			Gp = "integer",     # pointers to row groups of Zt
                        dims = "integer",   # dimensions and indicators
                        ## slots that vary during optimization
			ST = "list", # list of TSST' rep of rel. var. mats
			L = "CHMfactor",    # Cholesky factor of P(V'V + I)P'
			Vt = "dgCMatrix",   # sparse form of V'=(ZTS)'
			deviance = "numeric", # ML and REML deviance and components
			fixef = "numeric",  # fixed effects (length p)
			ranef = "numeric",  # random effects (length q)
                        uvec = "numeric",   # orthogonal random effects (q)
                        eta = "numeric",    # linear predictor
                        mu = "numeric",     # fitted values at current beta and b
                        resid = "numeric",  # raw residuals at current beta and b
                        "VIRTUAL"),
         validity = function(object) .Call(mer_validate, object))

setClass("lmer", ## linear mixed models
	 representation(## original data
                        ZtXy = "matrix",  # dense form of Z'[X:y]
                        XytXy = "matrix", # dense form of [X:y]'[X:y]
                        ## slots that vary during optimization
                        RVXy = "matrix",  # dense sol. to L RVXy = S T'ZtXy
                        RXy = "matrix"),  # Cholesky factor of downdated XytXy
         contains = "mer")

setClass("glmer", ## generalized linear mixed models
	 representation(## original data
                        env = "environment", # evaluation env for family
                        famName = "character"), # name of GLM family and link
         contains = "mer")

setClass("nlmer", ## nonlinear mixed models
	 representation(## original data
                        env = "environment", # evaluation environment for model
                        model = "call",    # nonlinear model
                        pnames = "character", # parameter names for nonlinear model
                        ## slots that vary during optimization
                        Mt = "dgCMatrix"), # transpose of gradient matrix d mu/d u
         contains = "mer")

setClass("summary.mer",                 # Additional slots in a summary object
	 representation(           
			methTitle = "character",
			logLik= "logLik",
			ngrps = "integer",
			sigma = "numeric", # scale, non-negative number
			coefs = "matrix",
			vcov = "dpoMatrix",
			REmat = "matrix",
			AICtab= "data.frame",
                        "VIRTUAL"))

setClass("summary.lmer", # the "lmer" result ``enhanced'' :
	 contains = c("lmer", "summary.mer"))

setClass("summary.glmer",          # the "glmer" result ``enhanced'' :
	 contains = c("glmer", "summary.mer"))

setClass("summary.nlmer",          # the "glmer" result ``enhanced'' :
	 contains = c("nlmer", "summary.mer"))

setClass("ranef.mer", contains = "list")

setClass("coef.lmer", contains = "list")

setClass("pedigree", representation =
	 list(sire = "integer", dam = "integer", label = "character"),
	 validity = function(object) {
	     n <- length(sire <- object@sire)
	     if (length(dam <- object@dam) != n)
		 return("sire and dam slots must be the same length")
	     if (length(object@label) != n)
		 return("'label' slot must have the same length as 'sire' and 'dam'")
	     if(n == 0) return(TRUE)
	     animal <- 1:n
	     snmiss <- !is.na(sire)
	     dnmiss <- !is.na(dam)
	     if (any(sire[snmiss] >= animal[snmiss]) ||
		 any(dam[dnmiss] >= animal[dnmiss]))
		 return("the sire and dam must precede the offspring")
             if (any(sire[snmiss] < 1 | sire[snmiss] > n) |
                 any(dam[dnmiss] < 1 | dam[dnmiss] > n))
                 return(paste("Non-missing sire or dam must be in [1,",
                              n, "]", sep = ''))
	     TRUE
	 })
