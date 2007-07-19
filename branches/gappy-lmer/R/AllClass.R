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

setClass("mer", ## Slots common to all three types of mixed models
	 representation(## original data
                        frame = "data.frame", # model frame or empty frame
                        call = "call", # matched call to model-fitting function
                        terms = "terms",    # terms for fixed-effects
			flist = "list",     # list of grouping factors
                        X = "matrix", # model matrix for fixed effects (may have 0 rows in lmer only)
			Zt = "dgCMatrix",   # sparse form of Z'
			weights = "numeric", # can be of length 0 for constant wts
                        offset = "numeric", # can be length 0 (for no offset)
                        y = "numeric",      # response vector
			cnames = "list", # row/column names of matrices in ST
			Gp = "integer", # pointers to groups of rows in Zt
                        dims = "integer", # dimensions and indicators
                        ## slots that vary during optimization
			ST = "list", # list of TSST' rep of rel. var. mats
			L = "CHMfactor", # sparse Cholesky factor of V'V + I
			Vt = "dgCMatrix",   # sparse form of V'=(ZTS)'
			deviance = "numeric", # ML and REML deviance and components
			fixef = "numeric", # fixed effects coefficients (length p)
			ranef = "numeric", # random effects (length q)
                        uvec = "numeric",  # orthogonal random effects (length q)
                        "VIRTUAL"),
         validity = function(object) .Call(mer_validate, object))

setClass("lmer", ## linear mixed models
	 representation(## original data
                        ZtXy = "matrix",  # dense form of Z'[X:y]
                        XytXy = "matrix", # dense form of [X:y]'[X:y]
                        ## slots that vary during optimization
                        RXy = "matrix", # dense Cholesky factor of downdated XytXy
                        RVXy = "matrix"), # dense solution to L RVXy = S T'ZtXy
         contains = "mer",
         validity = function(object) .Call(lmer_validate, object))

setClass("glmer", ## generalized linear mixed models
	 representation(## original data
                        env = "environment", # evaluation env for family
                        famName = "character"), # name of GLM family and link
         contains = "mer",
         validity = function(object) .Call(glmer_validate, object))

setClass("nlmer", ## nonlinear mixed models
	 representation(## original data
                        env = "environment", # evaluation environment for model
                        model = "call",      # nonlinear model
                        pnames = "character", # parameter names for nonlinear model
                        ## slots that vary during optimization
                        mu = "numeric", # fitted values at current values of beta and b
                        Mt = "dgCMatrix"), # transpose of gradient matrix d mu/d u
         contains = "mer",
         validity = function(object) .Call(nlmer_validate, object))

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

setClass("ranef.lmer", contains = "list")

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
