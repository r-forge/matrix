# Ensure that the stats4 and lattice packages are available
.onLoad <- function(lib, pkg)
{
    if ("package:nlme" %in% search()) {
        stop(paste("Package lme4 conflicts with package nlme.\n",
                   "To attach lme4 you must restart R without package nlme."))
    }
    require("methods", quietly = TRUE)
    require("stats4", quietly = TRUE)
    require("lattice", quietly = TRUE)
    require("Matrix", quietly = TRUE)
    cat(paste(" This package is in development.  For production work use\n",
              "lme from package nlme or glmmPQL from package MASS.\n"))
}

setClass("summary.pdMat", representation(cor = "corrmatrix",
                                         structName = "character",
                                         noCorrelation = "logical",
                                         formula = "formula"),
         prototype=list(structName="", formula=formula(NULL)))


setClass("lme", representation(call = "call",
                               facs = "list",
                               x = "list",
                               model = "data.frame",
                               REML = "logical",
                               rep = "ssclme",
                               fitted = "numeric",
                               residuals = "numeric"))

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("VarCorr",
         representation(scale="numeric",
                        reSumry="list",
                        useScale="logical"),
         prototype = list(scale = 1.0, useScale = TRUE))

setClass("groupedData",
         representation(data = "data.frame",
                        formula = "formula",
                        outer = "formula",
                        inner = "formula",
                        labels = "list",
                        units = "list"))

setClass("summary.ssclme",
         representation(coefficients="matrix",
                        scale="numeric",
                        denomDF="integer",
                        REML="logical",
                        ngrps="integer",
                        nobs="integer",
                        corFixed="corrmatrix",
                        VarCorr="VarCorr",
                        useScale="logical",
                        showCorrelation="logical"
                        ))

setOldClass("logLik")

setClass("summary.lme",
         representation(call = "call",
                        logLik = "logLik",
                        re = "summary.ssclme",
                        residuals = "numeric"))

setOldClass("family")

