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

#setClass("summary.glmm", representation(method="character",
#                                        family="character",
#                                        link="character"),
#         contains="summary.lme")

setClass("lme", representation(call = "call",
                               facs = "list",
                               mmats = "list",
                               model = "data.frame",
                               REML = "logical",
                               rep = "ssclme"))

## This is needed for the family slot of glmmStruct
setOldClass("family")

## Structure for fitting glmm classes
setClass("glmm",
         representation(family="family", # The glm family
                        origy="numeric",
                        n="numeric",
                        prior.weights="numeric",
                        init.weights="numeric",
                        init.y="numeric",
                        method="character"),
         contains = "lme")

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

setClass("VarCorr",
         representation(scale="numeric",
                        reSumry="list",
                        useScale="logical"),
         prototype = list(scale = 1.0, useScale = TRUE))

## Deepayan experimenting with possible groupedData structures:

setClass("groupedData",
         representation(data = "data.frame",
                        formula = "formula",
                        outer = "formula",
                        inner = "formula",
                        labels = "list",
                        units = "list"))

