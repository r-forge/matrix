library(lme4)
data(guImmun)
# system.time(fm <-
#             sparse GLMM(immun ~ kid2p + mom25p + ord + ethn +
#                         momEd + husEd + momWork + rural + pcInd81,
#                         data = guImmun, family = binomial,
#                         random = list(mom = ~ 1,comm = ~1))
# summary(fm)
q("no")
