library(lme4)
data(ScotsSec, package = "Matrix")
fm1 = lme4:::ssclme(attain ~ verbal, ScotsSec,
   list(primary = ~1, second = ~1))
#stopifnot(all.equal(c(unlist(tapply(ScotsSec$attain, ScotsSec$primary,
#                                    FUN = sum))), ZtX[1:148, 3]))
#stopifnot(all.equal(c(unlist(tapply(ScotsSec$attain, ScotsSec$second,
#                                    FUN = sum))), ZtX[-(1:148), 3]))
#stopifnot(all.equal(c(unlist(tapply(ScotsSec$verbal, ScotsSec$primary,
#                                    FUN = sum))), ZtX[1:148, 2]))
#stopifnot(all.equal(c(unlist(tapply(ScotsSec$verbal, ScotsSec$second,
#                                    FUN = sum))), ZtX[-(1:148), 2]))
str(fm1)
deviance(fm1, TRUE)
coef(fm1)
system.time(.Call("ssclme_EMstepsGets", fm1, 20, TRUE, PACKAGE = "Matrix"))
deviance(fm1, TRUE)
coef(fm1)
