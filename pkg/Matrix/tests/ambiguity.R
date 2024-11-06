library(Matrix)
requireNamespace("SparseM")
requireNamespace("graph")

ns <- asNamespace("Matrix")
generics <- getGenerics(ns)
l <- sapply(generics, testInheritedMethods, simplify = FALSE)

(namb <- vapply(l, function(x) length(x@target), 0L))
(todo <- names(namb)[namb > 0L])

## Make sure that we are removing from and never adding to this
## list of generic functions with ambiguous method dispatch :
todo. <- c("%%", "%/%", "&", "*", "+", "-", "/",
           "Arith", "Compare", "Logic", "Ops",
           "[<-", "^", "crossprod", "tcrossprod")
stopifnot(identical(todo, todo.))
