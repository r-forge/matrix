facshuffle = function(sslm, facs)       # unexported utility
{
    if (!length(sslm[[2]])) return(facs)
    s1 = sslm[[1]]
    s2 = sslm[[2]]
    lens = diff(s1@Gp)
    ff = vector("list", length(facs))
    for (i in seq(along = lens)) {
        sq = seq(lens[i])
        perm = 1 + s2[sq]
        s2 = s2[-sq] - lens[i]
        fi = facs[[i]]
        fip = factor(perm[as.integer(fi)])
        levels(fip)[perm] = levels(fi)
        ff[[i]] = fip
    }
    ff
}

ssclme =
    function(formula, data, random, weights, subset, na.action, method)
{
    random = as(random, "list")
#        lapply(random,
#               function(x) if(inherits(x, "formula")) pdLogChol(x) else x)
    method = if (missing(method)) "REML" else {
        match.arg(method, c("REML", "ML"))
    }
    mCall <- match.call(expand.dots = FALSE)
    mCall[[1]] <- as.name("model.frame")
    names(mCall)[2] <- "formula"
    mCall$random <- mCall$correlation <- mCall$method <-
        mCall$control <- NULL
    form <- formula
    form[[3]] <- (~a+b)[[2]]
    form[[3]][[2]] <- formula[[3]]
    form[[3]][[3]] <-
        as.formula((parse(text=paste("~",
                          paste(names(random),
                                collapse = "+")))[[1]]))[[2]]
    for (pdm in random) {
        tmp <- form
        tmp[[3]] <- (~a+b)[[2]]
        tmp[[3]][[2]] <- form[[3]]
        tmp[[3]][[3]] <- formula(pdm)[[2]]
        form <- tmp
    }
    environment(form) <- environment(formula)
    mCall$formula <- form
    mCall$drop.unused.levels <- TRUE
    data <- eval(mCall, parent.frame())
    grps <- lapply(names(random), function(x) eval(as.name(x), envir = data))
    names(grps) <- names(random)
    mmats <- c(lapply(random,
                      function(x) model.matrix(formula(x), data = data)),
               list(.Xy = cbind(model.matrix(formula, data = data),
                    .response = model.response(data))))
    obj = .Call("ssclme_create", grps, unlist(lapply(mmats, ncol)),
                as.integer(2e5), PACKAGE = "Matrix")
    facs = facshuffle(obj, grps)
    obj = obj[[1]]
    .Call("ssclme_update_mm", obj, facs, mmats, PACKAGE="Matrix")
    .Call("ssclme_initial", obj, PACKAGE="Matrix")
    .Call("ssclme_factor", obj, PACKAGE = "Matrix")
    obj
}

setMethod("deviance", signature(object = "ssclme"),
          function(object, REML = FALSE, ...) {
              .Call("ssclme_factor", object, PACKAGE = "Matrix")
              object@deviance[ifelse(REML, 2, 1)]
          })

setMethod("coef", signature(object = "ssclme"),
          function(object, ...) {
              .Call("ssclme_coef", object, PACKAGE = "Matrix")
          })

setMethod("ranef", signature(object = "ssclme"),
          function(object, ...) {
              .Call("ssclme_ranef", object, PACKAGE = "Matrix")
          })

setMethod("fixef", signature(object = "ssclme"),
          function(object, ...) {
              .Call("ssclme_fixef", object, PACKAGE = "Matrix")
          })
