bbCrosstab <- function(flist) {
    ind <- function(i,j) ((i-1) * i)/2 + j # index in rowwise lower triangle
    flist <- lapply(as.list(flist), function (x) as(x, "factor")[drop = TRUE])
    nfac <- length(flist)
    if (nfac < 1) return(new("bbSparseSy", x = list(), uplo = 'L'))
    nlev <- unlist(lapply(flist, function(x) length(levels(x))))
    ord <- rev(order(nlev))
    if (any(ord != seq(along = ord))) {
        nlev <- nlev[ord]
        flist <- flist[ord]
    }
    nobs <- length(flist[[1]])
    if (any(lapply(flist, length) != nobs))
        stop("all factors in flist must have the same length")
    ones <- rep(1, nobs)
    lst <- vector("list", choose(nfac + 1, 2))
    if (is.null(nms <- names(flist))) nms <- paste("V", 1:nfac, sep = "")
    nmat <- outer(nms, nms, FUN = "paste", sep = ":")
    nms <- vector("character", length(lst))
    zb <- lapply(flist, function(x) as.integer(x) - 1:1) # zero-based indices
    for (j in 1:nfac) {
        for (i in j:nfac) {
            lst[[ind(i,j)]] <- as(as(new("tripletMatrix",
                                       i = zb[[i]], j = zb[[j]],
                                       x = ones, Dim = nlev[c(i,j)]),
                                   "cscMatrix"),
                                "cscBlocked")
            nms[ind(i,j)] <- nmat[i, j]
        }
    }
    names(lst) <- nms
    new("bbCrosstab", x = lst, uplo = "L")
}

setMethod("isNested", signature(object = "bbCrosstab"),
          function(object, ...)
          all(unlist(lapply(object@x, function(x) any(diff(x@p) > 1)))))

bbCrossNfac <- function(x) {            # determine the number of factors
    x <- as(x, "bbCrosstab")
    as.integer((-1 + sqrt(1 + 8*length(x@x)))/2)
}

bbCrossProject <- function(x) {         # Project the first column onto the others
    ind <- function(i,j) ((i-1) * i)/2 + j # index in compressed lower triangle
    x <- as(x, "bbCrosstab")
    nf <- bbCrossNfac(x)
    if (nf < 2) stop("Number of factors must be > 1 to project")
    ans <- new("bbCrosstab", x = x@x[-(1:nf)], uplo = "L")
    for (i in 2:nf) {
        ans@x[[ind(i-1,i-1)]] <-
            as(tcrossprod(as(x@x[[i]], "cscMatrix")), "cscBlocked") 
    }
    ans
}

bbCrossPerms <- function(x) {
    ind <- function(i,j) ((i-1) * i)/2 + j
    needsPerm <- function(x) any(diff(as(x, "cscBlocked")@p) != 1)
    x <- as(x, "bbCrosstab")
    nf <- bbCrossNfac(x)
    perm <- vector("list", nf)
#    Lmat <- new("bbSparseTr", uplo = 'L', diag = 'U',
#                Linv = vector("list", nf), x = x@x)
    Linv <- vector("list", nf)
    for (j in 1:nf) {
        db <- x@x[[ind(j,j)]]             # diagonal block
        if (needsPerm(db)) {
            res <- .Call("sscMatrix_ldl_symbolic",
                         as(as(db, "cscMatrix"), "sscMatrix"), TRUE)
        } else {
            ncol <- length(db@p) - 1
            res <- list(rep(-(1:1), ncol), NULL, seq(ncol) - 1:1)
        }
        perm[[j]] <- res[[3]]
        Linv[[j]] <- as(.Call("Parent_inverse", res[[1]], TRUE),
                        "cscBlocked")    
    }
#    Lmat$Linv <- Linv
#    list(perm = perm, Lmat = Lmat)
    list(perm = perm, Linv = Linv) 
}

