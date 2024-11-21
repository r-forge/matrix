validDim <-
function(d)
    .Call(R_Dim_validate, d)

validDimGetsValue <-
function(value, dim, len = .Call(R_Dim_prod, as.integer(dim))) {
    if (mode(value) != "numeric")
        gettextf("assigned dimensions are not of type \"%s\" or \"%s\"",
                 "integer", "double")
    else if (length(value) != 2L)
        gettextf("assigned dimensions do not have length %d",
                 2L)
    else if (anyNA(value))
        gettext("assigned dimensions are NA")
    else if (min(value) <= -1L)
        gettext("assigned dimensions are negative")
    else if (is.double(value) &&
             max(value) - 1 >= .Machine[["integer.max"]])
        gettextf("assigned dimensions exceed %s",
                 "2^31-1")
    else if (!identical(p <- .Call(R_Dim_prod, as.integer(value)), len))
        gettextf("product of assigned dimensions [%.0f] is not equal to object length [%.0f]",
                 p, len)
    else as.integer(value)
}

validLengthGetsValue <-
function(value, max) {
    if (mode(value) != "numeric")
        gettextf("assigned length is not of type \"%s\" or \"%s\"",
                 "integer", "double")
    else if (length(value) != 1L)
        gettextf("assigned length does not have length %d",
                 2L)
    else if (is.na(value))
        gettext("assigned length is NA")
    else if (value <= -1L)
        gettext("assigned length is negative")
    else if (is.double(value) &&
             value - 1 >= max)
        gettextf("assigned length exceeds %.0f",
                 max)
    else (if (is.double(value)) trunc else as.integer)(value)
}

validDN <-
function(dn, d)
    .Call(R_DimNames_validate, dn, d)

fixupDN <-
function(dn)
    .Call(R_DimNames_fixup, dn)

fixupDN.if.valid <-
function(dn, d) {
    if (is.character(s <- validDim(d)) ||
        is.character(s <- validDN(dn, d)))
        stop(s, domain = NA)
    fixupDN(dn)
}

symDN <-
function(dn)
    .Call(R_symDN, dn)

symmetrizeDN <-
function(x) {
    if (isS4(x)) # assuming is(x, "Matrix")
        `dimnames<-`(x, symDN(x@Dimnames))
    else if (!is.null(dn <- dimnames(x))) # assuming list of length 2
        `dimnames<-`(x, symDN(dn))
    else x
}

isSymmetricDN <-
function(dn)
    .Call(R_DimNames_is_symmetric, dn)

is.null.DN <-
function(dn) {
    if (is.null(dn))
        return(TRUE)
    if (!is.null(names(dn)))
        names(dn) <- NULL
    ch0 <- character(0L)
    identical(dn, list(NULL, NULL)) ||
    identical(dn, list( ch0, NULL)) ||
    identical(dn, list(NULL,  ch0)) ||
    identical(dn, list( ch0,  ch0))
}


## METHODS FOR GENERIC: dim
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim", c(x = "Matrix"),
          function(x) x@Dim)

setMethod("dim", c(x = "MatrixFactorization"),
          function(x) x@Dim)

setMethod("dim", c(x = "sparseVector"),
          function(x) NULL)


## METHODS FOR GENERIC: dim<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim<-", c(x = "denseMatrix", value = "numeric"),
          function(x, value) {
              d <- x@Dim
              value <- validDimGetsValue(value, d)
              if (is.character(value))
                 stop(value, domain = NA)
              if (all(value == d))
                  return(x)
              r <- new(paste0(.M.kind(x), "geMatrix"))
              r@Dim <- value
              r@x <- .M2v(x)
              r
          })

setMethod("dim<-", c(x = "sparseMatrix", value = "numeric"),
          function(x, value) {
              d <- x@Dim
              value <- validDimGetsValue(value, d)
              if (is.character(value))
                 stop(value, domain = NA)
              if (all(value == d))
                  return(x)
              cl <- switch(.M.repr(x),
                           "C" = ".gC",
                           "R" = ".gR",
                           "T" = ".gT",
                           "d" = ".gC",
                           "i" = if (x@margin == 1L) ".gR" else ".gC")
              .V2sparse(.M2V(x), cl, nrow = value[1L], ncol = value[2L])
          })

setMethod("dim<-", c(x = "sparseVector", value = "numeric"),
          function(x, value) {
              value <- validDimGetsValue(value, len = length(x))
              if (is.character(value))
                 stop(value, domain = NA)
              .V2sparse(x, ".gC", nrow = value[1L], ncol = value[2L])
          })


## METHODS FOR GENERIC: length
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("Matrix", "MatrixFactorization"))
setMethod("length", c(x = .cl),
          function(x)
              .Call(R_Dim_prod, x@Dim))

setMethod("length", c(x = "sparseVector"),
          function(x) {
              r <- x@length
              if (is.integer(r))
                  r
              else if (r - 1 < .Machine[["integer.max"]])
                  as.integer(r)
              else trunc(r)
          })

rm(.cl)


## METHODS FOR GENERIC: length<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("length<-", c(x = "denseMatrix", value = "numeric"),
          function(x, value) {
              value <- validLengthGetsValue(value, 0x1p+52)
              if (is.character(value))
                  stop(value, domain = NA)
              mn <- length(x)
              if (value == mn)
                  x
              else `length<-`(.M2v(x), value)
          })

setMethod("length<-", c(x = "sparseMatrix", value = "numeric"),
          function(x, value) {
              value <- validLengthGetsValue(value, 0x1p+53)
              if (is.character(value))
                  stop(value, domain = NA)
              mn <- length(x)
              if (value == mn)
                  x
              else `length<-`(.M2V(x), value)
          })

setMethod("length<-", c(x = "sparseVector", value = "numeric"),
          function(x, value) {
              value <- validLengthGetsValue(value, 0x1p+53)
              if (is.character(value))
                  stop(value, domain = NA)
              mn <- length(x)
              if (value == mn)
                  return(x)
              i <- x@i
              kind <- .M.kind(x)
              y <- new(paste0(kind, "sparseVector"))
              y@length <- value
              if (mn < value) {
                  y@i <- c(i, seq.int(from = mn + 1L, to = value))
                  if (kind != "n")
                      y@x <- c(x@x, rep.int(x@x[NA_integer_], value - mn))
              } else if (length(i) > 0L) {
                  if (i[length(i)] - 1L < value) {
                  y@i <- i
                  if (kind != "n")
                      y@x <- x@x
                  } else {
                  k <- which(i - 1L < value)
                  y@i <- i[k]
                  if (kind != "n")
                      y@x <- x@x[k]
                  }
              }
              y
          })


## METHODS FOR GENERIC: dimnames
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("Matrix", "MatrixFactorization"))
setMethod("dimnames", c(x = .cl),
          function(x) x@Dimnames)

setMethod("dimnames", c(x = "symmetricMatrix"),
          function(x) symDN(x@Dimnames))

setMethod("dimnames", c(x = "sparseVector"),
          function(x) NULL)

rm(.cl)


## METHODS FOR GENERIC: dimnames<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("Matrix", "MatrixFactorization")) {
setMethod("dimnames<-", c(x = .cl, value = "NULL"),
          function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", c(x = .cl, value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })
}

for (.cl in c("generalMatrix", "symmetricMatrix")) {
setMethod("dimnames<-", c(x = .cl, value = "NULL"),
          function(x, value) {
              if (length(x@factors))
                  x@factors <- list()
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", c(x = .cl, value = "list"),
          function(x, value) {
              if (length(x@factors))
                  x@factors <- list()
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })
}

rm(.cl)


## METHODS FOR GENERIC: drop
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("drop", c(x = "Matrix"),
          function(x)
              if (any(x@Dim == 1L)) drop(.M2m(x)) else x)
