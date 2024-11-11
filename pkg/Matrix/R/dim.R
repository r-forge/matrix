validDim <-
function(d)
    .Call(R_Dim_validate, d)

validDimGetsValue <-
function(value, mn) {
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
             max(value <- trunc(value)) - 1 >= .Machine[["integer.max"]])
        gettextf("assigned dimensions exceed %s",
                 "2^31-1")
    else if ((p <- prod(value)) != mn)
        gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                 p, as.double(mn))
    else TRUE
}

validLengthGetsValue <-
function(value, mn) {
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
             (value <- trunc(value)) - 1 >= 0x1p+52)
        gettextf("assigned length exceeds %.0f",
                 0x1p+52)
    else TRUE
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
              s <- validDimGetsValue(value, prod(d))
              if (is.character(s))
                 stop(s, domain = NA)
              value <- as.integer(value)
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
              s <- validDimGetsValue(value, prod(d))
              if (is.character(s))
                 stop(s, domain = NA)
              value <- as.integer(value)
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
              s <- validDimGetsValue(value, length(x))
              if (is.character(s))
                 stop(s, domain = NA)
              value <- as.integer(value)
              .V2sparse(x, ".gC", nrow = value[1L], ncol = value[2L])
          })


## METHODS FOR GENERIC: length
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("Matrix", "MatrixFactorization"))
setMethod("length", c(x = .cl),
          function(x)
              if ((r <- prod(x@Dim)) > .Machine[["integer.max"]])
                  r
              else as.integer(r))

setMethod("length", c(x = "sparseVector"),
          function(x)
              if (is.integer(r <- x@length))
                  r
              else if (r - 1 < .Machine[["integer.max"]])
                  as.integer(r)
              else trunc(r))

rm(.cl)


## METHODS FOR GENERIC: length<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (.cl in c("denseMatrix", "sparseMatrix", "sparseVector"))
setMethod("length<-", c(x = .cl, value = "numeric"),
          function(x, value) {
              mn <- length(x)
              s <- validLengthGetsValue(value, mn)
              if (is.character(s))
                  stop(s, domain = NA)
              if (is.double(value))
                  value <- trunc(value)
              if (value == mn)
                  x
              else x[seq_len(value)]
          })

rm(.cl)


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
