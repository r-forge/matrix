validDim <- function(dim)
    .Call(R_Dim_validate, dim)

validDN <- function(dn, dim)
    .Call(R_DimNames_validate, dn, dim)

fixupDN <- function(dn)
    .Call(R_DimNames_fixup, dn)

fixupDN.if.valid <- function(dn, dim) {
    if(is.character(s <- validDim(dim)) || is.character(s <- validDN(dn, dim)))
        stop(s)
    fixupDN(dn)
}

symmDN <- function(dn)
    .Call(R_symmDN, dn)

symmetrizeDN <- function(x) {
    if(isS4(x)) # assuming is(x, "Matrix")
        `dimnames<-`(x, symmDN(x@Dimnames))
    else if(!is.null(dn <- dimnames(x))) # assuming list of length 2
        `dimnames<-`(x, symmDN(dn))
    else x
}

isSymmetricDN <- function(dn)
    .Call(R_DimNames_is_symmetric, dn)

is.null.DN <- function(dn) {
    if(is.null(dn))
        return(TRUE)
    if(!is.null(names(dn)))
        names(dn) <- NULL
    ch0 <- character(0L)
    identical(dn, list(NULL, NULL)) ||
    identical(dn, list( ch0, NULL)) ||
    identical(dn, list(NULL,  ch0)) ||
    identical(dn, list( ch0,  ch0))
}


## METHODS FOR GENERIC: dim
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim", signature(x = "Matrix"),
          function(x) x@Dim)

setMethod("dim", signature(x = "MatrixFactorization"),
          function(x) x@Dim)


## METHODS FOR GENERIC: dim<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dim<-", signature(x = "denseMatrix"),
          function(x, value) {
              if(is.double(value)) {
                  if(any(value >= .Machine$integer.max + 1))
                      stop("dimensions cannot exceed 2^31-1")
              } else if(!is.integer(value))
                  stop("dimensions must be numeric")
              if(length(value) != 2L)
                  stop("dimensions must have length 2")
              value <- as.integer(value)
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0L))
                  stop("dimensions cannot be negative")
              if(all(value == (d <- x@Dim)))
                  return(x)
              if((pv <- prod(value)) != (pd <- prod(d)))
                  stop(gettextf("assigned dimensions [product %.0f] do not match Matrix length [%.0f]",
                                pv, pd),
                       domain = NA)
              r <- .M2gen(x)
              r@Dim <- value
              if(length(r@factors))
                  r@factors <- list()
              r
          })

setMethod("dim<-", signature(x = "sparseMatrix"),
          function(x, value) {
              if(is.double(value)) {
                  if(any(value >= .Machine$integer.max + 1))
                      stop("dimensions cannot exceed 2^31-1")
              } else if(!is.integer(value))
                  stop("dimensions must be numeric")
              if(length(value) != 2L)
                  stop("dimensions must have length 2")
              value <- as.integer(value)
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0L))
                  stop("dimensions cannot be negative")
              if(all(value == (d <- x@Dim)))
                  return(x)
              if((pv <- prod(value)) != (pd <- prod(d)))
                  stop(gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                                pv, pd),
                       domain = NA)
              r <- spV2M(.M2V(x), nrow = value[1L], ncol = value[2L])
              switch(.M.repr(x), "C" = .M2C(r), "R" = .M2R(r), r)
          })

setMethod("dim<-", signature(x = "sparseVector"),
          function(x, value) {
              if(is.double(value)) {
                  if(any(value >= .Machine$integer.max + 1))
                      stop("dimensions cannot exceed 2^31-1")
              } else if(!is.integer(value))
                  stop("dimensions must be numeric")
              if(length(value) != 2L)
                  stop("dimensions must have length 2")
              value <- as.integer(value)
              if(anyNA(value))
                  stop("dimensions cannot contain NA")
              if(any(value < 0L))
                  stop("dimensions cannot be negative")
              if((pv <- prod(value)) != (pd <- x@length))
                  stop(gettextf("assigned dimensions [product %.0f] do not match object length [%.0f]",
                                pv, pd),
                       domain = NA)
              spV2M(x, nrow = value[1L], ncol = value[2L])
          })


## METHODS FOR GENERIC: length
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("length", "Matrix",
          function(x)
              if((r <- prod(x@Dim)) > .Machine$integer.max)
                  r
              else as.integer(r))

setMethod("length", "MatrixFactorization",
          function(x)
              if((r <- prod(x@Dim)) > .Machine$integer.max)
                  r
              else as.integer(r))

setMethod("length", "sparseVector",
          function(x)
              if(is.integer(r <- x@length) || r > .Machine$integer.max)
                  r
              else as.integer(r))


## METHODS FOR GENERIC: dimnames
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dimnames", signature(x = "Matrix"),
          function(x) x@Dimnames)

setMethod("dimnames", signature(x = "symmetricMatrix"),
          function(x) symmDN(x@Dimnames))

setMethod("dimnames", signature(x = "MatrixFactorization"),
          function(x) x@Dimnames)


## METHODS FOR GENERIC: dimnames<-
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("dimnames<-", signature(x = "Matrix", value = "NULL"),
          function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "NULL"),
          function(x, value) {
              if(length(x@factors))
                  x@factors <- list()
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "MatrixFactorization", value = "NULL"),
          function(x, value) {
              x@Dimnames <- list(NULL, NULL)
              x
          })

setMethod("dimnames<-", signature(x = "Matrix", value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "compMatrix", value = "list"),
          function(x, value) {
              if(length(x@factors))
                  x@factors <- list()
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })

setMethod("dimnames<-", signature(x = "MatrixFactorization", value = "list"),
          function(x, value) {
              x@Dimnames <- fixupDN.if.valid(value, x@Dim)
              x
          })


## METHODS FOR GENERIC: unname
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("unname", signature(obj = "Matrix"),
          function(obj, force = FALSE) {
              obj@Dimnames <- list(NULL, NULL)
              obj
          })

setMethod("unname", signature(obj = "MatrixFactorization"),
          function(obj, force = FALSE) {
              obj@Dimnames <- list(NULL, NULL)
              obj
          })


## METHODS FOR GENERIC: drop
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("drop", signature(x = "Matrix"),
          function(x) if(any(x@Dim == 1L)) drop(.M2m(x)) else x)
