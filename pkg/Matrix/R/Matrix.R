## METHODS FOR CLASS: Matrix (virtual)
## mother class containing all matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~ CONSTRUCTORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Matrix <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE,
                   dimnames = NULL, sparse = NULL,
                   doDiag = TRUE, forceCheck = FALSE)
{
    i.M <- i.sM <- i.dM <- i.sV <- i.m <- FALSE
    mnrow <- missing(nrow)
    mncol <- missing(ncol)
    if(isS4(data)) {
        cld <- getClassDef(class(data))
        i.M <- extends(cld, "Matrix")
        if(i.M) {
            i.sM <- extends(cld, "sparseMatrix")
            i.dM <- i.sM && extends(cld, "diagonalMatrix")
        } else if(extends(cld, "sparseVector")) {
            ## need to transmit missingness to 'spV2M'
            call. <- quote(spV2M(x = data, nrow =, ncol =, byrow = byrow))
            if(!mnrow)
                call.[[3L]] <- quote(nrow)
            if(!mncol)
                call.[[4L]] <- quote(ncol)
            data <- eval(call.)
            i.M <- i.sM <- i.sV <- forceCheck <- TRUE
        }
    } else {
        i.m <- is.matrix(data)
    }
    if(!i.M) {
        ## validate non-Matrix 'data', throwing type errors _early_
        if(is.object(data)) {
            if(i.m)
                class(data) <- NULL # retaining 'dim'
            else
                data <- as.vector(data)
        }
        mode. <- mode(data)
        kind <- switch(mode., numeric = "d", logical = "l",
                       stop("invalid 'data'"))
    }
    if(i.M || i.m) {
        ## 'data' is a Matrix or a numeric or logical matrix
        ## without a 'class' attribute
        if(!i.sV && !(mnrow && mncol && missing(byrow)))
            warning("'nrow', 'ncol', 'byrow' disregarded for [mM]atrix 'data'")
        if(!is.null(dimnames))
            dimnames(data) <- dimnames
        if(is.null(sparse))
            sparse <- sparseDefault(data)
        if(i.M) {
            ## return early in these cases:
            if(i.dM)
                ## !doDiag has been documented to result in a coercion to
                ## symmetricMatrix; we must use diag2*() below because the
                ## "usual" as(<diagonalMatrix>, "(Csparse|unpacked)Matrix")
                ## inherits from triangularMatrix, _not_ symmetricMatrix
                return(if(doDiag)
                           data
                       else if(sparse)
                           .diag2sparse(data, "s", "C", "U")
                       else .diag2dense(data, "s", FALSE, "U"))
            if(!forceCheck)
                return(if(i.sM == sparse)
                           data
                       else if(sparse)
                           as(data, "CsparseMatrix")
                       else as(data, "unpackedMatrix"))
        }
    } else {
        ## 'data' is a numeric or logical vector or non-matrix array
        ## without a 'class' attribute
        if(length(data) == 1L && !is.na(data) && data == 0 &&
           (is.null(sparse) || sparse)) {
            ## Matrix(0, ...): sparseMatrix unless sparse=FALSE
            ## MJ: we should _try_ to behave as R's do_matrix()
            ##     in the edge cases ... integer overflow is "OK"
            ##     since anyNA(Dim) is caught by validity methods
            if(mnrow == mncol) {
                nrow <- as.integer(nrow)
                ncol <- as.integer(ncol)
            } else if(mnrow) {
                ncol <- as.integer(ncol)
                if(ncol == 0L)
                    stop("data is too long")
                nrow <- as.integer(ceiling(1 / ncol))
            } else {
                nrow <- as.integer(nrow)
                if(nrow == 0L)
                    stop("data is too long")
                ncol <- as.integer(ceiling(1 / nrow))
            }
            square <- nrow == ncol
            if(is.null(dimnames))
                dimnames <- list(NULL, NULL)
            if(square && doDiag)
                return(new(paste0(kind, "diMatrix"),
                           Dim = c(nrow, ncol),
                           Dimnames = dimnames,
                           x = vector(mode., nrow)))
            data <- new(paste0(kind, if(square) "s" else "g", "CMatrix"),
                        Dim = c(nrow, ncol),
                        Dimnames = dimnames,
                        p = integer(ncol + 1))
            i.M <- i.sM <- sparse <- TRUE
        } else {
            ## usual case: vector|array->matrix
            data <- .External(Mmatrix,
                              data, nrow, ncol, byrow, dimnames, mnrow, mncol)
            if(is.null(sparse))
                sparse <- sparseDefault(data)
            i.m <- TRUE
        }
    }

    ## 'data' is a Matrix (but _not_ a diagonalMatrix) or a
    ## numeric or logical matrix without a 'class' attribute
    if(doDiag && isDiagonal(data))
        ## as(<[mM]atrix>, "diagonalMatrix") uses check = TRUE (a waste)
        return(forceDiagonal(data))
    if(i.m || i.sM != sparse) {
        data <- as(data, if(sparse) "CsparseMatrix" else "unpackedMatrix")
        if(i.m)
            ## as(<matrix>, "CsparseMatrix"), as(<matrix>, "unpackedMatrix")
            ## already check for symmetric, triangular structure
            return(data)
    }
    if(!is(data, "generalMatrix"))
        data
    else if(isSymmetric(data))
        forceSymmetric(data)
    else if(!(it <- isTriangular(data)))
        data
    else if(attr(it, "kind") == "U")
        triu(data)
    else tril(data)
}


## ~~~~ METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("diff", signature(x = "Matrix"),
          ## Mostly cut and paste of 'base::diff.default' :
          function(x, lag = 1L, differences = 1L, ...) {
              if(length(lag) != 1L || length(differences) > 1L ||
                  lag < 1L || differences < 1L)
                  stop("'lag' and 'differences' must be integers >= 1")
              if(lag * differences >= x@Dim[1L])
                  return(x[0L])
              i1 <- -seq_len(lag)
              for(i in seq_len(differences)) {
                  m <- x@Dim[1L]
                  x <- x[i1, , drop = FALSE] -
                      x[-m:-(m - lag + 1L), , drop = FALSE]
              }
              x
          })

if(FALSE) { ## still does not work for c(1, Matrix(2))
## For the same reason (and just in case) also do both S3 and S4 here:
c.Matrix <- function(...) unlist(lapply(list(...), as.vector))
## NB: Must use   signature  '(x, ..., recursive = FALSE)' :
setMethod("c", "Matrix", function(x, ..., recursive) c.Matrix(x, ...))
## The above is not sufficient for  c(NA, 3:2, <Matrix>, <matrix>)
setMethod("c", "numMatrixLike", function(x, ..., recursive) c.Matrix(x, ...))
}# not yet

## We want to use all.equal.numeric() *and* make sure that uses
## not just base::as.vector but the generic with our methods:
all.equal_num <- base::all.equal.numeric
##               ^^^^^^<R>/src/library/base/R/all.equal.R
environment(all.equal_num) <- environment() # our namespace

all.equal_Mat <- function(target, current, check.attributes = TRUE,
                          factorsCheck = FALSE, ...)
{
    msg <- attr.all_Mat(target, current, check.attributes=check.attributes,
                        factorsCheck=factorsCheck, ...)
    if(is.list(msg)) msg[[1]]
    else .a.e.comb(msg,
                   all.equal_num(as.vector(target), as.vector(current),
                                 check.attributes=check.attributes, ...))
}

## The all.equal() methods for dense matrices (and fallback):
setMethod("all.equal", c(target = "Matrix", current = "Matrix"),
          all.equal_Mat)
setMethod("all.equal", c(target = "Matrix", current = "ANY"),
          all.equal_Mat)
setMethod("all.equal", c(target = "ANY", current = "Matrix"),
          all.equal_Mat)
rm(all.equal_Mat)
## -> ./sparseMatrix.R, ./sparseVector.R  have specific methods
