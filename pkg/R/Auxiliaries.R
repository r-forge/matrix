#### "Namespace private" Auxiliaries  such as method functions
#### (called from more than one place --> need to be defined early)

## For %*% (M = Matrix; v = vector (double or integer {complex maybe?}):
.M.v <- function(x, y) callGeneric(x, as.matrix(y))
.v.M <- function(x, y) callGeneric(rbind(x), y)

.has.DN <- ## has non-trivial Dimnames slot?
    function(x) !identical(list(NULL,NULL), x@Dimnames)

## chol() via "dpoMatrix"
cholMat <- function(x, pivot, LINPACK) {
    px <- as(x, "dpoMatrix")
    if (isTRUE(validObject(px, test=TRUE))) chol(px)
    else stop("'x' is not positive definite -- chol() undefined.")
}

rowCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(da[1] != db[1])
	stop(gettextf("Matrices must have same number of rows in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE)
    ## return the common nrow()
    da[1]
}

colCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(da[2] != db[2])
	stop(gettextf("Matrices must have same number of columns in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE)
    ## return the common ncol()
    da[2]
}


prTriang <- function(x, digits = getOption("digits"),
                     justify = "none", right = TRUE)
{
    ## modeled along stats:::print.dist
    diag <- TRUE
    upper <- x@uplo == "U"

    m <- as(x, "matrix")
    cf <- format(m, digits = digits, justify = justify)
    if(upper)
        cf[row(cf) > col(cf)] <- "."
    else
        cf[row(cf) < col(cf)] <- "."
    print(cf, quote = FALSE, right = right)
    invisible(x)
}

prMatrix <- function(x, digits = getOption("digits")) {
    d <- dim(x)
    cl <- class(x)
    cat(sprintf('%d x %d Matrix of class "%s"\n', d[1], d[2], cl))
    maxp <- getOption("max.print")
    if(prod(d) <= maxp) {
        if(is(x, "triangularMatrix"))
            prTriang(x, digits = digits)
        else
            print(as(x, "matrix"), digits = digits)
    }
    else { ## d[1] > maxp / d[2] >= nr :
        m <- as(x, "matrix")
	nr <- maxp %/% d[2]
	n2 <- ceiling(nr / 2)
	print(head(m, max(1, n2)))
	cat("\n ..........\n\n")
	print(tail(m, max(1, nr - n2)))
    }
    ## DEBUG: cat("str(.):\n") ; str(x)
    invisible(x)# as print() S3 methods do
}

## For sparseness handling
non0ind <- function(x) {
    if(is.numeric(x))
        return(if((n <- length(x))) (0:(n-1))[x != 0] else integer(0))

    ## else return a (i,j) matrix of non-zero-indices

    stopifnot(is(x, "sparseMatrix"))
    if(is(x, "gTMatrix"))
        stop("'x' must be column- or row-compressed  'sparseMatrix'")
    isCol <- function(M) any("i" == slotNames(M))
    .Call("compressed_non_0_ij", x, isCol(x), PACKAGE = "Matrix")
}

### These are currently tests in ../tests/dgTMatrix.R !!!
uniq <- function(x) {
    if(is(x, "gTMatrix")) {
        ## Purpose: produce a *unique* triplet representation:
        ##		by having (i,j) sorted and unique
        ## -----------------------------------------------------------
        ## The following is *not* efficient {but easy to program}:
        if(is(x, "dgTMatrix")) as(as(x, "dgCMatrix"), "dgTMatrix")
        else if(is(x, "lgTMatrix")) as(as(x, "lgCMatrix"), "lgTMatrix")
        else stop("not implemented for class", class(x))

    } else x      # not 'gT' ; i.e. "uniquely" represented in any case
}

if(FALSE) ## try an "efficient" version
uniq_gT <- function(x)
{
    ## Purpose: produce a *unique* triplet representation:
    ##		by having (i,j) sorted and unique
    ## ----------------------------------------------------------------------
    ## Arguments: a "gT" Matrix
    stopifnot(is(x, "gTMatrix"))
    if((n <- length(x@i)) == 0) return(x)
    ii <- order(x@i, x@j)
    if(any(ii != 1:n)) {
        x@i <- x@i[ii]
        x@j <- x@j[ii]
        x@x <- x@x[ii]
    }
    ij <- x@i + nrow(x) * x@j
    if(any(dup <- duplicated(ij))) {

    }
    ### We should use a .Call() based utility for this!

}

