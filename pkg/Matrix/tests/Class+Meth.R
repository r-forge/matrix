library(Matrix)
source(system.file("test-tools.R", package = "Matrix"))# identical3(),
                                        # further  checkMatrix(), etc
if(interactive()) options(error = recover)
options(warn=1)# show as they happen
cat("doExtras:",doExtras,"\n")

setClass("myDGC", contains = "dgCMatrix")
(M <- new("myDGC", as(Matrix(c(-2:4, rep(0,9)), 4), "CsparseMatrix")))
stopifnot(exprs = {
    M[-4L, 2L] == 2:4
    MatrixClass(     "myDGC") == "dgCMatrix"
    MatrixClass( "dpoMatrix") == "dsyMatrix"
    MatrixClass( "dppMatrix") == "dspMatrix"
    MatrixClass( "corMatrix") == "dsyMatrix"
    MatrixClass( "copMatrix") == "dspMatrix"
    MatrixClass( "indMatrix") == "indMatrix"
    MatrixClass(   "pMatrix") == "indMatrix"
})

## [matrix-Bugs][6182] Coercion method doesn't work on child class
## Bugs item #6182, at 2015-09-01 17:49 by Vitalie Spinu
setClass("A", contains = "ngCMatrix")
ngc <- as(as(as(diag(3), "CsparseMatrix"), "generalMatrix"), "nMatrix")
validObject(dd <- as(ngc, "dMatrix")) # fine
A. <- as(ngc, "A")
stopifnot(identical(as(A., "dMatrix"), dd))
## as(.) coercion failed in Matrix <= 1.2.3
stopifnot(all( abs(A.)# failed too
              == diag(3)))

d <- Diagonal(3)
(dC <- as(d, "CsparseMatrix")) # "dtCMatrix" (unitriangular)
(dgC <- as(dC, "generalMatrix"))
stopifnot(exprs = {
    is(dgC, "dgCMatrix") # was wrong in Matrix 1.3.2
    ## identical(dgC, as(dC, "dgCMatrix")) # deprecated
    identical(dC , new("dtCMatrix", p = rep(0L, 4), Dim = c(3L, 3L), diag = "U"))
    identical(dC , diagN2U(as(dgC, "triangularMatrix")))
})

setClass("posDef", contains = "dspMatrix")
N <- as(crossprod(M) + Diagonal(4), "packedMatrix")
(N <- new("posDef", N))
stopifnot(is(N[1:2, 1:2], "symmetricMatrix"))

#### Automatically display the class inheritance structure
#### possibly augmented with methods

allCl <- getClasses("package:Matrix")
cat("actual and virtual classes:\n")
tt <- table( isVirt <- sapply(allCl, isVirtualClass) )
names(tt) <- c('"actual"', "virtual")
tt
## The "actual" Matrix classes:
aCl <- allCl[!isVirt]
(aMcl <- aCl[grep("Matrix$", aCl)]) # length 48
aMc2 <-  aCl[sapply(aCl, extends, class2 = "Matrix")]
stopifnot(all( aMcl %in% aMc2 ))
aMc2[!(aMc2 %in% aMcl)] ## only 4 : p?Cholesky & p?BunchKaufman

## Really nice would be to construct an inheritance graph and display
## it.  Following things are computational variations on the theme..

## We use a version of  canCoerce()  that works with two *classes* instead of
## canCoerce <- function (object, Class)
classCanCoerce <- function (class1, class2)
{
    extends(class1, class2) ||
    !is.null(selectMethod("coerce", optional = TRUE,
			  signature    = c(from = class1, to = class2),
			  useInherited = c(from = TRUE,	  to = FALSE)))
}
.dq <- function(ch) paste0('"', ch, '"')
.subclasses <- function(cnam) {
    cd <- getClass(cnam)
    unique(c(cd@className, unlist(lapply(names(cd@subclasses), .subclasses))))
}
for(n in allCl) {
    if(isVirtualClass(n))
        cat("Virtual class", .dq(n),"\n")
    else {
        cat("\"Actual\" class", .dq(n),":\n")
        x <- new(n)
        if(doExtras) for(m in allCl)
            if(classCanCoerce(n,m)) {
                ext <- extends(n, m)
                if(ext) {
                    cat(sprintf("   extends  %20s %20s \n", "", .dq(m)))
                } else {
                    cat(sprintf("   can coerce: %20s -> %20s: ", .dq(n), .dq(m)))
                    tt <- try(as(x, m), silent = TRUE)
                    if(inherits(tt, "try-error")) {
                        cat("\t *ERROR* !!\n")
                    } else {
                        cat("as() ok; validObject: ")
                        vo <- validObject(tt, test = TRUE)
                        cat(if(isTRUE(vo)) "ok" else paste("OOOOOOPS:", vo), "\n")
                    }
                }
            }
        cat("---\n")
    }
}

cat('Time elapsed: ', proc.time(),'\n') # for the above "part I"


if(doExtras && !interactive()) { # don't want to see on source()

cat("All classes in the 'Matrix' package:\n")
for(cln in allCl) {
    cat("\n-----\n\nClass", dQuote(cln),":\n      ",
        paste(rep("~",nchar(cln)),collapse=''),"\n")
    ## A smarter version would use  getClass() instead of showClass(),
    ## build the "graph" and only then display.
    ##
    showClass(cln)
}

cat("\n\n")

## One could extend the `display' by using (something smarter than)
## are the "coerce" methods showing more than the 'Extends' output above?
cat("All (S4) methods in the 'Matrix' package:\n")
showMethods(where="package:Matrix")

} # end{non-interactive}

## 1-indexing instead of 0-indexing for direct "dgT" should give error:
ii <- as.integer(c(1,2,2))
jj <- as.integer(c(1,1,3))
assertError(new("dgTMatrix",  i=ii, j=jj,        x= 10*(1:3), Dim=2:3))
assertError(new("dgTMatrix",  i=ii, j=jj - 1:1,  x= 10*(1:3), Dim=2:3))
assertError(new("dgTMatrix",  i=ii - 1:1, j=jj,  x= 10*(1:3), Dim=2:3))
(mm <- new("dgTMatrix",  i=ii - 1:1, j=jj - 1:1, x= 10*(1:3), Dim=2:3))
validObject(mm)

### Sparse Logical:
m <- Matrix(c(0,0,2:0), 3,5)
mT <- as(mC <- as(m, "CsparseMatrix"), "TsparseMatrix")
stopifnot(identical(as(mT,"CsparseMatrix"), mC))
(mC. <- as(mT[1:2, 2:3], "CsparseMatrix"))
(mlC <- as(mC. , "lMatrix"))
as(mlC, "triangularMatrix")

chk0 <- function(M, kind, shape2 = "ge") {
## try({
    stopifnot(identical(dim(M), c(0L,0L)))
    kM <- Matrix:::.M.kind(M)
    if(!isValid(M, (cl <- paste0(kind, shape2, "Matrix"))))
        stop(sprintf('*not* "%s", but "%s"', cl, class(M)))
## })
}

for(knd in c("d","l","n","i","z")) {
    show(m0 <- new(paste0(knd, "geMatrix")))
    dkind <- if(knd == "z") "z" else "d"
    nlknd <- if(knd == "n") "n" else "l"
    ## for _all_ `knd` :
    for(dm in list(m0^2, 0 + m0, m0 / 2, m0 %% 2)) # Arith
        chk0(dm, dkind)
    ## Compare 1 :  (FIXME ?)
    for(dm in list(m0 != m0, m0 == m0, m0 >= m0)) # Compare
        chk0(dm, nlknd)
    ## Compare 2 :  (FIXME ?)
    for(dm in list(m0 <= 0, 1 >= m0)) # Compare
        chk0(dm, "l", shape2 = "sy")
    for(dm in list(!!m0, m0 | m0, m0 & m0)) # Logic
        chk0(dm, nlknd)
    stopifnot(identical(!!m0, .M2kind(m0, nlknd)))
    switch(knd,# integer arithmetic __ mostly not yet __
           "l" = ,
           "n" = for(dm in list(-m0 #_not yet_, 0L + m0, m0 %/% 2L, m0 %% 2L)
                                )) chk0(dm, kind = "i"),
           "i" = for(dm in list(-m0, 0L + m0, m0 %/% 2L, m0 %% 2L)
                                 ) chk0(dm, kind = "i"),
           "d" =,
           "z" = invisible(), # nothing
           stop("invalid 'knd':", knd))
}


if(!doExtras && !interactive()) q("no") ## (saving testing time)

### Test all classes:  validObject(new( * )) should be fulfilled -----------

## need stoplist for now:
Rcl.struc <- c("gR", "sR", "tR")
(dR.classes <- paste0(paste0("d", Rcl.struc[Rcl.struc != "gR"]),   "Matrix"))
(.R.classes <- paste0(sort(outer(c("l", "n"), Rcl.struc, paste0)), "Matrix"))
                                        # have only stub implementation

Mat.MatFact <- c("Cholesky", "pCholesky",
                 "BunchKaufman", "pBunchKaufman")##, "LDL"
##FIXME maybe move to ../../MatrixModels/tests/ :
## (modmat.classes <- .subclasses("modelMatrix"))
no.t.etc <- c(.R.classes, dR.classes, Mat.MatFact)#, modmat.classes)
no.t.classes <- c(no.t.etc)     # no t() available
no.norm.classes <- no.t.classes
not.Ops      <- NULL            # "Ops", e.g. "+" fails
not.coerce1  <- no.t.etc        # not coercable from "dgeMatrix"
not.coerce2  <- no.t.etc        # not coercable from "matrix"

tstMatrixClass <-
    function(cl, mM = Matrix(c(2,1,1,2) + 0, 2,2,
                 dimnames=rep( list(c("A","B")), 2)), # dimnames: *symmetric*
             mm = as(mM, "matrix"), recursive = TRUE, offset = 0)
{
    ## Purpose: Test 'Matrix' class {and do this for all of them}
    ## ----------------------------------------------------------------------
    ## Arguments: cl: class object of a class that extends "Matrix"
    ##            mM: a "Matrix"-matrix which will be coerced to class 'cl'
    ##            mm: a S3-matrix       which will be coerced to class 'cl'
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler

    ## from pkg sfsmisc :
    bl.string <- function(no) sprintf("%*s", no, "")

    ## Compute a few things only once :
    mM <- as(as(as(mM, "unpackedMatrix"), "generalMatrix"), "dMatrix") # dge
    trm <- mm; trm[lower.tri(mm)] <- 0
    ## not yet used:
    ## summList <- lapply(getGroupMembers("Summary"), get,
    ##                    envir = asNamespace("Matrix"))
    if(recursive)
        cList <- character(0)

    extraValid <- function(m, cl = class(m)) {
        sN <- slotNames(cl)
        sN <- sN[sN != "factors"]
	for(nm in sN)
	    if(!is.null(a <- attributes(slot(m, nm))))
		stop(sprintf("slot '%s' with %d attributes, named: ",
			     nm, length(a)), paste(names(a), collapse=", "))
        invisible(TRUE)
    }


    ## This is the recursive function
    dotestMat <- function(cl, offset)
    {
        cat. <- function(...) cat(bl.string(offset), ...)

        clNam <- cl@subClass
        deprecated <- grepl("^[nlidz](ge|sy|sp|po|pp|tr|tp|[gspt][CRT])Matrix$", clNam)
	cat("\n==>")
        cat.(clNam)
        ##---------
        clD <- getClassDef(clNam)
        if(isVirtualClass(clD)) {
            cat(" - is virtual\n")
            if(recursive) {
                cat.("----- begin{class :", clNam, "}----new subclasses----\n")
                for(ccl in clD@subclasses) {
                    cclN <- ccl@subClass
                    if(cclN %in% cList)
                        cat.(cclN,": see above\n")
                    else {
                        cList <<- c(cList, cclN)
                        dotestMat(ccl, offset = offset + 3)
                    }
                }
                cat.("----- end{class :", clNam, "}---------------------\n")
            }
        } else { ## --- actual class ---
            genC <- extends(clD, "generalMatrix")
            symC <- extends(clD, "symmetricMatrix")
            triC <- extends(clD, "triangularMatrix")
            diaC <- extends(clD, "diagonalMatrix")
            indC <- extends(clD, "indMatrix")
            if(!(genC || symC || triC || diaC || indC))
                stop("does not extend one of 'general', 'symmetric', 'triangular', 'diagonal', 'ind'")
            sparseC <- extends(clD, "sparseMatrix")
            denseC  <- extends(clD, "denseMatrix")
            if(!(sparseC || denseC))
                stop("does not extend either 'sparse' or 'dense'")
	    cat("; new(*): ")
            m <- new(clNam) ; cat("ok; ")
	    m0 <- matrix(,0,0)
	    if(canCoerce(m0, clNam)) {
		cat("; canCoerce(matrix(,0,0), *) => as(m0, <.>): ")
                m0. <-
                    if(deprecated)
                        eval(Matrix:::.as.via.virtual(
                                          "matrix", clNam, quote(m0)))
                    else as(m0, clNam)
                if(.hasSlot(m, "diag") && .hasSlot(m0., "diag") &&
                   identical(m@diag, "N") && identical(m0.@diag, "U"))
                    ## tolerate as(0-by-0, .) formally having unit diagonal
                    m0.@diag <- "N"
		stopifnot(Qidentical(m, m0.)); cat("ok; ")
	    }
            is_p <- extends(clD, "indMatrix")
            is_cor <- extends(clD, "corMatrix") || extends(clD, "copMatrix")
            ## ^^^ has diagonal divided out
	    if(canCoerce(mm, clNam)) { ## replace 'm' by `non-empty' version
		cat("canCoerce(mm, *) ")
		m0 <- {
		    if(triC) trm
		    else if(is_p)
			mm == 1     # logical *and* "true" permutation
		    else mm
		}
		if(extends(clD, "lMatrix") ||
		   extends(clD, "nMatrix"))
		    storage.mode(m0) <- "logical"
		else if(extends(clD, "zMatrix"))
		    storage.mode(m0) <- "complex"
		validObject(m) ## validity of trivial 'm' before replacing
		m <-
                    if(deprecated)
                        eval(Matrix:::.as.via.virtual(
                                          "matrix", clNam, quote(m0)))
                    else as(m0, clNam)
		if(is_cor)
                    m0 <- cov2cor(m0)
	    } else {
                m0 <- vector(Matrix:::.type.kind[Matrix:::.M.kind(m)])
                dim(m0) <- c(0L,0L)
            }
	    ## m0 is the 'matrix' version of our 'Matrix' m
	    m. <- m0 ##m. <- if(is_p) as.integer(m0) else m0
            EQ <- if(is_cor) all.equal else identical
	    stopifnot(EQ(m0[FALSE], m[FALSE])
		      , EQ(m.[TRUE],  m[TRUE])
		      , if(length(m) >= 2) EQ(m.[2:1], m[2:1]) else TRUE)

	    if(all(dim(m) > 0)) { ## matrix(0,0,0)[FALSE,]  is invalid too
		m00 <- m[FALSE,FALSE]
		m.. <- m[TRUE , TRUE]
		stopifnot(dim(m00) == c(0L,0L),
			  dim(m..) == dim(m))
		## not yet  , class(m00) == clNam , identical(m.. , m)
	    }

            cat("valid: ", validObject(m), extraValid(m, clNam),"\n")

            ## This can only work as long as 'm' has no NAs :
            ## not yet -- have version in not.Ops below
            ## once we have is.na():
            ## 		stopifnot(all(m == m | is.na(m))) ## check all() and "==" [Compare]
            ## 		if(any(m != m && !is.na(m)))

	    show(m)
	    ## coerce to 'matrix'
	    m.m <- as(m, "matrix")
	    ##=========##
	    checkMatrix(m, m.m,
	    ##=========##
			do.t= !(clNam %in% no.t.classes),
			doNorm= !(clNam %in% no.norm.classes),
			doOps = all(clNam != not.Ops),
			doCoerce = all(clNam != not.coerce1),
			catFUN = cat.)

### FIXME: organize differently :
### 1) produce 'mM'  and 'mm' for the other cases,
### 2) use identical code for all cases

            if(is(m, "dMatrix") &&
               (is(m, "generalMatrix") || is(m, "symmetricMatrix"))) {
                if(any(clNam == not.coerce1))
                    cat.("not coercable_1\n")
                else if(canCoerce(mM, clNam)) {
                    m2 <-
                        if(deprecated)
                            eval(Matrix:::.as.via.virtual(
                                              class(mM), clNam, quote(mM)))
                        else as(mM, clNam)
                    cat("valid:", validObject(m2), "\n")
                    if(!is_cor) ## as.vector()
                        stopifnot(as.vector(m2) == as.vector(mM))
                    cat.("[cr]bind2():"); mm2 <- cbind2(m2,m2)
                    stopifnot(dim(rbind2(m2,m2)) == 2:1 * dim(mM)); cat(" ok")
                    if(genC && class(mm2) == clNam) ## non-square matrix when "allowed"
                        m2 <- mm2
                    dd <- diag(m2)
                    cat("; `diag<-` ")
                    diag(m2) <- 10*dd
                    stopifnot(is_cor || identical(dd, diag(mM)),
                              identical(10*dd, diag(m2))); cat("ok ")
                }
##                 if(all(clNam != not.coerce2)) {
                    if(canCoerce("matrix", clNam)) {
                        cat.("as(matrix, <class>): ")
                        m3 <-
                            if(deprecated)
                                eval(Matrix:::.as.via.virtual(
                                                  "matrix", clNam, quote(mm)))
                            else as(mm, clNam)
                        cat("valid:", validObject(m3), "\n")
                    } else cat.(" not coerceable from \"matrix\"\n")
##                 }
            }

            ## else { ... no happens in tstMatrix() above .. }

            ## if(is(m, "denseMatrix")) {
            ##     ## .........
            ##     cat.("as dsparse* ")
            ##     msp <- as(m, "dsparseMatrix")
            ##     cat.("; valid coercion: ", validObject(msp), "\n")
            ## } else if(is(m, "sparseMatrix")) {

            ## } else cat.("-- not dense nor sparse -- should not happen(!?)\n")

            if(is(m, "dsparseMatrix")) {
                if(any(clNam == not.coerce1))
                    cat.("not coercable_1\n")
                else {
                    ## make sure we can coerce to dgT* -- needed, e.g. for "image"
                    cat.("as dgT* ")
                    mgT <- eval(Matrix:::.as.via.virtual(
                                             class(m), "dgTMatrix", quote(m)))
                    cat(sprintf("; valid dgT* coercion: %s\n",
                                validObject(mgT)))
                }
            }
        }
    } # end{dotestMat}

    for(scl in getClass(cl)@subclasses)
        dotestMat(scl, offset + 1)
}
## in case we want to make progress:
## codetools::checkUsage(tstMatrixClass, all=TRUE)

tstMatrixClass("Matrix")
if(FALSE)## or just a sub class
tstMatrixClass("triangularMatrix")


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

if(!interactive()) warnings()
