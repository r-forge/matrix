## ~~~~ BACKWARDS COMPATIBILITY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.Rv <- getRversion() # remove after using!
if(.Rv < "4.0.0") {
    deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...)
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    ...length <- function()
        ## not equivalent to base::...length, as list(...) evaluates
        eval(quote(length(list(...))), sys.frame(-1L))
    tryInvokeRestart <- function(r, ...) {
        if(!isRestart(r))
            r <- findRestart(r)
        if(is.null(r))
            invisible(NULL)
        else .Internal(.invokeRestart(r, list(...)))
    }
}
rm(.Rv)


## ~~~~ PACKAGE ENVIRONMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Supporting once-per-session warnings, etc.
.MatrixEnv <- new.env(parent = emptyenv(), hash = FALSE)

## Storing settings from 'cholmod_common'
.chm_common <- new.env(parent = emptyenv())


## ~~~~ NAMESPACE HOOKS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.onLoad <- function(libname, pkgname) {
    verbose <-
        (!is.null(v <- getOption("Matrix.verbose")) && v >= 1L) ||
        (interactive() && identical(Sys.info()[["user"]], "maechler"))

    ## ambiguityNotes:
    ## show S4 method dispatch ambiguity notes if TRUE
    aN <-
        !is.null(getOption("ambiguousMethodSelection")) ||
        verbose || isTRUE(getOption("Matrix.ambiguityNotes"))
    if(!aN)
        options(ambiguousMethodSelection =
                    `environment<-`(function(cond) NULL, emptyenv()))
    assign("ambiguityNotes", aN, envir = .MatrixEnv)

    ## warnDeprecatedCoerce:
    ## <=0 ... no conditions signaled
    ##   1 ... persistent warning
    ## >=2 ... persistent error
    ##  NA ... one-time message { d(g.|.C)Matrix } or warning { others }
    wDC <-
        if(is.atomic(w <- getOption("Matrix.warnDeprecatedCoerce")) &&
           length(w) == 1L)
            as.integer(w)
        else if(verbose)
            1L
        else NA_integer_
    if(!is.null(w))
        options(Matrix.warnDeprecatedCoerce = wDC)
    assign("warnDeprecatedCoerce", wDC, envir = .MatrixEnv)

    .Call(CHM_set_common_env, .chm_common)
    NULL
}

.onUnload <- function(libpath) {
    library.dynam.unload("Matrix", libpath)
    if(!isTRUE(.MatrixEnv[["ambiguityNotes"]]))
	options(ambiguousMethodSelection = NULL)
    NULL
}


## ~~~~ DEPRECATED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Utility for Matrix.DeprecatedCoerce(); see below
.as.via.virtual <- function(Class1, Class2, from = quote(from)) {
    if(!isClassDef(Class1))
        Class1 <- getClassDef(Class1)
    if(!isClassDef(Class2))
        Class2 <- getClassDef(Class2)
    if(!grepl("^[dln](di|ge|tr|sy|tp|sp|[gts][CRT])Matrix$", Class2@className))
        stop("invalid 'Class2'")
    contains1 <- names(Class1@contains)
    contains2 <- names(Class2@contains)
    virtual <- list(c("dMatrix", "lMatrix", "nMatrix"),
                    c("generalMatrix", "triangularMatrix", "symmetricMatrix"),
                    c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix",
                      "diagonalMatrix", "unpackedMatrix", "packedMatrix"))
    to <- from
    for(v in virtual) {
        if(any(m <- match(v, contains2, 0L) > 0L)) {
            v1 <- v[m][1L]
            if(match(v1, contains1, 0L) == 0L)
                to <- call("as", to, v1)
        }
    }
    to
}

Matrix.DeprecatedCoerce <- function(Class1, Class2) {
    if(!isClassDef(Class1))
        Class1 <- getClassDef(Class1)
    if(!isClassDef(Class2))
        Class2 <- getClassDef(Class2)
    if(is.null(w <- getOption("Matrix.warnDeprecatedCoerce")))
        w <- .MatrixEnv[["warnDeprecatedCoerce"]]
    if(is.atomic(w) && length(w) == 1L &&
       ((w.na <- is.na(w <- as.integer(w))) || w > 0L)) {
        cln1 <- Class1@className
        cln2 <- Class2@className
        warning. <-
            if(w.na && grepl("d(g.|.C)Matrix", cln2))
                function(..., call., domain) message(..., domain = domain)
            else if(w.na || w == 1L)
                warning
            else stop
        warning.(gettextf("as(<%s>, \"%s\") is deprecated since Matrix 1.4-2; do %s instead",
                          cln1, cln2,
                          deparse1(.as.via.virtual(Class1, Class2, quote(.)))),
                 call. = FALSE, domain = NA)
        if(w.na)
            options(Matrix.warnDeprecatedCoerce = 0L)
    }
    invisible(NULL)
}

## "Granular" coercions available in Matrix 1.4-1,
## all candidates for deprecation in Matrix 1.4-2:
if(FALSE) {
stopifnot(packageVersion("Matrix") == "1.4.1")
dput(lapply(grep("to=\"[dln](di|ge|tr|sy|tp|sp|[gts][CRT])Matrix\"",
                 capture.output(showMethods("coerce")),
                 value = TRUE),
            function(s) unname(eval(str2lang(paste0("c(", s, ")"))))))
}
.from.to <- list(c("ddenseMatrix", "dgeMatrix"), c("ddiMatrix", "dgCMatrix"),
                 c("ddiMatrix", "dgeMatrix"), c("ddiMatrix", "dtCMatrix"),
                 c("dgCMatrix", "dgeMatrix"), c("dgCMatrix", "dgTMatrix"),
                 c("dgCMatrix", "dsCMatrix"), c("dgCMatrix", "dtCMatrix"),
                 c("dgCMatrix", "lgCMatrix"), c("dgCMatrix", "ngCMatrix"),
                 c("dgeMatrix", "dgCMatrix"), c("dgeMatrix", "dgTMatrix"),
                 c("dgeMatrix", "dspMatrix"), c("dgeMatrix", "dsTMatrix"),
                 c("dgeMatrix", "dsyMatrix"), c("dgeMatrix", "dtrMatrix"),
                 c("dgeMatrix", "lgeMatrix"), c("dgTMatrix", "dgCMatrix"),
                 c("dgTMatrix", "dgeMatrix"), c("dgTMatrix", "dsTMatrix"),
                 c("dgTMatrix", "dtCMatrix"), c("dgTMatrix", "dtTMatrix"),
                 c("dsCMatrix", "dgCMatrix"), c("dsCMatrix", "dgeMatrix"),
                 c("dsCMatrix", "dgTMatrix"), c("dsCMatrix", "dsRMatrix"),
                 c("dsCMatrix", "dsTMatrix"), c("dsCMatrix", "dsyMatrix"),
                 c("dsCMatrix", "lsCMatrix"), c("dsCMatrix", "nsCMatrix"),
                 c("dspMatrix", "dsyMatrix"), c("dspMatrix", "lspMatrix"),
                 c("dsTMatrix", "dgeMatrix"), c("dsTMatrix", "dgTMatrix"),
                 c("dsTMatrix", "dsCMatrix"), c("dsTMatrix", "dsyMatrix"),
                 c("dsTMatrix", "lsTMatrix"), c("dsyMatrix", "dsCMatrix"),
                 c("dsyMatrix", "dspMatrix"), c("dsyMatrix", "dsTMatrix"),
                 c("dsyMatrix", "lsyMatrix"), c("dtCMatrix", "dgCMatrix"),
                 c("dtCMatrix", "dgeMatrix"), c("dtCMatrix", "dgTMatrix"),
                 c("dtCMatrix", "dsCMatrix"), c("dtCMatrix", "dtrMatrix"),
                 c("dtCMatrix", "dtTMatrix"), c("dtCMatrix", "ltCMatrix"),
                 c("dtCMatrix", "ntCMatrix"), c("dtpMatrix", "dtrMatrix"),
                 c("dtpMatrix", "dtTMatrix"), c("dtpMatrix", "ltpMatrix"),
                 c("dtrMatrix", "dtpMatrix"), c("dtrMatrix", "ltrMatrix"),
                 c("dtTMatrix", "dgeMatrix"), c("dtTMatrix", "dgTMatrix"),
                 c("dtTMatrix", "dtCMatrix"), c("dtTMatrix", "dtrMatrix"),
                 c("indMatrix", "ngeMatrix"), c("indMatrix", "ngTMatrix"),
                 c("lgCMatrix", "dgCMatrix"), c("lgCMatrix", "lgeMatrix"),
                 c("lgCMatrix", "lgTMatrix"), c("lgCMatrix", "ltCMatrix"),
                 c("lgeMatrix", "dgeMatrix"), c("lgeMatrix", "lgCMatrix"),
                 c("lgeMatrix", "lgTMatrix"), c("lgeMatrix", "lspMatrix"),
                 c("lgeMatrix", "lsyMatrix"), c("lgeMatrix", "ltpMatrix"),
                 c("lgeMatrix", "ltrMatrix"), c("lgTMatrix", "dgTMatrix"),
                 c("lgTMatrix", "lgCMatrix"), c("lgTMatrix", "lgeMatrix"),
                 c("lgTMatrix", "lsCMatrix"), c("lgTMatrix", "ltTMatrix"),
                 c("lMatrix", "dgCMatrix"), c("lsCMatrix", "dsCMatrix"),
                 c("lsCMatrix", "lgCMatrix"), c("lsCMatrix", "lgTMatrix"),
                 c("lsCMatrix", "lsTMatrix"), c("lspMatrix", "dspMatrix"),
                 c("lspMatrix", "lgeMatrix"), c("lspMatrix", "lsyMatrix"),
                 c("lsTMatrix", "lgCMatrix"), c("lsTMatrix", "lgTMatrix"),
                 c("lsTMatrix", "lsCMatrix"), c("lsTMatrix", "lsyMatrix"),
                 c("lsyMatrix", "dsyMatrix"), c("lsyMatrix", "lgeMatrix"),
                 c("lsyMatrix", "lspMatrix"), c("ltCMatrix", "lgCMatrix"),
                 c("ltCMatrix", "ltTMatrix"), c("ltpMatrix", "dtpMatrix"),
                 c("ltpMatrix", "lgeMatrix"), c("ltpMatrix", "ltrMatrix"),
                 c("ltrMatrix", "dtrMatrix"), c("ltrMatrix", "lgeMatrix"),
                 c("ltrMatrix", "ltpMatrix"), c("ltTMatrix", "dtTMatrix"),
                 c("ltTMatrix", "lgCMatrix"), c("ltTMatrix", "lgTMatrix"),
                 c("ltTMatrix", "ltCMatrix"), c("ltTMatrix", "ltrMatrix"),
                 ## c("matrix.coo", "dgCMatrix"), c("matrix.coo", "dgTMatrix"),
                 ## c("matrix.csc", "dgCMatrix"), c("matrix.csr", "dgCMatrix"),
                 ## c("matrix.csr", "dgRMatrix"),
                 c("matrix", "dgCMatrix"),
                 c("matrix", "dgeMatrix"), c("matrix", "dgRMatrix"),
                 c("matrix", "dgTMatrix"), c("matrix", "dsCMatrix"),
                 c("matrix", "dspMatrix"), c("matrix", "dsTMatrix"),
                 c("matrix", "dsyMatrix"), c("matrix", "dtCMatrix"),
                 c("matrix", "dtpMatrix"), c("matrix", "dtrMatrix"),
                 c("matrix", "dtTMatrix"), c("matrix", "lgCMatrix"),
                 c("matrix", "lgeMatrix"), c("matrix", "lgTMatrix"),
                 c("matrix", "lsCMatrix"), c("matrix", "lspMatrix"),
                 c("matrix", "lsyMatrix"), c("matrix", "ltCMatrix"),
                 c("matrix", "ltpMatrix"), c("matrix", "ltrMatrix"),
                 c("matrix", "ltTMatrix"), c("matrix", "ngCMatrix"),
                 c("matrix", "ngeMatrix"), c("matrix", "ngTMatrix"),
                 c("matrix", "nspMatrix"), c("matrix", "nsyMatrix"),
                 c("matrix", "ntCMatrix"), c("matrix", "ntpMatrix"),
                 c("matrix", "ntrMatrix"), c("matrix", "ntTMatrix"),
                 c("ngCMatrix", "dgCMatrix"), c("ngCMatrix", "lgCMatrix"),
                 c("ngCMatrix", "ntCMatrix"), c("ngeMatrix", "dgeMatrix"),
                 c("ngeMatrix", "lgeMatrix"), c("ngeMatrix", "ngCMatrix"),
                 c("ngeMatrix", "ngTMatrix"), c("ngeMatrix", "nspMatrix"),
                 c("ngeMatrix", "nsyMatrix"), c("ngeMatrix", "ntpMatrix"),
                 c("ngeMatrix", "ntrMatrix"), c("ngTMatrix", "dgTMatrix"),
                 c("ngTMatrix", "lgeMatrix"), c("ngTMatrix", "lgTMatrix"),
                 c("ngTMatrix", "ngCMatrix"), c("ngTMatrix", "ngeMatrix"),
                 c("ngTMatrix", "ntTMatrix"), c("nsCMatrix", "dsCMatrix"),
                 c("nsCMatrix", "lsCMatrix"), c("nsCMatrix", "ngCMatrix"),
                 c("nsCMatrix", "nsTMatrix"), c("nspMatrix", "dspMatrix"),
                 c("nspMatrix", "lspMatrix"), c("nspMatrix", "ngeMatrix"),
                 c("nspMatrix", "nsyMatrix"), c("nsTMatrix", "dsTMatrix"),
                 c("nsTMatrix", "ngCMatrix"), c("nsTMatrix", "ngTMatrix"),
                 c("nsTMatrix", "nsCMatrix"), c("nsTMatrix", "nsyMatrix"),
                 c("nsyMatrix", "dsyMatrix"), c("nsyMatrix", "lsyMatrix"),
                 c("nsyMatrix", "ngeMatrix"), c("nsyMatrix", "nspMatrix"),
                 c("ntCMatrix", "dtCMatrix"), c("ntCMatrix", "ltCMatrix"),
                 c("ntCMatrix", "ngCMatrix"), c("ntpMatrix", "dtpMatrix"),
                 c("ntpMatrix", "ltpMatrix"), c("ntpMatrix", "ngeMatrix"),
                 c("ntpMatrix", "ntrMatrix"), c("ntrMatrix", "dtrMatrix"),
                 c("ntrMatrix", "ltrMatrix"), c("ntrMatrix", "ngeMatrix"),
                 c("ntrMatrix", "ntpMatrix"), c("ntTMatrix", "dtTMatrix"),
                 c("ntTMatrix", "ngCMatrix"), c("ntTMatrix", "ngTMatrix"),
                 c("ntTMatrix", "ntCMatrix"), c("ntTMatrix", "ntrMatrix"),
                 c("numLike", "dgeMatrix"), c("RsparseMatrix", "dgeMatrix"))

.def.template <- function(from) {
    cd1 <- getClassDef(.FROM)
    cd2 <- getClassDef(.TO)
    Matrix.DeprecatedCoerce(cd1, cd2);
    to <- .CALL
    if(identical(as.character(class(to)), .TO))
        return(to)
    ## Coercion via virtual generated a _subclass_ of the target class
    to.strict <- new(.TO)
    for (nm in slotNames(cd2))
        slot(to.strict, nm) <- slot(to, nm)
    to.strict
}
for (.f.t in .from.to) {
    .f <- .f.t[1L]
    .t <- .f.t[2L]
    .def <- .def.template
    .env <- list(.FROM = .f, .TO = .t, .CALL = .as.via.virtual(.f, .t))
    body(.def) <- do.call(substitute, list(body(.def), .env))
    setAs(.f, .t, .def)
}
rm(.from.to, .f.t, .f, .t, .def.template, .def, .env)


## ~~~~ DEFUNCT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cBind <- function(..., deparse.level = 1) {
    .Defunct(msg = "'cBind' is defunct; 'base::cbind' handles S4 objects since R 3.2.0")
}
rBind <- function(..., deparse.level = 1) {
    .Defunct(msg = "'rBind' is defunct; 'base::rbind' handles S4 objects since R 3.2.0")
}


## ~~~~ "MISCELLANEOUS" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.SuiteSparse_version <- function() {
    v <- .Call(get_SuiteSparse_version)
    package_version(list(major = v[1L], minor = paste(v[2:3], collapse = ".")))
}
