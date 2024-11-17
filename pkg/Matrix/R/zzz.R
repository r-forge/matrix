## ~~~~ VERSION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Matrix.Version <- function() {
    n <- .Call(R_Matrix_version)
    v <- .mapply(function(n, p, b, class) {
                     r <- integer(p)
                     while (p > 0L) {
                         r[p] <- tmp <- n %% b
                         n <- (n - tmp) %/% b
                         p <- p - 1L
                     }
                     v <- list(r)
                     class(v) <- c(class, "numeric_version")
                     v
                 },
                 list(n = n, p = c(3L, 1L, 3L), b = c(256L, 10L, 256L),
                      class = list("package_version", NULL, NULL)),
                 NULL)
    names(v) <- names(n)
    v
}


## ~~~~ ENVIRONMENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Recording default values of Matrix.* options
.MatrixEnv <- new.env(parent = emptyenv(), hash = FALSE)


## ~~~~ NAMESPACE HOOKS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.onLoad <- function(libname, pkgname) {
    ## For backwards compatibility with earlier versions of R,
    ## at least until x.y.z if we have Depends: R (>= x.y.z)
    Mns <- parent.env(environment())
    if(!environmentIsLocked(Mns)) {
        ## Namespace not locked yet, but being defensive here
        Rv <- getRversion()
        if(Rv < "4.4.0") {
        assign("%||%", envir = Mns, inherits = FALSE,
               function(x, y) if(is.null(x)) y else x)
        if(Rv < "4.1.3") {
        assign("...names", envir = Mns, inherits = FALSE,
               function() eval(quote(names(list(...))), sys.frame(-1L)))
        if(Rv < "4.0.0") {
        assign("deparse1", envir = Mns, inherits = FALSE,
               function(expr, collapse = " ", width.cutoff = 500L, ...)
                   paste(deparse(expr, width.cutoff, ...),
                         collapse = collapse))
        assign("sequence.default", envir = Mns, inherits = FALSE,
               function(nvec, from = 1L, by = 1L, ...) {
                   if(length(nvec) == 0L)
                       return(integer(0L))
                   else if(length(from) == 0L || length(by) == 0L)
                       stop(gettextf("'%s' has length 0 but '%s' does not",
                                     if(length(from) == 0L) "from" else "by", "nvec"),
                            domain = NA)
                   unlist(.mapply(seq.int,
                                  list(from = as.integer(from),
                                       by = as.integer(by),
                                       length.out = as.integer(nvec)),
                                  NULL),
                          recursive = FALSE, use.names = FALSE)
               })
        assign("tryInvokeRestart", envir = Mns, inherits = FALSE,
               function(r, ...)
                   tryCatch(invokeRestart(r, ...),
                            error = function(e) invisible(NULL)))
        } # Rv < "4.0.0"
        } # Rv < "4.1.3"
        } # Rv < "4.4.0"
    }

    ## verbose:
    ## logical/integer (but often supplied as double),
    ## deciding _if_ conditions are signaled
    v <- as.integer(Sys.getenv("R_MATRIX_VERBOSE", "0"))
    assign("verbose", if(is.na(v)) 0L else v, envir = .MatrixEnv)

    ## warn:
    ## logical/integer (but often supplied as double),
    ## deciding _what_ conditions are signaled
    ## (0=message, 1=warning, 2=error)
    w <- as.integer(Sys.getenv("R_MATRIX_WARN", "0"))
    assign("warn", if(is.na(w)) 0L else w, envir = .MatrixEnv)

    ## ambiguityNotes:
    ## show S4 method dispatch ambiguity notes if TRUE
    aN <- as.logical(Sys.getenv("R_MATRIX_AMBIGUITY_NOTES", "false"))
    aN <- (!is.na(aN) && aN) || !is.null(getOption("ambiguousMethodSelection"))
    assign("ambiguityNotes", aN, envir = .MatrixEnv)
    if(!aN)
        options(ambiguousMethodSelection = # ?methods::testInheritedMethods
                    `environment<-`(function(cond) NULL, emptyenv()))

    ## warnDeprecatedCoerce:
    ## <=0 ... no conditions signaled
    ##   1 ... persistent warning
    ## >=2 ... persistent error
    ##  NA ... one-time message { d(g.|.C)Matrix } or warning { others }
    wDC <- as.integer(Sys.getenv("R_MATRIX_WARN_DEPRECATED_COERCE", NA))
    assign("warnDeprecatedCoerce", wDC, envir = .MatrixEnv)

    ## warnSqrtDefault:
    ## <=0 ... no conditions signaled
    ##   1 ... persistent warning
    ## >=2 ... persistent error
    ##  NA ... one-time warning
    wSD <- as.integer(Sys.getenv("R_MATRIX_WARN_SQRT_DEFAULT", NA))
    assign("warnSqrtDefault", wSD, envir = .MatrixEnv)

    NULL
}

.onUnload <- function(libpath) {
    library.dynam.unload("Matrix", libpath)
    if(!.MatrixEnv[["ambiguityNotes"]])
        options(ambiguousMethodSelection = NULL)
    NULL
}


## ~~~~ DEPRECATED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

..2dge <- function(from) {
    .Deprecated(new = ".M2gen(*, \"d\") or .m2dense(*, \"dge\")", package = "Matrix")
    if(isS4(from))
        .M2gen(from, "d")
    else .m2dense(from, "dge")
}
.C2nC <- function(from, isTri) {
    .Deprecated(new = ".M2kind", package = "Matrix")
    .M2kind(from, "n")
}
.T2Cmat <- function(from, isTri) {
    .Deprecated(new = ".M2C", package = "Matrix")
    .M2C(from)
}
.asmatrix <- function(x) {
    .Deprecated(new = "as(., \"matrix\")", package = "Matrix")
    as(x, "matrix")
}
.dense2sy <- function(from, ...) {
    .Deprecated(new = ".M2sym", package = "Matrix")
    .M2sym(from, ...)
}
.diag2mat <- function(from) {
    .Deprecated(new = ".M2m", package = "Matrix")
    .M2m(from)
}
.diag2sT <- function(from, uplo = "U", kind = ".", drop0 = TRUE) {
    .Deprecated(new = ".diag2sparse", package = "Matrix")
    r <- .diag2sparse(from, kind, "s", "T", uplo)
    if(drop0)
        r <- .drop0(r)
    r
}
.diag2tT <- function(from, uplo = "U", kind = ".", drop0 = TRUE) {
    .Deprecated(new = ".diag2sparse", package = "Matrix")
    to <- .diag2sparse(from, kind, "t", "T", uplo)
    if(drop0)
        to <- .drop0(to)
    to
}
.dsy2dsp <- function(from) {
    .Deprecated(new = ".M2packed", package = "Matrix")
    .M2packed(from)
}
.dsy2mat <- function(from, keep.dimnames = TRUE) {
    .Deprecated(new = ".M2m", package = "Matrix")
    to <- .M2m(from)
    if(!keep.dimnames)
        dimnames(to) <- NULL
    to
}
.dxC2mat <- function(from, chkUdiag) {
    .Deprecated(new = ".M2m", package = "Matrix")
    .M2m(from)
}
.m2dgC <- function(from) {
    .Deprecated(new = ".m2sparse", package = "Matrix")
    .m2sparse(from, "dgC")
}
.m2lgC <- function(from) {
    .Deprecated(new = ".m2sparse", package = "Matrix")
    .m2sparse(from, "lgC")
}
.m2ngC <- function(from) {
    .Deprecated(new = ".m2sparse", package = "Matrix")
    if(anyNA(from))
        stop(gettextf("attempt to coerce matrix with NA to %s", "ngCMatrix"),
             domain = NA)
    .m2sparse(from, "ngC")
}
.m2ngCn <- function(from, na.is.not.0 = FALSE) {
    .Deprecated(new = ".m2sparse", package = "Matrix")
    if(!na.is.not.0 && anyNA(from))
        stop(gettextf("attempt to coerce matrix with NA to %s", "ngCMatrix"),
             domain = NA)
    .m2sparse(from, "ngC")
}
.m2ngTn <- function(from, na.is.not.0 = FALSE) {
    .Deprecated(new = ".m2sparse", package = "Matrix")
    if(!na.is.not.0 && anyNA(from))
        stop(gettextf("attempt to coerce matrix with NA to %s", "ngTMatrix"),
             domain = NA)
    .m2sparse(from, "ngT")
}
.n2dgT <- function(from) {
    .Deprecated(new = ".M2kind", package = "Matrix")
    .M2kind(from, "d")
}
.nC2d <- function(from) {
    .Deprecated(new = ".M2kind", package = "Matrix")
    .M2kind(from, "d")
}
.nC2l <- function(from) {
    .Deprecated(new = ".M2kind", package = "Matrix")
    .M2kind(from, "l")
}

.dense2m <- .sparse2m <- function(from) {
    if(FALSE) {
    .Deprecated(new = ".M2m", package = "Matrix")
    }
    .M2m(from)
}

.dense2v <- .sparse2v <- function(from) {
    if(FALSE) {
    .Deprecated(new = ".M2v", package = "Matrix")
    }
    .M2v(from)
}

.dense2kind <- function(from, kind) {
    if(FALSE) {
    .Deprecated(new = ".M2kind", package = "Matrix")
    }
    .M2kind(from, kind)
}

.sparse2kind <- function(from, kind, drop0 = FALSE) {
    if(FALSE) {
    .Deprecated(new = ".M2kind", package = "Matrix")
    }
    .M2kind(if(drop0) .drop0(from) else from, kind)
}

.dense2g <- .sparse2g <- function(from, kind = ".") {
    if(FALSE) {
    .Deprecated(new = ".M2gen", package = "Matrix")
    }
    .M2gen(from, kind)
}

.CR2RC <- function(from) {
    if(.M.repr(from) != "C") {
        if(FALSE) {
        .Deprecated(new = ".M2C", package = "Matrix")
        }
        .M2C(from)
    } else {
        if(FALSE) {
        .Deprecated(new = ".M2R", package = "Matrix")
        }
        .M2R(from)
    }
}

.CR2T <- function(from) {
    if(FALSE) {
    .Deprecated(new = ".M2T", package = "Matrix")
    }
    .M2T(from)
}

.T2CR <- function(from, Csparse = TRUE) {
    if(Csparse) {
        if(FALSE) {
        .Deprecated(new = ".M2C", package = "Matrix")
        }
        .M2C(from)
    } else {
        if(FALSE) {
        .Deprecated(new = ".M2R", package = "Matrix")
        }
        .M2R(from)
    }
}

.tCR2RC <- function(from) {
    if(FALSE) {
    .Deprecated(new = ".tCRT", package = "Matrix")
    }
    .tCRT(from)
}

uniqTsparse <- function(x, class.x = class(x)) {
    if(FALSE) {
    .Deprecated(new = "asUniqueT", package = "Matrix")
    }
    asUniqueT(x, isT = extends(class.x, "TsparseMatrix"))
}

.SuiteSparse_version <- function() {
    if(FALSE) {
    .Deprecated(new = "Matrix.Version", package = "Matrix")
    }
    Matrix.Version()[["suitesparse"]]
}

## Utility for Matrix.DeprecatedCoerce(); see below
.as.via.virtual <- function(Class1, Class2, from = quote(from)) {
    if(!isClassDef(Class1))
        Class1 <- getClassDef(Class1)
    if(!isClassDef(Class2))
        Class2 <- getClassDef(Class2)
    if(!grepl("^[nlidz](ge|sy|sp|po|pp|tr|tp|[gspt][CRT]|di)Matrix$", Class2@className))
        stop("invalid 'Class2'")
    contains1 <- names(Class1@contains)
    contains2 <- names(Class2@contains)
    virtual <- list(c("nMatrix", "lMatrix", "iMatrix",
                      "dMatrix", "zMatrix"),
                    c("triangularMatrix", "posdefMatrix",
                      "symmetricMatrix", "generalMatrix"),
                    c("diagonalMatrix", "CsparseMatrix",
                      "RsparseMatrix", "TsparseMatrix",
                      "unpackedMatrix", "packedMatrix"))
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
    w <- getOption("Matrix.warnDeprecatedCoerce",
                   .MatrixEnv[["warnDeprecatedCoerce"]])
    if(is.atomic(w) && length(w) == 1L &&
       ((w.na <- is.na(w <- as.integer(w))) || w > 0L)) {
        cln1 <- Class1@className
        cln2 <- Class2@className

        old <- sprintf("as(<%s>, \"%s\")", cln1, cln2)
        new <- deparse1(.as.via.virtual(Class1, Class2, quote(.)))

        if(w.na)
            on.exit(options(Matrix.warnDeprecatedCoerce = 0L))
        if(w.na && grepl("d(g.|.C)Matrix", cln2)) {
            cond <-
                tryCatch(.Deprecated(old = old, new = new, package = "Matrix"),
                         deprecatedWarning = identity)
            message(conditionMessage(cond), domain = NA)
        } else {
            if(!w.na && w > 1L) {
                oop <- options(warn = 2L)
                on.exit(options(oop))
            }
            .Deprecated(old = old, new = new, package = "Matrix")
        }
    }
    invisible(NULL)
}

if (FALSE)
dput(lapply(grep("to=\"[nlidz](ge|sy|sp|po|pp|tr|tp|[gspt][CRT]|di)Matrix\"",
                 capture.output(showMethods("coerce")),
                 value = TRUE),
            function(s) unname(eval(str2lang(paste0("c(", s, ")"))))))
.from.to <- list(c("Matrix", "dpoMatrix"),
                 c("Matrix", "dppMatrix"),
                 c("RsparseMatrix", "dgeMatrix"),
                 c("ddenseMatrix", "dgeMatrix"),
                 c("ddiMatrix", "dgCMatrix"),
                 c("ddiMatrix", "dgeMatrix"),
                 c("ddiMatrix", "dtCMatrix"),
                 c("dgCMatrix", "dgTMatrix"),
                 c("dgCMatrix", "dgeMatrix"),
                 c("dgCMatrix", "dsCMatrix"),
                 c("dgCMatrix", "dtCMatrix"),
                 c("dgCMatrix", "lgCMatrix"),
                 c("dgCMatrix", "ngCMatrix"),
                 c("dgTMatrix", "dgCMatrix"),
                 c("dgTMatrix", "dgeMatrix"),
                 c("dgTMatrix", "dsTMatrix"),
                 c("dgTMatrix", "dtCMatrix"),
                 c("dgTMatrix", "dtTMatrix"),
                 c("dgeMatrix", "dgCMatrix"),
                 c("dgeMatrix", "dgTMatrix"),
                 c("dgeMatrix", "dsTMatrix"),
                 c("dgeMatrix", "dspMatrix"),
                 c("dgeMatrix", "dsyMatrix"),
                 c("dgeMatrix", "dtrMatrix"),
                 c("dgeMatrix", "lgeMatrix"),
                 c("dpoMatrix", "dppMatrix"),
                 c("dppMatrix", "dpoMatrix"),
                 c("dsCMatrix", "dgCMatrix"),
                 c("dsCMatrix", "dgTMatrix"),
                 c("dsCMatrix", "dgeMatrix"),
                 c("dsCMatrix", "dsRMatrix"),
                 c("dsCMatrix", "dsTMatrix"),
                 c("dsCMatrix", "dsyMatrix"),
                 c("dsCMatrix", "lsCMatrix"),
                 c("dsCMatrix", "nsCMatrix"),
                 c("dsTMatrix", "dgTMatrix"),
                 c("dsTMatrix", "dgeMatrix"),
                 c("dsTMatrix", "dsCMatrix"),
                 c("dsTMatrix", "dsyMatrix"),
                 c("dsTMatrix", "lsTMatrix"),
                 c("dspMatrix", "dppMatrix"),
                 c("dspMatrix", "dsyMatrix"),
                 c("dspMatrix", "lspMatrix"),
                 c("dsyMatrix", "dpoMatrix"),
                 c("dsyMatrix", "dsCMatrix"),
                 c("dsyMatrix", "dsTMatrix"),
                 c("dsyMatrix", "dspMatrix"),
                 c("dsyMatrix", "lsyMatrix"),
                 c("dtCMatrix", "dgCMatrix"),
                 c("dtCMatrix", "dgTMatrix"),
                 c("dtCMatrix", "dgeMatrix"),
                 c("dtCMatrix", "dsCMatrix"),
                 c("dtCMatrix", "dtTMatrix"),
                 c("dtCMatrix", "dtrMatrix"),
                 c("dtCMatrix", "ltCMatrix"),
                 c("dtCMatrix", "ntCMatrix"),
                 c("dtTMatrix", "dgTMatrix"),
                 c("dtTMatrix", "dgeMatrix"),
                 c("dtTMatrix", "dtCMatrix"),
                 c("dtTMatrix", "dtrMatrix"),
                 c("dtpMatrix", "dtTMatrix"),
                 c("dtpMatrix", "dtrMatrix"),
                 c("dtpMatrix", "ltpMatrix"),
                 c("dtrMatrix", "dtpMatrix"),
                 c("dtrMatrix", "ltrMatrix"),
                 c("indMatrix", "ngTMatrix"),
                 c("indMatrix", "ngeMatrix"),
                 c("lMatrix", "dgCMatrix"),
                 c("lgCMatrix", "dgCMatrix"),
                 c("lgCMatrix", "lgTMatrix"),
                 c("lgCMatrix", "lgeMatrix"),
                 c("lgCMatrix", "ltCMatrix"),
                 c("lgTMatrix", "dgTMatrix"),
                 c("lgTMatrix", "lgCMatrix"),
                 c("lgTMatrix", "lgeMatrix"),
                 c("lgTMatrix", "lsCMatrix"),
                 c("lgTMatrix", "ltTMatrix"),
                 c("lgeMatrix", "dgeMatrix"),
                 c("lgeMatrix", "lgCMatrix"),
                 c("lgeMatrix", "lgTMatrix"),
                 c("lgeMatrix", "lspMatrix"),
                 c("lgeMatrix", "lsyMatrix"),
                 c("lgeMatrix", "ltpMatrix"),
                 c("lgeMatrix", "ltrMatrix"),
                 c("lsCMatrix", "dsCMatrix"),
                 c("lsCMatrix", "lgCMatrix"),
                 c("lsCMatrix", "lgTMatrix"),
                 c("lsCMatrix", "lsTMatrix"),
                 c("lsTMatrix", "lgCMatrix"),
                 c("lsTMatrix", "lgTMatrix"),
                 c("lsTMatrix", "lsCMatrix"),
                 c("lsTMatrix", "lsyMatrix"),
                 c("lspMatrix", "dspMatrix"),
                 c("lspMatrix", "lgeMatrix"),
                 c("lspMatrix", "lsyMatrix"),
                 c("lsyMatrix", "dsyMatrix"),
                 c("lsyMatrix", "lgeMatrix"),
                 c("lsyMatrix", "lspMatrix"),
                 c("ltCMatrix", "lgCMatrix"),
                 c("ltCMatrix", "ltTMatrix"),
                 c("ltTMatrix", "dtTMatrix"),
                 c("ltTMatrix", "lgCMatrix"),
                 c("ltTMatrix", "lgTMatrix"),
                 c("ltTMatrix", "ltCMatrix"),
                 c("ltTMatrix", "ltrMatrix"),
                 c("ltpMatrix", "dtpMatrix"),
                 c("ltpMatrix", "lgeMatrix"),
                 c("ltpMatrix", "ltrMatrix"),
                 c("ltrMatrix", "dtrMatrix"),
                 c("ltrMatrix", "lgeMatrix"),
                 c("ltrMatrix", "ltpMatrix"),
                 ## c("matrix", "dgCMatrix"),
                 c("matrix", "dgRMatrix"),
                 c("matrix", "dgTMatrix"),
                 c("matrix", "dgeMatrix"),
                 c("matrix", "dpoMatrix"),
                 c("matrix", "dppMatrix"),
                 c("matrix", "dsCMatrix"),
                 c("matrix", "dsTMatrix"),
                 c("matrix", "dspMatrix"),
                 c("matrix", "dsyMatrix"),
                 c("matrix", "dtCMatrix"),
                 c("matrix", "dtTMatrix"),
                 c("matrix", "dtpMatrix"),
                 c("matrix", "dtrMatrix"),
                 c("matrix", "lgCMatrix"),
                 c("matrix", "lgTMatrix"),
                 c("matrix", "lgeMatrix"),
                 c("matrix", "lsCMatrix"),
                 c("matrix", "lspMatrix"),
                 c("matrix", "lsyMatrix"),
                 c("matrix", "ltCMatrix"),
                 c("matrix", "ltTMatrix"),
                 c("matrix", "ltpMatrix"),
                 c("matrix", "ltrMatrix"),
                 c("matrix", "ngCMatrix"),
                 c("matrix", "ngTMatrix"),
                 c("matrix", "ngeMatrix"),
                 c("matrix", "nspMatrix"),
                 c("matrix", "nsyMatrix"),
                 c("matrix", "ntCMatrix"),
                 c("matrix", "ntTMatrix"),
                 c("matrix", "ntpMatrix"),
                 c("matrix", "ntrMatrix"),
                 ## c("matrix.coo", "dgCMatrix"),
                 ## c("matrix.coo", "dgTMatrix"),
                 ## c("matrix.csc", "dgCMatrix"),
                 ## c("matrix.csr", "dgCMatrix"),
                 ## c("matrix.csr", "dgRMatrix"),
                 c("ngCMatrix", "dgCMatrix"),
                 c("ngCMatrix", "lgCMatrix"),
                 c("ngCMatrix", "ntCMatrix"),
                 c("ngTMatrix", "dgTMatrix"),
                 c("ngTMatrix", "lgTMatrix"),
                 c("ngTMatrix", "lgeMatrix"),
                 c("ngTMatrix", "ngCMatrix"),
                 c("ngTMatrix", "ngeMatrix"),
                 c("ngTMatrix", "ntTMatrix"),
                 c("ngeMatrix", "dgeMatrix"),
                 c("ngeMatrix", "lgeMatrix"),
                 c("ngeMatrix", "ngCMatrix"),
                 c("ngeMatrix", "ngTMatrix"),
                 c("ngeMatrix", "nspMatrix"),
                 c("ngeMatrix", "nsyMatrix"),
                 c("ngeMatrix", "ntpMatrix"),
                 c("ngeMatrix", "ntrMatrix"),
                 c("nsCMatrix", "dsCMatrix"),
                 c("nsCMatrix", "lsCMatrix"),
                 c("nsCMatrix", "ngCMatrix"),
                 c("nsCMatrix", "nsTMatrix"),
                 c("nsTMatrix", "dsTMatrix"),
                 c("nsTMatrix", "ngCMatrix"),
                 c("nsTMatrix", "ngTMatrix"),
                 c("nsTMatrix", "nsCMatrix"),
                 c("nsTMatrix", "nsyMatrix"),
                 c("nspMatrix", "dspMatrix"),
                 c("nspMatrix", "lspMatrix"),
                 c("nspMatrix", "ngeMatrix"),
                 c("nspMatrix", "nsyMatrix"),
                 c("nsyMatrix", "dsyMatrix"),
                 c("nsyMatrix", "lsyMatrix"),
                 c("nsyMatrix", "ngeMatrix"),
                 c("nsyMatrix", "nspMatrix"),
                 c("ntCMatrix", "dtCMatrix"),
                 c("ntCMatrix", "ltCMatrix"),
                 c("ntCMatrix", "ngCMatrix"),
                 c("ntTMatrix", "dtTMatrix"),
                 c("ntTMatrix", "ngCMatrix"),
                 c("ntTMatrix", "ngTMatrix"),
                 c("ntTMatrix", "ntCMatrix"),
                 c("ntTMatrix", "ntrMatrix"),
                 c("ntpMatrix", "dtpMatrix"),
                 c("ntpMatrix", "ltpMatrix"),
                 c("ntpMatrix", "ngeMatrix"),
                 c("ntpMatrix", "ntrMatrix"),
                 c("ntrMatrix", "dtrMatrix"),
                 c("ntrMatrix", "ltrMatrix"),
                 c("ntrMatrix", "ngeMatrix"),
                 c("ntrMatrix", "ntpMatrix"),
                 c("numLike", "dgeMatrix"))

.def.template <-
function(from) {
    cd1 <- getClassDef(.FROM)
    cd2 <- getClassDef(.TO)
    Matrix.DeprecatedCoerce(cd1, cd2)
    to <- .CALL
    if(identical(as.character(class(to)), .TO))
        return(to)
    ## Coercion via virtual generated a _subclass_ of the target class
    to.strict <- new(.TO)
    for(nm in slotNames(cd2))
        slot(to.strict, nm) <- slot(to, nm)
    to.strict
}
for(.f.t in .from.to) {
    .f <- .f.t[1L]
    .t <- .f.t[2L]
    .def <- .def.template
    .env <- list(.FROM = .f, .TO = .t, .CALL = .as.via.virtual(.f, .t))
    body(.def) <- do.call(substitute, list(body(.def), .env))
    setAs(.f, .t, .def)
}
rm(.from.to, .f.t, .f, .t, .def.template, .def, .env)

setAs("sparseCholesky", "Matrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"Matrix\")",
                      new = "expand1(., \"L\")",
                      package = "Matrix")
          }
          expand1(from, "L")
      })

setAs("sparseCholesky", "dMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"dMatrix\")",
                      new = "expand1(., \"L\")",
                      package = "Matrix")
          }
          expand1(from, "L")
      })

setAs("sparseCholesky", "dsparseMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"dsparseMatrix\")",
                      new = "expand1(., \"L\")",
                      package = "Matrix")
          }
          expand1(from, "L")
      })

setAs("sparseCholesky", "sparseMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"sparseMatrix\")",
                      new = "expand1(., \"L\")",
                      package = "Matrix")
          }
          expand1(from, "L")
      })

setAs("sparseCholesky", "CsparseMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"CsparseMatrix\")",
                      new = "expand1(., \"L\")",
                      package = "Matrix")
          }
          expand1(from, "L")
      })

setAs("sparseCholesky", "RsparseMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"RsparseMatrix\")",
                      new = "as(expand1(., \"L\"), \"RsparseMatrix\")",
                      package = "Matrix")
          }
          as(expand1(from, "L"), "RsparseMatrix")
      })

setAs("sparseCholesky", "TsparseMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"TsparseMatrix\")",
                      new = "as(expand1(., \"L\"), \"TsparseMatrix\")",
                      package = "Matrix")
          }
          as(expand1(from, "L"), "TsparseMatrix")
      })

setAs("sparseCholesky", "triangularMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"triangularMatrix\")",
                      new = "as(expand1(., \"L\"), \"triangularMatrix\")",
                      package = "Matrix")
          }
          as(expand1(from, "L"), "triangularMatrix")
      })

setAs("sparseCholesky", "pMatrix",
      function(from) {
          if(FALSE) {
          .Deprecated(old = "as(<sparseCholesky>, \"pMatrix\")",
                      new = "expand1(., \"P1\")",
                      package = "Matrix")
          }
          expand1(from, "P1")
      })

setMethod("chol2inv", c(x = "sparseCholesky"),
          function(x, ...) {
              if(FALSE) {
              .Deprecated(old = "chol2inv(<sparseCholesky>)",
                          new = "solve(.)",
                          package = "Matrix")
              }
              solve(x)
          })


## ~~~~ DEFUNCT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cBind <- function(..., deparse.level = 1)
    .Defunct(new = "cbind", package = "Matrix")
rBind <- function(..., deparse.level = 1)
    .Defunct(msg = "rbind", package = "Matrix")
