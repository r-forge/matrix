## METHODS FOR CLASS: CsparseMatrix (virtual)
## sparse matrices in compressed sparse column (CSC) format
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.validateCsparse <- function(x, sort.if.needed = FALSE)
    .Call(Csparse_validate2, x, sort.if.needed)
##-> to be used in sparseMatrix(.), e.g. --- but is unused currently
## NB: 'sort.if.needed' is called 'maybe_modify' in C -- so be careful
## more useful:
.sortCsparse <- function(x) .Call(Csparse_sort, x) ## modifies 'x' !!

setMethod("writeMM", "CsparseMatrix",
	  function(obj, file, ...)
	  .Call(Csparse_MatrixMarket, obj, path.expand(as.character(file))))

dmperm <- function(x, nAns = 6L, seed = 0L) {
    stopifnot(length(nAns <- as.integer(nAns)) == 1L, nAns %in% c(2L, 4L, 6L),
              length(seed <- as.integer(seed)) == 1L, seed %in% -1:1)
    if(isS4(x)) {
        cld <- getClassDef(class(x))
        if(!extends(cld, "CsparseMatrix"))
            cld <- getClassDef(class(x <- as(x, "CsparseMatrix")))
        if(extends(cld, "symmetricMatrix"))
            cld <- getClassDef(class(x <- .M2gen(x)))
        if(!(extends(cld, "dMatrix") || extends(cld, "nMatrix")))
            x <- .M2kind(x, "d")
    } else { # typically a traditional matrix
        x <- .m2sparse(x, "dgC", NULL, NULL)
    }
    .Call(Csparse_dmperm, x, seed, nAns) # tolerating only [dn][gt]CMatrix 'x'
}
