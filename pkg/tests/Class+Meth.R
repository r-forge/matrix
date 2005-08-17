library(Matrix)

#### Automatically display the class inheritance structure
#### possibly augmented with methods

allCl <- getClasses("package:Matrix")

## Really nice would be to construct an inheritance graph and display
## it.  The following is just a cheap first step.

cat("All classes in the 'Matrix' package:\n")
for(cln in allCl) {
    cat("\n-----\n\nClass", dQuote(cln),":\n      ",
        paste(rep("~",nchar(cln)),collapse=''),"\n")
    ## A smarter version would use  getClass() instead of showClass(),
    ## build the "graph" and only then display.

    showClass(cln)
}

cat("\n\n")

## One could extend the `display' by using (something smarter than)
## are the "coerce" methods showing more than the 'Extends' output above?
cat("All (S4) methods in the 'Matrix' package:\n")
showMethods(where="package:Matrix")


### Show sparse logical matrices:
mM <- Matrix(1:4 >= 4, 2,2)
for(cl in getClass("lsparseMatrix")@subclasses) {
    clNam <- cl@subClass
    cat(clNam,"; new(..):")
    m <- new(clNam)
## FAILS e.g. for ''lgCMatrix'
    if(FALSE) ## FIXME
    cat("valid:", validObject(m), "as(*, dge*):")
    if(FALSE) { ## FIXME
        m2 <- as(mM, clNam)
        cat("valid:", validObject(m2))
    }
}
