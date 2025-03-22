## METHODS FOR GENERIC: Ops = {Arith, Compare, Logic} (group)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## > getGroupMembers("Ops")
## [1] "Arith"   "Compare" "Logic"
## > getGroupMembers("Arith")
## [1] "+"   "-"   "*"   "^"   "%%"  "%/%" "/"
## > getGroupMembers("Compare")
## [1] "==" ">"  "<"  "!=" "<=" ">="
## > getGroupMembers("Logic") # excluding unary "!" -> ./not.R
## [1] "&" "|"

.Ops.invalid <- function(x) {
    if(is.object(x))
        gettextf("invalid class \"%s\" in '%s' method", class(x)[1L], "Ops")
    else gettextf("invalid type \"%s\" in '%s' method", typeof(x), "Ops")
}


if(FALSE) {

## vvvv MJ: for _after_ 1.6-2, ditto ./(Arith|Compare|Logic).R
for(.cl in c("Matrix", "sparseVector")) {
setMethod("Ops", c(e1 = .cl, e2 = "ANY"),
          function(e1, e2)
              if(any(typeof(e2) == c("logical", "integer", "double"))) {
                  if(is.matrix(e2))
                      callGeneric(e1, unclass(e2))
                  else callGeneric(e2, as.vector(e2))
              } else stop(.Ops.invalid(e2), domain = NA))

setMethod("Ops", c(e1 = "ANY", e2 = .cl),
          function(e1, e2)
              if(any(typeof(e1) == c("logical", "integer", "double"))) {
                  if(is.matrix(e1))
                      callGeneric(unclass(e1), e2)
                  else callGeneric(as.vector(e1), e2)
              } else stop(.Ops.invalid(e1), domain = NA))

setMethod("Ops", c(e1 = .cl, e2 = "NULL"),
          function(e1, e2) callGeneric(e1, .NULL2vec(e1)))

setMethod("Ops", c(e1 = "NULL", e2 = .cl),
          function(e1, e2) callGeneric(.NULL2vec(e2), e2))

## MJ: OK, but I'd prefer to handle all "matrix" as ".geMatrix"
setMethod("Ops", c(e1 = .cl, e2 = "matrix"),
          function(e1, e2)
              if(any(typeof(e2) == c("logical", "integer", "double")))
                  callGeneric(e1, Matrix(e2))
              else stop(.Ops.invalid(e2), domain = NA))

## MJ: OK, but I'd prefer to handle all "matrix" as ".geMatrix"
setMethod("Ops", c(e1 = "matrix", e2 = .cl),
          function(e1, e2)
              if(any(typeof(e1) == c("logical", "integer", "double")))
                  callGeneric(Matrix(e1), e2)
              else stop(.Ops.invalid(e1), domain = NA))
}
rm(.cl)
## ^^^^ MJ: for _after_ 1.6-2, ditto ./(Arith|Compare|Logic).R
}

## "sparseVector" o "matrix"  is *not* dealt with correctly otherwise;
setMethod("Ops", c(e1 = "sparseVector", e2 = "matrix"),
          function(e1, e2)
              if(is.atomic(e2) && nzchar(kn <- .M.kind(e2)))
                  callGeneric(e1, as(e2, paste0(kn,"Matrix")))
              else stop(.Ops.invalid(e2), domain = NA))
setMethod("Ops", c(e1 = "matrix", e2 = "sparseVector"),
          function(e1, e2)
              if(is.atomic(e1) && nzchar(kn <- .M.kind(e1)))
                  callGeneric(as(e1, paste0(kn,"Matrix")), e2)
              else stop(.Ops.invalid(e1), domain = NA))



.Ops.checkDim <- function(d.a, d.b) {
    if(any(d.a != d.b))
        stop(gettextf("non-conformable matrix dimensions in %s",
                      deparse(sys.call(sys.parent()))),
             call. = FALSE, domain = NA)
    d.a
}

.Ops.checkDimNames <- function(dn.a, dn.b, useFirst = TRUE, check = FALSE) {
    ## behave as described in ?Arithmetic
    nullDN <- list(NULL, NULL)
    h.a <- !identical(nullDN, dn.a)
    h.b <- !identical(nullDN, dn.b)
    if(h.a || h.b) {
        if(useFirst) {
            if(!h.a) dn.b else dn.a
        } else {
            if (!h.b) dn.a
            else if(!h.a) dn.b
            else { ## both have non-trivial dimnames
                r <- dn.a # "default" result
                for(j in 1:2)
                    if(!is.null(dn <- dn.b[[j]])) {
                        if(is.null(r[[j]]))
                            r[[j]] <- dn
                        else if(check && !identical(r[[j]], dn))
                            warning(gettextf("dimnames [%d] mismatch in %s", j,
                                             deparse(sys.call(sys.parent()))),
                                    call. = FALSE, domain = NA)
                    }
                r
            }
        }
    } else nullDN
}

## cache them [rather in package 'methods' ??]
.ArithGenerics   <- getGroupMembers("Arith")
if(FALSE) { # unused
.CompareGenerics <- getGroupMembers("Compare")
.LogicGenerics   <- getGroupMembers("Logic")
}

## For now:
## MM: __FIXME__: treats vectors as matrices; not ok, e.g. for M + integer()
.fn1 <-
function(e1, e2)
    as(callGeneric(as(e1, "matrix"),    e2           ), "Matrix")
.fn2 <-
function(e1, e2)
    as(callGeneric(   e1           , as(e2, "matrix")), "Matrix")
.fn3 <-
function(e1, e2)
    as(callGeneric(as(e1, "matrix"), as(e2, "matrix")), "Matrix")
for(.cl in c("iMatrix", "zMatrix")) {
setMethod("Ops", c(e1 =      .cl, e2 =    "ANY"), .fn1)
setMethod("Ops", c(e1 =    "ANY", e2 =      .cl), .fn2)
setMethod("Ops", c(e1 =      .cl, e2 = "Matrix"), .fn3)
setMethod("Ops", c(e1 = "Matrix", e2 =      .cl), .fn3)
}
rm(.cl, .fn1, .fn2, .fn3)

### Design decision for *sparseMatrix*:
### work via Csparse  since Tsparse are not-unique (<-> slots not compatible)

### --  0 -- (not dense *or* sparse) -----------------------------------

##-------- originally from ./Matrix.R --------------------

## Some ``Univariate'' "Arith" (univariate := 2nd argument 'e2' is missing)
setMethod("+", c(e1 = "Matrix", e2 = "missing"), function(e1,e2) e1)

## -- unary `-` :  - <Matrix> ----------------------
## "fallback":
setMethod("-", c(e1 = "Matrix", e2 = "missing"),
          function(e1, e2) {
              warning("inefficient method used for \"- e1\"")
              0 - e1
          })
setMethod("-", c(e1 = "denseMatrix", e2 = "missing"),
          function(e1, e2) {
              e1@x <- -e1@x
              if(.hasSlot(e1, "factors") && length(e1@factors))
                  e1@factors <- list()
              e1
          })
## with these exceptions:
setMethod("-", c(e1 = "ndenseMatrix", e2 = "missing"), function(e1, e2) -.M2kind(e1, "i"))
setMethod("-", c(e1 = "ldenseMatrix", e2 = "missing"), function(e1, e2) -.M2kind(e1, "i"))

setMethod("-", c(e1 = "diagonalMatrix", e2 = "missing"),
          function(e1, e2) {
              kind <- .M.kind(e1)
              r <- new(if(kind == "z") "zdiMatrix" else "ddiMatrix")
              r@Dim <- d <- e1@Dim
              r@Dimnames <- e1@Dimnames
              r@x <-
                  if(e1@diag != "N")
                      rep.int(if(kind == "z") -1+0i else -1, d[1L])
                  else -(if(kind == "n") e1@x | is.na(e1@x) else e1@x)
              r
          })



## old-style matrices are made into new ones
setMethod("Ops", c(e1 = "Matrix", e2 = "matrix"),
          function(e1, e2) callGeneric(e1, Matrix(e2)))

setMethod("Ops", c(e1 = "matrix", e2 = "Matrix"),
          function(e1, e2) callGeneric(Matrix(e1), e2))
## Note: things like  callGeneric(Matrix(e1, sparse=is(e2,"sparseMatrix")), e2))
##   may *not* be better: e.g. Matrix(.) can give *diagonal* instead of sparse

##  Used to treat NULL as logical(0), now (+integer+complex) more careful:
setMethod("Ops", c(e1 = "Matrix", e2 = "NULL"),
          function(e1, e2) callGeneric(e1, .NULL2vec(e1)))
setMethod("Ops", c(e1 = "NULL", e2 = "Matrix"),
          function(e1, e2) callGeneric(.NULL2vec(e2), e2))

setMethod("Ops", c(e1 = "Matrix", e2 = "ANY"),
          function(e1, e2)
              if(is.object(e2) && is.matrix(e2))
                  callGeneric(e1, unclass(e2)) # e.g., for "table"
              else .bail.out.2(.Generic, class(e1), class(e2)))
setMethod("Ops", c(e1 = "ANY", e2 = "Matrix"),
          function(e1, e2)
              if(is.object(e1) && is.matrix(e1))
                  callGeneric(unclass(e1), e2) # e.g., for "table"
              else .bail.out.2(.Generic, class(e1), class(e2)))

## "General principle"
##  - - - - - - - - -
setMethod("Arith", c(e1 = "zMatrix", e2 = "Matrix"),
          function(e1, e2) callGeneric(e1, .M2kind(e2, "z")))
setMethod("Arith", c(e1 = "Matrix", e2 = "zMatrix"),
          function(e1, e2) callGeneric(.M2kind(e1, "z"), e2))
## use  iMatrix for {n,l,i} Arithmetic:
for(.cl in c("nMatrix", "lMatrix")) {
setMethod("Arith", c(e1 = "iMatrix", e2 = .cl),
          function(e1, e2) callGeneric(e1, .M2kind(e2, "i")))
setMethod("Arith", c(e1 = .cl, e2 = "iMatrix"),
          function(e1, e2) callGeneric(.M2kind(e1, "i"), e2))
}
rm(.cl)
## ==> Ops.x.x  for  <dM.>  o <dM.> etc

## if(FALSE)  <<<<< FIXME >>>>>>>>>
## just in case --- FIXME --should *not* be needed !!
## all other combinations
setMethod("Arith", c(e1 = "Matrix", e2 = "Matrix"),
          function(e1, e2) callGeneric(as(e1, "dMatrix"), as(e2, "dMatrix")))

## For "Compare", this would be feasible too, but is clearly suboptimal,
## particularly for "==" and "!="
## and for "lMatrix" and "nMatrix"  should not coerce at all
if(FALSE)
setMethod("Compare", c(e1 = "Matrix", e2 = "Matrix"),
          function(e1, e2) {
              if(is.na(match(.Generic, c("==", "!="))))
                  callGeneric(as(e1, "dMatrix"), as(e2, "dMatrix"))
              else { ## no coercion needed for "==" or "!="
                  ##
                  ## what now ?  <<<<<<<<<<< FIXME >>>>>>>>>
                  .bail.out.2(.Generic, class(e1), class(e2))
              }
          })


## Working entirely on "matching" x slot:
## can be done for matching-dim "*geMatrix", and also
## matching-{dim + uplo} for *packed* (only!) symmetric+triangular
.Ops.via.x <- function(e1, e2) {
    .Ops.checkDim(dim(e1), dim(e2))
    e1@x <- callGeneric(e1@x, e2@x)
    if(.hasSlot(e1, "factors") && length(e1@factors))
        e1@factors <- list()
    e1
}


###-------- originally from ./dMatrix.R --------------------
##
## Note that there extra methods for <sparse> o <sparse> !
##
## "Compare" -> returning  logical Matrices;  .Cmp.swap() is in ./Auxiliaries.R
setMethod("Compare", c(e1 = "numeric", e2 = "dMatrix"), .Cmp.swap)
setMethod("Compare", c(e1 = "logical", e2 = "dMatrix"), .Cmp.swap)
setMethod("Compare", c(e1 = "numeric", e2 = "lMatrix"), .Cmp.swap)
setMethod("Compare", c(e1 = "logical", e2 = "lMatrix"), .Cmp.swap)
setMethod("Compare", c(e1 = "numeric", e2 = "nMatrix"), .Cmp.swap)
setMethod("Compare", c(e1 = "logical", e2 = "nMatrix"), .Cmp.swap)


## This is parallel to Logic.Mat.atomic() below --->  __keep parallel__ !
Cmp.Mat.atomic <- function(e1, e2) { ## result will inherit from "lMatrix"
    n1 <- prod(d <- e1@Dim)
    cl <- class(e1)
    if(!(l2 <- length(e2)))
        return(if(n1 == 0) .M2kind(e1, "l") else as.logical(e2))
    ## else
    if(n1 && n1 < l2)
        stop(sprintf("dim [product %d] do not match the length of object [%d]",
                     n1, l2))
    cl1 <- getClassDef(cl)
    slots1 <- names(cl1@slots)
    has.x <- any("x" == slots1)# *fast* check for "x" slot presence
    if(l2 > 1 && has.x)
        return(if(n1 == 0) {
                   ## use `slots1`   sNms <-  .slotNames(e1)
                   r <- copyClass(e1, class2(cl, "l"),
                                  sNames = slots1[slots1 != "x"], check = FALSE)
                   r@x <- callGeneric(e1@x, e2)
                   r
               } else # cannot simply compare e2, e1@x -> use another method
                   callGeneric(e1, Matrix(e2, nrow=d[1], ncol=d[2]))
               )
    ## else
    isDiag <- extends(cl1, "diagonalMatrix")
    r0 <- callGeneric(0, e2)
    ar0 <- all(r0)
    if(isDiag && ar0) # all off-diagonals became TRUE ==> dense symmetric
        return(callGeneric(.diag2dense(e1, shape = "s"), e2))
    r <- callGeneric(if(has.x) e1@x else TRUE, e2)
    isTri <- extends(cl1, "triangularMatrix")
    ## Udg: append the diagonal at *end*, similar to diagU2N():
    Udg <- (isTri || isDiag) && e1@diag == "U"
    r. <- if(Udg) c(r, d.r <- callGeneric(..diag.x(e1), e2)) else r
    ## trivial case first (beware of NA)
    if(isTRUE(ar0 && all(r.))) { # |--> all TRUE
        r <- new(if(d[1] == d[2]) "lsyMatrix" else "lgeMatrix")
        r@Dim <- d
        r@Dimnames <- e1@Dimnames
        r@x <- rep.int(TRUE, n1)
    } else if(extends(cl1, "denseMatrix")) {
        full <- !.isPacked(e1) # << both "dtr" and "dsy" are 'full'
        if(full || allFalse(r0) || extends(cl1, "symmetricMatrix")) {
            if(sparseDefault(r, tol = 0.25)) # switch to sparse result
                return( callGeneric(as(e1, "sparseMatrix"), e2) )
            ## else
            ## FIXME? using copyClass() to copy "relevant" slots
            r <- new(class2(cl, "l"), x = r, Dim = d, Dimnames = e1@Dimnames)
            if(extends(cl1, "symmetricMatrix")) {
                r@uplo <- e1@uplo
            } else if(isTri) {
                r@uplo <- e1@uplo
                r@diag <- e1@diag
            }
        } else {
            ## packed matrix with structural 0 and r0 is not all FALSE:
            ##--> result cannot be packed anymore
            ## [dense & packed & not symmetric ] ==> must be "dtp*" :
            if(!extends(cl1, "dtpMatrix"))
                stop("internal bug in \"Compare\" method (Cmp.Mat.atomic); please report")
            rx <- rep_len(r0, n1)
            rx[indTri(d[1], upper = (e1@uplo == "U"), diag=TRUE)] <- r.
            r <- new("lgeMatrix", x = rx, Dim = d, Dimnames = e1@Dimnames)
        }
    } else if(isDiag) { ## && !ar0 := !all(r0)
        if(Udg)
            e1 <- .diagU2N(e1, cl = cl1, checkDense = TRUE)
        tx <- typeof(r.)
        if(tx != typeof(e1@x))
            e1 <- .M2kind(e1, .kind.type[[tx]])
        e1@x <- if(Udg) d.r else r
        return(if(Udg) .diagN2U(e1) else e1)
    } else {
        ##---- e1 is(. , sparseMatrix) -----------------
        remainSparse <- allFalse(r0) ## <==> things remain sparse
        if(Udg) {          # e1 *is* unit-diagonal (triangular sparse)
            r1 <- callGeneric(1, e2)
            Udg <- all(r1)       # maybe Unit-diagonal (sparse) result
            ## if(!remainSparse) we'll use non0ind() which *has* unit-diag. indices at end
            ##
            if(Udg && remainSparse) {
            } else { ## result will not be unit-diagonal sparse
                e1 <- .diagU2N(e1, cl = cl1) # otherwise, result is U-diag
                ## FIXME? rather
                ## if(extends1of(cl1, c("CsparseMatrix", "RsparseMatrix","TsparseMatrix")) {
                if(extends(cl1, "CsparseMatrix")) {
                    ## repeat computation if e1 has changed
                    r. <- callGeneric(if(has.x) e1@x else TRUE, e2)
                }
            }
        }

        if(remainSparse) {
            if(!anyNA(r) && ((Ar <- all(r)) || !any(r))) {
                lClass <- class2(cl, "l") # is "lsparse*"
                r <- new(lClass)
                r@Dim <- d
                r@Dimnames <- e1@Dimnames
                if(Ar) {       # 'TRUE' instead of 'x': same sparsity:
                    for(n in intersect(c("i","j","p","uplo","diag"), slots1))
                        slot(r, n) <- slot(e1, n)
                    n <- if(has.x)
                             length(e1@x)
                         else if(any("p" == slots1))
                             e1@p[d[2]+1L]
                         else length(e1@i)
                    r@x <- rep.int(TRUE, n)
                } else {
                    ## !any(r): all FALSE: keep empty 'r' matrix
                    ## but may need a valid 'pointer' slot:
                    if(extends(lClass, "CsparseMatrix"))
                        r@p <- rep.int(0L, 1+ncol(r))
                    else if(extends(lClass, "RsparseMatrix"))
                        r@p <- rep.int(0L, 1+nrow(r))
                }
            } else {  # some TRUE, FALSE, NA : go via unique 'Tsparse'
                M <- asUniqueT(e1)
                nCl <- class2(class(M), 'l') # logical Tsparse
                sN <- slotNames(nCl)
                ## copy "the other slots" (important for "tr"/"sym"):
                r <- copyClass(M, nCl, sNames = sN[is.na(match(sN, c("x", "factors")))])
                r@x <- callGeneric(if(has.x) M@x else 1, e2)
                if(extends(cl1, "CsparseMatrix"))
                    r <- .M2C(r)
                else if(extends(cl1, "RsparseMatrix"))
                    r <- .M2R(r)
            }
        } else {
            ## non sparse result; triangularity also gone, typically
            lClass <- if(extends(cl1, "symmetricMatrix"))
                          "lsyMatrix"
                      else "lgeMatrix"
            Matrix.message(sprintf("sparse to dense (%s) coercion in '%s' -> %s",
                               lClass, .Generic, "Cmp.Mat.atomic"),
                       .M.level = 2)
            rx <- rep_len(r0, n1)
            ## Here, we assume that 'r.' and the indices align (!)
            encI <- .Call(m_encodeInd,
                          non0ind(e1, cl1, uniqT=FALSE, xtendSymm=FALSE),
                          di = d, orig1=FALSE, checkBounds=FALSE)
            rx[1L + encI] <- r.
            r <- new(lClass, x = rx, Dim = d, Dimnames = e1@Dimnames)
        }
    }
    r
} # Cmp.Mat.atomic

## "Compare" -> returning  logical Matrices;
for(Mcl in c("dMatrix", "lMatrix", "nMatrix", "zMatrix", "iMatrix"))
    for(cl in c("numeric", "logical", "complex")) {
        setMethod("Compare", c(e1 = cl, e2 = Mcl), .Cmp.swap)# .Cmp.swap() --> ./Auxiliaries.R
        setMethod("Compare", c(e1 = Mcl, e2 = cl),  Cmp.Mat.atomic)
    }
rm(Mcl, cl)
rm(Cmp.Mat.atomic)

## "xMatrix <-> work with 'x' slot {was originally just for "Compare"}:
##  -------  {also used for "Arith"}:
Ops.x.x <- function(e1, e2)
{
    d <- .Ops.checkDim(dim(e1), dim(e2))
    if((dens1 <- extends(c1 <- class(e1), "denseMatrix")))
        gen1 <- extends(c1, "generalMatrix")
    if((dens2 <- extends(c2 <- class(e2), "denseMatrix")))
        gen2 <- extends(c2, "generalMatrix")
    if(dens1 && dens2) { ## both inherit from ddense*
        geM <- TRUE
        if(!gen1) {
            if(!gen2) { ## consider preserving "triangular" / "symmetric"
                geM <- FALSE
                le <- prod(d)
                isPacked <- function(x) length(x@x) < le
                Mclass <-
                    if(extends(c1, "symmetricMatrix") &&
                       extends(c2, "symmetricMatrix")) {
                      if(.M.kind(e1) != "z" || .M.kind(e2) != "z" || e1@trans == e2@trans) {
                        if(e1@uplo != e2@uplo)
                            ## one is upper, one is lower
                            e2 <- t(e2)
                        if((p1 <- isPacked(e1)) | (p2 <- isPacked(e2))) {
                            ## at least one is packed
                            if(p1 != p2) {
                                ## one is not packed --> *do* pack it:
                                if(p1)
                                    e2 <- pack(e2)
                                else
                                    e1 <- pack(e1)
                            }
                            "spMatrix"
                        } else "syMatrix"
                      }
                      else ## symmetric complex w/ differing `trans` slot
                          geM <- TRUE
                    } else if(extends(c1, "triangularMatrix") &&
                              extends(c2, "triangularMatrix")) {
                        geM <- e1@uplo != e2@uplo || isN0(callGeneric(0,0))
                        if(!geM) {
                            if(e1@diag == "U")
                                e1 <- ..diagU2N(e1)
                            if(e2@diag == "U")
                                e2 <- ..diagU2N(e2)
                            p1 <- isPacked(e1)
                            p2 <- isPacked(e2)
                            if(p1 || p2) { ## at least one is packed
                                if(p1 != p2) {
                                    ## one is not packed --> *do* pack it:
                                    if(p1) e2 <- pack(e2)
                                    else   e1 <- pack(e1)
                                }
                                "tpMatrix"
                            } else "trMatrix"
                        }
                    } else {
                        ## not symmetric, not triangular  ==> "general"
                        geM <- TRUE
                    }
                if(geM)
                    e2 <- .M2gen(e2)
            }
            if(geM)
                e1 <- .M2gen(e1) # was "dgeMatrix"
        } else { ## gen1
            if(!gen2)
                e2 <- .M2gen(e2)
        }
        ## now, in all cases @x should be matching & correct
        ## {only "uplo" part is used}
        r <- callGeneric(e1@x, e2@x)
        if(is.integer(r)) ## as "igeMatrix" does not yet exist!
            r <- as.double(r)
        kr <- .M.kind(r)
        if(geM)
            new(paste0(kr, "geMatrix"), x = r, Dim = d, Dimnames = e1@Dimnames)
        else
            new(paste0(kr, Mclass), x = r, Dim = d, Dimnames = e1@Dimnames,
                uplo = e1@uplo)
    } else { # at least one is sparse
        ## FIXME: e.g. if e1 is diagonal & e2 triangular / symmetric
        ##  coerce e1 to "like-e2" (in case of symmetric + z: w/ correct `trans` = ..
        if(!dens1 && .isDiagonal(e1)) # preserve symmetry
            e1 <- .diag2sparse(e1, shape = "s")
        if(!dens2 && .isDiagonal(e2))
            e2 <- .diag2sparse(e2, shape = "s")
        r <- if(!dens1 && !dens2) # _both_ are sparse
                 ## Now (new method dispatch, 2009-01) *does* happen
                 ## even though we have <sparse> o <sparse> methods
                 callGeneric(.M2C(e1), .M2C(e2))
             else if(dens1 && !dens2)
                 ## go to dense
                 callGeneric(e1, as(e2, "denseMatrix"))
             else # if(!dens1 && dens2)
                 callGeneric(as(e1, "denseMatrix"), e2)
        if(!.isSparse(r) && sparseDefault(r))
            as(r, "sparseMatrix")
        else r
    }
} # Ops.x.x

setMethod("Ops", c(e1 = "dMatrix", e2 = "dMatrix"), Ops.x.x)
setMethod("Ops", c(e1 = "lMatrix", e2 = "lMatrix"), Ops.x.x)
## n*: for "Arith" go via dMatrix, for "Logic" via "lMatrix"
setMethod("Compare", c(e1 = "nMatrix", e2 = "nMatrix"), Ops.x.x)
setMethod("Ops", c(e1 = "iMatrix", e2 = "iMatrix"), Ops.x.x)
setMethod("Ops", c(e1 = "zMatrix", e2 = "zMatrix"), Ops.x.x)

## l o d : depends on *kind* of Ops -- but Ops.x.x works on slots - correctly:
setMethod("Ops", c(e1="lMatrix", e2="dMatrix"), Ops.x.x)
setMethod("Ops", c(e1="dMatrix", e2="lMatrix"), Ops.x.x)

## lMatrix & nMatrix ... probably should also just use "Matrix" ?
##
## Hmm, the coercion should differ, depending on subgroup ("Logic", "Arith",..)
## --> try to get rid of these
## if(FALSE) {
setMethod("Ops", c(e1="lMatrix", e2="numeric"),
          function(e1, e2) callGeneric(as(e1,"dMatrix"), e2))
setMethod("Ops", c(e1="numeric", e2="lMatrix"),
          function(e1, e2) callGeneric(e1, as(e2,"dMatrix")))
setMethod("Ops", c(e1="nMatrix", e2="numeric"),
          function(e1, e2) callGeneric(as(e1,"dMatrix"), e2))
setMethod("Ops", c(e1="numeric", e2="nMatrix"),
          function(e1, e2) callGeneric(e1, as(e2,"dMatrix")))
## setMethod("Ops", c(e1="Matrix", e2="logical"),
##           function(e1,e2) callGeneric(as(e1,"lMatrix"), e2))
## setMethod("Ops", c(e1="logical", e2="Matrix"),
##           function(e1,e2) callGeneric(e1, as(e2,"lMatrix")))
## } # FALSE

### _______ FIXME ______  extend to  "posdefMatrix": all 10 of  [dz]p[opCRT]Matrix

## "dpoMatrix" / "dppMatrix" :
## Positive-definiteness is lost with all "Ops" but some "Arith" cases
for(cl in c("numeric", "logical")) { # "complex", "raw" : basically "replValue"

setMethod("Arith", c(e1 = cl, e2 = "dpoMatrix"),
          function(e1, e2) {
              if(!(l1 <- length(e1)))
                  double(0L)
              else if(l1 == 1 && any(.Generic == c("*","/","+")) && (e1 > 0)) {
                  e2@x <- callGeneric(e1, e2@x)
                  if(length(e2@factors))
                      e2@factors <- list()
                  e2 # remains "dpo"
              } else
                  callGeneric(e1, as(e2, "dsyMatrix"))
          })
setMethod("Arith", c(e1 = cl, e2 = "dppMatrix"),
          function(e1, e2) {
              if(!(l1 <- length(e1)))
                  double(0L)
              else if(l1 == 1 && any(.Generic == c("*","/","+")) && (e1 > 0)) {
                  e2@x <- callGeneric(e1, e2@x)
                  if(length(e2@factors))
                      e2@factors <- list()
                  e2 # remains "dpp"
              } else
                  callGeneric(e1, as(e2, "dspMatrix"))
          })
setMethod("Arith", c(e1 = "dpoMatrix", e2 = cl),
          function(e1, e2) {
              if(!(l2 <- length(e2)))
                  double(0L)
              else if(l2 == 1 && any(.Generic == c("*","/","+")) && (e2 > 0)) {
                  e1@x <- callGeneric(e1@x, e2)
                  if(length(e1@factors))
                      e1@factors <- list()
                  e1 # remains "dpo"
              } else
                  callGeneric(as(e1, "dsyMatrix"), e2)
          })
setMethod("Arith", c(e1 = "dppMatrix", e2 = cl),
          function(e1, e2) {
              if(!(l2 <- length(e2)))
                  double(0L)
              else if(l2 == 1 && any(.Generic == c("*","/","+")) && (e2 > 0)) {
                  e1@x <- callGeneric(e1@x, e2)
                  if(length(e1@factors))
                      e1@factors <- list()
                  e1 # remains "dpp"
              } else
                  callGeneric(as(e1, "dspMatrix"), e2)
          })

setMethod("Ops", c(e1 = cl, e2 = "dpoMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "dsyMatrix")))
setMethod("Ops", c(e1 = cl, e2 = "dppMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "dspMatrix")))
setMethod("Ops", c(e1 = "dpoMatrix", e2 = cl),
          function(e1, e2) callGeneric(as(e1, "dsyMatrix"), e2))
setMethod("Ops", c(e1 = "dppMatrix", e2 = cl),
          function(e1, e2) callGeneric(as(e1, "dspMatrix"), e2))
}# for(cl...)
rm(cl)

### --  I -- dense -----------------------------------------------------------

##-------- originally from ./dgeMatrix.R --------------------

for(Mcl in c(paste0(c("d","l","n", "i", "z"), "geMatrix")))
    setMethod("Arith", c(e1 = Mcl, e2 = Mcl),
              function(e1, e2) { ## M o M, M: [dlniz]geMatrix
                  ## NB:  triangular, symmetric, etc may need own method
                  .Ops.checkDim(d <- e1@Dim, e2@Dim)# now enforced by R's <matrix> too
                  dn <- .Ops.checkDimNames(dimnames(e1), dimnames(e2))
                  x <- callGeneric(e1@x, e2@x)
                  new(paste0(.kind.type[[typeof(x)]], "geMatrix"),
                      Dim = d, Dimnames = dn, x = x)
              })


A.M.n <- function(e1, e2) { ## M : [dlniz]geMatrix
    if(!is.null(dim(e2)))
        return(callGeneric(e1, .m2dense(e2)))
    d <- e1@Dim
    le <- length(e2)
    if(le == 0) {
        if(prod(d) == 0)
            new(class2(class(e1), .kind.type[[typeof(e2)]]), Dim = d, Dimnames = e1@Dimnames)
        else
            as.vector(e2)
    } else if(le == 1L || le == d[1] || any(prod(d) == c(le, 0))) { # matching dim
        tx <- typeof(x <- callGeneric(e1@x, as.vector(e2)))
        if(tx != typeof(e1@x))
            e1 <- .M2kind(e1, .kind.type[[tx]])
        e1@x <- x
        if(length(e1@factors))
            e1@factors <- list()
        e1
    } else
        stop("length of 2nd arg does not match dimension of first")
}

A.n.M <- function(e1, e2) { ## M : [dlniz]geMatrix
    if(!is.null(dim(e1)))
        return(callGeneric(.m2dense(e1), e2))
    d <- e2@Dim
    le <- length(e1)
    if(le == 0) {
        if(prod(d) == 0)
            new(class2(class(e2), .kind.type[[typeof(e1)]]), Dim = d, Dimnames = e2@Dimnames)
        else
            as.vector(e1)
    } else if(le == 1L || le == d[1] || any(prod(d) == c(le, 0))) { # matching dim
        tx <- typeof(x <- callGeneric(as.vector(e1), e2@x))
        if(tx != typeof(e2@x))
            e2 <- .M2kind(e2, .kind.type[[tx]])
        e2@x <- x
        if(length(e2@factors))
            e2@factors <- list()
        e2
    } else
        stop("length of 1st arg does not match dimension of 2nd")
}

for(cl in c("numeric", "logical", "complex", "sparseVector"))
    for(Mcl in c("dgeMatrix", "lgeMatrix", "ngeMatrix", "zgeMatrix", "igeMatrix")) {
        setMethod("Arith", c(e1 =  cl, e2 = Mcl), A.n.M)
        setMethod("Arith", c(e1 = Mcl, e2 =  cl), A.M.n)
    }
rm(cl, Mcl)
##
rm(A.M.n, A.n.M)

##-------- originally from ./ddenseMatrix.R --------------------

.Arith.denseM.atom <- function(e1, e2) {
    ## since e1 = "dgeMatrix" has its own method, we have
    ## either symmetric or triangular !
    n1 <- prod(d <- e1@Dim)
    le <- length(e2 <- as.vector(e2))
    if(n1 && n1 < le)
        stop(sprintf("dim [product %d] do not match the length of object [%d]",
                     n1, le))
    if(le == 0) {
        if(prod(d) == 0)
            new(class2(class(e1), .kind.type[[typeof(e2)]]), Dim = d, Dimnames = e1@Dimnames)
        else
            e2 # as.v.. above
    } else if(le == 1L || le == d[1] || any(prod(d) == c(le, 0))) { # matching dim
        if(is(e1, "triangularMatrix")) {
            r0 <- callGeneric(0, e2)
            if(all0(r0)) {              # result remains triangular
                if(e1@diag == "U" && !all(1 == callGeneric(1,e2)))
                    e1 <- diagU2N(e1)
                tx <- typeof(x <- callGeneric(e1@x, e2))
                if(tx != typeof(e1@x))
                    e1 <- .M2kind(e1, .kind.type[[tx]])
                e1@x <- x
                e1
            } else { ## result *general*
                callGeneric(.M2gen(e1), e2)
            }
        } else {                    ## symmetric
            if(le == 1) {           ## result remains symmetric
                tx <- typeof(x <- callGeneric(e1@x, e2))
                if(tx != typeof(e1@x))
                    e1 <- .M2kind(e1, .kind.type[[tx]])
                e1@x <- x
                if(length(e1@factors))
                    e1@factors <- list()
                e1
            } else { ## (le == d[1] || prod(d) == le)
                ## *might* remain symmetric, but 'x' may contain garbage
                ## *testing* for symmetry is also a bit expensive ==> simple:
                callGeneric(.M2gen(e1), e2)
            }
        }
    } else
        stop("length of 2nd arg does not match dimension of first")
}
## was only for "ddenseMatrix", now for all [zdiln]denseMatrix" :
setMethod("Arith", c(e1 = "denseMatrix", e2 =      "numeric"), .Arith.denseM.atom)
setMethod("Arith", c(e1 = "denseMatrix", e2 =      "logical"), .Arith.denseM.atom)
setMethod("Arith", c(e1 = "denseMatrix", e2 = "sparseVector"), .Arith.denseM.atom)
rm(.Arith.denseM.atom)

.Arith.atom.denseM <- function(e1, e2) {
    d <- e2@Dim
    ## note that e2 is either symmetric or triangular here
    le <- length(e1 <- as.vector(e1))
    if(le == 0) {
        if(prod(d) == 0)
            new(class2(class(e2), .kind.type[[typeof(e1)]]), Dim = d, Dimnames = e2@Dimnames)
        else
            e1 # as.v.. above
    } else if(le == 1L || le == d[1] || any(prod(d) == c(le, 0))) { # matching dim
        if(is(e2, "triangularMatrix")) {
            r0 <- callGeneric(e1, 0)
            if(all0(r0)) {		# result remains triangular
                if(e2@diag == "U" && !all(1 == callGeneric(e1,1)))
                    e2 <- diagU2N(e2)
                tx <- typeof(x <- callGeneric(e1, e2@x))
                if(tx != typeof(e2@x))
                    e2 <- .M2kind(e2, .kind.type[[tx]])
                e2@x <- x
                e2
            } else {			# result *general*
                callGeneric(e1, .M2gen(e2))
            }
        } else { ## symmetric
            if(le == 1) {		# result remains symmetric
                tx <- typeof(x <- callGeneric(e1, e2@x))
                if(tx != typeof(e2@x))
                    e2 <- .M2kind(e2, .kind.type[[tx]])
                e2@x <- x
                if(length(e2@factors))
                    e2@factors <- list()
                e2
            } else { ## (le == d[1] || prod(d) == le)
                ## *might* remain symmetric, but 'x' may contain garbage
                ## *testing* for symmetry is also a bit expensive ==> simple:
                callGeneric(e1, .M2gen(e2))
            }
        }
    } else
        stop("length of 1st arg does not match dimension of 2nd")
}
setMethod("Arith", c(e1 =      "numeric", e2 = "denseMatrix"), .Arith.atom.denseM)
setMethod("Arith", c(e1 =      "logical", e2 = "denseMatrix"), .Arith.atom.denseM)
setMethod("Arith", c(e1 = "sparseVector", e2 = "denseMatrix"), .Arith.atom.denseM)
rm(.Arith.atom.denseM)

## "Logic"
## -------

##-------- originally from ./ldenseMatrix.R --------------------

## These all had "Logic", now also for "Compare",
## but "Arith" differs: result will be "dgeMatrix' :
.Ops2dge.via.x <- function(e1,e2) {
    .Ops.checkDim(dim(e1), dim(e2))
    r <- copyClass(e1, "dgeMatrix", sNames = c("Dim","Dimnames"))
    r@x <- as.numeric(callGeneric(e1@x, e2@x))
    r
}

setMethod("Compare", c(e1="lgeMatrix", e2="lgeMatrix"), .Ops.via.x)
setMethod("Logic",   c(e1="lgeMatrix", e2="lgeMatrix"), .Ops.via.x)
setMethod("Arith",   c(e1="lgeMatrix", e2="lgeMatrix"), .Ops2dge.via.x)

setMethod("Compare", c(e1="ngeMatrix", e2="ngeMatrix"), .Ops.via.x)
setMethod("Logic",   c(e1="ngeMatrix", e2="ngeMatrix"), .Ops.via.x)
setMethod("Arith",   c(e1="ngeMatrix", e2="ngeMatrix"), .Ops2dge.via.x)

rm(.Ops.via.x, .Ops2dge.via.x)

## nMatrix -> lMatrix  conversions when "the other" is not nMatrix
## Use Ops.x.x unless both are sparse
setMethod("Ops", c(e1="dMatrix", e2="dMatrix"), Ops.x.x)
setMethod("Ops", c(e1="lMatrix", e2="lMatrix"), Ops.x.x)
setMethod("Ops", c(e1="nMatrix", e2="nMatrix"), Ops.x.x)
setMethod("Ops", c(e1="nMatrix", e2="lMatrix"), Ops.x.x)
setMethod("Ops", c(e1="lMatrix", e2="nMatrix"), Ops.x.x)
setMethod("Ops", c(e1="nMatrix", e2="dMatrix"), Ops.x.x)
setMethod("Ops", c(e1="dMatrix", e2="nMatrix"), Ops.x.x)
### __FIXME__  + zMatrix , iMatrix
rm(Ops.x.x)

## ... both are sparse: cannot use Ops.x.x
setMethod("Ops", c(e1="nsparseMatrix", e2="lsparseMatrix"),
          function(e1,e2) callGeneric(.M2kind(e1, "l"), e2))
setMethod("Ops", c(e1="lsparseMatrix", e2="nsparseMatrix"),
          function(e1,e2) callGeneric(e1, .M2kind(e2, "l")))
setMethod("Ops", c(e1="nsparseMatrix", e2="dsparseMatrix"),
          function(e1,e2) callGeneric(.M2kind(e1, "l"), e2))
setMethod("Ops", c(e1="dsparseMatrix", e2="nsparseMatrix"),
          function(e1,e2) callGeneric(e1, .M2kind(e2, "l")))
## For "Arith"  go to "d*", not "l*": {the above, replaced "l by "d : __FIXME?__ s/d/i/ :
setMethod("Arith", c(e1="nsparseMatrix", e2="lsparseMatrix"),
          function(e1,e2) callGeneric(.M2kind(e1, "d"), e2))
setMethod("Arith", c(e1="lsparseMatrix", e2="nsparseMatrix"),
          function(e1,e2) callGeneric(e1, .M2kind(e2, "d")))
setMethod("Arith", c(e1="nsparseMatrix", e2="dsparseMatrix"),
          function(e1,e2) callGeneric(.M2kind(e1, "d"), e2))
setMethod("Arith", c(e1="dsparseMatrix", e2="nsparseMatrix"),
          function(e1,e2) callGeneric(e1, .M2kind(e2, "d")))

if(FALSE) { ##-- not yet ---------
## New: for both "nsparseMatrix", *preserve* nsparse* -- via Tsp -- "nTsparseMatrix"
setMethod("Ops", c(e1 = "nsparseMatrix", e2 = "nsparseMatrix"),
          function(e1, e2) callGeneric(as(e1, "TsparseMatrix"),
                                       as(e2, "TsparseMatrix")))

Ops.nT.nT <- function(e1,e2) {
    d <- .Ops.checkDim(dim(e1), dim(e2))
    ## e1, e2 are nTsparse, i.e., inheriting from  "ngTMatrix", "ntTMatrix", "nsTMatrix"
    gen1 <- extends(cD1 <- getClassDef(class(e1)), "generalMatrix")
    gen2 <- extends(cD2 <- getClassDef(class(e2)), "generalMatrix")
    sym1 <- !gen1 && extends(cD1, "symmetricMatrix")
    sym2 <- !gen2 && extends(cD2, "symmetricMatrix")
    tri1 <- !gen1 && !sym1
    tri2 <- !gen2 && !sym2
    G <- gen1 && gen2
    S <- sym1 && sym2 && e1@uplo == e2@uplo
    T <- tri1 && tri2 && e1@uplo == e2@uplo
}

setMethod("Ops", c(e1="nTsparseMatrix", e2="nTsparseMatrix"), Ops.nT.nT)
}##--- not yet -------------



## Have this for "Ops" already above
## setMethod("Logic", c(e1 = "logical", e2 = "Matrix"),
##           function(e1, e2) callGeneric(e1, as(e2, "lMatrix")))
## setMethod("Logic", c(e1 = "Matrix", e2 = "logical"),
##           function(e1, e2) callGeneric(as(e1, "lMatrix"), e2))

.ll <- function(e1, e2) callGeneric(.M2kind(e1,"l"), .M2kind(e2, "l"))
setMethod("Logic", c(e1 = "nMatrix", e2 =  "Matrix"), .ll)
setMethod("Logic", c(e1 =  "Matrix", e2 = "nMatrix"), .ll)
setMethod("Logic", c(e1 = "nMatrix", e2 = "nMatrix"), .ll)
rm(.ll)

### "ANY" here means "any non-Matrix" (since "Ops"(ANY) has already bailout above):
setMethod("Logic", c(e1 = "ANY", e2 = "Matrix"),
          function(e1, e2) callGeneric(as.logical(e1), .M2kind(e2, "l")))
setMethod("Logic", c(e1 = "Matrix", e2 = "ANY"),
          function(e1, e2) callGeneric(as(e1, "lMatrix"), as.logical(e2)))

## "swap RHS and LHS" and use the method below -- can do this, since
## 	"Logic" := { "&" , "|" } and both are commutative
for(Mcl in c("lMatrix","nMatrix","dMatrix", "iMatrix", "zMatrix"))
    for(cl in c("logical", "numeric", "complex", "sparseVector"))
    setMethod("Logic", c(e1 = cl, e2 = Mcl),
              function(e1,e2) callGeneric(e2, e1))
## conceivably "numeric" could use  callGeneric(e2, as.logical(e1))
## but that's not useful at the moment, since Logic.Mat.atomic() does as.logical()

## This is parallel to Cmp.Mat.atomic() above --->  __keep parallel__ !
Logic.Mat.atomic <- function(e1, e2) { ## result will typically be "like" e1:
    l2 <- length(e2 <- as.logical(e2))
    n1 <- prod(d <- e1@Dim)
    if(n1 && n1 < l2)
        stop(sprintf("dim [product %d] do not match the length of object [%d]",
                     n1, l2))
    if(.Generic == "&" && l2 && allTrue (e2)) return(.M2kind(e1, "l"))
    if(.Generic == "|" && l2 && allFalse(e2)) return(.M2kind(e1, "l"))
    cl <- class(e1)
    if(l2 == 0)
        return(if(n1 == 0) .M2kind(e1, "l") else as.logical(e2))
    ## else
    cl1 <- getClassDef(cl)
    slots1 <- names(cl1@slots)
    has.x <- any("x" == slots1)# *fast* check for "x" slot presence
    if(l2 > 1 && has.x)
        return(if(n1 == 0) {
                   ## use `slots1`   sNms <-  .slotNames(e1)
                   r <- copyClass(e1, class2(cl, "l"),
                                  sNames = slots1[slots1 != "x"], check = FALSE)
                   r@x <- callGeneric(e1@x, e2)
                   r
               } else # cannot simply compare e2, e1@x -> use another method
                   callGeneric(e1, Matrix(e2, nrow=d[1], ncol=d[2]))
               )
    ## else
    isDiag <- extends(cl1, "diagonalMatrix")
    r0 <- callGeneric(0, e2)
    ar0 <- all(r0)
    if(isDiag && ar0) # all off-diagonals became TRUE ==> dense symmetric
        return(callGeneric(.diag2dense(e1, shape = "s"), e2))
    r <- callGeneric(if(has.x) e1@x else TRUE, e2)
    isTri <- extends(cl1, "triangularMatrix")
    ## Udg: append the diagonal at *end*, similar to diagU2N():
    Udg <- (isTri || isDiag) && e1@diag == "U"
    r. <- if(Udg) c(r, d.r <- callGeneric(..diag.x(e1), e2)) else r
    ## trivial case first (beware of NA)
    if(isTRUE(ar0 && all(r.))) { # |--> all TRUE
        r <- new(if(d[1] == d[2]) "lsyMatrix" else "lgeMatrix")
        r@Dim <- d
        r@Dimnames <- e1@Dimnames
        r@x <- rep.int(TRUE, n1)
    } else if(extends(cl1, "denseMatrix")) {
        full <- !.isPacked(e1) # << both "dtr" and "dsy" are 'full'
        if(full || allFalse(r0) || extends(cl1, "symmetricMatrix")) {
            if(sparseDefault(r, tol = 0.25)) # switch to sparse result
                return( callGeneric(as(e1, "sparseMatrix"), e2) )
            ## else
            ## FIXME? using copyClass() to copy "relevant" slots
            r <- new(class2(cl, "l"), x = r, Dim = d, Dimnames = e1@Dimnames)
            if(extends(cl1, "symmetricMatrix")) {
                r@uplo <- e1@uplo
            } else if(isTri) {
                r@uplo <- e1@uplo
                r@diag <- e1@diag
            }
        } else {
            ## packed matrix with structural 0 and r0 is not all FALSE:
            ##--> result cannot be packed anymore
            ## [dense & packed & not symmetric ] ==> must be "ltp*" :
            if(!extends(cl1, "ltpMatrix"))
                stop("internal bug in \"Logic\" method (Logic.Mat.atomic); please report")
            rx <- rep_len(r0, n1)
            rx[indTri(d[1], upper = (e1@uplo == "U"), diag=TRUE)] <- r.
            r <- new("lgeMatrix", x = rx, Dim = d, Dimnames = e1@Dimnames)
        }
    } else if(isDiag) { ## && !ar0 := !all(r0)
        if(Udg)
            e1 <- .diagU2N(e1, cl = cl1, checkDense = TRUE)
        tx <- typeof(r.)
        if(tx != typeof(e1@x))
            e1 <- .M2kind(e1, .kind.type[[tx]])
        e1@x <- if(Udg) d.r else r
        return(if(Udg) .diagN2U(e1) else e1)
    } else {
        ##---- e1 is(. , sparseMatrix) -----------------
        remainSparse <- allFalse(r0) ## <==> things remain sparse
        if(Udg) {          # e1 *is* unit-diagonal (triangular sparse)
            r1 <- callGeneric(1, e2)
            Udg <- all(r1)       # maybe Unit-diagonal (sparse) result
            ## if(!remainSparse) we'll use non0ind() which *has* unit-diag. indices at end
            ##
            if(Udg && remainSparse) {
            } else { ## result will not be unit-diagonal sparse
                e1 <- .diagU2N(e1, cl = cl1) # otherwise, result is U-diag
                ## FIXME? rather
                ## if(extends1of(cl1, c("CsparseMatrix", "RsparseMatrix","TsparseMatrix")) {
                if(extends(cl1, "CsparseMatrix")) {
                    ## repeat computation if e1 has changed
                    r. <- callGeneric(if(has.x) e1@x else TRUE, e2)
                }
            }
        }

        if(remainSparse) {
            if(!anyNA(r) && ((Ar <- all(r)) || !any(r))) {
                lClass <- class2(cl, "l") # is "lsparse*"
                r <- new(lClass)
                r@Dim <- d
                r@Dimnames <- e1@Dimnames
                if(Ar) {       # 'TRUE' instead of 'x': same sparsity:
                    for(n in intersect(c("i","j","p","uplo","diag"), slots1))
                        slot(r, n) <- slot(e1, n)
                    n <- if(has.x)
                             length(e1@x)
                         else if(any("p" == slots1))
                             e1@p[d[2]+1L]
                         else length(e1@i)
                    r@x <- rep.int(TRUE, n)
                } else {
                    ## !any(r): all FALSE: keep empty 'r' matrix
                    ## but may need a valid 'pointer' slot:
                    if(extends(lClass, "CsparseMatrix"))
                        r@p <- rep.int(0L, 1+ncol(r))
                    else if(extends(lClass, "RsparseMatrix"))
                        r@p <- rep.int(0L, 1+nrow(r))
                }
            } else {  # some TRUE, FALSE, NA : go via unique 'Tsparse'
                M <- asUniqueT(e1)
                nCl <- class2(class(M), 'l') # logical Tsparse
                sN <- slotNames(nCl)
                ## copy "the other slots" (important for "tr"/"sym"):
                r <- copyClass(M, nCl,
                               sNames = sN[is.na(match(sN, c("x","factors")))])
                r@x <- callGeneric(if(has.x) M@x else TRUE, e2)
                if(extends(cl1, "CsparseMatrix"))
                    r <- .M2C(r)
                else if(extends(cl1, "RsparseMatrix"))
                    r <- .M2R(r)
            }
        } else {
            ## non sparse result; triangularity also gone, typically
            lClass <- if(extends(cl1, "symmetricMatrix"))
                "lsyMatrix" else "lgeMatrix"
            Matrix.message(sprintf("sparse to dense (%s) coercion in '%s' -> %s",
                               lClass, .Generic, "Logic.Mat.atomic"),
                       .M.level = 2)
            rx <- rep_len(r0, n1)

            ## Here, we assume that 'r.' and the indices align (!)
            encI <- .Call(m_encodeInd,
                          non0ind(e1, cl1, uniqT=FALSE, xtendSymm=FALSE),
                          di = d, orig1=FALSE, checkBounds=FALSE)
            rx[1L + encI] <- r.
            r <- new(lClass, x = rx, Dim = d, Dimnames = e1@Dimnames)
        }
    }
    r
} # Logic.Mat.atomic

for(Mcl in c("lMatrix","nMatrix","dMatrix", "iMatrix", "zMatrix"))
    for(cl in c("logical", "numeric", "complex", "sparseVector"))
        setMethod("Logic", c(e1 = Mcl, e2 = cl), Logic.Mat.atomic)
rm(Logic.Mat.atomic, Mcl, cl)

### -- II -- sparse ----------------------------------------------------------

Ops.x.x.via.d <- function(e1, e2) callGeneric(.M2kind(e1, "d"), .M2kind(e2, "d"))


## Have lgC o lgC  and then lgT o lgT  Logic - quite similarly -
## also lsC o *  and ltC o * :

## Here's the common functionality
.do.Logic.lsparse <- function(e1,e2, d, dn, isOR, ij1, ij2) {

    ## NB non-diagonalMatrix := Union{ general, symmetric, triangular}
    gen1 <- extends(cD1 <- getClassDef(class(e1)), "generalMatrix")
    gen2 <- extends(cD2 <- getClassDef(class(e2)), "generalMatrix")
    sym1 <- !gen1 && extends(cD1, "symmetricMatrix")
    sym2 <- !gen2 && extends(cD2, "symmetricMatrix")
    tri1 <- !gen1 && !sym1
    tri2 <- !gen2 && !sym2
    G <- gen1 && gen2
    S <- sym1 && sym2 && e1@uplo == e2@uplo
    T <- tri1 && tri2 && e1@uplo == e2@uplo
    if(T && e1@diag != e2@diag) {
        ## one is "U" the other "N"
        if(e1@diag == "U")
            e1 <- diagU2N(e1)
        else ## (e2@diag == "U"
            e2 <- diagU2N(e2)
        shape <- "t"
    } else if(!G && !S && !T) {
        ## e.g. one symmetric, one general
        ## coerce to generalMatrix and go :
        if(!gen1) e1 <- .M2gen(e1)
        if(!gen2) e2 <- .M2gen(e2)
        shape <- "g"
    } else {
        shape <- if(T) "t" else if(S) "s" else "g"
    }

    ii <- WhichintersectInd(ij1, ij2, di=d)
    I1 <- ii[[1]] ; has1 <- length(I1) > 0
    I2 <- ii[[2]] ; has2 <- length(I2) > 0

    ## 1) common indices
    i <- ij1[I1, 1]
    j <- ij1[I1, 2]

    if(isOR) { ## i.e. .Generic == "|" i.e. not "&"
        x <- e1@x[I1] | e2@x[I2]

        ## 2) "e1 o  FALSE":
        x2 <- if(has1) e1@x[- I1] else e1@x # == callGeneric(e1@x[- I1], FALSE)
        ## 3) "0  o e1":
        x3 <- if(has2) e2@x[- I2] else e2@x # == callGeneric(FALSE, e2@x[- I2])
        i <- c(i,
               if(has1) ij1[-I1, 1] else ij1[, 1],
               if(has2) ij2[-I2, 1] else ij2[, 1])
        j <- c(j,
               if(has1) ij1[-I1, 2] else ij1[, 2],
               if(has2) ij2[-I2, 2] else ij2[, 2])
        x <- c(x, x2, x3)
    } else { ## AND
        x <- e1@x[I1] & e2@x[I2]
    }


    if(any(!(x. <- x | is.na(x)))) { ## drop 'FALSE's
        i <- i[x.]
        j <- j[x.]
        x <- x[x.]
    }
    if(shape == "g")
        new("lgTMatrix", Dim = d, Dimnames = dn, i = i, j = j, x = x)
    else
        new(paste0("l",shape,"TMatrix"), Dim = d, Dimnames = dn,
            i = i, j = j, x = x, uplo = e1@uplo)
}

Logic.lCMat <- function(e1, e2, isOR) {
    stopifnot(is.logical(isOR))
    d <- .Ops.checkDim(dim(e1), dim(e2))
    dn <- .Ops.checkDimNames(dimnames(e1), dimnames(e2))
    ## Very easy case first :
    if(identical(e1@i, e2@i) && identical(e1@p, e2@p)) {
        e1@x <- if(isOR) e1@x | e2@x else e1@x & e2@x
        if(.hasSlot(e1, "factors") && length(e1@factors))
            e1@factors <- list()
        return(e1)
    }
    ## else :

    .M2C(.do.Logic.lsparse(e1, e2, d = d, dn = dn, isOR = isOR,
                           ij1 = .Call(compressed_non_0_ij, e1, TRUE),
                           ij2 = .Call(compressed_non_0_ij, e2, TRUE)))
}

m.Logic.lCMat <- function(e1, e2) Logic.lCMat(e1, e2, isOR = .Generic == "|")

Logic.lTMat <- function(e1,e2) {
    d <- .Ops.checkDim(dim(e1), dim(e2))
    dn <- .Ops.checkDimNames(dimnames(e1), dimnames(e2))
    ## Very easy case first :
    if(identical(e1@i, e2@i) && identical(e1@j, e2@j)) {
        e1@x <- callGeneric(e1@x, e2@x)
        if(.hasSlot(e1, "factors") && length(e1@factors))
            e1@factors <- list()
        return(e1)
    }
    ## else :
    cld <- getClassDef(class(e1))
    .do.Logic.lsparse(e1, e2, d = d, dn = dn,
                      isOR = .Generic == "|",
                      ij1 = non0ind(e1, cld),
                      ij2 = non0ind(e2, cld))
}

setMethod("Logic", c(e1="lgCMatrix", e2="lgCMatrix"), m.Logic.lCMat)

setMethod("Logic", c(e1="lgTMatrix", e2="lgTMatrix"), Logic.lTMat)

rm(m.Logic.lCMat, Logic.lTMat)

setMethod("Logic", c(e1 = "lsCMatrix", e2 = "lsCMatrix"),
          function(e1, e2) {
              if(e1@uplo == e2@uplo)
                  Logic.lCMat(e1, e2, isOR = .Generic == "|")
              else Logic.lCMat(e1, t(e2), isOR = .Generic == "|")
          })

setMethod("Logic", c(e1 = "ltCMatrix", e2 = "ltCMatrix"),
          function(e1, e2) {
              isOR <- .Generic == "|"
              if(e1@uplo == e2@uplo) {
                  if(e1@diag == e2@diag) ## both "N" or both "U" (!)
                      Logic.lCMat(e1, e2, isOR=isOR)
                  else if(e1@diag == "U")
                      Logic.lCMat(diagU2N(e1), e2, isOR=isOR)
                  else ## e1@diag == "N"  *and*	 e2@diag == "U"
                      Logic.lCMat(e1, diagU2N(e2), isOR=isOR)
              } else {
                  ## differing triangle (upper *and* lower):
                  if(isOR) # both triangles => "general"
                      Logic.lCMat(.M2gen(e1), .M2gen(e2), isOR=TRUE)
                  else { ## have '&': all FALSE apart from diagonal
                      d <- .Ops.checkDim(dim(e1), dim(e2))
                      .diag2sparse(new("ldiMatrix", Dim = d,
                                       x = get(.Generic)(diag(e1), diag(e2))),
                                   kind = ".", shape = "t", repr = "C",
                                   uplo = e1@uplo)
                  }
              }
          })

## Now the other "Ops" for the "lgT" and "lgC" cases:
setMethod("Arith", c(e1="lgCMatrix", e2="lgCMatrix"), Ops.x.x.via.d)
setMethod("Arith", c(e1="lgTMatrix", e2="lgTMatrix"), Ops.x.x.via.d)
rm(Ops.x.x.via.d)

## More generally:  Arith: l* and n*  via  d* _or_ i* (future: "i")

A.M.lM <- function(e1, e2) {
    knd <- .Arith.kind(e1, e2, .Generic)
    callGeneric(.M2kind(e1, knd), .M2kind(e2, knd))
}
setMethod("Arith", c(e1="lsparseMatrix", e2="Matrix"), A.M.lM)
setMethod("Arith", c(e1="Matrix", e2="lsparseMatrix"), A.M.lM)
setMethod("Arith", c(e1="nsparseMatrix", e2="Matrix"), A.M.lM)
setMethod("Arith", c(e1="Matrix", e2="nsparseMatrix"), A.M.lM)
rm(A.M.lM)
##
for(cl in c("numeric", "logical")) # "complex", "raw" : basically "replValue"
  for(Mcl in c("lMatrix", "nMatrix")) {
      setMethod("Arith", c(e1=Mcl, e2=cl),
                function(e1, e2) callGeneric(.M2kind(e1, "d"), e2))
      setMethod("Arith", c(e1=cl, e2=Mcl),
                function(e1, e2) callGeneric(e1, .M2kind(e2,"d")))
  }
rm(cl, Mcl)

## FIXME: These are really too cheap: currently almost all go via dgC*() :
## setMethod("Compare", c(e1="lgCMatrix", e2="lgCMatrix"),
## setMethod("Compare", c(e1="lgTMatrix", e2="lgTMatrix"),
## setMethod("Compare", c(e1="lsparseMatrix", e2="lsparseMatrix"),
##           function(e1, e2) callGeneric(as(e1, "dgCMatrix"),
##                                        as(e2, "dgCMatrix")))
##. Have "Ops" below which only goes *conditionally* via Csparse
##.setMethod("Compare", c(e1="lsparseMatrix", e2="lsparseMatrix"),
##.          function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
##.                                       as(e2, "CsparseMatrix")))
## setMethod("Compare", c(e1="lgTMatrix", e2="lgTMatrix"),
##           function(e1, e2) callGeneric(as(e1, "dgCMatrix"),
##                                        as(e2, "dgCMatrix")))

###--- Sparse ... ----------


setMethod("Ops", c(e1="lsparseMatrix", e2="lsparseMatrix"),
          function(e1,e2) callGeneric(as(e1, "CsparseMatrix"),
                                      as(e2, "CsparseMatrix")))

setMethod("Logic", c(e1="lsparseMatrix", e2="ldenseMatrix"),
          function(e1,e2) callGeneric(as(e1, "generalMatrix"),
                                      as(e2, "sparseMatrix")))

setMethod("Logic", c(e1="ldenseMatrix", e2="lsparseMatrix"),
          function(e1,e2) callGeneric(as(e1, "sparseMatrix"),
                                      as(e2, "generalMatrix")))

setMethod("Logic", c(e1="lsparseMatrix", e2="lsparseMatrix"),
          function(e1,e2) {
              if(!is(e1,"generalMatrix"))
                  callGeneric(as(.M2gen(e1), "CsparseMatrix"), e2)
              else if(!is(e2,"generalMatrix"))
                  callGeneric(e1, as(.M2gen(e2), "CsparseMatrix"))
              else # both are general, i.e. lg[CRT]
                  callGeneric(as(e1, "CsparseMatrix"), as(e2, "CsparseMatrix"))
          })

## Via coercion to "generalMatrix":
## 1. symmetric  V
Mcl <- "zsCMatrix" # special because of 'trans':
setMethod("Arith", c(e1 = Mcl, e2 = Mcl),
          function(e1, e2) {
              Matrix.message(sprintf(
                  "suboptimal 'Arith' implementation of   %s  o  %s", class(e1), class(e2)))
              forceSymmetric(callGeneric(.M2gen(e1), .M2gen(e2)),
                             trans = if((t1 <- e1@trans) == e2@trans) t1 else "T")
          })
for(Mcl in c(paste0(c("d","l","n", "i"), "sCMatrix")))
    setMethod("Arith", c(e1 = Mcl, e2 = Mcl),
              function(e1, e2) {
                  Matrix.message(sprintf(
                      "suboptimal 'Arith' implementation of   %s  o  %s", class(e1), class(e2)))
                  forceSymmetric(callGeneric(.M2gen(e1), .M2gen(e2)))
              })

## 2.                              triangular  V
for(Mcl in c(paste0(c("d","l","n", "i", "z"), "tCMatrix")))
    setMethod("Arith", c(e1 = Mcl, e2 = Mcl),
              function(e1, e2) {
                  Matrix.message(sprintf(
                      "suboptimal 'Arith' implementation of   %s  o  %s", class(e1), class(e2)))
                  r <- callGeneric(.M2gen(e1), .M2gen(e2))
                  if(e1@uplo == e2@uplo) forceTriangular(r, uplo = e1@uplo) else r
          })
rm(Mcl)

##-------- originally from ./dgCMatrix.R --------------------
## for [dliz]gCMatrix (if triangular=FALSE),  or [dliz]tCMatrix (if triangular=TRUE)
.Arith.Csparse <- function(e1, e2, Generic, class.,
                           triangular = FALSE, check.dimnames = TRUE) {
    ## Generic is one of  "+", "-", "*", "^", "%%", "%/%", "/"

    ## triangular:  TRUE  iff e1,e2 are triangular  _and_  e1@uplo == e2@uplo
    d <- .Ops.checkDim(dim(e1), dim(e2))
    dn <- .Ops.checkDimNames(dimnames(e1), dimnames(e2),
                             check = check.dimnames)
    knd <- .Arith.kind(e1, e2, Generic)
    if(triangular) {
        ## need these for the 'x' slots in any case
        e1 <- .Call(R_sparse_diag_U2N, e1)
        e2 <- .Call(R_sparse_diag_U2N, e2)
        ## slightly more efficient than non0.i() or non0ind():
        ij1 <- .Call(compressed_non_0_ij, e1, isC=TRUE)
        ij2 <- .Call(compressed_non_0_ij, e2, isC=TRUE)
        newTMat <- function(i,j,x)
            new(paste0(knd,"tTMatrix"), Dim = d, Dimnames = dn, i = i, j = j,
                x = .v2kind(x, knd),
                uplo = e1@uplo)
    } else {
        cld <- getClassDef(class.)
        ij1 <- non0ind(e1, cld)
        ij2 <- non0ind(e2, cld)
        newTMat <- function(i,j,x)
            new(paste0(knd,"gTMatrix"), Dim = d, Dimnames = dn, i = i, j = j,
                x = .v2kind(x, knd))
    }

    ## NB.  knd != "n"
    if(!.hasSlot(e1, "x")) e1 <- .M2kind(e1, knd)
    if(!.hasSlot(e2, "x")) e2 <- .M2kind(e2, knd)
    switch(Generic,
       "+" = , "-" =
       {
           ## care for over-allocated 'x' slot:
           nc1 <- d[2] + 1L
           if((nz <- e1@p[nc1]) < length(e1@x)) e1@x <- e1@x[seq_len(nz)]
           if((nz <- e2@p[nc1]) < length(e2@x)) e2@x <- e2@x[seq_len(nz)]
           ## special "T" convention: repeated entries are *summed*
           .M2C(newTMat(i = c(ij1[,1], ij2[,1]),
                        j = c(ij1[,2], ij2[,2]),
                        x = if(Generic == "+")
                                 c(e1@x,   e2@x)
                            else c(e1@x, - e2@x)))
       },

       "*" =
       {
           ##  X * 0 == 0 * X == 0 --> keep common non-0 {treating NA * 0 ==> 0 !}
           ii <- WhichintersectInd(ij1, ij2, di=d)
           ij <- ij1[ii[[1L]], , drop = FALSE] # == ij2[[ii[[2]], ,drop=FALSE]
           .M2C(newTMat(i = ij[,1L],
                        j = ij[,2L],
                        x = e1@x[ii[[1L]]] * e2@x[ii[[2L]]]))
       },

       "^" =
       {
           ## Cases of  0 - ops :
           ## 1)  X^0 := 1  (even for X=0, and NA) ==> dense
           ## 2a) 0^NA := NA (!!)                      =====
           ## 2b) 0^<complex> := NaN + 1i*NaN (!!)
           ## 2c) 0^ Y := 0    for Y > 0
           ## 2d) 0^ Y := Inf  for Y < 0
           ## 3)  x^ y (both not 0):  dense
           ## ... {too complicated; rarely sparse, but e.g., M ^ diag(x = 1:n)
           r <- get(Generic)(.M2unpacked(e1), .M2m(e2))
           if(sparseDefault(r)) as(r, "sparseMatrix") else r
       },

       "%%" = , "%/%" = , "/" = ## 0 op 0	 |-> NaN => dense
       {
           ## FIXME? if(sparseDefault(r)) .... as above ?
           get(Generic)(.M2unpacked(e1), e2)
       }) # end{switch(..)}
} # .Arith.Csparse()


.Arith.CsparseCl <- function(e1,e2) .Arith.Csparse(e1,e2, .Generic, class.= class(e1))
for(cl in paste0(c("z","d","i","l"), "gCMatrix")) { # *not* "n" - have no @x
    setMethod("Arith", c(e1 = cl, e2 = cl), .Arith.CsparseCl)
}
rm(cl, .Arith.CsparseCl)


.Arith.CsparseTri <- function(e1, e2) {
    U1 <- e1@uplo
    ## will the result definitely be triangular?
    isTri <- U1 == e2@uplo && .Generic != "^"
    knd <- .Arith.kind(e1, e2, .Generic) # knd is *not* "n" -> have '@x'
    if(isTri)
        .Arith.Csparse(e1,e2, .Generic, class. = paste0(knd, "tCMatrix"),
                       triangular = TRUE)
    else # lowerTri  o  upperTri: |--> "all 0" {often} -- FIXME?
        .Arith.Csparse(.M2gen(e1), .M2gen(e2),
                       .Generic, class. = paste0(knd, "gCMatrix"))
}

for(cl in paste0(c("z","d","i","l","n"), "tCMatrix")) {
    setMethod("Arith", c(e1 = cl, e2 = cl), .Arith.CsparseTri)
}
rm(cl, .Arith.CsparseTri)



## TODO : Consider going a level up, and do this for all "Ops"
##
## NB: For "dgCMatrix" have special method ==> this is for dsC*, lgC*, ...
##     now also for Tsparse etc {*must* as method directly: "callGeneric()"}
.Arith.CM.atom <- function(e1, e2) {
    knd <- .Arith.kind(e1, e2, .Generic)
    if(length(e2) == 1L) { ## e.g.,  Mat ^ a
        if(is0(callGeneric(0, e2))) { ## remain sparse, symm., tri.,...
            e1 <- .M2kind(e1, knd)
            if(!extends(cld <- getClassDef(class(e1)), "CsparseMatrix"))
                cld <- getClassDef(class(e1 <- as(e1, "CsparseMatrix")))
            if(extends(cld, "triangularMatrix") &&
               e1@diag == "U" && !all(1 == callGeneric(1, e2)))
                e1 <- .diagU2N(e1, cld)
            e1@x <- callGeneric(e1@x, e2)
            if(.hasSlot(e1, "factors") && length(e1@factors))
                e1@factors <- list()
            return(.promote(e1))
        }
    }
    ## if(.M.shape(e1) == "g" && .M.kind(e1) == knd && .M.repr(e1) == "C") { # <k>gCMatrix already
    ##     e1@x <- callGeneric(e1@x, e2)
    ##     e1
    ## } else ## all other (potentially non-sparse) cases: give up symm, tri,..
        callGeneric(.M2C(.M2gen(.M2kind(e1, knd))), e2)
}
## The same,  e1 <-> e2 :
.Arith.atom.CM <- function(e1, e2) {
    knd <- .Arith.kind(e1, e2, .Generic)
    if(length(e1) == 1L) {
        if(is0(callGeneric(e1, 0))) {
            e2 <- .M2kind(e2, knd)
            if(!extends(cld <- getClassDef(class(e2)), "CsparseMatrix"))
                cld <- getClassDef(class(e2 <- as(e2, "CsparseMatrix")))
            if(extends(cld, "triangularMatrix") &&
               e2@diag == "U" && !all(1 == callGeneric(e1, 1)))
                e2 <- .diagU2N(e2, cld)
            e2@x <- callGeneric(e1, e2@x)
            if(.hasSlot(e2, "factors") && length(e2@factors))
                e2@factors <- list()
            return(.promote(e2))
        }
    }
    ## if(.M.shape(e2) == "g" && .M.kind(e2) == knd && .M.repr(e2) == "C") { # <k>gCMatrix already
    ##     e2@x <- callGeneric(e1, e2@x)
    ##     e2
    ## } else ## all other (potentially non-sparse) cases: give up symm, tri,..
        callGeneric(e1, .M2C(.M2gen(.M2kind(e2, knd))))
}

for(cl in c("logical", "numeric", "complex")) { # *not* "sparseVector"
  for(Mcl in c("CsparseMatrix",
               paste0(c("l","n","d", "i", "z"), "sparseMatrix"))) {
      setMethod("Arith", c(e1 = Mcl, e2 = cl), .Arith.CM.atom)
      setMethod("Arith", c(e1 = cl, e2 = Mcl), .Arith.atom.CM)
  }
}
##
rm(cl, Mcl, .Arith.atom.CM, .Arith.CM.atom)


##' compute indices for recycling <numeric> of length 'len'
##' to match sparseMatrix 'spM' (Csparse or Rsparse)
.Ops.recycle.ind <- function(spM, len, Csparse = TRUE) {
    n <- prod(d <- dim(spM))
    if(n && n < len) stop("vector too long in Matrix - vector operation")
    if(n %% len != 0) ## identical warning as in main/arithmetic.c
        warning("longer object length\n\tis not a multiple of shorter object length")
    ## TODO(speedup!): construction of value (1L + in0 %% len) via one .Call()
    in0 <- .Call(m_encodeInd, .Call(compressed_non_0_ij, spM, Csparse),
                 d, FALSE, FALSE)
    1L + in0 %% len
}

A.M.n <- function(e1, e2) { #  <sparseM> o <atomic>
    if(!(l2 <- length(e2))) # return 0-vector of correct kind, as matrix()+<0>
        return(if(length(e1)) vector(.type.kind[[.Arith.kind(e1, e2, .Generic)]]) else e1)
    if((rep1 <- .M.repr(e1)) == "T") { e1 <- .M2C(e1) ; rep1 <- "C" } # to avoid early inf.loop
    is0f <- is0(f0 <- callGeneric(0, e2)) #
    if((a0 <- all(is0f)) || mean(is0f) > 7/8) { ## result keeps sparseness structure of e1
                                        #--- 7/8 is *somewhat* arbitrary
        Udg <- (.M.shape(e1) == "t") && e1@diag == "U" # unit-diagonal
        if(Udg) e1 <- ..diagU2N(e1)
        if(l2 > 1L) {
            if(a0) # "recycle" e2 "carefully"
                e2 <- e2[.Ops.recycle.ind(e1, len = l2, Csparse = rep1 == "C")]
            else # not all 0, but many --> e2 is typically sparse
                return(callGeneric(if(Udg) ..diagN2U(e1, sparse=TRUE) else e1,
                                   as(e2, "sparseVector")))
        }
        if((Mknd <- .M.kind(e1)) == "n") e1 <- .M2kind(e1, (Mknd <- "l")) # to get 'x' slot
        knd <- .kind.type[[typeof(r <- callGeneric(e1@x, e2))]]
        i2d <- .Arith.i.d(e1, knd, Mknd)
        e1 <- .Arith.fixK(e1, knd, Mknd, i2d)
        e1@x <- if(i2d) as.double(r) else r
        if(.hasSlot(e1, "factors") && length(e1@factors))
            e1@factors <- list()
        if(Udg && all(diag(e1) == 1)) ..diagN2U(e1, sparse=TRUE) else e1
    } else { ## non-sparse result, as '0 o e2' is not (all) 0
        r <- as(e1, "matrix")
        knd <- .Arith.kind(e1, e2, .Generic)
        class <- paste0(knd, "ge")
        if(l2 == 1 && .M.kind(e1) != "n") {
            r[] <- f0
            r[non0ind(e1, getClassDef(paste0(knd, "gCMatrix"))) + 1L] <-
                callGeneric(e1@x, e2)
            .m2dense(r, class)
        } else {
            .m2dense(callGeneric(r, e2), class)
        }
    }
}

A.n.M <- function(e1, e2) { #  <atomic> o <sparseM>
    if(!(l1 <- length(e1))) ## return 0-vector of correct kind, as <0> + matrix()
        return(if(length(e2)) vector(.type.kind[[.Arith.kind(e1, e2, .Generic)]]) else e2)
    if((rep2 <- .M.repr(e2)) == "T") { e2 <- .M2C(e2) ; rep2 <- "C" } # to avoid early inf.loop
    is0f <- is0(f0 <- callGeneric(e1, 0))
    if((a0 <- all(is0f)) || mean(is0f) > 7/8) { ## result keeps sparseness structure of e1
                                        #--- 7/8 is *somewhat* arbitrary
        Udg <- (.M.shape(e2) == "t") && e2@diag == "U" # unit-diagonal
        if(Udg) e2 <- ..diagU2N(e2)
        if(l1 > 1L) {
            if(a0) # "recycle" e1 "carefully"
                e1 <- e1[.Ops.recycle.ind(e2, len = l1, Csparse = rep2 == "C")]
            else # not all 0, but many --> e1 is typically sparse
                return(callGeneric(as(e1, "sparseVector"),
                                   if(Udg) ..diagN2U(e2, sparse=TRUE) else e2))
        }
        if((Mknd <- .M.kind(e2)) == "n") e2 <- .M2kind(e2, (Mknd <- "l")) # to get 'x' slot
        knd <- .kind.type[[typeof(r <- callGeneric(e1, e2@x))]]
        i2d <- .Arith.i.d(e2, knd, Mknd)
        e2 <- .Arith.fixK(e2, knd, Mknd, i2d)
        e2@x <- if(i2d) as.double(r) else r
        if(.hasSlot(e2, "factors") && length(e2@factors))
            e2@factors <- list()
        if(Udg && all(diag(e2) == 1)) ..diagN2U(e2, sparse=TRUE) else e2
    } else { ## non-sparse result, as 'e1 o 0' is not (all) 0
        r <- as(e2, "matrix")
        knd <- .Arith.kind(e1, e2, .Generic)
        class <- paste0(knd, "ge")
        if(l1 == 1 && .M.kind(e2) != "n") {
            r[] <- f0
            r[non0ind(e2, getClassDef(paste0(knd, "gCMatrix"))) + 1L] <-
                callGeneric(e1, e2@x)
            .m2dense(r, class)
        } else {
            .m2dense(callGeneric(e1, r), class)
        }
    }
}
## coercing to "general*" / "dgC*"  would e.g. lose symmetry of  '3 * S'
for(kind in c("l","n","d", "i", "z")) {
    for(Mcl in paste0(kind, c("gC", "sparse" # "sparse" ?? -- have Csparsematrix above
                              ),
                      "Matrix"))
    for(cl in c("logical", "numeric", "complex")) {
        setMethod("Arith", c(e1 = Mcl, e2 = cl), A.M.n)
        setMethod("Arith", c(e1 = cl, e2 = Mcl), A.n.M)
    }
}
##
rm(kind, Mcl, cl, A.M.n, A.n.M)


##-------- originally from ./dgTMatrix.R --------------------

## Uses the triplet convention of *adding* entries with same (i,j):
setMethod("+", c(e1 = "dgTMatrix", e2 = "dgTMatrix"),
          function(e1, e2) {
              .Ops.checkDim(dim(e1), dim(e2))
              new("dgTMatrix", i = c(e1@i, e2@i), j = c(e1@j, e2@j),
                  x = c(e1@x, e2@x), Dim = e1@Dim, Dimnames = e1@Dimnames)
          })


##-------- originally from ./Csparse.R --------------------

setMethod("Arith", c(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
          function(e1, e2) {
              ## go via "symmetric" if both are symmetric, etc...
              s1 <- .M.shape(e1)
              s2 <- .M.shape(e2)
              ## as(*, "d.CMatrix") is deprecated:
              ## viaCl <- paste0("d", if(s1 == s2) s1 else "g", "CMatrix")
              knd <- .Arith.kind(e1,e2, .Generic)
              if(s1 != s2) ## go via "general"
                  callGeneric(.M2gen(.M2kind(e1, knd)),
                              .M2gen(.M2kind(e2, knd)))
              else
                  callGeneric(.M2kind(e1, knd), .M2kind(e2, knd))
          })

setMethod("Logic", c(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
          ## go via "symmetric" if both are symmetric, etc...
          function(e1, e2) {
              s1 <- .M.shape(e1)
              s2 <- .M.shape(e2)
              ## as(*, "d.CMatrix") is deprecated:
              ## viaCl <- paste0("l", if(s1 == s2) s1 else "g", "CMatrix")
              if(s1 != s2) ## go via "general"
                  callGeneric(.M2gen(.M2kind(e1, "l")),
                              .M2gen(.M2kind(e2, "l")))
              else
                  callGeneric(.M2kind(e1, "l"), .M2kind(e2, "l"))
          })

setMethod("Compare", c(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
          function(e1, e2) {
              ## How do the "0" or "FALSE" entries compare?
              ## Depends if we have an "EQuality RELation" or not:
              EQrel <- switch(.Generic,
                              "==" =, "<=" =, ">=" = TRUE,
                              "!=" =, "<"  =, ">"  = FALSE)
              if(EQrel) {
                  ## The (0 op 0) or  (FALSE op FALSE) comparison gives TRUE
                  ## -> result becomes *dense*; the following may be suboptimal
                  return(callGeneric(.sparse2dense(e1), .sparse2dense(e2)))
              }
              ## else: INequality:   0 op 0 gives FALSE ---> remain sparse!

              d <- .Ops.checkDim(dim(e1), dim(e2))
              cD1 <- getClassDef(class(e1))
              cD2 <- getClassDef(class(e2))
              Matrix.message(sprintf("Compare <Csparse> -- \"%s\" %s \"%s\" :\n",
                                 cD1@className, .Generic, cD2@className),
                         .M.level = 2)

              ## NB non-diagonalMatrix := Union{ general, symmetric, triangular}
              gen1 <- extends(cD1, "generalMatrix")
              gen2 <- extends(cD2, "generalMatrix")
              sym1 <- !gen1 && extends(cD1, "symmetricMatrix")
              sym2 <- !gen2 && extends(cD2, "symmetricMatrix")
              tri1 <- !gen1 && !sym1
              tri2 <- !gen2 && !sym2
              G <- gen1 && gen2
              S <- sym1 && sym2 && e1@uplo == e2@uplo
              T <- tri1 && tri2 && e1@uplo == e2@uplo

              if(T) { # can not use diag="U" -- "U" < != > "U" |--> FALSE on diagonal
                  if(e1@diag == "U")
                      e1 <- ..diagU2N(e1)
                  if(e2@diag == "U")
                      e2 <- ..diagU2N(e2)
                  shape <- "t"
              }
              else if(!G && !S && !T) {
                  ## e.g. one symmetric, one general
                  ## coerce to generalMatrix and go :
                  if(!gen1) e1 <- .M2gen(e1)
                  if(!gen2) e2 <- .M2gen(e2)
                  shape <- "g"
              } else { # not T but at least one of S, G :
                  shape <- if(S) "s" else "g"
              }

              dn <- .Ops.checkDimNames(dimnames(e1), dimnames(e2))
              ## ^^ FIXME: for 'S'; allow staying
              e1is.n <- extends(cD1, "nMatrix")
              e2is.n <- extends(cD2, "nMatrix")
              both.n <- e1is.n && e2is.n
              eq.pattern <- identical(e1@i, e2@i) && identical(e1@p, e2@p)
              ## the result object:
              ## FIXME: "n" result when e1 & e2 are "n", .. but for now only
              ## TODO? become "n" even whenever there are no NAs ??
              r.is.n <- both.n && eq.pattern
              newC <- sub("^.", if(r.is.n) "n" else "l", MatrixClass(class(e1)))
              ## For both.n && !eq.pattern, would need drop 0s from rx <- callGeneric(...)
              r <- if(S || T)
                        new(newC, uplo = e1@uplo) # diag = "N"  default
                   else new(newC)
              ## Easy case: identical sparsity pattern
              if(eq.pattern) {
                  if(e1is.n) {
                      if(e2is.n) # here have r.is.n
                          ## non-equality of identical pattern matrices: all FALSE
                          r@p <- rep.int(0L, d[2]+1L) # and r@i remains empty
                      else { # e1 pattern, e2@x
                          rx <- callGeneric(TRUE, e2@x)
                          if(allFalse(rx))
                              r@p <- rep.int(0L, d[2]+1L) # r@i, r@x remain empty
                          else {
                              r@x <- rx
                              r@i <- e2@i
                              r@p <- e2@p
                          }
                      }
                  } else if(e2is.n) { ## e1@x, e2 pattern
                      rx <- callGeneric(e1@x, TRUE)
                      if(allFalse(rx))
                          r@p <- rep.int(0L, d[2]+1L) # and r@i, r@x remain empty
                      else {
                          r@x <- rx
                          r@i <- e1@i
                          r@p <- e1@p
                      }
                  } else {		# both have 'x' slot
                      r@x <- callGeneric(e1@x, e2@x)
                      ## and all others are  '0 op 0' which give FALSE
                      r@i <- e1@i
                      r@p <- e1@p
                  }
                  r@Dim <- d
                  r@Dimnames <- dn
                  r
              } else {
                  ## now the 'x' slots ``match'' insofar as they are for the
                  ## same "space" (triangle for tri* and symm*; else rectangle)

                  ## not non0ind() which gives more;
                  ## want only those which correspond to 'x' slot
                  ij1 <- .Call(compressed_non_0_ij, e1, TRUE)
                  ij2 <- .Call(compressed_non_0_ij, e2, TRUE)
                  ii <- WhichintersectInd(ij1, ij2, di=d)
                  I1 <- ii[[1]]; has1 <- length(I1) > 0
                  I2 <- ii[[2]]; has2 <- length(I2) > 0

                  ## potentially could be faster for 'nsparse'
                  ## but this is simple:
                  e1x <- if(e1is.n) rep.int(1L, length(e1@i)) else e1@x
                  e2x <- if(e2is.n) rep.int(1L, length(e2@i)) else e2@x

                  ## 1) common
                  x <- callGeneric(e1x[I1],
                                   e2x[I2])
                  ## 2) "e1 o  0":
                  x2 <- callGeneric(if(has1) e1x[- I1] else e1x, 0)
                  ## 3) "0  o e2":
                  x3 <- callGeneric(0, if(has2) e2x[- I2] else e2x)

                  i <- c(ij1[I1, 1],
                         if(has1) ij1[-I1, 1] else ij1[, 1],
                         if(has2) ij2[-I2, 1] else ij2[, 1])
                  j <- c(ij1[I1, 2],
                         if(has1) ij1[-I1, 2] else ij1[, 2],
                         if(has2) ij2[-I2, 2] else ij2[, 2])
                  x <- c(x, x2, x3)
                  if(any(i0x <- is0(x))) { # drop 'FALSE's
                      n0 <- !i0x
                      i <- i[n0]
                      j <- j[n0]
                      x <- x[n0]
                  }
                  .M2C(if(e1is.n && e2is.n)
                           new(paste0("n",shape,"TMatrix"), Dim = d,
                               Dimnames = dn, i = i, j = j)
                       else if(!S && !T)
                           new(paste0("l",shape,"TMatrix"), Dim = d,
                               Dimnames = dn, i = i, j = j, x = x)
                       else # S or T
                           new(paste0("l",shape,"TMatrix"), Dim = d,
                               Dimnames = dn, i = i, j = j, x = x,
                               uplo = e1@uplo))
              }
          })


##-------- originally from ./sparseMatrix.R --------------------

## "Arith" short cuts / exceptions
setMethod("-", c(e1 = "sparseMatrix", e2 = "missing"),
          function(e1, e2) {
              e1 <- diagU2N(e1)
              e1@x <- -e1@x
              if(.hasSlot(e1, "factors") && length(e1@factors))
                  e1@factors <- list()
              e1
          })
## with the following exceptions (" unary `-` "):
##--- __FIXME__ use_iMatrix ==>  s/d/i/ :
setMethod("-", c(e1 = "nsparseMatrix", e2 = "missing"), function(e1, e2) - .M2kind(e1, "d"))
setMethod("-", c(e1 = "lsparseMatrix", e2 = "missing"), function(e1, e2) - .M2kind(e1, "d"))
setMethod("-", c(e1 = "indMatrix", e2 = "missing"), function(e1, e2) -as(e1, "dsparseMatrix"))

## Group method  "Arith"

## have CsparseMatrix methods above
## which may preserve "symmetric", "triangular" -- simply defer to those:

setMethod("Ops", c(e1 = "sparseMatrix", e2 = "nsparseMatrix"),
          function(e1, e2)
              callGeneric(as(e1, "CsparseMatrix"), .M2kind(e2, "l")))
setMethod("Ops", c(e1 = "nsparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2)
              callGeneric(.M2kind(e1, "l"), as(e2, "CsparseMatrix")))

## these were 'Arith', now generalized:
if(FALSE) { ## just shifts the ambiguity warnings ..
## <sparse> o <sparse> more complicated - against PITA disambiguation warnings:
setMethod("Ops", c(e1 = "TsparseMatrix", e2 = "TsparseMatrix"),
          function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
                                       as(e2, "CsparseMatrix")))
setMethod("Ops", c(e1 = "TsparseMatrix", e2 = "CsparseMatrix"),
          function(e1, e2) callGeneric(as(e1, "CsparseMatrix"), e2))
setMethod("Ops", c(e1 = "CsparseMatrix", e2 = "TsparseMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "CsparseMatrix")))
}
## catch the rest:  Rsparse*  and  T*  o  R*
setMethod("Ops", c(e1 = "sparseMatrix", e2 = "sparseMatrix"),
          function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
                                       as(e2, "CsparseMatrix")))

setMethod("Ops", c(e1 = "sparseMatrix", e2 = "numeric"),
          function(e1, e2) callGeneric(as(e1, "CsparseMatrix"), e2))
setMethod("Ops", c(e1 = "numeric", e2 = "sparseMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "CsparseMatrix")))

## setMethod("Compare", c(e1 = "sparseMatrix", e2 = "sparseMatrix"),
##           function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
##                                        as(e2, "CsparseMatrix")))


###-------- sparseVector -------------
###-------- ============ -------------

## Catch all remaining
setMethod("Ops", c(e1 = "sparseVector", e2 = "ANY"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))
setMethod("Ops", c(e1 = "ANY", e2 = "sparseVector"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))

## 1)  spVec  o  (sp)Vec : -------------

## FIXME:
##   2. <spVec>  o  <non-NA numeric>  should also happen directly and
##                               |-> sparse for o = {'*', "/", '&&', '==', ...

setMethod("Ops", c(e1 = "sparseVector", e2 = "vector"),
          function(e1, e2) {
              if(is.object(e2) || is.array(e2) || is.recursive(e2))
                  stop(gettextf("invalid class \"%s\" object in '%s' method",
                                data.class(e2), "Ops"),
                       domain = NA)
              if(length(e2) == 1) { ## scalar ------ special case - "fast"
                  if(all0(callGeneric(FALSE, e2))) { # result remains sparse
                      if(is(e1, "nsparseVector")) { # no 'x' slot, i.e. all TRUE
                          r <- callGeneric(TRUE, e2)
                          if(is.logical(r)) {
                              if(isTRUE(all(r))) # (could be NA)
                                  e1	# result unchanged
                              else
                                  newSpVec("lsparseVector", x = r, e1)
                          } else {
                              newSpVec(paste0(if(is.integer(r)) "i" else "d",
                                              "sparseVector"),
                                       x = r, e1)
                          }
                      } else { # has x slot
                          r <- callGeneric(e1@x, e2)
                          if(identical(class(r), class(e1@x))) {
                              e1@x <- r
                              e1
                          } else {
                              newSpVec(paste0(.M.kind(r), "sparseVector"),
                                       x = r, e1)
                          }
                      }
                  }
                  else ## non-sparse result
                      callGeneric(sp2vec(e1), e2)
              }
              else ## e2 is not scalar
                  callGeneric(e1, as(e2, "sparseVector"))
          })

setMethod("Ops", c(e1 = "vector", e2 = "sparseVector"),
          function(e1, e2) {
              if(is.object(e1) || is.array(e1) || is.recursive(e1))
                  stop(gettextf("invalid class \"%s\" object in '%s' method",
                                data.class(e1), "Ops"),
                       domain = NA)
              if(length(e1) == 1) { ## scalar ------ special case - "fast"
                  if(all0(callGeneric(e1, FALSE))) { # result remains sparse
                      if(is(e2, "nsparseVector")) { # no 'x' slot, i.e. all TRUE
                          r <- callGeneric(e1, TRUE)
                          if(is.logical(r)) {
                              if(isTRUE(all(r))) # (could be NA)
                                  e2	# result unchanged
                              else
                                  newSpVec("lsparseVector", x = r, e2)
                          } else {
                              newSpVec(paste0(if(is.integer(r)) "i" else "d",
                                              "sparseVector"),
                                       x = r, e2)
                          }
                      } else { # has x slot
                          r <- callGeneric(e1, e2@x)
                          if(identical(class(r), class(e2@x))) {
                              e2@x <- r
                              e2
                          } else {
                              newSpVec(paste0(.M.kind(r), "sparseVector"),
                                       x = r, e2)
                          }
                      }
                  }
                  else ## non-sparse result
                      callGeneric(e1, sp2vec(e2))
              }
              else ## e1 is not scalar
                  callGeneric(as(e1, "sparseVector"), e2)
          })


Ops.spV.spV <- function(e1, e2) {
    n1 <- e1@length
    n2 <- e2@length
    if(!n1 || !n2) ## return 0-length :
        return(if(is.na(match(.Generic, .ArithGenerics)))
                   logical()
               else numeric())
    ## else  n1, n2 >= 1 :
    if(n1 != n2) {
        if(n1 < n2) {
            n <- n1 ; N <- n2
        } else {
            n <- n2 ; N <- n1
        }
        if(n == 1L) { # simple case, do not really recycle
            if(n1 < n2)
                return(callGeneric(sp2vec(e1), e2))
            else return(callGeneric(e1, sp2vec(e2)))
        }
        ## else : 2 <= n < N
        if(N %% n != 0)
            warning("longer object length is not a multiple of shorter object length")
        ## recycle the shorter one
        if(n1 < n2) {
            e1 <- rep(e1, length.out = N)
        } else {
            e2 <- rep(e2, length.out = N)
        }
    } else { ## n1 == n2
        N <- n1
    }
    ## ---- e1 & e2 now are both of length 'N' ----

    ## First check the (0  o  0) result
    is1n <- extends(class(e1), "nsparseVector")
    is2n <- extends(class(e2), "nsparseVector")
    r00 <- callGeneric(if(is1n) FALSE else as0(e1@x),
                       if(is2n) FALSE else as0(e2@x))
    if(is0(r00)) { ## -> sparseVector
        e1x <- if(is1n) TRUE else e1@x
        e2x <- if(is2n) TRUE else e2@x
        sp <- .setparts(e1@i, e2@i)
        ## Idea: Modify 'e2' and return it :
        new.x <- c(callGeneric(e1x[sp[["ix.only"]]], 0), # e1-only
                   callGeneric(0, e2x[sp[["iy.only"]]]), # e2-only
                   callGeneric(e1x[sp[["my"]]],          # common to both
                               e2x[sp[["mx"]]]))
        i. <- c(sp[["x.only"]], sp[["y.only"]], sp[["int"]])
        cl2x <- typeof(e2x) ## atomic "class"es - can use in is(), as(), too:
        if(!is2n && is(new.x, cl2x)) {
            i. <- sort.int(i., method = "quick", index.return=TRUE)
            e2@x <- as(new.x, cl2x)[i.$ix]
            e2@i <- i.$x
            e2
        } else {
            newSpV(paste0(.kind.type[[typeof(new.x)]],"sparseVector"),
                   x = new.x, i = i., length = e2@length)
        }
    } else { ## 0 o 0  is NOT in {0 , FALSE} --> "dense" result
        callGeneric(sp2vec(e1),	 sp2vec(e2))
    }
} ## {Ops.spV.spV}

## try to use it in all cases
setMethod("Ops", c(e1 = "sparseVector", e2 = "sparseVector"),
          Ops.spV.spV)
## was    function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))

setMethod("Arith", c(e1 = "sparseVector", e2 = "sparseVector"),
          function(e1, e2) callGeneric(as(e1, "dsparseVector"),
                                       as(e2, "dsparseVector")))
setMethod("Arith", c(e1 = "dsparseVector", e2 = "dsparseVector"),
          Ops.spV.spV)

## "Arith"  exception (shortcut)
setMethod("-", c(e1 = "dsparseVector", e2 = "missing"),
          function(e1) { e1@x <- -e1@x ; e1 })


setMethod("Logic", c(e1 = "sparseVector", e2 = "sparseVector"),
          ## FIXME: this is suboptimal for "nsparseVector" !!
          function(e1, e2) callGeneric(as(e1, "lsparseVector"),
                                       as(e2, "lsparseVector")))

setMethod("Logic", c(e1 = "lsparseVector", e2 = "lsparseVector"),
          Ops.spV.spV)

## "nsparse" have no 'x' slot --> version of Ops.spV.spV..
## --------  but also for (nsp.. o lsp..) etc, when  lsp... has no NA
if(FALSE) ### FIXME
setMethod("Logic", c(e1 = "nsparseVector", e2 = "nsparseVector"),
          function(e1, e2) {
              .bail.out.2(.Generic, class(e1), class(e2))
          })

rm(Ops.spV.spV)

## 2)  spVec  o  [Mm]atrix : -------------

Ops.M.spV <- function(e1, e2) {
    d <- e1@Dim
    n1 <- prod(d)
    n2 <- e2@length
    if(n1 != n2) {
        if(n1 && n1 < n2) { # 0-extent matrix + vector is fine
            stop(sprintf("dim [product %d] do not match the length of object [%d]",
                         n1, n2))
        }
        ## else	 n1 > n2 [vector]
        N <- n1
        if(n2 == 1) ## simple case, do not really recycle
            return(callGeneric(e1, sp2vec(e2)))
        if(N %% n2 != 0)
            warning("longer object length is not a multiple of shorter object length")
        ## else : 2 <= n < N --- recycle the vector
        e2 <- rep(e2, length.out = N)
    } else { ## n1 == n2
        N <- n1
    }
    ## ---- e1 & e2 now are both of length 'N' ----
    dim(e2) <- d #-> sparseMatrix (!)
    callGeneric(e1, e2)
}## {Ops.M.spV}

Ops.spV.M <- function(e1, e2) {
    n1 <- e1@length
    d <- e2@Dim
    n2 <- prod(d)
    if(n2 != n1) {
        if(n2 && n2 < n1) { # vector + 0-extent matrix  is fine
            stop(sprintf("dim [product %d] do not match the length of object [%d]",
                         n2, n1))
        }
        ## else	 n2 > n1 [vector]
        N <- n2
        if(n1 == 1) ## simple case, do not really recycle
            return(callGeneric(sp2vec(e1), e2))
        if(N %% n1 != 0)
            warning("longer object length is not a multiple of shorter object length")
        ## else : 2 <= n < N --- recycle the vector
        e1 <- rep(e1, length.out = N)
    } else { ## n2 == n1
        N <- n2
    }
    ## ---- e2 & e1 now are both of length 'N' ----
    dim(e1) <- d #-> sparseMatrix (!)
    callGeneric(e1, e2)
}## {Ops.spV.M}

## try to use it in all cases
setMethod("Ops", c(e1 = "Matrix", e2 = "sparseVector"), Ops.M.spV)
setMethod("Ops", c(e1 = "sparseVector", e2 = "Matrix"), Ops.spV.M)

rm(Ops.M.spV, Ops.spV.M)


###---------------- diagonalMatrix ----------------------

.diag2tT.smart <- function(from, x, kind = ".") {
    shape <- .M.shape(x)
    uplo <- if(shape == "t") x@uplo else "U"
    .diag2sparse(from, kind, "t", "T", uplo)
}
.diag2T.smart <- function(from, x, kind = ".") {
    shape <- .M.shape(x)
    uplo <- if(shape == "s" || shape == "t") x@uplo else "U"
    .diag2sparse(from, kind, if(shape == "s") "s" else "t", "T", uplo)
}

 .diag.x <- function(m) if(m@diag != "N") rep.int(as1(m@x), m@Dim[1L]) else m@x
..diag.x <- function(m)                   rep.int(as1(m@x), m@Dim[1L])

## Use as S4 method for several signatures ==>  using callGeneric()
diagOdiag <- function(e1,e2) {
    ## result should also be diagonal _ if possible _
    r <- callGeneric(.diag.x(e1), .diag.x(e2)) # error if not "compatible"
    ## Check what happens with non-diagonals, i.e. (0 o 0), (FALSE o 0), ...:
    r00 <- callGeneric(as0(e1@x), as0(e2@x))
    if(is0(r00)) { ##  r00 == 0 or FALSE --- result *is* diagonal
        if(is.complex(r)) {
            if(is.complex(e2@x)) {
                e2@x <- if(is.complex(r)) r else as.complex(r)
                e2@diag <- "N"
                return(e2)
            }
            if(!is.complex(e1@x)) # e.g. logical;
                e1 <- .M2kind(e1, "z")
        }
        else if(is.double(r)) {
            if(is.double(e2@x)) {
                e2@x <- if(is.double(r)) r else as.double(r)
                e2@diag <- "N"
                return(e2)
            }
            if(!is.double(e1@x))
                e1 <- .M2kind(e1, "d")
        }
        else if(is.integer(r)) {
            if(is.integer(e2@x)) {
                e2@x <- if(is.integer(r)) r else as.integer(r)
                e2@diag <- "N"
                return(e2)
            }
            if(!is.integer(e1@x))
                e1 <- .M2kind(e1, "i")
        }
        else if(is.logical(r))
            e1 <- .M2kind(e1, "l")
        else stop(gettextf("intermediate 'r' is of type %s",
                           typeof(r)), domain=NA)
        e1@x <- r
        e1@diag <- "N"
        e1
    }
    else { ## result not diagonal, but at least symmetric
        ## e.g., m == m
        ## isInt <- (is.integer(r) && is.integer(r00)) # both (right ?)
        isNum <- (is.numeric(r) || is.numeric(r00))
        isLog <- (is.logical(r) || is.logical(r00))
        isClx <- (is.complex(r) || is.complex(r00))
        Matrix.message("exploding <diag> o <diag> into dense matrix", .M.level = 2)
        d <- e1@Dim
        n <- d[1L]
        stopifnot(length(r) == n)
        ## FIXME -- keep "iMatrix" eventually !
        if(isNum && !is.double(r) && !is.complex(r))
            r <- as.double(r)
        if(isClx && !is.complex(r))
            r <- as.complex(r)
        ## faster (?) than  m <- matrix(r00,n,n); diag(m) <- r ; as.vector(m)
        xx <- rbind(r, matrix(r00,n,n), deparse.level=0L)[seq_len(n*n)]
        newcl <-
            paste0(if(is.integer(xx)) "i"
                   else if(is.double(xx)) "d"
                   else if(is.complex(xx)) "z" ## TODO: symmetric or hermitian ?
                   else if(isLog) {
                       if(!anyNA(r) && !anyNA(r00)) "n" else "l"
                   } else stop("not yet implemented .. please report"), "syMatrix")

        if(isClx)
             new(newcl, Dim = e1@Dim, Dimnames = e1@Dimnames, x = xx, trans = "T")
        else new(newcl, Dim = e1@Dim, Dimnames = e1@Dimnames, x = xx)
    }
} # end  diagOdiag

### This would be *the* way, but we get tons of "ambiguous method dispatch"
diaTypes <- c("d","l","n", "i","z")
diCls <- paste0(diaTypes, "diMatrix")
if(FALSE) {
setMethod("Ops", c(e1 = "diagonalMatrix", e2 = "diagonalMatrix"),
          diagOdiag)
} else { ## These are just for method disambiguation:
    for(c1 in diCls)
        for(c2 in diCls)
            setMethod("Ops", c(e1 = c1, e2 = c2), diagOdiag)
    rm(c1, c2)
}
rm(diagOdiag)

## diagonal  o  triangular  |-->  triangular (upper/lower) ideally
## ~~~~~~~~ ~~~ ~~~~~~~~~~  when tri. is sparse, "work here",
##                          before conversion to sparse, since that loses "diagonality"
diagOtri <- function(e1, e2) {
    ## result to be triangular ("mostly")
    r <- callGeneric(d1 <- .diag.x(e1),
                     d2 <-  diag  (e2)) # error if not "compatible"
    ## Check what happens with non-diagonals, i.e. (0 o 0), (FALSE o 0), ...:
    r00 <- callGeneric(e1.0 <- as0(type = typeof(d1)),
                       e2.0 <- as0(type = typeof(d2)))
    rTri <- is0(r00) # r00 == 0 or FALSE  ==> result triangular !
    n <- e2@Dim[1L]
    upper <- e2@uplo == "U"
    knd <- .kind.type[[typeof(r00)]]
    if(!rTri) { # result not triangular but generalMatrix
        e2 <- .M2gen(e2, kind = knd)
    }
    diag(e2) <- r # the diagonal; now check what happens
    ## "in the off-diagonal (non-0) triangle";  e2.2 := as2(.) :
    e2.2 <- if(is.logical(e2.0)) TRUE else (t. <- as1(e2[0L])) + t.
    offD <- callGeneric(e1.0, e2.2) != e2.2 # off-diagonal values can change, too
    ## offD <- (r02 <- callGeneric(e1.0, e2.2)) != e2.2 # off-diagonal values can change, too
    if(rTri) { #
        if(offD) { # values "in triangle" can change, too:
            ## knd <- .kind.type[[typeof(r02)]] # never "n"
            if(knd != .M.kind(e2)) # e.g. "l" != "i";  (FIXME?  typeof(e2.2) instead of .M.kind(e2))
                e2 <- .M2kind(e2, knd) # no longer "n*" -> has 'x' slot
            if(.isPacked(e2)) {
                it <- indTri(n, upper = upper, diag = FALSE, packed = TRUE)
                e2@x[it] <- callGeneric(e1.0, e2@x[it])
            } else {
                it <- indTri(n, upper = upper, diag = FALSE, packed = FALSE)
                ## [<- is really calling into replTmat() -- does keep triang since 2025-03:
                e2[it] <- callGeneric(e1.0, e2[it])
            }
        }
    }
    else if(length(e2)) { # not triangular (nor trivial)
        ## non-zero triangle:
        it <- indTri(n, upper)
        e2[it] <- callGeneric(e1.0, e2[it])
        ## zero triangle off-diagonals:
        e2[indTri(n, !upper)] <- r00
    }
    e2
} # {diagOtri}

setMethod("Ops", c(e1 = "diagonalMatrix", e2 = "triangularMatrix"), diagOtri)
rm(diagOtri)

## For the reverse,  Ops == "Arith" | "Compare" | "Logic"  need to be split :
##   'Arith'  :=  '"+"', '"-"', '"*"', '"^"', '"%%"', '"%/%"', '"/"'
setMethod("Arith", c(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1, e2) { # e1 may be sparse _or_ dense
              switch(.Generic,
                     "+" = `diag<-`(e1, as.vector(diag(e1, names=FALSE) + .diag.x(e2))),
                     "-" = `diag<-`(e1, as.vector(diag(e1, names=FALSE) - .diag.x(e2))),
                     "*" = {
                         n <- e2@Dim[1L]
                         d2 <- if(e2@diag == "U") { # unit-diagonal
                                   d <- rep.int(as1(e2@x), n)
                                   e2@x <- d
                                   e2@diag <- "N"
                                   d
                               } else e2@x
                         e2@x <- diag(e1) * d2
                         e2
                     },
                     "^" = { ## will be dense ( as  <ANY> ^ 0 == 1 ):
                         e1 ^ .diag2dense(e2, ".", "g", FALSE)
                     },
                     ## otherwise:
                     callGeneric(e1, .diag2T.smart(e2, e1)))
          })

## Compare --> 'swap' (e.g.   e1 < e2   <==>  e2 > e1 ):
setMethod("Compare", c(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          .Cmp.swap)
## '&' and "|'  are commutative:
setMethod("Logic", c(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1, e2) callGeneric(e2, e1))

## For almost everything else, diag* shall be treated "as sparse" :
## These are cheap implementations via coercion

## For disambiguation --- define this for "sparseMatrix" , then for "ANY";
## and because we can save an .M.kind() call, we use this explicit
## "hack" for all diagonalMatrix *subclasses* instead of just "diagonalMatrix" :
##
## ddi*, ldi*, zdi*...
setMethod("Ops", c(e1 = "diagonalMatrix", e2 = "sparseMatrix"),
          function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = .M.kind(e1)), e2))
setMethod("Ops", c(e1 = "sparseMatrix", e2 = "diagonalMatrix"),
          function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind = .M.kind(e2))))

## Ops:	 Arith	--> numeric : "dMatrix"
##	 Compare --> logical
##	 Logic	 --> logical: "lMatrix"

atomV <- function(e1, e2, Generic, n = 0L) {
    eval(call(.type.kind[[.Arith.kind(e1,e2, Generic)]], n)) # double(0) etc
}
## TODO -- use below

## Other = "numeric" : stay diagonal if possible
## ddi*: Arith: result numeric, potentially ddiMatrix
for(arg2 in c("numeric","logical","complex"))
setMethod("Arith", c(e1 = "diagonalMatrix", e2 = arg2), # NB: *has* @x, even "ndiMatrix"
          function(e1,e2) {
              if(length(e2) == 0L)
                  return(if(length(e1)) vector(.type.kind[[.Arith.kind(e1,e2, .Generic)]], 0L)
                         else           .M2kind(e1,        .Arith.kind(e1,e2, .Generic)))
              f0 <- callGeneric(0, e2)
              if(all0(f0)) { # remain diagonal
                  Udg <- e1@diag == "U"
                  e1x <- if(Udg) as1(e1@x) else e1@x
                  n <- e1@Dim[1L]
                  if(n) { # n >= 1
                      r <- callGeneric(e1x, e2)
                      tx <- typeof(r)
                      if(tx != typeof(e1x))
                          e1 <- .M2kind(e1, .kind.type[[tx]])
                      if(Udg && any(r != 1)) {
                          e1@diag <- "N"
                          e1@x[seq_len(n)] <- r # possibly recycling r
                      }
                      else if(!Udg) {
                          L1 <- (le <- length(e2)) == 1L
                          e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                      } ## else: result = e1  with "U" diag
                  }
                  else { # n == 0
                      knd <- .Arith.kind(e1,e2, .Generic)
                      if(knd != .kind.type[[typeof(e1@x)]])
                          e1 <- .M2kind(e1, knd)
                  }
                  e1
              } else { # -> dense
                  r <- callGeneric(.M2unpacked(e1), e2)
                                        # .diag2tT.smart(e1, e2, kind = .M.kind(e1))
                  if(length(e2) <= 1L || .M.shape(r) == "s") .M2sym(r, trans="T") else r
              }
          })
rm(arg2)

for(arg1 in c("numeric","logical","complex"))
setMethod("Arith", c(e1 = arg1, e2 = "diagonalMatrix"),
          function(e1,e2) {
              if(length(e1) == 0L)
                  return(if(length(e2)) vector(.type.kind[[.Arith.kind(e1,e2, .Generic)]], 0L) # double(0) etc
                         else           .M2kind(e2,        .Arith.kind(e1,e2, .Generic)))
              f0 <- callGeneric(e1, 0)
              if(all0(f0)) { # remain diagonal
                  Udg <- e2@diag == "U"
                  e2x <- if(Udg) as1(e2@x) else e2@x
                  n <- e2@Dim[1L]
                  if(n) { # n >= 1
                      r <- callGeneric(e1, e2x)
                      tx <- typeof(r)
                      if(tx != typeof(e2x))
                          e2 <- .M2kind(e2, .kind.type[[tx]])
                      if(Udg && any(r != 1)) {
                          e2@diag <- "N"
                          e2@x[seq_len(n)] <- r # possibly recycling r
                      }
                      else if(!Udg) {
                          L1 <- (le <- length(e1)) == 1L
                          e2@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                      } ## else: result = e2  with "U" diag
                  }
                  else { # n == 0
                      knd <- .Arith.kind(e1,e2, .Generic)
                      if(knd != .kind.type[[typeof(e2@x)]])
                          e2 <- .M2kind(e2, knd)
                  }
                  e2
              } else { # -> dense
                  r <- callGeneric(e1, .M2unpacked(e2))
                                        # .diag2tT.smart(e1, e2, kind = .M.kind(e1))
                  if(length(e1) <= 1L || .M.shape(r) == "s") .M2sym(r, trans="T") else r
              }
          })
rm(arg1)

if(FALSE) { ## now covered above ??
## ldi* Arith --> result numeric, potentially ddiMatrix
for(arg2 in c("numeric","logical",))
setMethod("Arith", c(e1 = "ldiMatrix", e2 = arg2),
          function(e1,e2) {
              n <- e1@Dim[1L]
              if(length(e2) == 0L)
                  return(if(n) numeric()
                         else copyClass(e1, "ddiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE))
              f0 <- callGeneric(0, e2)
              if(all0(f0)) { # remain diagonal
                  E <- copyClass(e1, "ddiMatrix", c("diag", "Dim", "Dimnames"), check=FALSE)
                  ## storage.mode(E@x) <- "double"
                  if(e1@diag == "U") {
                      if(any((r <- callGeneric(1, e2)) != 1)) {
                          E@diag <- "N"
                          E@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = E  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e2)) == 1L
                      r <- callGeneric(e1@x, e2)
                      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "l"), e2)
          })
rm(arg2)

for(arg1 in c("numeric","logical"))
setMethod("Arith", c(e1 = arg1, e2 = "ldiMatrix"),
          function(e1,e2) {
              n <- e2@Dim[1L]
              if(length(e1) == 0L)
                  return(if(n) numeric()
                         else copyClass(e2, "ddiMatrix",
                                        c("diag", "Dim", "Dimnames"),
                                        check=FALSE))
              f0 <- callGeneric(e1, 0)
              if(all0(f0)) { # remain diagonal
                  E <- copyClass(e2, "ddiMatrix",
                                 c("diag", "Dim", "Dimnames"),
                                 check=FALSE)
                  ## storage.mode(E@x) <- "double"
                  if(e2@diag == "U") {
                      if(any((r <- callGeneric(e1, 1)) != 1)) {
                          E@diag <- "N"
                          E@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = E  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e1)) == 1L
                      r <- callGeneric(e1, e2@x)
                      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <-
                          if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(e1, .diag2tT.smart(e2, e1, kind = "l"))
          })
rm(arg1)
}## --- no longer; covered above

## ddi*: for "Ops" without "Arith": <Compare> or <Logic> --> result logical, potentially ldi
##
## Note that  ("numeric", "ddiMatrix")  is simply swapped, e.g.,
if(FALSE) {
    selectMethod("<", c("numeric","lMatrix"))# Compare
    selectMethod("&", c("numeric","lMatrix"))# Logic
} ## so we don't need to define a method here :

if(FALSE) { ## now covered above ??
## __ FIXME __ no longer ??? __ otherwise not only for "ddiMatrix" and then "ldiMatrix" !! ___
for(arg2 in c("numeric","logical"))
setMethod("Ops", c(e1 = "ddiMatrix", e2 = arg2),
          function(e1,e2) {
              n <- e1@Dim[1L]
              if(length(e2) == 0L)
                  return(if(n) logical()
                         else copyClass(e1, "ldiMatrix",
                                        c("diag", "Dim", "Dimnames"),
                                        check=FALSE))
              f0 <- callGeneric(0, e2)
              if(all0(f0)) { # remain diagonal
                  E <- copyClass(e1, "ldiMatrix",
                                 c("diag", "Dim", "Dimnames"),
                                 check=FALSE)
                  ## storage.mode(E@x) <- "logical"
                  if(e1@diag == "U") {
                      if(any((r <- callGeneric(1, e2)) != 1)) {
                          E@diag <- "N"
                          E@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = E  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e2)) == 1L
                      r <- callGeneric(e1@x, e2)
                      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <-
                          if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "d"), e2)
          })
rm(arg2)

## __ FIXME __ no longer ??? (see "ddiMatrix" above)
## ldi*: for "Ops" without "Arith": <Compare> or <Logic> --> result logical, potentially ldi
for(arg2 in c("numeric","logical"))
setMethod("Ops", c(e1 = "ldiMatrix", e2 = arg2),
          function(e1,e2) {
              n <- e1@Dim[1L]
              if(length(e2) == 0L)
                  return(if(n) logical() else e1)
              f0 <- callGeneric(FALSE, e2)
              if(all0(f0)) { # remain diagonal
                  if(e1@diag == "U") {
                      if(any((r <- callGeneric(TRUE, e2)) != 1)) {
                          e1@diag <- "N"
                          e1@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = e1  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e2)) == 1L
                      r <- callGeneric(e1@x, e2)
                      tx <- typeof(r)
                      if(tx != typeof(e1@x))
                          e1 <- .M2kind(e1, .kind.type[[tx]])
                      e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  e1
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "l"), e2)
          })
rm(arg2)
}## --- no longer; covered above


## Not {"sparseMatrix", "numeric} :  {"denseMatrix", "matrix", ... }
for(other in c("ANY", "Matrix", "dMatrix")) {
    setMethod("Ops", c(e1 = "diagonalMatrix", e2 = other),
              function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = .M.kind(e1)), e2))
    setMethod("Ops", c(e1 = other, e2 = "diagonalMatrix"),
              function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind = .M.kind(e2))))
}
rm(other)

## Direct subclasses of "denseMatrix":
## diaTypes <- c("d","l","n", "i","z")
dense.subCl <- paste0(diaTypes, "denseMatrix")
diCls       <- paste0(diaTypes, "diMatrix")
for(typ in diaTypes) { # "d" "l" "n" "i" "z"
    DI <- paste0(typ, "diMatrix")
    dMeth <- switch(typ,
             "d" = function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "d"), e2),
             "l" =,
             "n" = function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "l"), e2),
             "i" = function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "i"), e2),
             "z" = function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "z"), e2),
             stop("invalid 'typ'"))

    for(c2 in c(dense.subCl, "Matrix")) {
        for(Fun in c("*", "&")) {
            setMethod(Fun, c(e1 = DI, e2 = c2),
                      function(e1,e2) callGeneric(e1, Diagonal(x = diag(e2))))
            setMethod(Fun, c(e1 = c2, e2 = DI),
                      function(e1,e2) callGeneric(Diagonal(x = diag(e1)), e2))
        }
        setMethod("^", c(e1 = c2, e2 = DI),
                  function(e1,e2) callGeneric(Diagonal(x = diag(e1)), e2))
        for(Fun in c("%%", "%/%", "/")) ## 0 <op> 0 |--> NaN  for these.
            setMethod(Fun, c(e1 = DI, e2 = c2), dMeth)
    }
}
rm(diaTypes, dense.subCl, diCls, typ, DI, dMeth, c2, Fun)
