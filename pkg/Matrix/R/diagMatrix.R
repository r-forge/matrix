## METHODS FOR CLASS: diagonalMatrix (virtual)
## diagonal matrices
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## For group methods
.diag2tT.smart <- function(from, x, kind = ".") {
    shape <- .M.shape(x)
    uplo <- if(shape == "t") x@uplo else "U"
    .diag2sparse(.M2kind(from, kind), "t", "T", uplo)
}
.diag2T.smart <- function(from, x, kind = ".") {
    shape <- .M.shape(x)
    uplo <- if(shape == "s" || shape == "t") x@uplo else "U"
    .diag2sparse(.M2kind(from, kind), if(shape == "s") "s" else "t", "T", uplo)
}

 .diag.x <- function(m) if(m@diag != "N") rep.int(as1(m@x), m@Dim[1L]) else m@x
..diag.x <- function(m)                   rep.int(as1(m@x), m@Dim[1L])

setMethod("band", signature(x = "diagonalMatrix"),
          function(x, k1, k2, ...) {
              if(k1 <= 0L && k2 >= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("triu", signature(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if(k <= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("tril", signature(x = "diagonalMatrix"),
          function(x, k = 0L, ...) {
              if(k >= 0L)
                  return(x)
              r <- new(.M.nonvirtual(x))
              r@Dim <- d <- x@Dim
              r@Dimnames <- x@Dimnames
              r@x <- vector(typeof(x@x), d[1L])
              r
          })

setMethod("diag", signature(x = "diagonalMatrix"),
          function(x, nrow, ncol, names = TRUE) {
              r <- .diag.x(x)
              if(names &&
                 !any(vapply(dn <- x@Dimnames, is.null, NA)) &&
                 {
                     i <- seq_len(min(x@Dim))
                     identical(nms <- dn[[1L]][i], dn[[2L]][i])
                 })
                  names(r) <- nms
              r
          })

setMethod("diag<-", signature(x = "diagonalMatrix"),
          function(x, value) {
              n <- x@Dim[1L]
              nv <- length(value)
              if(nv != 1L && nv != n)
                  stop("replacement diagonal has wrong length")
              x@x <-
                  if(is.logical(x@x))
                      switch(typeof(value),
                             logical = rep_len(value, n),
                             integer =,
                             double =
                                 {
                                     x <- .M2kind(x, "d")
                                     rep_len(as.double(x), n)
                                 },
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"", typeof(value)),
                                  domain = NA))
                  else
                      switch(typeof(value),
                             logical =,
                             integer =,
                             double = rep_len(as.double(value), n),
                             stop(gettextf("replacement diagonal has incompatible type \"%s\"", typeof(value)),
                                  domain = NA))
              x@diag <- "N"
              x
          })

setMethod("t", signature(x = "diagonalMatrix"),
          function(x) { x@Dimnames <- x@Dimnames[2:1]; x })

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "missing"),
          function(x, uplo) .diag2sparse(x, "s", "C",  "U"))

setMethod("forceSymmetric", signature(x = "diagonalMatrix", uplo = "character"),
          function(x, uplo) .diag2sparse(x, "s", "C", uplo))

setMethod("symmpart", signature(x = "diagonalMatrix"),
          function(x) {
              kind <- .M.kind(x)
              r <- new(if(kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- x@Dim
              r@Dimnames <- symmDN(x@Dimnames)
              if(x@diag != "N")
                  r@diag <- "U"
              else
                  r@x <- switch(kind,
                                "n" = as.double(x@x | is.na(x@x)),
                                "l" = ,
                                "i" = ,
                                "d" = as.double(x@x),
                                "z" = complex(real = Re(x@x), imaginary = 0))
              r
          })

setMethod("skewpart", signature(x = "diagonalMatrix"),
          function(x) {
              kind <- .M.kind(x)
              r <- new(if(kind != "z") "ddiMatrix" else "zdiMatrix")
              r@Dim <- d <- x@Dim
              r@Dimnames <- symmDN(x@Dimnames)
              r@x <- if(kind != "z")
                         double(d[1L])
                     else if(x@diag != "N")
                         complex(d[1L])
                     else complex(real = 0, imaginary = Im(x@x))
              r
          })

setMethod("isDiagonal", signature(object = "diagonalMatrix"),
          function(object) TRUE)

setMethod("isTriangular", signature(object = "diagonalMatrix"),
          function(object, upper = NA, ...)
              if(is.na(upper)) `attr<-`(TRUE, "kind", "U") else TRUE)

setMethod("isSymmetric", signature(object = "diagonalMatrix"),
          function(object, checkDN = TRUE, ...) {
              if(checkDN) {
                  ca <- function(check.attributes = TRUE, ...)
                      check.attributes
                  if(ca(...) && !isSymmetricDN(object@Dimnames))
                      return(FALSE)
              }
              TRUE
          })


###---------------- <Ops> (<Arith>, <Logic>, <Compare> ) ----------------------

## Use as S4 method for several signatures ==>  using callGeneric()
diagOdiag <- function(e1,e2) {
    ## result should also be diagonal _ if possible _
    r <- callGeneric(.diag.x(e1), .diag.x(e2)) # error if not "compatible"
    ## Check what happens with non-diagonals, i.e. (0 o 0), (FALSE o 0), ...:
    r00 <- callGeneric(if(is.numeric(e1@x)) 0 else FALSE,
                       if(is.numeric(e2@x)) 0 else FALSE)
    if(is0(r00)) { ##  r00 == 0 or FALSE --- result *is* diagonal
        if(is.numeric(r)) { # "double" *or* "integer"
            if(!is.double(r))
                r <- as.double(r)
            if(is.double(e2@x)) {
                e2@x <- r
                e2@diag <- "N"
                return(e2)
            }
            if(!is.double(e1@x))
                ## e.g. e1, e2 are logical;
                e1 <- .M2kind(e1, "d")
        }
        else if(is.logical(r))
            e1 <- .M2kind(e1, "l")
        else stop(gettextf("intermediate 'r' is of type %s",
                           typeof(r)), domain=NA)
        e1@x <- r
        e1@diag <- "N"
        e1
    }
    else { ## result not diagonal, but at least symmetric:
        ## e.g., m == m
        isNum <- (is.numeric(r) || is.numeric(r00))
        isLog <- (is.logical(r) || is.logical(r00))
        Matrix.message("exploding <diag> o <diag> into dense matrix", .M.level = 2)
        d <- e1@Dim
        n <- d[1L]
        stopifnot(length(r) == n)
        if(isNum && !is.double(r))
            r <- as.double(r)
        ## faster (?) than  m <- matrix(r00,n,n); diag(m) <- r ; as.vector(m)
        xx <- rbind(r, matrix(r00,n,n), deparse.level=0L)[seq_len(n*n)]
        newcl <-
            paste0(if(isNum) "d"
                   else if(isLog) {
                       if(!anyNA(r) && !anyNA(r00)) "n" else "l"
                   } else stop("not yet implemented .. please report"), "syMatrix")

        new(newcl, Dim = e1@Dim, Dimnames = e1@Dimnames, x = xx)
    }
}

### This would be *the* way, but we get tons of "ambiguous method dispatch"
## we use this hack instead of signature  x = "diagonalMatrix" :
diCls <- names(getClassDef("diagonalMatrix")@subclasses)
if(FALSE) {
setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "diagonalMatrix"),
          diagOdiag)
} else { ## These are just for method disambiguation:
    for(c1 in diCls)
        for(c2 in diCls)
            setMethod("Ops", signature(e1 = c1, e2 = c2), diagOdiag)
    rm(c1, c2)
}
rm(diagOdiag)

## diagonal  o  triangular  |-->  triangular
## diagonal  o  symmetric   |-->  symmetric
##    {also when other is sparse: do these "here" --
##     before conversion to sparse, since that loses "diagonality"}
diagOtri <- function(e1,e2) {
    ## result must be triangular
    r <- callGeneric(d1 <- .diag.x(e1), diag(e2)) # error if not "compatible"
    ## Check what happens with non-diagonals, i.e. (0 o 0), (FALSE o 0), ...:
    e1.0 <- if(is.numeric(d1)) 0 else FALSE
    r00 <- callGeneric(e1.0, if(.n2 <- is.numeric(e2[0L])) 0 else FALSE)
    if(is0(r00)) { ##  r00 == 0 or FALSE --- result *is* triangular
        diag(e2) <- r
        ## check what happens "in the triangle"
        e2.2 <- if(.n2) 2 else TRUE
        if(!callGeneric(e1.0, e2.2) == e2.2) { # values "in triangle" can change:
            n <- dim(e2)[1L]
            it <- indTri(n, upper = (e2@uplo == "U"))
            e2[it] <- callGeneric(e1.0, e2[it])
        }
        e2
    }
    else { ## result not triangular ---> general
        rr <- as(e2, "generalMatrix")
        diag(rr) <- r
        rr
    }
}


setMethod("Ops", signature(e1 = "diagonalMatrix", e2 = "triangularMatrix"),
          diagOtri)
rm(diagOtri)

## For the reverse,  Ops == "Arith" | "Compare" | "Logic"
##   'Arith'  :=  '"+"', '"-"', '"*"', '"^"', '"%%"', '"%/%"', '"/"'
setMethod("Arith", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1, e2) { ## this must only trigger for *dense* e1
              switch(.Generic,
                     "+" = `diag<-`(e1, as.double(diag(e1, names=FALSE) + .diag.x(e2))),
                     "-" = `diag<-`(e1, as.double(diag(e1, names=FALSE) - .diag.x(e2))),
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
                         e1 ^ .diag2dense(e2, "g", FALSE)
                     },
                     ## otherwise:
                     callGeneric(e1, .diag2T.smart(e2, e1)))
          })

## Compare --> 'swap' (e.g.   e1 < e2   <==>  e2 > e1 ):
setMethod("Compare", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          .Cmp.swap)
## '&' and "|'  are commutative:
setMethod("Logic", signature(e1 = "triangularMatrix", e2 = "diagonalMatrix"),
          function(e1, e2) callGeneric(e2, e1))

## For almost everything else, diag* shall be treated "as sparse" :
## These are cheap implementations via coercion

## For disambiguation --- define this for "sparseMatrix" , then for "ANY";
## and because we can save an .M.kind() call, we use this explicit
## "hack" for all diagonalMatrix *subclasses* instead of just "diagonalMatrix" :
##
## ddi*:
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = "sparseMatrix"),
          function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "d"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "ddiMatrix"),
          function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind = "d")))
## ldi*
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = "sparseMatrix"),
          function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "l"), e2))
setMethod("Ops", signature(e1 = "sparseMatrix", e2 = "ldiMatrix"),
          function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind = "l")))

## Ops:	 Arith	--> numeric : "dMatrix"
##	 Compare --> logical
##	 Logic	 --> logical: "lMatrix"

## Other = "numeric" : stay diagonal if possible
## ddi*: Arith: result numeric, potentially ddiMatrix
for(arg2 in c("numeric","logical"))
setMethod("Arith", signature(e1 = "ddiMatrix", e2 = arg2),
          function(e1,e2) {
              n <- e1@Dim[1L]
              if(length(e2) == 0L)
                  return(if(n) numeric() else e1)
              f0 <- callGeneric(0, e2)
              if(all0(f0)) { # remain diagonal
                  if(e1@diag == "U") {
                      if(any((r <- callGeneric(1, e2)) != 1)) {
                          e1@diag <- "N"
                          e1@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = e1  (is "U" diag)
                  } else if(n) {
                      L1 <- (le <- length(e2)) == 1L
                      r <- callGeneric(e1@x, e2)
                      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
                      e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  e1
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "d"), e2)
          })
rm(arg2)

for(arg1 in c("numeric","logical"))
setMethod("Arith", signature(e1 = arg1, e2 = "ddiMatrix"),
          function(e1,e2) {
              n <- e2@Dim[1L]
              if(length(e1) == 0L)
                  return(if(n) numeric() else e2)
              f0 <- callGeneric(e1, 0)
              if(all0(f0)) { # remain diagonal
                  if(e2@diag == "U") {
                      if(any((r <- callGeneric(e1, 1)) != 1)) {
                          e2@diag <- "N"
                          e2@x[seq_len(n)] <- r # possibly recycling r
                      } ## else: result = e2  (is "U" diag)
                  } else {
                      L1 <- (le <- length(e1)) == 1L
                      r <- callGeneric(e1, e2@x)
                      ## "future fixme": if we have idiMatrix, and r is 'integer', use idiMatrix
                      e2@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  e2
              } else
                  callGeneric(e1, .diag2tT.smart(e2, e1, kind = "d"))
          })
rm(arg1)

## ldi* Arith --> result numeric, potentially ddiMatrix
for(arg2 in c("numeric","logical"))
setMethod("Arith", signature(e1 = "ldiMatrix", e2 = arg2),
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
setMethod("Arith", signature(e1 = arg1, e2 = "ldiMatrix"),
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
                      ## "future fixme": if we have idiMatrix,
                      ## and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <-
                          if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(e1, .diag2tT.smart(e2, e1, kind = "l"))
          })
rm(arg1)

## ddi*: for "Ops" without "Arith": <Compare> or <Logic> --> result logical, potentially ldi
##
## Note that  ("numeric", "ddiMatrix")  is simply swapped, e.g.,
if(FALSE) {
    selectMethod("<", c("numeric","lMatrix"))# Compare
    selectMethod("&", c("numeric","lMatrix"))# Logic
} ## so we don't need to define a method here :

for(arg2 in c("numeric","logical"))
setMethod("Ops", signature(e1 = "ddiMatrix", e2 = arg2),
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
                      ## "future fixme": if we have idiMatrix,
                      ### and r is 'integer', use idiMatrix
                      E@x[seq_len(n)] <-
                          if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  E
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "d"), e2)
          })
rm(arg2)

## ldi*: for "Ops" without "Arith": <Compare> or <Logic> --> result logical, potentially ldi
for(arg2 in c("numeric","logical"))
setMethod("Ops", signature(e1 = "ldiMatrix", e2 = arg2),
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
                      ## "future fixme": if we have idiMatrix,
                      ## and r is 'integer', use idiMatrix
                      e1@x[] <- if(L1) r else r[1L + ((n+1)*(0:(n-1L))) %% le]
                  }
                  e1
              } else
                  callGeneric(.diag2tT.smart(e1, e2, kind = "l"), e2)
          })
rm(arg2)

## Not {"sparseMatrix", "numeric} :  {"denseMatrix", "matrix", ... }
for(other in c("ANY", "Matrix", "dMatrix")) {
    ## ddi*:
    setMethod("Ops", signature(e1 = "ddiMatrix", e2 = other),
              function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind="d"), e2))
    setMethod("Ops", signature(e1 = other, e2 = "ddiMatrix"),
              function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind="d")))
    ## ldi*:
    setMethod("Ops", signature(e1 = "ldiMatrix", e2 = other),
              function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind="l"), e2))
    setMethod("Ops", signature(e1 = other, e2 = "ldiMatrix"),
              function(e1,e2) callGeneric(e1, .diag2T.smart(e2, e1, kind="l")))
}
rm(other)

## Direct subclasses of "denseMatrix": currently ddenseMatrix, ldense... :
if(FALSE) # now also contains "geMatrix"
dense.subCl <- local({ dM.scl <- getClassDef("denseMatrix")@subclasses
    names(dM.scl)[vapply(dM.scl, slot, 0, "distance") == 1] })
dense.subCl <- paste0(c("d","l","n"), "denseMatrix")
for(DI in diCls) {
    dMeth <-
        if(extends(DI, "dMatrix"))
            function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "d"), e2)
        else # "lMatrix", the only other kind for now
            function(e1,e2) callGeneric(.diag2T.smart(e1, e2, kind = "l"), e2)
    for(c2 in c(dense.subCl, "Matrix")) {
        for(Fun in c("*", "&")) {
            setMethod(Fun, signature(e1 = DI, e2 = c2),
                      function(e1,e2) callGeneric(e1, Diagonal(x = diag(e2))))
            setMethod(Fun, signature(e1 = c2, e2 = DI),
                      function(e1,e2) callGeneric(Diagonal(x = diag(e1)), e2))
        }
        setMethod("^", signature(e1 = c2, e2 = DI),
                  function(e1,e2) callGeneric(Diagonal(x = diag(e1)), e2))
        for(Fun in c("%%", "%/%", "/")) ## 0 <op> 0 |--> NaN  for these.
            setMethod(Fun, signature(e1 = DI, e2 = c2), dMeth)
    }
}
rm(dense.subCl, DI, dMeth, c2, Fun)

## Group methods "Math", "Math2" in			--> ./Math.R

### "Summary" : "max"   "min"   "range" "prod"  "sum"   "any"   "all"
### ----------   the last 4: separately here
for(cl in diCls) {
setMethod("any", cl,
          function (x, ..., na.rm) {
              if(any(x@Dim == 0)) FALSE
              else if(x@diag == "U") TRUE else any(x@x, ..., na.rm = na.rm)
          })
setMethod("all",  cl,
          function (x, ..., na.rm) {
              n <- x@Dim[1L]
              if(n >= 2) FALSE
              else if(n == 0 || x@diag == "U") TRUE
              else all(x@x, ..., na.rm = na.rm)
          })
setMethod("prod", cl,
          function (x, ..., na.rm) {
              n <- x@Dim[1L]
              if(n >= 2) 0
              else if(n == 0 || x@diag == "U") 1
              else ## n == 1, diag = "N" :
                  prod(x@x, ..., na.rm = na.rm)
          })
setMethod("sum", cl,
          function(x, ..., na.rm) {
              r <- sum(x@x, ..., na.rm = na.rm)# double or integer, correctly
              if(x@diag == "U" && !is.na(r)) r + x@Dim[1L] else r
          })
}
rm(cl, diCls)

## The remaining ones are  max, min, range :

setMethod("Summary", "ddiMatrix",
          function(x, ..., na.rm) {
              if(any(x@Dim == 0)) callGeneric(numeric(0), ..., na.rm=na.rm)
              else if(x@diag == "U")
                  callGeneric(x@x, 0, 1, ..., na.rm=na.rm)
              else callGeneric(x@x, 0, ..., na.rm=na.rm)
          })
setMethod("Summary", "ldiMatrix",
          function(x, ..., na.rm) {
              if(any(x@Dim == 0)) callGeneric(logical(0), ..., na.rm=na.rm)
              else if(x@diag == "U")
                  callGeneric(x@x, FALSE, TRUE, ..., na.rm=na.rm)
              else callGeneric(x@x, FALSE, ..., na.rm=na.rm)
          })
