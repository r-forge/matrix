#### ------- All "Ops"  group methods for all Matrix classes ------------
###               ===
### Note that the "Ops" group consists of
### sub-groups   "Arith", "Compare", and "Logic"
###               -----    -------        -----
### where 'Arith'   :=  '"+"', '"-"', '"*"', '"^"', '"%%"', '"%/%"', '"/"'
###       'Compare' := '"=="', '">"', '"<"', '"!="', '"<="', '">="'
###       'Logic'   :=  '"&"', '"|"'  (( but *not* '"!"' since that has
###			                 only one argument ))

## in shell, find them with
##    egrep 'Method\("(Ops|Compare|Arith|Logic)"' R/*R

### "Ops" ---- remember Ops = {Arith, Compare, Logic}  [Logic: since R 2.4.1]
### -----


### Note: diagonalMatrix are handled by special methods -> ./diagMatrix.R
###       --------------                                   ~~~~~~~~~~~~~~


### --  0 -- (not dense *or* sparse) -----------------------------------

##-------- originally from ./Matrix.R --------------------

## Some ``Univariate'' "Arith":
setMethod("+", signature(e1 = "Matrix", e2 = "missing"), function(e1) e1)
## "fallback":
setMethod("-", signature(e1 = "Matrix", e2 = "missing"),
          function(e1) {
              warning("inefficient method used for \"- e1\"")
              0-e1
          })

## old-style matrices are made into new ones
setMethod("Ops", signature(e1 = "Matrix", e2 = "matrix"),
	  function(e1, e2) callGeneric(e1, Matrix(e2)))
##	    callGeneric(e1, Matrix(e2, sparse=is(e1,"sparseMatrix"))))
setMethod("Ops", signature(e1 = "matrix", e2 = "Matrix"),
	  function(e1, e2) callGeneric(Matrix(e1), e2))

## bail-outs -- on highest possible level, hence "Ops", not "Compare"/"Arith" :
setMethod("Ops", signature(e1 = "Matrix", e2 = "Matrix"),
          function(e1, e2) {
              d <- dimCheck(e1,e2)
              .bail.out.2(.Generic, class(e1), class(e2))
          })
setMethod("Ops", signature(e1 = "Matrix", e2 = "ANY"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))
setMethod("Ops", signature(e1 = "ANY", e2 = "Matrix"),
          function(e1, e2) .bail.out.2(.Generic, class(e1), class(e2)))




##-------- originally from ./dMatrix.R --------------------

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  signature(e1 = "dMatrix", e2 = "dMatrix"),
          ## Going -> dense* (= ddense*) -> dgeMatrix
	  function(e1, e2) {
	      d <- dimCheck(e1,e2)
	      callGeneric(as(e1, "denseMatrix"),
			  as(e2, "denseMatrix"))
	  })

## "Compare" -> returning  logical Matrices
setMethod("Compare", signature(e1 = "numeric", e2 = "dMatrix"),
	  function(e1,e2) {
	      ## "swap RHS and LHS" and use the method below:
	      switch(.Generic,
		     "==" =, "!=" = callGeneric(e2, e1),
		     "<"  = e2 >  e1,
		     "<=" = e2 >= e1,
		     ">"  = e2 <  e1,
		     ">=" = e2 <= e1)
	  })

setMethod("Compare", signature(e1 = "dMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      lClass <- class2(class(e1), "l")
	      fullCl <- if(isSymmetric(e1)) "lsyMatrix" else "lgeMatrix"
	      ## Dbg cat("Compare", class(e1), "|-> ",lClass, "\n")
	      r	 <- callGeneric(e1@x, e2)
	      r0 <- callGeneric(0, e2)
              d <- e1@Dim
	      ## trivial case first
	      if(isTRUE(r0) && all(r)) {
		  r <- new(fullCl)
		  r@Dim <- d
		  r@Dimnames <- e1@Dimnames
		  r@x <- rep.int(TRUE, prod(d))
	      }
	      else if(is(e1, "denseMatrix")) {
		  full <- !isPacked(e1) # << both "dtr" and "dsy" are 'full'
		  if(full || identical(r0, FALSE) || is(e1, "symmetricMatrix"))
		      r <- new(lClass, x = r, Dim = d, Dimnames = dimnames(e1))
		  else { ## packed matrix with structural 0 and r0 is not FALSE:
		      ##--> result cannot be packed anymore
                      ## [dense & packed & not symmetric ] ==> must be "dtp*" :
                      if(!is(e1, "dtpMatrix"))
                          stop("internal bug in \"Compare\" method for \"dMatrix\"; please report")
                      rx <- rep.int(r0, d[1]*d[2])
                      rx[indTri(d[1], upper = (e1@uplo == "U"))] <- r
                      r <- new(fullCl, x = rx, Dim = d, Dimnames = dimnames(e1))
		  }

	      }
	      else { ## dsparseMatrix => lClass is "lsparse*"

		  if(identical(r0, FALSE)) { ## things remain sparse
		      if(!any(is.na(r)) && ((Ar <- all(r)) || !any(r))) {
			  r <- new(lClass)
			  r@Dim <- d
			  r@Dimnames <- dimnames(e1)
			  if(Ar) { # 'TRUE' instead of 'x': same sparsity:
			      r@x <- rep.int(TRUE, length(e1@x))
			      for(n in intersect(c("i","j","p"), slotNames(r)))
				  slot(r, n) <- slot(e1, n)
                          }
			  ## else: all FALSE: keep empty 'r' matrix
		      } else { # some TRUE, FALSE, NA : go via unique 'Tsparse'
			  M <- asTuniq(e1)
			  nCl <- class2(class(M), 'l') # logical Tsparse
			  r <- new(nCl)
			  r@x <- callGeneric(M@x, e2)
			  ## copy "the other slots" (important for "tr"/"sym"):
			  ## "%w/o%" <- function(x,y) x[is.na(match(x, y))]
			  sN <- slotNames(nCl)
			  for(n in sN[is.na(match(sN, "x"))])
			      slot(r, n) <- slot(M, n)
			  if(is(e1, "CsparseMatrix"))
			      r <- as(r, "CsparseMatrix")
			  else if(is(e1, "RsparseMatrix"))
			      r <- as(r, "RsparseMatrix")
		      }
		  } else {
		      ## non sparse result
		      message(sprintf("sparse to dense (%s) coercion in '%s'",
				      lClass, .Generic))
		      rx <- rep.int(r0, d[1]*d[2])
		      if(isTriangular(e1) && e1@diag == "U")
			  r <- c(r, rep.int(callGeneric(1, e2),d[1]))
		      rx[1:1 + encodeInd(non0ind(e1), nr = d[1])] <- r
		      r <- new(fullCl, x = rx, Dim = d, Dimnames = dimnames(e1))
		  }
	      }
	      r
	  })

## "dMatrix <-> work with 'x' slot
## FIXME? use 'Ops' and not just 'Compare' :
setMethod("Compare", signature(e1 = "dMatrix", e2 = "dMatrix"),
	  function(e1, e2) {
	      d <- dimCheck(e1,e2)
	      if((dens1 <- is(e1, "denseMatrix"))) gen1 <- is(e1, "generalMatrix")
	      if((dens2 <- is(e2, "denseMatrix"))) gen2 <- is(e2, "generalMatrix")

	      if(dens1 && dens2) { ## both inherit from ddense*

		  if(!gen1) e1 <- as(e1, "dgeMatrix")
		  if(!gen2) e2 <- as(e2, "dgeMatrix")
		  ## now, both are dge {ddense* & general*}

		  r <- new("lgeMatrix", x = callGeneric(e1@x, e2@x),
			   Dim = d, Dimnames = dimnames(e1))
	      }
	      else {
		  if(!dens1 && !dens2) {
		      ## both e1 _and_ e2 are sparse
		      ## should not happen since we have <sparse> o <sparse> methods
		      stop("Mistaken intended method dispatch -- please report to ",
			   packageDescription("Matrix")$Author)
		  }
		  ## else
		  if(dens1 && !dens2) ## go to dense
		      r <- callGeneric(e1, as(e2, "denseMatrix"))
		  else ## if(!dens1 && dens2)
		      r <- callGeneric(as(e1, "denseMatrix"), e2)

		  ## criterion "2 * nnz(.) < ." as in sparseDefault() in Matrix()  [./Matrix.R] :
		  if(2 * nnzero(r, na.counted = TRUE) < prod(d))
		      r <- as(r, "sparseMatrix")
	      }
	      r
	  })

### --  I -- dense -----------------------------------------------------------

##-------- originally from ./dgeMatrix.R --------------------

## ----- only work with NAMESPACE importFrom(methods, ..)

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  signature(e1 = "dgeMatrix", e2 = "dgeMatrix"),
	  function(e1, e2) {
	      ## NB:  triangular, symmetric, etc may need own method
	      d1 <- e1@Dim
	      d2 <- e2@Dim
	      eqD <- d1 == d2
	      if (!eqD[1])
		  stop("Matrices must have same number of rows for arithmetic")
	      same.dim <- eqD[2]
	      if (same.dim) {
		  d <- d1
		  dn <- dimNamesCheck(e1, e2)
	      }
	      else { # nrows differ ----> maybe recycling
		  if(d2[2] %% d1[2] == 0) { # nrow(e2) is a multiple
		      e1@x <- rep.int(e1@x, d2[2] %/% d1[2])
		      d <- d2
		      dn <- e2@Dimnames
		  } else if(d1[2] %% d2[2] == 0) { # nrow(e1) is a multiple
		      e2@x <- rep.int(e2@x, d1[2] %/% d2[2])
		      d <- d1
		      dn <- e1@Dimnames
		  } else
		      stop("number of rows are not compatible for ", .Generic)
	      }

	      ## be smart and preserve, e.g., triangular, or symmetric
	      ## but this sucks: For these,
	      ## 'uplo' and 'diag' also must coincide or be dealt with properly

	      ## ==> triangular, symmetric, etc may need own method
	      ##     also since their @x is `non-typical'

##		 if(same.dim) {
##		     if(extends(class(e1), class(e2))) {
##			 e2@x <- callGeneric(e1@x, e2@x)
##			 e2@Dimnames <- dn
##			 e2
##		     }
##		     else if(extends(class(e2), class(e1))) {
##			 e1@x <- callGeneric(e1@x, e2@x)
##			 e1@Dimnames <- dn
##			 e1
##		     }
##		 }
##		 else
		  new("dgeMatrix", Dim = d, Dimnames = dn,
		      x = callGeneric(e1@x, e2@x))
	  })

setMethod("Arith",
	  signature(e1 = "dgeMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      d <- e1@Dim
	      le <- length(e2)
	      if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
		  e1@x <- callGeneric(e1@x, as.vector(e2))
		  e1
	      } else stop ("length of 2nd arg does not match dimension of first")
	  })

setMethod("Arith",
	  signature(e1 = "numeric", e2 = "dgeMatrix"),
	  function(e1, e2) {
	      d <- e2@Dim
	      le <- length(e1)
	      if(le == 1 || le == d[1] || prod(d) == le) { # matching dim
		  e2@x <- callGeneric(as.vector(e1), e2@x)
		  e2
	      } else stop ("length of 1st arg does not match dimension of 2nd")
	  })

##-------- originally from ./ddenseMatrix.R --------------------

## Cheap version: work via "dgeMatrix" and use the group methods there:
## FIXME(?): try to preserve "symmetric", "triangular", ...
setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
          signature(e1 = "ddenseMatrix", e2 = "ddenseMatrix"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"),
                                       as(e2, "dgeMatrix")))
setMethod("Arith",
          signature(e1 = "ddenseMatrix", e2 = "numeric"),
          function(e1, e2) callGeneric(as(e1, "dgeMatrix"), e2))
setMethod("Arith",
          signature(e1 = "numeric", e2 = "ddenseMatrix"),
          function(e1, e2) callGeneric(e1, as(e2, "dgeMatrix")))

## "Logic"
## -------

##-------- originally from ./ldenseMatrix.R --------------------

setMethod("Logic", signature(e1="lgeMatrix", e2="lgeMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      e1@x <- callGeneric(e1@x, e2@x)
	      e1
	  })

setMethod("Logic", signature(e1="ldenseMatrix", e2="ldenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      callGeneric(as(e1, "lgeMatrix"), as(e2, "lgeMatrix"))
	  })

##-------- originally from ./ndenseMatrix.R --------------------

setMethod("Logic", signature(e1="ngeMatrix", e2="ngeMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      e1@x <- callGeneric(e1@x, e2@x)
	      e1
	  })

setMethod("Logic", signature(e1="ndenseMatrix", e2="ndenseMatrix"),
	  function(e1,e2) {
	      d <- dimCheck(e1, e2)
	      callGeneric(as(e1, "ngeMatrix"), as(e2, "ngeMatrix"))
	  })



### -- II -- sparse ----------------------------------------------------------

##-------- originally from ./dgCMatrix.R --------------------

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  signature(e1 = "dgCMatrix", e2 = "dgCMatrix"),
	  function(e1, e2) {
	      d <- dimCheck(e1, e2)
	      dn <- dimNamesCheck(e1, e2)
	      ij1 <- non0ind(e1)
	      ij2 <- non0ind(e2)
	      switch(.Generic,
		     "+" = , "-" =
		     ## special "T" convention: repeated entries are *summed*
		     as(new("dgTMatrix", Dim = d, Dimnames = dn,
			    i = c(ij1[,1], ij2[,1]),
			    j = c(ij1[,2], ij2[,2]),
			    x = c(callGeneric(e1@x, 0), callGeneric(0,e2@x))),
			"dgCMatrix"),

		     "*" =
		 { ##  X * 0 == 0 * X == 0 --> keep common non-0
		     ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
		     ij <- ij1[ii[[1]], , drop = FALSE]
		     as(new("dgTMatrix", Dim = d, Dimnames = dn,
			    i = ij[,1],
			    j = ij[,2],
			    x = e1@x[ii[[1]]] * e2@x[ii[[2]]]),
			"dgCMatrix")
		 },

		     "^" =
		 {
		     ii <- WhichintersectInd(ij1, ij2, nrow=d[1])
		     ## 3 cases:
		     ## 1) X^0 := 1  (even for X=0) ==> dense
		     ## 2) 0^Y := 0  for Y != 0		=====
		     ## 3) x^y :

		     ## FIXME:	dgeM[cbind(i,j)] <- V  is not yet possible
		     ##	    nor dgeM[ i_vect   ] <- V
		     ## r <- as(e2, "dgeMatrix")
		     ## ...
		     r <- as(e2, "matrix")
		     Yis0 <- is0(r)
		     r[complementInd(ij1, dim=d)] <- 0	    ## 2)
		     r[1:1 + ij2[ii[[2]], , drop=FALSE]] <-
			 e1@x[ii[[1]]] ^ e2@x[ii[[2]]]	    ## 3)
		     r[Yis0] <- 1			    ## 1)
		     as(r, "dgeMatrix")
		 },

		     "%%" = , "%/%" = , "/" = ## 0 op 0	 |-> NaN => dense
		     callGeneric(as(e1, "dgeMatrix"), e2)
		     )
	  })

setMethod("Arith",
	  signature(e1 = "dgCMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) { ## e.g.,  Mat ^ a
		  f0 <- callGeneric(0, e2)
		  if(is0(f0)) { # remain sparse
		      e1@x <- callGeneric(e1@x, e2)
		      e1
		  } else { ## non-sparse, since '0 o e2' is not 0

		      ## FIXME: dgeMatrix [cbind(i,j)] <- .. is not yet possible
		      ##		  r <- as(e1, "dgeMatrix")
		      ##		  r[] <- f0
		      ##		  r[non0ind(e1)] <- callGeneric(e1@x, e2)
		      r <- as(e1, "matrix")
		      r[] <- f0
		      r[non0ind(e1) + 1:1] <- callGeneric(e1@x, e2)
		      as(r, "dgeMatrix")
		  }
	      } else {
		  ## FIXME: maybe far from optimal:
		  warning("coercing sparse to dense matrix for arithmetic")
		  callGeneric(as(e1, "dgeMatrix"), e2)
	      }
	  })

setMethod("Arith",
	  signature(e1 = "numeric", e2 = "dgCMatrix"),
	  function(e1, e2) {
	      if(length(e1) == 1) {
		  f0 <- callGeneric(e1, 0)
		  if(is0(f0)) { # stay sparse, even "dgC"
		      e2@x <- callGeneric(e1, e2@x)
		      e2
		  } else {
		      ## FIXME: dgeMatrix [cbind(i,j)] <- .. is not yet possible
		      r <- as(e2, "matrix")
		      r[] <- f0
		      r[non0ind(e2) + 1:1] <- callGeneric(e1, e2@x)
		      as(r, "dgeMatrix")
		  }
	      } else {
		  ## FIXME: maybe far from optimal:
		  warning("coercing sparse to dense matrix for arithmetic")
		  callGeneric(e1, as(e2, "dgeMatrix"))
	      }
	  })

##-------- originally from ./Csparse.R --------------------
## TODO : Consider going a level up, and do this for all "Ops"

setMethod("Arith",
	  signature(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "dgCMatrix"),
				       as(e2, "dgCMatrix")))

setMethod("Arith",
	  signature(e1 = "CsparseMatrix", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e2) == 1) { ## e.g.,  Mat ^ a
		  f0 <- callGeneric(0, e2)
		  if(is0(f0)) { # remain sparse, symm., tri.,...
		      e1@x <- callGeneric(e1@x, e2)
		      return(e1)
		  }
	      }
	      ## all other (potentially non-sparse) cases: give up symm, tri,..
	      callGeneric(as(e1, paste(.M.kind(e1), "gCMatrix", sep='')), e2)
	  })

setMethod("Compare", signature(e1 = "CsparseMatrix", e2 = "CsparseMatrix"),
	  function(e1, e2) {
	      d <- dimCheck(e1,e2)

	      ## How do the "0" or "FALSE" entries compare?
	      ## Depends if we have an "EQuality RELation" or not:
	      EQrel <- switch(.Generic,
			      "==" =, "<=" =, ">=" = TRUE,
			      "!=" =, "<"  =, ">"  = FALSE)
	      if(EQrel) {
		  ## The (0 op 0) or  (FALSE op FALSE) comparison gives TRUE
		  ## -> result becomes *dense*; the following may be suboptimal
		  return( callGeneric(as(e1, "denseMatrix"),
				      as(e2, "denseMatrix")))
	      }

	      ## else: INequality:   0 op 0 gives FALSE ---> remain sparse!

	      ## NB non-diagonalMatrix := Union{ general, symmetric, triangular}
	      gen1 <- is(e1, "generalMatrix")
	      gen2 <- is(e2, "generalMatrix")
	      sym1 <- !gen1 && is(e1, "symmetricMatrix")
	      sym2 <- !gen2 && is(e2, "symmetricMatrix")
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
	      }
	      else if(!G && !S && !T) { ## coerce to generalMatrix and go
		  message("*** sparseMatrix comparison -- *unusual* case")
		  if(!gen1) e1 <- as(e1, "generalMatrix", strict = FALSE)
		  if(!gen2) e2 <- as(e2, "generalMatrix", strict = FALSE)
	      }

	      ## now the 'x' slots *should* match

	      newC <- sub("^.", "l", class(e1))
	      r <- new(newC)
	      r@x <- callGeneric(e1@x, e2@x)
	      for(sn in c("Dim", "Dimnames", "i", "p"))
		  slot(r, sn) <- slot(e1, sn)
	      r
	  })

## The same,  e1 <-> e2 :
setMethod("Arith",
	  signature(e1 = "numeric", e2 = "CsparseMatrix"),
	  function(e1, e2) {
	      if(length(e1) == 1) {
		  f0 <- callGeneric(e1, 0)
		  if(is0(f0)) {
		      e2@x <- callGeneric(e1, e2@x)
		      return(e2)
		  }
	      }
	      callGeneric(e1, as(e2, paste(.M.kind(e2), "gCMatrix", sep='')))
	  })

##-------- originally from ./sparseMatrix.R --------------------

## "Arith" short cuts / exceptions
setMethod("-", signature(e1 = "sparseMatrix", e2 = "missing"),
          function(e1) { e1@x <- -e1@x ; e1 })
## with the following exceptions:
setMethod("-", signature(e1 = "nsparseMatrix", e2 = "missing"),
          function(e1) callGeneric(as(e1, "dgCMatrix")))
setMethod("-", signature(e1 = "pMatrix", e2 = "missing"),
          function(e1) callGeneric(as(e1, "ngTMatrix")))

## Group method  "Arith"

## have CsparseMatrix methods (-> ./Csparse.R )
## which may preserve "symmetric", "triangular" -- simply defer to those:

setMethod("Arith", ##  "+", "-", "*", "^", "%%", "%/%", "/"
	  signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
				       as(e2, "CsparseMatrix")))
setMethod("Arith",
	  signature(e1 = "sparseMatrix", e2 = "numeric"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"), e2))
setMethod("Arith",
	  signature(e1 = "numeric", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(e1, as(e2, "CsparseMatrix")))

setMethod("Math",
	  signature(x = "sparseMatrix"),
	  function(x) callGeneric(as(x, "CsparseMatrix")))

setMethod("Compare", signature(e1 = "sparseMatrix", e2 = "sparseMatrix"),
	  function(e1, e2) callGeneric(as(e1, "CsparseMatrix"),
				       as(e2, "CsparseMatrix")))

