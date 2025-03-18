\name{packedMatrix-class}
\title{Virtual Class \sQuote{packedMatrix} of Packed Dense Format Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{packedMatrix-class}
%
\alias{nspMatrix-class}
\alias{lspMatrix-class}
\alias{ispMatrix-class}
\alias{dspMatrix-class}
\alias{zspMatrix-class}
%
\alias{dppMatrix-class}
\alias{zppMatrix-class}
%
\alias{ntpMatrix-class}
\alias{ltpMatrix-class}
\alias{itpMatrix-class}
\alias{dtpMatrix-class}
\alias{ztpMatrix-class}
%
\description{
  \code{packedMatrix} is a virtual subclass of
  \code{\linkS4class{denseMatrix}} representing packed dense format
  matrices, as defined by LAPACK:
  \url{https://netlib.org/lapack/lug/node121.html}.  To represent an
  \code{n}-by-\code{n} symmetric, Hermitian, or triangular matrix
  \code{A}, this format stores a vector of length
  \code{choose(n+1, 2) == n*(n+1)/2} (the number of entries in the upper
  or lower triangle), of which \code{n*(n+1)/2} or fewer elements are
  actually used (fewer if there are additional formal constraints on the
  structure of \code{A}).
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.  \code{Dim[1]} and \code{Dim[2]} must
      be equal as \code{A} is square by definition.}
    \item{\code{uplo}}{a character string, either \code{"U"} or
      \code{"L"}, indicating that \code{x} stores upper or lower
      triangular entries of \code{A}.}
    \item{\code{x}}{a vector of length \code{n*(n+1)/2},
      \code{n=Dim[1]}, storing upper or lower triangular (depending on
      \code{uplo}) entries of \code{A} in column-major order.
      Authors of subclasses of \code{packedMatrix} may declare that only
      a subset \code{k} of the elements of \code{x} is used to determine
      \code{A}, defining implicitly a function \code{F} such that
      \code{A = F(x[k])}; see \sQuote{Subclasses} for examples in
      package \pkg{Matrix}.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{packedMatrix} defined in package
  \pkg{Matrix} are listed below with links to virtual superclasses
  defining their data type and structure.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Structure\cr
    \code{nspMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{lspMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{ispMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dspMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{zspMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dppMatrix} \tab \code{\linkS4class{dMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{zppMatrix} \tab \code{\linkS4class{zMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{ntpMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ltpMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{itpMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{dtpMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ztpMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{triangularMatrix}}
  }
  The subclasses of \code{iMatrix} and \code{zMatrix} were not defined
  until \pkg{Matrix} version 1.8-0.  The other subclasses predate
  \code{packedMatrix} itself, which was not defined until \pkg{Matrix}
  version 1.4-1; previously, those subclasses were direct extensions of
  \code{denseMatrix}.

  The mapping from the \code{x} slot of length \code{n*(n+1)/2} to the
  \code{n}-by-\code{n} represented matrix \code{A} is given by:
\preformatted{
                                                                        A[i, j]                                        
                               |            i < j           |            i = j            |            i > j           
-------------------------------|----------------------------|-----------------------------|----------------------------
.spMatrix, uplo="U"            |      x[i+      j*(j-1)/2]  |    x[i+      j*(j-1)/2]     |      x[j+      i*(i-1)/2]  
zspMatrix, uplo="U", trans="C" |      x[i+      j*(j-1)/2]  | Re(x[i+      j*(j-1)/2])+0i | Conj(x[j+      i*(i-1)/2]) 
.spMatrix, uplo="L"            |      x[j+(2*n-i)*(i-1)/2]  |    x[i+(2*n-j)*(j-1)/2]     |      x[i+(2*n-j)*(j-1)/2]  
zspMatrix, uplo="L", trans="C" | Conj(x[j+(2*n-i)*(i-1)/2]) | Re(x[i+(2*n-j)*(j-1)/2])+0i |      x[i+(2*n-j)*(j-1)/2]  
.tpMatrix, uplo="U",  diag="N" |      x[i+      j*(j-1)/2]  |    x[i+      j*(j-1)/2]     |                         0  
.tpMatrix, uplo="U",  diag="U" |      x[i+      j*(j-1)/2]  |                       1     |                         0  
.tpMatrix, uplo="L",  diag="N" |                         0  |    x[i+(2*n-j)*(j-1)/2]     |      x[i+(2*n-j)*(j-1)/2]  
.tpMatrix, uplo="L",  diag="U" |                         0  |                       1     |      x[i+(2*n-j)*(j-1)/2]  
}

  Classes \code{.spMatrix} represent Hermitian or symmetric \code{A}.
  Objects have \code{x} and \code{uplo} slots.  If \code{uplo="U"}, then
  \code{x} stores the \code{n*(n+1)/2} upper triangular entries of
  \code{A}.  If \code{uplo="L"}, then \code{x} stores the
  \code{n*(n+1)/2} lower triangular entries of \code{A}.  Objects of
  class \code{zspMatrix} have an additional \code{trans} slot that
  distinguishes Hermitian and symmetric \code{A}, which use
  \code{trans="C"} and \code{"T"}, respectively.  If \code{trans="C"},
  then the imaginary parts of \code{x[du]} (\code{uplo="U"}) or
  \code{x[dl]} (\code{uplo="L"}) are not referenced, where \code{du=d},
  \code{dl=length(x)+1-d}, and \code{d=cumsum(seq_len(n))}.  (The
  diagonal entries of Hermitian matrices are real by definition.)

  Classes \code{dppMatrix} and \code{zppMatrix} extend \code{dspMatrix}
  and \code{zspMatrix} and represent positive semidefinite \code{A}.
  Positive semidefiniteness of \code{A} induces a constraint on the
  \code{x} slot.  However, this constraint is not strictly enforced by
  the validity methods, which test only that the diagonal entries of
  \code{A} are nonnegative (a necessary, not sufficient condition).
  For \code{zppMatrix}, only \code{trans="C"} is valid.

  Classes \code{.tpMatrix} represent triangular \code{A}.  Objects have
  \code{x}, \code{uplo}, and \code{diag} slots.  If \code{uplo="U"},
  then \code{A} is upper triangular and \code{x[-du]} stores the
  \code{n*(n-1)/2} strictly upper triangular entries.  If
  \code{uplo="L"}, then \code{A} is lower triangular and \code{x[-dl]}
  stores the \code{n*(n-1)/2} strictly lower triangular entries.  If
  \code{diag="N"}, then \code{x[du]} (\code{uplo="U"}) or \code{x[dl]}
  (\code{uplo="L"}) stores the \code{n} diagonal entries of \code{A}.
  If \code{diag="U"}, then \code{A} is unit triangular and those
  elements of \code{x} are unreferenced.
}
\seealso{
  \dQuote{Complementary} class \code{\linkS4class{unpackedMatrix}}
  representing conventional (or \dQuote{unpacked}) dense format
  matrices.  Generic function \code{\linkS4methods{unpack}}, an
  alternative to explicit coercion to \code{unpackedMatrix}.
}
\examples{
showClass("packedMatrix")
showMethods(classes = "packedMatrix")

%% moved from ./nsyMatrix-class.Rd :
(s0 <- new("nsyMatrix"))

(M2 <- Matrix(c(TRUE, NA, FALSE, FALSE), 2, 2)) # logical dense (ltr)
(sM <- M2 & t(M2))                       # -> "lge"
class(sM <- as(sM, "nMatrix"))           # -> "nge"
     (sM <- as(sM, "symmetricMatrix"))   # -> "nsy"
str(sM <- as(sM, "packedMatrix")) # -> "nsp", i.e., packed symmetric

%% moved from ./lsyMatrix-class.Rd :
(M2 <- Matrix(c(TRUE, NA, FALSE, FALSE), 2, 2)) # logical dense (ltr)
str(M2)
# can
(sM <- M2 | t(M2)) # "lge"
as(sM, "symmetricMatrix")
str(sM <- as(sM, "packedMatrix")) # packed symmetric

%% moved from ./dsyMatrix-class.Rd :
## Only upper triangular part matters (when uplo == "U" as per default)
(sy2 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, NA,32,77)))
str(t(sy2)) # uplo = "L", and the lower tri. (i.e. NA is replaced).

chol(sy2) #-> "Cholesky" matrix
(sp2 <- pack(sy2)) # a "dspMatrix"

## Coercing to dpoMatrix gives invalid object:
sy3 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, -1, 2, -7))
try(as(sy3, "dpoMatrix")) # -> error: not positive definite
\dontshow{
tr <- try(as(sy3, "dpoMatrix"), silent=TRUE)
stopifnot(1 == grep("not a positive definite matrix", as.character(tr)),
          is(sp2, "dspMatrix"))
}

## 4x4 example
m <- matrix(0,4,4); m[upper.tri(m)] <- 1:6
(sym <- m+t(m)+diag(11:14, 4))
(S1 <- pack(sym))
(S2 <- t(S1))
stopifnot(all(S1 == S2)) # equal "seen as matrix", but differ internally :
str(S1)
S2@x

%% moved from ./ntrMatrix-class.Rd :
showClass("ntrMatrix")

str(new("ntpMatrix"))
(nutr <- as(upper.tri(matrix(, 4, 4)), "ndenseMatrix"))
str(nutp <- pack(nutr)) # packed matrix: only 10 = 4*(4+1)/2 entries
!nutp # the logical negation (is *not* logical triangular !)
## but this one is:
stopifnot(all.equal(nutp, pack(!!nutp)))

%% moved from ./ltrMatrix-class.Rd :
showClass("ltrMatrix")

str(new("ltpMatrix"))
(lutr <- as(upper.tri(matrix(, 4, 4)), "ldenseMatrix"))
str(lutp <- pack(lutr)) # packed matrix: only 10 = 4*(4+1)/2 entries
!lutp # the logical negation (is *not* logical triangular !)
## but this one is:
stopifnot(all.equal(lutp, pack(!!lutp)))

%% moved from ./dtrMatrix-class-dense.Rd :
(m <- rbind(2:3, 0:-1)+0)
(M <- as(m, "generalMatrix"))

(T <- as(M, "triangularMatrix")) # formally upper triangular
(T2 <- as(t(M), "triangularMatrix"))
stopifnot(T@uplo == "U", T2@uplo == "L", identical(T2, t(T)))

m <- matrix(0,4,4); m[upper.tri(m)] <- 1:6
(t1 <- Matrix(m+diag(,4)))
str(t1p <- pack(t1))
(t1pu <- diagN2U(t1p))
stopifnot(exprs = {
   inherits(t1 , "dtrMatrix"); validObject(t1)
   inherits(t1p, "dtpMatrix"); validObject(t1p)
   inherits(t1pu,"dtCMatrix"); validObject(t1pu)
   t1pu@x == 1:6
   all(t1pu == t1p)
   identical((t1pu - t1)@x, numeric())# sparse all-0
})

%% moved from ./dtpMatrix-class.Rd
showClass("dtpMatrix")

(p1 <- pack(T2))
str(p1)
(pp <- pack(T))
ip1 <- solve(p1)
stopifnot(length(p1@x) == 3, length(pp@x) == 3,
          p1 @ uplo == T2 @ uplo, pp @ uplo == T @ uplo,
	  identical(t(pp), p1), identical(t(p1), pp),
	  all((l.d <- p1 - T2) == 0), is(l.d, "dtpMatrix"),
	  all((u.d <- pp - T ) == 0), is(u.d, "dtpMatrix"),
	  l.d@uplo == T2@uplo, u.d@uplo == T@uplo,
	  identical(t(ip1), solve(pp)), is(ip1, "dtpMatrix"),
	  all.equal(as(solve(p1,p1), "diagonalMatrix"), Diagonal(2)))
}
