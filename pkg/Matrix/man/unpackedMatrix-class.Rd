\name{unpackedMatrix-class}
\title{Virtual Class \sQuote{unpackedMatrix} of Conventional Dense Format Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{unpackedMatrix-class}
%
\alias{ngeMatrix-class}
\alias{lgeMatrix-class}
\alias{igeMatrix-class}
\alias{dgeMatrix-class}
\alias{zgeMatrix-class}
%
\alias{nsyMatrix-class}
\alias{lsyMatrix-class}
\alias{isyMatrix-class}
\alias{dsyMatrix-class}
\alias{zsyMatrix-class}
%
\alias{dpoMatrix-class}
\alias{zpoMatrix-class}
%
\alias{ntrMatrix-class}
\alias{ltrMatrix-class}
\alias{itrMatrix-class}
\alias{dtrMatrix-class}
\alias{ztrMatrix-class}
%
\description{
  \code{unpackedMatrix} is a virtual subclass of
  \code{\linkS4class{denseMatrix}} representing conventional dense
  format matrices, as defined by LAPACK:
  \url{https://netlib.org/lapack/lug/node121.html}.  To represent an
  \code{m}-by-\code{n} matrix \code{A}, this format stores a vector of
  length \code{m*n}, of which \code{m*n} \emph{or fewer} elements are
  actually used (fewer if there are formal constraints on the structure
  of \code{A}).  Within package \pkg{Matrix}, this format is called
  \dQuote{unpacked} to emphasize that it is complementary to the packed
  format also defined by LAPACK.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{x}}{a vector of length \code{prod(Dim)} storing entries
      of \code{A} in column-major order.  Authors of subclasses of
      \code{unpackedMatrix} may declare that only a subset \code{k} of
      the elements of \code{x} is used to determine \code{A}, defining
      implicitly a function \code{F} such that \code{A = F(x[k])}; see
      \sQuote{Subclasses} for examples in package \pkg{Matrix}.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{unpackedMatrix} defined in package
  \pkg{Matrix} are listed below with links to virtual superclasses
  defining their data type and structure.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Structure\cr
    \code{ngeMatrix} \tab \code{\linkS4class{nMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{lgeMatrix} \tab \code{\linkS4class{lMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{igeMatrix} \tab \code{\linkS4class{iMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{dgeMatrix} \tab \code{\linkS4class{dMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{zgeMatrix} \tab \code{\linkS4class{zMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{nsyMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{lsyMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{isyMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dsyMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{zsyMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dpoMatrix} \tab \code{\linkS4class{dMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{zpoMatrix} \tab \code{\linkS4class{zMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{ntrMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ltrMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{itrMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{dtrMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ztrMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{triangularMatrix}}
  }
  The subclasses of \code{iMatrix} and \code{zMatrix} were not defined
  until \pkg{Matrix} version 1.8-0.  The other subclasses predate
  \code{unpackedMatrix} itself, which was not defined until \pkg{Matrix}
  version 1.5-0; previously, those subclasses were direct extensions of
  \code{denseMatrix}.

  The mapping from the \code{x} slot of length \code{m*n} to the
  \code{m}-by-\code{n} represented matrix \code{A} is given by:
\preformatted{
                                                             A[i, j]                            
                               |        i < j       |         i = j        |        i > j       
                               |--------------------|----------------------|--------------------
.geMatrix                      |      x[(j-1)*m+i]  |     x[(j-1)*m+i]     |      x[(j-1)*m+i]  
.syMatrix, uplo="U"            |      x[(j-1)*m+i]  |     x[(j-1)*m+i]     |      x[(i-1)*m+j]  
zsyMatrix, uplo="U", trans="C" |      x[(j-1)*m+i]  |  Re(x[(j-1)*m+i])+0i | Conj(x[(i-1)*m+j]) 
.syMatrix, uplo="L"            |      x[(i-1)*m+j]  |     x[(j-1)*m+i]     |      x[(j-1)*m+i]  
zsyMatrix, uplo="L", trans="C" | Conj(x[(i-1)*m+j]) |  Re(x[(j-1)*m+i])+0i |      x[(j-1)*m+i]  
.trMatrix, uplo="U",  diag="N" |      x[(j-1)*m+i]  |     x[(j-1)*m+i]     |                 0  
.trMatrix, uplo="U",  diag="U" |      x[(j-1)*m+i]  |                1     |                 0  
.trMatrix, uplo="L",  diag="N" |                 0  |     x[(j-1)*m+i]     |      x[(j-1)*m+j]  
.trMatrix, uplo="L",  diag="U" |                 0  |                1     |      x[(j-1)*m+j]  
}

  Classes \code{.geMatrix} represent general (that is, unstructured)
  \code{A}.  Objects have an \code{x} slot that stores all \code{m*n}
  entries of \code{A}.  These classes closely approximate traditional
  matrices.  For example, a traditional matrix \code{z} of type
  \code{integer} and its analogue \code{Z} of formal class
  \code{igeMatrix} are related by:
  \preformatted{	z <- array(Z@x, dim = Z@Dim, dimnames = Z@Dimnames)
	Z <- new("igeMatrix", Dim = dim(z), Dimnames = dimnames(z), x = as.vector(z))}

  Classes \code{.syMatrix} represent Hermitian or symmetric \code{A}.
  Objects have \code{x} and \code{uplo} slots.  If \code{uplo="U"}, then
  \code{x[row(A) <= col(A)]} stores the \code{n*(n+1)/2} upper
  triangular entries of \code{A} and \code{x[row(A) > col(A)]} is
  unreferenced.  If \code{uplo="L"}, then \code{x[row(A) >= col(A)]}
  stores the \code{n*(n+1)/2} lower triangular entries of \code{A} and
  \code{x[row(A) < col(A)]} is unreferenced.  Objects of class
  \code{zsyMatrix} have an additional \code{trans} slot that
  distinguishes Hermitian and symmetric \code{A}, which use
  \code{trans="C"} and \code{"T"}, respectively.  If \code{trans="C"},
  then the imaginary parts of \code{x[row(A) == col(A)]} are not
  referenced.  (The diagonal entries of Hermitian matrices are real by
  definition.)

  Classes \code{dpoMatrix} and \code{zpoMatrix} extend \code{dsyMatrix}
  and \code{zsyMatrix} and represent positive semidefinite \code{A}.
  Positive semidefiniteness of \code{A} induces a constraint on the
  \code{x} slot.  However, this constraint is not strictly enforced by
  the validity methods, which test only that the diagonal entries of
  \code{A} are nonnegative (a necessary, not sufficient condition).
  For \code{zpoMatrix}, only \code{trans="C"} is valid.

  Classes \code{.trMatrix} represent triangular \code{A}.  Objects have
  \code{x}, \code{uplo}, and \code{diag} slots.  If \code{uplo="U"},
  then \code{A} is upper triangular, \code{x[row(A) < col(A)]} stores
  the \code{n*(n-1)/2} strictly upper triangular entries, and
  \code{x[row(A) > col(A)]} is unreferenced.  If \code{uplo="L"}, then
  \code{A} is lower triangular, \code{x[row(A) > col(A)]} stores the
  \code{n*(n-1)/2} strictly lower triangular entries, and
  \code{x[row(A) < col(A)]} is unreferenced.  If \code{diag="N"}, then
  \code{x[row(A) == col(A)]} stores the \code{n} diagonal entries of
  \code{A}.  If \code{diag="U"}, then \code{A} is unit triangular and
  those elements of \code{x} are unreferenced.
}
\seealso{
  \dQuote{Complementary} class \code{\linkS4class{packedMatrix}}
  representing packed dense format matrices.  Generic function
  \code{\linkS4methods{pack}}, an alternative to explicit coercion to
  \code{packedMatrix}.
}
\examples{
showClass("unpackedMatrix")
showMethods(classes = "unpackedMatrix")

%% moved from ./ngeMatrix-class.Rd :
showClass("ngeMatrix")
## "lgeMatrix" is really more relevant

%% moved from ./lgeMatrix-class.Rd :
showClass("lgeMatrix")
str(new("lgeMatrix"))
set.seed(1)
(lM <- Matrix(matrix(rnorm(28), 4,7) > 0))# a simple random lgeMatrix
set.seed(11)
(lC <- Matrix(matrix(rnorm(28), 4,7) > 0))# a simple random lgCMatrix
as(lM, "CsparseMatrix")

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

%% moved from ./dtrMatrix-class.Rd :
(m <- rbind(2:3, 0:-1))
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
}
