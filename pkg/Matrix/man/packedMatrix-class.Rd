\name{packedMatrix-class}
\title{Virtual Class \code{"packedMatrix"} of Packed Dense Matrices}
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
\alias{ntpMatrix-class}
\alias{ltpMatrix-class}
\alias{itpMatrix-class}
\alias{dtpMatrix-class}
\alias{ztpMatrix-class}
%
\alias{coerce,matrix,packedMatrix-method}
%
\description{
  Class \code{"packedMatrix"} is the \emph{virtual} class of dense
  symmetric or triangular matrices in "packed" format, storing only
  the \code{choose(n+1,2) == n*(n+1)/2} elements of the upper or
  lower triangle of an \code{n}-by-\code{n} matrix. It is used to
  define common methods for efficient subsetting, transposing, etc.
  of its \emph{proper} subclasses: currently \code{"[dln]spMatrix"}
  (packed symmetric), \code{"[dln]tpMatrix"} (packed triangular),
  and subclasses of these, such as
  \code{"\linkS4class{dppMatrix}"}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{\code{"character"}; either "U", for upper triangular, and "L", for lower.}
    \item{\code{Dim}, \code{Dimnames}:}{as all \code{\linkS4class{Matrix}} objects.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{denseMatrix}"}, directly.
  Class \code{"\linkS4class{Matrix}"}, by class \code{"denseMatrix"},
  distance 2.
}
\section{Methods}{
  \describe{
    \item{pack}{\code{signature(x = "packedMatrix")}: ... }
    \item{unpack}{\code{signature(x = "packedMatrix")}: ... }
    \item{isSymmetric}{\code{signature(object = "packedMatrix")}: ... }
    \item{isTriangular}{\code{signature(object = "packedMatrix")}: ... }
    \item{isDiagonal}{\code{signature(object = "packedMatrix")}: ... }
    \item{t}{\code{signature(x = "packedMatrix")}: ... }
    \item{diag}{\code{signature(x = "packedMatrix")}: ... }
    \item{diag<-}{\code{signature(x = "packedMatrix")}: ... }
  }
}
%% \references{
%% }
\author{Mikael Jagan}
%% \note{
%% }
\seealso{
  \code{\link{pack}} and \code{\link{unpack}}; its virtual "complement"
  \code{"\linkS4class{unpackedMatrix}"}; its proper subclasses
  \code{"\linkS4class{dspMatrix}"}, \code{"\linkS4class{ltpMatrix}"}, etc.
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
stopifnot(1 == grep("not a positive definite matrix",
                    as.character(tr)),
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
