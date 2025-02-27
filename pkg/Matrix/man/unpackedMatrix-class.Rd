\name{unpackedMatrix-class}
\title{Virtual Class \code{"unpackedMatrix"} of Unpacked Dense Matrices}
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
\alias{ntrMatrix-class}
\alias{ltrMatrix-class}
\alias{itrMatrix-class}
\alias{dtrMatrix-class}
\alias{ztrMatrix-class}
%
\description{
  Class \code{"unpackedMatrix"} is the \emph{virtual} class of dense
  matrices in "unpacked" format, storing all \code{m*n} elements of
  an \code{m}-by-\code{n} matrix. It is used to define common methods
  for efficient subsetting, transposing, etc. of its \emph{proper}
  subclasses: currently \code{"[dln]geMatrix"} (unpacked general),
  \code{"[dln]syMatrix"} (unpacked symmetric), \code{"[dln]trMatrix"}
  (unpacked triangular), and subclasses of these, such as
  \code{"\linkS4class{dpoMatrix}"}.
}
\section{Slots}{
  \describe{
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
    \item{pack}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{unpack}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{isSymmetric}{\code{signature(object = "unpackedMatrix")}: ... }
    \item{isTriangular}{\code{signature(object = "unpackedMatrix")}: ... }
    \item{isDiagonal}{\code{signature(object = "unpackedMatrix")}: ... }
    \item{t}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{diag}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{diag<-}{\code{signature(x = "unpackedMatrix")}: ... }
  }
}
%% \references{
%% }
\author{Mikael Jagan}
%% \note{
%% }
\seealso{
  \code{\link{pack}} and \code{\link{unpack}}; its virtual "complement"
  \code{"\linkS4class{packedMatrix}"}; its proper subclasses
  \code{"\linkS4class{dsyMatrix}"}, \code{"\linkS4class{ltrMatrix}"}, etc.
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

% moved from ./ntrMatrix-class.Rd :
showClass("ntrMatrix")

str(new("ntpMatrix"))
(nutr <- as(upper.tri(matrix(, 4, 4)), "ndenseMatrix"))
str(nutp <- pack(nutr)) # packed matrix: only 10 = 4*(4+1)/2 entries
!nutp # the logical negation (is *not* logical triangular !)
## but this one is:
stopifnot(all.equal(nutp, pack(!!nutp)))

% moved from ./ltrMatrix-class.Rd :
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
