\name{diagU2N}
\title{Transform Triangular/Diagonal Matrices from Unit to General Triangular/Diagonal and Back}
%
\keyword{array}
\keyword{attribute}
\keyword{utilities}
%
\alias{diagU2N}
\alias{diagN2U}
\alias{.diagU2N}
\alias{.diagN2U}
%
\description{
  Transform a triangular or diagonal matrix \code{x}, i.e., of \code{\link{class}}
  \code{\linkS4class{triangularMatrix}} or \code{\linkS4class{diagonalMatrix}},
  from (internally!) unit triangular (\dQuote{unitriangular}) to
  \dQuote{general} triangular or diagonal, respectively, (\code{diagU2N(x)})
  or back (\code{diagN2U(x)}).
  Note that the latter, \code{diagN2U(x)}, also sets the diagonal to one
  in cases where \code{diag(x)} was not all one.

  \code{.diagU2N(x)} and \code{.diagN2U(x)} assume \emph{without}
  checking that \code{x} is a \code{\linkS4class{triangularMatrix}} with
  suitable \code{diag} slot (\code{"U"} and \code{"N"}, respectively),
  hence they should be used with care.
}
\usage{
 diagU2N(x, cl = getClassDef(class(x)), checkDense = TRUE)
 diagN2U(x, cl = getClassDef(class(x)), checkDense = TRUE)

.diagU2N(x, cl = getClassDef(class(x)), checkDense = FALSE)
.diagN2U(x, cl = getClassDef(class(x)), checkDense = FALSE)
}
\arguments{
  \item{x}{a \code{\linkS4class{triangularMatrix}}, often sparse.}
  \item{cl}{(optional, for speedup only:) class (definition) of \code{x}.}
  \item{checkDense}{logical indicating if dense (see
    \code{\linkS4class{denseMatrix}}) matrices should be considered at
    all.  I.e., when false, the result will be sparse even when \code{x} is dense.}
}
\details{
  The concept of unit triangular matrices with a \code{diag} slot of
  \code{"U"} stems from LAPACK.
}
\note{Such internal storage details should rarely be of relevance to the
  user.  Hence, these functions really are rather \emph{internal}
  utilities.
}
\value{
  a triangular matrix of the same \code{\link{class}} but with a
  different \code{diag} slot.  For \code{diagU2N} (semantically) with
  identical entries as \code{x}, whereas in \code{diagN2U(x)}, the
  off-diagonal entries are unchanged and the diagonal is set to all
  \code{1} even if it was not previously.
}
\seealso{
  \code{"\linkS4class{diagonalMatrix}"},
  \code{"\linkS4class{triangularMatrix}"},
  \code{"\linkS4class{dtCMatrix}"}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL %
library(stats, pos = "package:base", verbose = FALSE)
}
D  <- Diagonal(7)
Dn <- diagU2N(D)
(iD <- Diagonal(7, 1L, int2dbl = FALSE))
iDn <- diagU2N(iD)
stopifnot(D @diag == "U", identical( D @x, numeric()),
          Dn@diag == "N", identical( Dn@x, rep(1,  7))
        ,iDn@diag == "N", identical(iDn@x, rep(1L, 7))
        , identical( D, diagN2U( Dn))
        , identical(iD, diagN2U(iDn))
          )

set.seed(64); (Tr <- Matrix(rpois(49, 1/4), 7, 7)) # sparse enough
(T <- D + triu(Tr, k = 1)) # sparse triangular
(uT <- diagN2U(T)) # "unitriangular"
(t.u <- diagN2U(10*T))# changes the diagonal!
stopifnot(all(T == uT), diag(t.u) == 1, T@x == 1,
          identical(T, diagU2N(uT)))
T[upper.tri(T)] <- 5 # still "dtC"
T <- diagN2U(as(T,"triangularMatrix"))
inherits(T, "sparseMatrix") # TRUE
## checkDense = T / F  :
dT <- as(T, "denseMatrix") # (unitriangular)
sT.n <-  diagU2N(dT, checkDense = TRUE)  # default
dT.n <- .diagU2N(dT, checkDense = FALSE) # default
stopifnot(is(dT.n, "denseMatrix"), is(sT.n, "sparseMatrix"),
          dT@diag == "U", dT.n@diag == "N", sT.n@diag == "N",
          all(dT == dT.n), all(dT == sT.n))
}
