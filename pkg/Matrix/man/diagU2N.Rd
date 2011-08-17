\name{diagU2N}
\alias{diagU2N}
\alias{.diagU2N}
\title{Transform Triangular Matrices from Unidiagonal to General and back}
\description{
  Transform a triangular matrix \code{x}, i.e., of \code{\link{class}}
  \code{"\linkS4class{triangularMatrix}"},
  from \dQuote{unidiagonal} to \dQuote{general} (\code{diagU2N(x)}) or
  back (\code{diagN2U(x)}).
}
\usage{
 diagN2U(x, cl = getClassDef(class(x)))
 diagU2N(x, cl = getClassDef(class(x)), checkDense = FALSE)
.diagU2N(x, cl, checkDense = FALSE)
}
\arguments{
  \item{x}{a \code{\linkS4class{triangularMatrix}}, often sparse.}
  \item{cl}{(optional, for speedup only:) class (definition) of \code{x}.}
  \item{checkDense}{logical indicating if dense (see
    \code{\linkS4class{denseMatrix}}) matrices should be considered at all.}
}
% \details{
% }
\value{
  a triangular matrix of the same \code{\link{class}} as \code{x},
  but with a different \code{diag} slot.
}
\seealso{
  \code{"\linkS4class{triangularMatrix}"},
  \code{"\linkS4class{dtCMatrix}"}.
}
\examples{
(T <- Diagonal(7) + triu(Matrix(rpois(49, 1/4), 7,7), k = 1))
(uT <- diagN2U(T)) # "unitriangular"
stopifnot(all(T == uT),
          identical(T, diagU2N(uT)))
}
\keyword{utilities}
\keyword{classes}
