\name{expand}
\title{Expand a (Matrix) Decomposition into Factors}
%
\alias{expand}
\alias{expand2}
%
\alias{expand,CHMfactor-method}
\alias{expand,denseLU-method}
\alias{expand,sparseLU-method}
%
\alias{expand2,CHMfactor-method}
\alias{expand2,Cholesky-method}
\alias{expand2,Schur-method}
\alias{expand2,denseLU-method}
\alias{expand2,pCholesky-method}
\alias{expand2,sparseLU-method}
%
\description{
  Expands decompositions stored in compact form into factors.
}
\usage{
expand(x, \dots)
}
\arguments{
  \item{x}{a matrix decomposition.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  The expanded decomposition, typically a list of matrix factors.
}
\details{
  This is a generic function with special methods for different types
  of decompositions, see \code{\link{showMethods}(expand)} to list them
  all.
}
\note{
  Factors for decompositions such as \code{lu} and \code{qr} can be
  stored in a compact form.
  The function \code{expand} allows all factors to be fully expanded.
}
\seealso{
  The LU \code{\link{lu}}, and the \code{\link{Cholesky}}
  decompositions which have \code{expand} methods;
  \code{\link{facmul}}.
}
\examples{
(x <- Matrix(round(rnorm(9),2), 3, 3))
(ex <- expand(lux <- lu(x)))
}
\keyword{algebra}
