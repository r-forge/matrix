\name{USCounties}
\title{USCounties Contiguity Matrix}
%
\docType{data}
\keyword{datasets}
%
\alias{USCounties}
%
\description{
  This matrix represents the contiguities of 3111 US counties using the
  Queen criterion of at least a single shared boundary point. The
  representation is as a row standardised spatial weights matrix
  transformed to a symmetric matrix (see Ord (1975), p. 125).
}
\usage{data(USCounties)}
\format{
  A \eqn{3111 ^2} symmetric sparse matrix of class
  \code{\linkS4class{dsCMatrix}} with 9101 non-zero entries.
}
\details{
  The data were read into \R using \code{\link[spdep]{read.gal}}, and
  row-standardised and transformed to symmetry using
  \code{\link[spdep]{nb2listw}} and \code{\link[spdep]{similar.listw}}.
  This spatial weights object was converted to class
  \code{\linkS4class{dsCMatrix}} using
  \code{\link[spdep]{as_dsTMatrix_listw}} and coercion.
}
\source{ %% no \url anymore, as broken (Jan. 2015)
  The data were retrieved from \code{http://sal.uiuc.edu/weights/zips/usc.zip},
  files \dQuote{usc.txt} and \dQuote{usc_q.GAL}, with permission for use
  and distribution from Luc Anselin (in early 2008).
}
\references{
  Ord, J. K. (1975)
  Estimation methods for models of spatial interaction;
  \emph{Journal of the American Statistical Association} \bold{70}, 120--126.
}
\examples{
data(USCounties)
(n <- ncol(USCounties))
I <- .symDiagonal(n)

set.seed(1)
r <- 50L
rho <- 1 / runif(r, 0, 0.5)

system.time(MJ0 <- sapply(rho, function(mult)
    determinant(USCounties + mult * I, logarithm = TRUE)$modulus))

## Can be done faster by updating the Cholesky factor:

C1 <- Cholesky(USCounties, Imult = 2)
system.time(MJ1 <- sapply(rho, function(mult)
    c(determinant(update(C1, USCounties, mult))$modulus)))
stopifnot(all.equal(MJ0, MJ1))

C2 <- Cholesky(USCounties, super = TRUE, Imult = 2)
system.time(MJ2 <- sapply(rho, function(mult)
    c(determinant(update(C2, USCounties, mult))$modulus)))
stopifnot(all.equal(MJ0, MJ2))
}
