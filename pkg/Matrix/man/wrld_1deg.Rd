\name{wrld_1deg}
\title{World 1-degree grid contiguity matrix}
%
\docType{data}
\keyword{datasets}
%
\alias{wrld_1deg}
%
\description{
  This matrix represents the distance-based contiguities of 15260
  one-degree grid cells of land areas. The representation is as a
  row standardised spatial weights matrix transformed to a symmetric
  matrix (see Ord (1975), p. 125).
}
\usage{data(wrld_1deg)}
\format{
  A \eqn{15260 ^2} symmetric sparse matrix of class
  \code{\linkS4class{dsCMatrix}} with 55973 non-zero entries.
}
\details{
  The data were created into \R using the coordinates of a
  \sQuote{SpatialPixels} object containing approximately one-degree grid
  cells for land areas only (world excluding Antarctica), using package
  \CRANpkg{spdep}'s
  \code{\link[spdep]{dnearneigh}} with a cutoff distance of \code{sqrt(2)},
  and row-standardised and transformed to symmetry using
  \code{\link[spdep]{nb2listw}} and \code{\link[spdep]{similar.listw}}.
  This spatial weights object was converted to a \code{\linkS4class{dsTMatrix}}
  using \code{\link[spdep]{as_dsTMatrix_listw}} and then coerced
  (column-compressed).
}
\source{
  The shoreline data was read into \R using \code{maptools::Rgshhs}
  from the GSHHS coarse shoreline database distributed with the
  former maptools package, omitting Antarctica.  A matching approximately
  one-degree grid was generated using
  \code{maptools::Sobj_SpatialGrid}, and the grids on land were
  found using the appropriate \code{\link[sp]{over}}% was {overlay} till 2016-02
  method for the \sQuote{SpatialPolygons} and \sQuote{SpatialGrid} objects,
  yielding a \sQuote{SpatialPixels} one containing only the grid cells with
  centres on land.
}
\references{
  Ord, J. K. (1975)
  Estimation methods for models of spatial interaction;
  \emph{Journal of the American Statistical Association} \bold{70}, 120--126.
}
\examples{
data(wrld_1deg)
(n <- ncol(wrld_1deg))
I <- .symDiagonal(n)

doExtras <- interactive() || nzchar(Sys.getenv("R_MATRIX_CHECK_EXTRA"))
set.seed(1)
r <- if(doExtras) 20L else 3L
rho <- 1 / runif(r, 0, 0.5)

system.time(MJ0 <- sapply(rho, function(mult)
    determinant(wrld_1deg + mult * I, logarithm = TRUE)$modulus))

## Can be done faster by updating the Cholesky factor:

C1 <- Cholesky(wrld_1deg, Imult = 2)
system.time(MJ1 <- sapply(rho, function(mult)
    c(determinant(update(C1, wrld_1deg, mult))$modulus)))
stopifnot(all.equal(MJ0, MJ1))

C2 <- Cholesky(wrld_1deg, super = TRUE, Imult = 2)
system.time(MJ2 <- sapply(rho, function(mult)
    c(determinant(update(C2, wrld_1deg, mult))$modulus)))
stopifnot(all.equal(MJ0, MJ2))
}
